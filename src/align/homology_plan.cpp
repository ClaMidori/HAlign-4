#include "align/homology_plan.hpp"

#include <algorithm>
#include <bit>
#include <cstdlib>
#include <array>
#include <cctype>
#include <fstream>
#include <limits>
#include <map>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <utility>

namespace homology_plan {
namespace {
struct Occurrence { std::uint32_t pos = 0; unsigned count = 0; };
using Index = std::unordered_map<std::uint64_t, Occurrence>;

void validateSequence(const Sequence& s) {
    if (s.id.empty() || s.id.find_first_of("\t\r\n") != std::string::npos)
        throw std::invalid_argument("FASTA identifier is empty or contains tabs/newlines");
    if (s.bases.empty()) throw std::invalid_argument("empty FASTA sequence: " + s.id);
    for (unsigned char c : s.bases) {
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T' && c != 'N')
            throw std::invalid_argument("homology planner accepts only A/C/G/T/N, not gaps or ambiguity codes: " + s.id);
    }
}
void validateInput(const std::vector<Sequence>& sequences, const Limits& limits) {
    if (sequences.size() < 2 || sequences.size() > limits.max_sequences)
        throw std::invalid_argument("diagnostic planner requires 2..8 sequences (bounded)");
    std::size_t total = 0;
    std::unordered_set<std::string> ids;
    for (const auto& s : sequences) {
        validateSequence(s);
        if (s.bases.size() > limits.max_per_sequence ||
            s.bases.size() > limits.max_total_bases - total)
            throw std::invalid_argument("planner input exceeds bounded 1Mb per sequence / 4Mb total budget");
        total += s.bases.size();
        if (!ids.insert(s.id).second)
            throw std::invalid_argument("duplicate FASTA identifier: " + s.id);
    }
}

// A rolling 2-bit 11-mer code. N and low-complexity seeds are excluded.
// Both strands are *not* merged; strand flips are unresolved in PR3D.1.
template<class Func>
void forKmers(const std::string& bases, unsigned k, Func fn) {
    const std::uint64_t mask = (std::uint64_t{1} << (2 * k)) - 1;
    std::uint64_t code = 0;
    unsigned valid = 0;
    for (std::size_t i = 0; i < bases.size(); ++i) {
        int b = -1;
        switch (bases[i]) {
            case 'A': b = 0; break; case 'C': b = 1; break;
            case 'G': b = 2; break; case 'T': b = 3; break;
            default: break;
        }
        if (b < 0) { code = 0; valid = 0; continue; }
        code = ((code << 2) | static_cast<std::uint64_t>(b)) & mask;
        if (valid < k) ++valid;
        if (valid < k) continue;
        std::uint64_t v = code;
        unsigned distinct = 0;
        for (unsigned j = 0; j < k; ++j) { distinct |= 1U << (v & 3U); v >>= 2; }
        if (std::popcount(distinct) < 3) continue;
        fn(code, static_cast<std::uint32_t>(i + 1 - k));
    }
}
Index makeIndex(const std::string& bases, unsigned k) {
    Index index;
    index.reserve(bases.size() / 2 + 1);
    forKmers(bases, k, [&](std::uint64_t code, std::uint32_t pos) {
        auto [it, inserted] = index.try_emplace(code, Occurrence{pos, 1});
        if (!inserted) it->second.count = 2; // 2 also represents >2
    });
    return index;
}
std::int64_t delta(const Anchor& a) {
    return static_cast<std::int64_t>(a.query) - a.reference;
}
std::int64_t floorDivide(std::int64_t value, std::int64_t divisor) {
    return value >= 0 ? value / divisor : -1 - (-value - 1) / divisor;
}

// A matching seed alone is insufficient: random and repeated backgrounds can
// produce dense diagonal coincidences, especially at megabase scale. Require
// sequence-context support in at least one direction, excluding the seed.
// This is NOT a local alignment and does not establish biological homology.
bool hasContext(const Anchor& a, const std::string& ref, const std::string& query,
                unsigned k) {
    constexpr std::size_t flank = 128;
    constexpr std::size_t minimum = 64;
    constexpr unsigned minimum_identity_percent = 42;
    const std::size_t r = a.reference;
    const std::size_t q = a.query;
    for (int direction : {-1, 1}) {
        unsigned matched = 0;
        std::size_t count = 0;
        for (std::size_t d = 1; d <= flank; ++d) {
            std::size_t rp = 0, qp = 0;
            if (direction < 0) {
                if (r < d || q < d) break;
                rp = r - d; qp = q - d;
            } else {
                if (r + k + d - 1 >= ref.size() || q + k + d - 1 >= query.size()) break;
                rp = r + k + d - 1; qp = q + k + d - 1;
            }
            ++count;
            if (ref[rp] == query[qp] && ref[rp] != 'N') ++matched;
        }
        if (count >= minimum && matched * 100U >= minimum_identity_percent * count)
            return true;
    }
    return false;
}

// Reject singleton matches and distant/repetitive diagonal coincidences.
// At least 4 nonoverlapping 11mers must be locally collinear on a narrow
// diagonal, then a deterministic LIS enforces global forward order.
std::vector<Anchor> evidenceFilter(std::vector<Anchor> candidates, const Config& c,
                                   const std::string& ref, const std::string& query) {
    candidates.erase(std::remove_if(candidates.begin(), candidates.end(),
        [&](const Anchor& a){return !hasContext(a, ref, query, c.k);}), candidates.end());
    std::sort(candidates.begin(), candidates.end(), [](const Anchor& a, const Anchor& b){
        return a.reference != b.reference ? a.reference < b.reference : a.query < b.query;
    });
    constexpr std::int64_t bin_width = 64;
    std::map<std::int64_t, std::vector<std::size_t>> bins;
    for (std::size_t i = 0; i < candidates.size(); ++i)
        bins[floorDivide(delta(candidates[i]), bin_width)].push_back(i);
    std::vector<Anchor> filtered;
    filtered.reserve(candidates.size());
    for (std::size_t i = 0; i < candidates.size(); ++i) {
        const auto& a = candidates[i];
        unsigned evidence = 1;
        std::array<std::size_t, 4> accepted = {i, i, i};
        const auto center_bin = floorDivide(delta(a), bin_width);
        for (std::int64_t b = center_bin - 1; b <= center_bin + 1 && evidence < c.min_independent_seeds; ++b) {
            const auto found = bins.find(b);
            if (found == bins.end()) continue;
            const auto& indexes = found->second;
            auto pivot = std::lower_bound(indexes.begin(), indexes.end(), a.reference,
                [&](std::size_t j, std::uint32_t pos){return candidates[j].reference < pos;});
            for (int direction : {-1, 1}) {
                auto iter = pivot;
                unsigned inspected = 0;
                while (inspected < 48 && evidence < c.min_independent_seeds) {
                    if (direction < 0) {
                        if (iter == indexes.begin()) break;
                        --iter;
                    } else {
                        if (iter == indexes.end()) break;
                    }
                    const auto j = *iter;
                    if (direction > 0) ++iter;
                    ++inspected;
                    if (j == i) continue;
                    const auto& z = candidates[j];
                    const auto dr = std::llabs(static_cast<std::int64_t>(z.reference) - a.reference);
                    if (dr > static_cast<std::int64_t>(c.local_span)) break;
                    if (dr < static_cast<std::int64_t>(c.k) ||
                        std::llabs(static_cast<std::int64_t>(z.query) - a.query) < static_cast<std::int64_t>(c.k) ||
                        std::llabs(delta(z) - delta(a)) > static_cast<std::int64_t>(c.diagonal_tolerance)) continue;
                    const bool forward = (z.reference < a.reference) == (z.query < a.query);
                    if (!forward) continue;
                    bool independent = true;
                    for (unsigned t = 1; t < evidence; ++t) {
                        const auto& prior = candidates[accepted[t]];
                        if (std::llabs(static_cast<std::int64_t>(z.reference) - prior.reference) < static_cast<std::int64_t>(c.k) ||
                            std::llabs(static_cast<std::int64_t>(z.query) - prior.query) < static_cast<std::int64_t>(c.k))
                            independent = false;
                    }
                    if (independent) accepted[evidence++] = j;
                }
            }
        }
        if (evidence >= c.min_independent_seeds) filtered.push_back(a);
    }
    // Strictly-increasing subsequence in query coordinates. Stable tie-breaking
    // uses ascending reference/query input and std::lower_bound replacement.
    if (filtered.empty()) return filtered;
    std::vector<std::size_t> tails, previous(filtered.size(), filtered.size());
    for (std::size_t i = 0; i < filtered.size(); ++i) {
        auto it = std::lower_bound(tails.begin(), tails.end(), filtered[i].query,
            [&](std::size_t j, std::uint32_t q){return filtered[j].query < q;});
        const auto position = static_cast<std::size_t>(it - tails.begin());
        if (position) previous[i] = tails[position - 1];
        if (it == tails.end()) tails.push_back(i);
        else *it = i;
    }
    std::vector<Anchor> chain;
    chain.reserve(tails.size());
    for (std::size_t p = tails.back(); p < filtered.size(); p = previous[p]) {
        chain.push_back(filtered[p]);
        if (previous[p] == filtered.size()) break;
    }
    std::reverse(chain.begin(), chain.end());
    return chain;
}
std::optional<Projection> nearest(const std::vector<Anchor>& anchors, std::uint32_t ref,
                                  std::uint32_t query_length, std::size_t radius) {
    auto it = std::lower_bound(anchors.begin(), anchors.end(), ref,
        [](const Anchor& a, std::uint32_t r){return a.reference < r;});
    const Anchor* found = nullptr;
    if (it != anchors.end()) found = &*it;
    if (it != anchors.begin()) {
        const Anchor* prev = &*std::prev(it);
        if (!found || std::llabs(static_cast<std::int64_t>(prev->reference) - ref) <=
                          std::llabs(static_cast<std::int64_t>(found->reference) - ref)) found = prev;
    }
    if (!found) return std::nullopt;
    const auto distance = std::llabs(static_cast<std::int64_t>(found->reference) - ref);
    if (distance > static_cast<std::int64_t>(radius)) return std::nullopt;
    const auto estimated = static_cast<std::int64_t>(found->query) +
                           (static_cast<std::int64_t>(ref) - found->reference);
    if (estimated < 0 || estimated > query_length) return std::nullopt;
    // If neighboring anchors straddle an insertion/deletion, the two nearest
    // projections can disagree. Reject that coordinate instead of silently
    // selecting one side and creating a misleading boundary.
    if (it != anchors.end() && it != anchors.begin()) {
        const auto& before = *std::prev(it);
        const auto& after = *it;
        if (ref >= before.reference && ref <= after.reference &&
            ref - before.reference <= radius && after.reference - ref <= radius) {
            const auto left = static_cast<std::int64_t>(before.query) + ref - before.reference;
            const auto right = static_cast<std::int64_t>(after.query) + ref - after.reference;
            if (std::llabs(left - right) > 64) return std::nullopt;
        }
    }
    return Projection{static_cast<std::uint32_t>(estimated), found->reference, found->query,
                      static_cast<std::uint32_t>(distance)};
}
Boundary makeBoundary(std::uint32_t ref, const Plan& p) {
    Boundary b;
    b.reference = ref;
    b.projections.reserve(p.supported_anchors.size());
    for (std::size_t j = 0; j < p.supported_anchors.size(); ++j) {
        const auto val = nearest(p.supported_anchors[j], ref,
                                 static_cast<std::uint32_t>(p.lengths[j + 1]),
                                 p.config.boundary_radius);
        if (val) ++b.supporting_queries;
        b.projections.push_back(val);
    }
    return b;
}
} // namespace

std::vector<Sequence> readFasta(const std::filesystem::path& input, Limits limits) {
    std::ifstream file(input, std::ios::binary);
    if (!file) throw std::runtime_error("cannot open input FASTA: " + input.string());
    // Scan a byte at a time: unlike getline, this cannot buffer an unbounded
    // malicious FASTA line before checking input budgets.
    std::vector<Sequence> reads;
    std::size_t total = 0, header_bytes = 0;
    bool line_start = true, in_header = false, id_complete = false;
    char ch = 0;
    while (file.get(ch)) {
        if (ch == '\r') continue;
        if (ch == '\n') { line_start = true; in_header = false; continue; }
        if (line_start) {
            line_start = false;
            if (ch == '>') {
                if (reads.size() >= limits.max_sequences)
                    throw std::invalid_argument("diagnostic planner limited to 8 sequences");
                reads.push_back({{}, {}});
                in_header = true;
                id_complete = false;
                header_bytes = 0;
                continue;
            }
            if (reads.empty()) throw std::invalid_argument("FASTA sequence before header");
        }
        if (in_header) {
            if (++header_bytes > 256) throw std::invalid_argument("FASTA header longer than 256 bytes");
            if (ch == ' ' || ch == '\t') id_complete = true;
            else if (!id_complete) reads.back().id.push_back(ch);
            continue;
        }
        if (reads.back().bases.size() >= limits.max_per_sequence ||
            total >= limits.max_total_bases)
            throw std::invalid_argument("FASTA input exceeds diagnostic size budget");
        const auto c = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T' && c != 'N')
            throw std::invalid_argument("FASTA contains unsupported base/gap");
        reads.back().bases.push_back(c);
        ++total;
    }
    if (file.bad()) throw std::runtime_error("I/O failure reading FASTA");
    validateInput(reads, limits);
    return reads;
}

Plan build(const std::vector<Sequence>& sequences, Config config, Limits limits) {
    validateInput(sequences, limits);
    if (config.k != 11 || config.target_block_bases < 5000 || config.target_block_bases > 100000 ||
        config.boundary_radius == 0 || config.boundary_radius > 1000 ||
        config.local_span == 0 || config.local_span > 20000 ||
        config.diagonal_tolerance == 0 || config.diagonal_tolerance > 64 ||
        config.min_independent_seeds != 4)
        throw std::invalid_argument("unsupported PR3D.1 experimental parameters");
    Plan plan;
    plan.config = config;
    for (const auto& s : sequences) {plan.ids.push_back(s.id); plan.lengths.push_back(s.bases.size());}
    const auto reference_index = makeIndex(sequences[0].bases, config.k);
    for (std::size_t i = 1; i < sequences.size(); ++i) {
        const auto query_index = makeIndex(sequences[i].bases, config.k);
        std::vector<Anchor> candidates;
        candidates.reserve(query_index.size() / 8 + 1);
        for (const auto& [key, q] : query_index) {
            const auto it = reference_index.find(key);
            if (it != reference_index.end() && it->second.count == 1 && q.count == 1)
                candidates.push_back({it->second.pos, q.pos});
        }
        plan.supported_anchors.push_back(evidenceFilter(std::move(candidates), config, sequences[0].bases, sequences[i].bases));
    }
    const auto ref_len = static_cast<std::uint32_t>(sequences[0].bases.size());
    plan.boundaries.push_back(makeBoundary(0, plan));
    const auto min_support = sequences.size() == 2 ? 1U : 2U;
    const auto min_span = std::max<std::size_t>(2000, config.target_block_bases / 4);
    const auto search_window = config.target_block_bases / 3;
    for (std::size_t target = config.target_block_bases;
         target + min_span < ref_len; target += config.target_block_bases) {
        const auto from = target > search_window ? target - search_window : 1;
        const auto to = std::min<std::size_t>(ref_len - min_span, target + search_window);
        std::vector<std::uint32_t> positions;
        for (const auto& chain : plan.supported_anchors) {
            auto it = std::lower_bound(chain.begin(), chain.end(), from,
                [](const Anchor& a, std::size_t x){return a.reference < x;});
            for (; it != chain.end() && it->reference <= to; ++it)
                positions.push_back(it->reference);
        }
        std::sort(positions.begin(), positions.end());
        positions.erase(std::unique(positions.begin(), positions.end()), positions.end());
        bool found = false;
        Boundary best;
        std::size_t closest = std::numeric_limits<std::size_t>::max();
        for (const auto r : positions) {
            if (r < plan.boundaries.back().reference + min_span || ref_len - r < min_span) continue;
            auto candidate = makeBoundary(r, plan);
            if (candidate.supporting_queries < min_support) continue;
            const auto distance = static_cast<std::size_t>(std::llabs(static_cast<std::int64_t>(r) - target));
            if (!found || candidate.supporting_queries > best.supporting_queries ||
                (candidate.supporting_queries == best.supporting_queries &&
                 (distance < closest || (distance == closest && r < best.reference)))) {
                found = true;
                best = std::move(candidate);
                closest = distance;
            }
        }
        if (found) plan.boundaries.push_back(std::move(best));
    }
    plan.boundaries.push_back(makeBoundary(ref_len, plan));
    return plan;
}

std::string asTsv(const Plan& plan) {
    if (plan.ids.size() < 2 || plan.boundaries.size() < 2)
        throw std::invalid_argument("invalid diagnostic plan");
    std::ostringstream out;
    out << "# PR3D.1 diagnostic only; NO alignment and NO verified homologous coordinates\n"
        << "# reference=" << plan.ids[0] << "\tlength=" << plan.lengths[0]
        << "\tk=" << plan.config.k << "\ttarget_bp=" << plan.config.target_block_bases << "\n"
        << "# coordinates: 0-based half-open; query boundaries are approximate projections from nearby unique-seed evidence\n"
        << "block\tread_id\tref_start\tref_end\tquery_start_est\tquery_end_est\tleft_anchor_distance\tright_anchor_distance\tsupported_seed_count\tstatus\n";
    for (std::size_t i = 1; i < plan.boundaries.size(); ++i) {
        const auto& left = plan.boundaries[i-1];
        const auto& right = plan.boundaries[i];
        for (std::size_t s = 0; s < plan.ids.size(); ++s) {
            out << i << '\t' << plan.ids[s] << '\t' << left.reference << '\t' << right.reference << '\t';
            if (s == 0) {
                out << left.reference << '\t' << right.reference << "\t0\t0\t0\tREFERENCE_ONLY\n";
                continue;
            }
            const auto& l = left.projections[s-1];
            const auto& r = right.projections[s-1];
            if (l) out << l->estimated_query; else out << "NA";
            out << '\t';
            if (r) out << r->estimated_query; else out << "NA";
            out << '\t';
            if (l) out << l->distance_to_anchor; else out << "NA";
            out << '\t';
            if (r) out << r->distance_to_anchor; else out << "NA";
            const auto& chain = plan.supported_anchors[s-1];
            const auto a = std::lower_bound(chain.begin(), chain.end(), left.reference,
                [](const Anchor& x, std::uint32_t p){return x.reference < p;});
            const auto b = std::lower_bound(chain.begin(), chain.end(), right.reference,
                [](const Anchor& x, std::uint32_t p){return x.reference < p;});
            const auto count = static_cast<std::size_t>(b - a);
            out << '\t' << count << '\t';
            if (l && r && l->estimated_query < r->estimated_query && count >= 3 &&
                right.reference - left.reference <= 2 * plan.config.target_block_bases)
                out << "CANDIDATE_ONLY\n";
            else out << "UNRESOLVED\n";
        }
    }
    return out.str();
}
void writeTsv(const Plan& plan, const std::filesystem::path& output) {
    if (output.empty()) throw std::invalid_argument("output TSV path is empty");
    if (std::filesystem::exists(output)) throw std::runtime_error("refusing to overwrite existing plan TSV");
    std::ofstream file(output, std::ios::binary | std::ios::out | std::ios::trunc);
    if (!file) throw std::runtime_error("cannot create plan TSV: " + output.string());
    file << asTsv(plan);
    file.flush();
    if (!file) throw std::runtime_error("failed writing plan TSV");
}
} // namespace homology_plan
