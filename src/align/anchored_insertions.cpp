#include "align/anchored_insertions.hpp"
#include "align.h"

#include <algorithm>
#include <limits>
#include <optional>
#include <unordered_set>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace align::anchored_insertions {
namespace {

struct Mapping {
    std::string anchors;
    std::map<std::size_t, std::string> inserts;
};

bool validDNA(const std::string& seq) {
    return !seq.empty() &&
           std::all_of(seq.begin(), seq.end(), [](char b) {
               return b == 'A' || b == 'C' || b == 'G' || b == 'T' || b == 'N';
           });
}

Mapping parseGlobalCigar(const std::string& reference,
                         const Query& query) {
    const auto operations = cigar::stringToCigar(query.cigar);
    if (operations.empty()) {
        throw std::runtime_error("Anchored insertions require a mapped CIGAR: " + query.id);
    }
    Mapping mapping;
    mapping.anchors.assign(reference.size(), '-');
    std::size_t r = 0, q = 0;
    for (const auto operation : operations) {
        char code{};
        std::uint32_t n{};
        cigar::intToCigar(operation, code, n);
        if (n == 0) throw std::runtime_error("Zero-length CIGAR operation");
        if (code == 'M' || code == '=' || code == 'X') {
            if (n > reference.size() - r || n > query.bases.size() - q)
                throw std::runtime_error("CIGAR match exceeds sequence bounds: " + query.id);
            mapping.anchors.replace(r, n, query.bases, q, n);
            r += n;
            q += n;
        } else if (code == 'D' || code == 'N') {
            if (n > reference.size() - r)
                throw std::runtime_error("CIGAR deletion exceeds reference: " + query.id);
            r += n;
        } else if (code == 'I') {
            if (n > query.bases.size() - q)
                throw std::runtime_error("CIGAR insertion exceeds query: " + query.id);
            mapping.inserts[r] += query.bases.substr(q, n);
            q += n;
        } else {
            throw std::runtime_error("Unsupported operation in global CIGAR: " + query.id);
        }
    }
    if (r != reference.size() || q != query.bases.size()) {
        throw std::runtime_error("Non-global or incomplete CIGAR: " + query.id);
    }
    return mapping;
}

std::string ungap(const std::string& aligned) {
    std::string result;
    result.reserve(aligned.size());
    for (char c : aligned) if (c != '-') result.push_back(c);
    return result;
}

// PR3B.5B: low-support reference blocks may contain competing lineage-specific
// insertions that were misclassified as reference-consuming CIGAR matches.
// This conservative heuristic uses only observed CIGARs and DNA, never an oracle.
struct SoftBlock {
    std::size_t start = 0, end = 0; // [start,end) in ungapped center
};
struct SoftHaplotype {
    std::size_t query_index;
    std::string bases;
    std::unordered_set<std::string> kmers;
};
struct SoftSplit {
    SoftBlock block;
    std::vector<std::vector<std::size_t>> groups; // indexes into haplotypes
    std::vector<SoftHaplotype> haplotypes;
    std::size_t center_group = 0;
    bool recruited = false; // split created by optional V3 fallback
};

std::unordered_set<std::string> kmerSet(const std::string& bases,
                                       std::size_t k = 7) {
    std::unordered_set<std::string> words;
    if (bases.empty()) return words;
    if (bases.size() < k) {
        words.insert(bases);
    } else {
        words.reserve(bases.size() - k + 1);
        for (std::size_t i = 0; i + k <= bases.size(); ++i)
            words.insert(bases.substr(i, k));
    }
    return words;
}

// Jaccard over SETS of 7-mers, exactly the benchmark's experimental criterion.
double kmerSimilarity(const std::unordered_set<std::string>& a,
                      const std::unordered_set<std::string>& b) {
    if (a.empty() && b.empty()) return 1.0;
    if (a.empty() || b.empty()) return 0.0;
    const auto& small = (a.size() <= b.size()) ? a : b;
    const auto& large = (a.size() <= b.size()) ? b : a;
    std::size_t overlap = 0;
    for (const auto& word : small) overlap += large.count(word);
    return static_cast<double>(overlap) / (a.size() + b.size() - overlap);
}

bool softCompatible(const SoftHaplotype& a, const SoftHaplotype& b) {
    if (a.bases.empty() || b.bases.empty()) return false;
    const auto shortest = std::min(a.bases.size(), b.bases.size());
    const auto longest = std::max(a.bases.size(), b.bases.size());
    if (static_cast<double>(shortest) / longest < 0.75) return false;
    return kmerSimilarity(a.kmers, b.kmers) >= 0.65;
}

std::string extractSoftHaplotype(const Mapping& mapping, SoftBlock block) {
    std::string result;
    for (std::size_t site = block.start; site < block.end; ++site) {
        const auto it = mapping.inserts.find(site);
        if (it != mapping.inserts.end()) result += it->second;
        if (mapping.anchors[site] != '-') result.push_back(mapping.anchors[site]);
    }
    const auto tail = mapping.inserts.find(block.end);
    if (tail != mapping.inserts.end()) result += tail->second;
    return result;
}

std::vector<SoftBlock> discoverSoftBlocks(const std::string& center,
                                          const std::vector<Mapping>& mappings) {
    // A candidate position has at most 70% exact center-base support and at
    // least 25% deletion/gap calls. Bridge up to 3 noncandidate columns.
    std::vector<bool> suspect(center.size(), false);
    for (std::size_t site = 0; site < center.size(); ++site) {
        std::size_t exact = 0, gap = 0;
        for (const auto& m : mappings) {
            exact += (m.anchors[site] == center[site]);
            gap += (m.anchors[site] == '-');
        }
        const double depth = static_cast<double>(mappings.size());
        suspect[site] = exact/depth <= 0.70 && gap/depth >= 0.25;
    }
    std::vector<SoftBlock> blocks;
    for (std::size_t site = 0; site < suspect.size();) {
        if (!suspect[site]) { ++site; continue; }
        const auto start = site;
        auto last = site;
        auto cursor = site + 1;
        while (cursor < suspect.size()) {
            if (suspect[cursor]) last = cursor;
            else if (cursor - last > 3) break;
            ++cursor;
        }
        if (last + 1 - start >= 20) blocks.push_back({start, last + 1});
        site = last + 1;
    }
    return blocks;
}

std::optional<SoftSplit> evaluateSoftBlock(const std::string& center,
                                           const std::vector<Mapping>& mappings,
                                           SoftBlock block) {
    SoftSplit split;
    split.block = block;
    for (std::size_t i = 0; i < mappings.size(); ++i) {
        std::string bases = extractSoftHaplotype(mappings[i], block);
        if (bases.empty()) continue;
        auto words = kmerSet(bases);
        split.haplotypes.push_back({i, std::move(bases), std::move(words)});
    }
    if (split.haplotypes.empty()) return std::nullopt;
    // Greedy complete linkage, preserving original FASTA order. The cluster
    // is selected by the first eligible group, matching the verified pilot.
    for (std::size_t i = 0; i < split.haplotypes.size(); ++i) {
        bool placed = false;
        for (auto& group : split.groups) {
            if (std::all_of(group.begin(), group.end(), [&](std::size_t member) {
                return softCompatible(split.haplotypes[i], split.haplotypes[member]);
            })) {
                group.push_back(i);
                placed = true;
                break;
            }
        }
        if (!placed) split.groups.push_back({i});
    }
    std::vector<std::size_t> core;
    std::size_t residual_rows = 0, residual_groups = 0;
    for (std::size_t gi = 0; gi < split.groups.size(); ++gi) {
        const auto& group = split.groups[gi];
        const bool long_enough = std::all_of(group.begin(), group.end(),
            [&](std::size_t h) { return split.haplotypes[h].bases.size() >= 20; });
        if (group.size() >= 2 && long_enough) core.push_back(gi);
        else { residual_rows += group.size(); ++residual_groups; }
    }
    if (core.size() < 2 || core.size() > 3 || residual_groups > 3 ||
        split.groups.size() > 6 ||
        residual_rows > std::max<std::size_t>(1, split.haplotypes.size()/5))
        return std::nullopt;

    const std::string fragment = center.substr(block.start, block.end - block.start);
    const auto fragment_words = kmerSet(fragment);
    std::vector<std::pair<double,std::size_t>> candidates;
    for (auto gi : core) {
        double lowest_score = 1.0;
        double lowest_ratio = 1.0;
        for (auto hi : split.groups[gi]) {
            const auto& h = split.haplotypes[hi];
            lowest_score = std::min(lowest_score, kmerSimilarity(fragment_words, h.kmers));
            lowest_ratio = std::min(lowest_ratio, static_cast<double>(
                std::min(fragment.size(), h.bases.size())) /
                std::max(fragment.size(), h.bases.size()));
        }
        if (lowest_score >= 0.65 && lowest_ratio >= 0.75)
            candidates.emplace_back(lowest_score, gi);
    }
    if (candidates.empty()) return std::nullopt;
    std::sort(candidates.begin(), candidates.end(),
        [](const auto& a, const auto& b) {
            if (a.first != b.first) return a.first > b.first;
            return a.second < b.second;
        });
    if (candidates.size() > 1 && candidates[0].first - candidates[1].first < 0.15)
        return std::nullopt;
    split.center_group = candidates[0].second;
    // Reject low-complexity cores: sparse 7-mer spectra can make two
    // homologous repetitive alleles look artificially unrelated.
    for (auto group_index : core) {
        for (auto hap_index : split.groups[group_index]) {
            if (split.haplotypes[hap_index].kmers.size() < 10)
                return std::nullopt;
        }
    }
    double max_cross = 0.0;
    for (std::size_t a = 0; a < core.size(); ++a) {
        for (std::size_t b = a + 1; b < core.size(); ++b) {
            for (const auto hi : split.groups[core[a]]) {
                for (const auto hj : split.groups[core[b]]) {
                    max_cross = std::max(max_cross, kmerSimilarity(
                        split.haplotypes[hi].kmers, split.haplotypes[hj].kmers));
                }
            }
        }
    }
    if (max_cross > 0.20) return std::nullopt;
    return split;
}


// PR3B.5C: when the strict V2 complete-linkage clusters fragment mutated
// copies into singletons, try a second, independently audited family model.
// No evolutionary labels or oracle are used at alignment time. Deliberately
// preserve the original FASTA order throughout union-find and local abPOA.
std::optional<SoftSplit> evaluateRecruitBlock(const std::string& center,
                                               const std::vector<Mapping>& mappings,
                                               SoftBlock block) {
    constexpr std::size_t k = 5;
    constexpr double edge = 0.25, min_internal = 0.20;
    constexpr double max_between5 = 0.15, max_between7 = 0.10;
    constexpr double min_center = 0.25, center_margin = 0.12;
    constexpr double min_ratio = 0.65;
    constexpr std::size_t min_core = 2, min_hap_len = 20;
    constexpr std::size_t min_unique_kmers = 12;

    SoftSplit split;
    split.block = block;
    split.recruited = true;
    for (std::size_t i = 0; i < mappings.size(); ++i) {
        auto bases = extractSoftHaplotype(mappings[i], block);
        if (!bases.empty())
            split.haplotypes.push_back({i, std::move(bases), {}});
    }
    const auto n = split.haplotypes.size();
    if (n < 2 * min_core) return std::nullopt;
    std::vector<std::unordered_set<std::string>> words5;
    std::vector<std::unordered_set<std::string>> words7;
    words5.reserve(n); words7.reserve(n);
    for (auto& h : split.haplotypes) {
        words5.push_back(kmerSet(h.bases, k));
        words7.push_back(kmerSet(h.bases, 7));
    }
    // Match the V3 pilot's stable connected components, but apply complete
    // pairwise audits before accepting any cluster as genuinely cohesive.
    std::vector<std::size_t> parent(n);
    for (std::size_t i = 0; i < n; ++i) parent[i] = i;
    auto rootOf = [&](std::size_t i) {
        while (parent[i] != i) {
            parent[i] = parent[parent[i]];
            i = parent[i];
        }
        return i;
    };
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = i + 1; j < n; ++j) {
            const auto& a = split.haplotypes[i].bases;
            const auto& b = split.haplotypes[j].bases;
            if (static_cast<double>(std::min(a.size(),b.size())) /
                    std::max(a.size(),b.size()) < min_ratio ||
                kmerSimilarity(words5[i],words5[j]) < edge) continue;
            auto ri = rootOf(i), rj = rootOf(j);
            if (ri != rj) parent[rj] = ri;
        }
    }
    // No unordered container iteration: the earliest observed haplotype
    // determines component order, exactly as in the Python validation.
    std::map<std::size_t, std::size_t> component_to_group;
    for (std::size_t i = 0; i < n; ++i) {
        auto root = rootOf(i);
        auto it = component_to_group.find(root);
        if (it == component_to_group.end()) {
            component_to_group[root] = split.groups.size();
            split.groups.emplace_back();
            it = component_to_group.find(root);
        }
        split.groups[it->second].push_back(i);
    }
    std::vector<std::size_t> core;
    std::size_t residual_rows = 0, residual_groups = 0;
    for (std::size_t gi = 0; gi < split.groups.size(); ++gi) {
        const auto& group = split.groups[gi];
        bool supported = group.size() >= min_core;
        for (auto h : group)
            supported = supported && split.haplotypes[h].bases.size() >= min_hap_len
                      && words5[h].size() >= min_unique_kmers;
        if (supported) core.push_back(gi);
        else { residual_rows += group.size(); ++residual_groups; }
    }
    if (core.size() < 2 || core.size() > 3 || residual_groups > 3 ||
        split.groups.size() > 6 || residual_rows > std::max<std::size_t>(1,n/5))
        return std::nullopt;

    for (auto gi : core) {
        const auto& group = split.groups[gi];
        for (std::size_t a = 0; a < group.size(); ++a)
            for (std::size_t b = a + 1; b < group.size(); ++b)
                if (kmerSimilarity(words5[group[a]],words5[group[b]]) < min_internal)
                    return std::nullopt;
    }
    for (std::size_t a = 0; a < core.size(); ++a) {
        for (std::size_t b = a + 1; b < core.size(); ++b) {
            for (auto ia : split.groups[core[a]]) {
                for (auto ib : split.groups[core[b]]) {
                    if (kmerSimilarity(words5[ia],words5[ib]) > max_between5 ||
                        kmerSimilarity(words7[ia],words7[ib]) > max_between7)
                        return std::nullopt;
                }
            }
        }
    }
    const auto fragment = center.substr(block.start,block.end-block.start);
    const auto fragment5 = kmerSet(fragment,k);
    std::vector<std::pair<double,std::size_t>> candidates;
    for (auto gi : core) {
        double lowest_score = 1.0, lowest_ratio = 1.0;
        for (auto h : split.groups[gi]) {
            lowest_score = std::min(lowest_score,kmerSimilarity(fragment5,words5[h]));
            lowest_ratio = std::min(lowest_ratio,static_cast<double>(
                std::min(fragment.size(),split.haplotypes[h].bases.size())) /
                std::max(fragment.size(),split.haplotypes[h].bases.size()));
        }
        if (lowest_score >= min_center && lowest_ratio >= 0.75)
            candidates.emplace_back(lowest_score,gi);
    }
    if (candidates.empty()) return std::nullopt;
    std::sort(candidates.begin(),candidates.end(),[](const auto& a,const auto& b) {
        if (a.first != b.first) return a.first > b.first;
        return a.second < b.second;
    });
    if (candidates.size() > 1 &&
        candidates[0].first-candidates[1].first < center_margin)
        return std::nullopt;
    split.center_group = candidates[0].second;
    return split;
}

std::vector<SoftSplit> selectSoftSplits(const std::string& center,
                                        const std::vector<Mapping>& mappings,
                                        std::size_t& candidate_count,
                                        std::size_t max_soft_blocks,
                                        bool recruit) {
    auto candidates = discoverSoftBlocks(center, mappings);
    candidate_count = candidates.size();
    // Reject the whole experimental soft path when there are too many
    // candidates. Keep the old anchored algorithm as a safe fallback.
    if (candidates.size() > max_soft_blocks) return {};
    std::vector<SoftSplit> chosen;
    for (const auto block : candidates) {
        if (auto split = evaluateSoftBlock(center, mappings, block))
            chosen.push_back(std::move(*split));
        else if (recruit) {
            if (auto alternate = evaluateRecruitBlock(center, mappings, block))
                chosen.push_back(std::move(*alternate));
        }
    }
    return chosen;
}

// Assemble only accepted soft blocks in separate haplotype lanes. For each
// other CIGAR insertion site, retain the PR3B.5A local MSA behavior unchanged.
Alignment assembleSoft(const std::string& center,
                       const std::vector<Query>& queries,
                       const std::vector<Mapping>& mappings,
                       const std::map<std::size_t, std::vector<std::size_t>>& sites,
                       const std::vector<SoftSplit>& splits,
                       const Options& limits,
                       std::size_t candidate_count) {
    Alignment result;
    result.soft_anchor_candidates = candidate_count;
    result.soft_anchor_splits = splits.size();
    result.soft_anchor_recruits = std::count_if(
        splits.begin(),splits.end(),[](const SoftSplit& s){return s.recruited;});
    result.occupied_sites = sites.size();
    const std::size_t row_count = queries.size()+1;
    result.rows.reserve(row_count);
    result.rows.push_back({"consensus", ""});
    for (const auto& query : queries) result.rows.push_back({query.id,""});
    // A maximum width of (center bases + all original query bases) is a
    // conservative upper bound, even when split blocks add separate lanes.
    const std::size_t width_limit = std::min(limits.max_output_cells/row_count,
                                             limits.max_total_bases);
    std::map<std::size_t, std::size_t> split_at;
    std::vector<bool> used_inserts(center.size()+1, false);
    for (std::size_t i = 0; i < splits.size(); ++i) {
        split_at.emplace(splits[i].block.start, i);
        for (std::size_t site = splits[i].block.start;
             site <= splits[i].block.end; ++site) {
            if (used_inserts[site])
                throw std::runtime_error("Overlapping soft-anchor insertion site");
            used_inserts[site] = true;
        }
    }
    auto append = [&](const std::vector<std::pair<std::size_t,std::string>>& items) {
        if (items.empty()) return;
        std::vector<native_msa::Sequence> input;
        input.reserve(items.size());
        bool identical = true;
        for (const auto& [row,bases] : items) {
            if (bases.empty() || row >= row_count)
                throw std::runtime_error("Invalid empty soft-anchor haplotype");
            if (!input.empty() && bases != input[0].bases) identical = false;
            input.push_back({result.rows[row].id,bases});
        }
        native_msa::Alignment aligned;
        if (input.size() == 1 || identical) {
            aligned.rows = std::move(input);
            aligned.columns = aligned.rows[0].bases.size();
        } else {
            native_msa::Options local_limits;
            aligned = native_msa::align(input, local_limits);
            ++result.local_msas;
        }
        if (aligned.columns == 0 || aligned.rows.size() != items.size())
            throw std::runtime_error("Soft-anchor local alignment has invalid shape");
        if (result.rows[0].bases.size() > width_limit ||
            aligned.columns > width_limit - result.rows[0].bases.size())
            throw std::invalid_argument("Soft-anchor output exceeds bounded cell budget");
        std::vector<std::string> local(row_count, std::string(aligned.columns, '-'));
        for (std::size_t j = 0; j < items.size(); ++j) {
            const auto row = items[j].first;
            if (aligned.rows[j].id != result.rows[row].id ||
                aligned.rows[j].bases.size() != aligned.columns ||
                ungap(aligned.rows[j].bases) != items[j].second)
                throw std::runtime_error("Soft-anchor local MSA changed query order or bases");
            local[row] = std::move(aligned.rows[j].bases);
        }
        for (std::size_t r = 0; r < row_count; ++r)
            result.rows[r].bases += local[r];
    };
    for (std::size_t site = 0; site <= center.size();) {
        const auto split_it = split_at.find(site);
        if (split_it != split_at.end()) {
            const auto& split = splits[split_it->second];
            for (std::size_t group_index = 0; group_index < split.groups.size(); ++group_index) {
                std::vector<std::pair<std::size_t,std::string>> items;
                if (group_index == split.center_group)
                    items.emplace_back(0, center.substr(split.block.start,
                                                         split.block.end - split.block.start));
                for (const auto member : split.groups[group_index]) {
                    const auto& h = split.haplotypes[member];
                    items.emplace_back(h.query_index + 1,h.bases);
                }
                append(items);
            }
            site = split.block.end;
            continue;
        }
        if (!used_inserts[site]) {
            const auto it = sites.find(site);
            if (it != sites.end()) {
                std::vector<std::pair<std::size_t,std::string>> items;
                items.reserve(it->second.size());
                for (const auto qi : it->second)
                    items.emplace_back(qi + 1,mappings[qi].inserts.at(site));
                append(items);
            }
        }
        if (site < center.size()) {
            if (result.rows[0].bases.size() >= width_limit)
                throw std::invalid_argument("Soft-anchor output exceeds bounded cell budget");
            result.rows[0].bases.push_back(center[site]);
            for (std::size_t qi = 0; qi < queries.size(); ++qi)
                result.rows[qi+1].bases.push_back(mappings[qi].anchors[site]);
        }
        ++site;
    }
    result.columns = result.rows[0].bases.size();
    for (std::size_t qi = 0; qi < queries.size(); ++qi) {
        if (result.rows[qi+1].bases.size() != result.columns ||
            ungap(result.rows[qi+1].bases) != queries[qi].bases)
            throw std::runtime_error("Soft-anchor MSA altered query bases");
    }
    if (ungap(result.rows[0].bases) != center)
        throw std::runtime_error("Soft-anchor MSA altered center bases");
    return result;
}

} // namespace

Alignment align(const std::string& consensus,
                const std::vector<Query>& queries,
                const Options& limits) {
    if (!validDNA(consensus) || consensus.size() > limits.max_center_bases)
        throw std::invalid_argument("Anchored center is invalid or exceeds 25kb budget");
    if (queries.empty() || queries.size() > limits.max_sequences)
        throw std::invalid_argument("Anchored query count exceeds the bounded budget");
    if (limits.max_total_bases < consensus.size())
        throw std::invalid_argument("Anchored total-base budget is smaller than the center");
    if (queries.size() == std::numeric_limits<std::size_t>::max())
        throw std::invalid_argument("Anchored row count overflow");

    std::unordered_set<std::string> names;
    std::vector<Mapping> mappings;
    mappings.reserve(queries.size());
    std::map<std::size_t, std::vector<std::size_t>> sites;
    std::size_t total_bases = consensus.size();
    std::size_t total_insertion_bases = 0;
    for (std::size_t i = 0; i < queries.size(); ++i) {
        const Query& query = queries[i];
        if (query.ordinal != i || query.id.empty() || query.id == "consensus" ||
            !names.insert(query.id).second)
            throw std::invalid_argument("Anchored queries need unique IDs and consecutive ordinals");
        if (!validDNA(query.bases) || query.bases.size() > limits.max_center_bases ||
            query.bases.size() > limits.max_total_bases - total_bases)
            throw std::invalid_argument("Anchored query DNA or sequence-size budget invalid");
        total_bases += query.bases.size();
        mappings.push_back(parseGlobalCigar(consensus, query));
        for (const auto& [site, bases] : mappings.back().inserts) {
            sites[site].push_back(i);
            if (bases.size() > limits.max_total_bases - total_insertion_bases)
                throw std::invalid_argument("Anchored insertion budget exceeded");
            total_insertion_bases += bases.size();
        }
    }

    // The soft mode is strictly opt-in; if no block has two distinct,
    // well-supported cores, preserve the PR3B.5A assembly byte-for-byte.
    std::size_t soft_candidates = 0;
    if (limits.soft_anchors) {
        const auto splits = selectSoftSplits(consensus, mappings, soft_candidates,
                                             limits.max_soft_blocks,
                                             limits.soft_anchor_recruit);
        if (!splits.empty())
            return assembleSoft(consensus, queries, mappings, sites,
                                splits, limits, soft_candidates);
    }

    // A local MSA width cannot exceed the sum of its input sequence lengths.
    // This is a conservative bound on the entire output profile allocation.
    const std::size_t max_width = consensus.size() + total_insertion_bases;
    const std::size_t row_count = queries.size() + 1;
    if (max_width > limits.max_output_cells / row_count)
        throw std::invalid_argument("Anchored final MSA exceeds cell budget");

    Alignment result;
    result.soft_anchor_candidates = soft_candidates;
    result.occupied_sites = sites.size();
    result.rows.reserve(row_count);
    result.rows.push_back({"consensus", ""});
    for (const auto& query : queries) result.rows.push_back({query.id, ""});
    for (auto& row : result.rows) row.bases.reserve(max_width);

    native_msa::Options local_limits;
    // Use the same bounded DNA backend as PR2B, never bypass its budgets.
    for (std::size_t site = 0; site <= consensus.size(); ++site) {
        auto it = sites.find(site);
        if (it != sites.end()) {
            std::vector<native_msa::Sequence> input;
            input.reserve(it->second.size());
            bool identical = true;
            std::string first;
            for (std::size_t qi : it->second) {
                const std::string& bases = mappings[qi].inserts.at(site);
                if (first.empty()) first = bases;
                if (bases != first) identical = false;
                input.push_back({queries[qi].id, bases});
            }
            native_msa::Alignment aligned;
            if (input.size() == 1 || identical) {
                aligned.rows = std::move(input);
                aligned.columns = first.size();
            } else {
                if (input.size() > local_limits.max_sequences)
                    throw std::invalid_argument("Anchored local insertion group exceeds NativeMSA budget");
                // Local sequences are ungapped; NativeMSA checks its own
                // 25kb/sequence, 100kb/total, 64-record limits before abPOA.
                aligned = native_msa::align(input, local_limits);
                ++result.local_msas;
            }
            if (aligned.rows.size() != it->second.size() || aligned.columns == 0)
                throw std::runtime_error("Anchored local MSA returned invalid dimensions");
            // NativeMSA keeps input order; verify IDs and bases for all paths.
            std::vector<std::string> inserted(queries.size());
            for (std::size_t j = 0; j < it->second.size(); ++j) {
                const std::size_t qi = it->second[j];
                if (aligned.rows[j].id != queries[qi].id ||
                    aligned.rows[j].bases.size() != aligned.columns ||
                    ungap(aligned.rows[j].bases) != mappings[qi].inserts.at(site))
                    throw std::runtime_error("Anchored local MSA corrupted input order or bases");
                inserted[qi] = std::move(aligned.rows[j].bases);
            }
            result.rows[0].bases.append(aligned.columns, '-');
            for (std::size_t qi = 0; qi < queries.size(); ++qi) {
                if (inserted[qi].empty())
                    result.rows[qi + 1].bases.append(aligned.columns, '-');
                else
                    result.rows[qi + 1].bases += inserted[qi];
            }
        }
        if (site < consensus.size()) {
            result.rows[0].bases.push_back(consensus[site]);
            for (std::size_t qi = 0; qi < queries.size(); ++qi)
                result.rows[qi + 1].bases.push_back(mappings[qi].anchors[site]);
        }
    }

    result.columns = result.rows[0].bases.size();
    for (std::size_t qi = 0; qi < queries.size(); ++qi) {
        if (result.rows[qi + 1].bases.size() != result.columns ||
            ungap(result.rows[qi + 1].bases) != queries[qi].bases)
            throw std::runtime_error("Anchored final MSA altered or lost query bases");
    }
    if (ungap(result.rows[0].bases) != consensus)
        throw std::runtime_error("Anchored final MSA altered consensus bases");
    return result;
}

} // namespace align::anchored_insertions
