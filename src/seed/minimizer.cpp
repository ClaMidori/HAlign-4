#include "seed.h"

#include <cstdint>
#include <string>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <cmath>
#include <limits>
#include "anchor.h"

namespace minimizer
{
    // splitmix64: widely-used 64-bit mixer with good speed and distribution.
    // minimap2 also applies hash64 mixing to k-mer encoding with same concept.
    static inline constexpr std::uint64_t splitmix64(std::uint64_t x) noexcept
    {
        x += 0x9e3779b97f4a7c15ULL;
        x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
        x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
        return x ^ (x >> 31);
    }

    // Cand: candidate (hash, pos) in window
    // pos is k-mer start position (0-based)
    struct Cand
    {
        std::uint64_t h;
        std::uint32_t pos;
    };

    // =========================================================
    // RingMinQueue: fixed-capacity ring monotonic queue (maintains window minimum)
    // =========================================================
    class RingMinQueue
    {
    public:
        explicit RingMinQueue(std::uint32_t capacity)
            : buf_(capacity), cap_(capacity)
        {
        }

        void clear() noexcept
        {
            head_ = 0;
            size_ = 0;
        }

        bool empty() const noexcept { return size_ == 0; }

        void push(std::uint64_t h, std::uint32_t pos) noexcept
        {
            const Cand c{h, pos};
            // Pop all elements >= current hash from back to maintain monotonic order
            while (size_ && back().h >= c.h) pop_back();
            buf_[idx(size_)] = c;
            ++size_;
        }

        void popExpired(std::uint32_t win_start) noexcept
        {
            while (size_ && front().pos < win_start) pop_front();
        }

        std::uint64_t minHash() const noexcept { return front().h; }
        std::uint32_t minPos() const noexcept { return front().pos; }

    private:
        std::vector<Cand> buf_;
        std::uint32_t cap_{0};
        std::uint32_t head_{0};
        std::uint32_t size_{0};

        std::uint32_t idx(std::uint32_t off) const noexcept
        {
            std::uint32_t i = head_ + off;
            if (i >= cap_) i -= cap_;
            return i;
        }

        Cand& front() noexcept { return buf_[head_]; }
        const Cand& front() const noexcept { return buf_[head_]; }

        Cand& back() noexcept { return buf_[idx(size_ - 1)]; }
        const Cand& back() const noexcept { return buf_[idx(size_ - 1)]; }

        void pop_front() noexcept
        {
            head_ = idx(1);
            --size_;
        }

        void pop_back() noexcept
        {
            --size_;
        }
    };


    // =============================================================
    // Goal: extract minimizers (with position hits) from single sequence
    // =============================================================
    MinimizerHits extractMinimizer(const std::string& seq,
                                   std::size_t k,
                                   std::size_t w,
                                   bool non_canonical)
    {
        MinimizerHits out;

        const std::uint32_t n = static_cast<std::uint32_t>(seq.size());
        if (k == 0 || w == 0 || n < k) return out;
        if (k > 31) return out;
        if (w >= 256) return out;

        const auto& nt4_table_ref = minimizer::nt4_table;

        const std::uint32_t total_kmer = n - static_cast<std::uint32_t>(k) + 1;
        const std::uint32_t win = std::min<std::uint32_t>(static_cast<std::uint32_t>(w), total_kmer);
        if (win == 0) return out;

        out.reserve(std::max<std::uint32_t>(1, n / win));

        const std::uint64_t mask = (1ULL << (2 * k)) - 1ULL;
        const std::uint64_t shift = 2ULL * (k - 1);

        std::uint64_t fwd = 0;
        std::uint64_t rev = 0;
        std::uint32_t valid = 0;

        RingMinQueue q(win);

        Cand last_out{0, 0};
        bool has_last = false;

        for (std::uint32_t i = 0; i < n; ++i) {
            const std::uint8_t c = nt4_table_ref[static_cast<std::uint8_t>(seq[i])];
            if (c >= 4) {
                fwd = rev = 0;
                valid = 0;
                q.clear();
                has_last = false;
                continue;
            }

            fwd = ((fwd << 2) | c) & mask;
            rev = (rev >> 2) | (std::uint64_t(3U ^ c) << shift);

            if (valid < k) ++valid;
            if (valid < k) continue;

            const std::uint32_t pos = i + 1 - static_cast<std::uint32_t>(k);

            const std::uint64_t code = non_canonical ? fwd : std::min(fwd, rev);
            const std::uint64_t h64 = splitmix64(code);
            const std::uint64_t h56 = (h64 >> 8); // Extract high 56 bits, reserve low 8 bits for span

            q.push(h56, pos);

            if (pos + 1 < win) continue;
            const std::uint32_t win_start = pos + 1 - win;
            q.popExpired(win_start);

            if (!q.empty()) {
                const Cand cur{q.minHash(), q.minPos()};
                if (!has_last || cur.h != last_out.h || cur.pos != last_out.pos) {
                    out.emplace_back(cur.h, cur.pos,
                                     /*rid*/ 0,
                                     non_canonical ? true : (fwd <= rev),
                                     static_cast<std::uint8_t>(k));
                    last_out = cur;
                    has_last = true;
                }
            }
        }

        return out;
    }


// ------------------------------------------------------------------
// Function: collect_anchors
// ------------------------------------------------------------------
// Purpose:
// Collect anchor list (Anchor) from ref_hits and qry_hits.
//
// Anchor semantics here:
// - One anchor represents: ref and query match on a minimizer hash;
// - Since same hash may appear multiple times in ref/query, one qry_hit may expand to multiple anchors.
// - Anchors are input for subsequent chaining.
//
// Key design points aligned with minimap2 (very important):
// 1) Filter first, then expand (occurrence expansion)
//    - minimap2's -f/-U/--q-occ-frac/-e strategies suppress repeat regions "before expanding occurrences".
//    - Without this, repeat regions would expand high-frequency minimizer to O(occ_ref * occ_qry) anchors,
//      directly causing memory/time explosion.
//
// 2) This implementation's filter parameters come from anchor::SeedFilterParams (defaults mimic minimap2 CLI):
//    - f_top_frac (-f): ignore most frequent top fraction minimizers in reference (by distinct count)
//    - u_floor/u_ceil (-U): bound occurrence threshold with upper/lower limits
//    - q_occ_frac (--q-occ-frac): discard if query-side too frequent (and exceeds reference threshold)
//    - sample_every_bp (-e): sparse position sampling for high-frequency minimizers (not full expansion)
//
// 3) Output anchors are unsorted:
//    - minimap2 later sorts by (rid, strand, diagonal, ref_pos, qry_pos) before DP chaining.
//    - This only collects; sorting is done by caller (anchor::sortAnchorsByDiagonal etc).
//
// Input:
// - ref_hits: reference minimizer hits (may be from full reference or portion)
// - qry_hits: query minimizer hits
//
// Output:
// - anchor::Anchors: anchor list (each records ref/qry rid/pos/span/is_rev)
//
// Complexity:
// - Sort ref_hits: O(R log R)
// - Count qry occurrences: O(Q)
// - Generate anchors: O(Q * avg_occ_ref_for_hash)
//   (repeat regions filtered/sparse-sampled, avoids worst-case explosion)
// ------------------------------------------------------------------
anchor::Anchors collect_anchors(const MinimizerHits& ref_hits, const MinimizerHits& qry_hits, anchor::SeedFilterParams params)
{
    anchor::Anchors anchors;

    // Boundary condition: if either side is empty, cannot produce anchors
    if (ref_hits.empty() || qry_hits.empty()) {
        return anchors;
    }


    // ------------------------------------------------------------------
    // Step 1: sort ref_hits + build hash -> (start, count) index
    // ------------------------------------------------------------------
    // Note:
    // - Same hash may appear multiple times in ref_hits (different positions/rids).
    // - After sorting by (hash, rid, pos, strand), hits of same hash become contiguous.
    // - We use unordered_map to store contiguous interval (start, count) for each hash;
    //   when query finds a hash, can expand in O(occ_ref) time.
    std::vector<minimizer::MinimizerHit> sorted_ref = ref_hits;
    std::sort(sorted_ref.begin(), sorted_ref.end()); // Use SeedHitBase's operator<

    std::unordered_map<hash_t, anchor::HashIndex> hash_index;
    hash_index.reserve(sorted_ref.size());

    // ref_occs: record occurrence count of each distinct hash in reference,
    // used later to estimate filter threshold based on -f (top fraction)
    std::vector<std::size_t> ref_occs;
    ref_occs.reserve(sorted_ref.size() / 2 + 1);

    if (!sorted_ref.empty()) {
        std::size_t start = 0;
        hash_t current_hash = sorted_ref[0].hash();

        for (std::size_t i = 1; i <= sorted_ref.size(); ++i) {
            // When hash changes or reach end => finalize interval for previous hash
            if (i == sorted_ref.size() || sorted_ref[i].hash() != current_hash) {
                const std::size_t occ = i - start;
                hash_index[current_hash] = anchor::HashIndex{start, occ};
                ref_occs.push_back(occ);
                if (i < sorted_ref.size()) {
                    start = i;
                    current_hash = sorted_ref[i].hash();
                }
            }
        }
    }

    // ------------------------------------------------------------------
    // Step 2: compute "reference-side high-frequency threshold" from ref_occs + (-f/-U)
    // ------------------------------------------------------------------
    // This threshold determines if a hash belongs to repeat region (high-frequency).
    // Note: this is key point of "filtering happens before occurrence expansion".
    const std::size_t ref_occ_thr = anchor::compute_ref_occ_threshold(ref_occs, params);

    // ------------------------------------------------------------------
    // Step 3: count occurrence of each hash on query side (for --q-occ-frac)
    // ------------------------------------------------------------------
    // minimap2's intuition:
    // - If a hash is extremely high-frequency in query, likely from low-complexity/repeat region
    // - These seeds don't help much with localization but produce many anchors
    std::unordered_map<hash_t, std::size_t> qry_occ;
    qry_occ.reserve(qry_hits.size());
    for (const auto& qh : qry_hits) {
        ++qry_occ[qh.hash()];
    }

    // q_occ_limit: convert --q-occ-frac from "ratio" to "count threshold"
    const double q_occ_limit = params.q_occ_frac > 0.0
        ? (params.q_occ_frac * static_cast<double>(qry_hits.size()))
        : std::numeric_limits<double>::infinity();

    // Estimate: anchor count usually same order as qry_hits (repeats suppressed)
    anchors.reserve(qry_hits.size());

    // ------------------------------------------------------------------
    // Step 4: iterate qry_hits, query ref index and generate anchors
    // ------------------------------------------------------------------
    for (const auto& qry_hit : qry_hits) {
        const hash_t qry_hash = qry_hit.hash();

        // Hash doesn't exist on ref side => cannot form anchor
        auto it = hash_index.find(qry_hash);
        if (it == hash_index.end()) continue;

        const anchor::HashIndex& idx = it->second;
        const std::size_t ref_occ = idx.count;

        // ---- 4.1 --q-occ-frac: discard if query-side too high-frequency (reduce explosion risk)
        // minimap2 semantics: when query hash occurrence exceeds threshold, discard that hash
        // This independent from ref-side filtering, handles low-complexity regions in query
        if (params.q_occ_frac > 0.0) {
            const std::size_t qocc = qry_occ[qry_hash];
            if (static_cast<double>(qocc) > q_occ_limit) {
                continue;
            }
        }

        // ---- 4.2 -f/-U + -e: reference-side high-frequency minimizers use sparse sampling
        // Semantics:
        // - If ref_occ > ref_occ_thr, this hash is very "repeated" in reference.
        // - minimap2 sparsely samples high-frequency minimizers instead of expanding all occurrences.
        // - We use simplest "modulo sampling by query position" to approximate -e behavior:
        //   only expand when qry_hit.pos % sample_every_bp == 0.
        //
        // Note: not unique implementation, but maintains key property of "hotpath not exploding".
        if (ref_occ > ref_occ_thr) {
            if (params.sample_every_bp == 0) continue;
            if ((static_cast<std::size_t>(qry_hit.pos()) % params.sample_every_bp) != 0) {
                continue;
            }
        }

        // ---- 4.3 expand reference occurrences, generate anchors
        // Generation only copies fields/light computation, avoids extra allocation
        for (std::size_t i = idx.start; i < idx.start + idx.count; ++i) {
            const auto& ref_hit = sorted_ref[i];

            anchor::Anchor anchor;
            anchor.hash = qry_hash;
            anchor.rid_ref = ref_hit.rid();
            anchor.pos_ref = ref_hit.pos();
            anchor.rid_qry = qry_hit.rid();
            anchor.pos_qry = qry_hit.pos();

            // span takes smaller of two: conservative estimate of "usable match length"
            // (downstream chaining's gap/penalty model typically needs just one span magnitude)
            anchor.span = std::min(ref_hit.span(), qry_hit.span());

            // Direction: ref XOR qry (consistent with minimap2's rev semantics)
            anchor.is_rev = (ref_hit.strand() != qry_hit.strand());

            anchors.emplace_back(anchor);
        }
    }

    return anchors;
}

} // namespace minimizer

