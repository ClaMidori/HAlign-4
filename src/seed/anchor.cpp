// ==================================================================
// chain.cpp - minimap2-style anchor collection (seeding) and preprocessing
// ==================================================================
//
// Functionality description:
// This file implements anchor collection for sequence alignment, following minimap2's design.
// Anchors are seed match positions shared between ref and query sequences, used
// for subsequent chaining and precise alignment.
//
// References:
// - minimap2/seed.c: mm_collect_matches() - collect seed matches
// - minimap2/hit.c: anchor processing and filtering
// - minimap2/lchain.c: chaining algorithm (DP dynamic programming)
//
// Core data flow (minimap2-style, see README Algorithm overview):
// 1) Build an index for the reference sequence minimizers (hash -> occurrences)
// 2) For each query minimizer, search the index; if the reference minimizer is not
//    in the top -f most frequent set, collect all its positions in the reference
// 3) seeds/anchors are then sorted by reference coordinate and DP-chained
//
// Important:
// - Minimap2's frequency filtering (-f/-U/--q-occ-frac/-e) occurs before occurrences are expanded.
//   "expanding occurrences"; otherwise repeats cause anchor explosion.
//
// Design notes:
// - Use C++ standard library (std::sort, std::unordered_map) instead of minimap2's ksort/khash
// - Keep the same semantics as minimap2: hash matches indicate potential homologous regions
// - Template design: supports any seed type inheriting from SeedHitBase
//
// Performance optimizations:
// - Pre-sort ref_hits by hash
// - Avoid redundant hash calculations (use cached hash values)
// - Pre-allocate output containers to reduce reallocation
// ==================================================================

#include "seed.h"
#include <algorithm>
#include <unordered_map>
#include <cstdint>
#include <vector>
#include <cmath>
#include <limits>

namespace anchor {

// ==================================================================

// ------------------------------------------------------------------
// Function: compute_occ_cutoff_top_frac
// Purpose: Calculate threshold for top f_top_frac frequent minimizers.
//
// Input:
// - occs: occurrence count of each distinct minimizer in reference
// - f_top_frac: top fraction (e.g., 2e-4 means skip most frequent 0.02% minimizers)
//
// Output:
// - Returns occ_cutoff: minimizers with occurrence >= occ_cutoff are considered high-frequency
//   (upper layer decides whether to filter/sparse sample).
//
// Key semantics (align with minimap2):
// - f_top_frac==0: no filter based on top fraction => return +inf
// - occs very small or n_skip==0: means "no filtering" => return +inf
// - f_top_frac>=1: extreme case, filter almost all => return 1
//
// Implementation details:
// - We need the "n_skip-th largest occurrence"; use nth_element for O(N) expected selection.
// - Use std::greater<std::size_t>() to get the n_skip-th largest in descending order.
//
// Complexity:
// - Time: O(N) expected (nth_element), worst O(N log N)
// - Space: O(N) (copy tmp to avoid modifying input)
// ------------------------------------------------------------------
std::size_t compute_occ_cutoff_top_frac(const std::vector<std::size_t>& occs,
                                        double f_top_frac)
{
    if (occs.empty()) return std::numeric_limits<std::size_t>::max();
    if (f_top_frac <= 0.0) return std::numeric_limits<std::size_t>::max();
    if (f_top_frac >= 1.0) return 1; // Extreme case: filter almost all

    // top f fraction of DISTINCT minimizers
    const std::size_t n = occs.size();

    // n_skip represents number of distinct minimizers to skip
    // Example: n=10000, f=2e-4 => n_skip=f*n=2 => skip top 2 most frequent minimizers
    const std::size_t n_skip = static_cast<std::size_t>(std::floor(f_top_frac * static_cast<double>(n)));
    if (n_skip == 0) return std::numeric_limits<std::size_t>::max();

    // nth_element rearranges tmp so that tmp[n_skip-1] is the "n_skip-th largest" element
    std::vector<std::size_t> tmp = occs;
    std::nth_element(tmp.begin(), tmp.begin() + static_cast<std::ptrdiff_t>(n_skip - 1), tmp.end(),
                     std::greater<std::size_t>());
    return tmp[n_skip - 1];
}

// ------------------------------------------------------------------
// Function: compute_ref_occ_threshold
// Purpose: Calculate the final "reference-side occurrence threshold".
//
// Common logic in minimap2 can be understood as:
// - First use -f to estimate a threshold f_cutoff (skip high-frequency minimizers in top fraction)
// - Then use -U to set bounds, constraining threshold to [u_floor, u_ceil] range
// - Final threshold = max(u_floor, min(u_ceil, f_cutoff))
//
// Benefits of this approach:
// - When reference is large and repetitive, -f can adaptively ignore most frequent minimizers
// - Meanwhile -U ensures threshold is not too small (killing normal seeds) or too large (causing repeat explosion)
//
// Complexity:
// - Primarily from compute_occ_cutoff_top_frac's nth_element
// ------------------------------------------------------------------
std::size_t compute_ref_occ_threshold(const std::vector<std::size_t>& occs,
                                      const SeedFilterParams& p)
{
    const std::size_t f_cutoff = compute_occ_cutoff_top_frac(occs, p.f_top_frac);
    const std::size_t capped = std::min(p.u_ceil, f_cutoff);
    return std::max(p.u_floor, capped);
}

// ==================================================================
// Helper function: sort anchors by diagonal (minimap2 chaining preprocessing)
// ==================================================================
//
// Why sort by diagonal?
// - Chaining links anchors "near same diagonal" to form chains.
// - Diagonal (ref_pos - qry_pos) for forward direction.
//   Homologous regions typically cluster on similar diagonals.
//
// Note: reverse chains (is_rev==true) have different diagonal definition:
// - minimap2 internally maps reverse query coords to reverse-complement system.
// - Without query length qlen, cannot directly compute (ref - q_rc).
// - But sorting needs only "monotonically equivalent" keys:
//   q_rc = qlen - (q + span)  => ref - q_rc = ref + q + span - qlen
//   qlen is constant, ignored, so use (ref + q + span) as reverse diagonal key.
//
// Sort keys:
// 1) rid_ref (separate different references)
// 2) is_rev (separate forward/reverse, avoid mixing)
// 3) diagonal key (see above)
// 4) pos_ref / pos_qry (for stability)
//
// Complexity: O(A log A), A=anchors.size()
// ------------------------------------------------------------------

// Seed selection point
// Apply similar selection algorithm
void sortAnchorsByDiagonal(Anchors& anchors)
{
    std::sort(anchors.begin(), anchors.end(),
        [](const Anchor& a, const Anchor& b) {
            if (a.rid_ref != b.rid_ref) return a.rid_ref < b.rid_ref;
            if (a.is_rev != b.is_rev)   return a.is_rev < b.is_rev;

            // Forward: diag = r - q
            // Reverse: equivalent to r - q_rc, where q_rc = qlen - (q + span)
            //          qlen is constant; sorting can use (r + q + span) instead
            const int64_t diag_a = a.is_rev
                ? (static_cast<int64_t>(a.pos_ref) + static_cast<int64_t>(a.pos_qry) + static_cast<int64_t>(a.span))
                : (static_cast<int64_t>(a.pos_ref) - static_cast<int64_t>(a.pos_qry));
            const int64_t diag_b = b.is_rev
                ? (static_cast<int64_t>(b.pos_ref) + static_cast<int64_t>(b.pos_qry) + static_cast<int64_t>(b.span))
                : (static_cast<int64_t>(b.pos_ref) - static_cast<int64_t>(b.pos_qry));
            if (diag_a != diag_b) return diag_a < diag_b;

            if (a.pos_ref != b.pos_ref) return a.pos_ref < b.pos_ref;
            return a.pos_qry < b.pos_qry;
        });
}

// ==================================================================
// Helper function: sort anchors by position
// ==================================================================
void sortAnchorsByPosition(Anchors& anchors)
{
    std::sort(anchors.begin(), anchors.end(),
        [](const Anchor& a, const Anchor& b) {
            if (a.rid_ref != b.rid_ref) return a.rid_ref < b.rid_ref;
            if (a.is_rev != b.is_rev)   return a.is_rev < b.is_rev;
            if (a.pos_ref != b.pos_ref) return a.pos_ref < b.pos_ref;
            return a.pos_qry < b.pos_qry;
        });
}

// ==================================================================
// Helper function: filter high-frequency anchors (post-filtering; not equivalent to minimap2's -f/-U)
// ==================================================================
void filterHighFrequencyAnchors(Anchors& anchors, std::size_t max_occ)
{
    if (anchors.empty() || max_occ == 0) return;

    std::unordered_map<hash_t, std::size_t> hash_count;
    for (const auto& anchor : anchors) {
        hash_count[anchor.hash]++;
    }

    // Note: this is "post-filtering", semantically different from minimap2's -f/-U (which filters before expanding occurrences).
    auto new_end = std::remove_if(anchors.begin(), anchors.end(),
        [&hash_count, max_occ](const Anchor& anchor) {
            return hash_count[anchor.hash] > max_occ;
        });

    anchors.erase(new_end, anchors.end());
}

// ==================================================================
// Chaining-related function implementations
// ==================================================================
//
// References minimap2/lchain.c's mg_lchain_dp function.
//
// Core ideas of minimap2's chaining algorithm:
// 1. Sort anchors by (target_pos, query_pos)
// 2. Use DP to compute the optimal cumulative score f[i] for reaching each anchor
// 3. For each anchor i, consider all possible predecessor anchors j (j < i)
// 4. Conditions for linking two anchors:
//    - On the same reference sequence (rid_ref same)
//    - Same chain direction (is_rev same)
//    - Reference and query distances both within allowed range
//    - Diagonal deviation within bandwidth
// 5. Chain score = base score - penalty
// 6. Use max_skip and max_iter for pruning optimization
// 7. Backtrack to extract all qualifying chains
// ==================================================================

// ------------------------------------------------------------------
// Helper function: compute log2 (used in penalty, see minimap2's mg_log2)
// ------------------------------------------------------------------
static inline float mg_log2(float x) {
    // Use standard library log2; return 0 for small values
    return x >= 1.0f ? std::log2(x) : 0.0f;
}

// ------------------------------------------------------------------
// chainScoreSimple - compute chaining score between two anchors
// ------------------------------------------------------------------
// Implementation details (corresponds to minimap2/lchain.c's comput_sc_simple):
//
// Input convention:
// - ai is the "later" anchor (larger position)
// - aj is the "earlier" anchor (smaller position)
//
// Linkability checks:
// - Must be on the same reference sequence: ai.rid_ref == aj.rid_ref
// - Chain direction must be same: ai.is_rev == aj.is_rev
// - Query-side gap (dq) must be > 0 and <= max_dist_x
// - Reference-side gap (dr) must be > 0 (unless special reverse case)
// - Diagonal deviation |dr - dq| must be <= bw
//
// Score calculation:
// - Base score = min(span_j, min(dq, dr))
// - Diagonal penalty = gap_penalty * |dr - dq| + skip_penalty * min(dr, dq)
// - Log penalty = 0.5 * log2(|dr - dq| + 1)
// - Final score = base score - diagonal penalty -log penalty
// ------------------------------------------------------------------
std::int32_t chainScoreSimple(const Anchor& ai, const Anchor& aj, const ChainParams& params)
{
    // Check if on the same reference sequence and chain direction
    if (ai.rid_ref != aj.rid_ref) return INT32_MIN;
    if (ai.is_rev != aj.is_rev)   return INT32_MIN;

    // Compute gap on query side
    const std::int32_t dq = static_cast<std::int32_t>(ai.pos_qry) - static_cast<std::int32_t>(aj.pos_qry);

    // dq must be > 0 (ai must be after aj) and within allowed range
    if (dq <= 0 || dq > params.max_dist_x) return INT32_MIN;

    // Compute gap on reference side
    const std::int32_t dr = static_cast<std::int32_t>(ai.pos_ref) - static_cast<std::int32_t>(aj.pos_ref);

    // For forward chains, dr must also be > 0
    // For reverse chains, the rules may differ due to coordinate system
    // Simplified here: require dr in reasonable range
    if (dr <= 0 || dr > params.max_dist_y) return INT32_MIN;

    // Compute diagonal deviation
    const std::int32_t dd = (dr > dq) ? (dr - dq) : (dq - dr);

    // Diagonal deviation must be within bandwidth
    if (dd > params.bw) return INT32_MIN;

    // Compute the smaller gap
    const std::int32_t dg = (dr < dq) ? dr : dq;

    // Base score: min(span_j, dg)
    // span represents the length covered by the previous anchor
    const std::int32_t q_span = static_cast<std::int32_t>(aj.span);
    std::int32_t sc = (q_span < dg) ? q_span : dg;

    // Penalty term (only penalize when there is gap deviation or gap > span)
    if (dd > 0 || dg > q_span) {
        // Linear penalty
        const float lin_pen = params.gap_penalty * static_cast<float>(dd)
                            + params.skip_penalty * static_cast<float>(dg);
        // Log penalty (take log of gap deviation)
        const float log_pen = (dd >= 1) ? mg_log2(static_cast<float>(dd + 1)) : 0.0f;

        // Total penalty = linear penalty + 0.5 * log penalty
        sc -= static_cast<std::int32_t>(lin_pen + 0.5f * log_pen);
    }

    return sc;
}

// ------------------------------------------------------------------
// chainAnchors - main chaining function (DP algorithm)
// ------------------------------------------------------------------
// Implementation details (references minimap2/lchain.c mg_lchain_dp):
//
// Algorithm flow:
// 1. Sort anchors by position: (rid_ref, is_rev, pos_ref, pos_qry)
// 2. Initialize DP arrays:
//    - f[i]: maximum chain score ending at anchor i
//    - p[i]: best predecessor for anchor i (-1 if none)
//    - v[i]: peak score before anchor i (for backtracking)
// 3. DP transition:
//    For each anchor i, traverse anchors j before it, compute chaining scores
//    f[i] = max{f[j] + score(j, i)} for all valid j
// 4. Backtrack:
//    Find highest-scoring anchor, backtrack to extract all anchors in chain
// 5. Return best chain's anchors (if meets min_cnt and min_score requirements)
//
// Optimization strategies:
// - max_iter: limit number of predecessors considered per anchor
// - max_skip: terminate early after skipping threshold consecutive anchors
//
// Complexity:
// - Worst O(N^2), but with max_iter/max_skip pruning typically O(N * max_iter)
//
// Notes:
// - Returns best chain's Anchors directly (not Chains structure)
// - Simplified logic, avoids building/managing multiple chains
// ------------------------------------------------------------------
Anchors chainAnchors(Anchors& anchors, const ChainParams& params)
{
    const std::int64_t n = static_cast<std::int64_t>(anchors.size());

    if (n == 0) return Anchors{};

    // ------------------------------------------------------------------
    // Step 1: Sort anchors by position
    // Sort key: (rid_ref, is_rev, pos_ref, pos_qry)
    // This clusters anchors in the same reference region, facilitating DP processing
    // ------------------------------------------------------------------
    std::sort(anchors.begin(), anchors.end(),
        [](const Anchor& a, const Anchor& b) {
            if (a.rid_ref != b.rid_ref) return a.rid_ref < b.rid_ref;
            if (a.is_rev != b.is_rev)   return a.is_rev < b.is_rev;
            if (a.pos_ref != b.pos_ref) return a.pos_ref < b.pos_ref;
            return a.pos_qry < b.pos_qry;
        });

    // ------------------------------------------------------------------
    // Step 2: Allocate DP arrays
    // ------------------------------------------------------------------
    std::vector<std::int32_t> f(n);      // f[i] = maximum chain score ending at anchor i
    std::vector<std::int64_t> p(n);      // p[i] = index of best predecessor for anchor i (-1 if none)
    std::vector<std::int32_t> t(n, 0);   // t[i] = mark array (used for max_skip optimization)

    // ------------------------------------------------------------------
    // Step 3: Fill DP table
    // ------------------------------------------------------------------
    // st: sliding window start position (exclude anchors that are too far)
    std::int64_t st = 0;
    std::int32_t max_f_global = 0;       // global maximum score
    std::int64_t best_end_idx = -1;      // ending anchor index of highest-scoring chain

    for (std::int64_t i = 0; i < n; ++i) {
        const Anchor& ai = anchors[static_cast<std::size_t>(i)];

        // Initialize: f[i] = span (score of anchor alone)
        std::int32_t max_f = static_cast<std::int32_t>(ai.span);
        std::int64_t max_j = -1;
        std::int32_t n_skip = 0;

        // Move window start: exclude anchors that are too far or in different reference sequence
        while (st < i) {
            const Anchor& ast = anchors[static_cast<std::size_t>(st)];
            // Different reference sequence or reference distance exceeds limit, then move window
            if (ast.rid_ref != ai.rid_ref ||
                ast.is_rev != ai.is_rev ||
                static_cast<std::int32_t>(ai.pos_ref - ast.pos_ref) > params.max_dist_x) {
                ++st;
            } else {
                break;
            }
        }

        // Limit iteration count
        std::int64_t iter_start = (i - st > params.max_iter) ? (i - params.max_iter) : st;

        // Iterate backwards through candidate predecessors
        for (std::int64_t j = i - 1; j >= iter_start; --j) {
            const Anchor& aj = anchors[static_cast<std::size_t>(j)];

            // Compute chaining score
            const std::int32_t sc = chainScoreSimple(ai, aj, params);
            if (sc == INT32_MIN) continue;

            // Accumulate predecessor's score
            const std::int32_t total_sc = f[static_cast<std::size_t>(j)] + sc;

            if (total_sc > max_f) {
                max_f = total_sc;
                max_j = j;
                // Found better predecessor, reset skip count
                if (n_skip > 0) --n_skip;
            } else if (t[static_cast<std::size_t>(j)] == static_cast<std::int32_t>(i)) {
                // This anchor has already been visited in current i's traversal
                // If skipped too many times consecutively, terminate early
                if (++n_skip > params.max_skip) break;
            }

            // Mark j's predecessor as visited (used for max_skip optimization)
            if (p[static_cast<std::size_t>(j)] >= 0) {
                t[static_cast<std::size_t>(p[static_cast<std::size_t>(j)])] = static_cast<std::int32_t>(i);
            }
        }

        // Record DP result
        f[static_cast<std::size_t>(i)] = max_f;
        p[static_cast<std::size_t>(i)] = max_j;

        // Update global maximum score and corresponding end index
        if (max_f > max_f_global) {
            max_f_global = max_f;
            best_end_idx = i;
        }
    }

    // ------------------------------------------------------------------
    // Step 4: Backtrack to extract best chain
    // ------------------------------------------------------------------
    // If global maximum score is below threshold, return empty
    if (max_f_global < params.min_score || best_end_idx < 0) {
        return Anchors{};
    }

    // Backtrack to collect anchor indices in chain
    std::vector<std::int64_t> chain_indices;
    std::int64_t cur = best_end_idx;
    while (cur >= 0) {
        chain_indices.push_back(cur);
        cur = p[static_cast<std::size_t>(cur)];
    }

    // Check if chain length meets requirements
    if (static_cast<std::int32_t>(chain_indices.size()) < params.min_cnt) {
        return Anchors{};
    }

    // Reverse so anchors are in position order (front to back)
    std::reverse(chain_indices.begin(), chain_indices.end());

    // Extract anchors to build result
    Anchors result;
    result.reserve(chain_indices.size());
    for (std::int64_t idx : chain_indices) {
        result.push_back(anchors[static_cast<std::size_t>(idx)]);
    }

    return result;
}


} // namespace anchor
