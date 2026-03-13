#ifndef HALIGN4_ANCHOR_H
#define HALIGN4_ANCHOR_H

#include <cstdint>
#include <type_traits>
#include <string_view>
#include <stdexcept>
#include <utility>
#include <vector>
#include <string>
#include <cstddef>
#include <algorithm>
#include <limits>
#include <cmath>
#include "hash.h"

// ================================================================
// anchor namespace: Data structures and utility functions related to anchors
// ================================================================
// Design motivation:
// - Anchor and its associated filtering/sorting functions are essentially public data structures for the "anchor/chaining stage",
//   and should not be tightly coupled with the seed abstraction interface (SeedHitBase/traits) in the same namespace.
// - Splitting them into anchor:: allows seed:: to focus only on "abstraction and extraction of seed/hit",
//   while anchor:: focuses on "how to form anchors from hits, how to sort/filter".
//
// Important:
// - This still uses the global hash_t (from include/hash.h) to stay consistent with existing code.
// - Currently, minimizer::collect_anchors returns anchor::Anchors; in the future, if syncmer/strobemer are supported,
//   the same set of anchor utility functions can be reused.
// ================================================================
namespace anchor
{
    // ------------------------------------------------------------------
    // Struct: Anchor
    // ------------------------------------------------------------------
    // Semantics: Describes an "anchor match" between ref and query on a certain seed/hash.
    // Downstream chaining will string a group of anchors into a chain (candidate alignment region).
    // ------------------------------------------------------------------
    struct Anchor
    {
        hash_t hash{};              // seed hash (expected to be the same for ref/query)
        std::uint32_t rid_ref{};    // ref sequence id
        std::uint32_t pos_ref{};    // position on ref (0-based)
        std::uint32_t rid_qry{};    // query sequence id (often fixed to 0 in many scenarios)
        std::uint32_t pos_qry{};    // position on query (0-based, forward coordinate system)
        std::uint32_t span{};       // coverage length (can use min(ref.span, qry.span))
        bool is_rev{};              // whether ref/query are on "opposite strands"
        // Optional: pre-cache diagonal, commonly used in chaining
        // int32_t diag{}; // (int32_t)pos_ref - (int32_t)pos_qry
    };

    using Anchors = std::vector<Anchor>;

    // ==================================================================
    // Internal helper structure: hash index for ref_hits
    // ==================================================================
    struct HashIndex {
        std::size_t start;  // Starting index in the sorted array
        std::size_t count;  // Number of elements with this hash
    };

    // ==================================================================
    // minimap2-style seeding filter parameters (default values same as minimap2 CLI)
    // ==================================================================
    struct SeedFilterParams {
        double f_top_frac = 2e-4;                 // -f
        std::size_t u_floor = 10;                 // -U lower
        std::size_t u_ceil  = 1000000;            // -U upper
        double q_occ_frac  = 0.01;                // --q-occ-frac
        std::size_t sample_every_bp = 500;        // -e
    };

    static inline SeedFilterParams default_mm2_params() {
        return SeedFilterParams{};
    }

    // Compute the occurrence threshold for -f (fraction): ignore the top f_top_frac most frequent minimizers
    // Return value: occ_cutoff (>=1). If there are very few distinct minimizers or f_top_frac==0, return +inf.
    std::size_t compute_occ_cutoff_top_frac(const std::vector<std::size_t>& occs,
                                            double f_top_frac);

    // Compute the final reference occurrence threshold: max{u_floor, min{u_ceil, -f}}
    std::size_t compute_ref_occ_threshold(const std::vector<std::size_t>& occs,
                                          const SeedFilterParams& p);

    // =====================================================================
    // sortAnchorsByDiagonal - Sort anchors by diagonal (for chaining algorithm)
    // =====================================================================
    void sortAnchorsByDiagonal(Anchors& anchors);

    // =====================================================================
    // sortAnchorsByPosition - Sort anchors by position
    // =====================================================================
    void sortAnchorsByPosition(Anchors& anchors);

    // =====================================================================
    // filterHighFrequencyAnchors - Filter high-frequency anchors (refer to minimap2)
    // =====================================================================
    void filterHighFrequencyAnchors(Anchors& anchors, std::size_t max_occ = 500);

    // =====================================================================
    // Chaining-related data structures and functions
    // =====================================================================
    // Refer to minimap2/lchain.c implementation, use dynamic programming (DP) to string anchors into a chain.
    //
    // Core idea:
    // - After sorting anchors by reference position, use DP to compute the optimal score to reach each anchor
    // - The score between two anchors = min(gap_ref, gap_qry, span) - penalty
    // - The penalty considers gap difference (diagonal offset) and gap size
    // - Backtrack to find the highest scoring chain
    // =====================================================================

    // ------------------------------------------------------------------
    // Chaining parameters (refer to minimap2 default values)
    // ------------------------------------------------------------------
    struct ChainParams {
        std::int32_t max_dist_x = 5000;       // Maximum distance in the direction of the reference sequence
        std::int32_t max_dist_y = 5000;       // Maximum distance in query sequence direction
        std::int32_t bw = 500;                // Bandwidth (diagonal offset tolerance)
        std::int32_t max_skip = 25;           // Maximum skip count (optimization)
        std::int32_t max_iter = 5000;         // Maximum iteration count (optimization)
        std::int32_t min_cnt = 3;             // Minimum anchor count for a chain
        std::int32_t min_score = 40;          // Minimum score for a chain
        float gap_penalty = 0.01f;            // gap penalty coefficient
        float skip_penalty = 0.01f;           // skip (gap size) penalty coefficient
    };

    // Return default chaining parameters
    inline ChainParams default_chain_params() {
        return ChainParams{};
    }

    // ------------------------------------------------------------------
    // chainScoreSimple - Compute the chaining score between two anchors (internal helper function)
    // ------------------------------------------------------------------
    // Input:
    //   ai    : current anchor (larger position)
    //   aj    : previous anchor (smaller position)
    //   params: chaining parameters
    //
    // Output:
    //   Returns chaining score, returns INT32_MIN if the two anchors cannot be linked
    //
    // Score calculation (refer to minimap2):
    //   base score = min(gap_ref, gap_qry, span)
    //   penalty = gap_penalty * |gap_ref - gap_qry| + skip_penalty * min(gap_ref, gap_qry)
    //   final score = base score - penalty - 0.5 * log2(|gap_ref - gap_qry| + 1)
    //
    // Note: This function is for internal use by chainAnchors only
    // ------------------------------------------------------------------
    std::int32_t chainScoreSimple(const Anchor& ai, const Anchor& aj, const ChainParams& params);

    // ------------------------------------------------------------------
    // chainAnchors - Use DP to chain anchors and return the best chain
    // ------------------------------------------------------------------
    // Function:
    // Uses dynamic programming algorithm to chain anchors and find the highest scoring anchor chain.
    // Refer to minimap2/lchain.c's mg_lchain_dp implementation.
    //
    // Input:
    //   anchors : list of anchors (will be sorted and modified)
    //   params  : chaining parameters
    //
    // Output:
    //   Returns the list of anchors in the best chain (sorted by position)
    //   If no chain meets the criteria, returns empty Anchors
    //
    // Algorithm steps:
    // 1. Sort anchors by (rid_ref, is_rev, pos_ref, pos_qry)
    // 2. Use DP to compute the optimal link score for each anchor
    // 3. Backtrack to extract the highest scoring chain
    // 4. Check if the chain meets min_cnt and min_score conditions
    // 5. Return the list of anchors in the best chain
    //
    // Note:
    //   - Input anchors will be sorted
    //   - Returned anchors are sorted by their position in the chain (from start to end)
    //   - If there are multiple chains with the same score, returns the first one
    // ------------------------------------------------------------------
    Anchors chainAnchors(Anchors& anchors, const ChainParams& params = default_chain_params());

} // namespace anchor


#endif //HALIGN4_ANCHOR_H