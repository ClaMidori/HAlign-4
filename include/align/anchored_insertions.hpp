#pragma once

#include "align/native_msa.hpp"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace align::anchored_insertions {

// This is an opt-in, bounded experimental path for short hierarchical DNA.
// Input SAM alignments must be global to the UNGAPPED consensus, and queries
// must be in original FASTA order. No evolutionary truth is used.
struct Query {
    std::string id;
    std::string bases;
    std::string cigar;
    std::uint64_t ordinal = 0;
};

struct Options {
    std::size_t max_sequences = 64;
    std::size_t max_center_bases = 25000;
    std::size_t max_total_bases = 100000;
    std::size_t max_output_cells = 8000000;
    bool soft_anchors = false;           // opt-in, conservative center-block rescue
    bool soft_anchor_recruit = false;   // opt-in V3 fallback when V2 rejects
    std::size_t max_soft_blocks = 8;     // fall back when too many candidates
};

struct Alignment {
    std::vector<native_msa::Sequence> rows; // center first, then queries in FASTA order
    std::size_t columns = 0;
    std::size_t occupied_sites = 0;
    std::size_t local_msas = 0;
    std::size_t soft_anchor_candidates = 0;
    std::size_t soft_anchor_splits = 0;
    std::size_t soft_anchor_recruits = 0;
};

// Each CIGAR's M/=/X and D/N consume reference bases, I consumes query bases.
// Site s is immediately before consensus base s (including the trailing site).
// All other SAM operations and partial/non-global CIGARs are rejected.
// Local insertion MSAs are performed with the pinned, in-process abPOA library.
Alignment align(const std::string& consensus,
                const std::vector<Query>& queries,
                const Options& limits = {});

} // namespace align::anchored_insertions
