#pragma once

#include <cstdint>
#include <filesystem>
#include <string>

namespace align_strategy {

// Statistics are computed using a fixed-size streaming buffer.
// No entire FASTA record is retained by the scanner.
struct FastaStats {
    std::uint64_t sequence_count = 0;
    std::uint64_t total_bases = 0;
    std::uint64_t longest_bases = 0;
    std::uint64_t shortest_bases = 0;
    bool contains_gaps = false;
};

struct Budget {
    std::uint64_t max_sequences = 64;
    std::uint64_t max_sequence_bases = 25000;
    std::uint64_t max_total_bases = 100000;
};

struct Selection {
    bool native_direct = false;
    std::string explanation;
};

FastaStats scanFasta(const std::filesystem::path& path);
// Throws if any FASTA row differs in aligned column count. Counts gaps as columns.
FastaStats scanRectangularFasta(const std::filesystem::path& path);

// PR3C.2. Opt-in bounded auto router. No I/O, allocations or heuristics.
// Its input comes from the existing fixed-buffer FASTA scanner, before
// preprocessing materializes the Top-K sequences in memory.
enum class AutoBackend { NativeDirect, NativeHierarchical, Reject };
struct AutoSelection {
    AutoBackend backend = AutoBackend::Reject;
    std::string explanation;
};

inline AutoSelection chooseBoundedAuto(const FastaStats& s, const Budget& b,
                                       std::uint64_t cons_n) {
    // The anchored hierarchical backend shares the native per-input limits.
    // Never substitute a different backend when a budget is exceeded.
    if (!s.sequence_count) return {AutoBackend::Reject, "no input sequences"};
    if (cons_n == 0) return {AutoBackend::Reject, "--cons-n must be positive"};
    if (s.contains_gaps) return {AutoBackend::Reject, "aligned gaps are not supported in bounded auto"};
    if (s.sequence_count > b.max_sequences || s.sequence_count > 64)
        return {AutoBackend::Reject, "sequence count exceeds 64-query budget"};
    if (s.longest_bases > b.max_sequence_bases || s.longest_bases > 25000)
        return {AutoBackend::Reject, "per-sequence length exceeds 25kb budget"};
    if (s.total_bases > b.max_total_bases || s.total_bases > 100000)
        return {AutoBackend::Reject, "input bases exceed 100kb total budget"};
    if (s.sequence_count <= cons_n)
        return {AutoBackend::NativeDirect, "all inputs fit native direct MSA"};
    if (cons_n < 2 || cons_n > 64)
        return {AutoBackend::Reject, "hierarchical --cons-n must be 2..64"};
    return {AutoBackend::NativeHierarchical,
            "bounded Top-K native center plus anchored insertions and Soft Anchors V3"};
}

Selection chooseNativeDirect(const FastaStats& stats, const Budget& budget,
                             std::uint64_t cons_n, bool keep_length,
                             bool has_center_path);

} // namespace align_strategy
