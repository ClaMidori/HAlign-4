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
Selection chooseNativeDirect(const FastaStats& stats, const Budget& budget,
                             std::uint64_t cons_n, bool keep_length,
                             bool has_center_path);

} // namespace align_strategy
