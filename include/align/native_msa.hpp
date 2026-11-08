#pragma once

#include <cstddef>
#include <string>
#include <vector>

namespace native_msa {

struct Sequence {
    std::string id;
    std::string bases;
};

struct Options {
    // Corresponds to abpoa -M 2 -X 2 -S (global, seeded).
    int match = 2;
    int mismatch = 2;
    bool seeding = true;
    // Initial short-sequence-only safety budget. The future block planner will
    // control each call independently, rather than submitting full ultra-long DNA.
    std::size_t max_sequence_bases = 25000;
    std::size_t max_total_bases = 100000;
    std::size_t max_sequences = 64;
};

struct Alignment {
    std::vector<Sequence> rows;
    std::size_t columns = 0;
};

// Build a local row-column MSA in-process through the abPOA C API.
// Only uppercase A/C/G/T/N are accepted in this first integration.
// Throws std::invalid_argument / std::runtime_error on invalid input or output.
Alignment align(const std::vector<Sequence>& sequences, const Options& options = {});

}  // namespace native_msa
