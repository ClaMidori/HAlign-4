#ifndef HALIGN4_MSA_VALIDATION_HPP
#define HALIGN4_MSA_VALIDATION_HPP

#include <cstdint>
#include <filesystem>

namespace msa {
struct ValidationStats {
    std::uint64_t records = 0;
    std::uint64_t columns = 0;
    std::uint64_t bases = 0;
};

// Validate an aligned FASTA against its original (cleaned) unaligned FASTA.
// Checks IDs, duplicates, equal column count and per-record ungapped
// length + two independent 64-bit fingerprints. Uses O(number of IDs)
// auxiliary memory rather than storing entire input sequences.
// Throws std::runtime_error if invalid. This is an integrity guard, not
// a biological homology/accuracy metric; fingerprints are not cryptographic.
ValidationStats validateFastaAlignment(
    const std::filesystem::path& input,
    const std::filesystem::path& aligned);
}  // namespace msa

#endif
