#pragma once
// PR3D.3A: audit-only, lossless FASTA tiling of bounded PR3D.2C TSV plans.
// Candidate projections are NOT biologically verified homologies. No MSA.
#include <cstddef>
#include <filesystem>
#include <string>
#include <string_view>

namespace lossless_extractor {
// Digest primitives reused by the PR3D.3B integrity audit.
std::string sha256Digest(std::string_view bytes);
std::string sha256FileDigest(const std::filesystem::path& path);
struct Summary {
    std::size_t reads = 0;
    std::size_t blocks = 0;
    std::size_t accepted_query_candidates = 0;
    std::size_t rejected_conflicting_candidates = 0;
    std::size_t tiles = 0;
    std::size_t unresolved_query_bases = 0;
};

// Input FASTA and planner TSV are never changed. Output must not exist.
// All outputs are staged in a sibling directory and published by rename.
Summary extract(const std::filesystem::path& fasta,
                const std::filesystem::path& plan_tsv,
                const std::filesystem::path& output_dir);
} // namespace lossless_extractor
