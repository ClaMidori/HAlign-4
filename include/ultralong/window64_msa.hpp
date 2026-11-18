#pragma once
#include <cstddef>
#include <cstdint>
#include <filesystem>

namespace ultralong {
struct WindowMsaConfig {
    std::uint64_t max_total_bases = 800000000ULL;
    std::size_t max_reads = 8;
    std::size_t local_max_bases = 12000;
    std::size_t local_max_total_bases = 100000;
    double min_reference_identity = 0.55;
    bool rescue_unresolved = false; // explicit opt-in; baseline unchanged
};
struct WindowMsaStats {
    std::uint64_t original_bases=0, windows=0, candidate_query_windows=0;
    std::uint64_t local_candidate_rows=0, aligned_blocks=0, aligned_reference_bases=0;
    std::uint64_t rejected_postcheck=0, skipped_budget=0, private_bases=0;
    std::uint64_t rescued_window_candidates=0, rescued_local_rows=0;
    std::uint64_t output_columns=0;
};
// Actual out-of-core MSA: indexed FASTA windows -> local homology planner
// -> bounded real abPOA -> on-disk row spools -> validated FASTA.
// No synthetic ground truth is ever loaded.
WindowMsaStats runWindowMsa(const std::filesystem::path& input,
                           const std::filesystem::path& output,
                           const std::filesystem::path& report,
                           WindowMsaConfig cfg = {});
} // namespace ultralong
