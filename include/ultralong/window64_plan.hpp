#pragma once
// PR3D.4D.2: bounded out-of-core, 64-bit window candidate geometry.
// NOT an MSA. No oracle. No assertion of true biological homology.
#include <cstddef>
#include <cstdint>
#include <filesystem>
namespace ultralong {
struct Window64Config {
    std::uint64_t window_bases = 65536;
    std::uint64_t flank_bases = 8192;
    std::uint64_t checkpoint_stride = 65536;
    std::size_t max_reads = 8;
    std::uint64_t max_total_bases = 800000000;
    bool rescue_unresolved = false; // opt-in k=10,9 rerun only on unresolved windows
};
struct Window64Stats {
    std::uint64_t total_bases = 0;
    std::uint64_t windows = 0;
    std::uint64_t candidate_pairs = 0;
    std::uint64_t unresolved_pairs = 0;
    std::uint64_t rescued_candidate_pairs = 0;
    std::uint64_t anchored_ref_bases = 0; // summed across candidate query rows
};
// Writes ONLY a TSV of potential windows, never a sequence alignment.
// Refuses overwrites of any output; leaves no published partial output on failure.
Window64Stats planWindows(const std::filesystem::path& input,
                          const std::filesystem::path& output,
                          Window64Config config = {});
} // namespace ultralong
