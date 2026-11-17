#pragma once

#include <cstddef>
#include <filesystem>

namespace ultralong_native {
struct Config {
    std::size_t max_total_bases = 8000000;
    std::size_t max_local_bases = 12000;
    double min_reference_identity = 0.55;
};
struct Summary {
    std::size_t original_bases = 0;
    std::size_t output_columns = 0;
    std::size_t accepted_blocks = 0;
    std::size_t aligned_reference_bases = 0;
    std::size_t candidate_query_rows = 0;
    std::size_t skipped_budget = 0;
    std::size_t rejected_postcheck = 0;
    std::size_t private_bases = 0;
};
// Only for short, memory-resident, colinear DNA. This is NOT out-of-core.
// Publishes new output+sidecar only after reconstructing all source bases.
Summary run(const std::filesystem::path& input,
            const std::filesystem::path& output,
            const std::filesystem::path& report,
            Config cfg = {});
} // namespace ultralong_native
