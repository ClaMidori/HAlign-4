#pragma once
// PR3D.4D.5B. Strict, conservative sequence-only rescue geometry.
// Synthetic truth is intentionally absent from all signatures.
#include "align/homology_plan.hpp"
#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>
namespace ultralong {
struct RescueRange { std::uint32_t query_start=0, query_end=0; std::size_t support=0; double sampled_identity=0; };
// Explicitly fail closed for sparse, inconsistent, repetitive or divergent evidence.
std::optional<RescueRange> projectRescueRange(
    const std::vector<homology_plan::Anchor>& anchors,
    const std::string& ref, const std::string& query,
    std::uint32_t ref_start, std::uint32_t ref_end,
    double min_sample_identity=0.53);
} // namespace ultralong
