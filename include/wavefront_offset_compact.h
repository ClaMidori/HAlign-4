#ifndef HALIGN4_WAVEFRONT_OFFSET_COMPACT_H
#define HALIGN4_WAVEFRONT_OFFSET_COMPACT_H

#include <cstdint>
#include <vector>
#include <stdexcept>
#include <spdlog/spdlog.h>

// ========================================================================
// Wavefront Offset Compaction: Reduce memory footprint from int32 to int16
// ========================================================================
// Strategy:
//   - Store first offset as int32 (baseline)
//   - Store subsequent offsets as int16 deltas (relative to previous)
//   - Savings: 50% reduction (4 bytes → 2 bytes per offset after first)
//
// Correctness:
//   - Offsets between adjacent diagonals differ by ~O(1) in typical alignments
//   - int16 range [-32768, 32767] is sufficient for all practical sequences
// ========================================================================

namespace wfa_compact {

using wf_offset_t = int32_t;           // Original type
using wf_offset_delta_t = int16_t;     // Compact delta type

// ========================================================================
// Convert absolute offsets to compact relative representation
// ========================================================================
// Input: vector of int32 offsets (lo to hi)
// Output: baseline (first offset) and vector of int16 deltas
// ========================================================================
inline void offsets_to_compact(
    const wf_offset_t* offsets,
    int num_offsets,
    wf_offset_t& out_baseline,
    std::vector<wf_offset_delta_t>& out_deltas)
{
    if (num_offsets <= 0) {
        out_baseline = 0;
        out_deltas.clear();
        return;
    }

    out_baseline = offsets[0];
    out_deltas.resize(num_offsets - 1);

    // Compute deltas and validate range
    for (int i = 1; i < num_offsets; ++i) {
        wf_offset_t delta_full = offsets[i] - offsets[i - 1];

        // Check if delta fits in int16
        if (delta_full < std::numeric_limits<wf_offset_delta_t>::min() ||
            delta_full > std::numeric_limits<wf_offset_delta_t>::max()) {
            spdlog::error(
                "offsets_to_compact: delta {} at offset {} exceeds int16 range [{}, {}]",
                delta_full, i,
                std::numeric_limits<wf_offset_delta_t>::min(),
                std::numeric_limits<wf_offset_delta_t>::max());
            throw std::runtime_error("Offset delta exceeds int16 range");
        }

        out_deltas[i - 1] = static_cast<wf_offset_delta_t>(delta_full);
    }
}

// ========================================================================
// Convert compact relative representation back to absolute offsets
// ========================================================================
// Input: baseline and vector of int16 deltas
// Output: reconstructed vector of int32 offsets
// ========================================================================
inline std::vector<wf_offset_t> compact_to_offsets(
    wf_offset_t baseline,
    const std::vector<wf_offset_delta_t>& deltas)
{
    std::vector<wf_offset_t> result;
    
    // Always include baseline as first element (even if no deltas)
    result.reserve(deltas.size() + 1);
    result.push_back(baseline);

    wf_offset_t current = baseline;
    for (wf_offset_delta_t delta : deltas) {
        current += static_cast<wf_offset_t>(delta);
        result.push_back(current);
    }

    return result;
}

// ========================================================================
// Validation: Round-trip test
// ========================================================================
inline bool validate_round_trip(
    const std::vector<wf_offset_t>& original)
{
    if (original.empty()) return true;

    wf_offset_t baseline;
    std::vector<wf_offset_delta_t> deltas;

    offsets_to_compact(original.data(), original.size(), baseline, deltas);
    auto recovered = compact_to_offsets(baseline, deltas);

    if (recovered.size() != original.size()) {
        spdlog::error("Round-trip validation failed: size mismatch {} vs {}",
                     recovered.size(), original.size());
        return false;
    }

    for (size_t i = 0; i < original.size(); ++i) {
        if (recovered[i] != original[i]) {
            spdlog::error("Round-trip validation failed at offset {}: {} vs {}",
                         i, recovered[i], original[i]);
            return false;
        }
    }

    return true;
}

// ========================================================================
// Estimate memory savings
// ========================================================================
// Returns (original_bytes, compact_bytes, savings_percent)
// ========================================================================
inline void estimate_savings(
    int num_offsets,
    double& out_original_bytes,
    double& out_compact_bytes,
    double& out_savings_percent)
{
    if (num_offsets <= 0) {
        out_original_bytes = 0;
        out_compact_bytes = 0;
        out_savings_percent = 0;
        return;
    }

    // Original: all int32
    out_original_bytes = num_offsets * sizeof(wf_offset_t);

    // Compact: first offset (int32) + rest (int16)
    out_compact_bytes = sizeof(wf_offset_t) + (num_offsets - 1) * sizeof(wf_offset_delta_t);

    // Savings
    if (out_original_bytes > 0) {
        out_savings_percent = 100.0 * (1.0 - out_compact_bytes / out_original_bytes);
    } else {
        out_savings_percent = 0;
    }
}

} // namespace wfa_compact

#endif // HALIGN4_WAVEFRONT_OFFSET_COMPACT_H
