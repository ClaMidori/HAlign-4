#pragma once
// PR3D.1: diagnostic, bounded DNA homology-block planning.
// This API NEVER produces an MSA, uses no oracle/truth labels, and never
// promotes estimated query coordinates to verified homology coordinates.
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

namespace homology_plan {
struct Limits {
    std::size_t max_sequences = 8;
    std::size_t max_per_sequence = 1000000;
    std::size_t max_total_bases = 4000000;
};
struct Config {
    unsigned k = 11;
    std::size_t target_block_bases = 20000;
    std::size_t boundary_radius = 500;
    std::size_t local_span = 10000;
    std::size_t diagonal_tolerance = 32;
    std::size_t min_independent_seeds = 4;
};
struct Sequence {
    std::string id;
    std::string bases;
};
struct Anchor {
    std::uint32_t reference = 0; // zero-based offset in first FASTA sequence
    std::uint32_t query = 0;     // zero-based offset in named query
};
struct Projection {
    std::uint32_t estimated_query = 0;
    std::uint32_t anchor_reference = 0;
    std::uint32_t anchor_query = 0;
    std::uint32_t distance_to_anchor = 0;
};
struct Boundary {
    std::uint32_t reference = 0;
    // Missing means NO confident coordinate estimate at this boundary.
    std::vector<std::optional<Projection>> projections;
    std::size_t supporting_queries = 0;
};
struct Plan {
    Config config;
    std::vector<std::string> ids; // first is deterministic coordinate reference
    std::vector<std::size_t> lengths;
    std::vector<std::vector<Anchor>> supported_anchors; // one vector per query
    std::vector<Boundary> boundaries; // always includes reference 0 and length
};

std::vector<Sequence> readFasta(const std::filesystem::path& input, Limits limits = {});
Plan build(const std::vector<Sequence>& sequences, Config config = {}, Limits limits = {});
std::string asTsv(const Plan& plan);
void writeTsv(const Plan& plan, const std::filesystem::path& output);
} // namespace homology_plan
