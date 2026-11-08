#pragma once

#include <cstddef>
#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

namespace align::profile_projection {

// Two MSAs are anchored on the same ungapped consensus, but may put gaps
// at different columns. The plan embeds both profiles in one column system.
struct Plan {
    std::size_t columns = 0;
    std::vector<std::size_t> reference_to_common;
    std::vector<std::size_t> insertion_to_common;
};

// Fails rather than guessing if the two anchors disagree at any residue.
Plan makePlan(std::string_view reference_anchor,
              std::string_view insertion_anchor,
              std::string_view expected_ungapped_anchor = {});

// Preserves all characters and their order. Newly introduced columns are '-'.
std::string projectRow(std::string_view row,
                       const std::vector<std::size_t>& old_to_common,
                       std::size_t common_columns);

// Verify the two rectangular FASTA profiles, then project each into the same
// coordinate system. The first record in both files must be the consensus.
// Output files are new paths; input files are not modified.
// The reference first-row ID can vary; its ungapped bases must match the
// expected consensus. The insertion first row must be named consensus.
std::size_t reconcileFastaProfiles(
    const std::filesystem::path& reference_fasta,
    const std::filesystem::path& insertion_fasta,
    const std::filesystem::path& projected_reference_fasta,
    const std::filesystem::path& projected_insertion_fasta,
    std::string_view expected_ungapped_consensus);

}  // namespace align::profile_projection
