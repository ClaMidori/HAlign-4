#pragma once
// PR 3D.3B: experimental abPOA MSAs *within* candidate-only PR3D.3A tiles.
// Never produces a global MSA and never treats a candidate as verified homology.
#include <cstddef>
#include <filesystem>

namespace block_local_msa {
struct Summary {
    std::size_t blocks = 0;
    std::size_t aligned = 0;
    std::size_t skipped_no_query = 0;
    std::size_t skipped_budget = 0;
    std::size_t aligned_rows = 0;
    std::size_t unresolved_tiles = 0;
};
// original_fasta must be exactly the source of the existing extraction.
// extraction_dir contains the five PR3D.3A files; output_dir must not exist.
// The aligned output is block-local only. UNRESOLVED tiles are never aligned.
Summary align(const std::filesystem::path& original_fasta,
              const std::filesystem::path& extraction_dir,
              const std::filesystem::path& output_dir);
}
