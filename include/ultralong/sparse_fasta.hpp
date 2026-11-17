#pragma once
// PR3D.4D.1: sparse 64-bit out-of-core FASTA index, no alignment or homology claims.
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

namespace ultralong {
struct Checkpoint {
    std::uint64_t nucleotide = 0;    // 0-based coordinate in the read
    std::uint64_t byte_offset = 0;   // 0-based byte offset in source FASTA
};
struct ReadIndex {
    std::string id;
    std::uint64_t bases = 0;
    std::vector<Checkpoint> checkpoints;
};
class SparseFasta {
public:
    static SparseFasta build(const std::filesystem::path& path,
                             std::uint64_t stride = 65536,
                             std::size_t max_reads = 1000);
    const std::vector<ReadIndex>& reads() const noexcept {return reads_;}
    std::uint64_t totalBases() const noexcept {return total_bases_;}
    std::uint64_t checkpointStride() const noexcept {return stride_;}
    const std::filesystem::path& path() const noexcept {return path_;}
    // Per-call DNA allocation <= max_bytes and independent of total FASTA length.
    std::string slice(std::size_t read, std::uint64_t start,
                      std::uint64_t count, std::uint64_t max_bytes = 1048576) const;
private:
    std::filesystem::path path_;
    std::uint64_t stride_ = 0;
    std::uint64_t source_bytes_ = 0;
    std::filesystem::file_time_type source_write_time_;
    std::uint64_t total_bases_ = 0;
    std::vector<ReadIndex> reads_;
    void verifySource() const;
};
} // namespace ultralong
