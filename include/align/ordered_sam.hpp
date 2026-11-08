#pragma once

#include "utils.h"

#include <cstdint>
#include <filesystem>
#include <memory>
#include <optional>
#include <vector>

namespace align::ordered_sam {

// Private SAM tag XI:Z: stores the input FASTA ordinal, not the QNAME.
// Using a string-valued SAM tag permits full uint64_t ordinals.
void addOrdinal(seq_io::SamRecord& record, std::uint64_t ordinal);
std::uint64_t getOrdinal(const seq_io::SamRecord& record);

// K-way merge: every thread-specific SAM stream must be ordered by input
// ordinal. At most one SAM record per source is resident in the heap.
class Reader {
public:
    explicit Reader(const std::vector<std::filesystem::path>& paths);
    Reader(const Reader&) = delete;
    Reader& operator=(const Reader&) = delete;

    bool next(seq_io::SamRecord& record);

private:
    struct Head {
        std::uint64_t ordinal;
        std::size_t source;
        seq_io::SamRecord record;
    };
    struct Greater {
        bool operator()(const Head& a, const Head& b) const noexcept {
            return a.ordinal > b.ordinal;
        }
    };

    void refill(std::size_t source);
    std::vector<std::unique_ptr<seq_io::SamReader>> readers_;
    std::vector<std::optional<std::uint64_t>> last_in_source_;
    std::optional<std::uint64_t> last_emitted_;
    std::vector<Head> heap_;
};

} // namespace align::ordered_sam
