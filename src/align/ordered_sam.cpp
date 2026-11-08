#include "align/ordered_sam.hpp"

#include <algorithm>
#include <charconv>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>

namespace align::ordered_sam {
namespace {
constexpr std::string_view kTag = "XI:Z:";
}

void addOrdinal(seq_io::SamRecord& record, std::uint64_t ordinal) {
    if (!record.opt.empty()) record.opt += '\t';
    record.opt += std::string(kTag) + std::to_string(ordinal);
}

std::uint64_t getOrdinal(const seq_io::SamRecord& record) {
    const std::string_view tags(record.opt);
    std::size_t pos = 0;
    std::optional<std::uint64_t> result;
    while (pos < tags.size()) {
        const std::size_t end = tags.find('\t', pos);
        const auto part = tags.substr(pos, end == std::string_view::npos ? end : end - pos);
        if (part.starts_with(kTag)) {
            if (result) throw std::runtime_error("Duplicate XI:Z ordinal SAM tag on " + record.qname);
            const auto value = part.substr(kTag.size());
            std::uint64_t ordinal{};
            const auto parsed = std::from_chars(value.data(), value.data() + value.size(), ordinal);
            if (value.empty() || parsed.ec != std::errc{} || parsed.ptr != value.data() + value.size()) {
                throw std::runtime_error("Invalid XI:Z ordinal SAM tag on " + record.qname);
            }
            result = ordinal;
        }
        if (end == std::string_view::npos) break;
        pos = end + 1;
    }
    if (!result) throw std::runtime_error("Missing XI:Z ordinal SAM tag on " + record.qname);
    return *result;
}

Reader::Reader(const std::vector<std::filesystem::path>& paths) {
    readers_.reserve(paths.size());
    last_in_source_.resize(paths.size());
    heap_.reserve(paths.size());
    for (const auto& path : paths) {
        if (!std::filesystem::exists(path)) {
            throw std::runtime_error("Missing thread SAM file: " + path.string());
        }
        readers_.push_back(std::make_unique<seq_io::SamReader>(path));
    }
    for (std::size_t i = 0; i < readers_.size(); ++i) refill(i);
}

void Reader::refill(std::size_t source) {
    seq_io::SamRecord record;
    if (!readers_[source]->next(record)) return;
    const auto ordinal = getOrdinal(record);
    if (last_in_source_[source] && ordinal <= *last_in_source_[source]) {
        throw std::runtime_error("Thread SAM ordinals not increasing: " + record.qname);
    }
    last_in_source_[source] = ordinal;
    heap_.push_back(Head{ordinal, source, std::move(record)});
    std::push_heap(heap_.begin(), heap_.end(), Greater{});
}

bool Reader::next(seq_io::SamRecord& record) {
    if (heap_.empty()) return false;
    std::pop_heap(heap_.begin(), heap_.end(), Greater{});
    Head head = std::move(heap_.back());
    heap_.pop_back();
    if (last_emitted_ && head.ordinal <= *last_emitted_) {
        throw std::runtime_error("Duplicate or descending global SAM ordinal: " +
                                 std::to_string(head.ordinal));
    }
    last_emitted_ = head.ordinal;
    record = std::move(head.record);
    refill(head.source);
    return true;
}

} // namespace align::ordered_sam
