#include "align/strategy.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>

#if defined(HALIGN4_HAVE_ZLIB) && HALIGN4_HAVE_ZLIB
#include <zlib.h>
#endif

namespace align_strategy {
namespace {
// Unlike the kseq-based scanner in PR3A, this parser never retains a DNA
// record. Input memory is bounded by the buffer, regardless of record length.
constexpr std::size_t kBufferBytes = 65536;
constexpr std::size_t kMaxIdentifierBytes = 1048576;

class InputStream {
public:
    explicit InputStream(const std::filesystem::path& path) : path_(path.string()) {
#if defined(HALIGN4_HAVE_ZLIB) && HALIGN4_HAVE_ZLIB
        file_ = gzopen(path_.c_str(), "rb");
        if (!file_) throw std::runtime_error("Cannot open FASTA for scanning: " + path_);
#else
        file_.open(path, std::ios::binary);
        if (!file_) throw std::runtime_error("Cannot open FASTA for scanning: " + path_);
        char magic[2] = {};
        file_.read(magic, 2);
        if (file_.gcount() == 2 && static_cast<unsigned char>(magic[0]) == 0x1f &&
            static_cast<unsigned char>(magic[1]) == 0x8b) {
            throw std::runtime_error("Compressed FASTA requires HAlign-4 built with ZLIB");
        }
        file_.clear();
        file_.seekg(0);
#endif
    }
    InputStream(const InputStream&) = delete;
    InputStream& operator=(const InputStream&) = delete;
    ~InputStream() {
#if defined(HALIGN4_HAVE_ZLIB) && HALIGN4_HAVE_ZLIB
        if (file_) gzclose(file_);
#endif
    }
    std::size_t read(char* dest, std::size_t count) {
#if defined(HALIGN4_HAVE_ZLIB) && HALIGN4_HAVE_ZLIB
        const int n = gzread(file_, dest, static_cast<unsigned int>(count));
        if (n < 0) {
            int error_number = 0;
            const char* message = gzerror(file_, &error_number);
            throw std::runtime_error("FASTA read error: " + path_ + ": " +
                                     (message ? message : "unknown gzip error"));
        }
        return static_cast<std::size_t>(n);
#else
        file_.read(dest, static_cast<std::streamsize>(count));
        const auto n = file_.gcount();
        if (file_.bad()) throw std::runtime_error("FASTA read error: " + path_);
        return static_cast<std::size_t>(n);
#endif
    }
private:
    std::string path_;
#if defined(HALIGN4_HAVE_ZLIB) && HALIGN4_HAVE_ZLIB
    gzFile file_ = nullptr;
#else
    std::ifstream file_;
#endif
};

class FastaCounter {
public:
    explicit FastaCounter(bool require_rectangular = false)
        : require_rectangular_(require_rectangular) {}

    void accept(char value) {
        if (value == '\r' || value == '\n') {
            line_start_ = true;
            header_ = false;
            return;
        }
        if (line_start_) {
            line_start_ = false;
            if (value == '>') {
                if (seen_) finishRecord();
                seen_ = true;
                header_ = true;
                identifier_chars_ = 0;
                in_identifier_ = true;
                current_bases_ = 0;
                return;
            }
            if (!seen_) {
                if (std::isspace(static_cast<unsigned char>(value))) return;
                throw std::runtime_error("FASTA sequence data precedes first header");
            }
        }
        if (header_) {
            if (in_identifier_) {
                if (std::isspace(static_cast<unsigned char>(value))) {
                    in_identifier_ = false;
                } else {
                    if (++identifier_chars_ > kMaxIdentifierBytes)
                        throw std::runtime_error("FASTA identifier is excessively long");
                }
            }
            return;
        }
        if (!std::isspace(static_cast<unsigned char>(value))) {
            if (value == '-') stats_.contains_gaps = true;
            if (current_bases_ == std::numeric_limits<std::uint64_t>::max())
                throw std::overflow_error("FASTA sequence length exceeds 64-bit limit");
            ++current_bases_;
        }
    }
    FastaStats result() {
        if (!seen_) throw std::runtime_error("Automatic alignment requires non-empty FASTA input");
        finishRecord();
        return stats_;
    }
private:
    void finishRecord() {
        if (identifier_chars_ == 0 || current_bases_ == 0)
            throw std::runtime_error("Automatic alignment requires non-empty FASTA IDs and sequences");
        if (stats_.sequence_count == std::numeric_limits<std::uint64_t>::max() ||
            current_bases_ > std::numeric_limits<std::uint64_t>::max() - stats_.total_bases)
            throw std::overflow_error("FASTA statistics exceed 64-bit limits");
        if (require_rectangular_ && stats_.sequence_count > 0 &&
            current_bases_ != first_record_bases_) {
            throw std::runtime_error("Non-rectangular aligned FASTA: row " +
                                     std::to_string(stats_.sequence_count + 1) +
                                     " has " + std::to_string(current_bases_) +
                                     " columns, expected " +
                                     std::to_string(first_record_bases_));
        }
        if (stats_.sequence_count == 0) first_record_bases_ = current_bases_;
        ++stats_.sequence_count;
        stats_.total_bases += current_bases_;
        stats_.longest_bases = std::max(stats_.longest_bases, current_bases_);
        stats_.shortest_bases = stats_.sequence_count == 1
            ? current_bases_ : std::min(stats_.shortest_bases, current_bases_);
    }
    FastaStats stats_{};
    bool seen_ = false;
    bool line_start_ = true;
    bool header_ = false;
    bool in_identifier_ = false;
    std::size_t identifier_chars_ = 0;
    std::uint64_t current_bases_ = 0;
    std::uint64_t first_record_bases_ = 0;
    bool require_rectangular_ = false;
};
} // namespace

namespace {
FastaStats scanFastaImpl(const std::filesystem::path& path, bool require_rectangular) {
    InputStream stream(path);
    FastaCounter scanner(require_rectangular);
    std::array<char, kBufferBytes> buffer{};
    for (;;) {
        const auto n = stream.read(buffer.data(), buffer.size());
        if (!n) break;
        for (std::size_t i = 0; i < n; ++i) scanner.accept(buffer[i]);
    }
    return scanner.result();
}
} // namespace

FastaStats scanFasta(const std::filesystem::path& path) {
    return scanFastaImpl(path, false);
}

FastaStats scanRectangularFasta(const std::filesystem::path& path) {
    return scanFastaImpl(path, true);
}

Selection chooseNativeDirect(const FastaStats& s, const Budget& b,
                             std::uint64_t cons_n, bool keep_length,
                             bool has_center_path) {
    if (keep_length) return {false, "--keep-length requires the hierarchical path"};
    if (has_center_path) return {false, "--center-path requires the hierarchical path"};
    if (!s.sequence_count) return {false, "no input sequences"};
    if (s.sequence_count > b.max_sequences) return {false, "sequence count exceeds native MSA limit"};
    if (s.sequence_count > cons_n) return {false, "sequence count exceeds --cons-n; native preprocessing would omit sequences"};
    if (s.longest_bases > b.max_sequence_bases) return {false, "a sequence exceeds the native per-sequence limit"};
    if (s.total_bases > b.max_total_bases) return {false, "total bases exceed the native MSA limit"};
    return {true, "all input sequences fit the bounded native direct MSA"};
}

} // namespace align_strategy
