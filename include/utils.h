#ifndef HALIGN4_UTILS_H
#define HALIGN4_UTILS_H

// ================================================================
// utils.h - File and sequence I/O public declarations (detailed English comments)
//
// This header file contains two main functional blocks:
//  1) file_io: A set of utility functions related to file/path/network fetching (copy/download/fetch/prep directory, etc.),
//     The purpose is to centralize all file system interaction logic for unified error handling and testing.
//  2) seq_io: A set of FASTA/FASTQ read/write abstractions based on kseq (KseqReader / FastaWriter / SeqRecord),
//     And provides sequence cleaning functions (cleanSequence) and convenient reader/writer constructors.
// Only interfaces are declared here; specific implementations are located in src/utils/file_io.cpp and src/utils/seq_io.cpp.
// ================================================================

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <system_error>
#include <regex>
#include <array>
#include <cstdint>

#if __has_include(<curl/curl.h>)
#include <curl/curl.h>
#define HALIGN4_HAVE_LIBCURL 1
#endif

namespace fs = std::filesystem;
using FilePath = std::filesystem::path;


namespace file_io
{
    // Below are the function declarations provided in file_io, with Chinese explanations added before each function to help users understand semantics and boundaries.

    // Format filesystem error information into a readable string (with path and error message)
    // Usage: Call when throwing exceptions or logging, to unify error output format
    std::string formatFsError(std::string_view msg,
                                     const FilePath& p,
                                     const std::error_code& ec);

    // Require path to exist (any type); throw runtime_error if not exist or check fails
    void requireExists(const FilePath& p, std::string_view what);

    // Require path to be a regular file; throw runtime_error if not
    void requireRegularFile(const FilePath& p, std::string_view what);

    // Require path to be a directory; throw runtime_error if not
    void requireDirectory(const FilePath& p, std::string_view what);

    // Ensure directory exists; create if not exist (including parent directories). Throw exception on error.
    void ensureDirectoryExists(const FilePath& p, std::string_view what = "directory");

    // Determine if path is empty (applicable to files)
    bool isEmpty(const FilePath& p);

    // Core: Prepare working directory (create if not exist; require to be directory if exist; optionally require to be empty)
    // - workdir: Directory path to prepare
    // - must_be_empty: If true, directory must be empty when exists, otherwise throw error
    void prepareEmptydir(const FilePath& workdir, bool must_be_empty);

    // Optional: Ensure parent directory of output file exists (call before writing file)
    void ensureParentDirExists(const FilePath& out_file);

    // Determine if given FilePath represents a remote URL (e.g., http://, https://, ftp://, s3://, etc.)
    // Return true means should download via network (remote link); return false means local path.
    // Note: This judgment is heuristic; for file:// or other edge schemes, caller should handle additionally as needed.
    bool isUrl(const FilePath& p);

    // copyFile: Copy file and overwrite destination; prefer std::filesystem::copy_file.
    // If cross-device error (EXDEV) encountered, fallback to stream-based copy to ensure cross-partition copying.
    // Throw exception on error.
    void copyFile(const FilePath& src, const FilePath& dst);

    // downloadFile: Download remote URL to local path dst
    // - If libcurl detected at compile time, use its API (more reliable, controllable)
    // - Otherwise fallback to command line tools (curl/wget), this method depends on runtime environment
    // Throw exception on error and try to delete partial residual files.
    void downloadFile(const std::string& url, const FilePath& dst);

    // fetchFile: If srcOrUrl is URL, call downloadFile, otherwise call copyFile
    // This is unified interface for upper layer calls, shielding local/remote differences
    void fetchFile(const FilePath& srcOrUrl, const FilePath& dst);

    // Recursively delete path (file or directory), throw exception on failure. Use with caution!
    void removeAll(const FilePath& p);

    // Read entire small file into std::string (only applicable when file size is moderate)
    std::string readFileToString(const FilePath& p);

}


namespace seq_io
{
#if defined(FilePath)
    using FilePath = ::FilePath;
#else
    using FilePath = std::filesystem::path;
#endif

    // SeqRecord: Used to represent a sequence record in memory (FASTA/FASTQ)
    // Field semantics:
    //  - id: header name, usually the token after '>' before the first space
    //  - desc: remaining header description information (optional)
    //  - seq: sequence string (A/C/G/T/U/-/N etc.)
    //  - qual: FASTQ quality string (can be filled if FASTQ; leave empty for FASTA)
    struct SeqRecord
    {
        std::string id;     // header name (before space)
        std::string desc;   // remaining header part
        std::string seq;    // sequence
        std::string qual;   // FASTQ quality string (optional), implementation needs to access this field
    };
    using SeqRecords     = std::vector<seq_io::SeqRecord>;

    // ------------------------------------------------------------------
    // SamRecord: Used to represent a SAM alignment record in memory
    // Field semantics (following SAM format standard):
    //  - qname: Query name (QNAME column)
    //  - flag: SAM flag bit flags (FLAG column, e.g., 0=forward unpaired, 16=reverse strand, etc.)
    //  - rname: Reference sequence name (RNAME column, "*" means unmapped)
    //  - pos: 1-based alignment start position (POS column, 0 or "*" means unknown)
    //  - mapq: Mapping quality (MAPQ column, 0-255, 255 means unknown)
    //  - cigar: CIGAR string (CIGAR column, e.g., "100M5I95M", "*" means unknown)
    //  - rnext: Mate/next reference name (RNEXT column, "=" means same, "*" means unknown)
    //  - pnext: Mate 1-based position (PNEXT column)
    //  - tlen: Template length (TLEN column, signed integer, 0 means unknown)
    //  - seq: Query sequence (SEQ column)
    //  - qual: Query quality string (QUAL column, "*" means unknown)
    //  - opt: Optional TAG fields (OPT column, all content starting from column 12)
    //
    // Design notes:
    // - All string fields use std::string (own field ownership, no lifecycle issues)
    // - Suitable for SamReader/SamWriter read/write and long-term storage
    // - Field default values follow SAM specification: "*" means unknown/missing
    // - Performance notes: Using std::string has string copy overhead, but gains simpler lifecycle management
    // ------------------------------------------------------------------
    struct SamRecord
    {
        std::string qname;
        std::uint16_t flag = 0;
        std::string rname = "*";
        std::uint32_t pos = 0;     // 1-based; 0 means '*'
        std::uint8_t mapq = 0;
        std::string cigar = "*";
        std::string rnext = "*";
        std::uint32_t pnext = 0;
        std::int32_t tlen = 0;
        std::string seq = "*";
        std::string qual = "*";
        std::string opt;  // Optional TAG fields (may or may not include leading '\t')
    };

    // makeCleanTable / clean_table:
    // - Used to map any character to a restricted base character set ('A','C','G','T','U','N','-'),
    // - Implemented as a 256-sized lookup table (for unsigned char), avoids runtime branches/conditional judgments, greatly improves batch cleaning performance;
    // - cleanSequence provides in-place cleaning interface, can reduce temporary allocations when reusing on large sequences.
    constexpr std::array<std::uint8_t, 256> makeCleanTable()
    {
        std::array<std::uint8_t, 256> table{};

        for (std::size_t i = 0; i < table.size(); ++i) {
            table[i] = static_cast<std::uint8_t>('N');
        }

        table[static_cast<unsigned char>('A')] = 'A';
        table[static_cast<unsigned char>('a')] = 'A';
        table[static_cast<unsigned char>('C')] = 'C';
        table[static_cast<unsigned char>('c')] = 'C';
        table[static_cast<unsigned char>('G')] = 'G';
        table[static_cast<unsigned char>('g')] = 'G';
        table[static_cast<unsigned char>('T')] = 'T';
        table[static_cast<unsigned char>('t')] = 'T';
        table[static_cast<unsigned char>('U')] = 'U';
        table[static_cast<unsigned char>('u')] = 'U';
        table[static_cast<unsigned char>('-')] = '-';

        return table;
    }

    // "Single definition" constant table in header file: inline constexpr is recommended practice
    inline constexpr auto clean_table = makeCleanTable();
    inline void cleanSequence(std::string& seq)
    {
        for (char& ch : seq) {
            ch = static_cast<char>(clean_table[static_cast<unsigned char>(ch)]);
        }
    }
    inline void cleanSequence(SeqRecord& seq)
    {
        for (char& ch : seq.seq) {
            ch = static_cast<char>(clean_table[static_cast<unsigned char>(ch)]);
        }
    }

    // ISequenceReader: Abstract reader interface, supports next(SeqRecord&) returning false to indicate EOF
    class ISequenceReader
    {
    public:
        virtual ~ISequenceReader() = default;
        virtual bool next(SeqRecord& rec) = 0; // false => EOF
    };

    // KseqReader: Efficient sequence reader based on kseq implementation
    // Note: Implementation located in src/utils/seq_io.cpp, constructor will choose gzopen/gzread or fopen/fread based on whether zlib is enabled
    class KseqReader final : public ISequenceReader
    {
    public:
        explicit KseqReader(const FilePath& file_path);
        ~KseqReader() override;

        KseqReader(const KseqReader&) = delete;
        KseqReader& operator=(const KseqReader&) = delete;

        KseqReader(KseqReader&&) noexcept;
        KseqReader& operator=(KseqReader&&) noexcept;

        bool next(SeqRecord& rec) override;

    private:
        struct Impl;
        std::unique_ptr<Impl> impl_;
    };

    // SeqWriter: Universal sequence/alignment result writer
    //
    // Supported formats:
    // - FASTA: Write SeqRecord (id/desc/seq)
    // - SAM  : Write minimal viable SAM (header + alignment records)
    //
    // Note:
    // - Uses 8MiB internal aggregation buffer by default to minimize disk write calls as much as possible;
    // - For compatibility with old usage, retain write(const SeqRecord&) as default entry for FASTA writing.
    class SeqWriter
    {
    public:
        enum class Format : std::uint8_t
        {
            fasta = 0,
            sam   = 1,
        };

        // FASTA (default)
        explicit SeqWriter(const FilePath& file_path, std::size_t line_width = 80);
        explicit SeqWriter(const FilePath& file_path, Format fmt, std::size_t line_width, std::size_t buffer_threshold_bytes);

        SeqWriter(const FilePath& file_path, std::size_t line_width, std::size_t buffer_threshold_bytes);

        SeqWriter(const SeqWriter&) = delete;
        SeqWriter& operator=(const SeqWriter&) = delete;

        SeqWriter(SeqWriter&&) noexcept = default;
        SeqWriter& operator=(SeqWriter&&) noexcept = default;


        // SAM constructor (factory function)
        static SeqWriter Sam(const FilePath& file_path, std::size_t buffer_threshold_bytes = 8ULL * 1024ULL * 1024ULL);

        // Current mode
        Format format() const noexcept { return format_; }

        // ------------------------- FASTA -------------------------
        void writeFasta(const SeqRecord& rec);

        // Compatible: default write() writes FASTA
        void write(const SeqRecord& rec) { writeFasta(rec); }

        // ------------------------- SAM -------------------------
        void writeSamHeader(std::string_view header_text);

        void writeSam(const SamRecord& r);

        // flush():
        // - First flush internal aggregation buffer, then flush ofstream
        void flush();

        ~SeqWriter();

    private:

        std::ofstream out_;
        Format format_{Format::fasta};
        std::size_t line_width_{80};

        std::string buffer_;
        std::size_t buffer_threshold_bytes_{8ULL * 1024ULL * 1024ULL};

        bool sam_header_written_{false};

        void flushBuffer_();
        void appendOrFlush_(std::string_view s);
        static void appendWrapped_(std::string& dst, std::string_view s, std::size_t width);
    };

    inline std::unique_ptr<ISequenceReader> openKseqReader(const FilePath& file_path)
    {
        return std::make_unique<KseqReader>(file_path);
    }

    // ------------------------------------------------------------------
    // Class: SamReader
    // Function: Read SAM file and parse into SamRecord
    //
    // Description:
    // 1. SAM format is tab-separated text file, each line contains one alignment record
    // 2. SAM file may contain header lines starting with '@', this class will automatically skip
    // 3. Parse all required SAM fields when reading (QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, SEQ, QUAL, etc.)
    // 4. Support efficient reading of large files (optimized with large buffer)
    //
    // Performance notes:
    // - Use std::getline for line-by-line reading, with large buffer to reduce system calls
    // - Use string_view for field splitting, avoid unnecessary string copies
    // - Default buffer size 8MiB, suitable for large-scale SAM files
    // - Internally save field strings to ensure string_view lifecycle safety
    // ------------------------------------------------------------------
    class SamReader
    {
    public:
        // Constructor: Open SAM file and initialize reader
        // Parameters:
        //   - file_path: SAM file path
        //   - buffer_size: Input buffer size (default 8MiB)
        explicit SamReader(const FilePath& file_path, std::size_t buffer_size = 8ULL * 1024ULL * 1024ULL);

        ~SamReader();

        // Disable copy, allow move
        SamReader(const SamReader&) = delete;
        SamReader& operator=(const SamReader&) = delete;
        SamReader(SamReader&&) noexcept;
        SamReader& operator=(SamReader&&) noexcept;


        bool next(SamRecord& rec);

    private:
        struct Impl;
        std::unique_ptr<Impl> impl_;
    };

    // ------------------------------------------------------------------
    // Function: convertSamToFasta
    // Function: Read SAM file and convert to FASTA file
    //
    // Description:
    // 1. Extract QNAME and SEQ fields from SAM file, write in FASTA format
    // 2. Automatically skip SAM header lines (starting with '@')
    // 3. Preserve query sequence name and sequence content, discard alignment information and quality values
    // 4. Support efficient conversion of large files (optimized with buffer)
    //
    // Parameters:
    //   - sam_path: Input SAM file path
    //   - fasta_path: Output FASTA file path
    //   - line_width: FASTA line width (default 80, 0 means no line wrapping)
    //
    // Performance:
    // - Use 8MiB input/output buffer to reduce I/O system calls
    // - Streaming processing, memory usage independent of file size
    //
    // Exceptions: Throw std::runtime_error when file opening fails or parsing error
    // ------------------------------------------------------------------
    void convertSamToFasta(const FilePath& sam_path, const FilePath& fasta_path, std::size_t line_width = 80);


    // ------------------------------------------------------------------
    // Function: makeSamRecord
    // Function: Build SAM record from SeqRecord and alignment information
    //
    // Parameters:
    //   - query: SeqRecord of query sequence
    //   - ref_name: Reference sequence name
    //   - cigar_str: CIGAR string (e.g., "100M5I95M")
    //   - pos: 1-based alignment start position (default 1)
    //   - mapq: mapping quality (default 60)
    //   - flag: SAM flag (default 0)
    //
    // Return: SamRecord (independent type under seq_io namespace)
    //
    // Description:
    //   1. This is a convenience function that encapsulates the conversion logic from SeqRecord to SamRecord
    //   2. Quality value handling: If query.qual is empty, use SAM default value "*"
    //   3. Lifecycle of all string_view fields is guaranteed by caller (query and cigar_str must be valid during use)
    //   4. Simplified SAM record, can be extended to add more fields later (AS, NM, etc.)
    //
    // Performance: O(1), only assignment operations, no memory allocation
    // ------------------------------------------------------------------
    SamRecord makeSamRecord(
        const SeqRecord& query,
        std::string_view ref_name,
        std::string_view cigar_str,
        std::uint32_t pos = 1,
        std::uint8_t mapq = 60,
        std::uint16_t flag = 0
    );

    // ------------------------------------------------------------------
    // Function: samRecordToSeqRecord
    // Function: Convert sequence information (QNAME/SEQ/QUAL) in SamRecord to SeqRecord
    //
    // Why need it:
    // - SamReader::next() returns SamRecord (many fields, contains all SAM columns)
    // - But many subsequent processing (e.g., write FASTA, cleanSequence, interface with existing code) are more accustomed to using SeqRecord
    // - This function centralizes the repetitive logic of "extracting sequence from SamRecord" in one place, avoiding copy code everywhere
    //
    // Input:
    // - sam_rec: SamRecord (usually from SamReader::next), where qname/seq/qual are std::string
    // - keep_qual: Whether to copy qual to SeqRecord (FASTA scenarios usually don't need qual, suggest false)
    //
    // Output:
    // - Return a SeqRecord:
    //   - id   = sam_rec.qname
    //   - desc = "" (SAM doesn't have FASTA header's desc concept)
    //   - seq  = sam_rec.seq
    //   - qual = (keep_qual ? sam_rec.qual : "")
    //
    // Performance and correctness:
    // - Since SamRecord already owns all strings, string copy will occur here
    // - Complexity O(|qname| + |seq| + |qual|), for large files IO is usually the bottleneck
    // ------------------------------------------------------------------
    // Output:
    // - Return a SeqRecord:
    //   - id   = sam_rec.qname
    //   - desc = "" (SAM doesn't have FASTA header's desc concept)
    //   - seq  = sam_rec.seq
    //   - qual = (keep_qual ? sam_rec.qual : "")
    //
    // Performance and correctness:
    // - Since SamRecord already owns all strings, string copy will occur here
    // - Complexity O(|qname| + |seq| + |qual|), for large files IO is usually the bottleneck
    // ------------------------------------------------------------------
    inline SeqRecord samRecordToSeqRecord(const SamRecord& sam_rec, bool keep_qual = false)
    {
        SeqRecord rec;
        rec.id = sam_rec.qname;
        rec.desc.clear();
        rec.seq = sam_rec.seq;

        if (keep_qual) {
            // SAM's qual may be "*", caller can filter again if needed.
            rec.qual = sam_rec.qual;
        } else {
            rec.qual.clear();
        }

        return rec;
    }

} // namespace seq_io

namespace cmd
{
    struct BuildOptions
    {
        bool quiet = true;                  // Whether to append silent redirection
        bool close_stdin = true;            // Whether to append "< /dev/null"
        bool detect_stdout_redirect = true; // Whether to decide silent strategy based on '>'
    };

    // cmd_template must contain {input} and {output}; optionally contain {thread}
    // Linux-only rules:
    // - If '>' detected in template (stdout already redirected): append "2>/dev/null"
    // - Otherwise append "> /dev/null 2>&1"
    // - Optionally append "< /dev/null" to close stdin
    std::string buildCommand(std::string cmd_template,
                             const std::string& input_path,
                             const std::string& output_path,
                             int thread,
                             const BuildOptions& opt = BuildOptions{});

    // Execute command and return real exit code:
    // - Normal exit: return exit code
    // - Signal termination: return 128 + signal
    // - system() itself failure: return -1
    int runCommand(const std::string& command);

    // Self-check with tiny.fasta, run template command and check if output file is generated
    bool testCommandTemplate(const std::string& cmd_template, const FilePath& workdir, int thread = 1);

} // namespace cmd


// ================================================================
// ProgressBar - Progress bar utility class (single line refresh, no need to know total in advance)
//
// Design goals:
// 1) Unify implementation of all progress bars in the project, avoid code redundancy
// 2) Support scenarios where total is not known in advance (only show processed count + rate)
// 3) Thread safety consideration: This class itself does not lock, suitable for serial area calls
//    If used in parallel areas, caller should ensure serial access (e.g., OpenMP serial area)
//
// Usage example:
//   ProgressBar progress("align");
//   for (...) {
//       // ... processing work ...
//       progress.tick();  // process one
//   }
//   progress.done();  // force output and line break
// ================================================================
class ProgressBar {
public:
    // Constructor
    // - label: Progress bar prefix label (e.g., "preprocess", "align", "merge")
    // - report_interval: How many items to process before refreshing (default 1000)
    explicit ProgressBar(std::string label, std::size_t report_interval = 1000) noexcept
        : label_(std::move(label))
        , report_interval_(report_interval)
        , count_(0)
        , next_report_(report_interval)
        , start_time_(std::chrono::steady_clock::now())
    {}

    // Disable copy (avoid accidental copying causing state confusion)
    ProgressBar(const ProgressBar&) = delete;
    ProgressBar& operator=(const ProgressBar&) = delete;

    // Allow move
    ProgressBar(ProgressBar&&) noexcept = default;
    ProgressBar& operator=(ProgressBar&&) noexcept = default;

    // tick: Call after processing one record, internally automatically determines if refresh is needed
    void tick() noexcept {
        ++count_;
        if (count_ >= next_report_) {
            print(false);
            next_report_ = ((count_ / report_interval_) + 1) * report_interval_;
        }
    }

    // tick: Call after processing multiple records (batch update, reduce call overhead)
    void tick(std::size_t n) noexcept {
        count_ += n;
        if (count_ >= next_report_) {
            print(false);
            next_report_ = ((count_ / report_interval_) + 1) * report_interval_;
        }
    }

    // done: Force output final status and line break (usually called after loop ends)
    void done() noexcept {
        print(true);
        std::fprintf(stderr, "\n");
    }

    // Get current processed count
    [[nodiscard]] std::size_t count() const noexcept { return count_; }

    // Reset progress bar (for reusing same instance to handle multiple stages)
    void reset() noexcept {
        count_ = 0;
        next_report_ = report_interval_;
        start_time_ = std::chrono::steady_clock::now();
    }

    // Reset progress bar and change label
    void reset(std::string new_label) noexcept {
        label_ = std::move(new_label);
        reset();
    }

private:
    // print: Output progress information to stderr (use \r to return to line start and overwrite)
    // - force: If true then force output, otherwise only output when threshold reached
    void print(bool force) const noexcept {
        if (!force && count_ < next_report_) {
            return;
        }

        const auto now = std::chrono::steady_clock::now();
        const double elapsed_sec = std::chrono::duration_cast<std::chrono::duration<double>>(
            now - start_time_).count();
        const double rate = (elapsed_sec > 1e-6)
            ? (static_cast<double>(count_) / elapsed_sec)
            : 0.0;

        // ANSI color code: green = \033[32m, reset = \033[0m
        // Trailing space used to overwrite possible remaining old characters
        std::fprintf(stderr,
                     "\r\033[32m[%s] processed=%zu  elapsed=%.1fs  rate=%.0f seq/s\033[0m   ",
                     label_.c_str(), count_, elapsed_sec, rate);
        std::fflush(stderr);
    }

    std::string label_;                                  // Progress bar prefix label
    std::size_t report_interval_;                        // Refresh interval (how many items to process before refreshing)
    std::size_t count_;                                  // Number of processed records
    std::size_t next_report_;                            // Next refresh threshold
    std::chrono::steady_clock::time_point start_time_;   // Start time
};

#endif //HALIGN4_UTILS_H

