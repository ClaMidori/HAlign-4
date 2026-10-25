#include "preprocess.h"
#include <chrono>
// Add a header file for calling external commands, constructing strings, and checking file existence.
#include <cstdlib>
#include <sstream>
#include <filesystem>
#include <stdexcept>
#include "msa_validation.hpp"

// This document contains the preprocessing logic for the FASTA input:
// - The input file (which can be a local path or a URL) is retrieved to the data/raw directory in the working directory;
// - Perform simple cleaning on the sequence (capitalize, replace non-AGCTU characters with N);
// - Write the cleaned sequences to the data/clean directory;
// - Maintain a Top-K selector to choose the longest K sequences (stably preserving earlier occurring sequences);
// - Finally, output the cleaned data and selected consensus to the specified file (this file only handles selection and writing interfaces, with actual writing handled by the caller).

// Explanation: The FilePath used here is file_io::FilePath (i.e., an alias for std::filesystem::path), which provides convenient path manipulation and is used throughout the codebase for file handling. The file_io namespace contains utility functions for file operations, such as checking if a path is a URL, ensuring directories exist, and handling file I/O errors. The preprocessInputFasta function will utilize these utilities to manage input files and directories effectively while performing the necessary preprocessing steps on the FASTA sequences.
// The seq_io namespace encapsulates the details of reading/writing FASTA (openKseqReader / SeqWriter / SeqRecord, etc.).

uint_t preprocessInputFasta(const std::string input_path, const std::string workdir, const int cons_n) {
    // Parameter description:
    // - input_path: Enter the path or URL (string) of FASTA.
    // - workdir: The working directory path (which should already be prepared or guaranteed by the caller). Subdirectories such as data/raw and data/clean will be created under this directory.
    // - cons_n: The number of longest sequences to retain (Top-K).

    // Timing is used for logging and helps with performance analysis.
    const auto t_start = std::chrono::steady_clock::now();

    spdlog::info("Preprocessing input FASTA file: {}", input_path);
    spdlog::info("Working directory: {}", workdir);

    // ---------- Catalog preparation ----------
    // 1) Ensure a "data" directory exists under the working directory.
    //    This directory stores both raw and cleaned data: data/raw and data/clean.
    FilePath data_dir = FilePath(workdir) / WORKDIR_DATA;
    file_io::ensureDirectoryExists(data_dir);
    spdlog::info("Ensured data directory exists: {}", data_dir.string());

    // 2) Create the raw_data directory under data for storing the original (uncleaned) input.
    FilePath raw_data_dir = data_dir / DATA_RAW;
    file_io::ensureDirectoryExists(raw_data_dir);
    spdlog::info("Ensured raw data directory exists: {}", raw_data_dir.string());

    // 3) Create the clean_data directory under data for storing cleaned output.
    FilePath clean_data_dir = data_dir / DATA_CLEAN;
    file_io::ensureDirectoryExists(clean_data_dir);
    spdlog::info("Ensured clean data directory exists: {}", clean_data_dir.string());

    // ========================================================================
    // Fetch input file (supports local/remote with performance optimization)
    // ========================================================================
    // Key optimizations:
    // - Remote file (URL): download into raw_data directory (requires local cache)
    // - Local file: read directly from original path (avoids unnecessary copy overhead)
    //
    // Performance impact:
    // - Large local files (e.g., 10GB FASTA): reduces minutes of copy time to 0 seconds
    // - Remote files: retains existing logic (must be downloaded locally)
    // - Memory usage: unaffected (streamed reads in all cases)
    // ========================================================================
    FilePath input_file = FilePath(input_path);
    FilePath actual_input_file;  // actual file path used for reading

    if (file_io::isUrl(input_file)) {
        // Remote file: download into raw_data directory
        FilePath raw_dest_file = raw_data_dir / input_file.filename();

        spdlog::info("Detected remote URL, downloading to: {} -> {}",
                     input_file.string(), raw_dest_file.string());

        // fetchFile internally calls downloadFile
        file_io::fetchFile(input_file, raw_dest_file);

        spdlog::info("Download completed: {}", raw_dest_file.string());
        actual_input_file = raw_dest_file;
    } else {
        // Local file: use original path directly without copying
        spdlog::info("Detected local file, reading directly from: {}", input_file.string());

        // Verify file exists and is readable (throws on failure)
        file_io::requireRegularFile(input_file, "input file");

        actual_input_file = input_file;

        spdlog::info("Local file verified, no copy needed (performance optimization)");
    }

    // ---------- Prepare output file names ----------
    // 5) Open the raw input and read records one-by-one; clean each sequence (cleanSequence) and write to clean_data.
    //    At the same time maintain a TopKLongestSelector to pick the longest cons_n sequences (for later consensus generation).
    // handle input filenames like `sample.fasta.gz` -> `sample.fasta`
    FilePath in_fname = input_file.filename();
    std::string in_name = in_fname.string();
    const std::string comp_ext = ".gz";
    if (in_name.size() > comp_ext.size() &&
        in_name.compare(in_name.size() - comp_ext.size(), comp_ext.size(), comp_ext) == 0) {
        in_name.resize(in_name.size() - comp_ext.size());
        spdlog::info("Detected compressed input; using output name: {}", in_name);
    }
    FilePath clean_dest_file = clean_data_dir / FilePath(in_name);
    FilePath consensus_file = clean_data_dir / CLEAN_CONS_UNALIGNED;
    spdlog::info("Clean output: {} ; Consensus output: {}", clean_dest_file.string(), consensus_file.string());

    // ---------- Open reader/writer and TopK selector ----------
    // seq_io::openKseqReader returns an abstract reader pointer for reading sequences one-by-one;
    // seq_io::SeqWriter is used to write cleaned sequences to the target file.
    //
    // Note: actual_input_file may be:
    // 1. A remote file: the downloaded file under raw_data
    // 2. A local file: the original path provided by the user (no copy needed, performance optimized)
    auto reader = seq_io::openKseqReader(actual_input_file);
    seq_io::SeqWriter clean_writer(clean_dest_file);
    TopKLongestSelector selector(cons_n);

    // Processing loop: read -> clean -> write -> submit to TopK selector
    seq_io::SeqRecord rec;
    std::size_t total_records = 0;

    // Use ProgressBar class instead of manual progress printing to reduce boilerplate
    ProgressBar progress("preprocess");

    // Important notes (performance and correctness):
    // - This loop is the hot path for preprocessing. For large inputs (tens of thousands / millions of records), watch IO and memory usage.
    // - Performance optimizations:
    //    * Use seq_io's KseqReader (based on fread/gzread) with a large buffer to significantly improve read throughput.
    //    * SeqWriter::write buffers a record's header and wrapped sequence into a temporary string and writes it in one shot,
    //      avoiding per-character writes and the associated system-call overhead; this is important for large files.
    //    * TopKLongestSelector should be implemented as a min-heap of size K, with insert/replace cost O(log K), which is suitable when K is much smaller than total record count.
    // - Memory trade-off: the TopK implementation keeps K full records (using O(K * avg_len) memory); if K is large, watch memory usage.
    while (reader->next(rec)) {
        ++total_records;
        // Normalize and clean the sequence: e.g., uppercase letters and replace non-AGCTU characters with N (implemented by seq_io::cleanSequence).
        // cleanSequence modifies rec.seq in place, avoiding repeated copies to save memory bandwidth.
        seq_io::cleanSequence(rec.seq);

        // Write the cleaned record to the clean_data directory
        // Note: FastaWriter::write already optimizes by concatenating and writing in one shot, which is performance-friendly.
        clean_writer.write(rec);

        // Submit the current record to the TopK selector (internally maintains a heap for O(log K) replacement cost)
        // Note: selector.consider should copy or take ownership of needed fields (e.g., id/seq) to avoid data corruption when rec is reused/overwritten.
        selector.consider(rec);

        // Update progress via ProgressBar::tick(); it decides internally whether to refresh.
        progress.tick();
    }
    progress.done();  // force final status output and newline

    // ---------- Write TopK results to consensus input file ----------
    // Note: takeSortedDesc returns records sorted by length in descending order (typically used to select the longest N sequences for consensus generation)
    seq_io::SeqWriter cons_writer(consensus_file);
    auto cons_seqs = selector.takeSortedDesc();

    // Write the selected sequences to the consensus file; keep consistent line-wrapping rules, which SeqWriter handles.
    for (const auto& cons_rec : cons_seqs) {
        // The cons_rec written here should be a deep-copied SeqRecord (returned by selector to ensure safety); if not, make a copy in selector.
        cons_writer.write(cons_rec);
    }

    // ---------- Stats and return value ----------
    const auto t_end = std::chrono::steady_clock::now();
    const double elapsed_s = std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start).count();

    spdlog::info("Preprocessing completed. Total records processed: {}. Selected top {} sequences: {}. Elapsed: {:.2f} s",
                 total_records, cons_n, cons_seqs.size(), elapsed_s);
    // Convert size_t to the project-level uint_t (defined in config.hpp); truncate to U_MAX to prevent overflow
    if (total_records > static_cast<std::size_t>(U_MAX)) {
        spdlog::warn("Processed records ({}) exceed U_MAX ({}); truncating to U_MAX", total_records, U_MAX);
        return static_cast<uint_t>(U_MAX);
    }

    return static_cast<uint_t>(total_records);
}


void alignConsensusSequence(const FilePath& input_file, const FilePath& output_file,
                            const std::string& msa_cmd, int threads)
{
    // Do not permit a failed/stale external command to supply a previous result.
    file_io::requireRegularFile(input_file, "unaligned consensus FASTA");
    if (std::filesystem::exists(output_file) &&
        std::filesystem::equivalent(input_file, output_file)) {
        throw std::runtime_error("input and output FASTA must be distinct files");
    }
    file_io::ensureParentDirExists(output_file);
    std::filesystem::remove(output_file);  // Clear stale output before launching.

    const auto start = std::chrono::steady_clock::now();
    const auto command = cmd::buildCommand(
        msa_cmd, input_file.string(), output_file.string(), threads, cmd::BuildOptions{});
    spdlog::info("Starting MSA command: {}", command);
    try {
        const int rc = cmd::runCommand(command);
        if (rc != 0) {
            throw std::runtime_error("MSA command failed (exit code " +
                                     std::to_string(rc) + "): " + command);
        }

        const auto stats = msa::validateFastaAlignment(input_file, output_file);
        const auto elapsed = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start).count();
        spdlog::info("MSA verified: {} sequences, {} columns, {} bases; elapsed {:.3f}s",
                     stats.records, stats.columns, stats.bases, elapsed);
    } catch (...) {
        // Remove partial or invalid output. Let main() propagate a nonzero exit status.
        std::error_code ignored;
        std::filesystem::remove(output_file, ignored);
        throw;
    }
}
