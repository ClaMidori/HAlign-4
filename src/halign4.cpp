#include <config.hpp>
#include <utils.h>
#include "preprocess.h"
#include "consensus.h"

#include "align.h"

// Program entry: command line parsing -> preprocessing -> consensus alignment -> sequence alignment -> result merging -> cleanup working directory

// Parameter validation and working directory preparation
static void checkOption(Options& opt) {
    // File verification
    file_io::requireRegularFile(opt.input, "input");
    if (!opt.center_path.empty()) {
        file_io::requireRegularFile(opt.center_path, "center_path");
    }

    // Numerical verification
    if (opt.threads <= 0) throw std::runtime_error("threads must be > 0");
    if (opt.kmer_size <= 0) throw std::runtime_error("kmer_size must be > 0");
    if (opt.kmer_window <= 0) throw std::runtime_error("kmer_window must be > 0");
    if (opt.cons_n <= 0) throw std::runtime_error("cons_n must be > 0");
    if (opt.kmer_size > 31) throw std::runtime_error("kmer_size too large (must be <= 31)");
    if (opt.kmer_window >= 256) {
        spdlog::warn("kmer_window >= 256 may be slow; current value: {}", opt.kmer_window);
    }

    // workdir preparation
#ifdef _DEBUG
    constexpr bool must_be_empty = false;
#else
    constexpr bool must_be_empty = true;
#endif
    file_io::prepareEmptydir(opt.workdir, must_be_empty);

    // MSA Command Template Parsing and Self-Test
    const std::string msa_cmd_str = resolveMsaCmdTemplate(opt.msa_cmd);
    if (cmd::testCommandTemplate(msa_cmd_str, opt.workdir, opt.threads)) {
        spdlog::info("msa_cmd template test passed.");
    } else {
        throw std::runtime_error("msa_cmd template test failed.");
    }
    opt.msa_cmd = msa_cmd_str;
}

// Clean working directory
static void cleanupWorkdir(const Options& opt) {
    if (!opt.save_workdir) {
        try {
            spdlog::info("Removing working directory: {}", opt.workdir);
            file_io::removeAll(FilePath(opt.workdir));
            spdlog::info("Working directory removed successfully");
        } catch (const std::exception& e) {
            spdlog::warn("Failed to remove working directory: {}", e.what());
        }
    } else {
        spdlog::info("Keeping working directory: {}", opt.workdir);
    }
}

int main(int argc, char** argv) {
    try
    {
        // Initialization logic: set up thread pool for spdlog, then parse CLI options, set up logger, and log the parsed options.
        spdlog::init_thread_pool(8192, 1);
        setupLogger();

        Options opt;
        CLI::App app{"halign4"};
        setupCli(app, opt);
        CLI11_PARSE(app, argc, argv);

        // Set default working directory if not provided, and log the final options for user confirmation.
        if (opt.workdir.empty()) {
            opt.workdir = makeDefaultWorkdir();
            spdlog::info("--workdir not provided, using default: {}", opt.workdir);
        }

        // Print parameters
        logParsedOptions(opt);
        spdlog::info("Starting halign4 version {}...", VERSION);

        // Check parameters
        checkOption(opt);
        setupLoggerWithFile(opt.workdir);

        // preprocessing
        const uint_t preproc_count = preprocessInputFasta(opt.input, opt.workdir, opt.cons_n);
        spdlog::info("Preprocessing produced {} records", preproc_count);

        // File path definitions
        const FilePath consensus_unaligned_file = FilePath(opt.workdir) / WORKDIR_DATA / DATA_CLEAN / CLEAN_CONS_UNALIGNED;
        const FilePath consensus_aligned_file = FilePath(opt.workdir) / WORKDIR_DATA / DATA_CLEAN / CLEAN_CONS_ALIGNED;
        const FilePath consensus_file = FilePath(opt.workdir) / WORKDIR_DATA / DATA_CLEAN / CLEAN_CONS_FASTA;
        const FilePath consensus_json_file = FilePath(opt.workdir) / WORKDIR_DATA / DATA_CLEAN / CLEAN_CONS_JSON;

        // Process consensus sequence: if user provided center_path, copy it to consensus_unaligned_file; otherwise generate consensus sequence from preprocessed data. If the number of sequences is small and length is not kept, skip consensus generation and directly align to get final output.
        if (!opt.center_path.empty())
        {
            spdlog::info("Using user-specified center sequence: {}", opt.center_path);
            if (std::filesystem::exists(consensus_unaligned_file)) {
                file_io::removeAll(consensus_unaligned_file);
            }
            file_io::copyFile(FilePath(opt.center_path), consensus_unaligned_file);
            spdlog::info("Center sequence copied to: {}", consensus_unaligned_file.string());
        }

        // Fast path: Output directly when the number of sequences is less than or equal to cons_n and length is not retained.
        if (preproc_count <= opt.cons_n && opt.keep_length == false)
        {
            alignConsensusSequence(consensus_unaligned_file, consensus_aligned_file, opt.msa_cmd, opt.threads);
            file_io::copyFile(consensus_aligned_file, FilePath(opt.output));
            spdlog::info("All sequences processed; final output written to {}", opt.output);

            cleanupWorkdir(opt);
            spdlog::info("halign4 End!");
            return 0;
        }
        else if (opt.center_path.empty())
        {
            // Generate consensus sequence
            alignConsensusSequence(consensus_unaligned_file, consensus_aligned_file, opt.msa_cmd, opt.threads);

            const std::string consensus_string = consensus::generateConsensusSequence(
                consensus_aligned_file,
                consensus_file,
                consensus_json_file,
                opt.cons_n,
                opt.threads,
                4096
            );

            spdlog::info("Consensus sequence generated with length {}", consensus_string.size());
        }

        // comparison stage: align sequences to consensus and merge results
        const FilePath ref_path = opt.center_path.empty() ? consensus_file : FilePath(opt.center_path);
        align::RefAligner ref_aligner(opt, ref_path);
        ref_aligner.alignQueryToRef(opt.input);
        ref_aligner.mergeAlignedResults(opt.output, 25600);

        cleanupWorkdir(opt);

        spdlog::info("halign4 End!");
        return 0;
    } catch (const std::exception &e) {
        spdlog::error("Fatal error: {}", e.what());
        spdlog::error("halign4 End!");
        return 1;
    } catch (...) {
        spdlog::error("Fatal error: unknown exception");
        spdlog::error("halign4 End!");
        return 1;
    }
}
