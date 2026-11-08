#include <config.hpp>
#include <utils.h>
#include "preprocess.h"
#include "consensus.h"

#include "align.h"
#include "align/strategy.hpp"
#include "align/native_msa.hpp"

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
    if (msa_cmd_str != "native") {
        if (cmd::testCommandTemplate(msa_cmd_str, opt.workdir, opt.threads)) {
            spdlog::info("msa_cmd template test passed.");
        } else {
            throw std::runtime_error("msa_cmd template test failed.");
        }
    } else {
        spdlog::info("Using in-process NativeMSA; no external command self-test needed.");
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

        // PR3C.2 is deliberately opt-in: the historic --mode auto and
        // --mode direct routing remain untouched. Conflicting flags must be
        // rejected before the scanner, workdir creation or external execution.
        if (opt.auto_hierarchical) {
            if (opt.mode != "auto") {
                throw std::invalid_argument("--auto-hierarchical requires --mode auto");
            }
            if (opt.keep_length || !opt.center_path.empty() ||
                opt.anchored_insertions || opt.native_center ||
                opt.soft_anchors || opt.soft_anchor_recruit) {
                throw std::invalid_argument(
                    "--auto-hierarchical manages backend options automatically: "
                    "do not pass --keep-length, --center-path, --native-center, "
                    "--anchored-insertions, --soft-anchors or --soft-anchor-recruit");
            }
        }

        // The experimental hierarchical center is a separate opt-in. Never
        // allow the direct/auto router to silently reinterpret its parameters.
        if (opt.native_center && opt.mode != "legacy") {
            throw std::invalid_argument("--native-center requires --mode legacy");
        }

        // PR3A: make the routing decision *before* preprocessing (which keeps
        // the Top-K full sequences in RAM). Explicit -p remains available in
        // the legacy mode; auto/direct cannot silently override it.
        if (opt.mode != "legacy") {
            if (!opt.msa_cmd.empty() && resolveMsaCmdTemplate(opt.msa_cmd) != "native") {
                throw std::invalid_argument("--mode auto/direct conflicts with -p/--msa-cmd; "
                                            "omit -p or use -p native");
            }
            const auto stats = align_strategy::scanFasta(opt.input);
            const native_msa::Options native_limits;
            const align_strategy::Budget budget{
                native_limits.max_sequences,
                native_limits.max_sequence_bases,
                native_limits.max_total_bases
            };
            spdlog::info("Alignment input: {} sequences, {} bases, min={} max={}",
                         stats.sequence_count, stats.total_bases,
                         stats.shortest_bases, stats.longest_bases);
            if (opt.auto_hierarchical) {
                if (opt.cons_n <= 0) {
                    throw std::invalid_argument("--auto-hierarchical requires positive --cons-n");
                }
                const auto route = align_strategy::chooseBoundedAuto(
                    stats, budget, static_cast<std::uint64_t>(opt.cons_n));
                if (route.backend == align_strategy::AutoBackend::Reject) {
                    throw std::runtime_error("--auto-hierarchical rejected input: " +
                                             route.explanation);
                }
                if (route.backend == align_strategy::AutoBackend::NativeHierarchical) {
                    // Reuse the fully validated PR3C.1 explicit backend, not
                    // a second implementation. Switch only after routing.
                    // This is the internal selection of the legacy *pipeline*,
                    // never the legacy external-MSA backend.
                    opt.mode = "legacy";
                    opt.native_center = true;
                    opt.anchored_insertions = true;
                    opt.soft_anchors = true;
                    opt.soft_anchor_recruit = true;
                    spdlog::info("Selected bounded auto hierarchical: {}", route.explanation);
                } else {
                    spdlog::info("Selected bounded auto direct: {}", route.explanation);
                }
                opt.msa_cmd = "native";
            } else {
                // Preserve the original --mode auto/direct behavior exactly.
                const auto choice = align_strategy::chooseNativeDirect(
                    stats, budget, static_cast<std::uint64_t>(opt.cons_n),
                    opt.keep_length, !opt.center_path.empty());
                if (!choice.native_direct) {
                    throw std::runtime_error("--mode " + opt.mode + ": " + choice.explanation +
                        ". PR3A intentionally does not auto-route large input to the "
                        "legacy Python chunker. Use --mode legacy explicitly, or wait "
                        "for native hierarchical/ultralong routing.");
                }
                spdlog::info("Selected bounded NativeMSA: {}", choice.explanation);
                opt.msa_cmd = "native";
            }
        }

        // PR3C.1: opt-in, bounded in-process Top-K center. The same native
        // backend is already used by direct MSA and local insertion MSAs.
        // Crucially this is NOT the native direct routing: Top-K must remain
        // longest-first, and all queries still pass through the anchored
        // hierarchical alignment below.
        if (opt.native_center) {
            if (!opt.anchored_insertions || opt.keep_length ||
                !opt.center_path.empty() || opt.cons_n < 2 || opt.cons_n > 64) {
                throw std::invalid_argument(
                    "--native-center requires --anchored-insertions, --keep-length off, "
                    "no --center-path, and --cons-n between 2 and 64");
            }
            // An explicit external -p is contradictory, rather than silently
            // ignored. -p native is accepted as redundant for convenience.
            if (!opt.msa_cmd.empty() && resolveMsaCmdTemplate(opt.msa_cmd) != "native") {
                throw std::invalid_argument(
                    "--native-center conflicts with external -p/--msa-cmd; omit -p");
            }
            opt.msa_cmd = "native";
            spdlog::info("Experimental native hierarchical center enabled (Top-{} longest)", opt.cons_n);
        }

        if (opt.soft_anchors && !opt.anchored_insertions) {
            throw std::invalid_argument("--soft-anchors requires --anchored-insertions");
        }
        if (opt.soft_anchor_recruit && !opt.soft_anchors) {
            throw std::invalid_argument("--soft-anchor-recruit requires --soft-anchors");
        }

        // PR3B.5A is intentionally opt-in and bounded; the default legacy
        // pipeline and PR3A direct/auto routing remain unchanged.
        if (opt.anchored_insertions &&
            (opt.mode != "legacy" || opt.keep_length ||
             (resolveMsaCmdTemplate(opt.msa_cmd) == "native" && !opt.native_center))) {
            throw std::invalid_argument("--anchored-insertions requires --mode legacy, "
                                        "--keep-length off; native center requires --native-center");
        }

        if (opt.anchored_insertions) {
            // Reject oversized inputs BEFORE preprocessing or alignment buffers
            // materialize a full batch. This is experimental short-DNA only.
            const auto stats = align_strategy::scanFasta(opt.input);
            if (stats.sequence_count < 2 || stats.sequence_count > 64 ||
                stats.longest_bases > 25000 || stats.total_bases > 100000 ||
                stats.contains_gaps) {
                throw std::invalid_argument(
                    "--anchored-insertions supports 2..64 ungapped DNA queries, "
                    "at most 25kb/sequence and 100kb input total");
            }
            if (!opt.center_path.empty()) {
                const auto center = align_strategy::scanFasta(opt.center_path);
                if (center.sequence_count != 1 || center.longest_bases > 25000 ||
                    center.contains_gaps) {
                    throw std::invalid_argument(
                        "--anchored-insertions requires exactly one ungapped center <=25kb");
                }
            }
        }

        // Check parameters
        checkOption(opt);
        setupLoggerWithFile(opt.workdir);

        // preprocessing
        const uint_t preproc_count = preprocessInputFasta(
            opt.input, opt.workdir, opt.cons_n,
            opt.msa_cmd == "native" && !opt.native_center);
        spdlog::info("Preprocessing produced {} records", preproc_count);

        if (opt.anchored_insertions && preproc_count <= static_cast<uint_t>(opt.cons_n)) {
            throw std::invalid_argument("--anchored-insertions requires hierarchical input count > --cons-n");
        }

        // PR 2B: deliberately expose only safe, small direct MSA. Do not
        // silently enter the hierarchical path with size-limited NativeMSA.
        // PR 3 will add budget-based automatic backend selection.
        if (opt.msa_cmd == "native" && !opt.native_center &&
            (preproc_count > static_cast<uint_t>(opt.cons_n) || preproc_count > 64 ||
             opt.keep_length || !opt.center_path.empty())) {
            throw std::runtime_error(
                "-p native currently supports only direct MSA of <=64 sequences "
                "with --keep-length disabled and without --center-path; "
                "use the legacy backend for hierarchical/large inputs");
        }

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
