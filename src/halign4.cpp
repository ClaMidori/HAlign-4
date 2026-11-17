#include <config.hpp>
#include <utils.h>
#include "preprocess.h"
#include "consensus.h"

#include "align.h"
#include "align/strategy.hpp"
#include "align/homology_plan.hpp"
#include "align/lossless_extractor.hpp"
#include "align/block_local_msa.hpp"
#include "align/native_msa.hpp"
#include "ultralong/native_pipeline.hpp"
#include "ultralong/window64_msa.hpp"
#include <chrono>

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

        // PR3D.4C-V2: explicit native dispatch BEFORE legacy workdir creation,
        // preprocessing, external MSA commands or other diagnostic modes.
        if (!opt.ultralong && !opt.ultralong_stream && (!opt.ultralong_report.empty()))
            throw std::invalid_argument("--ultralong-report requires --ultralong or --ultralong-stream");
        // PR3D.4D.3: 64-bit out-of-core MSA dispatcher. Dispatch before legacy
        // preprocessing and before the <=1Mb memory-resident --ultralong path.
        if (opt.ultralong_stream) {
            if (opt.ultralong || opt.mode != "legacy" || opt.keep_length || !opt.center_path.empty() ||
                !opt.msa_cmd.empty() || opt.homology_plan || opt.strict_homology_plan ||
                opt.recover_homology_plan || opt.homology_extract || opt.homology_block_align ||
                !opt.homology_plan_tsv.empty() || !opt.homology_extract_dir.empty() ||
                opt.anchored_insertions || opt.native_center || opt.auto_hierarchical ||
                opt.soft_anchors || opt.soft_anchor_recruit || !opt.workdir.empty() ||
                opt.save_workdir) {
                throw std::invalid_argument("--ultralong-stream cannot be combined with legacy or other MSA modes");
            }
            if (opt.threads < 1) throw std::invalid_argument("--threads must be positive");
            ultralong::WindowMsaConfig cfg;
            cfg.local_max_bases = opt.ultralong_max_local_bases;
            cfg.min_reference_identity = opt.ultralong_min_identity;
            // The bounded --ultralong mode keeps its original 8M default.
            // The streamed mode defaults to a 800M aggregate limit.
            cfg.max_total_bases = app.get_option("--ultralong-max-total-bases")->count()
                ? static_cast<std::uint64_t>(opt.ultralong_max_total_bases)
                : std::uint64_t{800000000};
            const auto sidecar = opt.ultralong_report.empty()
                ? opt.output + ".blocks.tsv" : opt.ultralong_report;
            const auto started = std::chrono::steady_clock::now();
            const auto stats = ultralong::runWindowMsa(opt.input,opt.output,sidecar,cfg);
            const auto elapsed = std::chrono::duration<double>(
                std::chrono::steady_clock::now()-started).count();
            spdlog::info("PR3D.4D.3 out-of-core: bases={} windows={} local_blocks={} "
                         "anchored_reference_bp={} private_bp={} msa_columns={} elapsed_sec={:.3f}",
                         stats.original_bases,stats.windows,stats.aligned_blocks,
                         stats.aligned_reference_bases,stats.private_bases,stats.output_columns,elapsed);
            spdlog::warn("Homology quality NOT evaluated in-process; review sidecar and independent pairwise F1");
            return 0;
        }
        if (opt.ultralong) {
            if (opt.mode != "legacy" || opt.keep_length || !opt.center_path.empty() ||
                !opt.msa_cmd.empty() || opt.homology_plan || opt.strict_homology_plan ||
                opt.recover_homology_plan || opt.homology_extract || opt.homology_block_align ||
                !opt.homology_plan_tsv.empty() || !opt.homology_extract_dir.empty() ||
                opt.anchored_insertions || opt.native_center || opt.auto_hierarchical ||
                opt.soft_anchors || opt.soft_anchor_recruit || !opt.workdir.empty() ||
                opt.save_workdir) {
                throw std::invalid_argument("--ultralong is standalone; cannot combine with legacy, "
                                            "planner/extractor, center or workdir options");
            }
            if (opt.threads < 1) throw std::invalid_argument("--threads must be positive");
            ultralong_native::Config cfg;
            cfg.max_local_bases = opt.ultralong_max_local_bases;
            cfg.max_total_bases = opt.ultralong_max_total_bases;
            cfg.min_reference_identity = opt.ultralong_min_identity;
            const auto sidecar = opt.ultralong_report.empty()
                ? opt.output + ".blocks.tsv" : opt.ultralong_report;
            const auto start = std::chrono::steady_clock::now();
            const auto stats = ultralong_native::run(opt.input,opt.output,sidecar,cfg);
            const double elapsed = std::chrono::duration<double>(
                std::chrono::steady_clock::now()-start).count();
            spdlog::info("PR3D.4C-V2 native preview: total_bases={} columns={} "
                         "accepted_blocks={} accepted_ref_bases={} private_bases={} elapsed_sec={:.3f}",
                         stats.original_bases,stats.output_columns,stats.accepted_blocks,
                         stats.aligned_reference_bases,stats.private_bases,elapsed);
            spdlog::warn("PR3D.4C-V2 native preview: 1Mb/read memory-resident planner, NOT out-of-core");
            return 0;
        }

        // Set default working directory if not provided, and log the final options for user confirmation.
        if (opt.workdir.empty()) {
            opt.workdir = makeDefaultWorkdir();
            spdlog::info("--workdir not provided, using default: {}", opt.workdir);
        }

        // Print parameters
        logParsedOptions(opt);
        spdlog::info("Starting halign4 version {}...", VERSION);

        // PR3D.3B: opt-in candidate-only local MSA; never touches the legacy
        // or global aligner, and never promotes UNRESOLVED to homology.
        if (!opt.homology_block_align && !opt.homology_extract_dir.empty())
            throw std::invalid_argument("--homology-extract-dir requires --homology-block-align");
        if (opt.homology_block_align) {
            if (opt.homology_extract_dir.empty() || opt.homology_extract ||
                !opt.homology_plan_tsv.empty() || opt.homology_plan ||
                opt.strict_homology_plan || opt.recover_homology_plan ||
                opt.mode != "legacy" || opt.keep_length || !opt.center_path.empty() ||
                !opt.msa_cmd.empty() || opt.anchored_insertions || opt.native_center ||
                opt.auto_hierarchical || opt.soft_anchors || opt.soft_anchor_recruit) {
                throw std::invalid_argument(
                    "--homology-block-align requires --homology-extract-dir and conflicts "
                    "with planning, extraction, routing, and global MSA flags");
            }
            const auto result = block_local_msa::align(
                opt.input, opt.homology_extract_dir, opt.output);
            spdlog::info("PR3D.3B local candidate MSAs: blocks={} aligned={} "
                         "skip_no_query={} skip_budget={} rows={} unresolved_tiles={} "
                         "(no global MSA)",
                         result.blocks, result.aligned, result.skipped_no_query,
                         result.skipped_budget, result.aligned_rows,
                         result.unresolved_tiles);
            return 0;
        }

        // PR3D.3A: standalone diagnostic extraction from an existing TSV.
        // Keep the planner and MSA codepaths unchanged; -o is a NEW directory.
        if (!opt.homology_extract && !opt.homology_plan_tsv.empty())
            throw std::invalid_argument("--homology-plan-tsv requires --homology-extract");
        if (opt.homology_extract) {
            if (opt.homology_plan_tsv.empty() || opt.homology_plan ||
                opt.strict_homology_plan || opt.recover_homology_plan ||
                opt.mode != "legacy" || opt.keep_length || !opt.center_path.empty() ||
                !opt.msa_cmd.empty() || opt.anchored_insertions || opt.native_center ||
                opt.auto_hierarchical || opt.soft_anchors || opt.soft_anchor_recruit) {
                throw std::invalid_argument(
                    "--homology-extract requires --homology-plan-tsv, and conflicts with "
                    "MSA/center/planner/routing flags (leave --mode legacy)");
            }
            const auto summary = lossless_extractor::extract(
                opt.input, opt.homology_plan_tsv, opt.output);
            spdlog::info("PR3D.3A diagnostic extraction: reads={} blocks={} candidates={} "
                         "rejected={} tiles={} unresolved_bases={} (no MSA)",
                         summary.reads, summary.blocks, summary.accepted_query_candidates,
                         summary.rejected_conflicting_candidates, summary.tiles,
                         summary.unresolved_query_bases);
            return 0;
        }

        // PR3D.1 audit-only fast exit. Do not preprocess, route, align, create
        // a work directory, launch abPOA, or write any FASTA/MSA results.
        // This branch runs before PR3C.2 routing and the old checkOption().
        if (opt.strict_homology_plan && !opt.homology_plan)
            throw std::invalid_argument("--homology-plan-strict requires --homology-plan");
        if (opt.recover_homology_plan && !opt.strict_homology_plan)
            throw std::invalid_argument("--homology-plan-recover requires --homology-plan-strict");
        if (opt.homology_plan) {
            if (opt.mode != "legacy" || opt.keep_length || !opt.center_path.empty() ||
                !opt.msa_cmd.empty() || opt.anchored_insertions || opt.native_center ||
                opt.auto_hierarchical || opt.soft_anchors || opt.soft_anchor_recruit) {
                throw std::invalid_argument(
                    "--homology-plan is standalone: omit MSA/center/routing flags "
                    "(leave --mode legacy, its default)");
            }
            const auto reads = homology_plan::readFasta(opt.input);
            homology_plan::Config planner_config;
            planner_config.strict_candidate_guard = opt.strict_homology_plan;
            planner_config.recover_indel_flanks = opt.recover_homology_plan;
            const auto plan = homology_plan::build(reads, planner_config);
            homology_plan::writeTsv(plan, opt.output);
            spdlog::info("PR3D.1 diagnostic saved {} reference blocks for {} sequences to {} "
                         "(no MSA generated)",
                         plan.boundaries.size() - 1, plan.ids.size(), opt.output);
            return 0;
        }

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
