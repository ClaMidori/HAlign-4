#include "consensus.h"

#include <algorithm>
#include <fstream>
#include <stdexcept>
#include <thread>
#include <cereal/archives/json.hpp>

#if __has_include(<omp.h>)
    #include <omp.h>
#endif

namespace consensus
{
    // Select consensus base: choose only among A/C/G/T/U (never returns N or '-').
    // Note: if all five counts are 0 (e.g., this position in input is all N or gap),
    // then fall back to 'A' by priority (can adjust to throw exception or other strategy as needed).
    char pickConsensusChar(const SiteCount& sc)
    {
        // Fixed priority (used when equal): A > C > G > T > U
        std::uint32_t best = sc.a;
        char best_ch = 'A';

        auto upd = [&](std::uint32_t v, char ch) {
            if (v > best) { best = v; best_ch = ch; }
        };

        // Only compare A/C/G/T/U, avoid selecting N or gap
        upd(sc.c, 'C');
        upd(sc.g, 'G');
        upd(sc.t, 'T');
        upd(sc.u, 'U');

        return best_ch;
    }

    // Write consensus sequence in FASTA format, line width fixed at 80
    // Note: ensures parent directory exists before writing
    void writeConsensusFasta(const FilePath& out_fasta, const std::string& seq)
    {
        file_io::ensureParentDirExists(out_fasta);

        std::ofstream ofs(out_fasta, std::ios::binary);
        if (!ofs) {
            throw std::runtime_error("failed to open fasta output: " + out_fasta.string());
        }

        ofs << ">consensus\n";
        constexpr std::size_t width = 80;
        for (std::size_t i = 0; i < seq.size(); i += width) {
            const std::size_t n = std::min(width, seq.size() - i);
            ofs.write(seq.data() + (std::streamoff)i, (std::streamsize)n);
            ofs.put('\n');
        }
    }

    // Use cereal to write counts as JSON (project prefers cereal)
    void writeCountsJson(const FilePath& out_json, const ConsensusJson& cj)
    {
        file_io::ensureParentDirExists(out_json);

        std::ofstream ofs(out_json, std::ios::binary);
        if (!ofs) {
            throw std::runtime_error("failed to open json output: " + out_json.string());
        }

        cereal::JSONOutputArchive ar(ofs);
        ar(cereal::make_nvp("consensus", cj));
    }



    // Extract single-sequence per-column counting logic into independent function.
    // This function only accumulates per-site counts from single sequence into cj.counts,
    // not responsible for updating cj.num_seqs (incremented by caller in safe context).
    //
    // Parallelization strategy: for single sequence, parallelize on position dimension (each thread handles different column index),
    // because each column's SiteCount is independent; multiple threads never write same index, so no atomic ops needed.
    // Note: if parallelizing multiple sequences (at call level), must ensure different threads don't concurrently write same column.
    //
    // Detailed optimization notes (for maintainers' reference)
    //
    // Overall goals:
    // - When processing a large number of aligned sequences, the bottleneck path of "counting bases by column" should be parallelized and vectorized as much as possible.
    //   while avoiding expensive atomic ops or frequent lock contention, maximizing cache locality and SIMD instructions.
    //
    // Main optimization techniques (visible in this file):
    // 1) Branch minimization
    //    - Traditional switch/case branches per character, misprediction costs high.
    //    - Use (idx == constant) boolean addition to eliminate branches, compiler more easily generates vector instructions.
    //
    // 2) Loop unrolling and prefetch
    //    - Via 4-way (or adjustable) loop unrolling reduce loop control overhead, increase instruction-level parallelism (ILP).
    //    - Use __builtin_prefetch to preload cache lines to be written, reduce cache miss latency (effective for long sequences/large worksets).
    //
    // 3) Thread-local accumulation + batching
    //    - Direct concurrent writes to global counts cause cache-line conflicts (false sharing). To avoid, use per-thread local count array,
    //      merge by column (reduce) after processing batch. This tradeoff of memory reduces concurrent write conflicts significantly,
    //      usually achieves best throughput with moderate thread count and alignment length.
    //    - Batch size needs tuning: too small causes scheduling/sync overhead, too large uses more memory and increases latency.
    //
    // 4) SoA (Structure of Arrays) vs AoS (Array of Structures)
    //    - AoS (SiteCount in cj.counts[i]) inconvenient for merging/vectorization, fields interleave within each SiteCount.
    //    - SoA puts each base's count in separate array (A[], C[], G[]...), facilitates sequential reads/writes per base column,
    //      more amenable to SIMD and prefetch. Implementation provides SoA batch path for maximum performance.
    //
    // 5) OpenMP + simd directives (#pragma omp parallel / #pragma omp simd / ivdep)
    //    - Use OpenMP for thread-level parallelism (assign work by position or sequence); add simd hints in inner loops to help compiler vectorize.
    //
    // 6) Compiler and compilation options
    //    - Enable -O3, -march=native, optional -flto in CMake Release mode, helps generate efficient vector instructions and inlining.
    //
    // Use cases and tradeoffs:
    // - When aln_len (alignment length) is large (thousands to tens of thousands) and batch size is large, SoA + thread-local merge usually fastest.
    // - When aln_len small (e.g., < 256), thread parallelization overhead may outweigh gains, should degrade to single thread or fewer threads.
    // - Local counts use extra memory: large thread count and alignment length cause noticeable memory usage (T * aln_len * sizeof(count)).
    // - Merge phase still scans by column; merge cost requires tradeoff between batch size and memory.
    //
    // Practical recommendations:
    // - Run small-scale benchmarks on target machine (different batch_size, threads, aln_len) to select optimal parameters.
    // - Combine with profiling tools (perf / VTune / likwid) to observe cache misses, memory bandwidth, and branch misses.
    // - For extreme performance, can further tune SoA path with AVX2/AVX512 intrinsics.
    static void processSequenceParallel(const std::string& s, ConsensusJson& cj, int thread)
    {
        const std::size_t aln_len = (std::size_t)cj.aln_len;
        if (s.size() != aln_len) {
            throw std::runtime_error("alignment length mismatch: expect " + std::to_string(aln_len) +
                                     ", got " + std::to_string(s.size()));
        }

        // Prefetch constants to reduce in-loop lookup overhead
        const unsigned char* data = reinterpret_cast<const unsigned char*>(s.data());
        std::uint8_t* base_map = const_cast<std::uint8_t*>(k_base_map.data());
        SiteCount* counts = cj.counts.data();

        const std::size_t limit = (aln_len / 4) * 4;

#if __has_include(<omp.h>)
        // Use parallel for simd, let OpenMP distribute threads and enable vectorization
        #pragma omp parallel for schedule(static) num_threads(thread)
#endif
        for (std::size_t i = 0; i < limit; i += 4) {
            // Light prefetch of future cache line, distance adjustable (16~64 bytes => ~16 positions)
            __builtin_prefetch(&counts[i + 16]);

            const std::uint8_t idx0 = base_map[data[i + 0]];
            const std::uint8_t idx1 = base_map[data[i + 1]];
            const std::uint8_t idx2 = base_map[data[i + 2]];
            const std::uint8_t idx3 = base_map[data[i + 3]];

            SiteCount& sc0 = counts[i + 0];
            SiteCount& sc1 = counts[i + 1];
            SiteCount& sc2 = counts[i + 2];
            SiteCount& sc3 = counts[i + 3];

            sc0.a += static_cast<std::uint32_t>(idx0 == 0);
            sc0.c += static_cast<std::uint32_t>(idx0 == 1);
            sc0.g += static_cast<std::uint32_t>(idx0 == 2);
            sc0.t += static_cast<std::uint32_t>(idx0 == 3);
            sc0.u += static_cast<std::uint32_t>(idx0 == 4);
            sc0.n += static_cast<std::uint32_t>(idx0 == 5);
            sc0.dash += static_cast<std::uint32_t>(idx0 == 6);

            sc1.a += static_cast<std::uint32_t>(idx1 == 0);
            sc1.c += static_cast<std::uint32_t>(idx1 == 1);
            sc1.g += static_cast<std::uint32_t>(idx1 == 2);
            sc1.t += static_cast<std::uint32_t>(idx1 == 3);
            sc1.u += static_cast<std::uint32_t>(idx1 == 4);
            sc1.n += static_cast<std::uint32_t>(idx1 == 5);
            sc1.dash += static_cast<std::uint32_t>(idx1 == 6);

            sc2.a += static_cast<std::uint32_t>(idx2 == 0);
            sc2.c += static_cast<std::uint32_t>(idx2 == 1);
            sc2.g += static_cast<std::uint32_t>(idx2 == 2);
            sc2.t += static_cast<std::uint32_t>(idx2 == 3);
            sc2.u += static_cast<std::uint32_t>(idx2 == 4);
            sc2.n += static_cast<std::uint32_t>(idx2 == 5);
            sc2.dash += static_cast<std::uint32_t>(idx2 == 6);

            sc3.a += static_cast<std::uint32_t>(idx3 == 0);
            sc3.c += static_cast<std::uint32_t>(idx3 == 1);
            sc3.g += static_cast<std::uint32_t>(idx3 == 2);
            sc3.t += static_cast<std::uint32_t>(idx3 == 3);
            sc3.u += static_cast<std::uint32_t>(idx3 == 4);
            sc3.n += static_cast<std::uint32_t>(idx3 == 5);
            sc3.dash += static_cast<std::uint32_t>(idx3 == 6);
        }

#if __has_include(<omp.h>)
        #pragma omp parallel for schedule(static) num_threads(thread)
#endif
        for (std::size_t i = limit; i < aln_len; ++i) {
            const std::uint8_t idx = base_map[data[i]];
            SiteCount& sc = counts[i];
            sc.a += static_cast<std::uint32_t>(idx == 0);
            sc.c += static_cast<std::uint32_t>(idx == 1);
            sc.g += static_cast<std::uint32_t>(idx == 2);
            sc.t += static_cast<std::uint32_t>(idx == 3);
            sc.u += static_cast<std::uint32_t>(idx == 4);
            sc.n += static_cast<std::uint32_t>(idx == 5);
            sc.dash += static_cast<std::uint32_t>(idx == 6);
        }
    }

    // Batch parallel processing: each thread maintains its own local counts array to accumulate multiple sequences,
    // merge local counts to global cj.counts after batch. This significantly reduces concurrent writes to global counts,
    // reduces false sharing and improves cache locality. Suitable for large aln_len and batch_size.
    static void processBatchParallel(const std::vector<std::string>& seqs, ConsensusJson& cj, int threads)
    {
        const std::size_t aln_len = (std::size_t)cj.aln_len;
        if (aln_len == 0 || seqs.empty()) return;

        int T = threads;
#if __has_include(<omp.h>)
        if (T <= 0) T = omp_get_max_threads();
        if (T <= 0) T = 1;
#else
        T = 1;
#endif

        // First verify all sequence lengths consistent in main thread, avoid unsafe exception in parallel region
        for (const auto& s : seqs) {
            if (s.size() != aln_len) {
                throw std::runtime_error("alignment length mismatch in batch processing");
            }
        }

        // Flatten local counts: one contiguous block, size T * aln_len
        std::vector<SiteCount> locals;
        try {
            locals.assign((std::size_t)T * aln_len, SiteCount{});
        } catch (...) {
            throw std::runtime_error("failed to allocate thread-local counts");
        }

        const std::uint8_t* base_map = k_base_map.data();

#if __has_include(<omp.h>)
        #pragma omp parallel for schedule(static) num_threads(T)
        for (std::size_t s = 0; s < seqs.size(); ++s) {
            const int tid = omp_get_thread_num();
#else
        for (std::size_t s = 0; s < seqs.size(); ++s) {
            const int tid = 0;
#endif
            const std::string& str = seqs[s];
            const unsigned char* data = reinterpret_cast<const unsigned char*>(str.data());
            SiteCount* local = locals.data() + (std::size_t)tid * aln_len;

            // Accumulate single sequence by position to thread-local count; use branch-minimizing boolean addition
            for (std::size_t i = 0; i < aln_len; ++i) {
                const std::uint8_t idx = base_map[data[i]];
                SiteCount& sc = local[i];
                sc.a += static_cast<std::uint32_t>(idx == 0);
                sc.c += static_cast<std::uint32_t>(idx == 1);
                sc.g += static_cast<std::uint32_t>(idx == 2);
                sc.t += static_cast<std::uint32_t>(idx == 3);
                sc.u += static_cast<std::uint32_t>(idx == 4);
                sc.n += static_cast<std::uint32_t>(idx == 5);
                sc.dash += static_cast<std::uint32_t>(idx == 6);
            }
        }

        // Merge local counts to global cj.counts (parallelize by column)
#if __has_include(<omp.h>)
        #pragma omp parallel for schedule(static) num_threads(T)
#endif
        for (std::size_t i = 0; i < aln_len; ++i) {
            // Use 64-bit temp accumulation to avoid overflow (even though SiteCount is uint32)
            std::uint64_t sa = 0, sc_ = 0, sg = 0, st = 0, su = 0, sn = 0, sd = 0;
            for (int t = 0; t < T; ++t) {
                const SiteCount& ls = locals[(std::size_t)t * aln_len + i];
                sa += ls.a; sc_ += ls.c; sg += ls.g; st += ls.t; su += ls.u; sn += ls.n; sd += ls.dash;
            }
            SiteCount& dst = cj.counts[i];
            dst.a += static_cast<std::uint32_t>(sa);
            dst.c += static_cast<std::uint32_t>(sc_);
            dst.g += static_cast<std::uint32_t>(sg);
            dst.t += static_cast<std::uint32_t>(st);
            dst.u += static_cast<std::uint32_t>(su);
            dst.n += static_cast<std::uint32_t>(sn);
            dst.dash += static_cast<std::uint32_t>(sd);
        }
    }

    /*
     Detailed notes on batch processing (thread-local accumulation + merge) functions:

     processBatchParallelWithLocals:
     - Each thread owns contiguous SiteCount array (locals), size aln_len.
     - Each thread accumulates its sequences to local locals (avoids concurrent writes to global cj.counts).
     - Merge phase parallelize by column: each thread merges some columns, adds local results to global counts.

     Advantages:
     - Greatly reduce false sharing (different threads don't frequently write same cache line).
     - Simple implementation, easy to understand and debug.

     Limitations:
     - Per-thread locals use more memory (T * aln_len * sizeof(SiteCount)).
     - SiteCount fields are interior AoS, merge accesses fields interleaved, vectorization limited.
    */

    // Batch processing (use externally allocated locals buffer to avoid per-batch allocation)
    static void processBatchParallelWithLocals(const std::vector<std::string>& seqs, ConsensusJson& cj, int threads, std::vector<SiteCount>& locals)
    {
        const std::size_t aln_len = (std::size_t)cj.aln_len;
        if (aln_len == 0 || seqs.empty()) return;

        int T = threads;
#if __has_include(<omp.h>)
        if (T <= 0) T = omp_get_max_threads();
        if (T <= 0) T = 1;
#else
        T = 1;
#endif

        const std::uint8_t* base_map = k_base_map.data();
        SiteCount* locals_ptr = locals.data();

#if __has_include(<omp.h>)
        #pragma omp parallel for schedule(static) num_threads(T)
        for (std::size_t s = 0; s < seqs.size(); ++s) {
            const int tid = omp_get_thread_num();
#else
        for (std::size_t s = 0; s < seqs.size(); ++s) {
            const int tid = 0;
#endif
            const std::string& str = seqs[s];
            const unsigned char* data = reinterpret_cast<const unsigned char*>(str.data());
            SiteCount* local = locals_ptr + (std::size_t)tid * aln_len;

            // More efficient per-position update: 4-way loop unroll using boolean addition, minimize branch
            std::size_t i = 0;
            const std::size_t limit = (aln_len / 4) * 4;
            for (; i < limit; i += 4) {
                __builtin_prefetch(&local[i + 16]);
                const std::uint8_t idx0 = base_map[data[i + 0]];
                const std::uint8_t idx1 = base_map[data[i + 1]];
                const std::uint8_t idx2 = base_map[data[i + 2]];
                const std::uint8_t idx3 = base_map[data[i + 3]];

                SiteCount& sc0 = local[i + 0];
                SiteCount& sc1 = local[i + 1];
                SiteCount& sc2 = local[i + 2];
                SiteCount& sc3 = local[i + 3];

                sc0.a += static_cast<std::uint32_t>(idx0 == 0);
                sc0.c += static_cast<std::uint32_t>(idx0 == 1);
                sc0.g += static_cast<std::uint32_t>(idx0 == 2);
                sc0.t += static_cast<std::uint32_t>(idx0 == 3);
                sc0.u += static_cast<std::uint32_t>(idx0 == 4);
                sc0.n += static_cast<std::uint32_t>(idx0 == 5);
                sc0.dash += static_cast<std::uint32_t>(idx0 == 6);

                sc1.a += static_cast<std::uint32_t>(idx1 == 0);
                sc1.c += static_cast<std::uint32_t>(idx1 == 1);
                sc1.g += static_cast<std::uint32_t>(idx1 == 2);
                sc1.t += static_cast<std::uint32_t>(idx1 == 3);
                sc1.u += static_cast<std::uint32_t>(idx1 == 4);
                sc1.n += static_cast<std::uint32_t>(idx1 == 5);
                sc1.dash += static_cast<std::uint32_t>(idx1 == 6);

                sc2.a += static_cast<std::uint32_t>(idx2 == 0);
                sc2.c += static_cast<std::uint32_t>(idx2 == 1);
                sc2.g += static_cast<std::uint32_t>(idx2 == 2);
                sc2.t += static_cast<std::uint32_t>(idx2 == 3);
                sc2.u += static_cast<std::uint32_t>(idx2 == 4);
                sc2.n += static_cast<std::uint32_t>(idx2 == 5);
                sc2.dash += static_cast<std::uint32_t>(idx2 == 6);

                sc3.a += static_cast<std::uint32_t>(idx3 == 0);
                sc3.c += static_cast<std::uint32_t>(idx3 == 1);
                sc3.g += static_cast<std::uint32_t>(idx3 == 2);
                sc3.t += static_cast<std::uint32_t>(idx3 == 3);
                sc3.u += static_cast<std::uint32_t>(idx3 == 4);
                sc3.n += static_cast<std::uint32_t>(idx3 == 5);
                sc3.dash += static_cast<std::uint32_t>(idx3 == 6);
            }
            for (; i < aln_len; ++i) {
                const std::uint8_t idx = base_map[data[i]];
                SiteCount& sc = local[i];
                sc.a += static_cast<std::uint32_t>(idx == 0);
                sc.c += static_cast<std::uint32_t>(idx == 1);
                sc.g += static_cast<std::uint32_t>(idx == 2);
                sc.t += static_cast<std::uint32_t>(idx == 3);
                sc.u += static_cast<std::uint32_t>(idx == 4);
                sc.n += static_cast<std::uint32_t>(idx == 5);
                sc.dash += static_cast<std::uint32_t>(idx == 6);
            }
        }

        // Merge local counts to global cj.counts (parallelize by column)
#if __has_include(<omp.h>)
        #pragma omp parallel for schedule(static) num_threads(T)
#endif
        for (std::size_t i = 0; i < aln_len; ++i) {
            std::uint64_t sa = 0, sc_ = 0, sg = 0, st = 0, su = 0, sn = 0, sd = 0;
            for (int t = 0; t < T; ++t) {
                const SiteCount& ls = locals_ptr[(std::size_t)t * aln_len + i];
                sa += ls.a; sc_ += ls.c; sg += ls.g; st += ls.t; su += ls.u; sn += ls.n; sd += ls.dash;
            }
            SiteCount& dst = cj.counts[i];
            dst.a += static_cast<std::uint32_t>(sa);
            dst.c += static_cast<std::uint32_t>(sc_);
            dst.g += static_cast<std::uint32_t>(sg);
            dst.t += static_cast<std::uint32_t>(st);
            dst.u += static_cast<std::uint32_t>(su);
            dst.n += static_cast<std::uint32_t>(sn);
            dst.dash += static_cast<std::uint32_t>(sd);
        }
    }

    /*
     processBatchParallelWithSoA (SoA version):
     - Maintain separate array per base (A[], C[], G[], T[], U[], N[], Dash[]), allocate contiguous segment per thread (T * aln_len).
     - In accumulation phase, each thread accumulates by position in own segment: a_ptr[i]++ etc. Same base counts are consecutive,
       read/write pattern more friendly to CPU vectorization and cache prefetch strategy.
     - Merge phase sums per-thread corresponding positions by column (i) and writes back to cj.counts[i].
     *
     * Advantages (why faster):
     * 1) Vectorization-friendly: SoA allows compiler to apply one instruction to consecutive a_ptr[] elements, generate SIMD instructions (e.g., AVX2).
     * 2) Reduce memory bandwidth waste: when processing base, only touch that base's array, reduce unnecessary cache writebacks.
     * 3) In merge step, sum consecutive memory across threads by column, simple access pattern, aids prefetch and hw merge.
     *
     * Cost:
     * - Extra memory overhead (7 uint32_t arrays * T * aln_len), but usually cheaper than atomics/locks/frequent cache sync.
     * - Slightly increased implementation complexity, but worth it for performance-sensitive scenarios.
    */

    static void processBatchParallelWithSoA(const std::vector<std::string>& seqs, ConsensusJson& cj, int threads,
                                            std::vector<std::uint32_t>& localsA,
                                            std::vector<std::uint32_t>& localsC,
                                            std::vector<std::uint32_t>& localsG,
                                            std::vector<std::uint32_t>& localsT,
                                            std::vector<std::uint32_t>& localsU,
                                            std::vector<std::uint32_t>& localsN,
                                            std::vector<std::uint32_t>& localsDash)
    {
        const std::size_t aln_len = (std::size_t)cj.aln_len;
        if (aln_len == 0 || seqs.empty()) return;

        int T = threads;
#if __has_include(<omp.h>)
        if (T <= 0) T = omp_get_max_threads();
        if (T <= 0) T = 1;
#else
        T = 1;
#endif

        const std::uint8_t* base_map = k_base_map.data();

#if __has_include(<omp.h>)
        #pragma omp parallel for schedule(static) num_threads(T)
        for (std::size_t s = 0; s < seqs.size(); ++s) {
            const int tid = omp_get_thread_num();
#else
        for (std::size_t s = 0; s < seqs.size(); ++s) {
            const int tid = 0;
#endif
            const std::string& str = seqs[s];
            const unsigned char* data = reinterpret_cast<const unsigned char*>(str.data());

            const std::size_t base_off = (std::size_t)tid * aln_len;
            std::uint32_t* a_ptr = localsA.data() + base_off;
            std::uint32_t* c_ptr = localsC.data() + base_off;
            std::uint32_t* g_ptr = localsG.data() + base_off;
            std::uint32_t* t_ptr = localsT.data() + base_off;
            std::uint32_t* u_ptr = localsU.data() + base_off;
            std::uint32_t* n_ptr = localsN.data() + base_off;
            std::uint32_t* d_ptr = localsDash.data() + base_off;

            const std::size_t limit = (aln_len / 4) * 4;
            std::size_t i = 0;
            for (; i < limit; i += 4) {
                __builtin_prefetch(a_ptr + i + 16);
                const std::uint8_t idx0 = base_map[data[i + 0]];
                const std::uint8_t idx1 = base_map[data[i + 1]];
                const std::uint8_t idx2 = base_map[data[i + 2]];
                const std::uint8_t idx3 = base_map[data[i + 3]];

                a_ptr[i + 0] += static_cast<std::uint32_t>(idx0 == 0);
                c_ptr[i + 0] += static_cast<std::uint32_t>(idx0 == 1);
                g_ptr[i + 0] += static_cast<std::uint32_t>(idx0 == 2);
                t_ptr[i + 0] += static_cast<std::uint32_t>(idx0 == 3);
                u_ptr[i + 0] += static_cast<std::uint32_t>(idx0 == 4);
                n_ptr[i + 0] += static_cast<std::uint32_t>(idx0 == 5);
                d_ptr[i + 0] += static_cast<std::uint32_t>(idx0 == 6);

                a_ptr[i + 1] += static_cast<std::uint32_t>(idx1 == 0);
                c_ptr[i + 1] += static_cast<std::uint32_t>(idx1 == 1);
                g_ptr[i + 1] += static_cast<std::uint32_t>(idx1 == 2);
                t_ptr[i + 1] += static_cast<std::uint32_t>(idx1 == 3);
                u_ptr[i + 1] += static_cast<std::uint32_t>(idx1 == 4);
                n_ptr[i + 1] += static_cast<std::uint32_t>(idx1 == 5);
                d_ptr[i + 1] += static_cast<std::uint32_t>(idx1 == 6);

                a_ptr[i + 2] += static_cast<std::uint32_t>(idx2 == 0);
                c_ptr[i + 2] += static_cast<std::uint32_t>(idx2 == 1);
                g_ptr[i + 2] += static_cast<std::uint32_t>(idx2 == 2);
                t_ptr[i + 2] += static_cast<std::uint32_t>(idx2 == 3);
                u_ptr[i + 2] += static_cast<std::uint32_t>(idx2 == 4);
                n_ptr[i + 2] += static_cast<std::uint32_t>(idx2 == 5);
                d_ptr[i + 2] += static_cast<std::uint32_t>(idx2 == 6);

                a_ptr[i + 3] += static_cast<std::uint32_t>(idx3 == 0);
                c_ptr[i + 3] += static_cast<std::uint32_t>(idx3 == 1);
                g_ptr[i + 3] += static_cast<std::uint32_t>(idx3 == 2);
                t_ptr[i + 3] += static_cast<std::uint32_t>(idx3 == 3);
                u_ptr[i + 3] += static_cast<std::uint32_t>(idx3 == 4);
                n_ptr[i + 3] += static_cast<std::uint32_t>(idx3 == 5);
                d_ptr[i + 3] += static_cast<std::uint32_t>(idx3 == 6);
            }
            for (; i < aln_len; ++i) {
                const std::uint8_t idx = base_map[data[i]];
                a_ptr[i] += static_cast<std::uint32_t>(idx == 0);
                c_ptr[i] += static_cast<std::uint32_t>(idx == 1);
                g_ptr[i] += static_cast<std::uint32_t>(idx == 2);
                t_ptr[i] += static_cast<std::uint32_t>(idx == 3);
                u_ptr[i] += static_cast<std::uint32_t>(idx == 4);
                n_ptr[i] += static_cast<std::uint32_t>(idx == 5);
                d_ptr[i] += static_cast<std::uint32_t>(idx == 6);
            }
        }

        // Merge: read locals arrays by column and write to cj.counts
#if __has_include(<omp.h>)
        #pragma omp parallel for schedule(static) num_threads(T)
#endif
        for (std::size_t i = 0; i < aln_len; ++i) {
            std::uint64_t sa = 0, sc_ = 0, sg = 0, st = 0, su = 0, sn = 0, sd = 0;
            const std::size_t stride = aln_len;
            for (int t = 0; t < T; ++t) {
                const std::size_t off = (std::size_t)t * stride + i;
                sa += localsA[off]; sc_ += localsC[off]; sg += localsG[off]; st += localsT[off]; su += localsU[off]; sn += localsN[off]; sd += localsDash[off];
            }
            SiteCount& dst = cj.counts[i];
            dst.a += static_cast<std::uint32_t>(sa);
            dst.c += static_cast<std::uint32_t>(sc_);
            dst.g += static_cast<std::uint32_t>(sg);
            dst.t += static_cast<std::uint32_t>(st);
            dst.u += static_cast<std::uint32_t>(su);
            dst.n += static_cast<std::uint32_t>(sn);
            dst.dash += static_cast<std::uint32_t>(sd);
        }
    }

    /*
     Deep performance optimization notes (supplementary)

     Notes for engineers wanting further optimization or debugging, contain specific advice and hardware considerations:

     1) False sharing and cache line alignment
        - False sharing occurs when multiple threads frequently write same cache line (e.g., SiteCount fields or adjacent indices in same 64B line).
        - Strategies to reduce false sharing:
          a) Thread-local accumulation (used in this file): each thread writes own local array, merge at end;
          b) Add padding to SiteCount to occupy whole cache line (increased memory);
             only consider when SiteCount updates extremely frequent and T very large; example: alignas(64) or add uint8_t pad[...] after SiteCount
          c) Let OpenMP allocate large contiguous chunk per thread (static schedule with large chunk_size), reduce interleaved writes.

     2) Vectorization and memory layout (SoA vs AoS)
        - AoS (SiteCount struct array) unfriendly to sequential single-field access, vectorization limited.
        - SoA (Structure of Arrays) puts A/C/G/T/U in separate arrays, facilitates large-range addition and SIMD generation.
        - This implementation provides SoA path: sequential writes to a_ptr/c_ptr/g_ptr in accumulation; then read by column and write back in merge.

     3) Batch size and thread count selection advice
        - First estimate max acceptable batch_size based on memory limit: batch_memory ≈ T * aln_len * sizeof(counter).
        - Typical strategy: keep batch_size large enough to amortize thread scheduling and merge cost, e.g., 1k~10k (depends on aln_len);
          but when aln_len very large (tens of thousands or million columns) batch_size can be smaller.
        - Prefer experiments: run suite on target machine (threads ∈ {1,2,4,8,...}, batch_size ∈ {64,256,1024,4096}) measure bases/sec.

     4) NUMA and affinity
        - On NUMA systems, pin threads to local NUMA node and allocate input/locals in same node (use numactl or pthread_setaffinity_np).
        - If memory bandwidth is bottleneck, consider allocate large batches to different NUMA nodes and merge separately to reduce cross-node traffic.

     5) Prefetch and warmup
        - __builtin_prefetch significantly reduces cache misses on some platforms, prefetch distance (e.g., i+16) needs tuning per cache line and access pattern.
        - Can try different prefetch distances in micro-benchmarks to find sweet spot.

     6) Memory and overflow safety
        - This implementation uses 32-bit counts (SiteCount fields); if per-column count may exceed 2^32 (very many sequences), switch to 64-bit (uint64_t) to avoid overflow.
        - In merge use uint64_t temp accumulation to avoid transient overflow, truncate back to 32-bit when writing SiteCount (or change SiteCount to uint64_t).

     7) Compiler-generated code inspection
        - To confirm vectorization works, view compiler-generated assembly with -O3 -march=native (-S or objdump), search for AVX/AVX2/AVX512 instructions.
        - Or use clang's -Rpass=loop-vectorize to have compiler report successful vectorizations.

     8) Performance regression and testing
        - After any optimization, use regression tests to ensure output consistency (project has test suite).
        - Also keep baseline (unoptimized version) for performance comparison.
    */

    // Supplementary notes on SoA batch processing: recommendations for tuning parameters and troubleshooting
    // - If merge phase takes much time, try chunking merge (e.g., merge one column range at a time for locality),
    //   or fix which columns each thread handles in merge to reduce memory thrashing.
    // - If memory bandwidth bottleneck, try reduce thread count or increase batch_size for more work per merge-writeback.

    /*
     Parallelization strategy and recommendations for generateConsensusSequence:
     - For small aln_len (e.g., < 256), prefer single thread or light parallelism; thread/merge overhead may outweigh benefits.
     - For large aln_len, use SoA + thread-local accumulation (processBatchParallelWithSoA) usually achieves best throughput.
     - batch_size value affects performance:
         * small batch: more real-time, low memory, but high sync overhead;
         * large batch: higher throughput but needs more memory and longer latency.
     - Before deployment, benchmark batch_size and threads on target machine (choose combo maximizing bases/sec).
     *
     * Practical debugging/performance collection tips:
     * - Use perf/top/htop to observe CPU utilization and memory bandwidth.
     * - Use perf record/report or VTune to see cache miss and branch miss hotspots.
     * - Compare different optimization flags (-O2 vs -O3, -march=native), confirm vectorization enabled (check compiler asm output).
    */
    std::string generateConsensusSequence(const FilePath& aligned_fasta,
                                                       const FilePath& out_fasta,
                                                       const FilePath& out_json,
                                                       std::uint64_t seq_limit,
                                                       int thread,
                                                       size_t batch_size)
    {
        file_io::requireRegularFile(aligned_fasta, "aligned_fasta");

        seq_io::KseqReader reader(aligned_fasta);

        // Read first to determine aln_len
        seq_io::SeqRecord rec;
        if (!reader.next(rec)) {
            throw std::runtime_error("aligned fasta is empty: " + aligned_fasta.string());
        }

        const std::size_t aln_len = rec.seq.size();
        if (aln_len == 0) {
            throw std::runtime_error("first sequence length is 0: " + aligned_fasta.string());
        }

        ConsensusJson cj;
        cj.aln_len = (std::uint64_t)aln_len;
        cj.counts.assign(aln_len, SiteCount{});

        std::uint64_t num_seqs = 0;

        // Batch processing params: batch_size tunable (heuristic), enlarging within memory limit reduces scheduling overhead
        std::vector<std::string> batch;
        batch.reserve(batch_size + 1);

        // Put first record in batch (and check length)
        if (rec.seq.size() != aln_len) {
            throw std::runtime_error("alignment length mismatch: first record length changed");
        }
        batch.push_back(std::move(rec.seq));
        ++num_seqs;

        // Pre-allocate and reuse thread-local buffers
        int T = thread;
#if __has_include(<omp.h>)
        if (T <= 0) T = omp_get_max_threads();
        if (T <= 0) T = 1;
#else
        if (T <= 0) T = 1;
#endif
        std::vector<SiteCount> locals;
        try {
            locals.assign((std::size_t)T * aln_len, SiteCount{});
        } catch (...) {
            throw std::runtime_error("failed to allocate thread-local counts");
        }

        // Read and process by batch
        while ((seq_limit == 0 || num_seqs < seq_limit)) {
            // Pad batch
            while (batch.size() < batch_size && (seq_limit == 0 || num_seqs < seq_limit)) {
                if (!reader.next(rec)) break;
                if (rec.seq.size() != aln_len) throw std::runtime_error("alignment length mismatch when reading");
                batch.push_back(std::move(rec.seq));
                ++num_seqs;
            }

            // Process current batch
            if (!batch.empty()) {
                // Clear locals in one shot (fast)
                std::memset(locals.data(), 0, locals.size() * sizeof(SiteCount));
                processBatchParallelWithLocals(batch, cj, thread, locals);
                batch.clear();
            }

            if (!reader.next(rec)) break; // EOF
        }


        if (cj.num_seqs == 0) cj.num_seqs = num_seqs;

        if (num_seqs == 0) {
            throw std::runtime_error("no sequences processed");
        }

        // Generate consensus sequence (single thread pick majority)
        std::string consensus_seq(aln_len, 'N');
        for (std::size_t i = 0; i < aln_len; ++i) {
            consensus_seq[i] = pickConsensusChar(cj.counts[i]);
        }

        writeConsensusFasta(out_fasta, consensus_seq);
        writeCountsJson(out_json, cj);

        return consensus_seq;
    }


} // namespace consensus
