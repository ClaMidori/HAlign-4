// align.h - Main interface of the HAlign-4 sequence alignment module
// Includes: CIGAR operations, alignment algorithms (KSW2/WFA2/MM2), reference sequence aligner

#ifndef HALIGN4_ALIGN_H
#define HALIGN4_ALIGN_H
#include "utils.h"
#include "mash.h"
#include "seed.h"
#include "ksw2.h"
#include "psw.h"
#include <unordered_map>
#include <filesystem>
#include <string>
#include <vector>
#include <functional>
#include "config.hpp"
#include "consensus.h"
#include "preprocess.h"

// CIGAR Operations: Encoding, Parsing, Sequence Projection
// - Compression Format: uint32_t (28-bit length on the top + 4-bit operator on the bottom)
// - Operation Codes: 0=M, 1=I, 2=D, 3=N, 4=S, 5=H, 6=P, 7==, 8=X

// - Operation codes: 0=M, 1=I, 2=D, 3=N, 4=S, 5=H, 6=P, 7==, 8=X
namespace cigar{
    using CigarUnit = uint32_t;  // Single CIGAR operation
    using Cigar_t = std::vector<CigarUnit>;  // Multiple CIGAR operations

    // Encoding/decoding
    CigarUnit cigarToInt(char operation, uint32_t len);
    void intToCigar(CigarUnit cigar, char& operation, uint32_t& len);

    // Detect and convert
    bool hasInsertion(const Cigar_t& cigar);
    std::string cigarToString(const Cigar_t& cigar);
    Cigar_t stringToCigar(const std::string& cigar_str);

    // Sequence projection alignment
    void padQueryToRefByCigar(std::string& query, const Cigar_t& cigar);
    void delQueryToRefByCigar(std::string& query, const Cigar_t& cigar);

    // CIGAR operation
    void appendCigar(Cigar_t& result, const Cigar_t& cigar_to_add);
    std::size_t getRefLength(const Cigar_t& cigar);
    std::size_t getQueryLength(const Cigar_t& cigar);
}

// Sequence alignment: KSW2, WFA2, anchor segmentation (MM2)
namespace align {
    // ------------------------------------------------------------------
    // Alternative type names: Seed and Seed Hit
    // ------------------------------------------------------------------
    // Explanation:
    // - The seed is a short k-mer or minimizer in the sequence, used to quickly locate potentially homogeneous regions.
    // - The seed hit indicates the seed's position in the reference and query sequences.
    // - The current implementation uses a minimizer as an initialization strategy (efficient and with low memory consumption).
    // ------------------------------------------------------------------
    using SeedHit = minimizer::MinimizerHit;   // Single seed settlement: (ref_pos, query_pos, hash)
    using SeedHits = std::vector<SeedHit>;     // List of anchor points (used for positioning reference points)
    static constexpr seed::SeedKind kSeedKind = seed::SeedKind::minimizer;  // Initial strategy: minimizer

    // ------------------------------------------------------------------
    // DNA character mapping table to index (compile-time constant)
    // ------------------------------------------------------------------
    // Note: Maps DNA characters (A/C/G/T/N, case-insensitive) to indices 0-4.
    // - 'A'/'a' -> 0
    // - 'C'/'c' -> 1
    // - 'G'/'g' -> 2
    // - 'T'/'t' -> 3
    // - 'N'/'n' Or another -> 4 (unknown base)
    // Purpose: KSW2 requires that sequences be encoded as arrays of integers before alignment can be performed.
    // ------------------------------------------------------------------
// 序列比对接口：KSW2 / WFA2 / 锚点分段（MM2）
namespace align {
    // 种子命中类型（当前统一使用 minimizer）
    using SeedHit = minimizer::MinimizerHit;   // (ref_pos, query_pos, hash)
    using SeedHits = std::vector<SeedHit>;
    static constexpr seed::SeedKind kSeedKind = seed::SeedKind::minimizer;

    typedef struct ProfileMatrix{
        int len;
        int dim;
        int depth;              /* profile 总序列数 */
        std::vector<uint32_t> prof;   /* 每列 dim 个计数；前 m 个通常是 residue/base 计数 */

        ProfileMatrix() : len(0), dim(5), depth(0), prof() {}

        // 从单条序列构造 profile：每列仅一个碱基计数为 1，其余为 0。
        // 约定 A/C/G/T/N -> 0/1/2/3/4，非法字符按 N 处理，保持与项目 DNA5 语义一致。
        explicit ProfileMatrix(const std::string& seq) : len(static_cast<int>(seq.size())), dim(5), depth(seq.empty() ? 0 : 1), prof(static_cast<std::size_t>(len) * 5, 0U) {
            for (int i = 0; i < len; ++i) {
                const char ch = seq[static_cast<std::size_t>(i)];
                int idx = 4;
                switch (ch) {
                case 'A': case 'a': idx = 0; break;
                case 'C': case 'c': idx = 1; break;
                case 'G': case 'g': idx = 2; break;
                case 'T': case 't': idx = 3; break;
                case 'U': case 'u': idx = 3; break; // RNA/U 按 T 处理
                case 'N': case 'n': idx = 4; break;
                default: idx = 4; break;
                }
                prof[static_cast<std::size_t>(i) * 5 + static_cast<std::size_t>(idx)] = 1U;
            }
        }
    };




    // DNA 字符映射到 0..4（A/C/G/T/N，大小写不敏感；其他字符按 N）
    static constexpr uint8_t ScoreChar2Idx[256] = {
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,  // 0-15
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,  // 16-31
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,  // 32-47   (spaces, etc.)
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,  // 48-63   (numbers)
        4,0,4,1,4,4,4,2,4,4,4,4,4,4,4,4,  // 64-79   (@,A,B,C,D,E,F,G,H,I,J,K,L,M,N,O)
        4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,  // 80-95   (P,Q,R,S,T,U,V,W,X,Y,Z,...)
        4,0,4,1,4,4,4,2,4,4,4,4,4,4,4,4,  // 96-111  (`,a,b,c,d,e,f,g,h,i,j,k,l,m,n,o)
        4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,  // 112-127 (p,q,r,s,t,u,v,w,x,y,z,...)
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,  // 128-143
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,  // 144-159
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,  // 160-175
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,  // 176-191
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,  // 192-207
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,  // 208-223
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,  // 224-239
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4   // 240-255
    };

    // ------------------------------------------------------------------
    // DNA5 substitution matrix (compile-time constant, used for global/extended KSW2 alignment)
    // ------------------------------------------------------------------
    // Illustrating:
    // 1. Matrix dimension: 5x5, indices 0-3 are A/C/G/T, index 4 is N (unknown base/wildcard)
    // 2. One-dimensional expansion: mat[i*5+j] corresponds to two-dimensional mat[i][j]
    // 3. Scoring rules:
    //    - Exact match (A-A, C-C, G-G, T-T): +5 points (reward for correct pairing)
    //    - Mismatch (A-C, A-G, etc.): -4 points (penalty for base mismatch)
    //    - Involving N (unknown base): 0 points (neither reward nor penalty)
    // 4. Parameter balance:
    //    - The match/mismatch ratio is approximately 5:4, consistent with common Blast parameters (+1/-1 or +5/-4).
    //    - With gap_open=6 and gap_extend=2, the indel cost is approximately 6+2k.
    //    - The cost of a single mismatch is 9 (decreasing from +5 to -4), and the cost of a single 1bp gap is 8, slightly favoring gaps.
    // 5. Usage requirements:
    //    - **Must be used in conjunction with the KSW_EZ_GENERIC_SC flag (to enable the full replacement matrix)**
    //    - Without this flag, KSW2 only supports simple match/mismatch modes and will ignore this matrix.
    // ------------------------------------------------------------------
    static constexpr int8_t dna5_simd_mat[25] = {
        // A   C   G   T   N
        5, -4, -4, -4,  0,  // A (i=0)
       -4,  5, -4, -4,  0,  // C (i=1)
       -4, -4,  5, -4,  0,  // G (i=2)
       -4, -4, -4,  5,  0,  // T (i=3)
        0,  0,  0,  0,  0   // N (i=4)
 };

    // ==================================================================
    // KSW2AlignConfig: Configuration structure of the KSW2 alignment algorithm
    // ==================================================================
    // Illustrating:
    // KSW2 is a high-performance sequence alignment algorithm (SSE/AVX accelerated) that supports global, local, and extended alignment.
    // This structure encapsulates all the parameters required by KSW2.
    //
    // ==================================================================
    struct KSW2AlignConfig {
        const int8_t* mat = dna5_simd_mat;  // Replacement matrix (flattened 5x5, A/C/G/T/N)
        int alphabet_size = 5;              // Alphabet size (DNA5=5)
        int gap_open = 6;                   // Gap open penalty points (positive number, internally treated as penalty points)
        int gap_extend = 2;                 // Gap extend (penalty points, positive number)
        int end_bonus = 0;                  // End-of-line reward (used for extended comparison, can be 0)
        int zdrop = -1;                     // Z-drop (-1 indicates using the library's default settings/disabling)
        int band_width = -1;                // Bandwidth (-1 indicates the entire matrix, with no bandwidth limit)
        int flag = KSW_EZ_GENERIC_SC | KSW_EZ_RIGHT; // KSW2 flag: Enable full replacement matrix, etc.
    };

    // ==================================================================
    // Function: auto_band - Automatically estimate KSW2 bandwidth parameter
    // ==================================================================
    // Function:
    // Automatically estimates a suitable bandwidth based on sequence length and expected indel (insertion/deletion) rate.
    // Used to accelerate KSW2 alignment for long sequences (limits DP matrix calculation range).
    //
    // Parameters:
    // @param qlen - query sequence length
    // @param tlen - target (reference) sequence length
    // @param indel_rate - expected indel rate (default 0.1 = 10%)
    //                     - 0.05: high similarity sequences (e.g., within the same genome)
    //                     - 0.10: moderate similarity (e.g., closely related species)
    //                     - 0.20: low similarity (e.g., distantly related species)
    // @param margin - safety margin (default 200)
    //                 - extra buffer to prevent real indels from exceeding bandwidth and causing errors
    //
    // Returns:
    // Suggested bandwidth value (integer)
    //
    // ==================================================================
    //------------------------------------------- Bandwidth estimation
    inline int auto_band(int qlen, int tlen,
        double indel_rate = 0.1,
        int    margin = 200)
    {
        // If the difference between qlen and tlen is too large, banded DP is not suitable: return -1 to indicate "disable bandwidth limit"
        if ((double)std::abs(qlen - tlen) / (double)std::max(qlen, tlen) > 0.5)
        {
            return -1;
        }

        // Bandwidth ≈ expected indel size + safety margin (empirical formula, does not affect correctness, only speed)
        return margin + static_cast<int>(indel_rate * (qlen + tlen / 2));

    }

    // ==================================================================
    // Sequence alignment algorithm interface
    // ==================================================================
    // Note:
    // The following functions provide a unified sequence alignment interface, all returning results in CIGAR format.
    // Input sequences can contain A/C/G/T/N (case-insensitive), other characters are treated as N
    // ==================================================================

    // ------------------------------------------------------------------
    // Function: globalAlignKSW2 - KSW2 global alignment
    // ------------------------------------------------------------------
    // Function:
    // Performs global sequence alignment using the KSW2 algorithm (Needleman-Wunsch mode)
    // Application: Complete alignment of two sequences, considering all bases at both ends
    //
    // Parameters:
    // @param ref - reference sequence (string, A/C/G/T/N)
    // @param query - query sequence (string, A/C/G/T/N)
    //
    // Returns:
    // CIGAR operation sequence (cigar::Cigar_t), describes how query aligns to ref
    //
    // Configuration:
    // - Uses dna5_simd_mat (5x5 substitution matrix)
    // - gap_open=6, gap_extend=2
    // - Full matrix (no bandwidth limit)
    //
    // Performance:
    // - Time complexity: O(M * N), M and N are the lengths of the two sequences
    // - Space complexity: O(M * N) (DP matrix)
    // - Acceleration: SSE/AVX instruction set acceleration (4-8x faster than ordinary DP)
    //
    // Example:
    // globalAlignKSW2("ACGT", "AGT") -> "1M1D2M" (query is missing C from ref)
    // ------------------------------------------------------------------
    cigar::Cigar_t globalAlignKSW2(const std::string& ref, const std::string& query);

    cigar::Cigar_t globalAlignKSW2(const std::string& ref, const std::string& query, align::AlignConfig cfg);

    // ------------------------------------------------------------------
    // Function: extendAlignKSW2 - KSW2 extension alignment
    // ------------------------------------------------------------------
    // Function:
    // Performs extension alignment using the KSW2 algorithm
    // Application: Extend from seed position to both sides until score drops too much or end of sequence is reached
    //
    // Parameters:
    // @param ref - reference sequence (string, A/C/G/T/N)
    // @param query - query sequence (string, A/C/G/T/N)
    // @param zdrop - Z-drop pruning threshold (default 200)
    //                - Stop extension when alignment score drops more than zdrop
    //                - Larger value, more complete extension but slower
    //
    // Returns:
    // CIGAR operation sequence (cigar::Cigar_t), describes alignment of the extended region
    //
    // Configuration:
    // - Uses dna5_simd_mat (5x5 substitution matrix)
    // - gap_open=6, gap_extend=2
    // - end_bonus=5 (reward for reaching sequence end)
    // - Automatically estimate bandwidth (auto_band)
    //
    // Performance:
    // - Usually faster than global alignment (due to pruning and bandwidth limit)
    // - Suitable for local alignment of long sequences (e.g., overlap detection in genome assembly)
    //
    // Usage suggestion:
    // - For highly similar sequences, lower zdrop (e.g., 100) to increase speed
    // - For complete alignment, use globalAlignKSW2 instead of extendAlignKSW2
    // ------------------------------------------------------------------
    cigar::Cigar_t extendAlignKSW2(const std::string& ref, const std::string& query, int zdrop = 200);

    // ------------------------------------------------------------------
    // Function: globalAlignWFA2 - WFA2 global alignment
    // ------------------------------------------------------------------
    // Function:
    // Performs global sequence alignment using the WFA2 (Wavefront Alignment) algorithm
    // WFA2 feature: For highly similar sequences (few indels), much faster than traditional DP algorithms
    //
    // Parameters:
    // @param ref - reference sequence (string, A/C/G/T/N)
    // @param query - query sequence (string, A/C/G/T/N)
    //
    // Returns:
    // CIGAR operation sequence (cigar::Cigar_t), describes how query aligns to ref
    //
    // Performance:
    // - Time complexity: O(s * N), s is edit distance, N is sequence length
    //   - For highly similar sequences (s << N), WFA2 is much faster than KSW2
    //   - For low similarity (s ≈ N), WFA2 may be slower than KSW2
    // - Space complexity: O(s^2) (wavefront matrix)
    //
    // Usage suggestion:
    // - Recommended for highly similar sequences (e.g., same genome, closely related species)
    // - For low or unknown similarity, use MinHash to estimate similarity first
    // - If estimated similarity > 90%, use WFA2; otherwise use KSW2
    //
    // Example:
    // globalAlignWFA2("ACGT", "AGT") -> "1M1D2M"
    // ------------------------------------------------------------------
    cigar::Cigar_t globalAlignWFA2(const std::string& ref, const std::string& query);

    cigar::Cigar_t globalAlignPSW(const ProfileMatrix& ref, const std::string& query, align::AlignConfig cfg);

    // ------------------------------------------------------------------
    // Type alias: alignment function type
    // ------------------------------------------------------------------
    // Used for functions like globalAlignMM2, can pass different alignment algorithms (KSW2 or WFA2)
    // Function signature: accepts two sequences (ref, query), returns CIGAR
    using AlignFunc = std::function<cigar::Cigar_t(const std::string&, const std::string&)>;

    // ------------------------------------------------------------------
    // Function: globalAlignMM2
    // Function: piecewise global alignment based on anchors
    // Parameters:
    //   - ref: reference sequence
    //   - query: query sequence
    //   - anchors: anchor collection (used for chaining and segmentation)
    //   - align_func: optional alignment function (default uses globalAlignKSW2)
    //                 can pass globalAlignWFA2 or other functions matching the signature
    // Returns: complete CIGAR sequence
    // Notes:
    //   - Use anchors to decompose sequences into multiple small fragments, align separately then concatenate
    //   - Through align_func parameter, can flexibly choose underlying alignment algorithm
    //   - If anchors are empty or invalid, falls back to align_func's global alignment
    // ------------------------------------------------------------------

    // 基于锚点的分段全局比对；锚点无效时退化为普通全局比对
    cigar::Cigar_t globalAlignSeq2Seq(const std::string& ref,
                                      const std::string& query,
                                      const anchor::Anchors& anchors);

    cigar::Cigar_t globalAlignSeq2Profile(const ProfileMatrix& ref,
                                  const std::string& query,
                                  const anchor::Anchors& anchors,
                                  int thread = 1);


    // ==================================================================
    // RefAligner class: High-performance reference sequence aligner (multiple sequence alignment MSA engine)
    // ==================================================================
    // Overview:
    // RefAligner is the core component of HAlign-4, responsible for aligning a large number of query sequences to the reference set
    // and generating the final multiple sequence alignment (MSA) result
    // Usage example:
    // ------------------------------------------------------------------
    // // 1. Create RefAligner
    // Options opt = parseCommandLine(argc, argv);
    // RefAligner aligner(opt, "refs.fasta");
    //
    // // 2. Align query sequences
    // aligner.alignQueryToRef("queries.fasta", 5120);
    //
    // // 3. Merge results to generate MSA
    // aligner.mergeAlignedResults("consensus_aligned.fasta", "mafft --auto");
    // ------------------------------------------------------------------
    // ==================================================================
    class RefAligner
    {
        public:
        // ------------------------------------------------------------------
        // Constructor 1: Direct parameter initialization
        //
        // Interface semantics (same as current implementation, and what the caller needs to understand):
        // - Reads the reference sequence set from ref_fasta_path and builds a sketch/minimizer index for each reference sequence;
        // - Generates a "central/consensus sequence" as the anchor for subsequent merge coordinate system:
        //   * keep_first_length=true  : central sequence is the first reference sequence (ref_sequences[0]);
        //   * keep_first_length=false : central sequence is generated by external MSA + generateConsensusSequence.
        //
        // Parameters:
        //   - work_dir: working directory (location for intermediate and output files)
        //   - ref_fasta_path: reference sequence FASTA
        //   - kmer_size/window_size/sketch_size: MinHash/Minimizer parameters
        //   - noncanonical: whether to enable reverse complement (true for DNA)
        //   - threads: OpenMP thread count (and for external MSA)
        //   - msa_cmd: external MSA command template (for consensus and insertion alignment)
        //   - keep_first_length/keep_all_length: final MSA column trimming strategy (see member comments)
        // ------------------------------------------------------------------
        RefAligner(const FilePath& work_dir, const FilePath& ref_fasta_path,
                   int kmer_size = 21, int window_size = 10,
                   int sketch_size = 2000, bool noncanonical = true,
                   int threads = 1, std::string msa_cmd = "",
                   bool keep_length = false,
                   bool enable_wfa = false);

        // ------------------------------------------------------------------
        // Constructor 2: Initialization based on Options struct (recommended)
        // ------------------------------------------------------------------
        // Parameters:
        //   @param opt - Options struct (contains all command line arguments and configuration)
        //                - Automatically extracts: threads, kmer_size, window_size, sketch_size,
        //                  noncanonical, msa_cmd, keep_first_length, keep_all_length
        //   @param ref_fasta_path - reference sequence FASTA file path
        //
        // Workflow:
        //   1. Read reference sequence file into ref_sequences
        //   2. Generate MinHash sketch for each reference sequence (for fast similarity estimation)
        //   3. Generate Minimizer index for each reference sequence (for fast homology region localization)
        //   4. If opt.consensus_string is not empty, use it directly; otherwise generate consensus sequence
        //   -------------------------
        RefAligner(const Options& opt, const FilePath& ref_fasta_path);

        // ==================================================================
        // Public methods: alignment and merging
        // ==================================================================

        // ------------------------------------------------------------------
        // Method: alignQueryToRef - Parallel alignment of query sequences to reference sequences
        // ------------------------------------------------------------------
        // Function:
        // Reads the query FASTA file and aligns each sequence to the most similar reference sequence
        // Outputs multiple SAM files (one per thread), recording the alignment results
        //
        // Parameters:
        //   @param qry_fasta_path - query sequence FASTA file path
        //   @param batch_size - number of sequences to read per batch (default 5120)
        //                       - Larger value, higher throughput, but more memory usage
        //                       - Suggestion: 1024-5120 for normal machines, 10000+ for high-performance servers
        //
        // ------------------------------------------------------------------
        // Note:
        // - threads <= 0 means using the default OpenMP runtime thread count (e.g., controlled by OMP_NUM_THREADS)
        // - batch_size controls the "streaming read" batch size; larger means higher throughput but more memory usage
        void alignQueryToRef(const FilePath& qry_fasta_path, std::size_t batch_size = 25600);

        // ------------------------------------------------------------------
        // Method: mergeAlignedResults - Merge all alignment results to generate the final MSA
        // ------------------------------------------------------------------
        // Function:
        // Merges the SAM files produced by multiple threads into a unified multiple sequence alignment (MSA) FASTA file
        // Includes secondary alignment of insertion sequences and coordinate unification
        //
        // Parameters:
        //   @param msa_cmd - external MSA tool command template (for aligning insertion sequences)
        //                    - Example: "mafft --auto {input} > {output}"
        //                    - Supported tools: MAFFT, Muscle, Clustal Omega
        //   @param batch_size - batch size: controls the number of sequences processed in parallel per batch
        //                       - Larger batch_size increases throughput but uses more memory
        //                       - Default value 1000 is suitable for most scenarios
        //   @param thread - number of parallel threads: for OpenMP parallel processing within the batch
        //                   - Default value 4, recommended to set to the number of CPU cores
        // ------------------------------------------------------------------
        void mergeAlignedResults(const FilePath output, std::size_t batch_size = 25600);

        // 全局比对统一入口（保留 similarity/minimizer 参数以兼容后续策略）
        cigar::Cigar_t Seq2SeqWithAnchor(const std::string& ref,
                                          const std::string& query,
                                          double similarity,
                                          const SeedHits* ref_minimizer = nullptr,
                                          const SeedHits* query_minimizer = nullptr) const;

        // ------------------------------------------------------------------
        // Helper function: removeRefGapColumns
        // Function: Delete columns where the reference is a gap according to ref_gap_pos (in-place modification)
        //
        // Usage scenario:
        // - The input sequence seq is usually a sequence that has already been MSA-aligned (contains gaps)
        // - ref_gap_pos records whether each column of the reference (first sequence) is a gap
        // - This function deletes all columns where ref_gap_pos[i]==true and keeps the rest
        //
        // Parameters:
        //   - seq: input/output sequence (aligned, contains gap '-') [in-place modification]
        //   - ref_gap_pos: for the reference (first sequence), whether each column is a gap; true means the column should be deleted
        // ------------------------------------------------------------------
        cigar::Cigar_t Seq2ProfileWithAnchor(const ProfileMatrix& ref,
                                    const std::string& ref_string,
                                   const std::string& query,
                                   double similarity,
                                   const SeedHits* ref_minimizer = nullptr,
                                   const SeedHits* query_minimizer = nullptr) const;

        // 删除“参考为 gap”的列（原地修改）
        static void removeRefGapColumns(
            std::string& seq,
            const std::vector<bool>& ref_gap_pos);


        private:
        // ------------------------------------------------------------------
        // Function: alignOneQueryToRef
        // Function: Align a single query and write to SAM file
        // Parameters:
        //   - q: query sequence to align
        //   - out: writer for normal output file (no insertion or no insertion after secondary alignment)
        //   - out_insertion: writer for insertion output file (still has insertion after secondary alignment)
        // Note:
        //   - This function does not modify shared reference data structures (thread-safe)
        //   - out and out_insertion are managed by the caller (per-thread)
        // ------------------------------------------------------------------
        void alignOneQueryToRef(const seq_io::SeqRecord& q,
                               seq_io::SeqWriter& out,
                               seq_io::SeqWriter& out_insertion) const;

        void alignOneQueryToProfile(const seq_io::SeqRecord& q,
                               seq_io::SeqWriter& out,
                               seq_io::SeqWriter& out_insertion) const;

        // Helper function: write SAM record (select correct reference name and output file)
        void writeSamRecord(const seq_io::SeqRecord& q, const cigar::Cigar_t& cigar,
                           std::string_view ref_name, seq_io::SeqWriter& out) const;

        // ------------------------------------------------------------------
        // Helper function: mergeConsensusAndSamToFasta
        // Function: Merge consensus_seq and a set of sequences from SAM files into a FASTA file.
        //
        // Meaning of keep (same as current implementation):
        // - keep=false: directly write query sequences from SAM to FASTA (no CIGAR projection)
        // - keep=true : project query to reference coordinates according to SAM CIGAR (currently "delete query insertions I",
        //              making its length closer to reference/consensus), for subsequent MSA/merge.
        //
        // Parameters:
        //   - sam_paths: list of input SAM file paths
        //   - fasta_path: output FASTA file path
        //   - keep: whether to project according to CIGAR (see above)
        //   - line_width: FASTA line width (default 80)
        // Returns: total number of sequences written (including consensus_seq at the beginning, at least 1)
        // ------------------------------------------------------------------
        std::size_t mergeConsensusAndSamToFasta(
            const std::vector<FilePath>& sam_paths,
            const FilePath& fasta_path,
            std::unordered_map<std::string, cigar::Cigar_t> ref_aligned_map,
            bool keep = false,
            std::size_t line_width = 80
            ) const;

        // ------------------------------------------------------------------
        // Helper function: convertSamToFastaRecord
        // Function: Convert a single SAM record to a FASTA record, and adjust sequence length according to CIGAR
        // ------------------------------------------------------------------
        void convertSamToFastaRecord(
            const seq_io::SamRecord& sam_rec,
            seq_io::SeqRecord& fasta_rec,
            const std::unordered_map<std::string, cigar::Cigar_t>& ref_aligned_map,
            std::size_t estimated_final_length) const;

        // ------------------------------------------------------------------
        // Helper function: processInsertionSequences
        // Function: Read insertion SAM, merge to FASTA, perform optional MSA
        // Returns: insertion sequence FASTA file path
        // ------------------------------------------------------------------
        FilePath processInsertionSequences(
            const FilePath& result_dir,
            const FilePath& aligned_insertion_fasta,
            std::unordered_map<std::string, cigar::Cigar_t>& ref_aligned_map) const;

        // ------------------------------------------------------------------
        // Helper function: writeConsensusAndReferences
        // Function: Read from alignment file and write consensus and reference sequences
        // Returns: number of sequences written
        // ------------------------------------------------------------------
        std::size_t writeConsensusAndReferences(
            seq_io::SeqWriter& final_writer,
            const FilePath& consensus_aligned_file,
            ProgressBar& progress) const;

        // ------------------------------------------------------------------
        // Helper function: writeInsertionSequences
        // Function: Read from aligned insertion file and write sequences (skip the first consensus)
        // Returns: number of sequences written
        // ------------------------------------------------------------------
        std::size_t writeInsertionSequences(
            seq_io::SeqWriter& final_writer,
            const FilePath& aligned_insertion_fasta,
            std::size_t& expected_length,
            bool& length_initialized,
            ProgressBar& progress) const;

        // ------------------------------------------------------------------
        // Helper function: processSamFileBatch
        // Function: Read SAM batch, convert to FASTA in parallel, write output serially
        // ------------------------------------------------------------------
        void processSamFileBatch(
            seq_io::SamReader& sam_reader,
            const std::size_t batch_size,
            seq_io::SeqWriter& final_writer,
            const std::unordered_map<std::string, cigar::Cigar_t>& ref_aligned_map,
            std::size_t estimated_final_length,
            std::size_t& expected_length,
            bool& length_initialized,
            std::size_t& seq_count,
            ProgressBar& progress) const;

        // ------------------------------------------------------------------
        // Helper function: parseAlignedReferencesToCigar
        // Function: Read the reference sequence file after MSA alignment and generate the CIGAR for each sequence (only M and D)
        //
        // Important changes (interface convention):
        // 1) No longer returns map as a return value, but outputs via parameter (avoids large object return/move, clearer for the caller)
        // 2) Added ref_gap_pos: marks whether each column of the aligned reference sequence is a gap ('-')
        //    - Here, "reference sequence" refers to the first sequence in the alignment file (usually consensus or central sequence)
        //    - ref_gap_pos[i] == true  means column i is a gap in the reference
        //    - ref_gap_pos[i] == false means column i is a base in the reference
        //
        // Note:
        // - For each subsequent sequence (reference/insertion), only its own character is encoded: base -> M, gap -> D
        // - No base consistency check (e.g., with the first sequence's column), keeps current logic equivalent
        // ------------------------------------------------------------------
        void parseAlignedReferencesToCigar(
            const FilePath& aligned_fasta_path,
            std::unordered_map<std::string, cigar::Cigar_t>& out_ref_aligned_map,
            std::vector<bool>& out_ref_gap_pos) const;


        // ==================================================================
        // Private member variables
        // ==================================================================

        // ------------------------------------------------------------------
        // Working directory and data paths
        // ------------------------------------------------------------------
        FilePath work_dir;  // Working directory (stores intermediate files and final results)

        // ------------------------------------------------------------------
        // Reference sequences and indexes
        // ------------------------------------------------------------------
        seq_io::SeqRecords ref_sequences;   // Reference sequence set (read from ref_fasta_path)
        mash::Sketches ref_sketch;          // MinHash sketch for each reference sequence (for fast similarity estimation)
        std::vector<SeedHits> ref_minimizers;  // Minimizer index for each reference sequence (for fast homology region localization)

        // ------------------------------------------------------------------
        // Consensus sequence
        // ------------------------------------------------------------------
        // Note:
        // - The consensus sequence is the "representative sequence" of all reference sequences
        // - If keep_first_length=false, generated by MSA tool (e.g., MAFFT)
        // - If keep_first_length=true, directly use ref_sequences[0]
        // - Used as the anchor for the coordinate system in mergeAlignedResults
        // ------------------------------------------------------------------
        seq_io::SeqRecord consensus_seq;  // Consensus sequence (stored as a member variable)
        mash::Sketch consensus_sketch;    // MinHash sketch of the consensus sequence (for similarity calculation in secondary alignment)
                          // Note: Initialized once in the constructor to avoid recalculating in every alignOneQueryToRef
        SeedHits consensus_minimizer;     // Minimizer index of the consensus sequence (for seed localization in secondary alignment)
                          // Note: Initialized once in the constructor to avoid recalculating in every alignOneQueryToRef

        // ------------------------------------------------------------------
        // MinHash and Minimizer parameters
        // ------------------------------------------------------------------
        int kmer_size = 21;      // k-mer size (shared by MinHash and Minimizer)
        int window_size = 10;    // Minimizer window size
        int sketch_size = 2000;  // MinHash sketch size (number of hashes)
        int random_seed = 42;    // Random seed (for MinHash hash function)

        // ------------------------------------------------------------------
        // Threads and external tool configuration
        // ------------------------------------------------------------------
        int threads = 1;            // OpenMP thread count (<=0 usually means let runtime decide; see .cpp for logic)
        std::string msa_cmd;        // External MSA command template (for consensus generation and insertion sequence MSA)

        bool keep_length = false; // true: trim columns where consensus is a gap, keep the original length of the central sequence

        // ------------------------------------------------------------------
        // MinHash calculation options
        // ------------------------------------------------------------------
        // Note:
        // - noncanonical: whether to use non-canonical k-mer (no strand distinction)
        //   * true: consider reverse complement of k-mer (for DNA sequences)
        //   * false: only consider original k-mer (for protein or single-stranded RNA)
        // ------------------------------------------------------------------
        bool noncanonical = true;   // Whether minimizer/sketch enables reverse complement (true for DNA)

        // ------------------------------------------------------------------
        // Output file paths
        // ------------------------------------------------------------------
        // Note:
        // - outs_path: normal SAM output file path for each thread
        //   * length = threads
        //   * outs_path[i] = {work_dir}/data/aligned_thread_{i}.sam
        // - outs_with_insertion_path: insertion SAM output file path for each thread
        //   * length = threads
        //   * outs_with_insertion_path[i] = {work_dir}/data/aligned_insertion_thread_{i}.sam
        // ------------------------------------------------------------------
        std::vector<FilePath> outs_path;
        std::vector<FilePath> outs_with_insertion_path;
    };

} // namespace align

#endif //HALIGN4_ALIGN_H
