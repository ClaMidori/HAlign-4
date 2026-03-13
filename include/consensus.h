#ifndef HALIGN4_CONSENSUS_H
#define HALIGN4_CONSENSUS_H

#include <cstddef>
#include "utils.h"
#include <cereal/cereal.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/string.hpp>


// ---------------------------------------------------------------------------
// TopKLongestSelector
// Notes:
// - Purpose: Maintain "longest K sequences" (Top-K by length) in a single streaming scan, used to select candidate set for building consensus.
// - Design goal: Low memory usage (only save K candidates), single pass (streaming) input to get result, time complexity O(N log K),
//   suitable for scenarios where N is large but K is relatively small.
// - Stability: When lengths are equal, prefer to keep earlier appearing sequences (stable tie-break). This is implemented through the order field (incrementing counter).
// - Concurrency: This class is not thread-safe; if concurrent insertion is needed in multi-threaded environment, please lock externally or use per-thread local TopK then merge strategy.
// - Memory note: Saved SeqRecord may contain long sequences, if K is large or sequences are long, pay attention to memory usage; alternatively, only save index or file offset.
// ---------------------------------------------------------------------------
class TopKLongestSelector
{
public:
    explicit TopKLongestSelector(std::size_t k = 0);

    void reset(std::size_t k);

    // Pass by value: Caller can choose to copy or std::move, internally will move to heap
    // Notes: Using pass-by-value semantics allows caller to decide whether to move SeqRecord, thus avoiding unnecessary copying.
    void consider(seq_io::SeqRecord rec);

    std::size_t size() const;
    std::size_t capacity() const;
    bool empty() const;

    // Take result (descending by length; same length ascending by input order), and clear internal state
    // Notes: Returns move semantics of SeqRecord, avoid copying large strings.
    std::vector<seq_io::SeqRecord> takeSortedDesc();

private:
    struct Item
    {
        std::size_t len{0};
        std::uint64_t order{0};   // Input order (for stabilizing tie-break)
        seq_io::SeqRecord rec;
    };

    // "Worse/smaller" judgment: shorter length is worse; same length then larger order (later appearance) is worse
    static bool worseThan(const Item& a, const Item& b);

    // Are the candidates better than the current worst (heap top)?
    static bool betterThan(const Item& cand, const Item& worst);

    void siftUp(std::size_t idx);
    void siftDown(std::size_t idx);

private:
    std::size_t k_{0};
    std::uint64_t order_counter_{0};

    // Self-implemented min-heap: heap_[0] is always the "worst" one (shortest/latest of the same length).
    std::vector<Item> heap_;
};

// ---------------------------------------------------------------------------
// consensus namespace: data structures and function declarations related to consensus sequence generation
// Contains: SiteCount (base counts for each site), ConsensusJson (for sequence count export),
// fast mapping table from character to index, and consensus generation/output interfaces.
// ---------------------------------------------------------------------------
namespace consensus
{
    // SiteCount: record counts of different categories for a site in a column
    // Field notes:
    // - a,c,g,t,u: corresponding base counts (U supports RNA/U),
    // - n: unknown base (N) count, dash: gap ('-' or '.') count
    // Performance tip: Currently uses uint32_t, if your dataset is extremely large (column counts may exceed 2^32),
    // can change type to uint64_t, but will increase memory usage.
    struct SiteCount
    {
        std::uint32_t a = 0;
        std::uint32_t c = 0;
        std::uint32_t g = 0;
        std::uint32_t t = 0;
        std::uint32_t u = 0;
        std::uint32_t n = 0;
        std::uint32_t dash = 0;

        template <class Archive>
        void serialize(Archive& ar)
        {
            ar(cereal::make_nvp("A", a),
               cereal::make_nvp("C", c),
               cereal::make_nvp("G", g),
               cereal::make_nvp("T", t),
               cereal::make_nvp("U", u),
               cereal::make_nvp("N", n),
               cereal::make_nvp("-", dash));
        }
    };

    // ConsensusJson: for serializing output consensus related statistics (e.g., write to JSON)
    // Field notes:
    // - num_seqs: total sequences participating in statistics
    // - aln_len: alignment length (length of count vector for each site)
    // - counts: SiteCount vector for each site, length is aln_len
    struct ConsensusJson
    {
        std::uint64_t num_seqs = 0;
        std::uint64_t aln_len = 0;
        std::vector<SiteCount> counts;

        template <class Archive>
        void serialize(Archive& ar)
        {
            ar(cereal::make_nvp("num_seqs", num_seqs),
               cereal::make_nvp("aln_len", aln_len),
               cereal::make_nvp("counts", counts));
        }
    };

    // -------------------- Character to index mapping table --------------------
    // Purpose: In statistics process, need to quickly map characters to 0..6 index (for array indexing and branch minimization),
    // use a 256 size lookup table (ASCII/unsigned char range) to complete mapping in constant time.
    // Mapping rules (idx): 0->A, 1->C, 2->G, 3->T, 4->U, 5->N, 6->gap('-' or '.')
    // Performance notes: Use inline constexpr initialization to ensure constant array generated at compile time, avoid runtime initialization overhead.
    inline constexpr std::array<std::uint8_t, 256> k_base_map = []() consteval {
        std::array<std::uint8_t, 256> m{};

        for (std::size_t i = 0; i < m.size(); ++i) m[i] = 5; // default N

        auto set = [&](unsigned char ch, std::uint8_t idx) { m[ch] = idx; };

        set((unsigned char)'A', 0); set((unsigned char)'a', 0);
        set((unsigned char)'C', 1); set((unsigned char)'c', 1);
        set((unsigned char)'G', 2); set((unsigned char)'g', 2);
        set((unsigned char)'T', 3); set((unsigned char)'t', 3);
        set((unsigned char)'U', 4); set((unsigned char)'u', 4);
        set((unsigned char)'N', 5); set((unsigned char)'n', 5);

        // gap: Both '-' and '.' are treated as gaps.
        set((unsigned char)'-', 6);
        set((unsigned char)'.', 6);

        return m;
    }();

    // Simple wrapper: Mapping a single character to an index
    inline std::uint8_t mapBase(char ch)
    {
        return k_base_map[(unsigned char)ch];
    }

    // -------------------- Consensus selection strategy --------------------
    // pickConsensusChar: given statistics count for a site, select a final consensus base
    // Strategy notes (current implementation suggestion):
    // - Only choose among A/C/G/T/U five (will not return N or gap), to ensure consensus sequence is as parseable nucleotide sequence as possible;
    // - When counts are equal, use fixed priority (e.g., A > C > G > T > U) to ensure reproducibility;
    // - Optional extension: if you want to return N at low coverage sites, modify strategy to return 'N' when total count is below threshold.
    char pickConsensusChar(const SiteCount& sc);

    // writeConsensusFasta: write final consensus sequence to FASTA (only write one >consensus)
    // Note: function implementation should ensure output directory exists (can use file_io::ensureParentDirExists before call)
    void writeConsensusFasta(const FilePath& out_fasta, const std::string& seq);

    // writeCountsJson: write ConsensusJson as JSON file using cereal
    void writeCountsJson(const FilePath& out_json, const ConsensusJson& cj);

    // -------------------- Consensus generation interface --------------------
    // generateConsensusSequence: given aligned FASTA file, count each site and generate consensus sequence.
    // Parameter notes:
    // - aligned_fasta: aligned FASTA (each sequence length should be consistent as aln_len)
    // - out_fasta: output consensus sequence FASTA file path
    // - out_json: output statistics count JSON file path
    // - seq_limit: if not 0, can limit processed sequence count (for debugging/sampling)
    // - thread: desired thread count (after passing, function will set OpenMP thread count internally or for thread pool size)
    // Return value: generated consensus sequence string (convenient for further processing in memory or test assertions)
    // Performance tip: implementation should support batch reading, thread local accumulation (avoid frequent atomic updates), SoA layout for vectorization, and parallelization in merge phase.
    // Also suggest recording time in implementation for benchmark analysis.
    std::string generateConsensusSequence(const FilePath& aligned_fasta,
                                                       const FilePath& out_fasta,
                                                       const FilePath& out_json,
                                                       std::uint64_t seq_limit,
                                                       int thread,
                                                       size_t batch_size = 4096);

} // namespace consensus

#endif // HALIGN4_CONSENSUS_H
