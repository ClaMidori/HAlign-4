// mash.cpp
// -----------------
// This file implements a simplified Mash-style sketch and related distance measure functions
// (Jaccard, Mash distance, ANI inference, etc).
// We document design considerations, parameter meanings, edge cases and implementation details
// in comments, to facilitate performance testing and maintenance.

#include "mash.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <queue>
#include <cstring> // for memcmp
#include "hash.h"
#include "robin_hood.h"
#include <unordered_map>
#include <unordered_set>
namespace mash
{
    // clamp01: helper function to clamp double value to [0,1] range
    // ensures numerical stability for returning probability/ratio-type values
    static inline double clamp01(double x) noexcept
    {
        if (x < 0.0) return 0.0;
        if (x > 1.0) return 1.0;
        return x;
    }

    Sketch sketchFromSequence(const std::string& seq,
                              std::size_t k,
                              std::size_t sketch_size,
                              bool noncanonical,
                              int seed)
    {
        Sketch sk;
        sk.k = k;
        sk.noncanonical = noncanonical;

        // --- Parameter/boundary checks ---
        // k=0 or sketch_size=0: meaningless, return empty sketch.
        // seq.size() < k: insufficient to form any k-mer.
        if (k == 0 || sketch_size == 0 || seq.size() < k) return sk;
        // This implementation uses 2-bit encoding to compress k-mer into 64 bits, thus limit k<=32.
        if (k > 32) return sk;

        // mask: retain low 2*k bits of rolling encoding
        const std::uint64_t mask = (1ULL << (2 * k)) - 1ULL;
        // shift: position for placing new base at highest two bits when maintaining reverse-complement rolling
        const std::uint64_t shift = 2ULL * (k - 1);

        // fwd: forward rolling 2-bit encoding
        // rev: reverse-complement rolling 2-bit encoding
        // valid: accumulated valid bases in current window (clears on N/invalid char)
        std::uint64_t fwd = 0;
        std::uint64_t rev = 0;
        std::size_t valid = 0;

        // max-heap (top is current set's "maximum" hash), maintains bottom-k (smallest k hashes)
        // complexity: each k-mer only O(log(sketch_size))
        std::priority_queue<hash_t> maxHeap;

        // seen: deduplication to avoid same hash entering heap multiple times (prevents bias)
        // here uses robin_hood::unordered_set (typically faster and more memory-efficient than std::unordered_set)
        // reserve *2: reduces rehash frequency
        std::unordered_set<hash_t> seen;
        seen.reserve(sketch_size * 2 + 1);

        // Main loop: roll through sequence, generate k-mer at each position (forward/reverse-complement),
        // then compute hash
        for (std::size_t i = 0; i < seq.size(); ++i)
        {
            // nt4_table: maps base to 0/1/2/3; other chars (e.g. N) map to 4
            const uint8_t c = nt4_table[static_cast<unsigned char>(seq[i])];
            if (c >= 4)
            {
                // Encountered N/invalid char: current rolling window becomes invalid, must re-accumulate from next char
                fwd = rev = 0;
                valid = 0;
                continue;
            }

            // forward rolling: left-shift 2 bits, add new base (low bits), then truncate to 2*k bits
            fwd = ((fwd << 2) | c) & mask;
            // reverse-complement rolling: right-shift 2 bits, put "complement base" at high bits
            // complement relation: A(0)<->T(3), C(1)<->G(2), corresponds to 3^c
            rev = (rev >> 2) | (std::uint64_t(3U ^ c) << shift);

            // valid count insufficient for k => cannot form complete k-mer yet
            if (valid < k) ++valid;
            if (valid < k)
            {
                continue;
            }

            // canonical selection:
            // - noncanonical=true: only take forward
            // - noncanonical=false: take min(fwd, rev) as canonical k-mer
            const std::uint64_t code = noncanonical ? fwd : std::min(fwd, rev);

            // Hash the 2-bit encoded k-mer directly to 64-bit
            // seed used for perturbing hash (different seed produces different sketch)
            const hash_t h = getHash2bit(code, static_cast<std::uint32_t>(seed));

            // Deduplication: skip if this hash was seen before
            if (!seen.insert(h).second)
            {
                continue;
            }

            if (maxHeap.size() < sketch_size)
            {
                // heap not full yet: insert directly
                maxHeap.push(h);
            }
            else if (h < maxHeap.top())
            {
                // heap full: only qualify for bottom-k if new hash is smaller
                // first pop current maximum and erase from seen; then insert new value
                seen.erase(maxHeap.top());
                maxHeap.pop();
                maxHeap.push(h);
            }
            // else: new hash not small enough, discard directly
            // (but note: seen already inserted it, must undo)
            else
            {
                seen.erase(h);
            }
        }

        // Export heap to vector (order not guaranteed at this point)
        sk.hashes.reserve(maxHeap.size());
        while (!maxHeap.empty())
        {
            sk.hashes.push_back(maxHeap.top());
            maxHeap.pop();
        }

        // Final sort is mandatory:
        // - intersection/jaccard require sorted unique
        // - other places also rely on sketch being sorted
        std::sort(sk.hashes.begin(), sk.hashes.end());

        return std::move(sk);
    }

    // jaccard: compute Jaccard similarity based on two sorted and unique sketch vectors
    // returns range [0,1], special case: two empty sets return 1.0
    // (definition choice: represents complete identity of empty sets)
    double jaccard(const Sketch& a, const Sketch& b)
    {
        if (a.k != b.k) throw std::invalid_argument("mash::jaccard: mismatched k");

        if (a.hashes.empty() && b.hashes.empty()) return 1.0;
        if (a.hashes.empty() || b.hashes.empty()) return 0.0;

        const std::size_t inter = intersectionSizeSortedUnique(a.hashes, b.hashes);
        const std::size_t uni = std::min(a.hashes.size(), b.hashes.size());
        if (uni == 0) return 1.0;
        return static_cast<double>(inter) / static_cast<double>(uni);
    }

    // mashDistanceFromJaccard: convert Jaccard similarity to Mash distance (based on Mash's mathematical derivation)
    // formula: x = (2*j)/(1+j) ; distance = -ln(x) / k
    // note edge cases: when j<=0 or x<=0 return +inf; when j>=1 return 0
    double mashDistanceFromJaccard(double j, std::size_t k)
    {
        if (k == 0) throw std::invalid_argument("mash::mashDistanceFromJaccard: k must be > 0");
        if (!(j > 0.0)) return std::numeric_limits<double>::infinity();
        if (j >= 1.0) return 0.0;

        const double x = (2.0 * j) / (1.0 + j);
        if (!(x > 0.0)) return std::numeric_limits<double>::infinity();
        return -std::log(x) / static_cast<double>(k);
    }

    // aniFromJaccard: estimate average nucleotide identity (ANI) from Jaccard using Mash's approximate relationship
    // formula: x = (2*j)/(1+j) ; ANI ~ x^(1/k)
    // return value clamped to [0,1]
    double aniFromJaccard(double j, std::size_t k)
    {
        if (k == 0) throw std::invalid_argument("mash::aniFromJaccard: k must be > 0");
        if (!(j > 0.0)) return 0.0;
        if (j >= 1.0) return 1.0;

        const double x = (2.0 * j) / (1.0 + j);
        if (!(x > 0.0)) return 0.0;

        return clamp01(std::pow(x, 1.0 / static_cast<double>(k)));
    }

    // aniFromMashDistance: infer ANI from Mash distance: ANI ~ exp(-d)
    // handles protections for infinity/non-finite inputs, limits output to [0,1]
    double aniFromMashDistance(double d)
    {
        if (!std::isfinite(d)) return 0.0;
        if (d <= 0.0) return 1.0;
        return clamp01(std::exp(-d));
    }

    // intersectionSizeSortedUnique
    // compute intersection size of two sorted, unique (no duplicates) hash vectors.
    // algorithm uses classic two-pointer linear scan.
    // requires input a,b already satisfy ascending order and no duplicates;
    // function doesn't modify input, only returns count of intersection elements.
    // time complexity: O(|a| + |b|). suitable for comparing sketch hash sets.
    // note: in this project hash_t = uint64_t, no need to branch on use64 etc
    std::size_t intersectionSizeSortedUnique(const std::vector<hash_t>& a,
                                             const std::vector<hash_t>& b) noexcept
    {
        std::size_t i = 0, j = 0, inter = 0;
        while (i < a.size() && j < b.size()) {
            const auto av = a[i];
            const auto bv = b[j];
            if (av == bv) {
                ++inter;
                ++i;
                ++j;
            } else if (av < bv) {
                ++i;
            } else {
                ++j;
            }
        }
        return inter;
    }

    bloom_filter filterFromSketch(const Sketch& sk, double false_positive_rate, int seed)
    {
        // ------------------------------------------------------------
        // Build Bloom Filter from sketch
        // ------------------------------------------------------------
        // Design goals:
        // 1) Make subsequent contains(hash) queries as fast as possible;
        // 2) false_positive_rate controllable (approximate);
        // 3) friendly to empty sketch: return default-constructed bloom_filter (operator!() is true).
        //
        // bloom_filter.hpp parameter model:
        // - projected_element_count: estimated number of elements to insert
        // - false_positive_probability: target false positive rate
        // - random_seed: affects salt generation, thus hash function family
        //
        // note: inserting hash_t (uint64_t) itself here, i.e., treating already-hashed values as key.
        // differs from directly inserting original k-mer sequences, but reasonable for set membership queries.
        // ------------------------------------------------------------

        bloom_filter bf;

        if (sk.hashes.empty())
        {
            return bf;
        }

        bloom_parameters p;
        p.projected_element_count = static_cast<unsigned long long>(sk.hashes.size());
        p.false_positive_probability = false_positive_rate;
        // Avoid random_seed falling into invalid range in bloom_filter.hpp (0 or all 1s)
        {
            const std::uint64_t s = static_cast<std::uint64_t>(static_cast<std::uint32_t>(seed));
            p.random_seed = (s == 0ULL) ? 0xA5A5A5A55A5A5A5AULL : (0xA5A5A5A55A5A5A5AULL ^ (s * 0x9E3779B97F4A7C15ULL));
            if (p.random_seed == 0ULL) p.random_seed = 0x1ULL;
            if (p.random_seed == 0xFFFFFFFFFFFFFFFFULL) p.random_seed = 0xFFFFFFFFFFFFFFFEULL;
        }

        // Compute optimal parameters for bloom filter (bit count and hash function count)
        if (!p.compute_optimal_parameters())
        {
            // Parameters invalid/computation failed => degrade to default-constructed (empty)
            return bloom_filter{};
        }

        bf = bloom_filter(p);

        // Insert all hash values
        for (const hash_t hv : sk.hashes)
        {
            bf.insert(hv);
        }

        return bf;
    }

    double jaccard(const bloom_filter& a, const Sketch& b)
    {
        // ------------------------------------------------------------
        // BloomFilter vs Sketch Jaccard (approximate)
        // ------------------------------------------------------------
        // In Sketch-Sketch implementation we use: inter / min(|A|,|B|)
        // This is approximate definition commonly used in Mash/MinHash scenarios (when both sides are bottom-k sets).
        //
        // Here a is BloomFilter (built from some sketch), it doesn't contain "set size information" itself,
        // therefore:
        // - intersection: estimated by contains query for each hv in b.hashes
        // - union size: still uses min(|A|,|B|) convention, but |A| takes a.element_count()
        //   (number of elements inserted during construction; note it may exceed true unique count,
        //
        // Edge cases:
        // - both empty => 1.0
        // - any empty => 0.0
        // ------------------------------------------------------------

        const std::size_t asz = static_cast<std::size_t>(a.element_count());
        const std::size_t bsz = b.hashes.size();

        if (asz == 0 && bsz == 0) return 1.0;
        if (asz == 0 || bsz == 0) return 0.0;

        std::size_t inter = 0;
        for (const hash_t hv : b.hashes)
        {
            if (a.contains(hv)) ++inter;
        }

        const std::size_t uni = std::min(asz, bsz);
        if (uni == 0) return 1.0;
        return static_cast<double>(inter) / static_cast<double>(uni);
    }

} // namespace mash
