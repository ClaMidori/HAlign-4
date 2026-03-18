#include <doctest/doctest.h>

#include <mash.h>

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <random>
#include <string>
#include <vector>

// ------------------------------------------------------------
// test_jaccard.cpp
// ------------------------------------------------------------
// This test file verifies correctness and performance of two Jaccard computations:
// 1) mash::jaccard(const Sketch&, const Sketch&)
//    - exact set intersection count (requires hashes to be sorted and unique).
// 2) mash::jaccard(const bloom_filter&, const Sketch&)
//    - uses a Bloom Filter for approximate membership queries to estimate intersection; affected by false positives.
//
// Performance evaluation requested:
// - fix a reference sketch
// - construct 10,000 query sketches (each size 2k)
// - measure time for 10k "query vs ref" Jaccard computations
//
// Note: perf cases are skipped by default to avoid slow CI/normal runs.
// Enable with env var: HALIGN4_RUN_PERF=1
// ------------------------------------------------------------

namespace {

// Generate a random DNA sequence (A/C/G/T).
static std::string random_dna(std::mt19937_64& rng, std::size_t len)
{
    static constexpr char bases[4] = {'A', 'C', 'G', 'T'};
    std::string s;
    s.reserve(len);
    for (std::size_t i = 0; i < len; ++i)
        s.push_back(bases[rng() & 3ULL]);
    return s;
}

} // namespace

TEST_SUITE("jaccard") {
// Note: in some environments doctest's TEST_SUITE macro is sensitive to newline/brace placement,
// so we use the same style as other test files in this project.

    TEST_CASE("Sketch-Sketch jaccard: basic correctness")
    {
        mash::Sketch a;
        mash::Sketch b;
        a.k = b.k = 21;

        SUBCASE("both empty")
        {
            CHECK(mash::jaccard(a, b) == doctest::Approx(1.0));
        }

        SUBCASE("one empty")
        {
            a.hashes = {1, 2, 3};
            CHECK(mash::jaccard(a, b) == doctest::Approx(0.0));
            CHECK(mash::jaccard(b, a) == doctest::Approx(0.0));
        }

        SUBCASE("identical")
        {
            a.hashes = {10, 20, 30};
            b.hashes = {10, 20, 30};
            CHECK(mash::jaccard(a, b) == doctest::Approx(1.0));
        }

        SUBCASE("partial overlap - uses min(|A|,|B|) denominator")
        {
            // intersection=2, min size=3 => 2/3
            a.hashes = {1, 2, 3};
            b.hashes = {2, 3, 4, 5};
            CHECK(mash::jaccard(a, b) == doctest::Approx(2.0 / 3.0));
        }
    }

    TEST_CASE("Sketch-Sketch jaccard: sketchFromSequence identical sequences should be 1")
    {
        // If you see "always 0", the most common causes are:
        // 1) sketch.hashes is not sorted/unique, breaking intersectionSizeSortedUnique
        // 2) the two sketches have different k (throws exception)
        // 3) seed/noncanonical parameters differ
        //
        // This case generates two sketches from the same sequence and parameters; jaccard must be near 1.
        const std::size_t k = 21;
        const std::size_t sketch_size = 2000;
        const int seed = 42;

        const std::string s = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        auto a = mash::sketchFromSequence(s, k, sketch_size, /*noncanonical*/true, seed);
        auto b = mash::sketchFromSequence(s, k, sketch_size, /*noncanonical*/true, seed);

        CHECK(a.k == k);
        CHECK(b.k == k);
        CHECK(std::is_sorted(a.hashes.begin(), a.hashes.end()));
        CHECK(std::is_sorted(b.hashes.begin(), b.hashes.end()));
        CHECK(std::adjacent_find(a.hashes.begin(), a.hashes.end()) == a.hashes.end());
        CHECK(std::adjacent_find(b.hashes.begin(), b.hashes.end()) == b.hashes.end());

        const auto inter = mash::intersectionSizeSortedUnique(a.hashes, b.hashes);
        const double j = mash::jaccard(a, b);
        MESSAGE("a.size=", a.size(), " b.size=", b.size(), " inter=", inter, " jaccard=", j);


        CHECK(j == doctest::Approx(1.0));
    }

    TEST_CASE("Sketch-Sketch jaccard: similar sequences should usually be > 0")
    {
        // This test is meant to catch the "often returns 0" behavior you described.
        // Note: for random sequences with k=21 and sketch_size=2000, the intersection is typically tiny,
        // so returning 0 is actually expected (MinHash intersection = 0 means estimated Jaccard is very low).
        // Therefore we use sequences that are similar but not identical, so they should share k-mers and jaccard should not be 0.

        const std::size_t k = 21;
        const std::size_t sketch_size = 2000;
        const int seed = 11;

        std::string s1 = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        std::string s2 = s1;
        // Modify a small region in the middle while keeping most k-mers the same
        if (s2.size() >= 10) {
            s2[5] = 'T';
            s2[6] = 'T';
            s2[7] = 'T';
        }

        auto a = mash::sketchFromSequence(s1, k, sketch_size, /*noncanonical*/true, seed);
        auto b = mash::sketchFromSequence(s2, k, sketch_size, /*noncanonical*/true, seed);

        const double j = mash::jaccard(a, b);
        MESSAGE("similar sequences jaccard=", j, " a.size=", a.size(), " b.size=", b.size());

        CHECK(j > 0.0);
    }

    TEST_CASE("BloomFilter-Sketch jaccard: exactness when FPP very small")
    {
        // Here we set a very low false positive rate.
        // At this small test scale, BloomFilter queries are almost equivalent to exact set membership.
        const std::size_t k = 21;

        mash::Sketch ref;
        ref.k = k;
        ref.hashes = {1, 2, 3, 100, 200, 500};

        mash::Sketch q1;
        q1.k = k;
        q1.hashes = {2, 3, 4, 5};

        // build a bloom filter (very low false positive rate)
        const double fpp = 1e-9;
        const int seed = 12345;
        auto bf = mash::filterFromSketch(ref, fpp, seed);

        // Since BloomFilter approximates intersection, we only require it to be close to the Sketch-Sketch result.
        const double j_exact = mash::jaccard(ref, q1);
        const double j_bf = mash::jaccard(bf, q1);

        // With a very small fpp, results should be very close (allowing tiny error).
        CHECK(j_bf == doctest::Approx(j_exact).epsilon(1e-6));
    }

    TEST_CASE("BloomFilter-Sketch jaccard: monotonic sanity (superset should not reduce) ")
    {
        // BloomFilter false positives only make contains larger (never turn true into false),
        // so for the same BloomFilter:
        // - if the query adds more elements (while keeping existing ones), estimated intersection should not decrease.
        mash::Sketch ref;
        ref.k = 21;
        ref.hashes = {10, 20, 30, 40, 50};

        auto bf = mash::filterFromSketch(ref, /*fpp*/1e-6, /*seed*/7);

        mash::Sketch q_small;
        q_small.k = 21;
        q_small.hashes = {10, 999};

        mash::Sketch q_large;
        q_large.k = 21;
        q_large.hashes = {10, 20, 999, 888};

        const double j1 = mash::jaccard(bf, q_small);
        const double j2 = mash::jaccard(bf, q_large);

        CHECK(j2 + 1e-12 >= j1);
    }

    TEST_CASE("jaccard_perf: fixed ref vs 10000 queries (sketch size=2k)")
    {
        const char* env = std::getenv("HALIGN4_RUN_PERF");
        if (!env || std::string(env) != "1")
        {
            DOCTEST_INFO("jaccard_perf skipped; set HALIGN4_RUN_PERF=1 to enable");
            return;
        }

        // Fixed params: sketch size 2k
        const std::size_t k = 21;
        const std::size_t sketch_size = 2000;
        const std::size_t seq_len = 300000;       // default 300k
        const std::size_t num_queries = 10000;   // 10k
        const int seed = 42;

        std::mt19937_64 rng(123456);

        // generate ref
        const std::string ref_seq = random_dna(rng, seq_len);
        mash::Sketch ref = mash::sketchFromSequence(ref_seq, k, sketch_size, /*noncanonical*/true, seed);
        REQUIRE(ref.size() == ref.hashes.size());

        // ref's bloom filter (used to accelerate approximate jaccard for query vs ref)
        bloom_filter ref_bf = mash::filterFromSketch(ref, /*fpp*/1e-8, /*seed*/seed);

        // generate query sketches (first generate sequences then sketch, to keep the test closer to real workload)
        std::vector<mash::Sketch> queries;
        queries.reserve(num_queries);
        for (std::size_t i = 0; i < num_queries; ++i)
        {
            std::string qseq = random_dna(rng, seq_len);
            queries.emplace_back(mash::sketchFromSequence(qseq, k, sketch_size, /*noncanonical*/true, seed));
        }

        // --- perf test 1: Sketch-Sketch ---
        {
            volatile double sink = 0.0;
            auto t0 = std::chrono::steady_clock::now();
            for (const auto& q : queries)
            {
                sink += mash::jaccard(ref, q);
            }
            auto t1 = std::chrono::steady_clock::now();
            const double sec = std::chrono::duration<double>(t1 - t0).count();
            MESSAGE("jaccard_perf Sketch-Sketch: queries=" << num_queries << " sketch_size=" << sketch_size << " took " << sec << " s" << " sink=" << sink);
        }

        // --- perf test 2: BloomFilter-Sketch ---
        {
            volatile double sink = 0.0;
            auto t0 = std::chrono::steady_clock::now();
            for (const auto& q : queries)
            {
                sink += mash::jaccard(ref_bf, q);
            }
            auto t1 = std::chrono::steady_clock::now();
            const double sec = std::chrono::duration<double>(t1 - t0).count();
            MESSAGE("jaccard_perf BloomFilter-Sketch: queries=" << num_queries << " sketch_size=" << sketch_size << " took " << sec << " s" << " sink=" << sink);
        }

        CHECK(true);
    }
}
