#include <doctest/doctest.h>
#include <chrono>
#include <random>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include "align.h"
#include "seed.h"  // used for minimizer extraction and anchor generation
#include "mash.h"  // used for computing Mash similarity

// ------------------------------------------------------------------
// Helper: generate random DNA sequence
// ------------------------------------------------------------------
static std::string generateRandomDNA(size_t length, unsigned seed = 42) {
    static const char bases[] = {'A', 'C', 'G', 'T'};
    std::mt19937 rng(seed);
    std::uniform_int_distribution<int> dist(0, 3);

    std::string seq;
    seq.reserve(length);
    for (size_t i = 0; i < length; ++i) {
        seq += bases[dist(rng)];
    }
    return seq;
}

// ------------------------------------------------------------------
// Helper: introduce random mutations (SNP + Indel) into a sequence
// ------------------------------------------------------------------
static std::string mutateSequence(const std::string& ref,
                                  double snp_rate = 0.01,     // 1% SNP
                                  double indel_rate = 0.005,  // 0.5% indel
                                  unsigned seed = 43) {
    static const char bases[] = {'A', 'C', 'G', 'T'};
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> prob(0.0, 1.0);
    std::uniform_int_distribution<int> base_dist(0, 3);
    std::uniform_int_distribution<int> indel_len(1, 5);  // indel length 1-5bp

    std::string mutated;
    mutated.reserve(ref.size() * 1.1);  // reserve capacity

    for (size_t i = 0; i < ref.size(); ++i) {
        double p = prob(rng);

        if (p < indel_rate) {
            // Random insertion or deletion
            if (prob(rng) < 0.5) {
                // Insertion
                int len = indel_len(rng);
                for (int j = 0; j < len; ++j) {
                    mutated += bases[base_dist(rng)];
                }
            } else {
                // Deletion: skip the current base
                continue;
            }
        } else if (p < indel_rate + snp_rate) {
            // SNP: replace with a different base
            char original = ref[i];
            char replacement;
            do {
                replacement = bases[base_dist(rng)];
            } while (replacement == original);
            mutated += replacement;
        } else {
            // Keep unchanged
            mutated += ref[i];
        }
    }

    return mutated;
}

// ------------------------------------------------------------------
// Helper: generate a sequence containing structural variants
// ------------------------------------------------------------------
// Supported SV types:
// - Large insertion (INS)
// - Large deletion (DEL)
// - Inversion (INV)
// - Tandem duplication (DUP)
static std::string generateSVSequence(const std::string& ref,
                                     const std::string& sv_type,
                                     size_t sv_pos,
                                     size_t sv_size,
                                     unsigned seed = 44) {
    std::string result;
    result.reserve(ref.size() + sv_size);

    if (sv_type == "INS") {
        // Insertion: insert random sequence at the specified position
        result = ref.substr(0, sv_pos);
        result += generateRandomDNA(sv_size, seed);
        result += ref.substr(sv_pos);
    }
    else if (sv_type == "DEL") {
        // Deletion: remove segment starting at the specified position
        result = ref.substr(0, sv_pos);
        if (sv_pos + sv_size < ref.size()) {
            result += ref.substr(sv_pos + sv_size);
        }
    }
    else if (sv_type == "INV") {
        // Inversion: reverse the specified region
        result = ref.substr(0, sv_pos);
        std::string inv_region = ref.substr(sv_pos, sv_size);
        std::reverse(inv_region.begin(), inv_region.end());
        // Reverse complement
        for (char& c : inv_region) {
            switch (c) {
                case 'A': c = 'T'; break;
                case 'T': c = 'A'; break;
                case 'C': c = 'G'; break;
                case 'G': c = 'C'; break;
            }
        }
        result += inv_region;
        if (sv_pos + sv_size < ref.size()) {
            result += ref.substr(sv_pos + sv_size);
        }
    }
    else if (sv_type == "DUP") {
        // Tandem duplication: duplicate the specified region and insert it after
        result = ref.substr(0, sv_pos + sv_size);
        result += ref.substr(sv_pos, sv_size);  // duplicate once
        if (sv_pos + sv_size < ref.size()) {
            result += ref.substr(sv_pos + sv_size);
        }
    }
    else {
        // Unknown type, return the original sequence
        result = ref;
    }

    return result;
}

// ------------------------------------------------------------------
// Helper: generate a complex sequence containing multiple SVs
// ------------------------------------------------------------------
struct SVEvent {
    std::string type;  // INS, DEL, INV, DUP
    size_t pos;        // position
    size_t size;       // size
};

static std::string generateComplexSVSequence(const std::string& ref,
                                            const std::vector<SVEvent>& events,
                                            unsigned seed = 45) {
    std::string result = ref;

    // Apply SVs from back to front to avoid index shift issues
    std::vector<SVEvent> sorted_events = events;
    std::sort(sorted_events.begin(), sorted_events.end(),
              [](const SVEvent& a, const SVEvent& b) { return a.pos > b.pos; });

    for (const auto& sv : sorted_events) {
        result = generateSVSequence(result, sv.type, sv.pos, sv.size, seed++);
    }

    return result;
}

// ------------------------------------------------------------------
// Helper: generate real anchors (based on minimizer matches)
// ------------------------------------------------------------------
// Note:
// Previous test cases used fixed-position anchors (assuming ref/query are aligned at the same positions),
// but after mutateSequence, due to indels, those anchors become completely inaccurate.
//
// This function uses minimizers to extract real shared k-mers and generate accurate anchors:
// 1) Extract minimizers from ref and query
// 2) Use collect_anchors to find matching minimizer hits
// 3) Return the list of real anchors
//
// Parameters:
// @param ref - reference sequence
// @param query - query sequence
// @param k - k-mer size (default 15)
// @param w - window size (default 10)
// @return real anchors list
// ------------------------------------------------------------------
static anchor::Anchors generateRealAnchors(const std::string& ref,
                                           const std::string& query,
                                           std::size_t k = 15,
                                           std::size_t w = 10) {
    // 1) Extract minimizers of ref
    minimizer::MinimizerHits ref_hits = minimizer::extractMinimizer(ref, k, w, false);

    // 2) Extract minimizers of query
    minimizer::MinimizerHits qry_hits = minimizer::extractMinimizer(query, k, w, false);

    // 3) Collect anchors (find matching minimizers)
    anchor::Anchors anchors = minimizer::collect_anchors(ref_hits, qry_hits);

    return anchors;
}

// ------------------------------------------------------------------
// Helper function: CIGAR to string (for debugging purposes)
// ------------------------------------------------------------------
static std::string cigarToString(const cigar::Cigar_t& cigar) {
    std::string result;
    for (auto unit : cigar) {
        char op;
        uint32_t len;
        cigar::intToCigar(unit, op, len);
        result += std::to_string(len) + op;
    }
    return result;
}

// ------------------------------------------------------------------
// Performance test helper: timer
// ------------------------------------------------------------------
struct Timer {
    std::chrono::high_resolution_clock::time_point start;

    Timer() : start(std::chrono::high_resolution_clock::now()) {}

    double elapsedMs() const {
        auto end = std::chrono::high_resolution_clock::now();
        return std::chrono::duration<double, std::milli>(end - start).count();
    }
};

// ------------------------------------------------------------------
// Test suite: correctness tests
// ------------------------------------------------------------------
TEST_SUITE("align") {

    TEST_CASE("globalAlignKSW2 - Exact match") {
        std::string seq = "ACGTACGTACGT";
        auto cigar = align::globalAlignKSW2(seq, seq);

        // Should be a perfect match; CIGAR should be 12M
        REQUIRE(cigar.size() >= 1);
        char op;
        uint32_t len;
        cigar::intToCigar(cigar[0], op, len);
        CHECK(op == 'M');
        CHECK(len == 12);
    }

    TEST_CASE("globalAlignKSW2 - Single mismatch") {
        std::string ref   = "ACGTACGTACGT";
        std::string query = "ACGTACCGTACGT";  // 6th position T->C mismatch

        auto cigar = align::globalAlignKSW2(ref, query);
        REQUIRE(cigar.size() > 0);

        // Verify CIGAR is not empty
        std::string cigar_str = cigarToString(cigar);
        MESSAGE("CIGAR: ", cigar_str);
        CHECK(!cigar_str.empty());
    }

    TEST_CASE("globalAlignKSW2 - Single insertion") {
        std::string ref   = "ACGTACGTACGT";
        std::string query = "ACGTAACGTACGT";  // Insert A after the 5th position

        auto cigar = align::globalAlignKSW2(ref, query);
        std::string cigar_str = cigarToString(cigar);
        MESSAGE("CIGAR: ", cigar_str);

        // Should contain insertion operations (I)
        bool has_insertion = cigar_str.find('I') != std::string::npos;
        CHECK(has_insertion);
    }

    TEST_CASE("extendAlignKSW2 - Basic extension") {
        std::string ref = generateRandomDNA(1000, 100);
        std::string query = ref.substr(100, 500);  // Extract the middle fragment

        auto cigar = align::extendAlignKSW2(ref, query, 200);
        REQUIRE(cigar.size() > 0);

        std::string cigar_str = cigarToString(cigar);
        MESSAGE("Extend CIGAR: ", cigar_str);
        CHECK(!cigar_str.empty());
    }

    TEST_CASE("globalAlignWFA2 - Exact match") {
        std::string seq = "ACGTACGTACGT";
        auto cigar = align::globalAlignWFA2(seq, seq);

        REQUIRE(cigar.size() >= 1);
        std::string cigar_str = cigarToString(cigar);
        MESSAGE("WFA2 CIGAR: ", cigar_str);
        CHECK(!cigar_str.empty());
    }

    TEST_CASE("Empty sequence boundary test") {
        std::string empty = "";
        std::string seq = "ACGT";

        // Verify empty sequences do not crash (behavior depends on implementation)
        CHECK_NOTHROW(align::globalAlignKSW2(empty, seq));
        CHECK_NOTHROW(align::globalAlignKSW2(seq, empty));
    }

    TEST_CASE("High-similarity sequence - 99% similarity") {
        // Generate a 1000bp sequence with 1% mismatches
        std::string ref = generateRandomDNA(1000, 12345);
        std::string query = mutateSequence(ref, 0.01, 0.0, 12346);  // only SNPs, no indels

        // All three methods should align correctly
        auto cigar_ksw2 = align::globalAlignKSW2(ref, query);
        auto cigar_extend = align::extendAlignKSW2(ref, query, 200);
        auto cigar_wfa2 = align::globalAlignWFA2(ref, query);

        // Verify all CIGARs are non-empty
        CHECK(cigar_ksw2.size() > 0);
        CHECK(cigar_extend.size() > 0);
        CHECK(cigar_wfa2.size() > 0);

        MESSAGE("KSW2 CIGAR: ", cigarToString(cigar_ksw2));
        MESSAGE("Extend CIGAR: ", cigarToString(cigar_extend));
        MESSAGE("WFA2 CIGAR: ", cigarToString(cigar_wfa2));
    }

    TEST_CASE("High-similarity sequence - 98% similarity (with indel)") {
        // Generate a sequence with 1% SNPs and 1% indels
        std::string ref = generateRandomDNA(500, 54321);
        std::string query = mutateSequence(ref, 0.01, 0.01, 54322);

        auto cigar_ksw2 = align::globalAlignKSW2(ref, query);
        auto cigar_wfa2 = align::globalAlignWFA2(ref, query);

        // Verify CIGAR contains different operation types
        std::string cigar_str_ksw2 = cigarToString(cigar_ksw2);
        std::string cigar_str_wfa2 = cigarToString(cigar_wfa2);

        MESSAGE("KSW2 CIGAR (with indels): ", cigar_str_ksw2);
        MESSAGE("WFA2 CIGAR (with indels): ", cigar_str_wfa2);

        CHECK(!cigar_str_ksw2.empty());
        CHECK(!cigar_str_wfa2.empty());

        // Should include match operations
        bool has_match_ksw2 = cigar_str_ksw2.find('M') != std::string::npos;
        bool has_match_wfa2 = cigar_str_wfa2.find('M') != std::string::npos;
        CHECK(has_match_ksw2);
        CHECK(has_match_wfa2);
    }

    TEST_CASE("Ultra-high similarity sequence - 99.9% similarity") {
        // Simulate sequencing errors: only 0.1% error rate
        std::string ref = generateRandomDNA(10000, 99999);
        std::string query = mutateSequence(ref, 0.0005, 0.0005, 100000);

        // At such high similarity, all methods should complete quickly
        Timer timer;
        auto cigar = align::globalAlignKSW2(ref, query);
        double elapsed = timer.elapsedMs();

        CHECK(cigar.size() > 0);
        MESSAGE("Ultra-high similarity 10k alignment elapsed: ", elapsed, " ms");

        // For 10k sequences with 99.9% similarity, it should finish in a reasonable time (<100ms)
        CHECK(elapsed < 100.0);
    }

    // ------------------------------------------------------------------
    // Test: cigarToString and stringToCigar are inverses
    // ------------------------------------------------------------------
    TEST_CASE("cigar::cigarToString and stringToCigar - Inverse operations") {
        // Test 1: standard CIGAR string
        SUBCASE("Standard CIGAR") {
            cigar::Cigar_t original;
            original.push_back(cigar::cigarToInt('M', 100));
            original.push_back(cigar::cigarToInt('I', 5));
            original.push_back(cigar::cigarToInt('M', 95));
            original.push_back(cigar::cigarToInt('D', 3));
            original.push_back(cigar::cigarToInt('M', 50));

            // Convert to string
            std::string cigar_str = cigar::cigarToString(original);
            CHECK(cigar_str == "100M5I95M3D50M");

            // Convert back to Cigar_t
            cigar::Cigar_t roundtrip = cigar::stringToCigar(cigar_str);

            // Verify inverse property
            REQUIRE(roundtrip.size() == original.size());
            for (size_t i = 0; i < original.size(); ++i) {
                CHECK(roundtrip[i] == original[i]);
            }
        }

        // Test 2: all CIGAR operators
        SUBCASE("All operators") {
            cigar::Cigar_t original;
            original.push_back(cigar::cigarToInt('M', 10));
            original.push_back(cigar::cigarToInt('I', 2));
            original.push_back(cigar::cigarToInt('D', 3));
            original.push_back(cigar::cigarToInt('N', 100));
            original.push_back(cigar::cigarToInt('S', 5));
            original.push_back(cigar::cigarToInt('H', 10));
            original.push_back(cigar::cigarToInt('P', 1));
            original.push_back(cigar::cigarToInt('=', 20));
            original.push_back(cigar::cigarToInt('X', 3));

            std::string cigar_str = cigar::cigarToString(original);
            cigar::Cigar_t roundtrip = cigar::stringToCigar(cigar_str);

            REQUIRE(roundtrip.size() == original.size());
            for (size_t i = 0; i < original.size(); ++i) {
                CHECK(roundtrip[i] == original[i]);
            }
        }

        // Test 3: special value "*"
        SUBCASE("Special value *") {
            cigar::Cigar_t empty_cigar = cigar::stringToCigar("*");
            CHECK(empty_cigar.empty());
        }

        // Test 4: empty string
        SUBCASE("Empty string") {
            cigar::Cigar_t empty_cigar = cigar::stringToCigar("");
            CHECK(empty_cigar.empty());
        }

        // Test 5: large numeric lengths
        SUBCASE("Large numeric lengths") {
            cigar::Cigar_t original;
            original.push_back(cigar::cigarToInt('M', 999999));
            original.push_back(cigar::cigarToInt('D', 123456));

            std::string cigar_str = cigar::cigarToString(original);
            CHECK(cigar_str == "999999M123456D");

            cigar::Cigar_t roundtrip = cigar::stringToCigar(cigar_str);
            REQUIRE(roundtrip.size() == original.size());
            CHECK(roundtrip[0] == original[0]);
            CHECK(roundtrip[1] == original[1]);
        }
    }

    // ------------------------------------------------------------------
// Test: stringToCigar error handling
// ------------------------------------------------------------------
TEST_CASE("cigar::stringToCigar - Error handling") {
    // Test 1: no number before operator
    SUBCASE("No number before operator") {
        CHECK_THROWS_AS(cigar::stringToCigar("M10"), std::runtime_error);
    }

    // Test 2: unknown operator
    SUBCASE("Unknown operator") {
        CHECK_THROWS_AS(cigar::stringToCigar("10Q"), std::runtime_error);
    }

    // Test 3: trailing number without operator
    SUBCASE("Trailing number without operator") {
        CHECK_THROWS_AS(cigar::stringToCigar("10M5"), std::runtime_error);
    }

    // Test 4: length is 0
            CHECK_THROWS_AS(cigar::stringToCigar("0M"), std::runtime_error);
        }

    // ------------------------------------------------------------------
// Test: stringToCigar robustness
// ------------------------------------------------------------------
TEST_CASE("cigar::stringToCigar - Robustness") {
    // Test 1: contains whitespace
            std::string cigar_with_spaces = " 10M 5I  3D ";
            cigar::Cigar_t result = cigar::stringToCigar(cigar_with_spaces);

            REQUIRE(result.size() == 3);
            char op;
            uint32_t len;

            cigar::intToCigar(result[0], op, len);
            CHECK(op == 'M');
            CHECK(len == 10);

            cigar::intToCigar(result[1], op, len);
            CHECK(op == 'I');
            CHECK(len == 5);

            cigar::intToCigar(result[2], op, len);
            CHECK(op == 'D');
            CHECK(len == 3);
        }
    }

// ------------------------------------------------------------------
// Performance test suite
// ------------------------------------------------------------------
TEST_SUITE("align_perf") {

    // ------------------------------------------------------------------
    // Performance test: short sequences (~100bp)
    // ------------------------------------------------------------------
    TEST_CASE("Performance - Short sequences (~100bp)") {
        constexpr int NUM_RUNS = 1000;
        constexpr size_t SEQ_LEN = 100;

        std::cout << "\n========== Short-sequence performance test (100bp, " << NUM_RUNS << " runs) ==========\n";

        // Generate test data
        std::vector<std::pair<std::string, std::string>> test_pairs;
        std::vector<anchor::Anchors> anchors_list;

        for (int i = 0; i < NUM_RUNS; ++i) {
            std::string ref = generateRandomDNA(SEQ_LEN, i * 2);
            std::string query = mutateSequence(ref, 0.02, 0.01, i * 2 + 1);
            test_pairs.emplace_back(ref, query);

            // Generate simulated anchors for MM2 (one every 30bp)
            anchor::Anchors anchors;
            for (size_t pos = 0; pos + 15 < SEQ_LEN; pos += 30) {
                anchor::Anchor a;
                a.hash = pos * 1000 + i;
                a.rid_ref = 0;
                a.pos_ref = static_cast<uint32_t>(pos);
                a.rid_qry = 0;
                a.pos_qry = static_cast<uint32_t>(pos);
                a.span = 15;
                a.is_rev = false;
                anchors.push_back(a);
            }
            anchors_list.push_back(anchors);
        }

        // Test globalAlignKSW2
        {
            Timer timer;
            for (const auto& [ref, query] : test_pairs) {
                auto cigar = align::globalAlignKSW2(ref, query);
                (void)cigar;  // prevent optimization removing it
            }
            double elapsed = timer.elapsedMs();
            std::cout << "  globalAlignKSW2:  " << std::fixed << std::setprecision(2)
                      << elapsed << " ms (" << (elapsed / NUM_RUNS) << " ms/run)\n";
        }

        // Test extendAlignKSW2
        {
            Timer timer;
            for (const auto& [ref, query] : test_pairs) {
                auto cigar = align::extendAlignKSW2(ref, query, 200);
                (void)cigar;
            }
            double elapsed = timer.elapsedMs();
            std::cout << "  extendAlignKSW2:  " << std::fixed << std::setprecision(2)
                      << elapsed << " ms (" << (elapsed / NUM_RUNS) << " ms/run)\n";
        }

        // Test globalAlignWFA2
        {
            Timer timer;
            for (const auto& [ref, query] : test_pairs) {
                auto cigar = align::globalAlignWFA2(ref, query);
                (void)cigar;
            }
            double elapsed = timer.elapsedMs();
            std::cout << "  globalAlignWFA2:  " << std::fixed << std::setprecision(2)
                      << elapsed << " ms (" << (elapsed / NUM_RUNS) << " ms/run)\n";
        }

        // Test globalAlignMM2 (with anchors)
        {
            Timer timer;
            for (size_t i = 0; i < test_pairs.size(); ++i) {
                const auto& [ref, query] = test_pairs[i];
                auto cigar = align::globalAlignMM2(ref, query, anchors_list[i]);
                (void)cigar;
            }
            double elapsed = timer.elapsedMs();
            std::cout << "  globalAlignMM2:   " << std::fixed << std::setprecision(2)
                      << elapsed << " ms (" << (elapsed / NUM_RUNS) << " ms/run)\n";
        }

        std::cout << "========================================================\n\n";
    }

    // ------------------------------------------------------------------
    // Performance test: medium sequences (~1000bp)
    // ------------------------------------------------------------------
    TEST_CASE("Performance - Medium sequences (~1000bp)") {
        constexpr int NUM_RUNS = 100;
        constexpr size_t SEQ_LEN = 1000;

        std::cout << "\n========== Medium-sequence performance test (1000bp, " << NUM_RUNS << " runs) ==========\n";

        std::vector<std::pair<std::string, std::string>> test_pairs;
        std::vector<anchor::Anchors> anchors_list;

        for (int i = 0; i < NUM_RUNS; ++i) {
            std::string ref = generateRandomDNA(SEQ_LEN, i * 2);
            std::string query = mutateSequence(ref, 0.02, 0.01, i * 2 + 1);
            test_pairs.emplace_back(ref, query);

            // Generate simulated anchors for MM2 (one every 150bp)
            anchor::Anchors anchors;
            for (size_t pos = 0; pos + 20 < SEQ_LEN; pos += 150) {
                anchor::Anchor a;
                a.hash = pos * 1000 + i;
                a.rid_ref = 0;
                a.pos_ref = static_cast<uint32_t>(pos);
                a.rid_qry = 0;
                a.pos_qry = static_cast<uint32_t>(pos);
                a.span = 20;
                a.is_rev = false;
                anchors.push_back(a);
            }
            anchors_list.push_back(anchors);
        }

        {
            Timer timer;
            for (const auto& [ref, query] : test_pairs) {
                auto cigar = align::globalAlignKSW2(ref, query);
                (void)cigar;
            }
            double elapsed = timer.elapsedMs();
            std::cout << "  globalAlignKSW2:  " << std::fixed << std::setprecision(2)
                      << elapsed << " ms (" << (elapsed / NUM_RUNS) << " ms/run)\n";
        }

        {
            Timer timer;
            for (const auto& [ref, query] : test_pairs) {
                auto cigar = align::extendAlignKSW2(ref, query, 200);
                (void)cigar;
            }
            double elapsed = timer.elapsedMs();
            std::cout << "  extendAlignKSW2:  " << std::fixed << std::setprecision(2)
                      << elapsed << " ms (" << (elapsed / NUM_RUNS) << " ms/run)\n";
        }

        {
            Timer timer;
            for (const auto& [ref, query] : test_pairs) {
                auto cigar = align::globalAlignWFA2(ref, query);
                (void)cigar;
            }
            double elapsed = timer.elapsedMs();
            std::cout << "  globalAlignWFA2:  " << std::fixed << std::setprecision(2)
                      << elapsed << " ms (" << (elapsed / NUM_RUNS) << " ms/run)\n";
        }

        {
            Timer timer;
            for (size_t i = 0; i < test_pairs.size(); ++i) {
                const auto& [ref, query] = test_pairs[i];
                auto cigar = align::globalAlignMM2(ref, query, anchors_list[i]);
                (void)cigar;
            }
            double elapsed = timer.elapsedMs();
            std::cout << "  globalAlignMM2:   " << std::fixed << std::setprecision(2)
                      << elapsed << " ms (" << (elapsed / NUM_RUNS) << " ms/run)\n";
        }

        std::cout << "========================================================\n\n";
    }

    // ------------------------------------------------------------------
    // Performance test: long sequences (~10000bp)
    // ------------------------------------------------------------------
    TEST_CASE("Performance - Long sequences (~10000bp)") {
        constexpr int NUM_RUNS = 10;
        constexpr size_t SEQ_LEN = 10000;

        std::cout << "\n========== Long-sequence performance test (10000bp, " << NUM_RUNS << " runs) ==========\n";

        std::vector<std::pair<std::string, std::string>> test_pairs;
        std::vector<anchor::Anchors> anchors_list;

        for (int i = 0; i < NUM_RUNS; ++i) {
            std::string ref = generateRandomDNA(SEQ_LEN, i * 2);
            std::string query = mutateSequence(ref, 0.02, 0.01, i * 2 + 1);
            test_pairs.emplace_back(ref, query);

            // Generate simulated anchors for MM2 (one every 500bp)
            anchor::Anchors anchors;
            for (size_t pos = 0; pos + 50 < SEQ_LEN; pos += 500) {
                anchor::Anchor a;
                a.hash = pos * 1000 + i;
                a.rid_ref = 0;
                a.pos_ref = static_cast<uint32_t>(pos);
                a.rid_qry = 0;
                a.pos_qry = static_cast<uint32_t>(pos);
                a.span = 50;
                a.is_rev = false;
                anchors.push_back(a);
            }
            anchors_list.push_back(anchors);
        }

        {
            Timer timer;
            for (const auto& [ref, query] : test_pairs) {
                auto cigar = align::globalAlignKSW2(ref, query);
                (void)cigar;
            }
            double elapsed = timer.elapsedMs();
            std::cout << "  globalAlignKSW2:  " << std::fixed << std::setprecision(2)
                      << elapsed << " ms (" << (elapsed / NUM_RUNS) << " ms/run)\n";
        }

        {
            Timer timer;
            for (const auto& [ref, query] : test_pairs) {
                auto cigar = align::extendAlignKSW2(ref, query, 200);
                (void)cigar;
            }
            double elapsed = timer.elapsedMs();
            std::cout << "  extendAlignKSW2:  " << std::fixed << std::setprecision(2)
                      << elapsed << " ms (" << (elapsed / NUM_RUNS) << " ms/run)\n";
        }

        {
            Timer timer;
            for (const auto& [ref, query] : test_pairs) {
                auto cigar = align::globalAlignWFA2(ref, query);
                (void)cigar;
            }
            double elapsed = timer.elapsedMs();
            std::cout << "  globalAlignWFA2:  " << std::fixed << std::setprecision(2)
                      << elapsed << " ms (" << (elapsed / NUM_RUNS) << " ms/run)\n";
        }

        {
            Timer timer;
            for (size_t i = 0; i < test_pairs.size(); ++i) {
                const auto& [ref, query] = test_pairs[i];
                auto cigar = align::globalAlignMM2(ref, query, anchors_list[i]);
                (void)cigar;
            }
            double elapsed = timer.elapsedMs();
            std::cout << "  globalAlignMM2:   " << std::fixed << std::setprecision(2)
                      << elapsed << " ms (" << (elapsed / NUM_RUNS) << " ms/run)\n";
        }

        std::cout << "========================================================\n\n";
    }


    // ------------------------------------------------------------------
    // Performance test: high similarity sequences (95%-99.9%, real sequencing scenario)
    // ------------------------------------------------------------------
    TEST_CASE("Performance - High similarity sequences (real-world)") {
        constexpr int NUM_RUNS = 200;
        constexpr size_t SEQ_LEN = 1000;

        std::cout << "\n========== High-similarity sequence performance test (1000bp, " << NUM_RUNS << " runs) ==========\n";
        std::cout << "Note: simulates real genome sequencing scenarios, sequence similarity >95%\n\n";

        // Test different similarity levels (by controlling SNP + indel rates)
        struct SimilarityLevel {
            double snp_rate;
            double indel_rate;
            const char* desc;
            double similarity;  // expected similarity
        };

        std::vector<SimilarityLevel> levels = {
            {0.0001, 0.0001, "Ultra-high similarity (>99.9%)", 99.98},
            {0.001,  0.0005, "Very high similarity (~99%)",   99.0},
            {0.005,  0.002,  "High similarity (~98%)",     98.0},
            {0.01,   0.005,  "Moderate similarity (~97%)",   97.0},
            {0.03,   0.01,   "Lower similarity (~95%)",   95.0},
            {0.05,   0.02,   "Low similarity (~90%)",     90.0},
            {0.08,   0.03,   "Even lower similarity (~85%)",   85.0},
            {0.12,   0.05,   "Very low similarity (~80%)",   80.0},
            {0.18,   0.07,   "Extremely low similarity (~70%)",   70.0}
        };

        for (const auto& level : levels) {
            std::cout << "---------- " << level.desc << " ----------\n";

            // Generate test data
            std::vector<std::pair<std::string, std::string>> test_pairs;
            std::vector<anchor::Anchors> anchors_list;

            for (int i = 0; i < NUM_RUNS; ++i) {
                std::string ref = generateRandomDNA(SEQ_LEN, i * 2);
                std::string query = mutateSequence(ref, level.snp_rate, level.indel_rate, i * 2 + 1);
                test_pairs.emplace_back(ref, query);

                // Adjust anchor point generation strategy based on similarity
                // High similarity: dense anchor points (per 150bp)
                // Low similarity: sparse anchor points (per 300bp) because there are fewer reliable anchor points.
                size_t anchor_interval = (level.similarity >= 90.0) ? 150 : 300;
                size_t anchor_span = (level.similarity >= 90.0) ? 20 : 15;

                anchor::Anchors anchors;
                for (size_t pos = 0; pos + anchor_span < SEQ_LEN; pos += anchor_interval) {
                    anchor::Anchor a;
                    a.hash = pos * 1000 + i;
                    a.rid_ref = 0;
                    a.pos_ref = static_cast<uint32_t>(pos);
                    a.rid_qry = 0;
                    a.pos_qry = static_cast<uint32_t>(pos);
                    a.span = static_cast<uint32_t>(anchor_span);
                    a.is_rev = false;
                    anchors.push_back(a);
                }
                anchors_list.push_back(anchors);
            }

            // ========== Mash similarity calculation ==========
            // Calculate the actual Jaccard similarity and ANI using Mash sketch.
            {
                constexpr std::size_t MASH_K = 21;          // k-mer size
                constexpr std::size_t MASH_SKETCH_SIZE = 2000; // sketch size

                double total_jaccard = 0.0;
                double total_ani = 0.0;
                int valid_count = 0;

                // Calculate Mash similarity for the first 20 sequence pairs (to avoid excessive computational overhead).
                int mash_sample_size = std::min(20, NUM_RUNS);
                for (int i = 0; i < mash_sample_size; ++i) {
                    const auto& [ref, query] = test_pairs[i];

                    // Generate sketch
                    auto sketch_ref = mash::sketchFromSequence(ref, MASH_K, MASH_SKETCH_SIZE);
                    auto sketch_query = mash::sketchFromSequence(query, MASH_K, MASH_SKETCH_SIZE);

                    // Calculate Jaccard similarity
                    double j = mash::jaccard(sketch_ref, sketch_query);

                    // Calculate ANI (Average Nucleotide Identity)
                    double ani = mash::aniFromJaccard(j, MASH_K);

                    total_jaccard += j;
                    total_ani += ani;
                    valid_count++;
                }

                if (valid_count > 0) {
                    double avg_jaccard = total_jaccard / valid_count;
                    double avg_ani = total_ani / valid_count;

                    std::cout << "  Mash similarity (k=" << MASH_K << ", s=" << MASH_SKETCH_SIZE << "):\n";
                    std::cout << "    Jaccard: " << std::fixed << std::setprecision(4) << avg_jaccard
                              << "  ANI: " << std::setprecision(2) << (avg_ani * 100.0) << "%\n";
                }
            }

            // ========== Alignment performance test ==========

            // KSW2 test
            {
                Timer timer;
                for (const auto& [ref, query] : test_pairs) {
                    auto cigar = align::globalAlignKSW2(ref, query);
                    (void)cigar;
                }
                double elapsed = timer.elapsedMs();
                double avg_ms = elapsed / NUM_RUNS;
                std::cout << "  KSW2 global alignment:  " << std::fixed << std::setprecision(3)
                          << avg_ms << " ms/run  (Throughput: " << (1000.0 / avg_ms)
                          << " runs/s)\n";
            }

            // KSW2 extend mode test
            {
                Timer timer;
                for (const auto& [ref, query] : test_pairs) {
                    auto cigar = align::extendAlignKSW2(ref, query, 200);
                    (void)cigar;
                }
                double elapsed = timer.elapsedMs();
                double avg_ms = elapsed / NUM_RUNS;
                std::cout << "  KSW2 extension mode:  " << std::fixed << std::setprecision(3)
                          << avg_ms << " ms/run  (Throughput: " << (1000.0 / avg_ms)
                          << " runs/s)\n";
            }

            // WFA2 test
            {
                Timer timer;
                for (const auto& [ref, query] : test_pairs) {
                    // Use the expected similarity (defined in SimilarityLevel).
                    auto cigar = align::globalAlignWFA2(ref, query);
                    (void)cigar;
                }
                double elapsed = timer.elapsedMs();
                double avg_ms = elapsed / NUM_RUNS;
                std::cout << "  WFA2 global alignment:  " << std::fixed << std::setprecision(3)
                          << avg_ms << " ms/run  (Throughput: " << (1000.0 / avg_ms)
                          << " runs/s)\n";
            }

            // MM2 test
            {
                Timer timer;
                for (size_t i = 0; i < test_pairs.size(); ++i) {
                    const auto& [ref, query] = test_pairs[i];
                    auto cigar = align::globalAlignMM2(ref, query, anchors_list[i]);
                    (void)cigar;
                }
                double elapsed = timer.elapsedMs();
                double avg_ms = elapsed / NUM_RUNS;
                std::cout << "  MM2 anchor-based alignment:   " << std::fixed << std::setprecision(3)
                          << avg_ms << " ms/run  (Throughput: " << (1000.0 / avg_ms)
                          << " runs/s)\n";
            }

            std::cout << "\n";
        }

        std::cout << "========================================================\n\n";
    }

    // ------------------------------------------------------------------
    // Performance test: length sensitivity (high similarity, varying lengths)
    // ------------------------------------------------------------------
    TEST_CASE("Performance - Length scaling (high similarity)") {
        constexpr int NUM_RUNS = 50;
        constexpr double SNP_RATE = 0.01;    // 1% SNP (~98% similarity)
        constexpr double INDEL_RATE = 0.005; // 0.5% indel

        std::cout << "\n========== Length scalability test (similarity ~98%, " << NUM_RUNS << " runs) ==========\n";

        std::vector<size_t> lengths = {100, 500, 1000, 5000, 10000};

        std::cout << std::setw(10) << "Length (bp)"
                  << std::setw(15) << "KSW2(ms)"
                  << std::setw(15) << "Extend(ms)"
                  << std::setw(15) << "WFA2(ms)"
                  << std::setw(15) << "MM2(ms)" << "\n";
        std::cout << std::string(70, '-') << "\n";

        for (size_t len : lengths) {
            // Generate test data
            std::vector<std::pair<std::string, std::string>> test_pairs;
            std::vector<anchor::Anchors> anchors_list;

            for (int i = 0; i < NUM_RUNS; ++i) {
                std::string ref = generateRandomDNA(len, i * 2);
                std::string query = mutateSequence(ref, SNP_RATE, INDEL_RATE, i * 2 + 1);
                test_pairs.emplace_back(ref, query);

                // Generate anchors (adjust density based on length)
                size_t anchor_interval = std::max(size_t(50), len / 10);  // about 10 anchors
                anchor::Anchors anchors;
                for (size_t pos = 0; pos + 20 < len; pos += anchor_interval) {
                    anchor::Anchor a;
                    a.hash = pos * 1000 + i;
                    a.rid_ref = 0;
                    a.pos_ref = static_cast<uint32_t>(pos);
                    a.rid_qry = 0;
                    a.pos_qry = static_cast<uint32_t>(pos);
                    a.span = 20;
                    a.is_rev = false;
                    anchors.push_back(a);
                }
                anchors_list.push_back(anchors);
            }

            double ksw2_time = 0.0;
            double extend_time = 0.0;
            double wfa2_time = 0.0;
            double mm2_time = 0.0;

            // KSW2
            {
                Timer timer;
                for (const auto& [ref, query] : test_pairs) {
                    auto cigar = align::globalAlignKSW2(ref, query);
                    (void)cigar;
                }
                ksw2_time = timer.elapsedMs() / NUM_RUNS;
            }

            // Extend
            {
                Timer timer;
                for (const auto& [ref, query] : test_pairs) {
                    auto cigar = align::extendAlignKSW2(ref, query, 200);
                    (void)cigar;
                }
                extend_time = timer.elapsedMs() / NUM_RUNS;
            }

            // WFA2
            {
                Timer timer;
                for (const auto& [ref, query] : test_pairs) {
                    auto cigar = align::globalAlignWFA2(ref, query);
                    (void)cigar;
                }
                wfa2_time = timer.elapsedMs() / NUM_RUNS;
            }

            // MM2
            {
                Timer timer;
                for (size_t i = 0; i < test_pairs.size(); ++i) {
                    const auto& [ref, query] = test_pairs[i];
                    auto cigar = align::globalAlignMM2(ref, query, anchors_list[i]);
                    (void)cigar;
                }
                mm2_time = timer.elapsedMs() / NUM_RUNS;
            }

            std::cout << std::setw(10) << len
                      << std::setw(15) << std::fixed << std::setprecision(3) << ksw2_time
                      << std::setw(15) << extend_time
                      << std::setw(15) << wfa2_time
                      << std::setw(15) << mm2_time << "\n";
        }

        std::cout << "========================================================\n\n";
    }

    // ------------------------------------------------------------------
    // Performance test: length difference scalability (fixed ref length, varying query length)
    // ------------------------------------------------------------------
    TEST_CASE("Performance - Length difference scalability (30k ref vs varying query)") {
        constexpr int NUM_RUNS = 30;
        constexpr size_t REF_LEN = 30000;  // fixed ref length 30k bp
        constexpr double SIMILARITY = 0.95; // 95% similarity
        constexpr double SNP_RATE = 0.03;   // 3% SNP
        constexpr double INDEL_RATE = 0.02; // 2% indel

        std::cout << "\n========== Length difference scalability test (ref=" << REF_LEN << "bp, similarity ~95%, " << NUM_RUNS << " runs) ==========\n";
        std::cout << "Note: simulate alignment performance between a long ref and queries of varying length (asymmetric length scenario)\n\n";

        std::vector<size_t> query_lengths = {100, 500, 1000, 5000, 10000, 15000, 20000, 25000};

        std::cout << std::setw(12) << "Query length"
                  << std::setw(15) << "KSW2(ms)"
                  << std::setw(15) << "Extend(ms)"
                  << std::setw(15) << "WFA2(ms)"
                  << std::setw(15) << "MM2(ms)" << "\n";
        std::cout << std::string(72, '-') << "\n";

        for (size_t query_len : query_lengths) {
            // Generate test data
            std::vector<std::pair<std::string, std::string>> test_pairs;
            std::vector<anchor::Anchors> anchors_list;

            for (int i = 0; i < NUM_RUNS; ++i) {
                // Generate a fixed 30k ref
                std::string ref = generateRandomDNA(REF_LEN, i * 7);

                // Select a contiguous segment from ref as the template for query (simulate real alignment)
                // Choose a position slightly before the middle to ensure enough space
                size_t start_pos = (REF_LEN - query_len) / 3;
                std::string ref_segment = ref.substr(start_pos, query_len);

                // Apply mutations to this segment to obtain the query (95% similarity)
                std::string query = mutateSequence(ref_segment, SNP_RATE, INDEL_RATE, i * 7 + 1);

                test_pairs.emplace_back(ref, query);

                // Generate anchors: in corresponding regions of ref and query
                // Adjust anchor interval based on query length
                size_t anchor_interval = std::max(size_t(100), query_len / 8);
                anchor::Anchors anchors;

                for (size_t offset = 0; offset + 20 < query_len; offset += anchor_interval) {
                    anchor::Anchor a;
                    a.hash = (start_pos + offset) * 1000 + i;
                    a.rid_ref = 0;
                    a.pos_ref = static_cast<uint32_t>(start_pos + offset);  // actual position on ref
                    a.rid_qry = 0;
                    a.pos_qry = static_cast<uint32_t>(offset);              // corresponding position on query
                    a.span = 20;
                    a.is_rev = false;
                    anchors.push_back(a);
                }
                anchors_list.push_back(anchors);
            }

            double ksw2_time = 0.0;
            double extend_time = 0.0;
            double wfa2_time = 0.0;
            double mm2_time = 0.0;

            // KSW2 global alignment
            {
                Timer timer;
                for (const auto& [ref, query] : test_pairs) {
                    auto cigar = align::globalAlignKSW2(ref, query);
                    (void)cigar;
                }
                ksw2_time = timer.elapsedMs() / NUM_RUNS;
            }

            // KSW2 extend mode (from anchors)
            {
                Timer timer;
                for (const auto& [ref, query] : test_pairs) {
                    auto cigar = align::extendAlignKSW2(ref, query, 200);
                    (void)cigar;
                }
                extend_time = timer.elapsedMs() / NUM_RUNS;
            }

            // WFA2 global alignment
            {
                Timer timer;
                for (const auto& [ref, query] : test_pairs) {
                    auto cigar = align::globalAlignWFA2(ref, query);
                    (void)cigar;
                }
                wfa2_time = timer.elapsedMs() / NUM_RUNS;
            }

            // MM2 anchor-assisted alignment
            {
                Timer timer;
                for (size_t i = 0; i < test_pairs.size(); ++i) {
                    const auto& [ref, query] = test_pairs[i];
                    auto cigar = align::globalAlignMM2(ref, query, anchors_list[i]);
                    (void)cigar;
                }
                mm2_time = timer.elapsedMs() / NUM_RUNS;
            }

            std::cout << std::setw(12) << query_len
                      << std::setw(15) << std::fixed << std::setprecision(3) << ksw2_time
                      << std::setw(15) << extend_time
                      << std::setw(15) << wfa2_time
                      << std::setw(15) << mm2_time << "\n";
        }

        std::cout << "\nNotes:\n";
        std::cout << "  - This test simulates a real scenario: aligning a long reference sequence with variable-length query sequences\n";
        std::cout << "  - Query is extracted from a region of Ref and mutated (95% similarity)\n";
        std::cout << "  - Anchor-assisted algorithms (MM2/Extend) should significantly reduce the search space in this scenario\n";
        std::cout << "  - KSW2 global alignment has O(m*n) complexity and incurs noticeable overhead when ref is very long\n";
        std::cout << "========================================================\n\n";
    }

    // ------------------------------------------------------------------
    // Performance test: low similarity sequences (70%-90%, challenging scenarios)
    // ------------------------------------------------------------------
    TEST_CASE("Performance - Low similarity sequences (70%-90%)") {
        constexpr int NUM_RUNS = 100;
        constexpr size_t SEQ_LEN = 1000;

        std::cout << "\n========== Low-similarity sequence performance test (1000bp, " << NUM_RUNS << " runs) ==========\n";
        std::cout << "Note: tests algorithm performance on high-mutation-rate sequences\n\n";

        struct LowSimilarityLevel {
            double snp_rate;
            double indel_rate;
            const char* desc;
        };

        std::vector<LowSimilarityLevel> levels = {
            {0.05,  0.02,  "90% similarity"},
            {0.10,  0.05,  "85% similarity"},
            {0.15,  0.08,  "75% similarity"},
            {0.20,  0.10,  "70% similarity"}
        };

        std::cout << std::setw(20) << "Similarity level"
                  << std::setw(15) << "KSW2(ms)"
                  << std::setw(15) << "Extend(ms)"
                  << std::setw(15) << "WFA2(ms)"
                  << std::setw(15) << "MM2(ms)" << "\n";
        std::cout << std::string(80, '-') << "\n";

        for (const auto& level : levels) {
            // Generate test data
            std::vector<std::pair<std::string, std::string>> test_pairs;
            std::vector<anchor::Anchors> anchors_list;

            for (int i = 0; i < NUM_RUNS; ++i) {
                std::string ref = generateRandomDNA(SEQ_LEN, i * 2);
                std::string query = mutateSequence(ref, level.snp_rate, level.indel_rate, i * 2 + 1);
                test_pairs.emplace_back(ref, query);

                // Low similarity scenario: anchors are sparser (one every 300bp)
                anchor::Anchors anchors;
                for (size_t pos = 0; pos + 15 < SEQ_LEN; pos += 300) {
                    anchor::Anchor a;
                    a.hash = pos * 1000 + i;
                    a.rid_ref = 0;
                    a.pos_ref = static_cast<uint32_t>(pos);
                    a.rid_qry = 0;
                    a.pos_qry = static_cast<uint32_t>(pos);
                    a.span = 15;
                    a.is_rev = false;
                    anchors.push_back(a);
                }
                anchors_list.push_back(anchors);
            }

            double ksw2_time = 0.0, extend_time = 0.0, wfa2_time = 0.0, mm2_time = 0.0;

            // KSW2
            {
                Timer timer;
                for (const auto& [ref, query] : test_pairs) {
                    auto cigar = align::globalAlignKSW2(ref, query);
                    (void)cigar;
                }
                ksw2_time = timer.elapsedMs() / NUM_RUNS;
            }

            // Extend
            {
                Timer timer;
                for (const auto& [ref, query] : test_pairs) {
                    auto cigar = align::extendAlignKSW2(ref, query, 200);
                    (void)cigar;
                }
                extend_time = timer.elapsedMs() / NUM_RUNS;
            }

            // WFA2 (using estimated similarity)
            {
                Timer timer;
                for (const auto& [ref, query] : test_pairs) {
                    double similarity = 1.0 - level.snp_rate - level.indel_rate;
                    auto cigar = align::globalAlignWFA2(ref, query);
                    (void)cigar;
                }
                wfa2_time = timer.elapsedMs() / NUM_RUNS;
            }

            // MM2
            {
                Timer timer;
                for (size_t i = 0; i < test_pairs.size(); ++i) {
                    const auto& [ref, query] = test_pairs[i];
                    auto cigar = align::globalAlignMM2(ref, query, anchors_list[i]);
                    (void)cigar;
                }
                mm2_time = timer.elapsedMs() / NUM_RUNS;
            }

            std::cout << std::setw(20) << level.desc
                      << std::setw(15) << std::fixed << std::setprecision(3) << ksw2_time
                      << std::setw(15) << extend_time
                      << std::setw(15) << wfa2_time
                      << std::setw(15) << mm2_time << "\n";
        }

        std::cout << "\nNote: at low similarity, MM2 anchors may be less reliable and performance gains can weaken\n";
        std::cout << "========================================================\n\n";
    }


    // ------------------------------------------------------------------
    // Performance test: globalAlignMM2 vs globalAlignKSW2 (with anchors vs without anchors)
    // ------------------------------------------------------------------
    TEST_CASE("Performance - globalAlignMM2 with anchors") {
        constexpr int NUM_RUNS = 100;
        constexpr size_t SEQ_LEN = 2000;

        std::cout << "\n========== globalAlignMM2 performance test (2000bp, " << NUM_RUNS << " runs) ==========\n";
        std::cout << "Note: compares anchor-assisted MM2 alignment vs pure global alignment\n\n";

        // Generate test data: high similarity sequences (98%)
        std::vector<std::pair<std::string, std::string>> test_pairs;
        std::vector<anchor::Anchors> anchors_list;

        for (int i = 0; i < NUM_RUNS; ++i) {
            std::string ref = generateRandomDNA(SEQ_LEN, i * 3);
            std::string query = mutateSequence(ref, 0.01, 0.01, i * 3 + 1);
            test_pairs.emplace_back(ref, query);

            // Generate simulated anchors for each sequence pair
            // In high similarity sequences, create one anchor every 200bp
            anchor::Anchors anchors;
            for (size_t pos = 0; pos + 50 < SEQ_LEN; pos += 200) {
                anchor::Anchor a;
                a.hash = pos * 100 + i;
                a.rid_ref = 0;
                a.pos_ref = static_cast<uint32_t>(pos);
                a.rid_qry = 0;
                a.pos_qry = static_cast<uint32_t>(pos);  // Assuming good alignment
                a.span = 50;
                a.is_rev = false;
                anchors.push_back(a);
            }
            anchors_list.push_back(anchors);
        }

        // Test globalAlignKSW2 (no anchors, pure global alignment)
        double ksw2_time = 0.0;
        {
            Timer timer;
            for (const auto& [ref, query] : test_pairs) {
                auto cigar = align::globalAlignKSW2(ref, query);
                (void)cigar;
            }
            ksw2_time = timer.elapsedMs();
        }

        // Test globalAlignMM2 (with anchors)
        double mm2_time = 0.0;
        {
            Timer timer;
            for (size_t i = 0; i < test_pairs.size(); ++i) {
                const auto& [ref, query] = test_pairs[i];
                auto cigar = align::globalAlignMM2(ref, query, anchors_list[i]);
                (void)cigar;
            }
            mm2_time = timer.elapsedMs();
        }

        // Test globalAlignMM2 (empty anchors, should degrade to KSW2)
        double mm2_empty_time = 0.0;
        {
            Timer timer;
            anchor::Anchors empty;
            for (const auto& [ref, query] : test_pairs) {
                auto cigar = align::globalAlignMM2(ref, query, empty);
                (void)cigar;
            }
            mm2_empty_time = timer.elapsedMs();
        }

        std::cout << "  globalAlignKSW2:              " << std::fixed << std::setprecision(2)
                  << ksw2_time << " ms (" << (ksw2_time / NUM_RUNS) << " ms/run)\n";
        std::cout << "  globalAlignMM2 (with anchors):      " << std::fixed << std::setprecision(2)
                  << mm2_time << " ms (" << (mm2_time / NUM_RUNS) << " ms/run)\n";
        std::cout << "  globalAlignMM2 (empty anchors):      " << std::fixed << std::setprecision(2)
                  << mm2_empty_time << " ms (" << (mm2_empty_time / NUM_RUNS) << " ms/run)\n";

        double speedup = ksw2_time / mm2_time;
        std::cout << "\n  Speedup (with anchors vs pure KSW2):    " << std::fixed << std::setprecision(2)
                  << speedup << "x\n";

        std::cout << "========================================================\n\n";
    }

    // ------------------------------------------------------------------
    // Performance test: globalAlignMM2 performance under different anchor densities
    // ------------------------------------------------------------------
    TEST_CASE("Performance - globalAlignMM2 anchor density") {
        constexpr int NUM_RUNS = 50;
        constexpr size_t SEQ_LEN = 3000;

        std::cout << "\n========== globalAlignMM2 anchor density test (3000bp, " << NUM_RUNS << " runs) ==========\n";

        // Generate test data
        std::vector<std::pair<std::string, std::string>> test_pairs;
        for (int i = 0; i < NUM_RUNS; ++i) {
            std::string ref = generateRandomDNA(SEQ_LEN, i * 4);
            std::string query = mutateSequence(ref, 0.01, 0.005, i * 4 + 1);
            test_pairs.emplace_back(ref, query);
        }

        // Test different anchor intervals (densities)
        std::vector<size_t> anchor_intervals = {500, 300, 200, 100, 50};

        std::cout << std::setw(15) << "Anchor interval (bp)"
                  << std::setw(15) << "Anchor count"
                  << std::setw(15) << "Elapsed (ms)" << "\n";
        std::cout << std::string(45, '-') << "\n";

        for (size_t interval : anchor_intervals) {
            // Generate anchors with the corresponding density
            std::vector<anchor::Anchors> anchors_list;
            size_t avg_anchor_count = 0;

            for (size_t run = 0; run < test_pairs.size(); ++run) {
                anchor::Anchors anchors;
                for (size_t pos = 0; pos + 50 < SEQ_LEN; pos += interval) {
                    anchor::Anchor a;
                    a.hash = pos * 1000 + run;
                    a.rid_ref = 0;
                    a.pos_ref = static_cast<uint32_t>(pos);
                    a.rid_qry = 0;
                    a.pos_qry = static_cast<uint32_t>(pos);
                    a.span = 50;
                    a.is_rev = false;
                    anchors.push_back(a);
                }
                avg_anchor_count += anchors.size();
                anchors_list.push_back(anchors);
            }
            avg_anchor_count /= test_pairs.size();

            // Measure performance
            Timer timer;
            for (size_t i = 0; i < test_pairs.size(); ++i) {
                const auto& [ref, query] = test_pairs[i];
                auto cigar = align::globalAlignMM2(ref, query, anchors_list[i]);
                (void)cigar;
            }
            double elapsed = timer.elapsedMs();

            std::cout << std::setw(15) << interval
                      << std::setw(15) << avg_anchor_count
                      << std::setw(15) << std::fixed << std::setprecision(2)
                      << (elapsed / NUM_RUNS) << "\n";
        }

        std::cout << "========================================================\n\n";
    }

    // ------------------------------------------------------------------
    // Test: globalAlignMM2 - anchor-based global alignment
    // ------------------------------------------------------------------
    TEST_CASE("globalAlignMM2 - Empty anchors degrade to global alignment") {
        std::string ref = "ACGTACGTACGT";
        std::string query = "ACGTACCGTACGT";  // The 6th position is mismatched.

        anchor::Anchors empty_anchors;
        auto cigar = align::globalAlignMM2(ref, query, empty_anchors);

        // It should degenerate into globalAlignKSW2
        REQUIRE(cigar.size() > 0);
        std::string cigar_str = cigarToString(cigar);
        MESSAGE("CIGAR (empty anchors): ", cigar_str);
        CHECK(!cigar_str.empty());
    }

    TEST_CASE("globalAlignMM2 - Single anchor") {
        std::string ref = "ACGTACGTACGT";
        std::string query = "ACGTACGTACGT";

        // Create an anchor point: at position 4, length 4.
        anchor::Anchors anchors;
        anchor::Anchor a;
        a.hash = 12345;
        a.rid_ref = 0;
        a.pos_ref = 4;
        a.rid_qry = 0;
        a.pos_qry = 4;
        a.span = 4;
        a.is_rev = false;
        anchors.push_back(a);

        auto cigar = align::globalAlignMM2(ref, query, anchors);
        REQUIRE(cigar.size() > 0);

        std::string cigar_str = cigarToString(cigar);
        MESSAGE("CIGAR (single anchor): ", cigar_str);

        // Verify the sequence length consumed by CIGAR
        std::size_t ref_len = cigar::getRefLength(cigar);
        std::size_t qry_len = cigar::getQueryLength(cigar);
        CHECK(ref_len == ref.size());
        CHECK(qry_len == query.size());
    }

    TEST_CASE("globalAlignMM2 - Multiple anchors form a chain") {
        std::string ref = "ACGTACGTACGTACGTACGT";  // 20bp
        std::string query = "ACGTACGTACGTACGTACGT"; // exact match

        // Create multiple anchor points
        anchor::Anchors anchors;

        // Anchor point 1: pos=0, span=4
        anchor::Anchor a1;
        a1.hash = 1001;
        a1.rid_ref = 0;
        a1.pos_ref = 0;
        a1.rid_qry = 0;
        a1.pos_qry = 0;
        a1.span = 4;
        a1.is_rev = false;
        anchors.push_back(a1);

        // Anchor point 2: pos=8, span=4
        anchor::Anchor a2;
        a2.hash = 1002;
        a2.rid_ref = 0;
        a2.pos_ref = 8;
        a2.rid_qry = 0;
        a2.pos_qry = 8;
        a2.span = 4;
        a2.is_rev = false;
        anchors.push_back(a2);

        // Anchor point 3: pos=16, span=4
        anchor::Anchor a3;
        a3.hash = 1003;
        a3.rid_ref = 0;
        a3.pos_ref = 16;
        a3.rid_qry = 0;
        a3.pos_qry = 16;
        a3.span = 4;
        a3.is_rev = false;
        anchors.push_back(a3);

        auto cigar = align::globalAlignMM2(ref, query, anchors);
        REQUIRE(cigar.size() > 0);

        std::string cigar_str = cigarToString(cigar);
        MESSAGE("CIGAR (multiple anchors): ", cigar_str);

        // Verify complete coverage
        std::size_t ref_len = cigar::getRefLength(cigar);
        std::size_t qry_len = cigar::getQueryLength(cigar);
        CHECK(ref_len == ref.size());
        CHECK(qry_len == query.size());
    }

    TEST_CASE("globalAlignMM2 - Gaps between anchors (small gaps)") {
        std::string ref =   "AAAA----CCCC----GGGG";  // 20bp (12bp after removing '-')
        std::string query = "AAAATTTTCCCCTTTTGGGG";  // 20bp

        // Actual sequence (without gaps)
        std::string ref_actual = "AAAACCCCGGGG";  // 12bp
        std::string query_actual = "AAAATTTTCCCCTTTTGGGG";  // 20bp

        // Create anchor points: only in the matched area
        anchor::Anchors anchors;

        // Anchor 1: AAAA (ref: 0-3, query: 0-3)
        anchor::Anchor a1;
        a1.hash = 2001;
        a1.rid_ref = 0;
        a1.pos_ref = 0;
        a1.rid_qry = 0;
        a1.pos_qry = 0;
        a1.span = 4;
        a1.is_rev = false;
        anchors.push_back(a1);

        // Anchor 2: CCCC (ref: 4-7, query: 8-11)
        anchor::Anchor a2;
        a2.hash = 2002;
        a2.rid_ref = 0;
        a2.pos_ref = 4;
        a2.rid_qry = 0;
        a2.pos_qry = 8;
        a2.span = 4;
        a2.is_rev = false;
        anchors.push_back(a2);

        // Anchor 3: GGGG (ref: 8-11, query: 16-19)
        anchor::Anchor a3;
        a3.hash = 2003;
        a3.rid_ref = 0;
        a3.pos_ref = 8;
        a3.rid_qry = 0;
        a3.pos_qry = 16;
        a3.span = 4;
        a3.is_rev = false;
        anchors.push_back(a3);

        auto cigar = align::globalAlignMM2(ref_actual, query_actual, anchors);
        REQUIRE(cigar.size() > 0);

        std::string cigar_str = cigarToString(cigar);
        MESSAGE("CIGAR (with small gaps): ", cigar_str);

        // Verify complete coverage
        std::size_t ref_len = cigar::getRefLength(cigar);
        std::size_t qry_len = cigar::getQueryLength(cigar);
        CHECK(ref_len == ref_actual.size());
        CHECK(qry_len == query_actual.size());

        // Should contain insertion operations (query is longer than ref)
        bool has_insertion = cigar_str.find('I') != std::string::npos;
        CHECK(has_insertion);
    }

    TEST_CASE("globalAlignMM2 - Large gaps between anchors (adaptive strategy test)") {
        // Create a 1000bp reference sequence
        std::string ref = generateRandomDNA(1000, 5000);
        std::string query = ref;  // Complete match first

        // Insert 150bp in the middle of the query (test large gap handling)
        query.insert(500, generateRandomDNA(150, 5001));

        // Create anchors covering regions before and after the insertion
        anchor::Anchors anchors;

        // Anchor 1: first half (ref: 0-99, query: 0-99)
        anchor::Anchor a1;
        a1.hash = 3001;
        a1.rid_ref = 0;
        a1.pos_ref = 0;
        a1.rid_qry = 0;
        a1.pos_qry = 0;
        a1.span = 100;
        a1.is_rev = false;
        anchors.push_back(a1);

        // Anchor 2: second half (ref: 500-599, query: 650-749)
        // query positions are shifted by 150bp (insertion length)
        anchor::Anchor a2;
        a2.hash = 3002;
        a2.rid_ref = 0;
        a2.pos_ref = 500;
        a2.rid_qry = 0;
        a2.pos_qry = 650;
        a2.span = 100;
        a2.is_rev = false;
        anchors.push_back(a2);

        auto cigar = align::globalAlignMM2(ref, query, anchors);
        REQUIRE(cigar.size() > 0);

        std::string cigar_str = cigarToString(cigar);
        MESSAGE("CIGAR (large gap): ", cigar_str);

        // Verify complete coverage
        std::size_t ref_len = cigar::getRefLength(cigar);
        std::size_t qry_len = cigar::getQueryLength(cigar);
        CHECK(ref_len == ref.size());
        CHECK(qry_len == query.size());
    }

    TEST_CASE("globalAlignMM2 - Result consistency with globalAlignKSW2 (no anchors)") {
        std::string ref = generateRandomDNA(500, 6000);
        std::string query = mutateSequence(ref, 0.02, 0.01, 6001);

        anchor::Anchors empty_anchors;
        auto cigar_mm2 = align::globalAlignMM2(ref, query, empty_anchors);
        auto cigar_ksw2 = align::globalAlignKSW2(ref, query);

        // Both should produce the same result (or at least be the same length).
        std::size_t mm2_ref_len = cigar::getRefLength(cigar_mm2);
        std::size_t mm2_qry_len = cigar::getQueryLength(cigar_mm2);
        std::size_t ksw2_ref_len = cigar::getRefLength(cigar_ksw2);
        std::size_t ksw2_qry_len = cigar::getQueryLength(cigar_ksw2);

        CHECK(mm2_ref_len == ksw2_ref_len);
        CHECK(mm2_qry_len == ksw2_qry_len);
        CHECK(mm2_ref_len == ref.size());
        CHECK(mm2_qry_len == query.size());
    }

    TEST_CASE("removeRefGapColumns - drop ref gap columns from aligned seq") {
        // Note: test the pure filtering behavior of "remove columns by ref_gap_pos" (in-place modification).
        // The input sequence is already aligned (contains gaps); this function only removes columns where the reference has gaps.

        // Input: aligned sequence "AC-GT" (length 5)
        std::string seq = "AC-GT";

        // ref_gap_pos: true indicates the reference has a gap at that column and it should be removed.
        // We remove column 2 (0-based), resulting in "ACGT"
        const std::vector<bool> ref_gap_pos = {false, false, true, false, false};

        align::RefAligner::removeRefGapColumns(seq, ref_gap_pos);
        CHECK(seq == "ACGT");
    }

    TEST_CASE("removeRefGapColumns - keep existing '-' as base when not in ref gap pos") {
        // Note: The input sequence itself contains '-', but it is preserved as long as the ref_gap_pos column is false.

        std::string seq = "A-CG"; // contains one original '-'
        const std::vector<bool> ref_gap_pos = {false, false, false, false}; // all columns are retained

        align::RefAligner::removeRefGapColumns(seq, ref_gap_pos);
        CHECK(seq == "A-CG");
    }


}

// ==================================================================
// Alignment accuracy test suite
// ==================================================================
TEST_SUITE("align") {

    // ------------------------------------------------------------------
    // Auxiliary function: Verify the correctness of CIGAR
    // ------------------------------------------------------------------
    static bool verifyCigar(const std::string& ref, const std::string& query,
                           const cigar::Cigar_t& cigar) {
        std::size_t ref_pos = 0;
        std::size_t qry_pos = 0;

        for (const auto& op : cigar) {
            char op_char;
            uint32_t len;
            cigar::intToCigar(op, op_char, len);

            switch (op_char) {
                case 'M':
                case '=':
                case 'X':
                    ref_pos += len;
                    qry_pos += len;
                    break;
                case 'I':
                    qry_pos += len;
                    break;
                case 'D':
                case 'N':
                    ref_pos += len;
                    break;
                case 'S':
                case 'H':
                    qry_pos += len;
                    break;
                default:
                    return false;
            }
        }

        return ref_pos == ref.size() && qry_pos == query.size();
    }

    // ------------------------------------------------------------------
    // Auxiliary function: Calculate the edit distance of CIGAR (simplified version)
    // ------------------------------------------------------------------
    static size_t getCigarEditDistance(const cigar::Cigar_t& cigar) {
        size_t edit_dist = 0;

        for (const auto& op : cigar) {
            char op_char;
            uint32_t len;
            cigar::intToCigar(op, op_char, len);

            if (op_char == 'I' || op_char == 'D') {
                edit_dist += len;
            }
        }

        return edit_dist;
    }

    // ------------------------------------------------------------------
    // Test: Accuracy of high similarity (95%-99%) comparison
    // ------------------------------------------------------------------
    TEST_CASE("Accuracy - High similarity (95%-99%)") {
        constexpr int NUM_TESTS = 50;
        constexpr size_t SEQ_LEN = 5000;

        std::cout << "\n========== High-similarity alignment accuracy test (500bp, " << NUM_TESTS << " runs) ==========\n";

        struct AccuracyStats {
            size_t total_tests = 0;
            size_t valid_cigars = 0;
            size_t perfect_match = 0;
            double avg_edit_dist = 0.0;
        };

        std::vector<double> similarity_levels = {0.99, 0.98, 0.97, 0.95};

        for (double similarity : similarity_levels) {
            double snp_rate = (1.0 - similarity) * 0.7;
            double indel_rate = (1.0 - similarity) * 0.3;

            std::map<std::string, AccuracyStats> algorithm_stats;
            algorithm_stats["KSW2"];
            algorithm_stats["WFA2"];
            algorithm_stats["MM2"];

            // Anchor statistics
            size_t total_anchors = 0;
            size_t min_anchors = std::numeric_limits<size_t>::max();
            size_t max_anchors = 0;

            for (int i = 0; i < NUM_TESTS; ++i) {
                std::string ref = generateRandomDNA(SEQ_LEN, i * 3);
                std::string query = mutateSequence(ref, snp_rate, indel_rate, i * 3 + 1);

                // Generate real anchors (based on minimizer matching)
                // k=15, w=10: Parameters suitable for high similarity scenarios
                anchor::Anchors anchors = generateRealAnchors(ref, query, 15, 10);

                // Count anchors
                total_anchors += anchors.size();
                min_anchors = std::min(min_anchors, anchors.size());
                max_anchors = std::max(max_anchors, anchors.size());

                // KSW2
                {
                    auto cigar = align::globalAlignKSW2(ref, query);
                    auto& stats = algorithm_stats["KSW2"];
                    stats.total_tests++;
                    if (verifyCigar(ref, query, cigar)) {
                        stats.valid_cigars++;
                        size_t edit_dist = getCigarEditDistance(cigar);
                        stats.avg_edit_dist += edit_dist;
                        if (edit_dist == 0) stats.perfect_match++;
                    }
                }

                // WFA2
                {
                    auto cigar = align::globalAlignWFA2(ref, query);
                    auto& stats = algorithm_stats["WFA2"];
                    stats.total_tests++;
                    if (verifyCigar(ref, query, cigar)) {
                        stats.valid_cigars++;
                        size_t edit_dist = getCigarEditDistance(cigar);
                        stats.avg_edit_dist += edit_dist;
                        if (edit_dist == 0) stats.perfect_match++;
                    }
                }

                // MM2
                {
                    auto cigar = align::globalAlignMM2(ref, query, anchors);
                    auto& stats = algorithm_stats["MM2"];
                    stats.total_tests++;
                    if (verifyCigar(ref, query, cigar)) {
                        stats.valid_cigars++;
                        size_t edit_dist = getCigarEditDistance(cigar);
                        stats.avg_edit_dist += edit_dist;
                        if (edit_dist == 0) stats.perfect_match++;
                    }
                }
            }

            std::cout << "\nSimilarity " << (similarity * 100) << "%:\n";
            std::cout << "  [Anchors] Avg=" << (total_anchors / NUM_TESTS)
                      << ", Min=" << (min_anchors == std::numeric_limits<size_t>::max() ? 0 : min_anchors)
                      << ", Max=" << max_anchors << "\n";
            for (const auto& [name, stats] : algorithm_stats) {
                double accuracy = 100.0 * stats.valid_cigars / stats.total_tests;
                double avg_dist = stats.avg_edit_dist / stats.total_tests;
                std::cout << "  " << std::setw(6) << name
                          << ": Accuracy=" << std::fixed << std::setprecision(1) << accuracy << "%"
                          << ", Avg edit distance=" << std::setprecision(2) << avg_dist
                          << ", Perfect matches=" << stats.perfect_match << "/" << stats.total_tests << "\n";
            }
        }

        std::cout << "========================================================\n\n";
    }

    // ------------------------------------------------------------------
    // Test: Accuracy of low similarity (70%-90%) comparison
    // ------------------------------------------------------------------
    TEST_CASE("Accuracy - Low similarity (70%-90%)") {
        constexpr int NUM_TESTS = 50;
        constexpr size_t SEQ_LEN = 5000;

        std::cout << "\n========== Low-similarity alignment accuracy test (500bp, " << NUM_TESTS << " runs) ==========\n";

        struct AccuracyStats {
            size_t total_tests = 0;
            size_t valid_cigars = 0;
            size_t failed = 0;
            double avg_edit_dist = 0.0;
        };

        std::vector<double> similarity_levels = {0.90, 0.85, 0.80, 0.75, 0.70};

        for (double similarity : similarity_levels) {
            double snp_rate = (1.0 - similarity) * 0.6;
            double indel_rate = (1.0 - similarity) * 0.4;

            std::map<std::string, AccuracyStats> algorithm_stats;
            algorithm_stats["KSW2"];
            algorithm_stats["WFA2"];
            algorithm_stats["MM2"];

            for (int i = 0; i < NUM_TESTS; ++i) {
                std::string ref = generateRandomDNA(SEQ_LEN, i * 3);
                std::string query = mutateSequence(ref, snp_rate, indel_rate, i * 3 + 1);

                // Generate real anchors (low similarity scenario: use smaller k and w)
                // k=13, w=8: at low similarity, use looser parameters to obtain enough anchors
                anchor::Anchors anchors = generateRealAnchors(ref, query, 13, 8);

                // KSW2
                {
                    auto cigar = align::globalAlignKSW2(ref, query);
                    auto& stats = algorithm_stats["KSW2"];
                    stats.total_tests++;
                    if (verifyCigar(ref, query, cigar)) {
                        stats.valid_cigars++;
                        stats.avg_edit_dist += getCigarEditDistance(cigar);
                    } else {
                        stats.failed++;
                    }
                }

                // WFA2
                {
                    auto cigar = align::globalAlignWFA2(ref, query);
                    auto& stats = algorithm_stats["WFA2"];
                    stats.total_tests++;
                    if (verifyCigar(ref, query, cigar)) {
                        stats.valid_cigars++;
                        stats.avg_edit_dist += getCigarEditDistance(cigar);
                    } else {
                        stats.failed++;
                    }
                }

                // MM2
                {
                    auto cigar = align::globalAlignMM2(ref, query, anchors);
                    auto& stats = algorithm_stats["MM2"];
                    stats.total_tests++;
                    if (verifyCigar(ref, query, cigar)) {
                        stats.valid_cigars++;
                        stats.avg_edit_dist += getCigarEditDistance(cigar);
                    } else {
                        stats.failed++;
                    }
                }
            }

            std::cout << "\nSimilarity " << (similarity * 100) << "%:\n";
            for (const auto& [name, stats] : algorithm_stats) {
                double accuracy = 100.0 * stats.valid_cigars / stats.total_tests;
                double avg_dist = stats.avg_edit_dist / (stats.valid_cigars > 0 ? stats.valid_cigars : 1);
                std::cout << "  " << std::setw(6) << name
                          << ": Accuracy=" << std::fixed << std::setprecision(1) << accuracy << "%"
                          << ", Avg edit distance=" << std::setprecision(2) << avg_dist
                          << ", Failures=" << stats.failed << "/" << stats.total_tests << "\n";
            }
        }

        std::cout << "\nNote: at low similarity, some algorithms may fail because of parameter limits\n";
        std::cout << "========================================================\n\n";
    }



}

// ------------------------------------------------------------------
// Quick guide: running performance tests
// ------------------------------------------------------------------
// Notes:
// 1. Correctness tests run by default to validate basic functionality
// 2. Performance tests are skipped by default (doctest::skip(true)); they must be manually enabled
// 3. Run performance tests:
//    ./halign4_tests -tc="*Performance*" --no-skip
// 4. Run only specific performance tests:
//    ./halign4_tests -tc="*Short*" --no-skip
// ------------------------------------------------------------------
