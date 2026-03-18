// ==================================================================
// test_anchor.cpp - unit tests for the anchor module
// ==================================================================
//
// Test coverage:
// 1. collect_anchors: anchor collection (based on minimizer hits)
// 2. chainAnchors: chaining algorithm (DP dynamic programming)
// 3. filter params: q_occ_frac, f_top_frac, sample_every_bp
// ==================================================================

#include <doctest/doctest.h>
#include "seed.h"
#include "anchor.h"

#include <algorithm>
#include <cstdint>
#include <string>
#include <vector>

using hash_t = std::uint64_t;

// ==================================================================
// Helper functions
// ==================================================================

// Create a MinimizerHit (for collect_anchors tests)
static minimizer::MinimizerHit makeHit(hash_t hash56, std::uint32_t pos,
                                       std::uint32_t rid = 0, bool strand = false, std::uint32_t span = 15)
{
    // Use convenient constructor: MinimizerHit(hash56, pos, rid, strand, span)
    return minimizer::MinimizerHit(hash56, pos, rid, strand, static_cast<std::uint8_t>(span));
}

// Create an Anchor (for chainAnchors tests)
static anchor::Anchor makeAnchor(hash_t hash, std::uint32_t pos_ref, std::uint32_t pos_qry,
                                  std::uint32_t rid_ref = 0, std::uint32_t rid_qry = 0,
                                  bool is_rev = false, std::uint32_t span = 15)
{
    anchor::Anchor a;
    a.hash = hash;
    a.rid_ref = rid_ref;
    a.pos_ref = pos_ref;
    a.rid_qry = rid_qry;
    a.pos_qry = pos_qry;
    a.span = span;
    a.is_rev = is_rev;
    return a;
}

// ==================================================================
// TEST SUITE: collect_anchors - anchor collection tests
// ==================================================================

TEST_SUITE("anchor")
{
    // ------------------------------------------------------------------
    // Basic functionality tests
    // ------------------------------------------------------------------

    TEST_CASE("collect_anchors - empty input returns empty")
    {
        minimizer::MinimizerHits ref_hits;
        minimizer::MinimizerHits qry_hits;

        auto anchors = minimizer::collect_anchors(ref_hits, qry_hits);
        CHECK(anchors.empty());
    }

    TEST_CASE("collect_anchors - no anchors when ref is empty")
    {
        minimizer::MinimizerHits ref_hits;
        minimizer::MinimizerHits qry_hits;
        qry_hits.push_back(makeHit(0x111111, 100));

        auto anchors = minimizer::collect_anchors(ref_hits, qry_hits);
        CHECK(anchors.empty());
    }

    TEST_CASE("collect_anchors - no anchors when qry is empty")
    {
        minimizer::MinimizerHits ref_hits;
        minimizer::MinimizerHits qry_hits;
        ref_hits.push_back(makeHit(0x111111, 100));

        auto anchors = minimizer::collect_anchors(ref_hits, qry_hits);
        CHECK(anchors.empty());
    }

    TEST_CASE("collect_anchors - single perfect match")
    {
        minimizer::MinimizerHits ref_hits;
        minimizer::MinimizerHits qry_hits;

        ref_hits.push_back(makeHit(0x123456, 100, 0, false, 20));
        qry_hits.push_back(makeHit(0x123456, 50, 0, false, 20));

        // Disable all filtering
        anchor::SeedFilterParams params;
        params.q_occ_frac = 0.0;
        params.f_top_frac = 0.0;

        auto anchors = minimizer::collect_anchors(ref_hits, qry_hits, params);

        REQUIRE(anchors.size() == 1);
        CHECK(anchors[0].hash == 0x123456);
        CHECK(anchors[0].pos_ref == 100);
        CHECK(anchors[0].pos_qry == 50);
        CHECK(anchors[0].span == 20);
        CHECK(anchors[0].is_rev == false);
    }

    TEST_CASE("collect_anchors - reverse match detection")
    {
        minimizer::MinimizerHits ref_hits;
        minimizer::MinimizerHits qry_hits;

        ref_hits.push_back(makeHit(0xABCDEF, 200, 0, false, 15));
        qry_hits.push_back(makeHit(0xABCDEF, 80, 0, true, 15));  // reverse

        anchor::SeedFilterParams params;
        params.q_occ_frac = 0.0;
        params.f_top_frac = 0.0;

        auto anchors = minimizer::collect_anchors(ref_hits, qry_hits, params);

        REQUIRE(anchors.size() == 1);
        CHECK(anchors[0].is_rev == true);  // ref XOR qry = false XOR true = true
    }

    TEST_CASE("collect_anchors - no anchors when no shared hash")
    {
        minimizer::MinimizerHits ref_hits;
        minimizer::MinimizerHits qry_hits;

        ref_hits.push_back(makeHit(0x111111, 100));
        ref_hits.push_back(makeHit(0x222222, 200));
        qry_hits.push_back(makeHit(0x333333, 50));
        qry_hits.push_back(makeHit(0x444444, 80));

        auto anchors = minimizer::collect_anchors(ref_hits, qry_hits);
        CHECK(anchors.empty());
    }

    TEST_CASE("collect_anchors - one-to-many expansion (occurrence expansion)")
    {
        minimizer::MinimizerHits ref_hits;
        minimizer::MinimizerHits qry_hits;

        // same hash appears 3 times on ref side
        ref_hits.push_back(makeHit(0x555555, 100));
        ref_hits.push_back(makeHit(0x555555, 200));
        ref_hits.push_back(makeHit(0x555555, 300));

        // hash appears once on qry side
        qry_hits.push_back(makeHit(0x555555, 50));

        anchor::SeedFilterParams params;
        params.q_occ_frac = 0.0;
        params.f_top_frac = 0.0;

        auto anchors = minimizer::collect_anchors(ref_hits, qry_hits, params);

        // should generate 3 anchors (1 qry × 3 ref)
        REQUIRE(anchors.size() == 3);

        // verify all anchors share the same qry position, but different ref positions
        CHECK(anchors[0].pos_qry == 50);
        CHECK(anchors[1].pos_qry == 50);
        CHECK(anchors[2].pos_qry == 50);

        // ref positions should be 100, 200, 300 (order may vary due to sorting)
        std::vector<std::uint32_t> ref_positions;
        for (const auto& a : anchors) {
            ref_positions.push_back(a.pos_ref);
        }
        std::sort(ref_positions.begin(), ref_positions.end());
        CHECK(ref_positions[0] == 100);
        CHECK(ref_positions[1] == 200);
        CHECK(ref_positions[2] == 300);
    }

    // ------------------------------------------------------------------
    // filter parameter tests
    // ------------------------------------------------------------------

    TEST_CASE("collect_anchors - q_occ_frac filters high-frequency hash")
    {
        minimizer::MinimizerHits ref_hits;
        minimizer::MinimizerHits qry_hits;

        ref_hits.push_back(makeHit(0x888888, 100));

        // qry has 100 hits, with same hash appearing 50 times
        for (std::uint32_t i = 0; i < 50; ++i) {
            qry_hits.push_back(makeHit(0x888888, i));  // high-frequency hash
        }
        for (std::uint32_t i = 0; i < 50; ++i) {
            qry_hits.push_back(makeHit(0x999900 + i, i));  // different hash
        }

        // set q_occ_frac = 10% (i.e., 10 hits); drop if exceeded
        anchor::SeedFilterParams params;
        params.q_occ_frac = 0.10;  // 10% of 100 = 10
        params.f_top_frac = 0.0;

        auto anchors = minimizer::collect_anchors(ref_hits, qry_hits, params);

        // high-frequency hash (0x888888) should be filtered out because 50 > 10
        CHECK(anchors.empty());
    }

    TEST_CASE("collect_anchors - span is min(ref.span, qry.span)")
    {
        minimizer::MinimizerHits ref_hits;
        minimizer::MinimizerHits qry_hits;

        ref_hits.push_back(makeHit(0xAAAAAA, 100, 0, false, 30));  // span=30
        qry_hits.push_back(makeHit(0xAAAAAA, 50, 0, false, 20));   // span=20

        anchor::SeedFilterParams params;
        params.q_occ_frac = 0.0;
        params.f_top_frac = 0.0;

        auto anchors = minimizer::collect_anchors(ref_hits, qry_hits, params);

        REQUIRE(anchors.size() == 1);
        CHECK(anchors[0].span == 20);  // min(30, 20) = 20
    }

    TEST_CASE("collect_anchors - multi-sequence rid passed through correctly")
    {
        minimizer::MinimizerHits ref_hits;
        minimizer::MinimizerHits qry_hits;

        ref_hits.push_back(makeHit(0xBBBBBB, 100, 2, false, 15));  // rid_ref=2
        qry_hits.push_back(makeHit(0xBBBBBB, 50, 3, false, 15));   // rid_qry=3

        anchor::SeedFilterParams params;
        params.q_occ_frac = 0.0;
        params.f_top_frac = 0.0;

        auto anchors = minimizer::collect_anchors(ref_hits, qry_hits, params);

        REQUIRE(anchors.size() == 1);
        CHECK(anchors[0].rid_ref == 2);
        CHECK(anchors[0].rid_qry == 3);
    }
}

// ==================================================================
// TEST SUITE: chainAnchors - chaining algorithm tests
// ==================================================================

TEST_SUITE("anchor")
{
    // ------------------------------------------------------------------
    // Basic functionality tests
    // ------------------------------------------------------------------

    TEST_CASE("chainAnchors - empty input returns empty chain")
    {
        anchor::Anchors anchors;
        auto best_chain = anchor::chainAnchors(anchors);
        CHECK(best_chain.empty());
    }

    TEST_CASE("chainAnchors - single anchor forms a single chain")
    {
        anchor::Anchors anchors;
        anchors.push_back(makeAnchor(0x111111, 100, 50, 0, 0, false, 20));

        anchor::ChainParams params;
        params.min_cnt = 1;  // allow single-anchor chains
        params.min_score = 10;

        auto best_chain = anchor::chainAnchors(anchors, params);

        REQUIRE(!best_chain.empty());
        CHECK(best_chain.size() == 1);
        CHECK(best_chain[0].span == 20);
    }

    TEST_CASE("chainAnchors - two linkable anchors form a single chain")
    {
        anchor::Anchors anchors;

        // two anchors, increasing positions, moderate gap
        anchors.push_back(makeAnchor(0x111111, 100, 50, 0, 0, false, 20));
        anchors.push_back(makeAnchor(0x222222, 150, 100, 0, 0, false, 20));

        anchor::ChainParams params;
        params.min_cnt = 2;
        params.min_score = 30;

        auto best_chain = anchor::chainAnchors(anchors, params);

        REQUIRE(!best_chain.empty());
        CHECK(best_chain.size() == 2);
        // verify anchors are sorted by position
        CHECK(best_chain[0].pos_ref < best_chain[1].pos_ref);
    }

    TEST_CASE("chainAnchors - only best chain returned across different reference IDs")
    {
        anchor::Anchors anchors;

        // Chain A: rid_ref=0, higher score
        anchors.push_back(makeAnchor(0x111111, 100, 50, 0, 0, false, 20));
        anchors.push_back(makeAnchor(0x222222, 150, 100, 0, 0, false, 20));
        anchors.push_back(makeAnchor(0x555555, 200, 150, 0, 0, false, 20));

        // Chain B: rid_ref=1, lower score
        anchors.push_back(makeAnchor(0x333333, 200, 150, 1, 0, false, 15));
        anchors.push_back(makeAnchor(0x444444, 250, 200, 1, 0, false, 15));

        anchor::ChainParams params;
        params.min_cnt = 2;
        params.min_score = 20;

        auto best_chain = anchor::chainAnchors(anchors, params);

        // should return chain A with higher score (rid_ref=0)
        REQUIRE(!best_chain.empty());
        CHECK(best_chain[0].rid_ref == 0);
        CHECK(best_chain.size() >= 2);
    }

    TEST_CASE("chainAnchors - only best chain returned for forward vs reverse anchors")
    {
        anchor::Anchors anchors;

        // forward chain (higher score)
        anchors.push_back(makeAnchor(0x111111, 100, 50, 0, 0, false, 20));
        anchors.push_back(makeAnchor(0x222222, 150, 100, 0, 0, false, 20));
        anchors.push_back(makeAnchor(0x666666, 200, 150, 0, 0, false, 20));

        // reverse chain (lower score)
        anchors.push_back(makeAnchor(0x333333, 300, 250, 0, 0, true, 15));
        anchors.push_back(makeAnchor(0x444444, 350, 300, 0, 0, true, 15));

        anchor::ChainParams params;
        params.min_cnt = 2;
        params.min_score = 20;

        auto best_chain = anchor::chainAnchors(anchors, params);

        // should return higher-scoring forward chain
        REQUIRE(!best_chain.empty());
        CHECK(best_chain[0].is_rev == false);
        CHECK(best_chain.size() >= 2);
    }

    TEST_CASE("chainAnchors - anchors too far apart do not chain")
    {
        anchor::Anchors anchors;

        // two anchors, ref distance exceeds max_dist_x (default 5000)
        anchors.push_back(makeAnchor(0x111111, 100, 50, 0, 0, false, 20));
        anchors.push_back(makeAnchor(0x222222, 6000, 5500, 0, 0, false, 20));

        anchor::ChainParams params;
        params.min_cnt = 1;  // allow single-anchor chains
        params.min_score = 10;
        params.max_dist_x = 5000;

        auto best_chain = anchor::chainAnchors(anchors, params);

        // should return only one anchor (too far to chain)
        REQUIRE(!best_chain.empty());
        CHECK(best_chain.size() == 1);
    }

    TEST_CASE("chainAnchors - do not chain when diagonal offset exceeds bandwidth")
    {
        anchor::Anchors anchors;

        // two anchors, diagonal offset = |dr - dq| = |150 - 100| = 50
        // if bw < 50, should not chain
        anchors.push_back(makeAnchor(0x111111, 100, 50, 0, 0, false, 20));
        anchors.push_back(makeAnchor(0x222222, 250, 150, 0, 0, false, 20));
        // dr = 250 - 100 = 150, dq = 150 - 50 = 100, dd = |150 - 100| = 50

        anchor::ChainParams params;
        params.min_cnt = 1;
        params.min_score = 10;
        params.bw = 30;  // bw < 50, should not chain

        auto best_chain = anchor::chainAnchors(anchors, params);

        // should return only one anchor (diagonal offset exceeds bandwidth)
        REQUIRE(!best_chain.empty());
        CHECK(best_chain.size() == 1);
    }

    TEST_CASE("chainAnchors - multiple chains sorted by descending score")
    {
        anchor::Anchors anchors;

        // Chain A: 3 anchors (higher score)
        anchors.push_back(makeAnchor(0x111111, 100, 50, 0, 0, false, 20));
        anchors.push_back(makeAnchor(0x222222, 150, 100, 0, 0, false, 20));
        anchors.push_back(makeAnchor(0x333333, 200, 150, 0, 0, false, 20));

        // Chain B: 2 anchors (lower score, and far from chain A)
        anchors.push_back(makeAnchor(0x444444, 10000, 5000, 0, 0, false, 15));
        anchors.push_back(makeAnchor(0x555555, 10100, 5100, 0, 0, false, 15));

        anchor::ChainParams params;
        params.min_cnt = 2;
        params.min_score = 20;

        auto best_chain = anchor::chainAnchors(anchors, params);

        REQUIRE(!best_chain.empty());

        // should return higher-scoring chain A (3 anchors)
        CHECK(best_chain.size() == 3);
    }

    // ------------------------------------------------------------------
    // parameter filter tests
    // ------------------------------------------------------------------

    TEST_CASE("chainAnchors - min_cnt filters short chains")
    {
        anchor::Anchors anchors;

        // only 2 anchors
        anchors.push_back(makeAnchor(0x111111, 100, 50, 0, 0, false, 15));
        anchors.push_back(makeAnchor(0x222222, 150, 100, 0, 0, false, 15));

        // require at least 3 anchors
        anchor::ChainParams params;
        params.min_cnt = 3;
        params.min_score = 10;

        auto best_chain = anchor::chainAnchors(anchors, params);

        // should be no chain (does not meet min_cnt)
        CHECK(best_chain.empty());
    }

    TEST_CASE("chainAnchors - min_score filters low score chains")
    {
        anchor::Anchors anchors;

        // 2 anchors with small span, low score
        anchors.push_back(makeAnchor(0x111111, 100, 50, 0, 0, false, 5));
        anchors.push_back(makeAnchor(0x222222, 150, 100, 0, 0, false, 5));

        anchor::ChainParams params;
        params.min_cnt = 2;
        params.min_score = 50;  // require score >= 50

        auto best_chain = anchor::chainAnchors(anchors, params);

        // score too low, should be filtered
        CHECK(best_chain.empty());
    }

    TEST_CASE("chainAnchors - returned anchors sorted by position")
    {
        anchor::Anchors anchors;

        anchors.push_back(makeAnchor(0x111111, 100, 50, 0, 0, false, 20));
        anchors.push_back(makeAnchor(0x222222, 200, 150, 0, 0, false, 25));

        anchor::ChainParams params;
        params.min_cnt = 2;
        params.min_score = 30;

        auto best_chain = anchor::chainAnchors(anchors, params);

        REQUIRE(!best_chain.empty());
        REQUIRE(best_chain.size() == 2);

        // verify anchors are in positional order
        CHECK(best_chain[0].pos_ref == 100);
        CHECK(best_chain[1].pos_ref == 200);
        CHECK(best_chain[0].pos_qry == 50);
        CHECK(best_chain[1].pos_qry == 150);
    }
}


