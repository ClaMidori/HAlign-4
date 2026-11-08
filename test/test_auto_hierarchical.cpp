#include "doctest/doctest.h"
#include "align/strategy.hpp"

TEST_SUITE("auto_hierarchical") {
    TEST_CASE("opt-in bounded auto route chooses exact direct/hierarchical safely") {
        using namespace align_strategy;
        Budget budget{};
        FastaStats s{3, 18, 8, 4};
        CHECK(chooseBoundedAuto(s, budget, 3).backend == AutoBackend::NativeDirect);
        CHECK(chooseBoundedAuto(s, budget, 4).backend == AutoBackend::NativeDirect);
        CHECK(chooseBoundedAuto(s, budget, 2).backend == AutoBackend::NativeHierarchical);
        CHECK(chooseBoundedAuto(s, budget, 1).backend == AutoBackend::Reject);
        CHECK(chooseBoundedAuto(s, budget, 0).backend == AutoBackend::Reject);
        s.contains_gaps = true;
        CHECK(chooseBoundedAuto(s, budget, 2).backend == AutoBackend::Reject);
        s.contains_gaps = false;
        s.sequence_count = 65;
        CHECK(chooseBoundedAuto(s, budget, 3).backend == AutoBackend::Reject);
        s.sequence_count = 3;
        s.longest_bases = 25001;
        CHECK(chooseBoundedAuto(s, budget, 3).backend == AutoBackend::Reject);
        s.longest_bases = 8;
        s.total_bases = 100001;
        CHECK(chooseBoundedAuto(s, budget, 2).backend == AutoBackend::Reject);
        s.total_bases = 18;
        CHECK(chooseBoundedAuto(s, budget, 65).backend == AutoBackend::NativeDirect);
        s.sequence_count = 1;
        CHECK(chooseBoundedAuto(s, budget, 1).backend == AutoBackend::NativeDirect);
        s.sequence_count = 0;
        CHECK(chooseBoundedAuto(s, budget, 3).backend == AutoBackend::Reject);
    }

    TEST_CASE("bounded auto routing is independent of thread count and never loosens budgets") {
        using namespace align_strategy;
        FastaStats s{24, 40000, 400, 250};
        Budget b{};
        const auto chosen = chooseBoundedAuto(s, b, 3);
        REQUIRE(chosen.backend == AutoBackend::NativeHierarchical);
        CHECK(chooseBoundedAuto(s, b, 3).backend == chosen.backend);
        b.max_total_bases = 39999;
        CHECK(chooseBoundedAuto(s, b, 3).backend == AutoBackend::Reject);
        b.max_total_bases = 100000;
        b.max_sequence_bases = 399;
        CHECK(chooseBoundedAuto(s, b, 3).backend == AutoBackend::Reject);
    }
}
