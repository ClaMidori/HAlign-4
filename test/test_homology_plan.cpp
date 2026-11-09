#include "doctest/doctest.h"
#include "align/homology_plan.hpp"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <random>
#include <string>
#include <vector>

namespace {
std::string randomDna(std::size_t n, std::uint32_t seed) {
    std::mt19937 rng(seed);
    std::string out;
    out.reserve(n);
    for (std::size_t i = 0; i < n; ++i) out.push_back("ACGT"[rng() % 4]);
    return out;
}
std::string mutated(std::string in, std::uint32_t seed, std::size_t every) {
    std::mt19937 rng(seed);
    for (std::size_t i = 1; i < in.size(); i += every) {
        auto old = in[i];
        auto next = "ACGT"[rng() % 4];
        in[i] = next == old ? "ACGT"[(std::string("ACGT").find(old) + 1) % 4] : next;
    }
    return in;
}
}
TEST_SUITE("homology_plan") {
    TEST_CASE("diagnostic anchors and boundaries are deterministic, forward and non-overlapping") {
        auto ref = randomDna(80000, 20261108);
        auto q1 = mutated(ref, 123, 30);
        auto q2 = mutated(ref, 456, 40);
        q1.insert(34000, randomDna(700, 77));
        const std::vector<homology_plan::Sequence> reads = {
            {"reference", ref}, {"carrier", q1}, {"second", q2}};
        const auto first = homology_plan::build(reads);
        const auto second = homology_plan::build(reads);
        REQUIRE(first.boundaries.size() >= 3);
        CHECK(homology_plan::asTsv(first) == homology_plan::asTsv(second));
        CHECK(first.supported_anchors.size() == 2);
        CHECK(first.supported_anchors[0].size() > 10);
        CHECK(first.supported_anchors[1].size() > 10);
        for (const auto& anchors : first.supported_anchors) {
            const auto bad = std::adjacent_find(anchors.begin(), anchors.end(),
                [](const auto& a, const auto& b) {
                    return a.reference >= b.reference || a.query >= b.query;
                });
            CHECK(bad == anchors.end());
        }
        CHECK(first.boundaries.front().reference == 0);
        CHECK(first.boundaries.back().reference == ref.size());
        for (std::size_t j = 1; j < first.boundaries.size(); ++j)
            CHECK(first.boundaries[j].reference > first.boundaries[j-1].reference);
        const auto tsv = homology_plan::asTsv(first);
        CHECK(tsv.find("CANDIDATE_ONLY") != std::string::npos);
        CHECK(tsv.find("NO verified homologous coordinates") != std::string::npos);
    }
    TEST_CASE("unrelated and low-complexity inputs are not forced into confident homologous blocks") {
        const auto a = randomDna(60000, 101);
        const auto b = randomDna(60000, 102);
        const auto random = homology_plan::build({{"a",a}, {"b",b}});
        CHECK(random.supported_anchors[0].size() < 5);
        CHECK(random.boundaries.size() == 2);
        const auto low = homology_plan::build({{"a",std::string(60000, 'A')}, {"b",std::string(60000,'A')}});
        CHECK(low.supported_anchors[0].empty());
        CHECK(low.boundaries.size() == 2);
        CHECK(homology_plan::asTsv(low).find("UNRESOLVED") != std::string::npos);
    }
    TEST_CASE("budget and DNA input checks reject unsafe or invalid data") {
        CHECK_THROWS_AS(homology_plan::build({{"single","ACGT"}}), std::invalid_argument);
        CHECK_THROWS_AS(homology_plan::build({{"a","ACGT"},{"a","ACGA"}}), std::invalid_argument);
        CHECK_THROWS_AS(homology_plan::build({{"a","ACGT-"},{"b","ACGT"}}), std::invalid_argument);
        CHECK_THROWS_AS(homology_plan::build({{"a","ACGY"},{"b","ACGT"}}), std::invalid_argument);
        CHECK_THROWS_AS(homology_plan::build({{"a",std::string(1000001, 'A')},{"b","ACGT"}}), std::invalid_argument);
        CHECK_THROWS_AS(homology_plan::build({{"a","ACGT"},{"b","ACGT"}}, {.k=9}), std::invalid_argument);
        const auto dir = std::filesystem::temp_directory_path();
        const auto file = dir / "halign4_pr3d1_doctest.fa";
        {
            std::ofstream out(file);
            out << ">a\nACGT\n>b\nACGN\n";
        }
        const auto reads = homology_plan::readFasta(file);
        REQUIRE(reads.size() == 2);
        CHECK(reads[1].bases == "ACGN");
        std::error_code ec;
        std::filesystem::remove(file, ec);
    }
}
