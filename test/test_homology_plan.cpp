#include "doctest/doctest.h"
#include "align/homology_plan.hpp"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <random>
#include <sstream>
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
    TEST_CASE("PR3D.2B opt-in guard vetoes discontinuous spans without changing default") {
        homology_plan::Plan p;
        p.ids = {"reference", "query"};
        p.lengths = {20000, 50000};
        p.supported_anchors = {{{100,100}, {5100,5100}, {10100,10100}, {15100,15100}}};
        homology_plan::Boundary left, right;
        left.reference = 0;
        left.projections = {homology_plan::Projection{0,0,0,0}};
        right.reference = 20000;
        right.projections = {homology_plan::Projection{50000,20000,50000,0}};
        p.boundaries = {left,right};
        const auto old = homology_plan::asTsv(p);
        CHECK(old.find("CANDIDATE_ONLY") != std::string::npos);
        p.config.strict_candidate_guard = true;
        const auto strict = homology_plan::asTsv(p);
        CHECK(strict.find("CANDIDATE_ONLY") == std::string::npos);
        CHECK(strict.find("UNRESOLVED") != std::string::npos);
        // Only a plausible span with distributed diagonal support may pass.
        p.lengths[1] = 20000;
        p.boundaries[1].projections[0] = homology_plan::Projection{20000,20000,20000,0};
        CHECK(homology_plan::asTsv(p).find("CANDIDATE_ONLY") != std::string::npos);
        p.config.strict_candidate_guard = false;
        CHECK(homology_plan::asTsv(p).find("PR3D.2B") == std::string::npos);
    }


    TEST_CASE("PR3D.2C flanks are rescued while suspected insertion core is not forced") {
        const auto ref = randomDna(120000, 20261119);
        auto q1 = mutated(ref, 1234, 30);
        auto q2 = mutated(ref, 5678, 40);
        q1.insert(60000, randomDna(30000, 4567));
        const std::vector<homology_plan::Sequence> reads = {
            {"reference", ref}, {"carrier", q1}, {"control", q2}};
        homology_plan::Config strict;
        strict.strict_candidate_guard = true;
        auto before = homology_plan::build(reads, strict);
        const auto previous = homology_plan::asTsv(before);
        strict.recover_indel_flanks = true;
        const auto after = homology_plan::build(reads, strict);
        const auto recovered = homology_plan::asTsv(after);
        CHECK(after.boundaries.size() > before.boundaries.size());
        CHECK(recovered.find("PR3D.2C: evidence-bracketed") != std::string::npos);
        CHECK(recovered.find("UNRESOLVED") != std::string::npos);
        CHECK(recovered.find("CANDIDATE_ONLY") != std::string::npos);
        const auto a = homology_plan::build(reads, strict);
        CHECK(homology_plan::asTsv(a) == recovered);
        const auto count_accepted = [](const std::string& tsv) {
            std::istringstream in(tsv);
            std::string line;
            std::uint64_t bp = 0;
            while (std::getline(in, line)) {
                if (line.empty() || line[0] == '#') continue;
                std::vector<std::string> fields;
                std::istringstream row(line);
                while (std::getline(row, line, '\t')) fields.push_back(line);
                if (fields.size() == 10 && fields[1] == "carrier" &&
                    fields[9] == "CANDIDATE_ONLY")
                    bp += std::stoull(fields[3]) - std::stoull(fields[2]);
            }
            return bp;
        };
        CHECK(count_accepted(recovered) > count_accepted(previous));
        CHECK(count_accepted(recovered) < ref.size());
        CHECK_THROWS_AS(homology_plan::build(reads, {.recover_indel_flanks=true}),
                        std::invalid_argument);
    }
    TEST_CASE("PR3D.2C disabled preserves exact strict and nonstrict diagnostics") {
        const auto ref = randomDna(50000, 20261120);
        const auto q = mutated(ref, 1234, 25);
        const std::vector<homology_plan::Sequence> reads = {{"r", ref}, {"q", q}};
        const auto old = homology_plan::asTsv(homology_plan::build(reads));
        CHECK(old.find("PR3D.2C") == std::string::npos);
        homology_plan::Config strict;
        strict.strict_candidate_guard = true;
        const auto old_strict = homology_plan::asTsv(homology_plan::build(reads, strict));
        CHECK(old_strict.find("PR3D.2C") == std::string::npos);
        CHECK(old_strict.find("PR3D.2B") != std::string::npos);
        strict.recover_indel_flanks = true;
        const auto recovery = homology_plan::asTsv(homology_plan::build(reads, strict));
        CHECK(recovery.find("PR3D.2C") != std::string::npos);
        CHECK(old_strict != recovery);
    }
}
