#include "doctest/doctest.h"
#include "align/profile_projection.hpp"
#include "align.h"

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <unordered_map>

namespace {
std::string ungapped(std::string_view row) {
    std::string result;
    for (const char c : row) if (c != '-') result += c;
    return result;
}
}

TEST_SUITE("profile_projection") {
    TEST_CASE("different consensus gap placements become one coordinate system") {
        const auto plan = align::profile_projection::makePlan(
            "A--CG-T", "A-C--GT", "ACGT");
        CHECK(plan.columns == 9);
        const auto ref = align::profile_projection::projectRow(
            "A--CG-T", plan.reference_to_common, plan.columns);
        const auto ins = align::profile_projection::projectRow(
            "A-C--GT", plan.insertion_to_common, plan.columns);
        CHECK(ref == ins);
        CHECK(ungapped(ref) == "ACGT");
        // Other rows retain every original base, even if their gap patterns differ.
        const auto ref_query = align::profile_projection::projectRow(
            "ATTCGGT", plan.reference_to_common, plan.columns);
        const auto insertion_query = align::profile_projection::projectRow(
            "ATCGGGT", plan.insertion_to_common, plan.columns);
        CHECK(ref_query.size() == plan.columns);
        CHECK(insertion_query.size() == plan.columns);
        CHECK(ungapped(ref_query) == "ATTCGGT");
        CHECK(ungapped(insertion_query) == "ATCGGGT");
    }

    TEST_CASE("consensus differences or incompatible row lengths are rejected") {
        CHECK_THROWS_AS(align::profile_projection::makePlan("ACGT", "ACGA"),
                        std::runtime_error);
        CHECK_THROWS_AS(align::profile_projection::makePlan("ACGT", "ACG"),
                        std::runtime_error);
        CHECK_THROWS_AS(align::profile_projection::makePlan("A-CGT", "ACGT", "ACGA"),
                        std::runtime_error);
        const auto plan = align::profile_projection::makePlan("A-CGT", "ACGT");
        CHECK_THROWS_AS(align::profile_projection::projectRow(
            "ACG", plan.reference_to_common, plan.columns), std::invalid_argument);
    }

    TEST_CASE("FASTA profiles are projected without losing DNA") {
        namespace fs = std::filesystem;
        const auto dir = fs::temp_directory_path() / ("halign4_profile_" +
            std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
        fs::create_directories(dir);
        const auto refs = dir / "refs.fa", ins = dir / "ins.fa";
        const auto projected_refs = dir / "refs_projected.fa";
        const auto projected_ins = dir / "ins_projected.fa";
        {
            std::ofstream stream(refs);
            stream << ">consensus\nA--CG-T\n>reference-1\nATTCGGT\n";
        }
        {
            std::ofstream stream(ins);
            stream << ">consensus\nA-C--GT\n>insertion-1\nATCGGGT\n";
        }
        const auto width = align::profile_projection::reconcileFastaProfiles(
            refs, ins, projected_refs, projected_ins, "ACGT");
        CHECK(width == 9);
        std::unordered_map<std::string, std::string> output;
        for (auto& path : {projected_refs, projected_ins}) {
            seq_io::KseqReader reader(path);
            seq_io::SeqRecord rec;
            while (reader.next(rec)) {
                CHECK(rec.seq.size() == width);
                if (rec.id == "consensus") {
                    if (!output.contains("consensus")) output["consensus"] = rec.seq;
                    CHECK(rec.seq == output["consensus"]);
                } else {
                    output[rec.id] = ungapped(rec.seq);
                }
            }
        }
        CHECK(output["reference-1"] == "ATTCGGT");
        CHECK(output["insertion-1"] == "ATCGGGT");
        fs::remove_all(dir);
    }

    TEST_CASE("user-named center is allowed only when its bases match the consensus") {
        namespace fs = std::filesystem;
        const auto dir = fs::temp_directory_path() / ("halign4_center_" +
            std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
        fs::create_directories(dir);
        const auto refs = dir / "refs.fa", ins = dir / "ins.fa";
        {
            std::ofstream stream(refs);
            stream << ">external-center\nACGT\n";
        }
        {
            std::ofstream stream(ins);
            stream << ">consensus\nAC-GT\n";
        }
        CHECK(align::profile_projection::reconcileFastaProfiles(
            refs, ins, dir / "out-r.fa", dir / "out-i.fa", "ACGT") == 5);
        {
            std::ofstream stream(refs);
            stream << ">external-center\nACGA\n";
        }
        CHECK_THROWS_AS(align::profile_projection::reconcileFastaProfiles(
            refs, ins, dir / "bad-r.fa", dir / "bad-i.fa", "ACGT"),
            std::runtime_error);
        fs::remove_all(dir);
    }

    TEST_CASE("non-rectangular source profile is not silently accepted") {
        namespace fs = std::filesystem;
        const auto dir = fs::temp_directory_path() / ("halign4_profile_invalid_" +
            std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
        fs::create_directories(dir);
        const auto refs = dir / "refs.fa", ins = dir / "ins.fa";
        {
            std::ofstream stream(refs);
            stream << ">consensus\nACGT\n>bad-row\nACG\n";
        }
        {
            std::ofstream stream(ins);
            stream << ">consensus\nACGT\n";
        }
        CHECK_THROWS_AS(align::profile_projection::reconcileFastaProfiles(
            refs, ins, dir / "out-r.fa", dir / "out-i.fa", "ACGT"),
            std::runtime_error);
        fs::remove_all(dir);
    }
}
