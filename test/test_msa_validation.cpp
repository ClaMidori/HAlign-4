#include <doctest/doctest.h>

#include "msa_validation.hpp"
#include "preprocess.h"

#include <chrono>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>

namespace {
namespace fs = std::filesystem;

struct Fixture {
    fs::path dir;
    Fixture() {
        const auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
        dir = fs::temp_directory_path() / ("halign4_msa_validation_" + std::to_string(stamp));
        fs::create_directories(dir);
    }
    ~Fixture() { std::error_code ec; fs::remove_all(dir, ec); }
    fs::path write(const std::string& filename, const std::string& content) const {
        auto path = dir / filename;
        std::ofstream f(path);
        f << content;
        f.close();
        return path;
    }
};
}

TEST_SUITE("msa_validation") {
    TEST_CASE("valid alignment, reordered FASTA IDs, gaps") {
        Fixture fx;
        auto in = fx.write("input.fa", ">S1\nACGT\n>S2\nAGT\n");
        auto out = fx.write("msa.fa", ">S2\nA-GT\n>S1\nACGT\n");
        auto stats = msa::validateFastaAlignment(in, out);
        CHECK(stats.records == 2);
        CHECK(stats.columns == 4);
        CHECK(stats.bases == 7);
    }
    TEST_CASE("rejects nonrectangular alignment") {
        Fixture fx;
        auto in = fx.write("input.fa", ">S1\nACGT\n>S2\nAGT\n");
        auto out = fx.write("msa.fa", ">S1\nACGT\n>S2\nAGT\n");
        CHECK_THROWS_AS(msa::validateFastaAlignment(in, out), std::runtime_error);
    }
    TEST_CASE("rejects changed bases, missing sequence and duplicate IDs") {
        Fixture fx;
        auto in = fx.write("input.fa", ">S1\nACGT\n>S2\nAGT\n");
        auto changed = fx.write("changed.fa", ">S1\nACGA\n>S2\nA-GT\n");
        auto missing = fx.write("missing.fa", ">S1\nACGT\n");
        auto duplicate = fx.write("duplicate.fa", ">S1\nACGT\n>S1\nACGT\n");
        CHECK_THROWS_AS(msa::validateFastaAlignment(in, changed), std::runtime_error);
        CHECK_THROWS_AS(msa::validateFastaAlignment(in, missing), std::runtime_error);
        CHECK_THROWS_AS(msa::validateFastaAlignment(in, duplicate), std::runtime_error);
    }
    TEST_CASE("rejects empty alignment") {
        Fixture fx;
        auto in = fx.write("input.fa", ">S1\nACGT\n");
        auto empty = fx.write("empty.fa", "");
        CHECK_THROWS_AS(msa::validateFastaAlignment(in, empty), std::runtime_error);
    }
    TEST_CASE("external command failures propagate and stale output is deleted") {
        Fixture fx;
        auto in = fx.write("input.fa", ">S1\nACGT\n>S2\nACGT\n");
        auto out = fx.write("output.fa", ">S1\nACGT\n>S2\nACGT\n");
        CHECK_THROWS_AS(
            alignConsensusSequence(in, out, "false # {input} {output}", 1),
            std::runtime_error);
        CHECK_FALSE(fs::exists(out));
        CHECK_THROWS_AS(
            alignConsensusSequence(in, out, "true # {input} {output}", 1),
            std::runtime_error);
        CHECK_FALSE(fs::exists(out));
    }
    TEST_CASE("external command valid output accepted") {
        Fixture fx;
        auto in = fx.write("input.fa", ">S1\nACGT\n>S2\nACGT\n");
        auto out = fx.dir / "output.fa";
        CHECK_NOTHROW(alignConsensusSequence(in, out, "cp {input} {output}", 1));
        CHECK(fs::exists(out));
    }
}
