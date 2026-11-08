#include "doctest/doctest.h"
#include "align/ordered_sam.hpp"

#include <array>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
namespace fs = std::filesystem;

struct Sandbox {
    fs::path dir;
    Sandbox() : dir(fs::temp_directory_path() / ("halign4_ordered_sam_" +
        std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()))) {
        fs::create_directories(dir);
    }
    ~Sandbox() {
        std::error_code ignored;
        fs::remove_all(dir, ignored);
    }
};

void writeSam(const fs::path& path, const std::vector<std::uint64_t>& ordinals) {
    auto writer = seq_io::SeqWriter::Sam(path);
    writer.writeSamHeader("@HD\tVN:1.6\tSO:unknown");
    for (const auto ordinal : ordinals) {
        seq_io::SamRecord r;
        r.qname = "Q" + std::to_string(ordinal);
        r.rname = "consensus";
        r.pos = 1;
        r.mapq = 60;
        r.cigar = "4M";
        r.seq = "ACGT";
        align::ordered_sam::addOrdinal(r, ordinal);
        writer.writeSam(r);
    }
    writer.flush();
}

std::vector<std::string> mergedNames(const std::vector<fs::path>& paths) {
    align::ordered_sam::Reader reader(paths);
    seq_io::SamRecord rec;
    std::vector<std::string> names;
    while (reader.next(rec)) names.push_back(rec.qname);
    return names;
}
}  // namespace

TEST_SUITE("ordered_sam") {
    TEST_CASE("SAM ordinal tags round-trip with 64-bit values") {
        seq_io::SamRecord rec;
        rec.qname = "SEQ";
        rec.opt = "NM:i:0";
        align::ordered_sam::addOrdinal(rec, 4294967298ULL);
        CHECK(align::ordered_sam::getOrdinal(rec) == 4294967298ULL);
        CHECK(rec.opt.find("NM:i:0") != std::string::npos);
        rec.opt = "NM:i:0";
        CHECK_THROWS_AS(align::ordered_sam::getOrdinal(rec), std::runtime_error);
        rec.opt = "XI:Z:abc";
        CHECK_THROWS_AS(align::ordered_sam::getOrdinal(rec), std::runtime_error);
        rec.opt = "XI:Z:1\tXI:Z:2";
        CHECK_THROWS_AS(align::ordered_sam::getOrdinal(rec), std::runtime_error);
        rec.opt = "XI:Z:18446744073709551616";
        CHECK_THROWS_AS(align::ordered_sam::getOrdinal(rec), std::runtime_error);
    }

    TEST_CASE("global order is identical for one, two and three SAM streams") {
        Sandbox d;
        const auto a = d.dir / "a.sam";
        const auto b = d.dir / "b.sam";
        const auto c = d.dir / "c.sam";
        writeSam(a, {0, 3, 6, 9});
        writeSam(b, {1, 4, 7});
        writeSam(c, {2, 5, 8});
        const auto expected = std::vector<std::string>{
            "Q0", "Q1", "Q2", "Q3", "Q4", "Q5", "Q6", "Q7", "Q8", "Q9"};
        CHECK(mergedNames({a, b, c}) == expected);
        writeSam(a, {0, 1, 2, 3, 4, 5, 6, 7, 8, 9});
        writeSam(b, {});
        CHECK(mergedNames({a, b}) == expected);
        CHECK(mergedNames({a}) == expected);
    }

    TEST_CASE("missing, duplicate and nonmonotone ordinals fail closed") {
        Sandbox d;
        const auto a = d.dir / "a.sam";
        const auto b = d.dir / "b.sam";
        CHECK_THROWS_AS(mergedNames({a}), std::runtime_error);
        writeSam(a, {0, 2});
        writeSam(b, {1, 2});
        CHECK_THROWS_AS(mergedNames({a, b}), std::runtime_error);
        writeSam(a, {0, 4, 3});
        writeSam(b, {});
        CHECK_THROWS_AS(mergedNames({a, b}), std::runtime_error);
    }
}
