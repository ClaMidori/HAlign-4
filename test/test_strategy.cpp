#include "doctest/doctest.h"
#include "align/strategy.hpp"

#include <chrono>
#include <filesystem>
#include <fstream>
#include <string>

#if defined(HALIGN4_HAVE_ZLIB) && HALIGN4_HAVE_ZLIB
#include <zlib.h>
#endif

TEST_SUITE("strategy") {
    TEST_CASE("scan FASTA statistics without loading the entire dataset") {
        namespace fs = std::filesystem;
        const auto file = fs::temp_directory_path() /
            ("halign4_stats_" + std::to_string(
                std::chrono::steady_clock::now().time_since_epoch().count()) + ".fa");
        {
            std::ofstream stream(file);
            stream << ">A\nACGT\n>B\nACGTACGT\n>C\nACGTAC\n";
        }
        const auto stats = align_strategy::scanFasta(file);
        CHECK(stats.sequence_count == 3);
        CHECK(stats.total_bases == 18);
        CHECK(stats.longest_bases == 8);
        CHECK(stats.shortest_bases == 4);
        std::error_code ec;
        fs::remove(file, ec);
    }

    TEST_CASE("large FASTA record is counted across fixed-size buffers") {
        namespace fs = std::filesystem;
        const auto file = fs::temp_directory_path() /
            ("halign4_stream_stats_" + std::to_string(
                std::chrono::steady_clock::now().time_since_epoch().count()) + ".fa");
        {
            std::ofstream out(file, std::ios::binary);
            out << ">long description\n";
            for (int i = 0; i < 1200; ++i) out << std::string(1000, 'A') << "\n";
            out << ">short\nACG\nTT\n";
        }
        const auto stats = align_strategy::scanFasta(file);
        CHECK(stats.sequence_count == 2);
        CHECK(stats.total_bases == 1200005);
        CHECK(stats.longest_bases == 1200000);
        CHECK(stats.shortest_bases == 5);
        std::error_code ec;
        fs::remove(file, ec);
    }

#if defined(HALIGN4_HAVE_ZLIB) && HALIGN4_HAVE_ZLIB
    TEST_CASE("compressed gzip FASTA has the same statistics") {
        namespace fs = std::filesystem;
        const auto file = fs::temp_directory_path() /
            ("halign4_gz_stats_" + std::to_string(
                std::chrono::steady_clock::now().time_since_epoch().count()) + ".fa.gz");
        gzFile out = gzopen(file.string().c_str(), "wb");
        REQUIRE(out != nullptr);
        const char* fasta = ">A\nACGT\n>B\nA\n";
        CHECK(gzwrite(out, fasta, 13) == 13);
        CHECK(gzclose(out) == Z_OK);
        const auto stats = align_strategy::scanFasta(file);
        CHECK(stats.sequence_count == 2);
        CHECK(stats.total_bases == 5);
        CHECK(stats.longest_bases == 4);
        CHECK(stats.shortest_bases == 1);
        std::error_code ec;
        fs::remove(file, ec);
    }
#endif

    TEST_CASE("scanner rejects empty DNA and invalid FASTA structure") {
        namespace fs = std::filesystem;
        const auto file = fs::temp_directory_path() /
            ("halign4_stream_invalid_" + std::to_string(
                std::chrono::steady_clock::now().time_since_epoch().count()) + ".fa");
        {
            std::ofstream out(file);
            out << ">A\nACGT\n>B\n";
        }
        CHECK_THROWS(align_strategy::scanFasta(file));
        {
            std::ofstream out(file);
            out << "ACGT\n>A\nACGT\n";
        }
        CHECK_THROWS(align_strategy::scanFasta(file));
        {
            std::ofstream out(file);
            out << ">\nACGT\n";
        }
        CHECK_THROWS(align_strategy::scanFasta(file));
        std::error_code ec;
        fs::remove(file, ec);
    }

    TEST_CASE("hierarchical merged FASTA scanner detects unequal row widths") {
        namespace fs = std::filesystem;
        const auto file = fs::temp_directory_path() /
            ("halign4_hierarchical_rect_" + std::to_string(
                std::chrono::steady_clock::now().time_since_epoch().count()) + ".fa");
        {
            std::ofstream out(file);
            out << ">reference\nACGT--\n>insertion\nACGTTT\n";
        }
        const auto valid = align_strategy::scanRectangularFasta(file);
        CHECK(valid.sequence_count == 2);
        CHECK(valid.longest_bases == 6);
        CHECK(valid.contains_gaps);
        {
            std::ofstream out(file);
            out << ">reference\nACGT\n>insertion\nACGTTT\n";
        }
        CHECK_THROWS_WITH_AS(align_strategy::scanRectangularFasta(file),
                             doctest::Contains("Non-rectangular"), std::runtime_error);
        const auto unaligned = align_strategy::scanFasta(file);
        CHECK(unaligned.sequence_count == 2);
        CHECK(unaligned.longest_bases == 6);
        std::error_code ec;
        fs::remove(file, ec);
    }

    TEST_CASE("native direct routing is bounded and will not drop FASTA records") {
        align_strategy::FastaStats stats{3, 18, 8, 4};
        const align_strategy::Budget budget{};
        CHECK(align_strategy::chooseNativeDirect(stats, budget, 1000, false, false).native_direct);
        CHECK_FALSE(align_strategy::chooseNativeDirect(stats, budget, 2, false, false).native_direct);
        CHECK_FALSE(align_strategy::chooseNativeDirect(stats, budget, 1000, true, false).native_direct);
        CHECK_FALSE(align_strategy::chooseNativeDirect(stats, budget, 1000, false, true).native_direct);
        stats.longest_bases = 25001;
        CHECK_FALSE(align_strategy::chooseNativeDirect(stats, budget, 1000, false, false).native_direct);
        stats.longest_bases = 8;
        stats.total_bases = 100001;
        CHECK_FALSE(align_strategy::chooseNativeDirect(stats, budget, 1000, false, false).native_direct);
        stats.total_bases = 18;
        stats.sequence_count = 65;
        CHECK_FALSE(align_strategy::chooseNativeDirect(stats, budget, 1000, false, false).native_direct);
    }
}
