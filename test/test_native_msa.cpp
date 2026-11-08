#include "doctest/doctest.h"
#include "align/native_msa.hpp"
#include "preprocess.h"
#include "msa_validation.hpp"

#include <stdexcept>
#include <filesystem>
#include <fstream>
#include <chrono>
#include <string>
#include <vector>

namespace {
std::string ungap(std::string text) {
    std::string result;
    for (const char base : text) if (base != '-') result.push_back(base);
    return result;
}
}

TEST_SUITE("native_msa") {
    TEST_CASE("identical short DNA, with in-process abPOA") {
        const std::vector<native_msa::Sequence> input = {
            {"seq1", "ACGTACGTACGT"}, {"seq2", "ACGTACGTACGT"},
            {"seq3", "ACGTACGTACGT"}
        };
        const auto output = native_msa::align(input);
        REQUIRE(output.rows.size() == input.size());
        CHECK(output.columns == input[0].bases.size());
        for (std::size_t i = 0; i < input.size(); ++i) {
            CHECK(output.rows[i].id == input[i].id);
            CHECK(output.rows[i].bases == input[i].bases);
        }
    }

    TEST_CASE("indels, N and original row order survive") {
        const std::vector<native_msa::Sequence> input = {
            {"A", "ACGTACGTACGTACGT"},
            {"B", "ACGTACGTTTACGTACGT"},
            {"C", "ACGTACGACGT"},
            {"D", "ACGTNCGTACGTACGT"}
        };
        const auto output = native_msa::align(input);
        REQUIRE(output.rows.size() == input.size());
        CHECK(output.columns >= input[1].bases.size());
        for (std::size_t i = 0; i < input.size(); ++i) {
            CHECK(output.rows[i].id == input[i].id);
            CHECK(output.rows[i].bases.size() == output.columns);
            CHECK(ungap(output.rows[i].bases) == input[i].bases);
        }
    }

    TEST_CASE("reject duplicate IDs and invalid symbols without abPOA") {
        CHECK_THROWS_AS(native_msa::align({{"A", "ACGT"}, {"A", "ACGA"}}),
                        std::invalid_argument);
        CHECK_THROWS_AS(native_msa::align({{"A", "ACGT"}, {"B", "AcGT"}}),
                        std::invalid_argument);
        CHECK_THROWS_AS(native_msa::align({{"A", "ACGT"}, {"B", "AC-GT"}}),
                        std::invalid_argument);
        CHECK_THROWS_AS(native_msa::align({{"A", ""}, {"B", "ACGT"}}),
                        std::invalid_argument);
    }

    TEST_CASE("refuse oversized direct alignment") {
        native_msa::Options options;
        options.max_sequence_bases = 5;
        CHECK_THROWS_AS(native_msa::align({{"A", "ACGTAC"}, {"B", "ACGT"}}, options),
                        std::invalid_argument);
    }
}

TEST_SUITE("native_msa") {
    TEST_CASE("NativeMSA FASTA route executes without external command") {
        namespace fs = std::filesystem;
        const auto tmp = fs::temp_directory_path() /
            ("halign4_native_route_" + std::to_string(
                std::chrono::steady_clock::now().time_since_epoch().count()));
        fs::create_directories(tmp);
        const auto input = tmp / "input.fa";
        const auto output = tmp / "output.fa";
        {
            std::ofstream f(input);
            f << ">A\nACGTACGTACGT\n>B\nACGTTACGTACGT\n>C\nACGTACGTACGT\n";
        }
        CHECK_NOTHROW(alignConsensusSequence(input, output, "native", 1));
        CHECK_NOTHROW(msa::validateFastaAlignment(input, output));
        std::error_code ec;
        fs::remove_all(tmp, ec);
    }

    TEST_CASE("NativeMSA FASTA route rejects oversize and removes stale output") {
        namespace fs = std::filesystem;
        const auto tmp = fs::temp_directory_path() /
            ("halign4_native_oversize_" + std::to_string(
                std::chrono::steady_clock::now().time_since_epoch().count()));
        fs::create_directories(tmp);
        const auto input = tmp / "input.fa";
        const auto output = tmp / "output.fa";
        {
            std::ofstream f(input);
            f << ">A\n" << std::string(25001, 'A') << "\n>B\nACGT\n";
        }
        {
            std::ofstream f(output);
            f << "stale\n";
        }
        CHECK_THROWS(alignConsensusSequence(input, output, "native", 1));
        CHECK_FALSE(fs::exists(output));
        std::error_code ec;
        fs::remove_all(tmp, ec);
    }
}

TEST_SUITE("native_msa") {
    TEST_CASE("native direct preprocessing preserves input order; hierarchical selection stays length-first") {
        namespace fs = std::filesystem;
        const auto tmp = fs::temp_directory_path() /
            ("halign4_input_order_" + std::to_string(
                std::chrono::steady_clock::now().time_since_epoch().count()));
        fs::create_directories(tmp);
        const auto input = tmp / "input.fasta";
        {
            std::ofstream f(input);
            f << ">A\nACGT\n>B\nACGTACGT\n>C\nACGTAC\n";
        }
        const auto get_ids = [](const fs::path& path) {
            std::vector<std::string> ids;
            auto reader = seq_io::openKseqReader(path);
            seq_io::SeqRecord record;
            while (reader->next(record)) ids.push_back(record.id);
            return ids;
        };

        const auto native_work = tmp / "native";
        REQUIRE(preprocessInputFasta(input.string(), native_work.string(), 10, true) == 3);
        const auto native_input = native_work / WORKDIR_DATA / DATA_CLEAN / CLEAN_CONS_UNALIGNED;
        const std::vector<std::string> expected_native = {"A", "B", "C"};
        CHECK(get_ids(native_input) == expected_native);
        const auto native_output = tmp / "native_output.fasta";
        CHECK_NOTHROW(alignConsensusSequence(native_input, native_output, "native", 1));
        CHECK(get_ids(native_output) == expected_native);
        CHECK_NOTHROW(msa::validateFastaAlignment(native_input, native_output));

        const auto hierarchical_work = tmp / "hierarchical";
        REQUIRE(preprocessInputFasta(input.string(), hierarchical_work.string(), 10, false) == 3);
        const auto hierarchical_input = hierarchical_work / WORKDIR_DATA / DATA_CLEAN / CLEAN_CONS_UNALIGNED;
        const std::vector<std::string> expected_hierarchical = {"B", "C", "A"};
        CHECK(get_ids(hierarchical_input) == expected_hierarchical);
        std::error_code ec;
        fs::remove_all(tmp, ec);
    }
}
