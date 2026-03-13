#include <doctest/doctest.h>

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <string>

#include "utils.h"   // The declaration of seq_io::KseqReader / SeqRecord is required.

namespace fs = std::filesystem;

// ------------------------- perf gating -------------------------

// perf switch: you already have this in test_consensus.cpp, but we put it here again for self-containment of this test file.
static bool perfEnabled() {
    const char* v = std::getenv("HALIGN4_RUN_PERF");
    return (v != nullptr) && (*v != '\0') && (std::string(v) != "0");
}
static bool shouldSkipPerf() { return !perfEnabled(); }

// ------------------------- helpers -------------------------

static fs::path makeTempDir(std::string_view name) {
    fs::path base;
    if (const char* p = std::getenv("HALIGN4_PERF_DIR"); p && *p) {
        base = fs::path(p);
    } else {
        base = fs::current_path(); // To avoid situations where /tmp is tmpfs or space is limited, please refer to the relevant documentation.
    }

    fs::path dir = base / std::string(name);
    std::error_code ec;
    fs::remove_all(dir, ec);
    fs::create_directories(dir, ec);
    REQUIRE_MESSAGE(!ec, "cannot create temp dir: " << dir.string() << " (" << ec.message() << ")");
    return dir;
}

// Writing FASTA: Each sequence is written on a single line to avoid unnecessary overhead; 1000x30000 is approximately 30MB in size.
static void writeAlignedFastaSingleLine(const fs::path& p, std::size_t n_seqs, std::size_t len) {
    std::ofstream ofs(p, std::ios::binary);
    REQUIRE_MESSAGE(ofs.good(), "cannot write fasta: " << p.string());

    std::string seq(len, 'A');
    seq.push_back('\n');

    for (std::size_t i = 0; i < n_seqs; ++i) {
        char hdr[64];
        int hn = std::snprintf(hdr, sizeof(hdr), ">s%zu\n", i);
        ofs.write(hdr, hn);
        ofs.write(seq.data(), (std::streamsize)seq.size());
        REQUIRE_MESSAGE(ofs.good(), "write failed at i=" << i << " (disk full?)");
    }
    ofs.flush();
}

static double toMiB(double bytes) { return bytes / (1024.0 * 1024.0); }

// ------------------------- tests -------------------------

TEST_SUITE("read_fasta")
{
    TEST_CASE("KseqReader - smoke read small fasta") {
        auto dir = makeTempDir("halign4_tests_read_smoke");
        fs::path in = dir / "small.fasta";

        {
            std::ofstream ofs(in, std::ios::binary);
            REQUIRE(ofs.good());
            ofs << ">a\nACGTNn\n>b\nAAAA\n>c\nTTNN\n";
        }

        seq_io::KseqReader r(in);
        seq_io::SeqRecord rec;

        std::size_t count = 0;
        while (r.next(rec)) {
            ++count;
            CHECK(!rec.id.empty());
            CHECK(!rec.seq.empty());

            if (rec.id == "a") {
                CHECK(rec.n_num == 2);
            } else if (rec.id == "b") {
                CHECK(rec.n_num == 0);
            } else if (rec.id == "c") {
                CHECK(rec.n_num == 2);
            }
        }
        CHECK(count == 3);
    }
}
TEST_SUITE("read_fasta")
{
    // Performance testing: Execution only occurs when HALIGN4_RUN_PERF=1.
    TEST_CASE("read time: n=10000 len=30000 (kseq)")
    {
        constexpr std::size_t N   = 10000;
        constexpr std::size_t LEN = 30000;

        auto dir = makeTempDir("halign4_tests_read_perf");
        fs::path in = dir / "aligned_1000_30000.fasta";
        fs::path out = dir / "aligned_1000_30000_out.fasta";
        seq_io::SeqWriter clean_writer(out);

        // 1) Preparing for input (excluding reading time)
        writeAlignedFastaSingleLine(in, N, LEN);

        std::error_code ec;
        const auto fsz = fs::file_size(in, ec);
        MESSAGE("input_ready: path=" << in.string()
                << " size_MiB=" << (ec ? 0.0 : toMiB((double)fsz)));

        // 2) Timing read
        auto t0 = std::chrono::steady_clock::now();

        seq_io::KseqReader r(in);
        seq_io::SeqRecord rec;

        std::size_t count = 0;
        std::size_t bad_len = 0;
        while (r.next(rec)) {
            ++count;
            if (rec.seq.size() != LEN) ++bad_len;
            seq_io::cleanSequence(rec.seq);
        }

        auto t1 = std::chrono::steady_clock::now();
        const double sec = std::chrono::duration<double>(t1 - t0).count();

        CHECK(count == N);
        CHECK(bad_len == 0);

        // 3) Output throughput
        const double mib = ec ? 0.0 : toMiB((double)fsz);
        const double mibps = (sec > 0.0) ? (mib / sec) : 0.0;

        MESSAGE("read_done: n=" << count
                << " len=" << LEN
                << " time_s=" << sec
                << " throughput_MiBps=" << mibps);

        // 4) Clean up (optional)
        fs::remove_all(dir, ec);
    }

} // TEST_SUITE(fasta_read)
