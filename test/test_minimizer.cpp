#include <doctest/doctest.h>

#include <seed.h>

#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <string>
#include <vector>

// ================================================================
// This file plugs extractMinimizer performance tests into doctest:
// - skipped by default (to avoid long CI/normal test runs)
// - explicitly enabled via env var HALIGN4_RUN_PERF=1
// - no strict perf assertions; it prints timing/throughput for manual comparison
// ================================================================

static bool perfEnabled() {
    const char* v = std::getenv("HALIGN4_RUN_PERF");
    return (v != nullptr) && (*v != '\0') && (std::string(v) != "0");
}
static bool shouldSkipPerf() { return !perfEnabled(); }

static std::string makeRandomDna(std::size_t len, std::uint32_t seed)
{
    static constexpr char bases[4] = {'A', 'C', 'G', 'T'};
    std::uint64_t x = seed;

    // Simple xorshift64* PRNG, lighter than mt19937 and sufficient for generating test data
    auto next = [&]() {
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        return x * 2685821657736338717ULL;
    };

    std::string s;
    s.resize(len);
    for (std::size_t i = 0; i < len; ++i) {
        s[i] = bases[static_cast<std::size_t>(next() & 3ULL)];
    }
    return s;
}

TEST_SUITE("minimizer" * doctest::skip(shouldSkipPerf()))
{
    TEST_CASE("extractMinimizer - throughput")
    {
        // Parameters can be overridden via env vars for easier tuning locally
        // HALIGN4_MINIMIZER_SEQ_LEN=10000
        // HALIGN4_MINIMIZER_NUM_SEQS=200
        // HALIGN4_MINIMIZER_ROUNDS=5
        auto getenv_u64 = [](const char* name, std::uint64_t defv) {
            if (const char* p = std::getenv(name); p && *p) return static_cast<std::uint64_t>(std::strtoull(p, nullptr, 10));
            return defv;
        };

        const std::size_t seq_len  = static_cast<std::size_t>(getenv_u64("HALIGN4_MINIMIZER_SEQ_LEN",  30000));
        const std::size_t num_seqs = static_cast<std::size_t>(getenv_u64("HALIGN4_MINIMIZER_NUM_SEQS", 100000));
        const std::size_t rounds   = static_cast<std::size_t>(getenv_u64("HALIGN4_MINIMIZER_ROUNDS",   1));

        // parameters: k/w
        const std::size_t k = static_cast<std::size_t>(getenv_u64("HALIGN4_MINIMIZER_K", 15));
        const std::size_t w = static_cast<std::size_t>(getenv_u64("HALIGN4_MINIMIZER_W", 10));

        MESSAGE("seq_len=" << seq_len << " num_seqs=" << num_seqs << " rounds=" << rounds << " k=" << k << " w=" << w);

        std::vector<std::string> seqs;
        seqs.reserve(num_seqs);
        for (std::size_t i = 0; i < num_seqs; ++i) {
            seqs.emplace_back(makeRandomDna(seq_len, static_cast<std::uint32_t>(1234 + i)));
        }

        // Warm-up: avoid first-run cache/branch predictor effects from dominating
        std::uint64_t checksum = 0;
        for (const auto& s : seqs) {
            auto mz = minimizer::extractMinimizer(s, k, w, false);
            checksum += static_cast<std::uint64_t>(mz.size());
            if (!mz.empty()) checksum ^= mz.front().hash();
        }

        const auto t0 = std::chrono::steady_clock::now();

        std::uint64_t total_minimizers = 0;
        for (std::size_t r = 0; r < rounds; ++r) {
            for (const auto& s : seqs) {
                auto mz = minimizer::extractMinimizer(s, k, w, false);
                total_minimizers += static_cast<std::uint64_t>(mz.size());
                if (!mz.empty()) checksum ^= mz.back().hash();
            }
        }

        const auto t1 = std::chrono::steady_clock::now();
        const double sec = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0).count();

        const double total_bp = static_cast<double>(seq_len) * static_cast<double>(num_seqs) * static_cast<double>(rounds);
        const double bp_per_s = total_bp / sec;
        const double seq_per_s = (static_cast<double>(num_seqs) * static_cast<double>(rounds)) / sec;

        MESSAGE("elapsed_s=" << sec);
        MESSAGE("throughput_bp_per_s=" << bp_per_s);
        MESSAGE("throughput_seq_per_s=" << seq_per_s);
        MESSAGE("total_minimizers=" << total_minimizers);
        MESSAGE("checksum=" << checksum);

        // Do not make hard perf assertions to avoid CI instability due to machine/load differences.
        CHECK(sec > 0.0);
    }
}
