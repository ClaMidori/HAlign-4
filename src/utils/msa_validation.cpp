#include "msa_validation.hpp"
#include "utils.h"

#include <cstdint>
#include <stdexcept>
#include <string>
#include <unordered_map>

namespace {
struct Signature {
    std::uint64_t bases = 0;
    std::uint64_t hash1 = 14695981039346656037ULL;
    std::uint64_t hash2 = 7809847782465536321ULL;
};

Signature signatureOf(const std::string& sequence) {
    Signature sig;
    for (const unsigned char ch : sequence) {
        if (ch == '-') continue;
        ++sig.bases;
        sig.hash1 = (sig.hash1 ^ ch) * 1099511628211ULL;
        sig.hash2 = (sig.hash2 ^ ch) * 14029467366897019727ULL;
    }
    return sig;
}

struct Expected {
    Signature signature;
    bool seen = false;
};
}

namespace msa {
ValidationStats validateFastaAlignment(const std::filesystem::path& input,
                                       const std::filesystem::path& aligned) {
    file_io::requireRegularFile(input, "unaligned FASTA");
    file_io::requireRegularFile(aligned, "aligned FASTA");
    if (file_io::isEmpty(aligned)) {
        throw std::runtime_error("aligned FASTA is empty: " + aligned.string());
    }

    std::unordered_map<std::string, Expected> expected;
    {
        seq_io::KseqReader reader(input);
        seq_io::SeqRecord record;
        while (reader.next(record)) {
            if (record.id.empty()) throw std::runtime_error("empty input FASTA ID");
            auto [it, added] = expected.emplace(record.id, Expected{signatureOf(record.seq), false});
            if (!added) throw std::runtime_error("duplicate input FASTA ID: " + record.id);
        }
    }
    if (expected.empty()) throw std::runtime_error("unaligned FASTA contains no sequences");

    ValidationStats stats;
    bool first = true;
    seq_io::KseqReader reader(aligned);
    seq_io::SeqRecord record;
    while (reader.next(record)) {
        const auto it = expected.find(record.id);
        if (it == expected.end()) {
            throw std::runtime_error("unexpected aligned FASTA ID: " + record.id);
        }
        if (it->second.seen) {
            throw std::runtime_error("duplicate aligned FASTA ID: " + record.id);
        }
        it->second.seen = true;

        const auto width = static_cast<std::uint64_t>(record.seq.size());
        if (first) {
            stats.columns = width;
            first = false;
        } else if (width != stats.columns) {
            throw std::runtime_error("non-rectangular MSA at ID: " + record.id);
        }

        const auto actual = signatureOf(record.seq);
        const auto& original = it->second.signature;
        if (actual.bases != original.bases ||
            actual.hash1 != original.hash1 ||
            actual.hash2 != original.hash2) {
            throw std::runtime_error("bases changed, missing or reordered at ID: " + record.id);
        }
        stats.bases += actual.bases;
        ++stats.records;
    }
    if (stats.records != expected.size()) {
        throw std::runtime_error("aligned FASTA has missing sequence IDs");
    }
    if (stats.columns == 0) {
        throw std::runtime_error("aligned FASTA contains no alignment columns");
    }
    return stats;
}
}  // namespace msa
