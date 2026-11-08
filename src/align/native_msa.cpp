#include "align/native_msa.hpp"

#include <abpoa.h>

#include <cstdint>
#include <cstdio>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

namespace native_msa {
namespace {

std::uint8_t encode(char base) {
    switch (base) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        case 'N': return 4;
        default: throw std::invalid_argument(
            "NativeMSA currently requires uppercase DNA A/C/G/T/N without gaps");
    }
}

char decode(std::uint8_t base) {
    switch (base) {
        case 0: return 'A';
        case 1: return 'C';
        case 2: return 'G';
        case 3: return 'T';
        case 4: return 'N';
        case 5: return '-';
        default: throw std::runtime_error("abPOA returned an invalid MSA base code");
    }
}

std::string ungap(const std::string& sequence) {
    std::string value;
    value.reserve(sequence.size());
    for (char c : sequence) {
        if (c != '-') value.push_back(c);
    }
    return value;
}

struct AbpoaDeleter {
    void operator()(abpoa_t* ptr) const noexcept { if (ptr) abpoa_free(ptr); }
};
struct ParaDeleter {
    void operator()(abpoa_para_t* ptr) const noexcept { if (ptr) abpoa_free_para(ptr); }
};
struct FileDeleter {
    void operator()(FILE* ptr) const noexcept { if (ptr) std::fclose(ptr); }
};

}  // namespace

Alignment align(const std::vector<Sequence>& sequences, const Options& options) {
    if (sequences.empty()) throw std::invalid_argument("NativeMSA requires sequences");
    if (options.max_sequences == 0 || options.max_sequence_bases == 0 ||
        options.max_total_bases == 0 || options.match <= 0 || options.mismatch <= 0) {
        throw std::invalid_argument("NativeMSA options must be positive");
    }
    if (sequences.size() > options.max_sequences ||
        sequences.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw std::invalid_argument("NativeMSA sequence count exceeds the budget or C API limits");
    }

    std::unordered_set<std::string> ids;
    std::size_t total = 0;
    std::vector<std::vector<std::uint8_t>> encoded;
    encoded.reserve(sequences.size());
    for (const auto& seq : sequences) {
        if (seq.id.empty() || !ids.insert(seq.id).second) {
            throw std::invalid_argument("NativeMSA requires non-empty, unique FASTA IDs");
        }
        if (seq.bases.empty() || seq.bases.size() > options.max_sequence_bases ||
            seq.bases.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
            seq.bases.size() > options.max_total_bases - total) {
            throw std::invalid_argument("NativeMSA sequence length exceeds configured budget");
        }
        total += seq.bases.size();
        encoded.emplace_back();
        auto& row = encoded.back();
        row.reserve(seq.bases.size());
        for (char base : seq.bases) row.push_back(encode(base));
    }

    if (sequences.size() == 1) return Alignment{sequences, sequences[0].bases.size()};

    std::vector<int> lengths;
    std::vector<std::uint8_t*> pointers;
    lengths.reserve(encoded.size());
    pointers.reserve(encoded.size());
    for (auto& row : encoded) {
        lengths.push_back(static_cast<int>(row.size()));
        pointers.push_back(row.data());
    }

    std::unique_ptr<abpoa_t, AbpoaDeleter> ab(abpoa_init());
    std::unique_ptr<abpoa_para_t, ParaDeleter> params(abpoa_init_para());
    if (!ab || !params) throw std::runtime_error("abPOA initialization failed");
    params->align_mode = ABPOA_GLOBAL_MODE;
    params->match = options.match;
    params->mismatch = options.mismatch;
    params->out_msa = 1;
    params->out_cons = 0;
    params->progressive_poa = 0;
    params->sort_input_seq = 0;
    // abPOA v1.5.6: CLI -S sets disable_seeding=0, not 1.
    params->disable_seeding = options.seeding ? 0 : 1;
    abpoa_post_set_para(params.get());

    // abpoa_msa() produces its in-memory MSA and ALSO writes to FILE* when
    // out_msa=1. Keep that output out of stdout, and remove the file on close.
    std::unique_ptr<FILE, FileDeleter> discard(std::tmpfile());
    if (!discard) throw std::runtime_error("NativeMSA failed to create temporary output stream");
    const int rc = abpoa_msa(ab.get(), params.get(), static_cast<int>(sequences.size()),
                             nullptr, lengths.data(), pointers.data(), nullptr, discard.get());
    if (rc != 0 || !ab->abc || ab->abc->n_seq != static_cast<int>(sequences.size()) ||
        ab->abc->msa_len <= 0 || !ab->abc->msa_base) {
        throw std::runtime_error("abPOA did not return a valid multiple alignment");
    }

    Alignment result;
    result.columns = static_cast<std::size_t>(ab->abc->msa_len);
    result.rows.reserve(sequences.size());
    for (std::size_t i = 0; i < sequences.size(); ++i) {
        if (!ab->abc->msa_base[i]) throw std::runtime_error("abPOA returned a null MSA row");
        std::string bases;
        bases.reserve(result.columns);
        for (std::size_t j = 0; j < result.columns; ++j) {
            bases.push_back(decode(ab->abc->msa_base[i][j]));
        }
        if (ungap(bases) != sequences[i].bases) {
            throw std::runtime_error("NativeMSA altered or reordered bases for: " + sequences[i].id);
        }
        result.rows.push_back(Sequence{sequences[i].id, std::move(bases)});
    }
    return result;
}

}  // namespace native_msa
