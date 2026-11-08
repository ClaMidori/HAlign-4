#include "align/profile_projection.hpp"
#include "utils.h"

#include <algorithm>
#include <cctype>
#include <limits>
#include <stdexcept>
#include <string>

namespace align::profile_projection {
namespace {

char normalized(char c) {
    return static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
}

void appendMappedColumns(std::vector<std::size_t>& mapping,
                         std::size_t n, std::size_t start) {
    if (n > std::numeric_limits<std::size_t>::max() - start)
        throw std::overflow_error("Profile projection exceeds addressable columns");
    for (std::size_t j = 0; j < n; ++j) mapping.push_back(start + j);
}

void projectProfile(const std::filesystem::path& source,
                    const std::filesystem::path& destination,
                    const std::vector<std::size_t>& map,
                    std::size_t common_columns) {
    seq_io::KseqReader reader(source);
    seq_io::SeqWriter writer(destination, 80);
    seq_io::SeqRecord record;
    std::size_t count = 0;
    while (reader.next(record)) {
        if (record.seq.size() != map.size())
            throw std::runtime_error("Non-rectangular profile: " + source.string() +
                                     " sequence " + record.id);
        record.seq = projectRow(record.seq, map, common_columns);
        writer.writeFasta(record);
        ++count;
    }
    if (count == 0) throw std::runtime_error("Empty MSA profile: " + source.string());
    writer.flush();
}

}  // namespace

Plan makePlan(std::string_view reference_anchor,
              std::string_view insertion_anchor,
              std::string_view expected_ungapped_anchor) {
    if (reference_anchor.empty() || insertion_anchor.empty())
        throw std::invalid_argument("Cannot merge MSA profiles with empty anchors");

    Plan plan;
    plan.reference_to_common.reserve(reference_anchor.size());
    plan.insertion_to_common.reserve(insertion_anchor.size());

    std::size_t r = 0, i = 0, emitted_bases = 0;
    while (r < reference_anchor.size() || i < insertion_anchor.size()) {
        const std::size_t r_start = r, i_start = i;
        while (r < reference_anchor.size() && reference_anchor[r] == '-') ++r;
        while (i < insertion_anchor.size() && insertion_anchor[i] == '-') ++i;

        const std::size_t rgaps = r - r_start;
        const std::size_t igaps = i - i_start;
        const std::size_t widest = std::max(rgaps, igaps);
        if (widest > std::numeric_limits<std::size_t>::max() - plan.columns)
            throw std::overflow_error("Merged MSA width overflows size_t");
        appendMappedColumns(plan.reference_to_common, rgaps, plan.columns);
        appendMappedColumns(plan.insertion_to_common, igaps, plan.columns);
        plan.columns += widest;

        // Both profiles must end at the same residue boundary.
        if (r == reference_anchor.size() || i == insertion_anchor.size()) {
            if (r != reference_anchor.size() || i != insertion_anchor.size())
                throw std::runtime_error("MSA profiles have different ungapped consensus lengths");
            break;
        }

        const char residue = normalized(reference_anchor[r]);
        if (residue != normalized(insertion_anchor[i]))
            throw std::runtime_error("MSA profiles disagree on consensus residue " +
                                     std::to_string(emitted_bases + 1));
        if (!expected_ungapped_anchor.empty() &&
            (emitted_bases >= expected_ungapped_anchor.size() ||
             residue != normalized(expected_ungapped_anchor[emitted_bases])))
            throw std::runtime_error("Profile anchor differs from RefAligner consensus");

        if (plan.columns == std::numeric_limits<std::size_t>::max())
            throw std::overflow_error("Merged MSA width overflows size_t");
        plan.reference_to_common.push_back(plan.columns);
        plan.insertion_to_common.push_back(plan.columns);
        ++plan.columns;
        ++r;
        ++i;
        ++emitted_bases;
    }
    if (emitted_bases == 0 ||
        (!expected_ungapped_anchor.empty() && emitted_bases != expected_ungapped_anchor.size()))
        throw std::runtime_error("Missing or truncated consensus in MSA profiles");
    if (plan.reference_to_common.size() != reference_anchor.size() ||
        plan.insertion_to_common.size() != insertion_anchor.size())
        throw std::logic_error("Profile projection map has invalid width");
    return plan;
}

std::string projectRow(std::string_view row,
                       const std::vector<std::size_t>& old_to_common,
                       std::size_t common_columns) {
    if (row.size() != old_to_common.size())
        throw std::invalid_argument("Profile row width does not match projection plan");
    std::string result(common_columns, '-');
    std::size_t previous = 0;
    for (std::size_t col = 0; col < row.size(); ++col) {
        const auto mapped = old_to_common[col];
        if (mapped >= common_columns || (col > 0 && mapped <= previous))
            throw std::invalid_argument("Invalid nonmonotone profile projection plan");
        result[mapped] = row[col];
        previous = mapped;
    }
    return result;
}

std::size_t reconcileFastaProfiles(
    const std::filesystem::path& reference_fasta,
    const std::filesystem::path& insertion_fasta,
    const std::filesystem::path& projected_reference_fasta,
    const std::filesystem::path& projected_insertion_fasta,
    std::string_view expected_ungapped_consensus) {
    if (reference_fasta == projected_reference_fasta ||
        insertion_fasta == projected_insertion_fasta ||
        reference_fasta == projected_insertion_fasta ||
        insertion_fasta == projected_reference_fasta ||
        projected_reference_fasta == projected_insertion_fasta)
        throw std::invalid_argument("Profile projection requires distinct input/output paths");

    seq_io::KseqReader reference_reader(reference_fasta);
    seq_io::KseqReader insertion_reader(insertion_fasta);
    seq_io::SeqRecord reference_first, insertion_first;
    if (!reference_reader.next(reference_first) || !insertion_reader.next(insertion_first))
        throw std::runtime_error("Cannot reconcile empty MSA profiles");
    // The reference FASTA can contain a user-named center. Its first row
    // need not be called 'consensus', but MUST match the ungapped consensus.
    // The insertion MSA always adds a synthetic 'consensus' as its first row.
    if (insertion_first.id != "consensus")
        throw std::runtime_error("Insertion profile must start with the 'consensus' record");

    const auto plan = makePlan(reference_first.seq, insertion_first.seq,
                               expected_ungapped_consensus);
    projectProfile(reference_fasta, projected_reference_fasta,
                   plan.reference_to_common, plan.columns);
    projectProfile(insertion_fasta, projected_insertion_fasta,
                   plan.insertion_to_common, plan.columns);
    return plan.columns;
}

}  // namespace align::profile_projection
