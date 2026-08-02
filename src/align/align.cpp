#include "align.h"
#include "ksw2.h"

extern "C" {
#include "alignment/cigar.h"
#include "wavefront/wavefront_align.h"
}

// Sequence alignment algorithm wrapper: KSW2 / WFA2
// - unified return cigar::Cigar_t (compressed format)
// - support global alignment, extension alignment, anchor-segment alignment

namespace align
{
    namespace {
        inline bool wfaCompactOffsetWouldOverflow(const std::string& ref, const std::string& query) {
#if HALIGN4_WFA_OFFSET_BITS == 16
            // Conservative bound: antidiagonal values in WFA kernels can approach ref_len + qry_len.
            // Keep a safety margin below int16 max to avoid undefined truncation behavior in compact mode.
            constexpr std::size_t kMaxSafeAntidiagonal = 30000;
            return (ref.size() + query.size()) > kMaxSafeAntidiagonal;
#else
            (void)ref;
            (void)query;
            return false;
#endif
        }
    } // namespace

    // KSW2 global alignment (end-to-end) - encode sequence and call KSW2
    cigar::Cigar_t globalAlignKSW2(const std::string& ref, const std::string& query)
    {
        align::KSW2AlignConfig cfg;
        return globalAlignKSW2(ref, query, cfg);
    }

    cigar::Cigar_t globalAlignKSW2(const std::string& ref, const std::string& query,
                                   align::KSW2AlignConfig cfg)
    {
        // Boundary: if any sequence empty, return pure I/D CIGAR
        if (ref.size() == 0 || query.size() == 0) {
            cigar::Cigar_t cigar;
            if (ref.size() == 0 && query.size() > 0) {
                cigar.push_back(cigar::cigarToInt('I', static_cast<uint32_t>(query.size())));
            } else if (query.size() == 0 && ref.size() > 0) {
                cigar.push_back(cigar::cigarToInt('D', static_cast<uint32_t>(ref.size())));
            }
            return cigar;
        }

        // Encode sequence: DNA5 (A/C/G/T/N -> 0..4)
        std::vector<uint8_t> ref_enc(ref.size());
        std::vector<uint8_t> qry_enc(query.size());

        for (size_t i = 0; i < ref.size(); ++i)
            ref_enc[i] = align::ScoreChar2Idx[static_cast<uint8_t>(ref[i])];
        for (size_t i = 0; i < query.size(); ++i)
            qry_enc[i] = align::ScoreChar2Idx[static_cast<uint8_t>(query[i])];

        cfg.band_width = align::auto_band(ref.size(), query.size());

        // Call KSW2
        ksw_extz_t ez{};
        ksw_extz2_sse(0,
            static_cast<int>(qry_enc.size()), qry_enc.data(),
            static_cast<int>(ref_enc.size()), ref_enc.data(),
            cfg.alphabet_size, cfg.mat,
            cfg.gap_open, cfg.gap_extend,
            cfg.band_width, cfg.zdrop, cfg.end_bonus,
            cfg.flag, &ez);

        // Copy and free CIGAR
        cigar::Cigar_t cigar;
        cigar.reserve(ez.n_cigar);
        for (int i = 0; i < ez.n_cigar; ++i)
            cigar.push_back(ez.cigar[i]);

        free(ez.cigar);
        return cigar;
    }

    // KSW2 extension alignment - use zdrop and EXTZ_ONLY flags
    cigar::Cigar_t extendAlignKSW2(const std::string& ref,
        const std::string& query,
        int zdrop)
    {
        // Encode sequence
        std::vector<uint8_t> ref_enc(ref.size());
        std::vector<uint8_t> qry_enc(query.size());
        for (size_t i = 0; i < ref.size(); ++i) ref_enc[i] = align::ScoreChar2Idx[(uint8_t)ref[i]];
        for (size_t i = 0; i < query.size(); ++i) qry_enc[i] = align::ScoreChar2Idx[(uint8_t)query[i]];

        // Configure params: EXTZ_ONLY + RIGHT + APPROX_DROP for extension
        align::KSW2AlignConfig cfg;
        cfg.mat = align::dna5_simd_mat;
        cfg.zdrop = zdrop;
        cfg.flag = KSW_EZ_EXTZ_ONLY | KSW_EZ_RIGHT | KSW_EZ_APPROX_DROP;
        cfg.end_bonus = 50;
        cfg.alphabet_size = 5;
        cfg.gap_open = 6;
        cfg.gap_extend = 2;
        cfg.band_width = align::auto_band(ref.size(), query.size());

        // Call KSW2
        ksw_extz_t ez{};
        ksw_extz2_sse(nullptr,
            static_cast<int>(qry_enc.size()), qry_enc.data(),
            static_cast<int>(ref_enc.size()), ref_enc.data(),
            cfg.alphabet_size, cfg.mat,
            cfg.gap_open, cfg.gap_extend,
            cfg.band_width, cfg.zdrop, cfg.end_bonus,
            cfg.flag, &ez);

        // Copy and free CIGAR
        cigar::Cigar_t cigar;
        cigar.reserve(ez.n_cigar);
        for (int i = 0; i < ez.n_cigar; ++i)
            cigar.push_back(ez.cigar[i]);

        free(ez.cigar);
        return cigar;
    }


    // WFA2 global alignment - use gap_affine mode
    cigar::Cigar_t globalAlignWFA2(const std::string& ref,
        const std::string& query)
    {
        if (wfaCompactOffsetWouldOverflow(ref, query)) {
            return globalAlignKSW2(ref, query);
        }

        // Build WFA2 attributes
        wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
        attributes.distance_metric = gap_affine;
        attributes.affine_penalties.mismatch = 3;
        attributes.affine_penalties.gap_opening = 4;
        attributes.affine_penalties.gap_extension = 1;
        attributes.memory_mode = wavefront_memory_high;

        // Create and run aligner
        wavefront_aligner_t* const wf_aligner = wavefront_aligner_new(&attributes);
        wavefront_align(wf_aligner, ref.c_str(), ref.length(), query.c_str(), query.length());

        // Extract CIGAR
        uint32_t* cigar_buffer = nullptr;
        int cigar_length = 0;
        cigar_get_CIGAR(wf_aligner->cigar, false, &cigar_buffer, &cigar_length);

        cigar::Cigar_t cigar;
        cigar.reserve(static_cast<std::size_t>(cigar_length));
        for (int i = 0; i < cigar_length; ++i)
            cigar.push_back(cigar_buffer[i]);

        wavefront_aligner_delete(wf_aligner);
        return cigar;
    }

    // ------------------------------------------------------------------
    // extendAlignWFA2: WFA2 extension alignment
    // ------------------------------------------------------------------
    // WFA2 extension alignment (ends-free extension):
    // - typically used to extend outward from seed position, quickly get local alignment
    // - configurable zdrop threshold, control early termination during extension
    //
    // This implementation uses high-memory mode of wavefront_aligner_attr_default:
    // - suits long sequence alignment, but may cause high memory usage
    // ------------------------------------------------------------------
    // cigar::Cigar_t extendAlignWFA2(const std::string& ref,
    //     const std::string& query, int zdrop)
    // {
    //     wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
    //     attributes.distance_metric = gap_affine;
    //     attributes.affine_penalties.mismatch = 2;      // X > 0
    //     attributes.affine_penalties.gap_opening = 3;   // O >= 0
    //     attributes.affine_penalties.gap_extension = 1; // E > 0
    //     attributes.memory_mode = wavefront_memory_high;
    //     attributes.heuristic.strategy = wf_heuristic_zdrop;
    //     attributes.heuristic.zdrop = zdrop;
    //     attributes.heuristic.steps_between_cutoffs = 1;
    //     //// Create a WFAligner
    //     //
    //     wavefront_aligner_t* const wf_aligner = wavefront_aligner_new(&attributes);
    //
    //     wavefront_align(wf_aligner, ref.c_str(), ref.length(), query.c_str(), query.length());
    //     /*wfa::WFAlignerGapAffine aligner(2, 3, 1, wfa::WFAligner::Alignment, wfa::WFAligner::MemoryUltralow);
    //
    //     aligner.alignEnd2End(ref, query);*/
    //
    //     uint32_t* cigar_buffer; // Buffer to hold the resulting CIGAR operations.
    //     int cigar_length = 0; // Length of the CIGAR string.
    //     // Retrieve the CIGAR string from the wavefront aligner.
    //     cigar_get_CIGAR(wf_aligner->cigar, true, &cigar_buffer, &cigar_length);
    //
    //     /* ---------- 4. Copy / Free CIGAR ---------- */
    //     cigar::Cigar_t cigar;
    //
    //     for (int i = 0; i < cigar_length; ++i)
    //         cigar.push_back(cigar_buffer[i]);
    //
    //     wavefront_aligner_delete(wf_aligner);
    //
    //     return cigar;
    // }

    // Anchor-based segment global alignment (minimap2 style)
    // - split by anchors into segments, align each segment globally, merge results
    cigar::Cigar_t globalAlignMM2(const std::string& ref,
                                  const std::string& query,
                                  const anchor::Anchors& anchors)
    {
        align::KSW2AlignConfig cfg;
        align::KSW2AlignConfig first_cfg;
        first_cfg.flag = KSW_EZ_GENERIC_SC;

        const std::size_t ref_len = ref.size();
        const std::size_t qry_len = query.size();

        // Chain anchors: get best chain
        anchor::Anchors sorted_anchors = anchors;
        anchor::ChainParams chain_params = anchor::default_chain_params();
        anchor::Anchors chain_anchors = anchor::chainAnchors(sorted_anchors, chain_params);
        if (chain_anchors.empty()) {
            return globalAlignKSW2(ref, query);
        }

        // Sort by query coordinate
        std::sort(chain_anchors.begin(), chain_anchors.end(),
                  [](const anchor::Anchor& a, const anchor::Anchor& b) {
                      if (a.pos_qry != b.pos_qry) return a.pos_qry < b.pos_qry;
                      return a.pos_ref < b.pos_ref;
                  });

        cigar::Cigar_t result;
        result.reserve(chain_anchors.size() * 2 + 2);

        std::size_t ref_pos = 0;
        std::size_t qry_pos = 0;

        auto append_segment = [&](std::size_t ref_start, std::size_t ref_end,
                                  std::size_t qry_start, std::size_t qry_end, align::KSW2AlignConfig seg_cfg) {
            // Boundary clipping
            ref_start = std::min(ref_start, ref_len);
            ref_end = std::min(ref_end, ref_len);
            qry_start = std::min(qry_start, qry_len);
            qry_end = std::min(qry_end, qry_len);

            if (ref_end < ref_start) ref_end = ref_start;
            if (qry_end < qry_start) qry_end = qry_start;

            const std::string seg_ref = ref.substr(ref_start, ref_end - ref_start);
            const std::string seg_qry = query.substr(qry_start, qry_end - qry_start);

            cigar::Cigar_t seg_cigar = globalAlignKSW2(seg_ref, seg_qry, seg_cfg);

            // Use CIGAR to infer consumed length
            const std::size_t seg_ref_len = seg_ref.size();
            const std::size_t seg_qry_len = seg_qry.size();
            const std::size_t c_ref = cigar::getRefLength(seg_cigar);
            const std::size_t c_qry = cigar::getQueryLength(seg_cigar);

            if (c_ref != seg_ref_len || c_qry != seg_qry_len) {
#ifdef _DEBUG
                spdlog::warn("globalAlignMM2(seg): segment cigar mismatch (expected ref:{}/qry:{}, got ref:{}/qry:{}); forcing robust fallback for this segment",
                             seg_ref_len, seg_qry_len, c_ref, c_qry);
#endif
                // Fallback strategy: Query all I, Ref all D
                cigar::Cigar_t forced_cigar;
                if (seg_qry_len > 0) {
                    forced_cigar.push_back(cigar::cigarToInt('I', static_cast<uint32_t>(seg_qry_len)));
                }
                if (seg_ref_len > 0) {
                    forced_cigar.push_back(cigar::cigarToInt('D', static_cast<uint32_t>(seg_ref_len)));
                }
                cigar::appendCigar(result, forced_cigar);

                ref_pos = ref_end;
                qry_pos = qry_end;
                return;
            }

            cigar::appendCigar(result, seg_cigar);

            ref_pos = ref_start + c_ref;
            qry_pos = qry_start + c_qry;
        };

        // Left end: start to first anchor
        {
            const auto& first = chain_anchors.front();
            append_segment(ref_pos, first.pos_ref, qry_pos, first.pos_qry, first_cfg);
        }

        // Per anchor: handle span and gap
        for (std::size_t i = 0; i < chain_anchors.size(); ++i) {
            const auto& a = chain_anchors[i];

            const std::size_t a_ref_start = static_cast<std::size_t>(a.pos_ref);
            const std::size_t a_qry_start = static_cast<std::size_t>(a.pos_qry);
            const std::size_t a_ref_end = a_ref_start + static_cast<std::size_t>(a.span);
            const std::size_t a_qry_end = a_qry_start + static_cast<std::size_t>(a.span);

            append_segment(ref_pos, a_ref_end, qry_pos, a_qry_end, cfg);

            if (i + 1 < chain_anchors.size()) {
                const auto& b = chain_anchors[i + 1];
                append_segment(ref_pos, b.pos_ref, qry_pos, b.pos_qry, cfg);
            }
        }

        // Right end: last anchor to end
        append_segment(ref_pos, ref_len, qry_pos, qry_len, cfg);

        // Final consistency check
        const std::size_t total_ref = cigar::getRefLength(result);
        const std::size_t total_qry = cigar::getQueryLength(result);
        if (total_ref != ref_len || total_qry != qry_len) {
            spdlog::error("globalAlignMM2: final cigar mismatch (ref:{}/{}, qry:{}/{}), fallback to global",
                         total_ref, ref_len, total_qry, qry_len);
            return globalAlignKSW2(ref, query);
        }

        return result;
    }

} // namespace align
