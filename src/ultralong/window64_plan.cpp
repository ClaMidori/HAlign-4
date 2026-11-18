#include "ultralong/window64_plan.hpp"
#include "ultralong/sparse_fasta.hpp"
#include "align/homology_plan.hpp"
#include "ultralong/rescue_geometry.hpp"
#include <algorithm>
#include <cstdlib>
#include <unistd.h>
#include <array>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace fs=std::filesystem;
namespace ultralong {
namespace {
struct Track {
    std::int64_t shift = 0; // q_global - r_global; updated ONLY on supported proposals
    std::uint64_t last_q_end = 0;
    bool initialized = false;
};
std::uint64_t clampedStart(std::int64_t predicted,std::uint64_t flank,std::uint64_t length) {
    const auto start=predicted <= 0 ? std::uint64_t{0} : static_cast<std::uint64_t>(predicted);
    return start < flank ? 0 : std::min(start-flank,length);
}
std::string statusFor(const homology_plan::Plan& p,
                      std::size_t reference_window,
                      std::uint64_t fetch_start, std::uint64_t full_query_length,
                      Track& track, std::uint64_t reference_start,
                      std::uint64_t reference_end,
                      std::uint64_t& qa, std::uint64_t& qb,
                      std::uint64_t& anchors_count) {
    const auto& anchors=p.supported_anchors.at(0);
    anchors_count=anchors.size();
    if (anchors.size()<4) return "UNRESOLVED_NO_SUPPORT";
    const auto& left=p.boundaries.front().projections.at(0);
    const auto& right=p.boundaries.back().projections.at(0);
    if (!left || !right) return "UNRESOLVED_EDGE";
    // Require independent evidence across every quarter. This rules out a
    // cluster at the endpoints masquerading as support throughout the window.
    std::array<bool,4> quarters{};
    for (const auto& a:anchors) {
        const auto quarter=std::min<std::size_t>(3, static_cast<std::uint64_t>(a.reference)*4/reference_window);
        quarters[quarter]=true;
    }
    if (!std::all_of(quarters.begin(),quarters.end(),[](bool x){return x;}))
        return "UNRESOLVED_PARTIAL_COVERAGE";
    if (left->estimated_query >= right->estimated_query) return "UNRESOLVED_NONMONOTONE";
    if (fetch_start>std::numeric_limits<std::uint64_t>::max()-right->estimated_query)
        throw std::overflow_error("query projected coordinate overflow");
    qa=fetch_start+left->estimated_query;
    qb=fetch_start+right->estimated_query;
    if (qb>full_query_length || qa>=qb) return "UNRESOLVED_RANGE";
    // Small per-window length skew only; do not claim whole long indels.
    const auto r_length=reference_end-reference_start;
    const auto q_length=qb-qa;
    const auto drift = r_length>q_length ? r_length-q_length : q_length-r_length;
    if (drift>std::max<std::uint64_t>(1024,r_length/8)) return "UNRESOLVED_LENGTH_SKEW";
    // The global query lane is monotone. Once a region has been proposed,
    // never silently place a subsequent window on an overlapping query region.
    if (track.initialized && qa<track.last_q_end) return "UNRESOLVED_OVERLAP";
    if (qa>static_cast<std::uint64_t>(INT64_MAX) || qb>static_cast<std::uint64_t>(INT64_MAX) ||
        reference_start>static_cast<std::uint64_t>(INT64_MAX) ||
        reference_end>static_cast<std::uint64_t>(INT64_MAX))
        throw std::overflow_error("signed drift conversion would overflow");
    track.shift=static_cast<std::int64_t>(qb)-static_cast<std::int64_t>(reference_end);
    track.last_q_end=qb;
    track.initialized=true;
    return "CANDIDATE_ONLY";
}
} // anonymous
Window64Stats planWindows(const fs::path& input,const fs::path& output,Window64Config cfg) {
    if (cfg.window_bases<8192 || cfg.window_bases>262144 ||
        cfg.flank_bases<1024 || cfg.flank_bases>32768 ||
        cfg.flank_bases*2>=cfg.window_bases ||
        cfg.max_reads<2 || cfg.max_reads>8 ||
        cfg.max_total_bases<2 || cfg.max_total_bases>800000000ULL)
        throw std::invalid_argument("invalid bounded window planner budget");
    const auto source=fs::absolute(input).lexically_normal();
    const auto dest=fs::absolute(output).lexically_normal();
    if(source==dest || fs::exists(dest) || fs::is_symlink(dest) || !fs::is_directory(dest.parent_path()))
        throw std::invalid_argument("planner output must be new and distinct from input");
    auto fasta=SparseFasta::build(source,cfg.checkpoint_stride,cfg.max_reads);
    if(fasta.reads().size()<2||fasta.reads().size()>cfg.max_reads ||
       fasta.totalBases()>cfg.max_total_bases)
        throw std::invalid_argument("input exceeds total/read count budget");
    const auto ref_size=fasta.reads().front().bases;
    if(ref_size>static_cast<std::uint64_t>(INT64_MAX/2))
        throw std::invalid_argument("reference too large for signed coordinate drift");
    Window64Stats stats;
    stats.total_bases=fasta.totalBases();
    const auto total_window=cfg.window_bases+2*cfg.flank_bases;
    if(total_window>1000000)throw std::invalid_argument("window too large for legacy planner");
    homology_plan::Limits local_limits;
    local_limits.max_sequences=2;
    local_limits.max_per_sequence=static_cast<std::size_t>(total_window);
    local_limits.max_total_bases=static_cast<std::size_t>(2*total_window);
    homology_plan::Config pc;
    pc.target_block_bases=5000;
    pc.boundary_radius=1000;
    pc.strict_candidate_guard=true;
    pc.recover_indel_flanks=true;
    std::vector<Track> tracks(fasta.reads().size());
    // Two-stage, no-overwrite publication: hold small TSV in a temporary file.
    // If runtime fails halfway, output does not appear.
    auto temp=dest.parent_path()/("."+dest.filename().string()+".stage.XXXXXX");
    auto pattern=temp.string();
    std::vector<char> buf(pattern.begin(),pattern.end());buf.push_back('\0');
    int fd=mkstemp(buf.data());
    if(fd<0)throw std::runtime_error("cannot create temporary report");
    close(fd);
    const fs::path stage(buf.data());
    struct Cleanup {fs::path p;~Cleanup(){std::error_code ec;fs::remove(p,ec);}} cleanup{stage};
    std::ofstream journal(stage,std::ios::binary|std::ios::trunc);
    if(!journal)throw std::runtime_error("cannot write staged plan");
    journal<<"# PR3D.4D.2 candidate-only, NO_MSA, NO_ORACLE, NO_VERIFIED_HOMOLOGY\n";
    journal<<"window\tread_id\tref_start\tref_end\tquery_start_est\tquery_end_est\tanchors\tstatus\n";
    for(std::uint64_t r0=0;r0<ref_size;){
        const auto rlen=std::min(cfg.window_bases,ref_size-r0);
        const auto r1=r0+rlen;
        // No global sequence is read into RAM. One bounded reference window
        // and one bounded query interval at a time, for every pair.
        const auto ref=fasta.slice(0,r0,rlen,total_window);
        ++stats.windows;
        for(std::size_t r=1;r<fasta.reads().size();++r){
            const auto& item=fasta.reads()[r];
            const auto& track=tracks[r];
            const auto predicted=static_cast<std::int64_t>(r0)+track.shift;
            const auto qs=clampedStart(predicted,cfg.flank_bases,item.bases);
            const auto qlen=std::min<std::uint64_t>(item.bases-qs,rlen+cfg.flank_bases*2);
            std::string status="UNRESOLVED_OUT_OF_RANGE";
            std::uint64_t qa=0,qb=0,anchor_count=0;
            if(qlen>=std::min<std::uint64_t>(rlen,4096) && rlen>=4096) {
                const auto query=fasta.slice(r,qs,qlen,total_window);
                std::vector<homology_plan::Sequence> local{{fasta.reads()[0].id,ref},{item.id,query}};
                // These *existing* PR3D.2C routines are strictly local and
                // continue to use 32-bit coordinates *within* a window only.
                const auto plan=homology_plan::build(local,pc,local_limits);
                status=statusFor(plan,rlen,qs,item.bases,tracks[r],r0,r1,qa,qb,anchor_count);
                // PR3D.4D.5B: opt-in bounded second pass ONLY for unresolved
                // windows. Never touches a primary candidate or its shift.
                if (cfg.rescue_unresolved && status!="CANDIDATE_ONLY") {
                    for (unsigned k : {10U,9U}) {
                        auto alternate=pc;
                        alternate.k=k;
                        alternate.allow_short_k=true;
                        const auto rescue_plan=homology_plan::build(local,alternate,local_limits);
                        if(rescue_plan.supported_anchors.empty())continue;
                        const auto proposed=projectRescueRange(
                            rescue_plan.supported_anchors[0],ref,query,0,
                            static_cast<std::uint32_t>(rlen));
                        if(!proposed)continue;
                        Track trial=tracks[r];
                        std::uint64_t a=0,b=0,c=0;
                        const auto trial_status=statusFor(rescue_plan,rlen,qs,item.bases,
                            trial,r0,r1,a,b,c);
                        if(trial_status!="CANDIDATE_ONLY")continue;
                        if(std::llabs(static_cast<std::int64_t>(a-qs)-proposed->query_start)>32 ||
                           std::llabs(static_cast<std::int64_t>(b-qs)-proposed->query_end)>32)continue;
                        tracks[r]=trial;
                        qa=a;qb=b;anchor_count=c;
                        status="RESCUE_CANDIDATE_ONLY";
                        ++stats.rescued_candidate_pairs;
                        break;
                    }
                }
            }
            journal<<stats.windows<<'\t'<<item.id<<'\t'<<r0<<'\t'<<r1<<'\t';
            if(status=="CANDIDATE_ONLY" || status=="RESCUE_CANDIDATE_ONLY") {
                journal<<qa<<'\t'<<qb<<'\t';
                ++stats.candidate_pairs;
                stats.anchored_ref_bases+=rlen;
            } else {
                journal<<"NA\tNA\t";
                ++stats.unresolved_pairs;
            }
            journal<<anchor_count<<'\t'<<status<<'\n';
        }
        r0=r1;
    }
    journal.flush();journal.close();
    if(!journal)throw std::runtime_error("staged window report write failed");
    // create_hard_link provides no-replace publication, even against races.
    fs::create_hard_link(stage,dest);
    return stats;
}
} // namespace ultralong
