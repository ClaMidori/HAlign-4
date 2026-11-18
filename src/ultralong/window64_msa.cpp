#include "ultralong/window64_msa.hpp"
#include "ultralong/window64_plan.hpp"
#include "ultralong/sparse_fasta.hpp"
#include "align/homology_plan.hpp"
#include "align/native_msa.hpp"
#include "ultralong/rescue_geometry.hpp"

#include <algorithm>
#include <array>
#include <cerrno>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <limits>
#include <memory>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <unistd.h>

namespace fs=std::filesystem;
namespace ultralong {
namespace {
using Pos=std::uint64_t;
struct Projection { Pos from=0,to=0; bool candidate=false, rescued=false; };
struct Window { Pos from=0,to=0; std::vector<Projection> query; };
struct Block { Pos from=0,to=0; std::vector<Projection> query; };
std::vector<std::string> tabs(const std::string& s) {
    std::vector<std::string> f;std::size_t start=0;
    while(true){const auto end=s.find('\t',start);f.push_back(s.substr(start,end==std::string::npos?end:end-start));if(end==std::string::npos)break;start=end+1;}
    return f;
}
Pos number(const std::string& s) {
    if(s.empty() || s.size()>20 || (s.size()>1 && s[0]=='0')) throw std::runtime_error("invalid noncanonical TSV number");
    Pos n=0;
    for(char c:s){if(c<'0'||c>'9'||n>(UINT64_MAX-static_cast<unsigned>(c-'0'))/10)throw std::runtime_error("invalid/overflow TSV number");n=n*10+static_cast<unsigned>(c-'0');}
    return n;
}
class PlanReader {
public:
    PlanReader(const fs::path& path,const SparseFasta& fasta): in_(path,std::ios::binary),fasta_(fasta) {
        if(!in_)throw std::runtime_error("cannot read staged window plan");
        std::string line;
        if(!std::getline(in_,line) || line!="# PR3D.4D.2 candidate-only, NO_MSA, NO_ORACLE, NO_VERIFIED_HOMOLOGY" ||
           !std::getline(in_,line) || line!="window\tread_id\tref_start\tref_end\tquery_start_est\tquery_end_est\tanchors\tstatus")
            throw std::runtime_error("unexpected window plan schema");
    }
    bool next(Window& w) {
        if(cursor_==fasta_.reads()[0].bases){
            std::string line;if(std::getline(in_,line))throw std::runtime_error("extra window plan records");
            if(in_.bad())throw std::runtime_error("window plan read error");
            return false;
        }
        if(!std::getline(in_,line_))throw std::runtime_error("truncated window plan");
        const auto first=tabs(line_);
        if(first.size()!=8 || first[1]!=fasta_.reads()[1].id)throw std::runtime_error("bad plan row order");
        const Pos win_idx=++window_index_;
        if(number(first[0])!=win_idx)throw std::runtime_error("window index discontinuity");
        const Pos from=number(first[2]),to=number(first[3]);
        if(from!=cursor_ || from>=to || to>fasta_.reads()[0].bases || to-from>65536)
            throw std::runtime_error("invalid window boundaries");
        w=Window{from,to,std::vector<Projection>(fasta_.reads().size()-1)};
        parseRow(first,1,w);
        for(std::size_t r=2;r<fasta_.reads().size();++r){
            if(!std::getline(in_,line_))throw std::runtime_error("window plan missing query row");
            const auto f=tabs(line_);
            if(f.size()!=8 || f[1]!=fasta_.reads()[r].id || number(f[0])!=win_idx)
                throw std::runtime_error("window plan read id/order mismatch");
            parseRow(f,r,w);
        }
        cursor_=to;
        return true;
    }
private:
    void parseRow(const std::vector<std::string>& f,std::size_t r,Window& w) {
        if(number(f[2])!=w.from || number(f[3])!=w.to)throw std::runtime_error("window plan interval mismatch");
        auto& q=w.query.at(r-1);
        if(f[7]=="CANDIDATE_ONLY" || f[7]=="RESCUE_CANDIDATE_ONLY"){
            q.from=number(f[4]);q.to=number(f[5]);
            if(q.from>=q.to||q.to>fasta_.reads()[r].bases || number(f[6])<4)
                throw std::runtime_error("unsafe candidate geometry in window plan");
            q.candidate=true; q.rescued=(f[7]=="RESCUE_CANDIDATE_ONLY");
        }else if(f[7].rfind("UNRESOLVED_",0)==0){
            if(f[4]!="NA"||f[5]!="NA")throw std::runtime_error("unresolved plan must have NA bounds");
        }else throw std::runtime_error("unknown window planner status");
    }
    std::ifstream in_;const SparseFasta& fasta_;
    std::string line_;Pos cursor_=0,window_index_=0;
};

// Validate every read/range returned by the pre-existing local planner.
std::vector<Block> localBlocks(const homology_plan::Plan& p,const std::vector<homology_plan::Sequence>& reads) {
    if(p.boundaries.size()<2)throw std::runtime_error("local plan has no blocks");
    std::vector<Block> out(p.boundaries.size()-1);
    std::vector<std::vector<bool>> seen(out.size(),std::vector<bool>(reads.size()));
    std::istringstream input(homology_plan::asTsv(p));
    std::string line;
    bool header=false;
    while(std::getline(input,line)){
        if(line.empty()||line[0]=='#')continue;
        if(!header){if(line!="block\tread_id\tref_start\tref_end\tquery_start_est\tquery_end_est\tleft_anchor_distance\tright_anchor_distance\tsupported_seed_count\tstatus")
            throw std::runtime_error("unexpected local planner schema");header=true;continue;}
        const auto f=tabs(line);
        if(f.size()!=10)throw std::runtime_error("bad local planner row");
        const Pos b=number(f[0]);
        if(b<1||b>out.size())throw std::runtime_error("local planner block out of range");
        std::size_t r=reads.size();
        for(std::size_t j=0;j<reads.size();++j)if(f[1]==reads[j].id){r=j;break;}
        if(r==reads.size() || seen[b-1][r])throw std::runtime_error("local planner duplicate/unknown read");
        seen[b-1][r]=true;
        auto& block=out[b-1];
        if(block.query.empty())block.query.resize(reads.size()-1);
        const Pos start=number(f[2]),end=number(f[3]);
        if(start!=p.boundaries[b-1].reference || end!=p.boundaries[b].reference || end<=start)
            throw std::runtime_error("inconsistent local planner reference bounds");
        block.from=start;block.to=end;
        if(r==0){
            if(f[9]!="REFERENCE_ONLY" || number(f[4])!=start || number(f[5])!=end)
                throw std::runtime_error("bad local reference row");
        }else if(f[9]=="CANDIDATE_ONLY"){
            auto& q=block.query[r-1];
            q.from=number(f[4]);q.to=number(f[5]);
            if(q.from>=q.to || q.to>reads[r].bases.size() || number(f[8])<3)
                throw std::runtime_error("out-of-bounds local projection");
            q.candidate=true;
        }else if(f[9]!="UNRESOLVED")throw std::runtime_error("unknown local status");
    }
    if(!header)throw std::runtime_error("local planner header missing");
    for(const auto& row:seen)for(bool found:row)if(!found)throw std::runtime_error("truncated local planner TSV");
    return out;
}

class Stage {
public:
    Stage(const fs::path& output,const fs::path& report,const SparseFasta& fasta):out_(output),report_(report),fasta_(fasta),cursor_(fasta.reads().size()) {
        if(out_==report_ || out_==fasta.path()||report_==fasta.path()||out_.parent_path()!=report_.parent_path()||
           fs::exists(out_)||fs::exists(report_)||!fs::is_directory(out_.parent_path()))
            throw std::invalid_argument("input, output, report must be distinct new files in a shared directory");
        auto temp=out_.parent_path()/(".halign4_outofcore_XXXXXX");
        auto pattern=temp.string();std::vector<char> buffer(pattern.begin(),pattern.end());buffer.push_back('\0');
        if(!mkdtemp(buffer.data()))throw std::runtime_error("cannot allocate private MSA staging directory");
        dir_=buffer.data();
        for(std::size_t r=0;r<fasta_.reads().size();++r){
            auto file=std::make_unique<std::ofstream>(dir_/("read_"+std::to_string(r)+".raw"),std::ios::binary|std::ios::trunc);
            if(!*file)throw std::runtime_error("cannot open MSA row spool");
            spools_.push_back(std::move(file));
        }
        journal_.open(dir_/"blocks.tsv",std::ios::binary|std::ios::trunc);
        if(!journal_)throw std::runtime_error("cannot open report spool");
        journal_<<"status\twindow\tblock\tread_id\tquery_start\tquery_end\tref_start\tref_end\tcolumns\tnote\n";
    }
    ~Stage(){if(!dir_.empty()){std::error_code ec;fs::remove_all(dir_,ec);}}
    const fs::path& dir()const{return dir_;}
    Pos cursor(std::size_t r)const{return cursor_.at(r);}
    Pos columns()const{return width_;}
    Pos privateBases()const{return private_bases_;}
    void privateRegion(std::size_t r,Pos to,Pos window,Pos block,const char* reason) {
        if(r>=cursor_.size()||to<cursor_[r]||to>fasta_.reads()[r].bases)throw std::runtime_error("private region out of bounds");
        const Pos begin=cursor_[r];
        static constexpr Pos CHUNK=1<<16;
        Pos pos=begin;
        while(pos<to){
            const auto n=std::min(CHUNK,to-pos);
            auto bases=fasta_.slice(r,pos,n,CHUNK);
            for(std::size_t j=0;j<cursor_.size();++j){
                if(j==r)write(j,bases);else gaps(j,n);
            }
            width_+=n;private_bases_+=n;pos+=n;
        }
        if(to!=begin) journal_<<"PRIVATE\t"<<window<<'\t'<<block<<'\t'<<fasta_.reads()[r].id<<'\t'<<begin<<'\t'<<to
                              <<"\tNA\tNA\t"<<to-begin<<'\t'<<reason<<'\n';
        cursor_[r]=to;
    }
    void accepted(Pos window,Pos block,const std::vector<std::size_t>& chosen,
                  const std::vector<Pos>& begin,const std::vector<Pos>& end,
                  const native_msa::Alignment& msa,Pos ref_from,Pos ref_to,
                  const std::vector<bool>& rescued){
        if(chosen.empty()||chosen.size()!=msa.rows.size()||chosen.size()!=rescued.size()||msa.columns==0)throw std::runtime_error("invalid local MSA shape");
        for(std::size_t i=0;i<chosen.size();++i){
            const auto r=chosen[i];
            if(r>=cursor_.size()||begin[i]!=cursor_[r]||end[i]<begin[i]||end[i]>fasta_.reads()[r].bases||
               msa.rows[i].id!=fasta_.reads()[r].id || msa.rows[i].bases.size()!=msa.columns)
                throw std::runtime_error("local MSA row ordering/cursor/width violation");
            std::string observed;
            for(char c:msa.rows[i].bases){
                if(c!='-'){
                    if(c!='A'&&c!='C'&&c!='G'&&c!='T'&&c!='N')throw std::runtime_error("local MSA returned unexpected base");
                    observed.push_back(c);
                }
            }
            if(observed!=fasta_.slice(r,begin[i],end[i]-begin[i],100000))
                throw std::runtime_error("local MSA DNA reconstruction failed");
        }
        for(std::size_t r=0;r<cursor_.size();++r){
            const auto it=std::find(chosen.begin(),chosen.end(),r);
            if(it==chosen.end())gaps(r,msa.columns);
            else write(r,msa.rows[static_cast<std::size_t>(it-chosen.begin())].bases);
        }
        for(std::size_t i=0;i<chosen.size();++i){
            const auto r=chosen[i];cursor_[r]=end[i];
            journal_<<"ANCHORED_LOCAL_ABPOA\t"<<window<<'\t'<<block<<'\t'<<fasta_.reads()[r].id<<'\t'
                    <<begin[i]<<'\t'<<end[i]<<'\t'<<ref_from<<'\t'<<ref_to<<'\t'<<msa.columns
                    <<(rescued[i] ? "\tRESCUE_K9_10_HEURISTIC_NOT_VERIFIED\n" :
                                      "\tHEURISTIC_NOT_VERIFIED\n");
        }
        width_+=msa.columns;
    }
    void publish() {
        for(std::size_t r=0;r<cursor_.size();++r)
            if(cursor_[r]!=fasta_.reads()[r].bases)throw std::runtime_error("not all original DNA was consumed");
        for(auto& o:spools_){o->flush();o->close();if(!*o)throw std::runtime_error("spool finalization failed");}
        journal_.flush();journal_.close();if(!journal_)throw std::runtime_error("report finalization failed");
        std::ofstream out(dir_/"final.fasta",std::ios::binary|std::ios::trunc);
        if(!out)throw std::runtime_error("cannot write staged aligned FASTA");
        static constexpr std::size_t SZ=1<<16;
        std::array<char,SZ> buf{};
        for(std::size_t r=0;r<cursor_.size();++r){
            const auto path=dir_/("read_"+std::to_string(r)+".raw");
            if(fs::file_size(path)!=width_)throw std::runtime_error("MSA row spool width mismatch");
            std::ifstream in(path,std::ios::binary);
            if(!in)throw std::runtime_error("cannot reopen MSA spool");
            out<<'>'<<fasta_.reads()[r].id<<'\n';
            Pos cursor=0, ref_start=0, ref_end=0;
            std::string ref_buffer;
            Pos col=0;
            while(in){
                in.read(buf.data(),static_cast<std::streamsize>(buf.size()));const auto count=in.gcount();
                for(std::streamsize i=0;i<count;++i){
                    const auto c=buf[static_cast<std::size_t>(i)];
                    if(c!='-'){
                        if(cursor>=ref_end){
                            if(cursor>=fasta_.reads()[r].bases)throw std::runtime_error("MSA has extra bases");
                            ref_start=cursor;ref_end=cursor+std::min<Pos>(SZ,fasta_.reads()[r].bases-cursor);
                            ref_buffer=fasta_.slice(r,ref_start,ref_end-ref_start,SZ);
                        }
                        if(c!=ref_buffer[static_cast<std::size_t>(cursor-ref_start)])
                            throw std::runtime_error("final streamed FASTA mismatches input DNA");
                        ++cursor;
                    }
                    out.put(c);++col;if(col%80==0)out.put('\n');
                }
            }
            if(in.bad()||col!=width_||cursor!=fasta_.reads()[r].bases)
                throw std::runtime_error("final FASTA not rectangular/lossless");
            if(col%80!=0)out.put('\n');
        }
        out.flush();out.close();if(!out)throw std::runtime_error("final FASTA write failed");
        fs::create_hard_link(dir_/"final.fasta",out_);
        try{fs::create_hard_link(dir_/"blocks.tsv",report_);}catch(...){std::error_code ec;fs::remove(out_,ec);throw;}
    }
    std::ofstream& report(){return journal_;}
private:
    void write(std::size_t r,const std::string& text) {
        spools_[r]->write(text.data(),static_cast<std::streamsize>(text.size()));
        if(!*spools_[r])throw std::runtime_error("MSA spool write error");
    }
    void gaps(std::size_t r,Pos n) {
        static const std::string chars(1<<16,'-');
        while(n){
            const auto take=static_cast<std::size_t>(std::min<Pos>(n,chars.size()));
            spools_[r]->write(chars.data(),static_cast<std::streamsize>(take));
            if(!*spools_[r])throw std::runtime_error("gap spool write failed");
            n-=take;
        }
    }
    fs::path dir_,out_,report_;
    const SparseFasta& fasta_;
    std::vector<Pos> cursor_;
    std::vector<std::unique_ptr<std::ofstream>> spools_;
    std::ofstream journal_;
    Pos width_=0,private_bases_=0;
};

bool plausible(const native_msa::Alignment& a,const std::vector<native_msa::Sequence>& local,double minimum){
    if(a.rows.size()!=local.size()||a.columns==0||a.columns>100000)return false;
    for(std::size_t r=0;r<local.size();++r){
        if(a.rows[r].id!=local[r].id||a.rows[r].bases.size()!=a.columns)return false;
    }
    for(std::size_t r=1;r<local.size();++r){
        std::size_t paired=0,match=0;
        for(std::size_t c=0;c<a.columns;++c){
            const auto ref=a.rows[0].bases[c],q=a.rows[r].bases[c];
            if(ref!='-'&&q!='-'){
                ++paired;if(ref==q&&ref!='N')++match;
            }
        }
        const auto minlen=std::min(local[0].bases.size(),local[r].bases.size());
        if(paired*100<minlen*65 || paired==0||static_cast<double>(match)/paired<minimum)return false;
    }
    return true;
}
} // anonymous namespace

WindowMsaStats runWindowMsa(const fs::path& in,const fs::path& output,const fs::path& report,WindowMsaConfig cfg) {
    if(cfg.max_reads<2 || cfg.max_reads>8 || cfg.max_total_bases<2 || cfg.max_total_bases>800000000ULL ||
       cfg.local_max_bases<1000 || cfg.local_max_bases>12000 ||
       cfg.local_max_total_bases<cfg.local_max_bases ||cfg.local_max_total_bases>100000 ||
       !(cfg.min_reference_identity>=0.5 && cfg.min_reference_identity<=1))
        throw std::invalid_argument("invalid 3D4D3 bounded MSA config");
    const auto source=fs::absolute(in).lexically_normal();
    const auto out=fs::absolute(output).lexically_normal();
    const auto side=fs::absolute(report).lexically_normal();
    if(source==out||source==side||out==side||!fs::is_directory(out.parent_path())||
       out.parent_path()!=side.parent_path()||fs::exists(out)||fs::exists(side))
        throw std::invalid_argument("MSA output/report must be distinct new paths in shared parent");
    const auto fasta=SparseFasta::build(source,65536,cfg.max_reads);
    if(fasta.reads().size()<2||fasta.totalBases()>cfg.max_total_bases)
        throw std::invalid_argument("source exceeds DNA/reads bounded budget");
    Stage stage(out,side,fasta);
    Window64Config wc;
    wc.max_reads=cfg.max_reads;wc.max_total_bases=cfg.max_total_bases;
    wc.rescue_unresolved=cfg.rescue_unresolved;
    const auto pstats=planWindows(source,stage.dir()/"window_plan.tsv",wc);
    PlanReader plans(stage.dir()/"window_plan.tsv",fasta);
    WindowMsaStats stats;
    stats.original_bases=fasta.totalBases();stats.windows=pstats.windows;
    stats.candidate_query_windows=pstats.candidate_pairs;
    stats.rescued_window_candidates=pstats.rescued_candidate_pairs;
    homology_plan::Config pc;
    pc.target_block_bases=5000;pc.boundary_radius=1000;
    pc.strict_candidate_guard=true;pc.recover_indel_flanks=true;
    Window w;Pos win_idx=0;
    while(plans.next(w)){
        ++win_idx;
        // Each window candidate interval is sliced ONCE and expressed locally.
        // The local homology planner introduces smaller evidence-based partitions.
        std::vector<homology_plan::Sequence> local;
        std::vector<std::size_t> index;
        std::vector<Pos> global_begin;
        local.push_back({fasta.reads()[0].id,fasta.slice(0,w.from,w.to-w.from,65536)});
        index.push_back(0);global_begin.push_back(w.from);
        for(std::size_t r=1;r<fasta.reads().size();++r){
            const auto q=w.query[r-1];
            if(!q.candidate || q.from<stage.cursor(r))continue;
            if(q.to-q.from>100000)continue;
            local.push_back({fasta.reads()[r].id,fasta.slice(r,q.from,q.to-q.from,100000)});
            index.push_back(r);global_begin.push_back(q.from);
        }
        if(local.size()<2){stage.privateRegion(0,w.to,win_idx,0,"WINDOW_NOT_SUPPORTED");continue;}
        homology_plan::Limits limits{local.size(),100000,static_cast<std::size_t>(local.size()*100000)};
        const auto plan=homology_plan::build(local,pc,limits);
        auto blocks=localBlocks(plan,local);
        // Primary reference partitions are preserved whenever they contain
        // at least one aligned candidate. Recovery never replaces accepted
        // primary sub-block boundaries in that case.
        std::vector<std::vector<homology_plan::Anchor>> alternate_lanes(local.size()-1);
        bool needs_rescue=false;
        if(cfg.rescue_unresolved){
            for(const auto& b:blocks){
                for(const auto& q:b.query)if(!q.candidate)needs_rescue=true;
            }
            if(needs_rescue){
                auto alt_cfg=pc;alt_cfg.k=9;alt_cfg.allow_short_k=true;
                const auto alt=homology_plan::build(local,alt_cfg,limits);
                alternate_lanes=alt.supported_anchors;
                // Only if ALL primary blocks are unsupported may we use
                // alternate partition boundaries. Otherwise preserve them.
                bool primary_has_candidate=false;
                for(const auto& b:blocks)
                    for(const auto& q:b.query)if(q.candidate)primary_has_candidate=true;
                if(!primary_has_candidate){
                    auto alternative_blocks=localBlocks(alt,local);
                    if(alternative_blocks.size()>blocks.size())
                        blocks=std::move(alternative_blocks);
                }
            }
        }
        for(std::size_t b=0;b<blocks.size();++b){
            const auto& block=blocks[b];
            const auto ref0=w.from+block.from,ref1=w.from+block.to;
            if(ref0!=stage.cursor(0) || ref1>w.to)throw std::runtime_error("local/global reference block boundary mismatch");
            std::vector<std::size_t> chosen{0};
            std::vector<Pos> begins{ref0}, ends{ref1};
            std::vector<bool> rescued{false};
            std::vector<native_msa::Sequence> dna{{local[0].id,local[0].bases.substr(block.from,block.to-block.from)}};
            std::size_t total=dna[0].bases.size();
            if(total<=cfg.local_max_bases){
                for(std::size_t i=1;i<local.size();++i){
                    auto q=block.query[i-1];
                    bool is_rescued=false;
                    if(!q.candidate && cfg.rescue_unresolved){
                        const auto proposal=projectRescueRange(
                            alternate_lanes[i-1],local[0].bases,local[i].bases,
                            static_cast<std::uint32_t>(block.from),
                            static_cast<std::uint32_t>(block.to));
                        if(proposal){
                            q.from=proposal->query_start;
                            q.to=proposal->query_end;
                            q.candidate=true;
                            is_rescued=true;
                        }
                    }
                    if(!q.candidate)continue;
                    ++stats.local_candidate_rows;
                    const auto r=index[i];
                    const auto start=global_begin[i]+q.from,stop=global_begin[i]+q.to;
                    if(start<stage.cursor(r)||stop<=start)continue;
                    if(stop-start>cfg.local_max_bases || stop-start>cfg.local_max_total_bases-total){
                        ++stats.skipped_budget;continue;
                    }
                    chosen.push_back(r);begins.push_back(start);ends.push_back(stop);
                    rescued.push_back(is_rescued || w.query[r-1].rescued);
                    dna.push_back({local[i].id,local[i].bases.substr(q.from,q.to-q.from)});
                    total+=q.to-q.from;
                }
            }else ++stats.skipped_budget;
            if(chosen.size()<2){stage.privateRegion(0,ref1,win_idx,b+1,"LOCAL_NO_SAFE_CANDIDATE");continue;}
            native_msa::Options options;
            options.max_sequences=cfg.max_reads;
            options.max_sequence_bases=cfg.local_max_bases;
            options.max_total_bases=cfg.local_max_total_bases;
            auto aligned=native_msa::align(dna,options);
            if(!plausible(aligned,dna,cfg.min_reference_identity)){
                // A rescued query must never veto an already supported
                // primary alignment: retry the primary-only read set.
                if(std::any_of(rescued.begin(),rescued.end(),[](bool x){return x;})){
                    std::vector<std::size_t> base_chosen;
                    std::vector<Pos> base_begins,base_ends;
                    std::vector<native_msa::Sequence> base_dna;
                    for(std::size_t j=0;j<chosen.size();++j){
                        if(rescued[j])continue;
                        base_chosen.push_back(chosen[j]);
                        base_begins.push_back(begins[j]);
                        base_ends.push_back(ends[j]);
                        base_dna.push_back(dna[j]);
                    }
                    if(base_dna.size()>=2){
                        auto base_aligned=native_msa::align(base_dna,options);
                        if(plausible(base_aligned,base_dna,cfg.min_reference_identity)){
                            chosen=std::move(base_chosen);
                            begins=std::move(base_begins);
                            ends=std::move(base_ends);
                            dna=std::move(base_dna);
                            rescued.assign(chosen.size(),false);
                            aligned=std::move(base_aligned);
                        }
                    }
                }
            }
            if(!plausible(aligned,dna,cfg.min_reference_identity)){
                ++stats.rejected_postcheck;
                stage.privateRegion(0,ref1,win_idx,b+1,"POSTCHECK_REJECTED");continue;
            }
            // Private query intervals BEFORE the accepted block, never promoted as homology.
            for(std::size_t i=1;i<chosen.size();++i)
                stage.privateRegion(chosen[i],begins[i],win_idx,b+1,"QUERY_PREFIX_UNRESOLVED");
            stage.accepted(win_idx,b+1,chosen,begins,ends,aligned,ref0,ref1,rescued);
            for(bool x:rescued)if(x)++stats.rescued_local_rows;
            ++stats.aligned_blocks;stats.aligned_reference_bases+=ref1-ref0;
        }
    }
    for(std::size_t r=0;r<fasta.reads().size();++r)
        stage.privateRegion(r,fasta.reads()[r].bases,win_idx+1,0,"READ_TAIL_UNRESOLVED");
    stats.private_bases=stage.privateBases();stats.output_columns=stage.columns();
    stage.publish();
    return stats;
}
} // namespace ultralong
