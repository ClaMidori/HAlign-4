#include "ultralong/native_pipeline.hpp"
#include "align/homology_plan.hpp"
#include "align/native_msa.hpp"

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <sstream>
#include <unordered_set>
#include <vector>

namespace ultralong_native {
namespace {
namespace fs = std::filesystem;
using Reads = std::vector<homology_plan::Sequence>;

struct Candidate {
    std::size_t qa = 0, qb = 0, seeds = 0;
    bool eligible = false;
};
struct Block {
    std::size_t ra=0,rb=0;
    std::vector<Candidate> query;
};

std::vector<std::string> tabs(const std::string& s) {
    std::vector<std::string> result;
    std::size_t from=0;
    while(true) {
        auto next=s.find('\t',from);
        result.push_back(s.substr(from,next==std::string::npos?next:next-from));
        if(next==std::string::npos)break;
        from=next+1;
    }
    return result;
}
std::size_t number(const std::string& s) {
    if(s.empty() || s.size()>18 || (s.size()>1 && s.front()=='0'))
        throw std::runtime_error("noncanonical numeric planner cell");
    std::size_t n=0;
    for(char c:s){
        if(c<'0'||c>'9'||n>(std::numeric_limits<std::size_t>::max()-static_cast<std::size_t>(c-'0'))/10)
            throw std::runtime_error("invalid numeric planner cell");
        n=n*10+static_cast<std::size_t>(c-'0');
    }
    return n;
}
std::vector<Block> blocksFromPlan(const homology_plan::Plan& plan, const Reads& reads) {
    const auto text=homology_plan::asTsv(plan); // applies pre-existing PR3D.2B candidate veto
    std::istringstream stream(text);
    std::string line;
    constexpr const char* schema = "block\tread_id\tref_start\tref_end\tquery_start_est\tquery_end_est\tleft_anchor_distance\tright_anchor_distance\tsupported_seed_count\tstatus";
    bool header=false;
    std::vector<Block> blocks(plan.boundaries.size()-1);
    std::vector<std::vector<bool>> found(blocks.size(),std::vector<bool>(reads.size(),false));
    while(std::getline(stream,line)) {
        if(line.empty()||line[0]=='#')continue;
        if(!header){if(line!=schema)throw std::runtime_error("unexpected planner TSV schema");header=true;continue;}
        auto fields=tabs(line);
        if(fields.size()!=10)throw std::runtime_error("invalid planner TSV row");
        auto b=number(fields[0]);
        if(b<1||b>blocks.size())throw std::runtime_error("planner block index invalid");
        auto& block=blocks[b-1];
        if(block.query.empty())block.query.resize(reads.size()-1);
        std::size_t read=reads.size();
        for(std::size_t j=0;j<reads.size();++j)if(fields[1]==reads[j].id){read=j;break;}
        if(read==reads.size() || found[b-1][read])throw std::runtime_error("planner unknown/duplicate read");
        found[b-1][read]=true;
        const auto ra=number(fields[2]),rb=number(fields[3]);
        if(ra!=plan.boundaries[b-1].reference || rb!=plan.boundaries[b].reference || ra>=rb)
            throw std::runtime_error("planner reference interval mismatch");
        block.ra=ra;block.rb=rb;
        if(read==0){
            if(fields[9]!="REFERENCE_ONLY" || number(fields[4])!=ra || number(fields[5])!=rb)
                throw std::runtime_error("planner reference row invalid");
        }else{
            if(fields[9]!="CANDIDATE_ONLY" && fields[9]!="UNRESOLVED")
                throw std::runtime_error("planner row status invalid");
            auto& q=block.query[read-1];
            if(fields[9]=="CANDIDATE_ONLY") {
                q.qa=number(fields[4]);q.qb=number(fields[5]);q.seeds=number(fields[8]);
                if(q.qa>=q.qb||q.qb>reads[read].bases.size()||q.seeds<3)
                    throw std::runtime_error("planner candidate exceeds read bounds");
                q.eligible=true;
            }
        }
    }
    if(!header)throw std::runtime_error("planner TSV header missing");
    for(const auto& row:found)for(bool v:row)if(!v)throw std::runtime_error("incomplete planner TSV");
    for(std::size_t i=1;i<blocks.size();++i)if(blocks[i].ra!=blocks[i-1].rb)
        throw std::runtime_error("planner reference blocks not contiguous");
    return blocks;
}

class Stage {
public:
    Stage(const fs::path& parent, std::size_t reads) {
        const auto stamp=std::chrono::steady_clock::now().time_since_epoch().count();
        for(std::size_t i=0;i<100;++i){
            auto candidate=parent/(".halign4_pr3d4c_"+std::to_string(stamp)+"_"+std::to_string(i));
            std::error_code ec;
            if(fs::create_directory(candidate,ec)){dir_=candidate;break;}
            if(ec)throw std::runtime_error("cannot create private staging directory: "+ec.message());
        }
        if(dir_.empty())throw std::runtime_error("cannot allocate unique staging directory");
        streams_.reserve(reads);
        paths_.reserve(reads);
        for(std::size_t i=0;i<reads;++i){
            auto path=dir_/("row_"+std::to_string(i));
            paths_.push_back(path);
            streams_.emplace_back(path,std::ios::binary);
            if(!streams_.back())throw std::runtime_error("cannot open temporary output spool");
        }
    }
    ~Stage(){std::error_code ec;fs::remove_all(dir_,ec);}
    Stage(const Stage&)=delete;Stage& operator=(const Stage&)=delete;
    void add(std::size_t row, const std::string& s){
        streams_.at(row).write(s.data(),static_cast<std::streamsize>(s.size()));
        if(!streams_[row])throw std::runtime_error("spool write failed");
    }
    void gaps(std::size_t row,std::size_t n){
        static const std::string eightk(8192,'-');
        while(n){auto part=std::min(n,eightk.size());streams_.at(row).write(eightk.data(),static_cast<std::streamsize>(part));n-=part;}
        if(!streams_[row])throw std::runtime_error("gap spool write failed");
    }
    void finish(){for(auto& f:streams_){f.flush();if(!f)throw std::runtime_error("spool flush failed");f.close();}}
    const fs::path& row(std::size_t i) const{return paths_.at(i);}
    const fs::path& dir()const{return dir_;}
private:
    fs::path dir_;
    std::vector<fs::path> paths_;
    std::vector<std::ofstream> streams_;
};

bool plausible(const native_msa::Alignment& msa,const std::vector<native_msa::Sequence>& local,double minimum) {
    if(msa.rows.size()!=local.size()||msa.columns==0)return false;
    for(std::size_t r=0;r<local.size();++r)
        if(msa.rows[r].id!=local[r].id||msa.rows[r].bases.size()!=msa.columns)return false;
    const auto& ref=msa.rows[0].bases;
    for(std::size_t j=1;j<local.size();++j){
        const auto& q=msa.rows[j].bases;
        std::size_t paired=0,match=0;
        for(std::size_t c=0;c<msa.columns;++c){
            if(ref[c]!='-'&&q[c]!='-'){
                ++paired;
                if(ref[c]==q[c]&&ref[c]!='N')++match;
            }
        }
        const auto shorter=std::min(local[0].bases.size(),local[j].bases.size());
        if(paired*100<shorter*65 || paired==0 || double(match)/paired<minimum)return false;
    }
    return true;
}

void makeFinal(const Stage& stage,const Reads& source,const fs::path& output,std::size_t expected_width){
    std::ofstream f(output,std::ios::binary|std::ios::trunc);
    if(!f)throw std::runtime_error("cannot create staged aligned FASTA");
    for(std::size_t r=0;r<source.size();++r){
        f<<'>'<<source[r].id<<'\n';
        std::ifstream row(stage.row(r),std::ios::binary);
        if(!row)throw std::runtime_error("cannot reopen aligned row spool");
        std::size_t width=0,bases=0,line=0;
        char c;
        while(row.get(c)){
            if(c!='-'){
                if(bases>=source[r].bases.size()||source[r].bases[bases]!=c)
                    throw std::runtime_error("FASTA reconstruction failed for "+source[r].id);
                ++bases;
            }
            f.put(c);++width;
            if(++line==80){f.put('\n');line=0;}
        }
        if(row.bad()||width!=expected_width||bases!=source[r].bases.size())
            throw std::runtime_error("aligned FASTA not lossless or rectangular: "+source[r].id);
        if(line!=0)f.put('\n');
    }
    f.flush();
    if(!f)throw std::runtime_error("aligned FASTA write failed");
}
} // namespace

Summary run(const fs::path& input,const fs::path& output,const fs::path& report,Config cfg){
    if(cfg.max_total_bases<2||cfg.max_total_bases>8000000||cfg.max_local_bases<1||
       cfg.max_local_bases>12000||!(cfg.min_reference_identity>0 && cfg.min_reference_identity<=1))
        throw std::invalid_argument("native preview budget/identity out of range");
    if(output.empty()||report.empty())throw std::invalid_argument("missing output/report");
    const auto out=fs::absolute(output).lexically_normal();
    const auto side=fs::absolute(report).lexically_normal();
    const auto src=fs::absolute(input).lexically_normal();
    if(out==side||out==src||side==src||out.parent_path()!=side.parent_path())
        throw std::invalid_argument("input/output/report conflict or output/report in different directories");
    if(!fs::is_directory(out.parent_path())||fs::exists(out)||fs::is_symlink(out)||
       fs::exists(side)||fs::is_symlink(side))
        throw std::invalid_argument("output and report must be new paths in an existing directory");
    const homology_plan::Limits limits{8,1000000,cfg.max_total_bases};
    const auto source=homology_plan::readFasta(input,limits);
    auto pc=homology_plan::Config{};
    pc.target_block_bases=5000;
    pc.strict_candidate_guard=true;
    pc.recover_indel_flanks=true;
    const auto plan=homology_plan::build(source,pc,limits);
    const auto blocks=blocksFromPlan(plan,source);
    Stage stage(out.parent_path(),source.size());
    std::ofstream journal(stage.dir()/"blocks.tsv",std::ios::binary);
    if(!journal)throw std::runtime_error("cannot write staged confidence TSV");
    journal<<"status\tblock\tread_id\tquery_start\tquery_end\tref_start\tref_end\tmsa_columns\tnote\n";
    Summary stats;
    for(const auto& s:source)stats.original_bases+=s.bases.size();
    std::vector<std::size_t> cursor(source.size(),0);
    auto privateSegment=[&](std::size_t read,std::size_t to){
        if(to<cursor[read]||to>source[read].bases.size())throw std::runtime_error("private interval invalid");
        const auto length=to-cursor[read];
        if(length==0)return;
        for(std::size_t i=0;i<source.size();++i){
            if(i==read)stage.add(i,source[i].bases.substr(cursor[i],length));
            else stage.gaps(i,length);
        }
        stats.private_bases+=length;
        stats.output_columns+=length;
        cursor[read]=to;
    };
    for(std::size_t index=0;index<blocks.size();++index){
        const auto& block=blocks[index];
        if(block.ra!=cursor[0]||block.rb>source[0].bases.size())throw std::runtime_error("reference cursor drift");
        std::vector<native_msa::Sequence> local{{source[0].id,source[0].bases.substr(block.ra,block.rb-block.ra)}};
        std::vector<std::size_t> selected{0};
        std::vector<std::size_t> begins{block.ra};
        std::vector<std::size_t> ends{block.rb};
        std::size_t local_bases=block.rb-block.ra;
        for(std::size_t r=1;r<source.size();++r){
            const auto& c=block.query.at(r-1);
            if(!c.eligible)continue;
            ++stats.candidate_query_rows;
            if(c.qa<cursor[r]||c.qb<=c.qa)continue;
            if(c.qb-c.qa>cfg.max_local_bases||block.rb-block.ra>cfg.max_local_bases||
               local_bases+c.qb-c.qa>100000){++stats.skipped_budget;continue;}
            local.push_back({source[r].id,source[r].bases.substr(c.qa,c.qb-c.qa)});
            selected.push_back(r);begins.push_back(c.qa);ends.push_back(c.qb);
            local_bases+=c.qb-c.qa;
        }
        if(local.size()==1 || local[0].bases.size()>cfg.max_local_bases){
            journal<<"UNRESOLVED_REF\t"<<index+1<<'\t'<<source[0].id<<"\tNA\tNA\t"<<block.ra<<'\t'<<block.rb<<"\t0\tNO_SAFE_LOCAL_CANDIDATE\n";
            privateSegment(0,block.rb);
            continue;
        }
        native_msa::Options opts;
        opts.max_sequences=source.size();opts.max_sequence_bases=cfg.max_local_bases;
        opts.max_total_bases=100000;
        auto msa=native_msa::align(local,opts);
        if(!plausible(msa,local,cfg.min_reference_identity)){
            ++stats.rejected_postcheck;
            journal<<"UNRESOLVED_REF\t"<<index+1<<'\t'<<source[0].id<<"\tNA\tNA\t"<<block.ra<<'\t'<<block.rb<<"\t0\tPOSTCHECK_REJECTED\n";
            privateSegment(0,block.rb);
            continue;
        }
        // Private query prefixes are emitted exclusively before a supported block.
        for(std::size_t i=1;i<selected.size();++i)privateSegment(selected[i],begins[i]);
        for(std::size_t r=0;r<source.size();++r){
            auto it=std::find(selected.begin(),selected.end(),r);
            if(it==selected.end())stage.gaps(r,msa.columns);
            else stage.add(r,msa.rows[static_cast<std::size_t>(it-selected.begin())].bases);
        }
        for(std::size_t i=0;i<selected.size();++i){
            cursor[selected[i]]=ends[i];
            journal<<"ANCHORED_LOCAL_BACKEND\t"<<index+1<<'\t'<<source[selected[i]].id<<'\t'
                   <<begins[i]<<'\t'<<ends[i]<<'\t'<<block.ra<<'\t'<<block.rb<<'\t'
                   <<msa.columns<<"\tHEURISTIC_NOT_VERIFIED\n";
        }
        stats.output_columns+=msa.columns;
        stats.accepted_blocks+=1;
        stats.aligned_reference_bases+=block.rb-block.ra;
    }
    for(std::size_t r=0;r<source.size();++r)privateSegment(r,source[r].bases.size());
    journal.flush();if(!journal)throw std::runtime_error("confidence sidecar flush failed");journal.close();
    stage.finish();
    makeFinal(stage,source,stage.dir()/"aligned.fasta",stats.output_columns);
    // Same filesystem + hard links: no-replace publication, including races.
    fs::create_hard_link(stage.dir()/"aligned.fasta",out);
    try {fs::create_hard_link(stage.dir()/"blocks.tsv",side);}
    catch(...){std::error_code ec;fs::remove(out,ec);throw;}
    return stats;
}
} // namespace ultralong_native
