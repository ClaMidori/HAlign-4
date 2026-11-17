#include "align/lossless_extractor.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cerrno>
#include <system_error>
#ifdef __linux__
#include <fcntl.h>
#include <linux/fs.h>
#include <sys/syscall.h>
#include <unistd.h>
#endif
#include <cctype>
#include <tuple>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

namespace lossless_extractor {
namespace {
namespace fs = std::filesystem;

// Local SHA-256, avoiding a shell command and an additional runtime dependency.
class Sha256 {
public:
    void update(std::string_view str) {
        for (unsigned char c : str) {
            buffer_[buffer_size_++] = c;
            if (buffer_size_ == 64) { process(); buffer_size_ = 0; }
            if (bytes_ == std::numeric_limits<std::uint64_t>::max())
                throw std::overflow_error("SHA-256 input too large");
            ++bytes_;
        }
    }
    std::string finish() {
        const std::uint64_t bits = bytes_ * 8U;
        buffer_[buffer_size_++] = 0x80;
        if (buffer_size_ > 56) {
            while (buffer_size_ < 64) buffer_[buffer_size_++] = 0;
            process(); buffer_size_ = 0;
        }
        while (buffer_size_ < 56) buffer_[buffer_size_++] = 0;
        for (int s = 56; s >= 0; s -= 8) buffer_[buffer_size_++] = static_cast<std::uint8_t>(bits >> s);
        process();
        std::ostringstream o;
        o << std::hex << std::setfill('0');
        for (const auto x : h_) o << std::setw(8) << x;
        return o.str();
    }
private:
    static constexpr std::array<std::uint32_t, 64> k_ = {
        0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
        0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
        0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
        0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
        0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
        0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
        0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
        0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2};
    std::array<std::uint32_t,8> h_ = {0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,
                                       0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19};
    std::array<std::uint8_t,64> buffer_{};
    std::size_t buffer_size_ = 0;
    std::uint64_t bytes_ = 0;
    static std::uint32_t rr(std::uint32_t x, unsigned n) {return (x >> n) | (x << (32U-n));}
    void process() {
        std::uint32_t w[64]{};
        for (int i=0;i<16;++i) {
            const auto j=4*i;
            w[i]=(std::uint32_t(buffer_[j])<<24)|(std::uint32_t(buffer_[j+1])<<16)|
                 (std::uint32_t(buffer_[j+2])<<8)|std::uint32_t(buffer_[j+3]);
        }
        for (int i=16;i<64;++i) {
            const auto s0=rr(w[i-15],7)^rr(w[i-15],18)^(w[i-15]>>3);
            const auto s1=rr(w[i-2],17)^rr(w[i-2],19)^(w[i-2]>>10);
            w[i]=w[i-16]+s0+w[i-7]+s1;
        }
        auto [a,b,c,d,e,f,g,h]=std::tuple{h_[0],h_[1],h_[2],h_[3],h_[4],h_[5],h_[6],h_[7]};
        for (int i=0;i<64;++i) {
            const auto s1=rr(e,6)^rr(e,11)^rr(e,25);
            const auto choice=(e&f)^(~e&g);
            const auto tmp1=h+s1+choice+k_[i]+w[i];
            const auto s0=rr(a,2)^rr(a,13)^rr(a,22);
            const auto maj=(a&b)^(a&c)^(b&c);
            const auto tmp2=s0+maj;
            h=g;g=f;f=e;e=d+tmp1;d=c;c=b;b=a;a=tmp1+tmp2;
        }
        h_[0]+=a;h_[1]+=b;h_[2]+=c;h_[3]+=d;
        h_[4]+=e;h_[5]+=f;h_[6]+=g;h_[7]+=h;
    }
};
std::string sha256(std::string_view s) {Sha256 h;h.update(s);return h.finish();}
std::string sha256File(const fs::path& path) {
    std::ifstream f(path,std::ios::binary);
    if(!f) throw std::runtime_error("cannot read input: "+path.string());
    Sha256 s;
    std::array<char,8192> buf{};
    while (f) {
        f.read(buf.data(),buf.size());
        if(f.gcount()>0) s.update(std::string_view(buf.data(),static_cast<std::size_t>(f.gcount())));
    }
    if(f.bad()) throw std::runtime_error("file read failed: "+path.string());
    return s.finish();
}

bool utf8Valid(std::string_view s) {
    for (std::size_t i=0;i<s.size();) {
        const auto c=static_cast<std::uint8_t>(s[i]);
        unsigned n=0; std::uint32_t cp=0;
        if(c<0x80){++i;continue;}
        if(c>=0xc2 && c<=0xdf){n=2;cp=c&0x1f;}
        else if(c>=0xe0 && c<=0xef){n=3;cp=c&0x0f;}
        else if(c>=0xf0 && c<=0xf4){n=4;cp=c&7;}
        else return false;
        if(i+n>s.size())return false;
        for(unsigned j=1;j<n;++j){const auto z=static_cast<std::uint8_t>(s[i+j]);if((z&0xc0)!=0x80)return false;cp=(cp<<6)|(z&0x3f);}
        if((n==2 && cp<0x80)||(n==3 && cp<0x800)||(n==4 && cp<0x10000)||
           cp>0x10ffff||(cp>=0xd800 && cp<=0xdfff)) return false;
        i+=n;
    }
    return true;
}
struct Read {std::string id,bases;};
std::vector<Read> parseFasta(const fs::path& path) {
    std::ifstream f(path,std::ios::binary);
    if(!f)throw std::runtime_error("cannot open input FASTA: "+path.string());
    std::vector<Read> reads;
    std::set<std::string> ids;
    std::string header;
    std::size_t total=0;
    bool start=true,in_header=false;
    char raw=0;
    auto finishHeader=[&]() {
        const auto begin=header.find_first_not_of(" \t\r\n\v\f");
        if(begin==std::string::npos)throw std::invalid_argument("empty FASTA ID");
        const auto end=header.find_first_of(" \t\r\n\v\f",begin);
        const auto name=header.substr(begin,end==std::string::npos?std::string::npos:end-begin);
        if(!utf8Valid(name)||name.empty())throw std::invalid_argument("invalid UTF-8 FASTA ID");
        if(!ids.insert(name).second)throw std::invalid_argument("duplicate FASTA ID: "+name);
        reads.back().id=name;
    };
    while(f.get(raw)) {
        if(raw=='\n'){
            if(in_header) {finishHeader();in_header=false;header.clear();}
            start=true;continue;
        }
        if(raw=='\r')continue;
        if(start){
            start=false;
            if(raw=='>'){
                if(in_header)throw std::invalid_argument("nested FASTA header");
                if(reads.size()>=8)throw std::invalid_argument("more than 8 FASTA records");
                reads.push_back({{}, {}});
                in_header=true;header.clear();continue;
            }
            if(reads.empty())throw std::invalid_argument("FASTA sequence before first header");
        }
        if(in_header){
            header.push_back(raw);
            if(header.size()>256)throw std::invalid_argument("FASTA header exceeds 256 bytes");
        }else{
            const char c=static_cast<char>(std::toupper(static_cast<unsigned char>(raw)));
            if(c!='A'&&c!='C'&&c!='G'&&c!='T'&&c!='N')
                throw std::invalid_argument("unsupported FASTA base");
            if(reads.back().bases.size()>=1000000||total>=4000000)
                throw std::invalid_argument("diagnostic FASTA size budget exceeded");
            reads.back().bases.push_back(c);++total;
        }
    }
    if(f.bad())throw std::runtime_error("FASTA I/O failure");
    if(in_header)finishHeader();
    if(reads.size()<2||reads.size()>8)throw std::invalid_argument("FASTA requires 2..8 reads");
    for(const auto& r:reads)if(r.bases.empty())throw std::invalid_argument("empty FASTA sequence: "+r.id);
    return reads;
}
std::vector<std::string> splitTabs(const std::string& line){
    std::vector<std::string> v;std::size_t at=0;
    while(true){const auto p=line.find('\t',at);v.push_back(line.substr(at,p==std::string::npos?p:p-at));
                if(p==std::string::npos)break;
                at=p+1;
    }
    return v;
}
std::size_t integer(const std::string& s,const std::string& label){
    if(s.empty()||(s.size()>1 && s[0]=='0'))throw std::invalid_argument("noncanonical integer for "+label);
    std::size_t value=0;
    for(char c:s){
        if(c<'0'||c>'9'||value>(std::numeric_limits<std::size_t>::max()-static_cast<unsigned>(c-'0'))/10)
            throw std::invalid_argument("invalid integer for "+label);
        value=value*10+static_cast<unsigned>(c-'0');
    }
    return value;
}
struct PlanRow {
    std::size_t a=0,b=0,qa=0,qb=0,seeds=0;
    bool candidate=false;
};
using Block=std::vector<PlanRow>;
std::vector<Block> parsePlan(const fs::path& path,const std::vector<Read>& reads) {
    if(fs::file_size(path)>16*1024*1024)
        throw std::invalid_argument("diagnostic TSV exceeds 16MiB safety budget");
    std::ifstream f(path,std::ios::binary);
    if(!f)throw std::runtime_error("cannot open diagnostic TSV: "+path.string());
    std::unordered_map<std::string,std::size_t> ids;
    for(std::size_t i=0;i<reads.size();++i)ids.emplace(reads[i].id,i);
    std::map<std::size_t,Block> by_block;
    std::map<std::size_t,std::vector<bool>> seen;
    const std::string schema="block\tread_id\tref_start\tref_end\tquery_start_est\tquery_end_est\tleft_anchor_distance\tright_anchor_distance\tsupported_seed_count\tstatus";
    std::string line;bool first=true,columns=false,reference_metadata=false;
    while(std::getline(f,line)){
        if(!line.empty()&&line.back()=='\r')line.pop_back();
        if(first){first=false;if(line.find("NO alignment and NO verified homologous coordinates")==std::string::npos)
            throw std::invalid_argument("not a diagnostic-only planner report");}
        if(line.empty())continue;
        if(line[0]=='#'){
            if(line.rfind("# reference=",0)==0){
                const auto fields=splitTabs(line.substr(2));
                if(fields.empty()||fields[0]!="reference="+reads[0].id ||
                    std::find(fields.begin(),fields.end(),"length="+std::to_string(reads[0].bases.size()))==fields.end())
                    throw std::invalid_argument("TSV reference metadata disagrees with FASTA");
                reference_metadata=true;
            }
            continue;
        }
        if(!columns){
            if(line!=schema)throw std::invalid_argument("unexpected 10-column diagnostic TSV schema");
            columns=true;continue;
        }
        const auto col=splitTabs(line);
        if(col.size()!=10)throw std::invalid_argument("malformed TSV row");
        const auto idx=integer(col[0],"block");
        if(idx<1||idx>100000)throw std::invalid_argument("invalid block number");
        auto [block_it,created]=by_block.try_emplace(idx,reads.size());
        if(created)seen.emplace(idx,std::vector<bool>(reads.size(),false));
        auto it=ids.find(col[1]);
        if(it==ids.end()||seen.at(idx)[it->second])
            throw std::invalid_argument("unknown/duplicate TSV read_id");
        auto& r=block_it->second[it->second];
        seen.at(idx)[it->second]=true;
        r.a=integer(col[2],"ref_start");r.b=integer(col[3],"ref_end");
        if(r.a>=r.b||r.b>reads[0].bases.size())throw std::invalid_argument("reference interval out of bounds");
        const auto j=it->second;
        if(j==0){
            if(col[9]!="REFERENCE_ONLY"||col[4]!=col[2]||col[5]!=col[3])
                throw std::invalid_argument("invalid reference row");
            r.qa=r.a;r.qb=r.b;
        }else{
            if(col[9]!="UNRESOLVED"&&col[9]!="CANDIDATE_ONLY")
                throw std::invalid_argument("unsupported query status");
            r.candidate=col[9]=="CANDIDATE_ONLY";
            if(col[4]=="NA"||col[5]=="NA"){
                if(r.candidate)throw std::invalid_argument("candidate lacks projection");
            } else {
                r.qa=integer(col[4],"query_start_est");r.qb=integer(col[5],"query_end_est");
                if(r.qa>reads[j].bases.size()||r.qb>reads[j].bases.size()||
                    (r.candidate&&r.qa>=r.qb))throw std::invalid_argument("query interval out of bounds");
            }
            r.seeds=integer(col[8],"supported_seed_count");
            if(r.candidate){
                if(r.seeds<3)throw std::invalid_argument("candidate lacks seed support");
                if(integer(col[6],"left_anchor_distance")>500||integer(col[7],"right_anchor_distance")>500)
                    throw std::invalid_argument("candidate lacks boundary support");
            }
        }
    }
    if(f.bad()||!columns||!reference_metadata||by_block.empty())
        throw std::invalid_argument("incomplete diagnostic TSV");
    if(by_block.rbegin()->first!=by_block.size())
        throw std::invalid_argument("block numbers must be contiguous from 1");
    std::vector<Block> blocks;
    blocks.reserve(by_block.size());
    std::size_t pos=0;
    for(const auto& [block_number,block]:by_block){
        for(std::size_t j=0;j<reads.size();++j){
            if(!seen.at(block_number)[j])throw std::invalid_argument("block missing sequence row");
            if(block[j].a!=block[0].a||block[j].b!=block[0].b)
                throw std::invalid_argument("inconsistent block reference boundaries");
        }
        if(block[0].a!=pos)throw std::invalid_argument("reference partition gap/overlap");
        pos=block[0].b;
        blocks.push_back(block);
    }
    if(pos!=reads[0].bases.size())throw std::invalid_argument("reference partition incomplete");
    return blocks;
}
struct Rejection {std::size_t block=0,read=0;std::set<std::string> why;};
struct Tile {std::size_t read=0,block=0,start=0,end=0;std::string kind,id,checksum;};
struct Pair {std::size_t block=0,read=0,a=0,b=0,qa=0,qb=0;};
std::string jsonEsc(std::string_view v){
    std::string s="\"";static constexpr char h[]="0123456789abcdef";
    for(unsigned char c:v){
        if(c=='\"'||c=='\\'){s+='\\';s+=static_cast<char>(c);}
        else if(c=='\b')s+="\\b";else if(c=='\f')s+="\\f";
        else if(c=='\n')s+="\\n";else if(c=='\r')s+="\\r";else if(c=='\t')s+="\\t";
        else if(c<32){s+="\\u00";s+=h[c>>4];s+=h[c&15];}
        else s+=static_cast<char>(c);
    }
    return s+"\"";
}
void checked(std::ofstream& f,const fs::path& p){f.flush();if(!f)throw std::runtime_error("failed writing "+p.string());}
void writeFasta(std::ofstream& f,std::string_view id,std::string_view bases){
    f<<'>'<<id<<'\n';
    for(std::size_t p=0;p<bases.size();p+=100)f<<bases.substr(p,100)<<'\n';
}
} // namespace

std::string sha256Digest(std::string_view bytes) { return sha256(bytes); }
std::string sha256FileDigest(const fs::path& path) { return sha256File(path); }

Summary extract(const fs::path& fasta,const fs::path& plan_tsv,const fs::path& output_dir){
    if(output_dir.empty())throw std::invalid_argument("empty output directory");
    const auto out=fs::absolute(output_dir).lexically_normal();
    if(fs::exists(out)||fs::is_symlink(out))
        throw std::invalid_argument("refusing existing/empty output directory");
    if(!fs::is_directory(out.parent_path()))
        throw std::invalid_argument("output parent directory does not exist");
    const auto reads=parseFasta(fasta);
    const auto blocks=parsePlan(plan_tsv,reads);
    std::vector<Pair> pairs;
    std::map<std::pair<std::size_t,std::size_t>,std::set<std::string>> rejected;
    std::vector<std::vector<std::pair<std::size_t,std::pair<std::size_t,std::size_t>>>> accepted(reads.size());
    for(std::size_t read=1;read<reads.size();++read){
        std::vector<std::size_t> candidates;
        for(std::size_t i=0;i<blocks.size();++i)if(blocks[i][read].candidate)candidates.push_back(i);
        for(std::size_t i=0;i<candidates.size();++i){
            const auto ai=candidates[i];const auto& a=blocks[ai][read];
            for(std::size_t j=i+1;j<candidates.size();++j){
                const auto bi=candidates[j];const auto& b=blocks[bi][read];
                const bool overlap=a.qa<b.qb && b.qa<a.qb;
                if(overlap){rejected[{ai,read}].insert("QUERY_OVERLAP");
                            rejected[{bi,read}].insert("QUERY_OVERLAP");}
                if(a.qa>=b.qa||a.qb>b.qa){
                    if(a.qb>b.qa && a.qa<b.qb)continue;
                    rejected[{ai,read}].insert("NON_COLLINEAR");
                    rejected[{bi,read}].insert("NON_COLLINEAR");
                }
            }
        }
        for(auto i:candidates)if(!rejected.contains({i,read})){
            const auto& r=blocks[i][read];
            accepted[read].push_back({i,{r.qa,r.qb}});
            const auto& ref=blocks[i][0];
            pairs.push_back({i+1,read,ref.a,ref.b,r.qa,r.qb});
        }
    }
    // Python writes associations in reference-block order, then FASTA query order.
    std::sort(pairs.begin(),pairs.end(),[](const auto& a,const auto& b){
        return a.block!=b.block?a.block<b.block:a.read<b.read;});
    std::vector<Tile> tiles;
    auto addTile=[&](std::size_t r,std::size_t a,std::size_t b,std::size_t block,const std::string& kind,std::size_t n){
        if(a>=b||b>reads[r].bases.size())throw std::logic_error("invalid tile coordinates");
        std::ostringstream id;id<<'t'<<std::setw(2)<<std::setfill('0')<<(r+1)<<'_'<<std::setw(5)<<n;
        tiles.push_back({r,block,a,b,kind,id.str(),sha256(std::string_view(reads[r].bases).substr(a,b-a))});
    };
    std::size_t unresolved=0;
    for(std::size_t r=0;r<reads.size();++r){
        std::vector<std::tuple<std::size_t,std::size_t,std::size_t,std::string>> spans;
        if(r==0){
            for(std::size_t i=0;i<blocks.size();++i)
                spans.emplace_back(blocks[i][0].a,blocks[i][0].b,i+1,"REFERENCE_TILE");
        }else{
            for(auto [idx,q]:accepted[r])spans.emplace_back(q.first,q.second,idx+1,"CANDIDATE_TILE");
            std::sort(spans.begin(),spans.end());
        }
        std::size_t pos=0,n=0,covered=0;
        for(const auto& [a,b,block,kind]:spans){
            if(a<pos||b>reads[r].bases.size())throw std::logic_error("overlapping accepted tiles");
            if(a>pos){addTile(r,pos,a,0,"UNRESOLVED",++n);unresolved+=a-pos;covered+=a-pos;}
            addTile(r,a,b,block,kind,++n);covered+=b-a;pos=b;
        }
        if(pos<reads[r].bases.size()){
            addTile(r,pos,reads[r].bases.size(),0,"UNRESOLVED",++n);
            unresolved+=reads[r].bases.size()-pos;
            covered+=reads[r].bases.size()-pos;
        }
        if(n==0||covered!=reads[r].bases.size())throw std::logic_error("lossless reconstruction failed");
    }
    Summary summary{reads.size(),blocks.size(),pairs.size(),rejected.size(),tiles.size(),unresolved};
    fs::path stage;
    // Create sibling staging directory; never overwrite output or source.
    for(unsigned i=0;i<64;++i){
        const auto ticks=std::chrono::steady_clock::now().time_since_epoch().count();
        stage=out.parent_path()/(".pr3d3a_native_"+std::to_string(ticks)+"_"+std::to_string(i));
        if(fs::create_directory(stage))break;
        if(i==63)throw std::runtime_error("cannot create unique staging directory");
    }
    try{
        auto create=[&](const char* name){
            const auto path=stage/name;
            std::ofstream f(path,std::ios::binary);
            if(!f)throw std::runtime_error("cannot write "+path.string());
            return f;
        };
        {
            auto f=create("tiles.fasta");
            for(const auto& t:tiles)
                writeFasta(f,t.id,std::string_view(reads[t.read].bases).substr(t.start,t.end-t.start));
            checked(f,stage/"tiles.fasta");
        }
        {
            auto f=create("tiles.tsv");
            f<<"tile_id\tread_id\tkind\tstart\tend\tblock\tlength\tsha256\n";
            for(const auto& t:tiles){
                f<<t.id<<'\t'<<reads[t.read].id<<'\t'<<t.kind<<'\t'<<t.start<<'\t'<<t.end<<'\t';
                if(t.block)f<<t.block;
                f<<'\t'<<(t.end-t.start)<<'\t'<<t.checksum<<'\n';
            }
            checked(f,stage/"tiles.tsv");
        }
        {
            auto f=create("candidate_pairs.tsv");
            f<<"block\tread_id\tref_start\tref_end\tquery_start\tquery_end\tstatus\n";
            for(const auto& p:pairs)f<<p.block<<'\t'<<reads[p.read].id<<'\t'<<p.a<<'\t'<<p.b
                <<'\t'<<p.qa<<'\t'<<p.qb<<"\tCANDIDATE_ONLY\n";
            checked(f,stage/"candidate_pairs.tsv");
        }
        {
            auto f=create("rejected_candidates.tsv");
            f<<"block\tread_id\treason\n";
            for(std::size_t i=0;i<blocks.size();++i)
                for(std::size_t r=1;r<reads.size();++r){
                    const auto it=rejected.find({i,r});if(it==rejected.end())continue;
                    f<<(i+1)<<'\t'<<reads[r].id<<'\t';
                    bool first=true;for(const auto& reason:it->second){if(!first)f<<';';f<<reason;first=false;}
                    f<<'\n';
                }
            checked(f,stage/"rejected_candidates.tsv");
        }
        {
            auto f=create("summary.json");
            f<<"{\n  \"schema\": \"pr3d3a-extraction-v1\",\n"
             <<"  \"diagnostic_only\": true,\n  \"no_msa\": true,\n"
             <<"  \"no_verified_homology\": true,\n"
             <<"  \"source_fasta_sha256\": \""<<sha256File(fasta)<<"\",\n"
             <<"  \"source_plan_sha256\": \""<<sha256File(plan_tsv)<<"\",\n"
             <<"  \"reference_id\": "<<jsonEsc(reads.front().id)<<",\n"
             <<"  \"reads\": [\n";
            for(std::size_t i=0;i<reads.size();++i){
                if(i)f<<",\n";
                f<<"    {\n      \"read_id\": "<<jsonEsc(reads[i].id)<<",\n"
                 <<"      \"length\": "<<reads[i].bases.size()<<",\n"
                 <<"      \"sha256\": \""<<sha256(reads[i].bases)<<"\"\n    }";
            }
            f<<"\n  ],\n  \"blocks\": "<<summary.blocks<<",\n"
             <<"  \"accepted_query_candidates\": "<<summary.accepted_query_candidates<<",\n"
             <<"  \"rejected_conflicting_candidates\": "<<summary.rejected_conflicting_candidates<<",\n"
             <<"  \"tiles\": "<<summary.tiles<<",\n"
             <<"  \"unresolved_query_bases\": "<<summary.unresolved_query_bases<<"\n}\n";
            checked(f,stage/"summary.json");
        }
        if(fs::exists(out)||fs::is_symlink(out))
            throw std::runtime_error("refusing output directory created concurrently");
#ifdef __linux__
        // Atomic no-overwrite publication. Unlike plain rename(), this refuses
        // a competing directory created after the initial existence check.
        if(::syscall(SYS_renameat2,AT_FDCWD,stage.c_str(),AT_FDCWD,
                     out.c_str(),RENAME_NOREPLACE)!=0)
            throw std::system_error(errno,std::generic_category(),
                                    "cannot atomically publish extracted directory");
#else
        fs::rename(stage,out);
#endif
    }catch(...){std::error_code ec;fs::remove_all(stage,ec);throw;}
    return summary;
}
} // namespace lossless_extractor
