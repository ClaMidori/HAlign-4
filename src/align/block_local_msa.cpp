#include "align/block_local_msa.hpp"
#include "align/homology_plan.hpp"
#include "align/lossless_extractor.hpp"
#include "align/native_msa.hpp"

#include <algorithm>
#include <array>
#include <cerrno>
#include <chrono>
#include <cctype>
#include <cstdint>
#include <fstream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>
#ifdef __linux__
#include <fcntl.h>
#include <linux/fs.h>
#include <sys/syscall.h>
#include <unistd.h>
#endif

namespace block_local_msa {
namespace {
namespace fs = std::filesystem;

struct Tile {
    std::string tile_id, read_id, kind, bases, hash;
    std::size_t start=0, end=0, block=0, read_index=0;
};
struct Candidate {
    std::size_t block=0, read_index=0, ref_start=0, ref_end=0, query_start=0, query_end=0;
};

std::vector<std::string> tabs(const std::string& line) {
    std::vector<std::string> v;
    std::size_t pos=0;
    while (true) {
        const auto next=line.find('\t',pos);
        v.push_back(line.substr(pos,next==std::string::npos?next:next-pos));
        if(next==std::string::npos)break;
        pos=next+1;
    }
    return v;
}
std::size_t num(const std::string& s) {
    if(s.empty()||(s.size()>1&&s[0]=='0'))throw std::invalid_argument("noncanonical TSV integer");
    std::size_t v=0;
    for(unsigned char c:s) {
        if(c<'0'||c>'9'||v>(std::numeric_limits<std::size_t>::max()-(c-'0'))/10)
            throw std::invalid_argument("malformed TSV integer");
        v=v*10+(c-'0');
    }
    return v;
}
std::string readSmall(const fs::path& p,std::uintmax_t limit) {
    if(!fs::is_regular_file(p)||fs::file_size(p)>limit)
        throw std::invalid_argument("missing or oversized extraction artifact: "+p.string());
    std::ifstream f(p,std::ios::binary);
    if(!f)throw std::runtime_error("cannot read extraction artifact: "+p.string());
    return {std::istreambuf_iterator<char>(f),std::istreambuf_iterator<char>()};
}
std::vector<std::string> lines(const std::string& data) {
    std::vector<std::string> out;
    std::istringstream s(data);
    std::string line;
    while(std::getline(s,line)) {
        if(!line.empty()&&line.back()=='\r')line.pop_back();
        if(line.empty())throw std::invalid_argument("blank lines forbidden in extraction TSV");
        out.push_back(std::move(line));
    }
    return out;
}
void writeFasta(std::ofstream& f,const std::string& id,const std::string& dna) {
    f<<'>'<<id<<'\n';
    for(std::size_t i=0;i<dna.size();i+=100)f<<dna.substr(i,100)<<'\n';
}
void flush(std::ofstream& f,const fs::path& p) {
    f.flush();if(!f)throw std::runtime_error("cannot write "+p.string());
}
void publish(const fs::path& stage,const fs::path& dest) {
    if(fs::exists(dest)||fs::is_symlink(dest))throw std::runtime_error("output directory exists");
#ifdef __linux__
    if(::syscall(SYS_renameat2,AT_FDCWD,stage.c_str(),AT_FDCWD,dest.c_str(),RENAME_NOREPLACE)!=0)
        throw std::system_error(errno,std::generic_category(),"atomic no-replace publication failed");
#else
    fs::rename(stage,dest);
#endif
}
std::string jsonEsc(std::string_view input) {
    std::string out="\"";
    constexpr char hex[]="0123456789abcdef";
    for(const unsigned char c:input) {
        if(c=='"'||c=='\\') {out+='\\';out+=static_cast<char>(c);}
        else if(c=='\n')out+="\\n";
        else if(c=='\r')out+="\\r";
        else if(c=='\t')out+="\\t";
        else if(c<32) {out+="\\u00";out+=hex[c>>4];out+=hex[c&15];}
        else out+=static_cast<char>(c);
    }
    return out+'"';
}
std::vector<Tile> readTiles(const fs::path& dir,const std::vector<homology_plan::Sequence>& source) {
    const auto table=lines(readSmall(dir/"tiles.tsv",16*1024*1024));
    if(table.empty()||table.front()!="tile_id\tread_id\tkind\tstart\tend\tblock\tlength\tsha256")
        throw std::invalid_argument("unexpected tiles.tsv schema");
    if(table.size()>100001)throw std::invalid_argument("too many extracted tiles");
    std::unordered_map<std::string,std::size_t> read_index;
    for(std::size_t i=0;i<source.size();++i)read_index.emplace(source[i].id,i);
    std::vector<Tile> tiles;
    std::set<std::string> ids;
    std::vector<std::size_t> next(source.size(),0);
    std::size_t last_read=0;
    for(std::size_t i=1;i<table.size();++i) {
        const auto c=tabs(table[i]);
        if(c.size()!=8)throw std::invalid_argument("malformed tiles.tsv row");
        const auto it=read_index.find(c[1]);
        if(it==read_index.end()||!ids.insert(c[0]).second||c[0].empty()||c[0].size()>32)
            throw std::invalid_argument("unknown read or duplicate/invalid tile ID");
        const auto rid=it->second;
        if(rid<last_read)throw std::invalid_argument("tile FASTA order disagrees with input FASTA");
        last_read=rid;
        const auto a=num(c[3]),b=num(c[4]);
        const auto block=c[5].empty()?0:num(c[5]);
        if(a!=next[rid]||b<=a||b>source[rid].bases.size()||num(c[6])!=b-a)
            throw std::invalid_argument("extracted tiles do not partition each read without loss");
        if((rid==0 && (c[2]!="REFERENCE_TILE"||block==0)) ||
           (rid>0 && c[2]=="REFERENCE_TILE") ||
           (c[2]=="UNRESOLVED" && block!=0) ||
           (c[2]=="CANDIDATE_TILE" && (rid==0||block==0)) ||
           (c[2]!="REFERENCE_TILE"&&c[2]!="CANDIDATE_TILE"&&c[2]!="UNRESOLVED") ||
           (rid>0&&c[2]!="UNRESOLVED"&&c[2]!="CANDIDATE_TILE"))
            throw std::invalid_argument("invalid extracted tile status or block");
        if(c[7].size()!=64||std::any_of(c[7].begin(),c[7].end(),[](char h){return !std::isxdigit(static_cast<unsigned char>(h));}))
            throw std::invalid_argument("invalid tile hash");
        next[rid]=b;
        tiles.push_back({c[0],c[1],c[2],{},c[7],a,b,block,rid});
    }
    for(std::size_t r=0;r<source.size();++r)
        if(next[r]!=source[r].bases.size())throw std::invalid_argument("missing trailing sequence tiles");

    // Read the FASTA incrementally, rejecting malformed sequence rows and IDs.
    const auto fasta=dir/"tiles.fasta";
    if(!fs::is_regular_file(fasta)||fs::file_size(fasta)>6*1024*1024)
        throw std::invalid_argument("missing or oversized tiles.fasta");
    std::ifstream f(fasta,std::ios::binary);
    if(!f)throw std::runtime_error("cannot open tiles.fasta");
    std::size_t idx=0;
    std::string line;
    while(std::getline(f,line)) {
        if(!line.empty()&&line.back()=='\r')line.pop_back();
        if(line.empty())throw std::invalid_argument("empty FASTA line");
        if(line.front()=='>') {
            if(idx>=tiles.size()||line.substr(1)!=tiles[idx].tile_id)
                throw std::invalid_argument("tile FASTA headers/order disagree with tiles.tsv");
            if(idx>0 && tiles[idx-1].bases.size()!=tiles[idx-1].end-tiles[idx-1].start)
                throw std::invalid_argument("incomplete tile FASTA");
            ++idx;
        } else {
            if(idx==0||line.size()>100||line.find_first_not_of("ACGTN")!=std::string::npos ||
               line.size()>tiles[idx-1].end-tiles[idx-1].start-tiles[idx-1].bases.size())
                throw std::invalid_argument("invalid tile FASTA sequence");
            tiles[idx-1].bases+=line;
        }
    }
    if(f.bad()||idx!=tiles.size())throw std::invalid_argument("incomplete tiles.fasta");
    for(const auto& t:tiles) {
        if(t.bases.size()!=t.end-t.start ||
           t.bases!=source[t.read_index].bases.substr(t.start,t.end-t.start) ||
           lossless_extractor::sha256Digest(t.bases)!=t.hash)
            throw std::invalid_argument("tile checksum or source FASTA mismatch: "+t.tile_id);
    }
    return tiles;
}
} // namespace

Summary align(const fs::path& original_fasta,const fs::path& extraction_dir,const fs::path& output_dir) {
    if(output_dir.empty()||!fs::is_directory(extraction_dir))
        throw std::invalid_argument("missing extracted tile directory or output");
    const auto out=fs::absolute(output_dir).lexically_normal();
    if(fs::exists(out)||fs::is_symlink(out)||!fs::is_directory(out.parent_path()))
        throw std::invalid_argument("output must be a new directory with an existing parent");
    const auto source=homology_plan::readFasta(original_fasta);
    const auto manifest=readSmall(extraction_dir/"summary.json",128*1024);
    const auto source_hash=lossless_extractor::sha256FileDigest(original_fasta);
    if(manifest.find("\"schema\": \"pr3d3a-extraction-v1\"")==std::string::npos ||
       manifest.find("\"diagnostic_only\": true")==std::string::npos ||
       manifest.find("\"no_msa\": true")==std::string::npos ||
       manifest.find("\"no_verified_homology\": true")==std::string::npos ||
       manifest.find("\"source_fasta_sha256\": \""+source_hash+"\"")==std::string::npos)
        throw std::invalid_argument("source FASTA differs from extraction or manifest is not diagnostic");

    auto tiles=readTiles(extraction_dir,source);
    std::map<std::size_t,const Tile*> refs;
    std::map<std::pair<std::size_t,std::size_t>,const Tile*> candidates;
    std::size_t unresolved=0;
    for(const auto& tile:tiles){
        if(tile.kind=="REFERENCE_TILE"){
            if(!refs.emplace(tile.block,&tile).second)throw std::invalid_argument("duplicate reference block tile");
        }else if(tile.kind=="CANDIDATE_TILE"){
            if(!candidates.emplace(std::make_pair(tile.block,tile.read_index),&tile).second)
                throw std::invalid_argument("duplicate query block tile");
        }else ++unresolved;
    }
    if(refs.empty()||refs.size()>100000 || refs.rbegin()->first!=refs.size())
        throw std::invalid_argument("reference blocks are not contiguous");
    std::size_t reference_pos=0;
    for(const auto& [block,tile]:refs){
        if(tile->start!=reference_pos)throw std::invalid_argument("reference block order/coordinates invalid");
        reference_pos=tile->end;
    }
    if(reference_pos!=source[0].bases.size())throw std::invalid_argument("reference coverage incomplete");
    for(const auto& [key,tile]:candidates)
        if(!refs.contains(key.first))throw std::invalid_argument("candidate references missing block");

    const auto pairlines=lines(readSmall(extraction_dir/"candidate_pairs.tsv",16*1024*1024));
    if(pairlines.empty()||pairlines[0]!="block\tread_id\tref_start\tref_end\tquery_start\tquery_end\tstatus")
        throw std::invalid_argument("unexpected candidate_pairs.tsv schema");
    std::map<std::pair<std::size_t,std::size_t>,bool> seen;
    for(std::size_t i=1;i<pairlines.size();++i){
        const auto c=tabs(pairlines[i]);
        if(c.size()!=7||c[6]!="CANDIDATE_ONLY")throw std::invalid_argument("invalid candidate pair");
        std::size_t read_index=source.size();
        for(std::size_t r=1;r<source.size();++r)if(source[r].id==c[1])read_index=r;
        const auto b=num(c[0]);
        const auto key=std::make_pair(b,read_index);
        if(read_index==source.size()||!candidates.contains(key)||!seen.emplace(key,true).second)
            throw std::invalid_argument("candidate association missing or duplicated");
        const auto* q=candidates.at(key);
        const auto* r=refs.at(b);
        if(num(c[2])!=r->start||num(c[3])!=r->end||num(c[4])!=q->start||num(c[5])!=q->end)
            throw std::invalid_argument("candidate coordinates do not match extracted tiles");
    }
    if(seen.size()!=candidates.size())throw std::invalid_argument("missing candidate pair row");
    const auto rejectlines=lines(readSmall(extraction_dir/"rejected_candidates.tsv",16*1024*1024));
    if(rejectlines.empty()||rejectlines[0]!="block\tread_id\treason")
        throw std::invalid_argument("unexpected rejected_candidates.tsv schema");
    std::set<std::pair<std::size_t,std::string>> rejects;
    for(std::size_t i=1;i<rejectlines.size();++i){
        const auto c=tabs(rejectlines[i]);
        if(c.size()!=3||c[2].empty()||!refs.contains(num(c[0]))||
           !rejects.insert({num(c[0]),c[1]}).second)
            throw std::invalid_argument("invalid rejected candidate row");
        std::size_t read_index=source.size();
        for(std::size_t r=1;r<source.size();++r)if(source[r].id==c[1])read_index=r;
        if(read_index==source.size()||candidates.contains({num(c[0]),read_index}))
            throw std::invalid_argument("rejected candidate was accepted or unknown read");
    }
    std::vector<std::tuple<std::size_t,std::string,std::string,std::size_t,std::size_t,std::size_t>> blockinfo;
    // blockinfo: block, status, reason, row count, total input bases, MSA columns.
    std::vector<std::tuple<std::size_t,const Tile*,std::string,std::string>> rowinfo;
    std::vector<std::pair<std::string,std::string>> fasta_rows;
    std::vector<std::tuple<std::size_t,std::string,std::string,std::string>> colops;
    Summary summary;
    summary.blocks=refs.size();summary.unresolved_tiles=unresolved;
    const native_msa::Options limits;
    for(const auto& [block,ref]:refs){
        std::vector<const Tile*> participating={ref};
        for(std::size_t r=1;r<source.size();++r){
            const auto it=candidates.find({block,r});
            if(it!=candidates.end())participating.push_back(it->second);
        }
        if(participating.size()==1){
            blockinfo.emplace_back(block,"SKIPPED","NO_QUERY_CANDIDATE",1,ref->bases.size(),0);
            ++summary.skipped_no_query;
            continue;
        }
        std::size_t total=0;bool too_big=participating.size()>limits.max_sequences;
        for(const auto* tile:participating){
            if(tile->bases.size()>limits.max_sequence_bases ||
               total+tile->bases.size()>limits.max_total_bases)too_big=true;
            total+=tile->bases.size(); // bounded by 4Mb from source FASTA
        }
        if(too_big){
            blockinfo.emplace_back(block,"SKIPPED","NATIVE_BUDGET",participating.size(),total,0);
            ++summary.skipped_budget;
            continue;
        }
        // PR3D.3B1: the original 25-kb / 100-kb length limits are NOT
        // an abPOA RAM budget. A modest number of long tiles can trigger
        // multi-gigabyte SIMD dynamic-programming matrix allocations.
        //
        // Upper-bound the pairwise length products across the progressively
        // inserted rows. This is a heuristic preflight (NOT an upper bound on
        // abPOA's actual bytes); refuse costly blocks rather than risk an
        // unchecked allocation that terminates the entire process.
        constexpr std::size_t kMaxEstimatedPairwiseCells = 4'000'000;
        std::size_t previous_bases = 0;
        std::size_t estimated_cells = 0;
        bool excessive_dp = false;
        for (const auto* tile : participating) {
            const auto length = tile->bases.size();
            if (previous_bases != 0 &&
                length > (kMaxEstimatedPairwiseCells - estimated_cells) / previous_bases) {
                excessive_dp = true;
                break;
            }
            estimated_cells += length * previous_bases; // division guard above
            previous_bases += length; // bounded by max_total_bases
        }
        if (excessive_dp) {
            blockinfo.emplace_back(block, "SKIPPED", "DP_CELL_BUDGET",
                                   participating.size(), total, 0);
            ++summary.skipped_budget;
            continue;
        }
        std::vector<native_msa::Sequence> input;
        for(const auto* tile:participating)input.push_back({tile->tile_id,tile->bases});
        const auto msa=native_msa::align(input,limits);
        if(msa.rows.size()!=input.size()||msa.columns==0||msa.columns>limits.max_total_bases)
            throw std::runtime_error("invalid abPOA local MSA dimensions");
        for(std::size_t i=0;i<input.size();++i){
            const auto& row=msa.rows[i];
            std::string recovered;
            for(const auto c:row.bases)if(c!='-')recovered+=c;
            if(row.id!=input[i].id||row.bases.size()!=msa.columns||recovered!=input[i].bases)
                throw std::runtime_error("abPOA changed tile bases, row order, or alignment width");
            const auto aligned_id="b"+std::to_string(block)+"|"+row.id;
            fasta_rows.emplace_back(aligned_id,row.bases);
            rowinfo.emplace_back(block,participating[i],aligned_id,
                                 lossless_extractor::sha256Digest(row.bases));
            ++summary.aligned_rows;
        }
        // Per-query operations are MSA-column oriented, NOT SAM CIGAR. P
        // explicitly records double-gap columns caused by other MSA rows.
        for(std::size_t r=1;r<msa.rows.size();++r){
            std::string ops;
            std::size_t ref_bases=0,query_bases=0;
            char prev='\0';std::size_t run=0;
            auto emit=[&](){if(run)ops+=std::to_string(run)+prev;};
            for(std::size_t c=0;c<msa.columns;++c){
                const bool a=msa.rows[0].bases[c]!='-';
                const bool b=msa.rows[r].bases[c]!='-';
                const char op=a?(b?'M':'D'):(b?'I':'P');
                ref_bases+=a;query_bases+=b;
                if(op!=prev){emit();prev=op;run=1;}else ++run;
            }
            emit();
            if(ref_bases!=participating[0]->bases.size()||
               query_bases!=participating[r]->bases.size())
                throw std::runtime_error("local MSA column operations lose bases");
            colops.emplace_back(block,participating[0]->tile_id,participating[r]->tile_id,ops);
        }
        ++summary.aligned;
        blockinfo.emplace_back(block,"TENTATIVE_LOCAL_MSA","CANDIDATE_ONLY",
                               participating.size(),total,msa.columns);
    }
    fs::path stage;
    for(unsigned i=0;i<64;++i){
        const auto ticks=std::chrono::steady_clock::now().time_since_epoch().count();
        stage=out.parent_path()/(".pr3d3b_stage_"+std::to_string(ticks)+"_"+std::to_string(i));
        if(fs::create_directory(stage))break;
        if(i==63)throw std::runtime_error("could not create staging directory");
    }
    try {
        {
            const auto p=stage/"local_msas.fasta";std::ofstream f(p,std::ios::binary);
            if(!f)throw std::runtime_error("cannot open local_msas.fasta");
            for(const auto& [id,row]:fasta_rows)writeFasta(f,id,row);
            flush(f,p);
        }
        {
            const auto p=stage/"blocks.tsv";std::ofstream f(p,std::ios::binary);
            if(!f)throw std::runtime_error("cannot open blocks.tsv");
            f<<"block\tstatus\treason\trows\tinput_bases\tcolumns\n";
            for(const auto& [block,status,reason,rows,bases,cols]:blockinfo)
                f<<block<<'\t'<<status<<'\t'<<reason<<'\t'<<rows<<'\t'<<bases<<'\t'<<cols<<'\n';
            flush(f,p);
        }
        {
            const auto p=stage/"rows.tsv";std::ofstream f(p,std::ios::binary);
            if(!f)throw std::runtime_error("cannot open rows.tsv");
            f<<"block\tread_id\ttile_id\tkind\tstart\tend\taligned_id\taligned_sha256\n";
            for(const auto& [block,tile,id,hash]:rowinfo)
                f<<block<<'\t'<<tile->read_id<<'\t'<<tile->tile_id<<'\t'<<tile->kind
                 <<'\t'<<tile->start<<'\t'<<tile->end<<'\t'<<id<<'\t'<<hash<<'\n';
            flush(f,p);
        }
        {
            const auto p=stage/"column_operations.tsv";std::ofstream f(p,std::ios::binary);
            if(!f)throw std::runtime_error("cannot open column_operations.tsv");
            f<<"block\tref_tile_id\tquery_tile_id\tmsa_column_ops\n";
            for(const auto& [block,ref,qry,ops]:colops)
                f<<block<<'\t'<<ref<<'\t'<<qry<<'\t'<<ops<<'\n';
            flush(f,p);
        }
        {
            const auto p=stage/"unaligned_tiles.tsv";std::ofstream f(p,std::ios::binary);
            if(!f)throw std::runtime_error("cannot open unaligned_tiles.tsv");
            f<<"tile_id\tread_id\tkind\tstart\tend\tblock\treason\n";
            for(const auto& tile:tiles){
                std::string reason;
                if(tile.kind=="UNRESOLVED")reason="UNRESOLVED";
                else if(std::get<1>(blockinfo[tile.block-1])!="TENTATIVE_LOCAL_MSA")reason="BLOCK_SKIPPED";
                if(!reason.empty())
                    f<<tile.tile_id<<'\t'<<tile.read_id<<'\t'<<tile.kind<<'\t'<<tile.start<<'\t'
                     <<tile.end<<'\t'<<tile.block<<'\t'<<reason<<'\n';
            }
            flush(f,p);
        }
        {
            const auto p=stage/"summary.json";std::ofstream f(p,std::ios::binary);
            if(!f)throw std::runtime_error("cannot open summary.json");
            f<<"{\n  \"schema\": \"pr3d3b-local-msas-v1\",\n"
             <<"  \"diagnostic_only\": true,\n  \"no_global_msa\": true,\n"
             <<"  \"candidate_homology_unverified\": true,\n"
             <<"  \"source_fasta_sha256\": \""<<source_hash<<"\",\n"
             <<"  \"source_extraction_summary_sha256\": \""
             <<lossless_extractor::sha256Digest(manifest)<<"\",\n"
             <<"  \"reference_id\": "<<jsonEsc(source[0].id)<<",\n"
             <<"  \"blocks\": "<<summary.blocks<<",\n"
             <<"  \"aligned_blocks\": "<<summary.aligned<<",\n"
             <<"  \"skipped_no_query\": "<<summary.skipped_no_query<<",\n"
             <<"  \"skipped_budget\": "<<summary.skipped_budget<<",\n"
             <<"  \"aligned_rows\": "<<summary.aligned_rows<<",\n"
             <<"  \"unresolved_tiles\": "<<summary.unresolved_tiles<<"\n}\n";
            flush(f,p);
        }
        publish(stage,out);
    }catch(...){std::error_code ec;fs::remove_all(stage,ec);throw;}
    return summary;
}
} // namespace block_local_msa
