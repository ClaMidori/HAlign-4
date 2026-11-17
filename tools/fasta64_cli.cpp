#include "ultralong/sparse_fasta.hpp"
#include <chrono>
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <cstdint>
#include <limits>

namespace fs=std::filesystem;
namespace {
std::uint64_t integer(const std::string& s) {
    if(s.empty()||s.find_first_not_of("0123456789")!=std::string::npos)
        throw std::invalid_argument("expected nonnegative integer: "+s);
    std::size_t consumed=0;
    const auto v=std::stoull(s,&consumed);
    if(consumed!=s.size())throw std::invalid_argument("bad integer: "+s);
    return v;
}
void link_new(const fs::path& staged,const fs::path& dst) {
    // Hard-link creation fails if destination already exists, including symlinks.
    fs::create_hard_link(staged,dst);
}
} // namespace
int main(int argc,char** argv) {
    fs::path source,report,fetch_out;
    std::string fetch;
    std::uint64_t stride=65536;
    try {
        for(int i=1;i<argc;++i) {
            const std::string arg=argv[i];
            if(arg=="--help"||arg=="-h"){
                std::cout<<"halign4_fasta64_index --input FASTA --output SUMMARY.tsv "
                          "[--stride 65536] [--fetch read_index:start:length --slice-output SLICE.fasta]\n"
                          "INDEX_ONLY, NOT AN MSA. Supports uncompressed FASTA with A/C/G/T/N.\n";
                return 0;
            }
            if(++i==argc)throw std::invalid_argument("missing argument value for "+arg);
            const std::string value=argv[i];
            if(arg=="--input")source=value;
            else if(arg=="--output")report=value;
            else if(arg=="--stride")stride=integer(value);
            else if(arg=="--fetch")fetch=value;
            else if(arg=="--slice-output")fetch_out=value;
            else throw std::invalid_argument("unknown option "+arg);
        }
        if(source.empty()||report.empty()||(!fetch.empty() && fetch_out.empty())||
           (fetch.empty() && !fetch_out.empty()))
            throw std::invalid_argument("missing required args or incomplete fetch pair");
        const auto src=fs::absolute(source).lexically_normal();
        const auto dst=fs::absolute(report).lexically_normal();
        const auto fdst=fetch.empty()?fs::path{}:fs::absolute(fetch_out).lexically_normal();
        if(src==dst||fs::exists(dst)||fs::is_symlink(dst)||
           !fs::exists(dst.parent_path()) || !fs::is_directory(dst.parent_path()))
            throw std::invalid_argument("output must be a new file in existing directory");
        if(!fetch.empty() && (src==fdst||dst==fdst||fs::exists(fdst)||fs::is_symlink(fdst)||
          dst.parent_path()!=fdst.parent_path()))
            throw std::invalid_argument("slice output must be distinct and new in report directory");
        const auto start=std::chrono::steady_clock::now();
        const auto index=ultralong::SparseFasta::build(src,stride);
        std::string slice;
        std::string slice_name;
        std::uint64_t slice_start=0;
        if(!fetch.empty()) {
            const auto p=fetch.find(':'),q=fetch.find(':',p==std::string::npos?p:p+1);
            if(p==std::string::npos||q==std::string::npos||fetch.find(':',q+1)!=std::string::npos)
                throw std::invalid_argument("fetch needs rowIndex:start:length");
            const auto row=integer(fetch.substr(0,p));
            slice_start=integer(fetch.substr(p+1,q-p-1));
            const auto length=integer(fetch.substr(q+1));
            if(row>=index.reads().size())throw std::out_of_range("fetch record outside index");
            slice_name=index.reads()[static_cast<std::size_t>(row)].id;
            slice=index.slice(static_cast<std::size_t>(row),slice_start,length);
        }
        const auto now=std::chrono::steady_clock::now().time_since_epoch().count();
        fs::path tmpdir;
        for(int k=0;k<128;++k) {
            const auto trial=dst.parent_path()/(".halign4_index_"+std::to_string(now)+"_"+std::to_string(k));
            std::error_code ec;
            if(fs::create_directory(trial,ec)){tmpdir=trial;break;}
            if(ec)throw std::runtime_error("cannot create staging directory: "+ec.message());
        }
        if(tmpdir.empty())throw std::runtime_error("cannot allocate staging directory");
        try {
            const auto intermediate=tmpdir/"summary.tsv";
            std::ofstream out(intermediate,std::ios::binary);
            if(!out)throw std::runtime_error("cannot create staged index report");
            out<<"read_index\tread_id\tbases\tcheckpoints\tstride\n";
            for(std::size_t r=0;r<index.reads().size();++r) {
                const auto& entry=index.reads()[r];
                out<<r<<'\t'<<entry.id<<'\t'<<entry.bases<<'\t'<<entry.checkpoints.size()
                   <<'\t'<<index.checkpointStride()<<'\n';
            }
            out.flush();out.close();
            if(!out)throw std::runtime_error("index summary write failure");
            if(!fetch.empty()) {
                std::ofstream x(tmpdir/"slice.fasta",std::ios::binary);
                if(!x)throw std::runtime_error("cannot create staged slice");
                x<<'>'<<slice_name<<"|start="<<slice_start<<"|length="<<slice.size()<<'\n';
                for(std::size_t i=0;i<slice.size();i+=80)
                    x.write(slice.data()+i,static_cast<std::streamsize>(std::min<std::size_t>(80,slice.size()-i)))<<'\n';
                x.flush();x.close();
                if(!x)throw std::runtime_error("slice output write failure");
            }
            link_new(intermediate,dst);
            try{if(!fetch.empty())link_new(tmpdir/"slice.fasta",fdst);}
            catch(...){std::error_code ec;fs::remove(dst,ec);throw;}
        }catch(...){std::error_code ec;fs::remove_all(tmpdir,ec);throw;}
        std::error_code ec;fs::remove_all(tmpdir,ec);
        const double elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        std::cout<<"INDEX64=PASS FASTA_READS="<<index.reads().size()
                 <<" DNA_BASES="<<index.totalBases()<<" STRIDE="<<index.checkpointStride()
                 <<" WALL_SEC="<<elapsed<<"\n"
                 <<"INFERENCE=NOT_RUN MSA=NOT_PRODUCED OUT_OF_CORE_FASTA_ACCESS=PASS\n";
        return 0;
    }catch(const std::exception& e){std::cerr<<"INDEX64_FAIL_CLOSED: "<<e.what()<<'\n';return 2;}
}
