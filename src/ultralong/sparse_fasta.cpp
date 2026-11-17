#include "ultralong/sparse_fasta.hpp"
#include <algorithm>
#include <array>
#include <cctype>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <unordered_set>

namespace fs = std::filesystem;
namespace ultralong {
namespace {
std::uint64_t fileLength(const fs::path& path) {
    const auto size=fs::file_size(path);
    if (size>std::numeric_limits<std::uint64_t>::max())
        throw std::runtime_error("FASTA source file too large for 64-bit indexing");
    return static_cast<std::uint64_t>(size);
}
char base(char ch) {
    const char upper=static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
    if (upper != 'A' && upper != 'C' && upper != 'G' && upper != 'T' && upper != 'N')
        throw std::invalid_argument("unsupported FASTA base: only A/C/G/T/N allowed; no gaps");
    return upper;
}
} // namespace
SparseFasta SparseFasta::build(const fs::path& path, std::uint64_t stride, std::size_t max_reads) {
    if (stride < 1024 || stride > (1ULL<<24) || max_reads < 1 || max_reads > 100000)
        throw std::invalid_argument("invalid FASTA checkpoint stride/max_reads");
    if (!fs::is_regular_file(path)) throw std::invalid_argument("FASTA must be a regular uncompressed file");
    SparseFasta result;
    result.path_=fs::absolute(path).lexically_normal();
    result.stride_=stride;
    result.source_bytes_=fileLength(result.path_);
    result.source_write_time_=fs::last_write_time(result.path_);
    std::ifstream in(result.path_,std::ios::binary);
    if(!in)throw std::runtime_error("cannot open FASTA for indexing");
    std::array<char,65536> buffer{};
    std::uint64_t offset=0;
    bool line_start=true, in_header=false, id_done=false;
    std::size_t header_length=0;
    std::unordered_set<std::string> ids;
    while(in) {
        in.read(buffer.data(),static_cast<std::streamsize>(buffer.size()));
        const auto got=in.gcount();
        for(std::streamsize p=0;p<got;++p) {
            const char ch=buffer[static_cast<std::size_t>(p)];
            if(ch=='\n') {
                if(in_header) {
                    if(result.reads_.back().id.empty() || !ids.insert(result.reads_.back().id).second)
                        throw std::invalid_argument("empty or duplicated FASTA identifier");
                }
                line_start=true;in_header=false;
            } else if(ch=='\r') {
                // CRLF is tolerated; the original parser also ignores CR.
            } else if(line_start) {
                line_start=false;
                if(ch=='>') {
                    if(!result.reads_.empty() && result.reads_.back().bases==0)
                        throw std::invalid_argument("empty FASTA read");
                    if(result.reads_.size()==max_reads)throw std::invalid_argument("too many FASTA records");
                    result.reads_.push_back(ReadIndex{});
                    in_header=true;id_done=false;header_length=0;
                } else {
                    if(result.reads_.empty())throw std::invalid_argument("sequence before FASTA header");
                    auto& r=result.reads_.back();
                    if(r.bases%stride==0)r.checkpoints.push_back({r.bases,offset});
                    (void)base(ch);
                    if(r.bases==UINT64_MAX||result.total_bases_==UINT64_MAX)
                        throw std::overflow_error("FASTA nucleotide coordinate overflow");
                    ++r.bases;++result.total_bases_;
                }
            } else if(in_header) {
                if(++header_length>256)throw std::invalid_argument("FASTA header longer than 256 bytes");
                if(ch==' '||ch=='\t')id_done=true;
                else if(!id_done){
                    if(static_cast<unsigned char>(ch)<33 || ch=='>')
                        throw std::invalid_argument("invalid FASTA read name");
                    result.reads_.back().id.push_back(ch);
                }
            } else {
                auto& r=result.reads_.back();
                if(r.bases%stride==0)r.checkpoints.push_back({r.bases,offset});
                (void)base(ch);
                if(r.bases==UINT64_MAX||result.total_bases_==UINT64_MAX)
                    throw std::overflow_error("FASTA nucleotide coordinate overflow");
                ++r.bases;++result.total_bases_;
            }
            if(offset==UINT64_MAX)throw std::overflow_error("FASTA byte offset overflow");
            ++offset;
        }
    }
    if(in.bad() || offset!=result.source_bytes_)throw std::runtime_error("FASTA read failed or source changed while indexing");
    if(in_header && (result.reads_.back().id.empty() || !ids.insert(result.reads_.back().id).second))
        throw std::invalid_argument("empty or duplicated FASTA identifier");
    if(result.reads_.empty() || result.reads_.back().bases==0)
        throw std::invalid_argument("missing records or final empty FASTA read");
    for(const auto& r:result.reads_)if(r.checkpoints.empty()||r.checkpoints.front().nucleotide!=0)
        throw std::runtime_error("index checkpoint invariant failed");
    result.verifySource();
    return result;
}
void SparseFasta::verifySource() const {
    if(!fs::is_regular_file(path_) || fileLength(path_)!=source_bytes_ ||
       fs::last_write_time(path_)!=source_write_time_)
        throw std::runtime_error("FASTA source changed after indexing; rebuild the index");
}
std::string SparseFasta::slice(std::size_t record,std::uint64_t start,
                                std::uint64_t count,std::uint64_t max_bytes) const {
    if(record>=reads_.size())throw std::out_of_range("unknown FASTA record index");
    const auto& r=reads_[record];
    if(start>r.bases || count>r.bases-start)throw std::out_of_range("FASTA slice out of range");
    if(count>max_bytes || max_bytes > (1ULL<<26))
        throw std::invalid_argument("FASTA slice exceeds bounded 64 MiB allocation policy");
    verifySource();
    if(!count)return {};
    const auto it=std::upper_bound(r.checkpoints.begin(),r.checkpoints.end(),start,
                       [](std::uint64_t pos,const Checkpoint& cp){return pos<cp.nucleotide;});
    if(it==r.checkpoints.begin())throw std::logic_error("missing FASTA checkpoint");
    const auto& cp=*(it-1);
    if(cp.byte_offset>static_cast<std::uint64_t>(std::numeric_limits<std::streamoff>::max()))
        throw std::overflow_error("FASTA file offset exceeds seekg limit");
    std::ifstream in(path_,std::ios::binary);
    if(!in)throw std::runtime_error("cannot open FASTA for slice");
    in.seekg(static_cast<std::streamoff>(cp.byte_offset));
    if(!in)throw std::runtime_error("FASTA seek failed");
    std::string result;
    result.reserve(static_cast<std::size_t>(count));
    std::uint64_t cursor=cp.nucleotide;
    // A checkpoint always points to a DNA byte. Newlines and CR are skipped.
    // A header before completion indicates source mutation or malformed index.
    char ch=0;
    bool at_line_start=false;
    while(cursor<start+count && in.get(ch)) {
        if(ch=='\n'){at_line_start=true;continue;}
        if(ch=='\r')continue;
        if(at_line_start && ch=='>')throw std::runtime_error("FASTA record ended unexpectedly during slice");
        at_line_start=false;
        const char b=base(ch);
        if(cursor>=start)result.push_back(b);
        ++cursor;
    }
    if(!in && !in.eof())throw std::runtime_error("I/O error while slicing FASTA");
    if(result.size()!=count)throw std::runtime_error("FASTA slice unexpectedly truncated");
    verifySource();
    return result;
}
} // namespace ultralong
