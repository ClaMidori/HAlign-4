#include "doctest/doctest.h"
#include "align/block_local_msa.hpp"
#include "align/lossless_extractor.hpp"

#include <chrono>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>

namespace {
namespace fs=std::filesystem;
struct Fixture {
    fs::path dir,fasta,plan,extract_dir,aligned_dir;
    Fixture(){
        dir=fs::temp_directory_path()/("halign4_pr3d3b_"+
            std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
        fs::create_directories(dir);
        fasta=dir/"input.fasta";plan=dir/"plan.tsv";
        extract_dir=dir/"extracted";aligned_dir=dir/"local";
        std::ofstream(fasta)<<">ref\nACGTACGT\n>q1\nACGTNNACGT\n>q2\nACGTACGT\n";
    }
    ~Fixture(){std::error_code ec;fs::remove_all(dir,ec);}
    void writePlan(bool no_candidates=false) {
        std::ofstream f(plan);
        f<<"# PR3D.1 diagnostic only; NO alignment and NO verified homologous coordinates\n"
         <<"# reference=ref\tlength=8\tk=11\ttarget_bp=20000\n"
         <<"block\tread_id\tref_start\tref_end\tquery_start_est\tquery_end_est\tleft_anchor_distance\tright_anchor_distance\tsupported_seed_count\tstatus\n";
        for(int b=0;b<2;++b){
            const int a=b*4,z=a+4;
            f<<b+1<<"\tref\t"<<a<<'\t'<<z<<'\t'<<a<<'\t'<<z<<"\t0\t0\t0\tREFERENCE_ONLY\n";
            if(no_candidates){
                f<<b+1<<"\tq1\t"<<a<<'\t'<<z<<"\tNA\tNA\tNA\tNA\t0\tUNRESOLVED\n";
                f<<b+1<<"\tq2\t"<<a<<'\t'<<z<<"\tNA\tNA\tNA\tNA\t0\tUNRESOLVED\n";
            } else {
                f<<b+1<<"\tq1\t"<<a<<'\t'<<z<<'\t'<<(b?6:0)<<'\t'<<(b?10:4)<<"\t1\t1\t9\tCANDIDATE_ONLY\n";
                f<<b+1<<"\tq2\t"<<a<<'\t'<<z<<'\t'<<a<<'\t'<<z<<"\t1\t1\t9\tCANDIDATE_ONLY\n";
            }
        }
    }
    void extract(bool no_candidates=false){writePlan(no_candidates);lossless_extractor::extract(fasta,plan,extract_dir);}
};
std::string read(const fs::path& path){
    std::ifstream f(path,std::ios::binary);
    return {std::istreambuf_iterator<char>(f),std::istreambuf_iterator<char>()};
}
}
TEST_SUITE("block_local_msa") {
    TEST_CASE("local abPOA MSAs preserve tile order, bases, and unresolved") {
        Fixture x;x.extract();
        const auto r=block_local_msa::align(x.fasta,x.extract_dir,x.aligned_dir);
        CHECK(r.blocks==2);CHECK(r.aligned==2);CHECK(r.aligned_rows==6);
        CHECK(r.skipped_no_query==0);CHECK(r.skipped_budget==0);
        CHECK(r.unresolved_tiles==1);
        const auto local=read(x.aligned_dir/"local_msas.fasta");
        CHECK(local.find(">b1|t01_00001\nACGT\n")!=std::string::npos);
        CHECK(local.find(">b2|t02_00003\nACGT\n")!=std::string::npos);
        CHECK(local.find("NN")==std::string::npos);
        CHECK(read(x.aligned_dir/"unaligned_tiles.tsv").find("UNRESOLVED")!=std::string::npos);
        CHECK(read(x.aligned_dir/"column_operations.tsv").find("4M")!=std::string::npos);
        CHECK(read(x.aligned_dir/"summary.json").find("\"no_global_msa\": true")!=std::string::npos);
    }
    TEST_CASE("no query candidates means no local alignment, never forcing unresolved") {
        Fixture x;x.extract(true);
        const auto r=block_local_msa::align(x.fasta,x.extract_dir,x.aligned_dir);
        CHECK(r.aligned==0);CHECK(r.skipped_no_query==2);
        CHECK(read(x.aligned_dir/"local_msas.fasta").empty());
        CHECK(read(x.aligned_dir/"blocks.tsv").find("NO_QUERY_CANDIDATE")!=std::string::npos);
    }
    TEST_CASE("source mismatch, forged candidate, and existing output fail closed") {
        Fixture x;x.extract();
        fs::create_directory(x.aligned_dir);
        CHECK_THROWS_AS(block_local_msa::align(x.fasta,x.extract_dir,x.aligned_dir),std::invalid_argument);
        fs::remove(x.aligned_dir);
        std::ofstream(x.fasta,std::ios::app)<<"N\n";
        CHECK_THROWS(block_local_msa::align(x.fasta,x.extract_dir,x.aligned_dir));
        CHECK_FALSE(fs::exists(x.aligned_dir));
        // Restore source, damage the cross-file pair coordinate.
        std::ofstream(x.fasta)<<">ref\nACGTACGT\n>q1\nACGTNNACGT\n>q2\nACGTACGT\n";
        auto pairs=read(x.extract_dir/"candidate_pairs.tsv");
        auto pos=pairs.find("\t0\t4\t0\t4\tCANDIDATE_ONLY");
        REQUIRE(pos!=std::string::npos);
        pairs.replace(pos,2,"\t1");
        std::ofstream(x.extract_dir/"candidate_pairs.tsv")<<pairs;
        CHECK_THROWS_AS(block_local_msa::align(x.fasta,x.extract_dir,x.aligned_dir),std::invalid_argument);
        CHECK_FALSE(fs::exists(x.aligned_dir));
    }
    TEST_CASE("nominally allowed 4-kb blocks are skipped when quadratic DP cost is too high") {
        Fixture x;
        const std::string dna(4096, 'A');
        std::ofstream(x.fasta) << ">ref\n" << dna << "\n>q1\n" << dna << "\n";
        {
            std::ofstream f(x.plan);
            f << "# PR3D.1 diagnostic only; NO alignment and NO verified homologous coordinates\n"
              << "# reference=ref\tlength=4096\tk=11\ttarget_bp=20000\n"
              << "block\tread_id\tref_start\tref_end\tquery_start_est\tquery_end_est\tleft_anchor_distance\tright_anchor_distance\tsupported_seed_count\tstatus\n"
              << "1\tref\t0\t4096\t0\t4096\t0\t0\t0\tREFERENCE_ONLY\n"
              << "1\tq1\t0\t4096\t0\t4096\t1\t1\t9\tCANDIDATE_ONLY\n";
        }
        lossless_extractor::extract(x.fasta, x.plan, x.extract_dir);
        const auto result = block_local_msa::align(x.fasta, x.extract_dir, x.aligned_dir);
        CHECK(result.blocks==1);
        CHECK(result.aligned==0);
        CHECK(result.skipped_budget==1);
        CHECK(read(x.aligned_dir/"blocks.tsv").find("DP_CELL_BUDGET") != std::string::npos);
        CHECK(read(x.aligned_dir/"local_msas.fasta").empty());
        CHECK(read(x.aligned_dir/"unaligned_tiles.tsv").find("BLOCK_SKIPPED") != std::string::npos);
    }
    TEST_CASE("over-budget projected region is deferred, not aligned") {
        Fixture x;
        std::ofstream(x.fasta)<<">ref\n"<<std::string(30001,'A')
            <<"\n>q1\n"<<std::string(30001,'A')<<"\n";
        std::ofstream f(x.plan);
        f<<"# PR3D.1 diagnostic only; NO alignment and NO verified homologous coordinates\n"
         <<"# reference=ref\tlength=30001\tk=11\ttarget_bp=20000\n"
         <<"block\tread_id\tref_start\tref_end\tquery_start_est\tquery_end_est\tleft_anchor_distance\tright_anchor_distance\tsupported_seed_count\tstatus\n"
         <<"1\tref\t0\t30001\t0\t30001\t0\t0\t0\tREFERENCE_ONLY\n"
         <<"1\tq1\t0\t30001\t0\t30001\t1\t1\t9\tCANDIDATE_ONLY\n";
        f.close();
        lossless_extractor::extract(x.fasta,x.plan,x.extract_dir);
        const auto s=block_local_msa::align(x.fasta,x.extract_dir,x.aligned_dir);
        CHECK(s.blocks==1);CHECK(s.aligned==0);CHECK(s.skipped_budget==1);
        CHECK(read(x.aligned_dir/"local_msas.fasta").empty());
        CHECK(read(x.aligned_dir/"blocks.tsv").find("NATIVE_BUDGET")!=std::string::npos);
    }
}
