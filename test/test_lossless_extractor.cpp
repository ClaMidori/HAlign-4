#include "doctest/doctest.h"
#include "align/lossless_extractor.hpp"

#include <chrono>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>

namespace {
namespace fs=std::filesystem;
struct Fixture {
    fs::path dir, fasta, plan, dest;
    Fixture() {
        dir=fs::temp_directory_path()/(
            "halign4_pr3d3a_unit_"+
            std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
        fs::create_directory(dir);
        fasta=dir/"input.fasta";plan=dir/"plan.tsv";dest=dir/"extracted";
        std::ofstream(fasta)<<">reference\nACGTACGT\n>query1\nACGTNNACGT\n";
    }
    ~Fixture(){std::error_code ec;fs::remove_all(dir,ec);}
    void writePlan(bool overlap=false, bool invalid=false) {
        std::ofstream f(plan);
        f<<"# PR3D.1 diagnostic only; NO alignment and NO verified homologous coordinates\n"
         <<"# reference=reference\tlength="<<(invalid?9:8)<<"\tk=11\ttarget_bp=20000\n"
         <<"block\tread_id\tref_start\tref_end\tquery_start_est\tquery_end_est\tleft_anchor_distance\tright_anchor_distance\tsupported_seed_count\tstatus\n"
         <<"1\treference\t0\t4\t0\t4\t0\t0\t0\tREFERENCE_ONLY\n"
         <<"1\tquery1\t0\t4\t0\t"<<(overlap?6:4)<<"\t2\t2\t9\tCANDIDATE_ONLY\n"
         <<"2\treference\t4\t8\t4\t8\t0\t0\t0\tREFERENCE_ONLY\n"
         <<"2\tquery1\t4\t8\t"<<(overlap?4:6)<<"\t10\t2\t2\t9\tCANDIDATE_ONLY\n";
    }
};
std::string read(const fs::path& file) {
    std::ifstream f(file,std::ios::binary);
    return {std::istreambuf_iterator<char>(f),std::istreambuf_iterator<char>()};
}
} // namespace
TEST_SUITE("lossless_extractor") {
    TEST_CASE("accepted blocks plus unresolved core reconstruct without loss") {
        Fixture x;x.writePlan();
        const auto s=lossless_extractor::extract(x.fasta,x.plan,x.dest);
        CHECK(s.reads==2);CHECK(s.blocks==2);
        CHECK(s.accepted_query_candidates==2);
        CHECK(s.rejected_conflicting_candidates==0);
        CHECK(s.tiles==5);CHECK(s.unresolved_query_bases==2);
        CHECK(read(x.dest/"tiles.tsv").find("UNRESOLVED\t4\t6\t\t2")!=std::string::npos);
        CHECK(read(x.dest/"tiles.fasta").find(">t02_00002\nNN\n")!=std::string::npos);
        CHECK(read(x.dest/"summary.json").find("\"no_msa\": true")!=std::string::npos);
        CHECK(read(x.dest/"candidate_pairs.tsv").find("CANDIDATE_ONLY")!=std::string::npos);
        CHECK(read(x.dest/"rejected_candidates.tsv")=="block\tread_id\treason\n");
    }
    TEST_CASE("overlap rejects BOTH candidates without deleting sequence bases") {
        Fixture x;x.writePlan(true);
        const auto s=lossless_extractor::extract(x.fasta,x.plan,x.dest);
        CHECK(s.accepted_query_candidates==0);
        CHECK(s.rejected_conflicting_candidates==2);
        CHECK(s.unresolved_query_bases==10);
        CHECK(s.tiles==3);
        CHECK(read(x.dest/"rejected_candidates.tsv").find("QUERY_OVERLAP")!=std::string::npos);
        CHECK(read(x.dest/"tiles.fasta").find(">t02_00001\nACGTNNACGT\n")!=std::string::npos);
    }
    TEST_CASE("bad metadata or pre-existing output is refused without changes") {
        Fixture x;x.writePlan(false,true);
        CHECK_THROWS_AS(lossless_extractor::extract(x.fasta,x.plan,x.dest),std::invalid_argument);
        CHECK_FALSE(fs::exists(x.dest));
        x.writePlan();fs::create_directory(x.dest);
        CHECK_THROWS_AS(lossless_extractor::extract(x.fasta,x.plan,x.dest),std::invalid_argument);
        CHECK(fs::is_empty(x.dest));
    }
}
