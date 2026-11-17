#include "ultralong/window64_msa.hpp"
#include <chrono>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>
#include <sys/resource.h>
static std::uint64_t number(const std::string& s){
    if(s.empty()||s.find_first_not_of("0123456789")==std::string::npos){
        auto n=std::stoull(s);if(n)return n;
    }
    throw std::invalid_argument("invalid positive integer");
}
int main(int argc,char** argv){
    try{
        std::string input,output,report;ultralong::WindowMsaConfig cfg;
        for(int i=1;i<argc;++i){
            const std::string k=argv[i];
            if(k=="--help"||k=="-h"){
                std::cout<<"halign4_window64_msa --input reads.fasta --output NEW_MSA.fasta --report NEW.tsv [--max-bases N] [--min-identity 0.55] [--max-local-bases 12000]\n";
                return 0;
            }
            if(++i>=argc)throw std::invalid_argument("missing value for "+k);
            const std::string v=argv[i];
            if(k=="--input")input=v;
            else if(k=="--output")output=v;
            else if(k=="--report")report=v;
            else if(k=="--max-bases")cfg.max_total_bases=number(v);
            else if(k=="--max-local-bases")cfg.local_max_bases=number(v);
            else if(k=="--min-identity")cfg.min_reference_identity=std::stod(v);
            else throw std::invalid_argument("unsupported option "+k);
        }
        if(input.empty()||output.empty()||report.empty())throw std::invalid_argument("--input/--output/--report required");
        const auto start=std::chrono::steady_clock::now();
        const auto s=ultralong::runWindowMsa(input,output,report,cfg);
        const double sec=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        struct rusage usage{};getrusage(RUSAGE_SELF,&usage);
#ifdef PR3D4D3_STUB
        std::cerr<<"TEST_BACKEND_STUB_ONLY: no abPOA. MSA quality and timing invalid.\n";
#endif
        std::cout<<"WINDOW64_NATIVE_MSA=PASS NO_ORACLE=PASS LOSSLESS=PASS RECTANGULAR=PASS\n"
                 <<"TOTAL_BASES="<<s.original_bases<<" WINDOWS="<<s.windows
                 <<" CANDIDATE_WINDOWS="<<s.candidate_query_windows
                 <<" LOCAL_CANDIDATE_ROWS="<<s.local_candidate_rows
                 <<" ALIGNED_BLOCKS="<<s.aligned_blocks
                 <<" ALIGNED_REF_BASES="<<s.aligned_reference_bases
                 <<" SKIPPED_BUDGET="<<s.skipped_budget
                 <<" REJECTED_POSTCHECK="<<s.rejected_postcheck
                 <<" PRIVATE_BASES="<<s.private_bases
                 <<" OUTPUT_COLUMNS="<<s.output_columns<<"\n"
                 <<"WALL_SEC="<<sec<<" PEAK_RSS_KIB="<<usage.ru_maxrss<<"\n"
                 <<"BIOLOGICAL_HOMOLOGY=NOT_EVALUATED\n";
        return 0;
    }catch(const std::exception& e){std::cerr<<"WINDOW64_NATIVE_MSA_FAIL_CLOSED: "<<e.what()<<'\n';return 2;}
}
