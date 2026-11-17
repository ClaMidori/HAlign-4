#include "ultralong/window64_plan.hpp"
#include <chrono>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>
#include <sys/resource.h>
static std::uint64_t num(const std::string& s){
    if(s.empty()||s.find_first_not_of("0123456789")!=std::string::npos)
        throw std::invalid_argument("positive integer required");
    return std::stoull(s);
}
int main(int argc,char** argv){
    try {
        std::string input,output; ultralong::Window64Config cfg;
        for(int i=1;i<argc;++i){
            const std::string opt=argv[i];
            if(opt=="--help"){
                std::cout<<"halign4_window64_plan --input input.fa --output NEW.tsv [--window 65536] [--flank 8192] [--max-bases 800000000]\n"
                         <<"This produces only 64-bit window homology CANDIDATES. No MSA and no oracle.\n";return 0;
            }
            if(i+1>=argc)throw std::invalid_argument("missing value for "+opt);
            const std::string value=argv[++i];
            if(opt=="--input")input=value;
            else if(opt=="--output")output=value;
            else if(opt=="--window")cfg.window_bases=num(value);
            else if(opt=="--flank")cfg.flank_bases=num(value);
            else if(opt=="--max-bases")cfg.max_total_bases=num(value);
            else throw std::invalid_argument("unsupported option: "+opt);
        }
        if(input.empty()||output.empty())throw std::invalid_argument("--input and --output required");
        auto start=std::chrono::steady_clock::now();
        const auto s=ultralong::planWindows(input,output,cfg);
        const auto sec=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        struct rusage usage{};getrusage(RUSAGE_SELF,&usage);
        std::cout<<"WINDOW64_PLAN=PASS NO_MSA=YES NO_ORACLE=YES NOT_VERIFIED_HOMOLOGY=YES\n"
                 <<"DNA_BASES="<<s.total_bases<<" WINDOWS="<<s.windows
                 <<" CANDIDATE_QUERY_WINDOWS="<<s.candidate_pairs
                 <<" UNRESOLVED_QUERY_WINDOWS="<<s.unresolved_pairs
                 <<" ANCHORED_REFERENCE_BP_SUM="<<s.anchored_ref_bases<<'\n'
                 <<"WALL_SEC="<<sec<<" PEAK_RSS_KIB="<<usage.ru_maxrss<<'\n';
        return 0;
    }catch(const std::exception& e){std::cerr<<"WINDOW64_FAIL_CLOSED: "<<e.what()<<'\n';return 2;}
}
