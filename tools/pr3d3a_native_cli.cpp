#include "align/lossless_extractor.hpp"
#include <iostream>
#include <stdexcept>
#include <string>
int main(int argc,char** argv) {
    try {
        std::string fasta,plan,out;
        for(int i=1;i<argc;++i){
            const std::string s=argv[i];
            if((s=="--fasta"||s=="--plan"||s=="--out-dir") && i+1<argc){
                auto& target=s=="--fasta"?fasta:s=="--plan"?plan:out;
                if(!target.empty())throw std::invalid_argument("duplicate option "+s);
                target=argv[++i];
            }else throw std::invalid_argument("usage: halign4_pr3d3a_extract --fasta INPUT --plan TSV --out-dir NEW_DIRECTORY");
        }
        if(fasta.empty()||plan.empty()||out.empty())throw std::invalid_argument("all arguments are required");
        const auto v=lossless_extractor::extract(fasta,plan,out);
        std::cout<<"PR3D3A PASS: reads="<<v.reads<<" blocks="<<v.blocks
                 <<" query_candidates="<<v.accepted_query_candidates
                 <<" conflict_rejections="<<v.rejected_conflicting_candidates
                 <<" tiles="<<v.tiles<<" unresolved_query_bases="<<v.unresolved_query_bases
                 <<" reconstructed="<<v.reads<<'\n';
        return 0;
    }catch(const std::exception& e){std::cerr<<"PR3D3A REFUSED: "<<e.what()<<'\n';return 2;}
}
