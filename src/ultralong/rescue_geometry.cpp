#include "ultralong/rescue_geometry.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <limits>

namespace ultralong {
namespace {
std::optional<std::int64_t> projected(const std::vector<homology_plan::Anchor>& a,
                                      std::uint32_t r, std::uint32_t radius) {
    const auto it=std::lower_bound(a.begin(),a.end(),r,
        [](const homology_plan::Anchor& v,std::uint32_t pos){return v.reference<pos;});
    const homology_plan::Anchor* best=nullptr;
    if(it!=a.end())best=&*it;
    if(it!=a.begin()){
        const auto& prev=*std::prev(it);
        if(!best || static_cast<std::uint64_t>(r-prev.reference)<=
                         static_cast<std::uint64_t>(std::llabs(std::int64_t(best->reference)-r))) best=&prev;
    }
    if(!best || std::llabs(std::int64_t(best->reference)-r)>radius)return std::nullopt;
    // Both sides of a nearby indel cannot imply distant, conflicting coordinates.
    if(it!=a.begin() && it!=a.end()){
        const auto& before=*std::prev(it);const auto& after=*it;
        if(std::llabs(std::int64_t(before.reference)-r)<=radius &&
           std::llabs(std::int64_t(after.reference)-r)<=radius &&
           std::llabs((std::int64_t(after.query)-after.reference) -
                      (std::int64_t(before.query)-before.reference))>32) return std::nullopt;
    }
    return static_cast<std::int64_t>(best->query) + std::int64_t(r)-best->reference;
}
}
std::optional<RescueRange> projectRescueRange(
    const std::vector<homology_plan::Anchor>& a,
    const std::string& ref,const std::string& query,
    std::uint32_t start,std::uint32_t end,double min_identity) {
    if(end<=start || end>ref.size() || start>ref.size() || query.size()>100000 ||
       ref.size()>100000 || a.size()<8 || min_identity<0.5 || min_identity>1.0) return std::nullopt;
    const auto span=static_cast<std::uint32_t>(end-start);
    if(span<2048)return std::nullopt;
    // Enforce the trusted monotone unique-coordinate contract at the API edge.
    for(std::size_t i=1;i<a.size();++i){
        if(a[i-1].reference>=a[i].reference || a[i-1].query>=a[i].query)
            return std::nullopt;
    }
    std::array<std::uint32_t,4> evidence{};
    std::size_t support=0;
    std::int64_t min_delta=INT64_MAX,max_delta=INT64_MIN;
    std::uint32_t last=0;
    bool have_last=false;
    for(const auto& v:a){
        if(v.reference<start || v.reference>=end || v.query>=query.size())continue;
        if(have_last && v.reference-last<9)continue;
        evidence[std::min<unsigned>(3,static_cast<unsigned>(4ull*(v.reference-start)/span))]++;
        last=v.reference;have_last=true;++support;
        const auto d=static_cast<std::int64_t>(v.query)-v.reference;
        min_delta=std::min(min_delta,d);
        max_delta=std::max(max_delta,d);
    }
    if(support<8 || max_delta-min_delta>32 ||
       std::any_of(evidence.begin(),evidence.end(),[](auto x){return x<1;}))return std::nullopt;
    const auto radius=std::min<std::uint32_t>(1000,span/4);
    const auto x=projected(a,start,radius),y=projected(a,end,radius);
    if(!x||!y||*x<0||*y<=*x||static_cast<std::uint64_t>(*y)>query.size())return std::nullopt;
    // Compensated large indels, repeats and abrupt offsets are NOT globally rescued.
    if(std::llabs((*y-*x)-span)>32)return std::nullopt;
    // Two physically separated context samples. Prevent identity decided by a small motif.
    std::uint64_t matches=0,total=0;
    const std::array<std::uint32_t,4> offsets{span/8,3*span/8,5*span/8,7*span/8};
    for(auto rdelta:offsets){
        const auto rs=start+rdelta;
        const auto q_est=projected(a,rs,radius);
        if(!q_est||*q_est<0)return std::nullopt;
        std::uint64_t match=0,count=0;
        for(std::int64_t d=-128;d<128;++d){
            const auto rp=std::int64_t(rs)+d, qp=*q_est+d;
            if(rp<static_cast<std::int64_t>(start)||rp>=end||qp<*x||qp>=*y)continue;
            if(qp>=static_cast<std::int64_t>(query.size()))return std::nullopt;
            const char r=ref[static_cast<std::size_t>(rp)], q=query[static_cast<std::size_t>(qp)];
            if(r=='N'||q=='N')continue;
            ++count;if(r==q)++match;
        }
        if(count<128 || static_cast<double>(match)/count<0.42)return std::nullopt;
        matches+=match;total+=count;
    }
    if(total<512 || static_cast<double>(matches)/total<min_identity)return std::nullopt;
    return RescueRange{static_cast<std::uint32_t>(*x),static_cast<std::uint32_t>(*y),
                       support,static_cast<double>(matches)/total};
}
} // namespace ultralong
