#include "doctest/doctest.h"
#include "align/anchored_insertions.hpp"

#include <stdexcept>
#include <string>
#include <vector>

namespace {
std::string ungap(const std::string& s) {
    std::string out;
    for (char c : s) if (c != '-') out.push_back(c);
    return out;
}
}

TEST_SUITE("anchored_insertions") {
    TEST_CASE("two separate loci preserve every consensus anchor and query base") {
        using namespace align::anchored_insertions;
        const std::vector<Query> queries = {
            {"q1", "ACTGT", "2M1I2M", 0},
            {"q2", "ACGTT", "4M1I", 1},
            {"q3", "AGT", "1M1D2M", 2}
        };
        const auto msa = ::align::anchored_insertions::align("ACGT", queries);
        REQUIRE(msa.rows.size() == 4);
        CHECK(msa.occupied_sites == 2);
        CHECK(msa.columns == 6);
        CHECK(msa.local_msas == 0); // each site has only one insertion
        CHECK(msa.rows[0].bases == "AC-GT-");
        CHECK(msa.rows[1].bases == "ACTGT-");
        CHECK(msa.rows[2].bases == "AC-GTT");
        CHECK(msa.rows[3].bases == "A--GT-");
    }

    TEST_CASE("same-locus insertions use bounded in-process abPOA") {
        using namespace align::anchored_insertions;
        const std::vector<Query> queries = {
            {"q1", "ACGT", "1M1I2M", 0}, // A C GT
            {"q2", "ACCGT", "1M2I2M", 1} // A CC GT
        };
        const auto msa = ::align::anchored_insertions::align("AGT", queries);
        REQUIRE(msa.rows.size() == 3);
        CHECK(msa.occupied_sites == 1);
        CHECK(msa.local_msas == 1);
        for (std::size_t i = 0; i < queries.size(); ++i) {
            CHECK(msa.rows[i+1].bases.size() == msa.columns);
            CHECK(ungap(msa.rows[i+1].bases) == queries[i].bases);
        }
        CHECK(ungap(msa.rows[0].bases) == "AGT");
    }

    TEST_CASE("invalid SAM CIGAR and duplicate IDs fail closed") {
        using namespace align::anchored_insertions;
        CHECK_THROWS(::align::anchored_insertions::align("ACGT", {{"x", "ACGT", "3M", 0}}));
        CHECK_THROWS(::align::anchored_insertions::align("ACGT", {{"x", "ACGT", "4S", 0}}));
        CHECK_THROWS(::align::anchored_insertions::align("ACGT", {{"x", "ACGT", "4M", 1}}));
        CHECK_THROWS(::align::anchored_insertions::align("ACGT", {{"x", "ACGT", "4M", 0},
                                    {"x", "ACGT", "4M", 1}}));
        CHECK_THROWS(::align::anchored_insertions::align("ACGT", {{"x", "ACGT", "1M1I3M", 0}}));
        CHECK_THROWS(::align::anchored_insertions::align("ACGT", {{"x", "acgt", "4M", 0}}));
    }

    TEST_CASE("bounded output and query counts prevent large allocations") {
        using namespace align::anchored_insertions;
        Options tight;
        tight.max_output_cells = 4;
        CHECK_THROWS(::align::anchored_insertions::align("ACGT", {{"x", "ACGT", "4M", 0}}, tight));
        tight = Options{};
        tight.max_sequences = 1;
        CHECK_THROWS(::align::anchored_insertions::align("ACGT", {{"x", "ACGT", "4M", 0},
                                    {"y", "ACGT", "4M", 1}}, tight));
    }
}


TEST_SUITE("soft_anchors") {
    TEST_CASE("conservative split separates unrelated center haplotypes and preserves all bases") {
        using namespace align::anchored_insertions;
        // Distinct deterministic pseudo-random DNA prevents misleadingly low
        // k-mer overlaps caused by a degenerate homopolymer test fixture.
        unsigned state = 37;
        auto dna = [&](std::size_t len) {
            std::string out;
            for (std::size_t i=0; i<len; ++i) {
                state = state*1664525u + 1013904223u;
                out += "ACGT"[(state >> 24)&3];
            }
            return out;
        };
        const auto left = dna(20), right = dna(20);
        const auto center_hap = dna(80), alternate = dna(60);
        const std::string center = left + center_hap + right;
        std::vector<Query> queries;
        for (std::size_t i=0; i<3; ++i) {
            queries.push_back({"alt"+std::to_string(i),left+alternate+right,
                               "80M20D20M",queries.size()});
            queries.push_back({"ref"+std::to_string(i),center,
                               "120M",queries.size()});
        }
        for (std::size_t i=0; i<3; ++i) {
            queries.push_back({"empty"+std::to_string(i),left+right,
                               "20M80D20M",queries.size()});
        }
        // Residual singleton with a partial center haplotype: tolerated only
        // because the two full-length cores each contain three sequences.
        const auto orphan_fragment = dna(20);
        queries.push_back({"orphan",left+orphan_fragment+right,
                           "20M50D20M10D20M",queries.size()});
        Options base;
        auto plain = ::align::anchored_insertions::align(center,queries,base);
        base.soft_anchors=true;
        auto rescued = ::align::anchored_insertions::align(center,queries,base);
        REQUIRE(rescued.rows.size()==queries.size()+1);
        CHECK(rescued.soft_anchor_candidates >= 1);
        CHECK(rescued.soft_anchor_splits == 1);
        CHECK(rescued.columns > plain.columns);
        for (std::size_t i=0;i<queries.size();++i)
            CHECK(ungap(rescued.rows[i+1].bases)==queries[i].bases);
        CHECK(ungap(rescued.rows[0].bases)==center);
        // The two independent haplotypes occupy disjoint lanes, preserving
        // shared left and right genomic anchors.
        const auto& alt = rescued.rows[1].bases;
        const auto& ref = rescued.rows[2].bases;
        std::size_t shared = 0;
        for (std::size_t i=20; i<rescued.columns-20; ++i)
            shared += (alt[i]!='-' && ref[i]!='-');
        CHECK(shared == 0);
        CHECK(alt.substr(0,20)==left);
        CHECK(ref.substr(rescued.columns-20)==right);
    }

    TEST_CASE("no supported competing groups falls back to byte-identical PR3B5A") {
        using namespace align::anchored_insertions;
        const std::vector<Query> queries = {
            {"q1", "ACGT", "1M1I2M", 0},
            {"q2", "ACCGT", "1M2I2M", 1}
        };
        const auto old = ::align::anchored_insertions::align("AGT",queries);
        Options opt;
        opt.soft_anchors = true;
        const auto candidate = ::align::anchored_insertions::align("AGT",queries,opt);
        CHECK(candidate.soft_anchor_splits == 0);
        CHECK(candidate.rows.size()==old.rows.size());
        for (std::size_t i=0; i<old.rows.size();++i)
            CHECK(candidate.rows[i].bases==old.rows[i].bases);
    }

    TEST_CASE("soft output cell budget is enforced") {
        using namespace align::anchored_insertions;
        Options opt;
        opt.soft_anchors = true;
        opt.max_output_cells = 1;
        CHECK_THROWS(::align::anchored_insertions::align("ACGT",{{"q","ACGT","4M",0}},opt));
    }

    TEST_CASE("low complexity competing haplotypes are not split by sparse k-mers") {
        using namespace align::anchored_insertions;
        const std::string left(20,'G'), right(20,'T'),
                          center_hap(80,'C'), alternative(60,'A');
        const std::string center=left+center_hap+right;
        std::vector<Query> q;
        for (int i=0;i<3;++i) {
            q.push_back({"a"+std::to_string(i),left+alternative+right,
                         "80M20D20M",q.size()});
            q.push_back({"c"+std::to_string(i),center,
                         "120M",q.size()});
            q.push_back({"n"+std::to_string(i),left+right,
                         "20M80D20M",q.size()});
        }
        Options soft; soft.soft_anchors=true;
        const auto old=::align::anchored_insertions::align(center,q);
        const auto got=::align::anchored_insertions::align(center,q,soft);
        CHECK(got.soft_anchor_splits==0);
        CHECK(got.rows.size()==old.rows.size());
        for (std::size_t i=0;i<old.rows.size();++i)
            CHECK(got.rows[i].bases == old.rows[i].bases);
    }

    TEST_CASE("homologous 60-base prefix of 85-base extension is not split") {
        using namespace align::anchored_insertions;
        unsigned state = 778899;
        auto dna = [&](std::size_t len) {
            std::string out;
            for (std::size_t i=0;i<len;++i) {
                state = state*1664525u+1013904223u;
                out += "ACGT"[(state>>24)&3];
            }
            return out;
        };
        const auto left=dna(20), right=dna(20);
        const auto shared=dna(60), extension=dna(25);
        const auto center=left+shared+extension+right;
        std::vector<Query> queries;
        for (std::size_t i=0;i<3;++i) {
            queries.push_back({"short"+std::to_string(i),left+shared+right,
                               "80M25D20M",queries.size()});
            queries.push_back({"long"+std::to_string(i),center,
                               "125M",queries.size()});
            queries.push_back({"none"+std::to_string(i),left+right,
                               "20M85D20M",queries.size()});
        }
        Options soft;soft.soft_anchors=true;
        const auto old=::align::anchored_insertions::align(center,queries);
        const auto got=::align::anchored_insertions::align(center,queries,soft);
        CHECK(got.soft_anchor_candidates>=1);
        CHECK(got.soft_anchor_splits==0);
        for (std::size_t i=0;i<old.rows.size();++i)
            CHECK(got.rows[i].bases==old.rows[i].bases);
    }


}

TEST_SUITE("soft_anchor_recruit") {
    TEST_CASE("V3 recovers mutated independent 60/85-base haplotypes when V2 refuses") {
        using namespace align::anchored_insertions;
        unsigned rng=0x12345678u;
        auto next = [&]() { rng=rng*1664525u+1013904223u; return rng; };
        auto dna = [&](std::size_t length) {
            std::string s;
            for (std::size_t i=0;i<length;++i) s += "ACGT"[(next()>>24)&3];
            return s;
        };
        auto mutate = [&](std::string s) {
            // 3% independent substitutions: 7-mer strict linkage tends to
            // fragment each family; 5-mer family recruitment can recover it.
            for (std::size_t i=0;i<s.size();++i) {
                if (next()%100 < 3) {
                    const char original=s[i];
                    const auto idx=std::string("ACGT").find(original);
                    s[i]="ACGT"[(idx+1+(next()%3))%4];
                }
            }
            return s;
        };
        const std::string left=dna(40), right=dna(40);
        const std::string center_hap=dna(85), alt_hap=dna(60);
        const std::string center=left+center_hap+right;
        std::vector<Query> queries;
        for (int i=0;i<9;++i) {
            const auto seq=left+mutate(alt_hap)+right;
            queries.push_back({"alt"+std::to_string(i),seq,
                               "100M25D40M",queries.size()});
        }
        for (int i=0;i<6;++i) {
            const auto seq=left+mutate(center_hap)+right;
            queries.push_back({"ref"+std::to_string(i),seq,
                               "165M",queries.size()});
        }
        for (int i=0;i<9;++i)
            queries.push_back({"empty"+std::to_string(i),left+right,
                               "40M85D40M",queries.size()});
        Options v2; v2.soft_anchors=true;
        const auto old=::align::anchored_insertions::align(center,queries,v2);
        Options v3=v2; v3.soft_anchor_recruit=true;
        const auto got=::align::anchored_insertions::align(center,queries,v3);
        CHECK(old.soft_anchor_recruits==0);
        CHECK(got.soft_anchor_candidates>=1);
        CHECK(got.soft_anchor_splits>=1);
        CHECK(got.soft_anchor_recruits>=1);
        CHECK(got.columns>old.columns);
        REQUIRE(got.rows.size()==queries.size()+1);
        CHECK(ungap(got.rows[0].bases)==center);
        for (std::size_t i=0;i<queries.size();++i) {
            CHECK(ungap(got.rows[i+1].bases)==queries[i].bases);
            CHECK(got.rows[i+1].bases.size()==got.columns);
        }
        // Two nonhomologous families must occupy disjoint block lanes.
        std::size_t shared=0;
        for (std::size_t col=40;col<got.columns-40;++col)
            shared += (got.rows[1].bases[col]!='-' &&
                       got.rows[10].bases[col]!='-');
        CHECK(shared==0);
    }

    TEST_CASE("V3 does not split a homologous 60-base prefix of an 85-base extension") {
        using namespace align::anchored_insertions;
        unsigned state=778899u;
        auto dna = [&](std::size_t n) {
            std::string s;
            for (std::size_t i=0;i<n;++i) {
                state=state*1664525u+1013904223u;
                s += "ACGT"[(state>>24)&3];
            }
            return s;
        };
        const auto left=dna(20),right=dna(20),shared=dna(60),tail=dna(25);
        const auto center=left+shared+tail+right;
        std::vector<Query> queries;
        for (int i=0;i<4;++i) {
            queries.push_back({"short"+std::to_string(i),left+shared+right,
                               "80M25D20M",queries.size()});
            queries.push_back({"long"+std::to_string(i),center,
                               "125M",queries.size()});
            queries.push_back({"none"+std::to_string(i),left+right,
                               "20M85D20M",queries.size()});
        }
        Options v2; v2.soft_anchors=true;
        const auto old=::align::anchored_insertions::align(center,queries,v2);
        Options v3=v2;v3.soft_anchor_recruit=true;
        const auto got=::align::anchored_insertions::align(center,queries,v3);
        CHECK(got.soft_anchor_splits==0);
        CHECK(got.soft_anchor_recruits==0);
        CHECK(got.rows.size()==old.rows.size());
        for (std::size_t i=0;i<old.rows.size();++i)
            CHECK(got.rows[i].bases==old.rows[i].bases);
    }


    TEST_CASE("V3 low-complexity competing homopolymers do not create spurious families") {
        using namespace align::anchored_insertions;
        const std::string left(30,'G'),right(30,'T'),
                          reference(85,'C'),alternative(60,'A');
        const auto center=left+reference+right;
        std::vector<Query> queries;
        for (int i=0;i<6;++i) {
            queries.push_back({"alt"+std::to_string(i),left+alternative+right,
                               "90M25D30M",queries.size()});
            queries.push_back({"ref"+std::to_string(i),center,
                               "145M",queries.size()});
            queries.push_back({"none"+std::to_string(i),left+right,
                               "30M85D30M",queries.size()});
        }
        Options v2;v2.soft_anchors=true;
        const auto before=::align::anchored_insertions::align(center,queries,v2);
        Options v3=v2;v3.soft_anchor_recruit=true;
        const auto after=::align::anchored_insertions::align(center,queries,v3);
        CHECK(after.soft_anchor_recruits==0);
        CHECK(after.soft_anchor_splits==before.soft_anchor_splits);
        REQUIRE(after.rows.size()==before.rows.size());
        for (std::size_t i=0;i<before.rows.size();++i)
            CHECK(after.rows[i].bases==before.rows[i].bases);
    }

    TEST_CASE("V3 leaves V2-accepted high-confidence splits byte-identical") {
        using namespace align::anchored_insertions;
        unsigned state=37;
        auto dna=[&](std::size_t n) {
            std::string s;
            for (std::size_t i=0;i<n;++i) {
                state=state*1664525u+1013904223u;
                s += "ACGT"[(state>>24)&3];
            }
            return s;
        };
        const auto left=dna(20),right=dna(20);
        const auto reference=dna(80),alternative=dna(60);
        const auto center=left+reference+right;
        std::vector<Query> queries;
        for (int i=0;i<4;++i) {
            queries.push_back({"alt"+std::to_string(i),left+alternative+right,
                               "80M20D20M",queries.size()});
            queries.push_back({"ref"+std::to_string(i),center,
                               "120M",queries.size()});
            queries.push_back({"none"+std::to_string(i),left+right,
                               "20M80D20M",queries.size()});
        }
        Options v2;v2.soft_anchors=true;
        const auto before=::align::anchored_insertions::align(center,queries,v2);
        REQUIRE(before.soft_anchor_splits==1);
        Options v3=v2;v3.soft_anchor_recruit=true;
        const auto after=::align::anchored_insertions::align(center,queries,v3);
        CHECK(after.soft_anchor_splits==1);
        CHECK(after.soft_anchor_recruits==0);
        REQUIRE(after.rows.size()==before.rows.size());
        for (std::size_t i=0;i<before.rows.size();++i)
            CHECK(after.rows[i].bases==before.rows[i].bases);
    }
}
