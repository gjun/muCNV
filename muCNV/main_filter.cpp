//
//  main_filter.cpp
//  muCNV
//
//  Created by Goo Jun on 11/24/18.
//  Copyright © 2018 Goo Jun. All rights reserved.
//

#include <stdio.h>
#include "muCNV.h"
// TCLAP headers
#include "tclap/CmdLine.h"
#include "tclap/Arg.h"

int main_filter(int argc, char** argv)
{
    string vcf_file;
    string interval_file;
    string out_filename;
    string supp_file;
    bool bMerge;
    bool bFilter;
    
    std::vector<string> sample_ids;
    std::vector<string> vcfs;
    std::vector<string> bam_names;
    
    // Parsing command-line arguments
    try
    {
        TCLAP::CmdLine cmd("Command description message", ' ', "0.06");
        
        TCLAP::ValueArg<string> argOut("o","out","Output filename",false,"muCNV.vcf","string");
        TCLAP::ValueArg<string> argVcf("v","vcf","VCF file containing candidate SVs",false,"","string");
        TCLAP::ValueArg<string> argInterval("i","interval", "Binary interval file containing candidate SVs", false, "", "string");
        TCLAP::ValueArg<string> argSupp("S","Support","Support VCF file containing suppporting info",false,"","string");
        TCLAP::SwitchArg switchFilter("f", "filter", "Filter candidate discovery set using supporting VCF", cmd, false);
        TCLAP::SwitchArg switchMerge("m", "merge", "Merge candidate discovery SVs based on RO", cmd, false);
        
        cmd.add(argOut);
        cmd.add(argVcf);
        cmd.add(argInterval);
        cmd.add(argSupp);

        cmd.add(switchFilter);
        cmd.add(switchMerge);
        
        cmd.parse(argc, argv);

        out_filename = argOut.getValue();
        vcf_file = argVcf.getValue();
        interval_file = argInterval.getValue();
        supp_file = argSupp.getValue();

        bFilter = switchFilter.getValue();
        bMerge = switchMerge.getValue();
    }
    catch (TCLAP::ArgException &e)
    {
        std::cerr << "Error: " << e.error() << " for arg " << e.argId() << std::endl;
        abort();
    }
    
    if (bMerge)
    {
        double RO_THRESHOLD = 0.8;
        
        std::vector<sv> dels;
        std::vector<sv> dups;
        std::vector<sv> invs;
        std::vector<sv> all;
        
        vcfs.push_back(vcf_file);
        
        read_intervals_from_vcf(sample_ids, vcfs, all);
        for(int i=0;i<(int)all.size();++i)
        {
            if (all[i].svtype == DEL)
            {
                dels.push_back(all[i]);
            }
            else if (all[i].svtype == DUP)
            {
                dups.push_back(all[i]);
            }
            else if (all[i].svtype == INV)
            {
                invs.push_back(all[i]);
            }
        }
        all.clear();
        
        std::sort(dels.begin(), dels.end());
        
        std::vector<sv> merged_candidates;
        merged_candidates.push_back(dels[0]);
        for(int i=1; i<(int)dels.size(); ++i)
        {
            if (RO(dels[i], merged_candidates.back() ) > RO_THRESHOLD)
            {
                merged_candidates.push_back(dels[i]);
            }
            else
            {
                sv best_sv;
                pick_sv_from_merged(best_sv, merged_candidates);
                all.push_back(best_sv);
                
                merged_candidates.clear();
                merged_candidates.push_back(dels[i]);
                
            }
        }
        if (merged_candidates.size() > 0)
        {
            sv best_sv;
            pick_sv_from_merged(best_sv, merged_candidates);
            all.push_back(best_sv);
        }
        merged_candidates.clear();
        
        merged_candidates.push_back(dups[0]);
        for(int i=1; i<(int)dups.size(); ++i)
        {
            if (RO(dups[i], merged_candidates.back() ) > RO_THRESHOLD)
            {
                merged_candidates.push_back(dups[i]);
            }
            else
            {
                sv best_sv;
                pick_sv_from_merged(best_sv, merged_candidates);
                all.push_back(best_sv);
                
                merged_candidates.clear();
                merged_candidates.push_back(dups[i]);
                
            }
        }
        if (merged_candidates.size() > 0)
        {
            sv best_sv;
            pick_sv_from_merged(best_sv, merged_candidates);
            all.push_back(best_sv);
        }
        merged_candidates.clear();
        
        merged_candidates.push_back(invs[0]);
        for(int i=1; i<(int)invs.size(); ++i)
        {
            if (RO(invs[i], merged_candidates.back() ) > RO_THRESHOLD)
            {
                merged_candidates.push_back(invs[i]);
            }
            else
            {
                sv best_sv;
                pick_sv_from_merged(best_sv, merged_candidates);
                all.push_back(best_sv);
                
                merged_candidates.clear();
                merged_candidates.push_back(invs[i]);
                
            }
        }
        if (merged_candidates.size() > 0)
        {
            sv best_sv;
            pick_sv_from_merged(best_sv, merged_candidates);
            all.push_back(best_sv);
        }
        merged_candidates.clear();
        
        std::sort(all.begin(), all.end());
        
        // Write to output file
        FILE *fp = fopen(out_filename.c_str(), "w");
        fprintf(fp, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
        for(int i=0; i<(int)all.size(); ++i)
        {
            fprintf(fp, "%d\t%d\tSV%d\t.\t<%s>\t.\tPASS\tSUPP=%d;SVTYPE=%s;END=%d\n", all[i].chrnum, all[i].pos, i+1, svTypeName(all[i].svtype).c_str(), all[i].supp, svTypeName(all[i].svtype).c_str(), all[i].end);
        }
        fclose(fp);
    }
    else if (bFilter)
    {
        std::vector<sv> dels;
        std::vector<sv> dups;
        std::vector<sv> invs;
        std::vector<sv> all;
        
        vcfs.push_back(supp_file);
        read_intervals_from_vcf(sample_ids, vcfs, all);
        FILE *fp = fopen(out_filename.c_str(), "w");
        
        fprintf(fp, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
        
        for(int i=0;i<(int)all.size();++i)
        {
            if (all[i].svtype == DEL)
            {
                dels.push_back(all[i]);
            }
            else if (all[i].svtype == DUP)
            {
                dups.push_back(all[i]);
            }
            else if (all[i].svtype == INV)
            {
                invs.push_back(all[i]);
            }
        }
        all.clear();
        
        std::ifstream vfile(vcf_file.c_str(), std::ios::in);
        int cnt = 1;
        
        while(vfile.good())
        {
            sv S;
            string suppvec;
            int idx = -1;
            
            if (read_candidate_vcf(vfile, S, suppvec)>0)
            {
                sv o_S;
                
                if (S.svtype == DEL)
                {
                    idx = find_overlap_sv(S, dels);
                    if (idx >=0)
                    {
                        o_S = dels[idx];
                    }
                }
                else if (S.svtype == DUP || S.svtype==CNV)
                {
                    idx = find_overlap_sv(S, dups);
                    if (idx >=0)
                    {
                        o_S = dups[idx];
                    }
                }
                else if (S.svtype == INV)
                {
                    idx = find_overlap_sv(S, invs);
                    if (idx >=0)
                    {
                        o_S = invs[idx];
                    }
                }
                //print record
                if (idx>=0)
                {
                    fprintf(fp, "%d\t%d\tSV%d\t.\t<%s>\t.\tPASS\tSVTYPE=%s;OVERLAP=%d:%d-%d;SUPP=%d;N=%d;END=%d;SUPP_VEC=%s\n", S.chrnum, S.pos, cnt++, svTypeName(S.svtype).c_str(), svTypeName(S.svtype).c_str(), o_S.chrnum, o_S.pos, o_S.end, o_S.supp, S.supp, S.end, suppvec.c_str());
                }
                else
                {
                    fprintf(fp, "%d\t%d\tSV%d\t.\t<%s>\t.\tPASS\tSVTYPE=%s;SUPP=1;N=%d;END=%d;SUPP_VEC=%s\n", S.chrnum, S.pos, cnt++, svTypeName(S.svtype).c_str(), svTypeName(S.svtype).c_str(), S.supp, S.end, suppvec.c_str());
                }
            }
        }
        fclose(fp);
        vfile.close();
    }
    return 0;
}
