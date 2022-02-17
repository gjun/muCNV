//
//  main_filter.cpp
//  muCNV
//
//  Created by Goo Jun on 11/24/18.
//  Copyright © 2018 Goo Jun. All rights reserved.
//

#include <stdio.h>
#include "in_vcf.h"
#include <fstream>

// TCLAP headers
#include "tclap/CmdLine.h"
#include "tclap/Arg.h"

double RO_THRESHOLD = 0.5;
double GENO_THRESHOLD = 0.9;

double concord(std::vector<int> &g1, std::vector<int> &g2)
{
    int N = 0;
    int N_match = 0;

    for(int i=0; i<g1.size(); ++i)
    {
        if (g1[i]>=0 && g2[i] >=0) // none of them is missing
        {
            if (g1[i]>0 && g2[i]>0)
            {
                N_match++;
            }
            if (g1[i]>0 || g2[i]>0)
            {
                N++;
            }
        }
    }

    return (double)N_match/(double)N;
}

void merge_overlaps(std::vector<sv> &all, std::vector<sv> &overlapping_svs, std::vector<std::vector <int> > &overlapping_genos)
{
    int N = overlapping_svs.size();
    std::vector< std::vector<double> > pair_RO(N);
    std::vector<sv> merged_candidates;

    for(int j=0; j<N; j++)
    {
        pair_RO[j] = std::vector<double>(N);
    }
    int clusters[N];
    // Merge based on RO within overlapping SV set
    for(int j=0;j<N-1;++j)
    {
        pair_RO[j][j] = 0;
        clusters[j] = j;
        for(int k=j+1;k<N; ++k)
        {
            pair_RO[j][k] = pair_RO[k][j] = mod_RO(overlapping_svs[j], overlapping_svs[k]);
        }
    }
    clusters[N-1] = N-1;
    pair_RO[N-1][N-1] = 0;

    for(int j=0; j<N-1;j++)
    {
        for(int k=j+1; k<N; k++)
        {
            if (pair_RO[j][k] >= RO_THRESHOLD && concord(overlapping_genos[j], overlapping_genos[k]) >= GENO_THRESHOLD)
            {
                clusters[k]=clusters[j];
            }
        }
    }
    for(int j=0;j<N;j++)
    {
        if (clusters[j] == j) // if not, it's already merged with a previous SV
        {
            merged_candidates.clear();
            merged_candidates.push_back(overlapping_svs[j]);
            for(int k=j+1;k<N;k++)
            {
                if (clusters[k] == j)
                {
                    merged_candidates.push_back(overlapping_svs[k]);
                }
            }
            if (merged_candidates.size()>1)
            {
                int best_idx = pick_sv_from_merged(merged_candidates);
                for(int k=0; k<merged_candidates.size(); ++k)
                {
                    if (k==best_idx)
                    {
                        merged_candidates[k].filter= "PASS";
                    }
                    else
                    {
                        merged_candidates[k].filter = "Collapsed";
                    }
                
                    all.push_back(merged_candidates[k]);
                }

            }
            else
            {
                merged_candidates[0].filter = "PASS";
                all.push_back(merged_candidates[0]);
            }
        }
    }
}

bool looks_good(sv &curr_sv, std::vector<int> &gt, std::vector< std::vector<double> > &dps, std::vector< std::vector<int> > &cnts)
{
    double sum_ref[4] = {0.0}, sum_het[4] = {0.0};
    int n_ref = 0, n_het=0;
//    double sumsq_ref[4], sumsq_het[4];

    for(int i=0; i<dps.size(); ++i)
    {
        dps[i][1] *= dps[i][0];
        dps[i][2] *= dps[i][3];

        if (gt[i] == 0)
        {
            for(int j=0; j<4; ++j)
            {
                sum_ref[j] += dps[i][j];
            }
            n_ref ++;
        }
        else if (gt[i] > 0)
        {
            for(int j=0; j<4; ++j)
            {
                sum_het[j] += dps[i][j];
            }
            n_het ++;
            
        }
    }
    double mean_ref[4], mean_het[4];
    for(int j=0; j<4; ++j)
    {
        if (n_ref > 0)
        {
            mean_ref[j] = sum_ref[j]/n_ref;
        }
        else
        {
            mean_ref[j] = 0;
        }
        if (n_het > 0)
        {
            mean_het[j] = sum_het[j]/n_het;
        }
        else
        {
            // std::cerr << " No allele n_ref " << n_ref << " n_het " << n_het << " ";
            curr_sv.filter = "AC0";
            return false; // No alleles, filter out
        }
    }
    if (mean_het[0] >1.3 || mean_het[3] > 1.3 || mean_het[0] < 0.7 || mean_het[3] < 0.7)
    {
       // std::cerr << " Peripheral low " << mean_ref[0] << "," << mean_ref[3] << " " << mean_het[0] << "," << mean_het[3] << " ";
        curr_sv.filter = "VarPrepostDepth";
        return false;  // Peripheral depth too low
    }
    else if (mean_ref[0] <0.7|| mean_ref[0] >1.3 || mean_ref[3] < 0.7 || mean_ref[3]>1.3)
    {
        curr_sv.filter = "RefPrepostDepth";
       // std::cerr << " Peripheral high " << mean_ref[0] << "," << mean_ref[3] << " " << mean_het[0] << "," << mean_het[3] << " ";
        return false;
    }
    if (curr_sv.multi_dp)
    {
        if (abs(mean_ref[1]-mean_ref[2]) > 0.5 || abs(mean_het[1]-mean_het[2])> 0.5 ) 
        {
         //   std::cerr << " MultiDP fail " << mean_ref[1] << "," << mean_ref[2] << " " << mean_het[1] << "," << mean_het[2] << " ";
            curr_sv.filter = "VarDepthDP2";
            return false;
        }
        else if (mean_ref[1] < 0.7 || mean_ref[2]<0.7 || mean_ref[1]>1.3 || mean_ref[2]>1.3)
        {
            curr_sv.filter = "RefDepthDP2";
          //  std::cerr << " Ref DP fail " << mean_ref[1] << "," << mean_ref[2] << " " << mean_het[1] << "," << mean_het[2] << " ";

            return false;
        }
    }
    

    return true;
}


void run_merging(std::vector<sv> &all, invcf &vcf)
{  
    std::vector<sv> overlapping_svs;
    std::vector<std::vector <int> > overlapping_genos;

    sv curr_sv;    
    std::vector<int> genos(vcf.n_sample); 
    std::vector<std::vector <double> > dps(vcf.n_sample);
    std::vector<std::vector <int> > cnts(vcf.n_sample);
    for(int i=0; i<vcf.n_sample; ++i)
    {
        dps[i] = std::vector<double>(4);
        cnts[i] = std::vector<int>(4);
    }

    bool b_fail = true;
    do
    {
        if (vcf.read_next(curr_sv) == 0)
        {
            if (curr_sv.svtype == DUP)
            {
                if (curr_sv.pre > 0.8 && curr_sv.pre < 1.2 && curr_sv.post > 0.8 && curr_sv.post <1.2)
                {
                    vcf.read_geno(curr_sv, genos);
                    vcf.read_dpcnt(curr_sv, dps, cnts);
                    if (looks_good(curr_sv, genos, dps, cnts))
                    {
                        b_fail = false;  
    //                    std::cerr <<  svTypeName(curr_sv.svtype) << " " << curr_sv.pos << "-" <<curr_sv.end << " supp " << curr_sv.supp << " multi_dp " << curr_sv.multi_dp;
    //                  std::cerr << " Pre " << curr_sv.pre << " Post " << curr_sv.post << std::endl;
                    }
                    else
                    {

                        all.push_back(curr_sv);
                    }
                //    else
                //   {
    //                  std::cerr <<  "FAIL: " << svTypeName(curr_sv.svtype) << " " << curr_sv.pos << "-" <<curr_sv.end << " supp " << curr_sv.supp << " multi_dp " << curr_sv.multi_dp;
    //                std::cerr << " Pre " << curr_sv.pre << " Post " << curr_sv.post << std::endl;
            //     }
                }
                else
                {
                    curr_sv.filter = "PrepostDepth";
                    all.push_back(curr_sv);
                }
            }
        }
        else
        {
            std::cerr << "Failed to read first valud SV" << std::endl;
            exit(-1);
        }
    } while (b_fail);
                    
    overlapping_svs.push_back(curr_sv);
    overlapping_genos.push_back(genos);
    
    int curr_end = curr_sv.end;
    int curr_chr = curr_sv.chrnum;
    
    while(vcf.read_next(curr_sv)==0)
    {
        if (curr_sv.svtype == DUP)
        {
            if (curr_sv.pre > 0.8 && curr_sv.pre < 1.2 && curr_sv.post > 0.8 && curr_sv.post <1.2)
            {
                vcf.read_geno(curr_sv, genos); 
                vcf.read_dpcnt(curr_sv, dps, cnts);
                if (looks_good(curr_sv, genos, dps, cnts))
                {
                    if (curr_sv.chrnum == curr_chr && curr_sv.pos < curr_end)
                    {
                        overlapping_svs.push_back(curr_sv);
                        overlapping_genos.push_back(genos);

                        if (curr_sv.end > curr_end)
                        {
                            curr_end = curr_sv.end;
                        }
                    }
                    else
                    {
                        if (overlapping_svs.size()>1)
                        {
                            merge_overlaps(all, overlapping_svs, overlapping_genos);
                        }
                        else
                        {
                            overlapping_svs[0].filter = "PASS";
                            all.push_back(overlapping_svs[0]);
                        }
                        overlapping_svs.clear();
                        overlapping_genos.clear();

                        curr_chr = curr_sv.chrnum;
                        curr_end = curr_sv.end;

                        overlapping_svs.push_back(curr_sv);
                        overlapping_genos.push_back(genos);
                    }
                }
                else
                {
                    all.push_back(curr_sv);
                }
            }
            else
            {
                curr_sv.filter = "PrepostDepth";
                all.push_back(curr_sv);
            }
        }
    }
    if (overlapping_svs.size()>1)
    {
        merge_overlaps(all, overlapping_svs, overlapping_genos);
    }
    else
    {
        overlapping_svs[0].filter = "PASS";
        all.push_back(overlapping_svs[0]);
    }
}

int main_filter(int argc, char** argv)
{
    std::string vcf_file;
    std::string interval_file;
    std::string out_filename;
    std::string supp_file;
    bool bMerge;
    bool bFilter;

    
    // Parsing command-line arguments
    try
    {
        TCLAP::CmdLine cmd("Command description message", ' ', "0.06");
        
        TCLAP::ValueArg<std::string> argOut("o","out","Output filename",false,"muCNV.vcf","string");
        TCLAP::ValueArg<std::string> argVcf("v","vcf","VCF file containing candidate SVs",false,"","string");
    //    TCLAP::ValueArg<std::string> argInterval("i","interval", "Binary interval file containing candidate SVs", false, "", "string");
    //    TCLAP::ValueArg<std::string> argSupp("S","Support","Support VCF file containing suppporting info",false,"","string");
    //    TCLAP::SwitchArg switchFilter("f", "filter", "Filter candidate discovery set using supporting VCF", cmd, false);
        TCLAP::SwitchArg switchMerge("m", "merge", "Merge candidate discovery SVs based on RO and stats", cmd, false);
        
        cmd.add(argOut);
        cmd.add(argVcf);
    //    cmd.add(argInterval);
    //    cmd.add(argSupp);

        cmd.parse(argc, argv);

        out_filename = argOut.getValue();
        vcf_file = argVcf.getValue();
    //    interval_file = argInterval.getValue();
    //    supp_file = argSupp.getValue();

    //    bFilter = switchFilter.getValue();
        bMerge = switchMerge.getValue();
    }
    catch (TCLAP::ArgException &e)
    {
        std::cerr << "Error: " << e.error() << " for arg " << e.argId() << std::endl;
        abort();
    }
        
    std::vector<std::string> sample_ids;
    std::vector<std::string> vcfs;
    std::vector<std::string> bam_names;

    if (bMerge)
    {
        // std::vector<sv> dels;
        // std::vector<sv> dups;
        // std::vector<sv> invs;
        std::vector<sv> all;
        
        
        // vcfs.push_back(vcf_file);        
        // read_intervals_from_vcf(sample_ids, vcfs, all);
        invcf myvcf;
        int n_sample = myvcf.open(vcf_file);
        std::cerr << "There are " << myvcf.n_sample << " samples." << std::endl;

        /*
        std::vector<int> genos(n_sample);
        sv tmp_sv;
        myvcf.read_next(tmp_sv, genos);
        return 0;

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
        std::sort(dups.begin(), dups.end());
        std::sort(invs.begin(), invs.end());

        run_merging(all, dels, RO_THRESHOLD);
        run_merging(all, dups, RO_THRESHOLD);
        run_merging(all, invs, RO_THRESHOLD);
        */
        run_merging(all, myvcf);
        std::sort(all.begin(), all.end());
        
        // Write to output file
        FILE *fp = fopen(out_filename.c_str(), "w");
        fprintf(fp, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
        for(int i=0; i<(int)all.size(); ++i)
        {
            fprintf(fp, "%d\t%d\t%s_%d:%d-%d\t.\t<%s>\t.\t%s\tSUPP=%d;SVTYPE=%s;END=%d\n", all[i].chrnum, all[i].pos, svTypeName(all[i].svtype).c_str(), all[i].chrnum, all[i].pos, all[i].end, svTypeName(all[i].svtype).c_str(), all[i].filter.c_str(), all[i].supp, svTypeName(all[i].svtype).c_str(), all[i].end);
        }
        fclose(fp);
    }
    /*
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
            std::string suppvec;
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
    }*/
    return 0;
}
