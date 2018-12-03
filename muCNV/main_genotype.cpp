//
//  main_genotype.cpp
//  muCNV
//
//  Created by Goo Jun on 11/23/18.
//  Copyright © 2018 Goo Jun. All rights reserved.
//

#include <stdio.h>
#include <map>

// TCLAP headers
#include "tclap/CmdLine.h"
#include "tclap/Arg.h"

#include "muCNV.h"
#include "gc_content.h"
#include "data_reader.h"

int main_genotype(int argc, char** argv)
{
    bool bNoHeader = false;
    bool bFail = false;
    int var_i = 0;
    
    string index_file;
    string vcf_file;
    string interval_file;
    string supp_file;
    string supp_id_file;
    string out_filename;
    string bam_file;
    string sChr;
    string gc_file;
    string sampID;
    string region;
    
    // Parsing command-line arguments
    try
    {
        TCLAP::CmdLine cmd("Command description message", ' ', "0.06");
        
        TCLAP::ValueArg<string> argOut("o","out","Output filename",false,"muCNV.vcf","string");
        TCLAP::ValueArg<string> argVcf("v","vcf","VCF file containing candidate SVs",false,"","string");
        TCLAP::ValueArg<string> argInterval("V","interVal", "Binary interval file containing candidate SVs", false, "", "string");
        TCLAP::ValueArg<int> argN("n","number", "N-th variant only", false, 466307, "integer");

//        TCLAP::ValueArg<string> argSupp("S","Support","Support VCF file containing suppporting info",false,"","string");
 //       TCLAP::ValueArg<string> argSuppID("I","IDinSupport","Sample ID list for support vectors",false,"","string");

        TCLAP::ValueArg<string> argIndex("i","index","List file containing list of intermediate pileups. Required with genotype option",false,"","string");
        TCLAP::ValueArg<string> argGcfile("f","gcFile","File containing GC content information",false, "GRCh38.gc", "string");

        TCLAP::ValueArg<string> argRegion("r", "region", "Genomic region (chr:start-end)", false, "", "string" );

        TCLAP::SwitchArg switchFail("a", "all", "Report filter failed variants", cmd, false);
        TCLAP::SwitchArg switchNoHeader("h", "headerless", "Do not print header in genoptyed VCF", cmd, false);
        
        cmd.add(argOut);
        cmd.add(argVcf);
        cmd.add(argInterval);
        cmd.add(argGcfile);
        cmd.add(argN);
        cmd.add(argIndex);
//        cmd.add(argSuppID);
 //       cmd.add(argSupp);

        cmd.add(argRegion);
        
        cmd.parse(argc, argv);
        
        var_i = argN.getValue();
        index_file = argIndex.getValue();
        out_filename = argOut.getValue();
        vcf_file = argVcf.getValue();
        interval_file = argInterval.getValue();
//        supp_file = argSupp.getValue();
//        supp_id_file = argSuppID.getValue();
        gc_file = argGcfile.getValue();
        bNoHeader = switchNoHeader.getValue();
        bFail = switchFail.getValue();

        region = argRegion.getValue();
        
        // TODO: other error checking...
        if (index_file == "")
        {
            std::cerr << "Error: list file is required for genotyping" << std::endl;
            exit(1);
        }
    }
    catch (TCLAP::ArgException &e)
    {
        std::cerr << "Error: " << e.error() << " for arg " << e.argId() << std::endl;
        abort();
    }
    
    // Multi-sample genotyping from summary VCFs
    

    
    std::vector<string> sample_ids;
    std::vector<string> vcfs;
    std::vector<string> pileup_names;
    std::vector<string> bam_names;
    
    GcContent gc;
    gc.initialize(gc_file);
    std::cerr << "GC content initialized" << std::endl;

    std::vector<breakpoint> vec_bp;
    std::vector<sv> vec_sv;
    
    if (vcf_file != "")
        read_svs_from_vcf(vcf_file, vec_bp, vec_sv);
    else if (interval_file != "")
        read_svs_from_intfile(interval_file, vec_bp, vec_sv);
    else
        std::cerr << "Error, VCF or Interval file is required." << std::endl;
    
    //vec_bp is not necessary for genotyping, but let's keep it for simplicity now
	std::cerr << vec_sv.size() << " SVs identified " << std::endl;
    
    read_list(index_file, pileup_names);
    int n_sample = 0;
    int n_pileup = (int) pileup_names.size();
    int n_var = (int) vec_sv.size();
    
    std::cerr << n_pileup << " pileup files identified from index file" << std::endl;
    
    DataReader reader;
    std::vector<SampleStat> stats;
    n_sample = reader.load(pileup_names, stats, gc);
    std::cerr << n_sample << " samples identified from pileup files\n" << std::endl;

//    for(int i=0; i<n_var; ++i)
    if (1)
    {
    //    int i=1000; // TEMPORARY FOR TESTING
//		int i=315268;
		int i=var_i; // a good one for testing
        
        vec_sv[i].print();

        std::vector< std::vector<double> > dvec_dp100 (n_sample);
        std::vector<double> var_dp (n_sample);
        
        int startpos = reader.read_depth100(vec_sv[i], dvec_dp100, gc);
        int endpos = startpos + (int)dvec_dp100[0].size() * 100;\
    
        string fname = "var" + std::to_string(i) + ".dp100.txt";
        
        FILE *fp = fopen(fname.c_str(), "wt");
        fprintf(fp, "index");
        for(int k=0; k<n_sample; ++k)
        {
            fprintf(fp, "\ts%d", k);
        }
        fprintf(fp, "\n");
        // TEMPORARY TEST --->
        for(int j=0; j<dvec_dp100[0].size(); ++j)
        {
            fprintf(fp, "%d",j);
            for(int k=0; k<n_sample; ++k)
            {
                fprintf(fp, "\t%f", dvec_dp100[k][j]);
            }
            fprintf(fp, "\n");

        }
        fclose(fp);
		/*
		std::cerr << reader.sample_ids[0];
		std::cerr << " " << vec_sv[i].chrnum;
		std::cerr << ":" << startpos << "-" << endpos << std::endl;
        fprintf(stderr, "%s: at chr%d:%d-%d\n", reader.sample_ids[0].c_str(), vec_sv[i].chrnum, startpos, endpos);
        
        for (int j=0; j<dvec_dp100[0].size(); ++j)
        {
            fprintf(stderr, "%f\n", dvec_dp100[0][j]);
        }
		*/
        // <---

		std::vector<ReadStat> rdstats (n_sample);
        reader.read_pair_split(vec_sv[i], rdstats, gc);

        reader.read_var_depth(i, var_dp);
        
        fname = "var" + std::to_string(i) + ".stat.txt";
        
        fp = fopen(fname.c_str(), "wt");
        fprintf(fp, "index\tprefr\tprerf\tpostfr\tpostrf\tprerpmiss\tpostrpmiss\tprespmiss\tpostspmiss\tpresplitout\tpresplitin\tpostsplitout\tpostsplitin\tdpvar\n");
        for(int j=0; j<n_sample; ++j)
        {
            fprintf(fp, "%d\t%d\t%d\t%d\t%d", j, rdstats[j].n_pre_FR, rdstats[j].n_pre_RF, rdstats[j].n_post_FR, rdstats[j].n_post_RF);
            fprintf(fp, "\t%d\t%d\t%d\t%d", rdstats[j].n_pre_rp_missing, rdstats[j].n_post_rp_missing, rdstats[j].n_pre_sp_missing, rdstats[j].n_post_sp_missing);
            fprintf(fp, "\t%d\t%d\t%d\t%d", rdstats[j].n_pre_split_out, rdstats[j].n_pre_split_in, rdstats[j].n_post_split_out, rdstats[j].n_post_split_in);
            fprintf(fp, "\t%f\n", var_dp[j]);
        }
        fclose(fp);
        // 0. Cluster var_dp 1-D
		// using var depth
        
        // 1. Cluster var_dp 2-D
		// Using dp100 
        
        // 2. Readpair / splitread counting
		// From individual genotyping, if RP/SP genotyped ones have expected avg. depth, they're genotyped
		// If not, consider them missing
        
        // 3. dp100 filtering / clustering (?)
		// median filtering and edge detection around breakpoint
    }
    
    return 0;
    
}

void obsolete()
{
    /*
     std::vector<double> avg_depths;
     std::vector<double> avg_isizes;
     std::vector<double> std_isizes;
     
     std::map<string, int> id_to_idx;
     
     
     invcfs V_list;
     
     V_list.initialize(vcfs, sample_ids, avg_depths, avg_isizes, std_isizes, region);
     
     n_sample = (int)sample_ids.size();
     
     for(int i=0;i<n_sample; ++i)
     {
     id_to_idx[sample_ids[i]] = i;
     }
     
     std::cerr << "Genotyping index loaded." << std::endl;
     std::cerr << n_sample << " samples identified." << std::endl;
     
     std::ifstream sfile;
     std::vector<int> id_map(n_sample, 0);
     
     if (supp_file != "")
     {
     if (supp_id_file != "")
     {
     std::cerr << "Supporting std::vector file and ID file provided." <<std::endl;
     }
     else
     {
     std::cerr << "Error. Supporting ID file is needed for supporting std::vector file." << std::endl;
     exit(1);
     }
     
     std::ifstream idfile(supp_id_file.c_str(), std::ios::in);
     for(int i=0;i<n_sample;++i)
     {
     string ln;
     getline(idfile, ln);
     if (id_to_idx.find(ln)!=id_to_idx.end())
     {
     id_map[i]  = id_to_idx[ln] ;
     }
     else
     {
     std::cerr << "Cannot find " << ln << " in ID" << std::endl;
     exit(1);
     }
     }
     idfile.close();
     }
     
     sfile.open(supp_file.c_str(), std::ios::in);
     
     outvcf vfile;
     
     vfile.open(out_filename);
     
     if (!bNoHeader)
     {
     vfile.write_header(sample_ids);
     }
     
     sv S;
     svdata D;
     svgeno G;
     
     D.set_size(n_sample);
     
     int val = 0;
     
     std::vector<sv> supp_list;
     std::vector<string> suppvec_list;
     
     while((val = V_list.read_interval_multi(S, D, region))>=0)
     {
     
     sv suppS;
     string suppvec;
     bool bMatch = false;
     
     
     if (!supp_list.empty())
     {
     if (supp_list[0].pos < S.pos)
     {
     supp_list.clear();
     suppvec_list.clear();
     }
     else
     {
     for(int i=0;i<(int)supp_list.size(); ++i)
     {
     if (supp_list[i].pos == S.pos && supp_list[i].end == S.end && supp_list[i].svtype == S.svtype)
     {
     bMatch = true;
     suppS = supp_list[i];
     suppvec = suppvec_list[i];
     break;
     }
     }
     }
     }
     if (!bMatch)
     {
     while(read_candidate_vcf(sfile, suppS, suppvec)>0)
     {
     if (suppS.pos == S.pos && !(suppS.svtype == S.svtype && suppS.end == S.end) )
     {
     supp_list.push_back(suppS);
     suppvec_list.push_back(suppvec);
     }
     else if (suppS.pos == S.pos && suppS.end == S.end && suppS.svtype == S.svtype)
     {
     break;
     }
     else if (suppS.pos > S.pos) // this should never happen
     {
     std::cerr << "Something wrong" <<std::endl;
     std::cerr << "Curren position is " << S.pos << "-" << S.end << " while supporting std::vector has passed " << suppS.pos << std::endl;
     exit(1);
     }
     }
     }
     std::vector<double> wt(n_sample, 1);
     double add_wt = 0;
     
     if (suppS.supp > 1 )
     {
     switch(suppS.supp)
     {
     case 2:
     add_wt = 2;
     break;
     case 3:
     add_wt = 5;
     break;
     case 4:
     add_wt = 10;
     break;
     case 5:
     add_wt = 20;
     break;
     }
     for(int i=0;i<n_sample; ++i)
     {
     if (suppvec[i] == '1')
     {
     // Give additional weights to samples found in supp_vec
     wt[id_map[i]] += add_wt;
     }
     }
     }
     
     if (val>0  && !in_centrome(S) )
     {
     gtype T;
     G.initialize(n_sample);
     D.normalize(S, avg_depths, avg_isizes);
     G.b_pass = false;
     
     if (S.svtype == DEL)
     {
     T.call_del(S, D, G, avg_isizes, std_isizes, wt);
     }
     else if (S.svtype == DUP || S.svtype == CNV)
     {
     T.call_cnv(S, D, G, avg_isizes, std_isizes, wt);
     }
     
     if (S.svtype != INV  && (bFail || G.b_pass))
     //                if (G.b_pass || bFail)
     {
     string ln;
     G.info += ";SUPP=" + std::to_string(suppS.supp);
     ln.reserve(n_sample*30);
     G.print(S, D, ln, wt);
     vfile.print(ln);
     }
     
     }
     
     }
     vfile.close();
     */
}
