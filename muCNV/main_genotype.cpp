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

#include "in_vcf.h"
#include "gc_content.h"
#include "data_reader.h"
#include "genotyper.h"

void write_varstat(sv& curr_sv, std::vector<SampleStat> &stats, std::vector<ReadStat> &rdstats, std::vector<double> &var_dp)
{
    int n_sample = (int) stats.size();
    
    std::string fname = svTypeName(curr_sv.svtype) + "_" + std::to_string(curr_sv.chrnum) + "_" + std::to_string(curr_sv.pos) + "-" + std::to_string(curr_sv.end) + ".stat.txt";
    
    FILE* fp = fopen(fname.c_str(), "wt");
    
    fprintf(fp, "index\tavgdp\tstddp\tavgisz\tstdisz\tprefr\tprerf\tpostfr\tpostrf\tprerpmiss\tpostrpmiss\tprespmiss\tpostspmiss\tpresplitout\tpresplitin\tpostsplitout\tpostsplitin\tdpvar\n");
    for(int j=0; j<n_sample; ++j)
    {
        fprintf(fp, "%d\t%f\t%f\t%f\t%f", j, stats[j].avg_dp, stats[j].std_dp, stats[j].avg_isize, stats[j].std_isize);
        fprintf(fp, "\t%d\t%d\t%d\t%d", rdstats[j].n_pre_FR, rdstats[j].n_pre_RF, rdstats[j].n_post_FR, rdstats[j].n_post_RF);
        fprintf(fp, "\t%d\t%d\t%d\t%d", rdstats[j].n_pre_rp_missing, rdstats[j].n_post_rp_missing, rdstats[j].n_pre_sp_missing, rdstats[j].n_post_sp_missing);
        fprintf(fp, "\t%d\t%d\t%d\t%d", rdstats[j].n_pre_split_out, rdstats[j].n_pre_split_in, rdstats[j].n_post_split_out, rdstats[j].n_post_split_in);
        fprintf(fp, "\t%f\n", var_dp[j]);
    }
    fclose(fp);
}

int main_genotype(int argc, char** argv)
{
    bool bNoHeader = false;
    bool bFail = false;
    
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
    string range;
    bool b_dumpstat;
    
    // Parsing command-line arguments
    try
    {
        TCLAP::CmdLine cmd("Command description message", ' ', "0.06");
        
        TCLAP::ValueArg<string> argOut("o","out","Output filename",false,"muCNV.vcf","string");
        TCLAP::ValueArg<string> argVcf("v","vcf","VCF file containing candidate SVs",false,"","string");
        TCLAP::ValueArg<string> argInterval("V","interVal", "Binary interval file containing candidate SVs", false, "", "string");
        TCLAP::ValueArg<string> argRange("n","numbers", "variants in range (from-to) only", false, "", "integer-integer");

        TCLAP::ValueArg<string> argIndex("i","index","List file containing list of intermediate pileups. Required with genotype option",false,"","string");
        TCLAP::ValueArg<string> argGcfile("f","gcFile","File containing GC content information",false, "GRCh38.gc", "string");

        TCLAP::ValueArg<string> argRegion("r", "region", "Genotype specific genomic region", false, "", "chr:startpos-endpos" );

        TCLAP::SwitchArg switchFail("a", "all", "Report filter failed variants", cmd, false);
        TCLAP::SwitchArg switchDumpstat("d", "dumpstat", "dump detailed statistics of variants (warning: large output)", cmd, false);
        TCLAP::SwitchArg switchNoHeader("l", "lessheader", "Do not print header in genoptyed VCF", cmd, false);
        
        cmd.add(argOut);
        cmd.add(argVcf);
        cmd.add(argInterval);
        cmd.add(argGcfile);
        cmd.add(argIndex);
        cmd.add(argRegion);
        cmd.add(argRange);
        cmd.parse(argc, argv);
        
        range = argRange.getValue();
        index_file = argIndex.getValue();
        out_filename = argOut.getValue();
        vcf_file = argVcf.getValue();
        interval_file = argInterval.getValue();
        gc_file = argGcfile.getValue();
        bNoHeader = switchNoHeader.getValue();
        bFail = switchFail.getValue();
        b_dumpstat = switchDumpstat.getValue();

        region = argRegion.getValue();
        
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
    
    int n_start = 0;
    int n_end = n_var;
    
    std::cerr << n_pileup << " pileup files identified from index file" << std::endl;
    
    DataReader reader;
    std::vector<SampleStat> stats;
    n_sample = reader.load(pileup_names, stats, gc);
    std::cerr << n_sample << " samples identified from pileup files\n" << std::endl;

    if (range != "")
    {
        std::string::size_type sz;   // alias of size_t
        
        n_start = std::stoi(range, &sz);
        n_end = std::stoi(range.substr(sz+1));
        std::cerr << "Variants index from " << n_start << " to " << n_end << " will be genotyepd" << std::endl;
    }
    
    for(int i=n_start; i<=n_end; ++i)
    {
        //vec_sv[i].print();

        SvGeno G(n_sample);
        SvData D(n_sample);
        Genotyper gtyper;
        
        reader.read_depth100(vec_sv[i], D.dp2, gc, b_dumpstat);
        std::vector<ReadStat> rdstats (n_sample);
        reader.read_pair_split(vec_sv[i], D.rdstats, gc);
        reader.read_var_depth(i, D.var_depth);
        
        
        for(int j=0; j<n_sample; ++j)
        {
            for(int k=0; k<D.dp2.size(); ++k)
            {
                D.dp2[k][j] /= (double)stats[j].avg_dp;
            }
            D.var_depth[j] /= (double)stats[j].avg_dp;
        }

        if (b_dumpstat)
        {
            write_varstat(vec_sv[i], stats, D.rdstats, D.var_depth);
        }
        
        gtyper.call(vec_sv[i], D, G);
        std::cout << gtyper.print(vec_sv[i], D, G);
    }
    
    return 0;
    
}
