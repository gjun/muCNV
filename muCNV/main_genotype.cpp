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
#include "out_vcf.h"
#include "gc_content.h"
#include "data_reader.h"
#include "genotyper.h"
#include "common.h"

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
	double max_p;
    bool b_dumpstat;
	bool b_kmeans;
	bool b_mahalanobis;
	int chr;
    
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
        TCLAP::ValueArg<double> argChr("c","chr","Chromosome number (1-24) if pileup contains only a single chromosome",false,0,"integer (1-24)");

        TCLAP::ValueArg<double> argPoverlap("p","pmax","Maximum overlap between depth clusters",false,0.2,"number(0-1.0)");

        TCLAP::SwitchArg switchFail("a", "all", "Report filter failed variants", cmd, false);
        TCLAP::SwitchArg switchDumpstat("d", "dumpstat", "dump detailed statistics of variants (warning: large output)", cmd, false);
        TCLAP::SwitchArg switchNoHeader("l", "lessheader", "Do not print header in genoptyed VCF", cmd, false);
        TCLAP::SwitchArg switchKmeans("k", "Kmeans", "Use K-Means instead of EM", cmd, false);
        TCLAP::SwitchArg switchMahalanobis("m","Mahalanobis", "Use Mahalanobis distance (use with K-means)", cmd, false);
        
        cmd.add(argOut);
        cmd.add(argVcf);
        cmd.add(argInterval);
        cmd.add(argGcfile);
        cmd.add(argIndex);
		cmd.add(argChr);
        cmd.add(argRegion);
        cmd.add(argRange);
		cmd.add(argPoverlap);
        cmd.parse(argc, argv);
        
        range = argRange.getValue();
        index_file = argIndex.getValue();
        out_filename = argOut.getValue();
        vcf_file = argVcf.getValue();
        interval_file = argInterval.getValue();
        gc_file = argGcfile.getValue();
        bNoHeader = switchNoHeader.getValue();
        bFail = switchFail.getValue();
		chr = argChr.getValue();
		max_p = argPoverlap.getValue();
        b_dumpstat = switchDumpstat.getValue();
		b_kmeans = switchKmeans.getValue();
		b_mahalanobis = switchMahalanobis.getValue();

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
    std::vector<string> pileup_names;
    
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
	std::cerr << vec_sv.size() << " SVs loaded from input file" << std::endl;
    

    read_list(index_file, pileup_names);
    int n_sample = 0;

    int n_pileup = (int) pileup_names.size();
    int n_var = (int) vec_sv.size();
	std::vector<int> n_vars (gc.num_chr+1, 0);
    
	for(int i=0; i<(int)vec_sv.size(); ++i)
	{
		n_vars[vec_sv[i].chrnum] ++;
	}

    int n_start = 0;
    int n_end = n_var-1;
    
	if (chr>0)
	{
		std::cerr << "Processing chromosome " << chr << " only" << std::endl;
	}
    std::cerr << n_pileup << " pileup files identified from index file" << std::endl;
    
    DataReader reader;
    std::vector<SampleStat> stats;

	if (chr>0)
	{
		for(int i=0; i<(int)pileup_names.size(); ++i)
		{
			pileup_names[i] += ".chr" + std::to_string(chr);
		}
	}
	
    n_sample = reader.load(pileup_names, stats, gc, chr);

    std::cerr << n_sample << " samples identified from pileup files" << std::endl;

    // TEMPORARY!
	/*
    if (chr>0)
        reader.adjust_gc_factor(gc, stats, chr);
		*/
/*
    printf("SAMPLE\tAVG_DP\tSTD_DP\tAVG_ISIZE\tSTD_ISIZE\n");
    for(int i=0 ;i<n_sample; ++i)
    {
        printf("%s\t%.2f\t%.2f\t%.2f\t%.2f\n", reader.sample_ids[i].c_str(), stats[i].avg_dp, stats[i].std_dp, stats[i].avg_isize, stats[i].std_isize);
    }
    exit (0);
    */

    if (region != "" && range != "")
    {
        std::cerr << "Region (chr:start-end) and Range (from-to) cannot be set together" << std::endl;
        exit(1);
    }

    int r_chr = 0;
    int r_start = 0;
    int r_end = 0;
    if (range != "")
    {
        std::string::size_type sz;   // alias of size_t
        
        n_start = std::stoi(range, &sz);
        n_end = std::stoi(range.substr(sz+1));
        std::cerr << "Variants index from " << n_start << " to " << n_end << " will be genotyped" << std::endl;
    }
    else if (region != "")
    {
        std::vector<std::string> C;
        split(region.c_str(), ":", C);
        r_chr = std::stoi(C[0]);
        if (r_chr < 1 || r_chr > gc.num_chr)
        {
            std::cerr << "Cannot parse region " << region << std::endl;
            exit(1);
        }
        std::vector<std::string> P;
        split(C[1].c_str(), "-", P);
        r_start = std::stoi(P[0]);
        r_end = std::stoi(P[1]);
        if (r_end <= r_start || r_start < 1 || r_end < 1)
        {
            std::cerr << "Cannot parse region " << region << std::endl;
            exit(1);
        }
    }


	OutVcf out_vcf;
	out_vcf.open(out_filename);
	if (!bNoHeader)
	{
		out_vcf.write_header(reader.sample_ids);
	}
    
	int vec_offset = 0;
    
	if (chr > 0)
	{
        while(vec_sv[vec_offset].chrnum < chr)
        {
            vec_offset++;
        }
        if (n_start == 0)
            n_start += vec_offset;
        if (n_end > vec_offset + n_vars[chr - 1])
            n_end = vec_offset + n_vars[chr] - 1;
        std::cerr << n_vars[chr] - 1 << " variants from " << n_start << ", ";
        vec_sv[n_start].print(stderr);
        std::cerr << " to " << n_end << ", ";
        vec_sv[n_end].print(stderr);
        std::cerr << " identified in the chromosome." << std::endl;
	}

    int del_count = 0;
    int dup_count = 0;
    int inv_count = 0;
    double min_GC = 0.225;
    double max_GC = 0.7;

    SvGeno G(n_sample);
    SvData D(n_sample);
    
    for(int i=n_start; i<=n_end; ++i)
    {
        // GC content of the candidate variant region.
        double sv_gc = gc.get_gc_content(vec_sv[i].chrnum, vec_sv[i].pos, vec_sv[i].end);

        // chr X and Y calling not supported yet
        if (((chr== 0 && vec_sv[i].chrnum < 23) || (chr>0 && vec_sv[i].chrnum == chr)) && (r_chr == 0 || (vec_sv[i].chrnum == r_chr && vec_sv[i].pos >= r_start && vec_sv[i].pos < r_end)) && !in_centrome(vec_sv[i]) && sv_gc > min_GC && sv_gc < max_GC)
        {

            if (vec_sv[i].svtype != INV || vec_sv[i].len > 100)
            {
                Genotyper gtyper;
                G.reset();
                D.reset();
                
                G.MAX_P_OVERLAP = max_p;
                if (vec_sv[i].svtype == DUP || vec_sv[i].svtype == CNV)
                {
                    G.MAX_P_OVERLAP *= 1.5;
                }

                reader.read_var_depth(i - vec_offset, D.dps[2]); // TODO: make read_var_depth to check whether first argument is in range

                // TODO: This is arbitrary threshold to filter out centromere region, add more systematic way to filter out problematic regions, by checking within-sample variance of regions
                if (average(D.dps[2]) < 150)
                {
                    reader.read_pair_split(vec_sv[i], D.rdstats, gc, D.all_rps, D.all_lclips, D.all_rclips);
                    if (vec_sv[i].svtype == DEL || vec_sv[i].svtype == DUP || vec_sv[i].svtype == CNV)
                    {
                        // var_depth gets GC-correction here
                        D.multi_dp = reader.read_depth100(vec_sv[i], D.dps, gc, b_dumpstat);
                    }
                    for(int j=0; j<n_sample; ++j)
                    {
                        for(int k=0; k<(int)D.dps.size(); ++k)
                        {
                            D.dps[k][j] /= (double)stats[j].avg_dp;
                        }
                    }

                    gtyper.call(vec_sv[i], D, G, b_kmeans, b_mahalanobis, stats);

                    G.info = "var" + std::to_string(i);

                    if (bFail || G.b_pass)
                    {
                        out_vcf.write_sv(vec_sv[i], D, G);

                        if (vec_sv[i].svtype == DEL)
                        {
                            del_count ++;
                        }
                        else if (vec_sv[i].svtype == DUP || vec_sv[i].svtype == CNV )
                        {
                            dup_count ++;
                        }
                        else if (vec_sv[i].svtype == INV)
                        {
                            inv_count ++;
                        }
                    }
                }
            }
        }
    }
    std::cerr<< del_count << " deletions, " << dup_count << " duplications, and " << inv_count << " inversions written to the output file" << std::endl;
	out_vcf.close();
    
    return 0;
    
}
