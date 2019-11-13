//
//  main_print_pileup.cpp
//  muCNV
//
//  Created by Goo Jun on 11/23/18.
//  Copyright © 2018 Goo Jun. All rights reserved.
//

#include <stdio.h>
#include <cmath>
// TCLAP headers
#include "tclap/CmdLine.h"
#include "tclap/Arg.h"

#include "in_vcf.h"
#include "gc_content.h"
#include "pileup.h"

int main_print_pileup(int argc, char** argv)
{
    std::string index_filename;
    std::string vcf_filename;
    std::string interval_filename;
    std::string gc_filename;
    std::string sampID;
    std::string pileupID;
    std::string out_filename;
//    std::string region;
	int chr;
    
    std::vector<std::string> sample_ids;
    
    // Parsing command-line arguments
    try
    {
        TCLAP::CmdLine cmd("Command description message", ' ', "0.06");
        
        TCLAP::ValueArg<std::string> argPileup("p","pileup","Prefix of the pileup file",false,"","string");
        TCLAP::ValueArg<std::string> argSample("s","sample","Sample ID of a single sample to be printed from the pileup",false,"","string");
        TCLAP::ValueArg<std::string> argVcf("v","vcf","VCF file containing candidate SVs",false,"","string");
        TCLAP::ValueArg<std::string> argInterval("V","interVal", "Binary interval file containing candidate SVs", false, "", "string");
        TCLAP::ValueArg<std::string> argGcfile("f","gcFile","File containing GC content information",false, "GRCh38.gc", "string");
        TCLAP::ValueArg<std::string> argOut("o","out","Output file, default: stdout",false,"","string");
//        TCLAP::ValueArg<std::string> argRegion("r", "region", "Genomic region (chr:start-end)", false, "", "string" );
        TCLAP::ValueArg<int> argChr("c", "chr", "Pileup contains single chromosome ( .chrN.pileup, as in merged pileups)", false, 0, "integer" );
        
        cmd.add(argVcf);
        cmd.add(argInterval);
        cmd.add(argGcfile);
        cmd.add(argSample);
        cmd.add(argPileup);
        cmd.add(argOut);
//        cmd.add(argRegion);
		cmd.add(argChr);
        
        cmd.parse(argc, argv);
        
        pileupID = argPileup.getValue();
        sampID = argSample.getValue();
        vcf_filename = argVcf.getValue();
        interval_filename = argInterval.getValue();
        gc_filename = argGcfile.getValue();
        out_filename = argOut.getValue();
//        region = argRegion.getValue();
		chr = argChr.getValue();
       
    }
    catch (TCLAP::ArgException &e)
    {
        std::cerr << "Error: " << e.error() << " for arg " << e.argId() << std::endl;
        abort();
    }
    
    int n_sample = 0;
    
    std::vector<sv> vec_sv;
    std::vector<breakpoint> vec_bp;
    
    // mvar, mpileup ?
    std::string pileup_name = pileupID + ".pileup";
    std::string varfile_name = pileupID + ".var";
    std::string idxfile_name = pileupID + ".idx";

    if (chr>0)
    {
        pileup_name = pileupID + ".chr" + std::to_string(chr) + ".pileup";
        varfile_name = pileupID +".chr" + std::to_string(chr) +  ".var";
        idxfile_name = pileupID +".chr" + std::to_string(chr) +  ".idx";
    }
    

    Pileup pup;
    Pileup var_file;
    BaseFile idx_file;
    
    pup.open(pileup_name, std::ios::in | std::ios::binary);

	if (!pup.good())
	{
		fprintf(stderr, "Cannot open %s \n", pileup_name.c_str());
		exit(1);
	}

    var_file.open(varfile_name, std::ios::in | std::ios::binary);

	if (!var_file.good())
	{
		fprintf(stderr, "Cannot open %s \n", varfile_name.c_str());
		exit(1);
	}
    idx_file.open(idxfile_name, std::ios::in | std::ios::binary);
    
	if (!idx_file.good())
	{
		fprintf(stderr, "Cannot open %s \n", idxfile_name.c_str());
		exit(1);
	}

    uint64_t curr_idx = 0;
    
    pup.read_int32(n_sample);
    printf("n_sample(pileup) : %d \n", n_sample);
    
    int sample_idx = -1;

    printf("SampleID:");
    for(int i=0; i<n_sample; ++i)
    {
        char buf[256];
        pup.read_sample_id(buf);
        std::string this_ID (buf);
        if (this_ID == sampID)
        {
            printf("Only sample %s at index %d will be printed \n", sampID.c_str(), i);
            sample_idx = i;
        }
        printf("\t%s", buf);
    }
    printf("\n");
    
    if (sampID != "" && sample_idx == -1)
    {
        // Single sample ID is given but not found in the pileup
        std::cerr << "Error, Sampel ID " << sampID << " cannot be found in pileup" << std::endl;
        exit(1);
    }
    
    for(int i=0; i<n_sample; ++i)
    {
        SampleStat s;
        pup.read_sample_stat(s);
        
        if (sample_idx < 0 || sample_idx == i)
        {
            printf("Sample %d, Mean DP: %f, Std DP: %f, Mean ISIZE: %f, Std ISIZE: %f \n", i, s.avg_dp, s.std_dp, s.avg_isize, s.std_isize);
        }
    }
    
    GcContent gc;
    gc.initialize_quick(gc_filename);
    
    for(int i=0; i<n_sample; ++i)
    {
        std::vector<double> gc_factor (gc.num_bin);
        pup.read_gc_factor(gc_factor, gc.num_bin);
    
        if (sample_idx < 0 || sample_idx == i)
        {
            printf("GC-factors for sample %d:\n", i);
            
            for(int j=0; j<gc.num_bin; ++j)
            {
                printf("GC-bin %d: %f\n", j, gc_factor[j]);
            }
        }
    }
	fflush(stdout);

    idx_file.read_uint64(curr_idx);
    printf("index position %d, tellg position %lu\n", (int)curr_idx, (unsigned long)pup.tellg());
    
    std::vector<uint64_t> dpsum (n_sample, 0);
    std::vector<uint64_t> n_dp (n_sample, 0);

    uint16_t *dp100 = (uint16_t*) malloc(sizeof(uint16_t) * n_sample);

	if (chr>0)
	{
        int N = gc.n_interval[chr];
        for(int j=0;j<N;++j)
        {
            pup.read_depth(dp100, n_sample);
            for(int i=0; i<n_sample; ++i)
            {
                dpsum[i] += dp100[i];
                n_dp[i] +=1;
            }
        }
	}
	else
	{
		for(int c=1; c<=gc.num_chr; ++c)
		{
			int N = gc.n_interval[c];
			for(int j=0;j<N;++j)
			{
				pup.read_depth(dp100, n_sample);
				for(int i=0; i<n_sample; ++i)
				{
					dpsum[i] += dp100[i];
					n_dp[i] +=1;
				}
			}
		}
	}
    delete [] dp100;
    for(int i=0; i<n_sample; ++i)
    {
        if (sample_idx < 0 || sample_idx == i)
        {
            printf("Sample %d, average DP100: %d\n", i, (int)round((double)dpsum[i]/n_dp[i]/32.0));
        }
    }
    
	if (chr>0)
	{
		int N = ceil((double)gc.chr_size[chr] / 10000.0) ;
		
		for(int j=1;j<=N;++j)
		{
			idx_file.read_uint64(curr_idx);
			printf("index position %d, tellg position %d\n", (int)curr_idx, (int)pup.tellg());
            printf("chr%d:%d-%d\n", chr, (j-1)*10000, j*10000);

			for(int i=0; i<n_sample; ++i)
			{
				uint32_t n_rp = 0;
				pup.read_uint32(n_rp);
                if (sample_idx < 0 || sample_idx == i)
                    printf("Sample %d, %d readpairs\n", i, n_rp);
                
				for(int k=0; k<(int)n_rp; ++k)
				{
					readpair rp;
					pup.read_readpair(rp);
                    if (sample_idx < 0 || sample_idx == i)
                    {
                        rp.chrnum = chr;

                        printf("\t%d\t%d\t%d\t%u\t%d\n", rp.chrnum, rp.selfpos, rp.matepos, rp.matequal, rp.pairstr);
                    }
				}
				
				uint32_t n_sp = 0;
				pup.read_uint32(n_sp);
                if (sample_idx < 0 || sample_idx == i)
                    printf("Sample %d, %d split reads\n", i, n_sp);
				for(int k=0; k<(int)n_sp; ++k)
				{
					splitread sp;
					pup.read_splitread(sp);
                    if (sample_idx < 0 || sample_idx == i)
                    {
                        sp.chrnum = chr;

                        printf("\t%d\t%d\t%d\t%d\t%d\n", sp.chrnum, sp.pos, sp.sapos, sp.firstclip, sp.secondclip);
                    }
				}
                
                uint32_t n_lclip = 0;
                pup.read_uint32(n_lclip);
                if (sample_idx < 0 || sample_idx == i)
                    printf("Sample %d, %d left clips\n", i, n_lclip);
                for(int k=0; k<(int)n_lclip; ++k)
                {
                    sclip myclip;
                    pup.read_softclip(myclip);
                    myclip.chrnum = chr;

                    if (sample_idx < 0 || sample_idx == i)
                        printf("\t%d\t%d\n", myclip.chrnum, myclip.pos);
                }

                uint32_t n_rclip = 0;
                pup.read_uint32(n_rclip);
                if (sample_idx < 0 || sample_idx == i)
                    printf("Sample %d, %d right clips\n", i, n_rclip);
                for(int k=0; k<(int)n_rclip; ++k)
                {
                    sclip myclip;
                    pup.read_softclip(myclip);
                    myclip.chrnum = chr;

                    if (sample_idx < 0 || sample_idx == i)
                        printf("\t%d\t%d\n", myclip.chrnum, myclip.pos);
                }

			}
		}
	}
	else
	{
		for(int c=1;c<=gc.num_chr; ++c)
		{
			int N = ceil((double)gc.chr_size[c] / 10000.0) ;
			
			for(int j=1;j<=N;++j)
			{
				idx_file.read_uint64(curr_idx);
				printf("index position %d, tellg position %d\n", (int)curr_idx, (int)pup.tellg());
                printf("chr%d:%d-%d\n", c, (j-1)*10000, j*10000);

				for(int i=0; i<n_sample; ++i)
				{
					uint32_t n_rp = 0;
					pup.read_uint32(n_rp);
                    if (sample_idx < 0 || sample_idx == i)
                        printf("Sample %d, %d readpairs\n", i, n_rp);
					for(int k=0; k<(int)n_rp; ++k)
					{
						readpair rp;
						pup.read_readpair(rp);
                        rp.chrnum = c;

                        if (sample_idx < 0 || sample_idx == i)
                            printf("\t%d\t%d\t%d\t%u\t%d\n", rp.chrnum, rp.selfpos, rp.matepos, rp.matequal, rp.pairstr);
					}
					
					uint32_t n_sp = 0;
					pup.read_uint32(n_sp);
                    if (sample_idx < 0 || sample_idx == i)
                        printf("Sample %d, %d split reads\n", i, n_sp);
					for(int k=0; k<(int)n_sp; ++k)
					{
						splitread sp;
						pup.read_splitread(sp);
                        sp.chrnum = c;
   
                        if (sample_idx < 0 || sample_idx == i)
                            printf("\t%d\t%d\t%d\t%d\t%d\n", sp.chrnum, sp.pos, sp.sapos, sp.firstclip, sp.secondclip);
					}
                    
                    uint32_t n_lclip = 0;
                    pup.read_uint32(n_lclip);
                    if (sample_idx < 0 || sample_idx == i)
                        printf("Sample %d, %d left clips\n", i, n_lclip);
                    for(int k=0; k<(int)n_lclip; ++k)
                    {
                        sclip myclip;
                        pup.read_softclip(myclip);
                        myclip.chrnum = c;
                        if (sample_idx < 0 || sample_idx == i)
                            printf("\t%d\t%d\n", myclip.chrnum, myclip.pos);
                    }
                    
                    uint32_t n_rclip = 0;
                    pup.read_uint32(n_rclip);
                    if (sample_idx < 0 || sample_idx == i)
                        printf("Sample %d, %d right clips\n", i, n_rclip);
                    for(int k=0; k<(int)n_rclip; ++k)
                    {
                        sclip myclip;
                        pup.read_softclip(myclip);
                        myclip.chrnum = c;
                        if (sample_idx < 0 || sample_idx == i)
                            printf("\t%d\t%d\n", myclip.chrnum, myclip.pos);
                    }

				}
			}
		}
	}

    pup.close();
    idx_file.close();
    
    // VAR file
    if (vcf_filename != "" && interval_filename == "")
    {
        read_svs_from_vcf(vcf_filename, vec_bp, vec_sv);
    }
    else if (vcf_filename == "" && interval_filename != "")
    {
        read_svs_from_intfile(interval_filename, vec_bp, vec_sv);
    }
    else
    {
        std::cerr << "Error, VCF file or Interval file is required (but not both)" << std::endl;
        exit (1);
    }
    
    
    int n_var = 0;
    
    var_file.read_int32(n_sample);
    var_file.read_int32(n_var);
    
    printf("Variant File, n_sample: %d, n_var : %d\n", n_sample, n_var);
    
    printf("Sample ID(s):");
    for(int i=0; i<n_sample; ++i)
    {
        char buf[256];

        var_file.read_sample_id(buf);
        if (sample_idx < 0 || sample_idx == i)
            printf("\t%s", buf);
    }
    printf("\n");
    
	int vec_offset = 0;
	if (chr>0)
	{
		while(vec_sv[vec_offset].chrnum < chr)
		{
			vec_offset++;
		}
		std::cerr << n_var << " variants from " << vec_offset  << std::endl;
	}
    for(int j=0;j<n_var;++j)
    {
        vec_sv[vec_offset+j].print(stdout);
        for(int i=0; i<n_sample; ++i)
        {
            uint16_t dp;
            var_file.read_depth(&dp, 1);
            if (sample_idx < 0 || sample_idx == i)
                printf("\t%f", (dp/32.0));
        }
        printf("\n");
    }
    var_file.close();
    
    return 0;
}
