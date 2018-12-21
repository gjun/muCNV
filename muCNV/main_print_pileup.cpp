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
    std::string index_file;
    std::string vcf_file;
    std::string interval_file;
    std::string gc_file;
    std::string sampID;
    std::string region;
	int chr;
    
    std::vector<std::string> sample_ids;
    
    // Parsing command-line arguments
    try
    {
        TCLAP::CmdLine cmd("Command description message", ' ', "0.06");
        

        TCLAP::ValueArg<std::string> argSampleID("s","sample","Sample ID",false,"","string");
        TCLAP::ValueArg<std::string> argVcf("v","vcf","VCF file containing candidate SVs",false,"","string");
        TCLAP::ValueArg<std::string> argInterval("V","interVal", "Binary interval file containing candidate SVs", false, "", "string");
        TCLAP::ValueArg<std::string> argGcfile("f","gcFile","File containing GC content information",false, "GRCh38.gc", "string");
        TCLAP::ValueArg<std::string> argRegion("r", "region", "Genomic region (chr:start-end)", false, "", "string" );
        TCLAP::ValueArg<int> argChr("c", "chr", "Pileup contains single chromosome", false, 0, "integer" );
        
        cmd.add(argVcf);
        cmd.add(argInterval);
        cmd.add(argGcfile);
        cmd.add(argSampleID);
        cmd.add(argRegion);
		cmd.add(argChr);
        
        cmd.parse(argc, argv);
        
        sampID = argSampleID.getValue();
        vcf_file = argVcf.getValue();
        interval_file = argInterval.getValue();
        gc_file = argGcfile.getValue();
        region = argRegion.getValue();
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
    std::string pileup_name = sampID + ".pileup";
    std::string varfile_name = sampID + ".var";
    std::string idxfile_name = sampID + ".idx";

    if (chr>0)
    {
        pileup_name = sampID + ".chr" + std::to_string(chr) + ".pileup";
        varfile_name = sampID +".chr" + std::to_string(chr) +  ".var";
        idxfile_name = sampID +".chr" + std::to_string(chr) +  ".idx";
    }
    
    // TODO: make this also work with VCF file
    // read out and print pileup info
    read_svs_from_intfile(interval_file, vec_bp, vec_sv);
    
    Pileup pup;
    Pileup var_file;
    BaseFile idx_file;
    
    pup.open(pileup_name, std::ios::in | std::ios::binary);
    var_file.open(varfile_name, std::ios::in | std::ios::binary);
    idx_file.open(idxfile_name, std::ios::in | std::ios::binary);
    
    uint64_t curr_idx = 0;
    
    pup.read_int32(n_sample);
    printf("n_sample(pileup) : %d \n", n_sample);

    printf("SampleID:");
    for(int i=0; i<n_sample; ++i)
    {
        char buf[256];
        pup.read_sample_id(buf);
        printf("\t%s", buf);
    }
    printf("\n");
    
    for(int i=0; i<n_sample; ++i)
    {
        SampleStat s;
        pup.read_sample_stat(s);
        printf("Sample %d, AVG DP: %f, STdev: %f, AVG ISIZE: %f, STdev: %f \n", i, s.avg_dp, s.std_dp, s.avg_isize, s.std_isize);
    }
    
    GcContent gc;
    gc.initialize(gc_file);
    
    for(int i=0; i<n_sample; ++i)
    {
        printf("GC-factors for sample %d:\n", i);
        std::vector<double> gc_factor (gc.num_bin);
        pup.read_gc_factor(gc_factor, gc.num_bin);
        for(int j=0; j<gc.num_bin; ++j)
        {
            printf("GC-bin %d: %f\n", j, gc_factor[j]);
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
        int N = ceil((double)gc.chr_size[chr] / 100.0) + 1;
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
			int N = ceil((double)gc.chr_size[c] / 100.0) + 1;
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
        printf("Sample %d, average DP100: %d\n", i, (int)round((double)dpsum[i]/n_dp[i]/32.0));
    }
    
	if (chr>0)
	{
		int N = ceil((double)gc.chr_size[chr] / 10000.0) ;
		
		for(int j=1;j<=N;++j)
		{
			idx_file.read_uint64(curr_idx);
			printf("index position %d, tellg position %d\n", (int)curr_idx, (int)pup.tellg());
			
			for(int i=0; i<n_sample; ++i)
			{
				uint32_t n_rp = 0;
				pup.read_uint32(n_rp);
				printf("Sample %d, %d readpairs\n", i, n_rp);
				for(int k=0; k<(int)n_rp; ++k)
				{
					readpair rp;
					pup.read_readpair(rp);
					printf("\t%d\t%d\t%d\t%u\t%d\n", rp.chrnum, rp.selfpos, rp.matepos, rp.matequal, rp.pairstr);
				}
				
				uint32_t n_sp = 0;
				pup.read_uint32(n_sp);
				printf("Sample %d, %d split reads\n", i, n_sp);
				for(int k=0; k<(int)n_sp; ++k)
				{
					splitread sp;
					pup.read_splitread(sp);
					printf("\t%d\t%d\t%d\t%d\t%d\n", sp.chrnum, sp.pos, sp.sapos, sp.firstclip, sp.secondclip);

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
				
				for(int i=0; i<n_sample; ++i)
				{
					uint32_t n_rp = 0;
					pup.read_uint32(n_rp);
					printf("Sample %d, %d readpairs\n", i, n_rp);
					for(int k=0; k<(int)n_rp; ++k)
					{
						readpair rp;
						pup.read_readpair(rp);
						printf("\t%d\t%d\t%d\t%u\t%d\n", rp.chrnum, rp.selfpos, rp.matepos, rp.matequal, rp.pairstr);
					}
					
					uint32_t n_sp = 0;
					pup.read_uint32(n_sp);
					printf("Sample %d, %d split reads\n", i, n_sp);
					for(int k=0; k<(int)n_sp; ++k)
					{
						splitread sp;
						pup.read_splitread(sp);
						printf("\t%d\t%d\t%d\t%d\t%d\n", sp.chrnum, sp.pos, sp.sapos, sp.firstclip, sp.secondclip);

					}

				}
			}
		}
	}

    pup.close();
    idx_file.close();
    
    int n_var = 0;
    
    var_file.read_int32(n_sample);
    var_file.read_int32(n_var);
    
    printf("Variant File, n_sample: %d, n_var : %d\n", n_sample, n_var);
    
    printf("Sample ID(s):");
    for(int i=0; i<n_sample; ++i)
    {
        char buf[256];

        var_file.read_sample_id(buf);
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
            printf("\t%f", (dp/32.0));
        }
        printf("\n");
    }
    var_file.close();
    
    return 0;
}
