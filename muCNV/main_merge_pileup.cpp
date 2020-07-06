//
//  main_merge_pileup.cpp
//  muCNV
//
//  Created by Goo Jun on 11/24/18.
//  Copyright © 2018 Goo Jun. All rights reserved.
//

#include <stdio.h>

#include "in_vcf.h"
#include "pileup.h"
#include <cmath>

// TCLAP headers
#include "tclap/CmdLine.h"
#include "tclap/Arg.h"


int main_merge_pileup(int argc, char** argv)
{
    std::string index_file;
    std::string output_name;
	std::string interval_file;
    std::string gc_file;
    
    // Parsing command-line arguments
    try
    {
        TCLAP::CmdLine cmd("Command description message", ' ', "0.06");
    
        TCLAP::ValueArg<std::string> argIndex("i","index","Text file containing list of pileup samples",true,"","string");
        TCLAP::ValueArg<std::string> argOut("o","output","Output base filename for merged pileup",true,"","string");
        TCLAP::ValueArg<std::string> argGcfile("f","gcFile","File containing GC content information",true, "GRCh38.gc", "string");
		TCLAP::ValueArg<std::string> argInterval("V","interVal", "Binary interval file containing candidate SVs", true, "", "string");
        
        cmd.add(argGcfile);
        cmd.add(argIndex);
        cmd.add(argOut);
		cmd.add(argInterval);
    
        cmd.parse(argc, argv);
        
        index_file = argIndex.getValue();
        output_name = argOut.getValue();
        gc_file = argGcfile.getValue();
		interval_file = argInterval.getValue();
    }
    catch (TCLAP::ArgException &e)
    {
        std::cerr << "Error: " << e.error() << " for arg " << e.argId() << std::endl;
        abort();
    }
    
    int n_sample = 0;
    
    std::vector<std::string> samples;
    read_list(index_file, samples);
    
    n_sample = (int) samples.size();
    
    GcContent gc;
    gc.initialize(gc_file);
    std::cerr << "GC content initialized" << std::endl;

    std::vector<Pileup> pileups(n_sample);
    std::cerr << "Pileup files are initialized" << std::endl;
    
    for(int i=0;i<n_sample; ++i)
    {
        std::string pileup_name = samples[i] + ".pileup";
        
        pileups[i].open(pileup_name, std::ios::in | std::ios::binary);
        
        if (!pileups[i].good() )
        {
            std::cerr << "Error: Cannot open " << pileup_name << std::endl;
            exit(1);
        }
        int n = 0;
        
        pileups[i].read_int32(n); // number of samples (should be all 1)
        
        if (n!=1)
        {
            std::cerr << "Error: Wrong number of samples in " << pileup_name << " : " << n << std::endl;
            exit(1);
        }
    }
    printf("n_sample(pileup) : %d \n", n_sample);

    std::vector<Pileup> mpups(gc.num_chr+1);
    std::vector<BaseFile> midxs(gc.num_chr+1);
    std::vector<size_t> curr_pos(gc.num_chr+1, 0);

    std::string filename;
    
	// Currently this 'merging' supports only merging from single-sample pilieups
	// It does not support merging already merged pileups

	// Currently merged pileups are automatically separaged by chromosomes

	for(int chr=1; chr<=gc.num_chr; ++chr)
	{
		filename = output_name + ".chr" + std::to_string(chr) + ".pileup";
		mpups[chr].open(filename, std::ios::out | std::ios::binary);
        if (!mpups[chr].good() )
        {
            std::cerr << "Error: Cannot open " << filename << " to write " << std::endl;
            exit(1);
        }
		filename = output_name + ".chr" + std::to_string(chr) +  ".idx";
		midxs[chr].open(filename, std::ios::out | std::ios::binary);
        if (!midxs[chr].good() )
        {
            std::cerr << "Error: Cannot open " << filename << " to write " << std::endl;
            exit(1);
        }
		curr_pos[chr] += mpups[chr].write_int32(n_sample);
	}
		
	for(int i=0; i<n_sample; ++i)
	{
		char buf[256];
		
		pileups[i].read_sample_id(buf);
		std::string sample_id(buf);

		
		for(int chr=1; chr<=gc.num_chr; ++chr)
		{
			curr_pos[chr] += mpups[chr].write_sample_id(sample_id);
		}
	}
		
	for(int i=0; i<n_sample; ++i)
	{
		SampleStat stat;
		pileups[i].read_sample_stat(stat);
		printf("Sample %d, AVG DP: %f, Stdev: %f, AVG ISIZE: %f, Stdev: %f \n", i, stat.avg_dp, stat.std_dp, stat.avg_isize, stat.std_isize);

		for(int chr=1; chr<=gc.num_chr; ++chr)
		{
			curr_pos[chr] += mpups[chr].write_sample_stat(stat);
		}
	}
		
		// Let's just write GC factor for everything...
	for(int i=0; i<n_sample; ++i)
	{
		std::vector<double> gc_factor;
		pileups[i].read_gc_factor(gc_factor, gc.num_bin);

		for(int chr=1; chr<=gc.num_chr; ++chr)
		{
			curr_pos[chr] += mpups[chr].write_gc_factor(gc_factor, gc.num_bin);
		}
	}

	// Write First Index Position
	for(int chr=1; chr<=gc.num_chr; ++chr)
	{
		midxs[chr].write_uint64(curr_pos[chr]);
	}
		

	for(int chr=1; chr<=gc.num_chr; ++chr)
	{
		std::vector<uint16_t*> dp100 (n_sample, NULL);
			
		int N = gc.n_interval[chr];

		for(int i=0; i<n_sample; ++i)
		{
			dp100[i] = (uint16_t *) calloc(N, sizeof(uint16_t));
			pileups[i].read_depth(dp100[i], N);
		}
		for(int j=0; j<N; ++j)
		{
			for(int i=0; i<n_sample; ++i)
			{
				curr_pos[chr] += mpups[chr].write_depth(&(dp100[i][j]), 1);
			}
		}
		for(int i=0;i<n_sample;++i)
		{
			delete [] dp100[i];
		}
	}
		

	for(int chr=1;chr<=gc.num_chr; ++chr)
	{
		int N = ceil((double)gc.chr_size[chr] / 10000.0) ;
		
		for(int j=1;j<=N;++j)
		{

			midxs[chr].write_uint64(curr_pos[chr]); // where each 10,000-bp interval starts;
		
			for(int i=0; i<n_sample; ++i)
			{
				uint32_t n_rp = 0;
				readpair rp;
				
				pileups[i].read_uint32(n_rp);
				curr_pos[chr] += mpups[chr].write_uint32(n_rp);
				
				for(int k=0; k<(int)n_rp; ++k)
				{
					pileups[i].read_readpair(rp);
					curr_pos[chr] += mpups[chr].write_readpair(rp);
				}
				
				uint32_t n_sp = 0;
				splitread sp;

				pileups[i].read_uint32(n_sp);
				curr_pos[chr] += mpups[chr].write_uint32(n_sp);
				
				for(int k=0; k<(int)n_sp; ++k)
				{
					pileups[i].read_splitread(sp);
					curr_pos[chr] += mpups[chr].write_splitread(sp);
				}
                
                uint32_t n_lclip = 0;
                sclip myclip;
                
                pileups[i].read_uint32(n_lclip);
                curr_pos[chr] += mpups[chr].write_uint32(n_lclip);
                
                for(int k=0; k<(int)n_lclip; ++k)
                {
                    pileups[i].read_softclip(myclip);
                    curr_pos[chr] += mpups[chr].write_softclip(myclip);
                }
                
                uint32_t n_rclip = 0;
                
                pileups[i].read_uint32(n_rclip);
                curr_pos[chr] += mpups[chr].write_uint32(n_rclip);
                
                for(int k=0; k<(int)n_rclip; ++k)
                {
                    pileups[i].read_softclip(myclip);
                    curr_pos[chr] += mpups[chr].write_softclip(myclip);
                }
                
			}

		}
	}

	for(int i=0; i<n_sample; ++i)
	{
		pileups[i].close();
	}
	for(int chr=1;chr<=gc.num_chr; ++chr)
	{
		mpups[chr].close();
		midxs[chr].close();
	}
	
	std::vector<Pileup> varfiles (n_sample);
	
	// TODO: change Pileup class to file-based class, create Varfile class, inherited from the same common file class
	int n_var_total = 0;
	
	for(int i=0;i<n_sample; ++i)
	{
		std::string var_name = samples[i] + ".var";
		varfiles[i].open(var_name, std::ios::in | std::ios::binary);
        if (!varfiles[i].good() )
        {
            std::cerr << "Error: Cannot open " << var_name << std::endl;
            exit(1);
        }
		int n = 0;
		
		varfiles[i].read_int32(n); // number of samples (should be all 1)
		if (n!=1)
		{
			std::cerr << "Wrong number of samples in " << var_name << " : " << n << std::endl;
			exit(1);
		}
		varfiles[i].read_int32(n);
		if (n_var_total == 0)
		{
			if (n==0)
			{
				std::cerr << "Error: n_var is zero in " << var_name << std::endl;
				exit(1);
			}
			n_var_total = n;
		}
		else
		{
			if (n_var_total != n)
			{
				std::cerr << "Error: n_var is " << n << " in " << var_name << " while previous n_var is " << n_var_total << std::endl;
				exit(1);
			}
		}
	}
	std::cerr << "Var files are initialized" << std::endl;

	std::vector<int> n_vars(gc.num_chr+1, 0);
	int sum_var = 0;
	if (1)
	{
		// to release memory vec_sv, vec_bp after loading number of variants only
		std::vector<sv> vec_sv;
		std::vector<breakpoint> vec_bp;
		read_svs_from_intfile(interval_file, vec_bp, vec_sv);
		std::cerr << vec_sv.size() << " SVs identified " << std::endl;

		for(int i=0; i<(int)vec_sv.size(); ++i)
		{
			n_vars[vec_sv[i].chrnum] ++;
		}

		for (int chr=1; chr<= gc.num_chr; ++chr)
		{
			std::cerr << "Number of variants in chr " << chr << " : " << n_vars[chr] << std::endl;
			sum_var += n_vars[chr];
		}
	}
	if (sum_var != n_var_total)
	{
		std::cerr << "Error: n_var_total is " << n_var_total << " while sum_var  is " << sum_var << std::endl;
		exit(1);
	}

	std::vector<Pileup> mvars (gc.num_chr+1);
	for(int chr=1; chr<=gc.num_chr; ++chr)	
	{
		filename = output_name + ".chr" + std::to_string(chr) + ".var";
		mvars[chr].open(filename, std::ios::out | std::ios::binary);
        if (!mvars[chr].good() )
        {
            std::cerr << "Error: Cannot open " << filename << " to write" << std::endl;
            exit(1);
        }
		mvars[chr].write_int32(n_sample);
		mvars[chr].write_int32(n_vars[chr]);
	}
	
	for(int i=0; i<n_sample; ++i)
	{
		char buf[256];
		varfiles[i].read_sample_id(buf);
		std::string sample_id(buf);

		for(int chr=1; chr<=gc.num_chr; ++chr)	
		{
			mvars[chr].write_sample_id(sample_id);
		}
	}
		
	for(int chr=1; chr<=gc.num_chr; ++chr)	
	{
		std::vector<uint16_t*> dp_var (n_sample, NULL);
		
		for(int i=0; i<n_sample; ++i)
		{
			dp_var[i] = (uint16_t *) calloc(n_vars[chr], sizeof(uint16_t));
			varfiles[i].read_depth(dp_var[i], n_vars[chr]);
		}
		
		for(int j=0; j<n_vars[chr]; ++j)
		{
			for(int i=0; i<n_sample; ++i)
			{
				mvars[chr].write_depth(&(dp_var[i][j]), 1);
			}
		}

		for(int i=0;i<n_sample;++i)
		{
			delete [] dp_var[i];
		}
		mvars[chr].close();
	}

	for(int i=0; i<n_sample; ++i)
	{
		varfiles[i].close();
	}

    return 0;
}
