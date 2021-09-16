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


int main_merge_var(int argc, char** argv)
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
//        TCLAP::ValueArg<int> argChr("c", "chr", "Single chromosome", false, 0, "integer" );
        cmd.add(argGcfile);
        cmd.add(argIndex);
        cmd.add(argOut);
//        cmd.add(argChr);
		cmd.add(argInterval);
    
        cmd.parse(argc, argv);
        
        index_file = argIndex.getValue();
        output_name = argOut.getValue();
        gc_file = argGcfile.getValue();
		interval_file = argInterval.getValue();
//        chr = argChr.getValue();
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
    std::vector<int> vec_chr;
	if (1)
	{
		// to release memory vec_sv, vec_bp after loading number of variants only
		std::vector<sv> vec_sv;
		std::vector<breakpoint> vec_bp;
		read_svs_from_intfile(interval_file, vec_bp, vec_sv);
		std::cerr << vec_sv.size() << " SVs identified " << std::endl;
        int prev_chr = -1;

		for(int i=0; i<(int)vec_sv.size(); ++i)
		{
			n_vars[vec_sv[i].chrnum] ++;
            if (vec_sv[i].chrnum != prev_chr)
            {
                vec_chr.push_back(vec_sv[i].chrnum);
                prev_chr = vec_sv[i].chrnum;
            }
		}
        for(int i=0; i<(int)vec_chr.size(); ++i)
        {
            int chr = vec_chr[i];
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
	for(int i=0; i<(int)vec_chr.size(); ++i)	
	{
        int chr = vec_chr[i];
		std::string filename = output_name + ".chr" + std::to_string(chr) + ".var";
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

		for(int i=0; i<(int)vec_chr.size(); ++i)	
		{
            int chr = vec_chr[i];
			mvars[chr].write_sample_id(sample_id);
		}
	}
		
	for(int i=0;i<(int)vec_chr.size();++i)
	{
        int chr = vec_chr[i];
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
