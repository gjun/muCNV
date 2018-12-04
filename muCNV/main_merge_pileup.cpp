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
    std::string gc_file;
    
    // Parsing command-line arguments
    try
    {
        TCLAP::CmdLine cmd("Command description message", ' ', "0.06");
    
        TCLAP::ValueArg<std::string> argIndex("i","index","Text file containing list of pileup samples",true,"","string");
        TCLAP::ValueArg<std::string> argOut("o","output","Output base filename for merged pileup",true,"","string");
        TCLAP::ValueArg<std::string> argGcfile("f","gcFile","File containing GC content information",false, "GRCh38.gc", "string");
        
        cmd.add(argGcfile);
        cmd.add(argIndex);
        cmd.add(argOut);
        
        cmd.parse(argc, argv);
        
        index_file = argIndex.getValue();
        output_name = argOut.getValue();
        gc_file = argGcfile.getValue();
        
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
        
        int n = 0;
        
        pileups[i].read_int32(n); // number of samples (should be all 1)
        
        if (n!=1)
        {
            std::cerr << "Wrong number of samples in " << pileup_name << " : " << n << std::endl;
            exit(1);
        }
    }
    printf("n_sample(pileup) : %d \n", n_sample);

    Pileup mpup;
    BaseFile midx;
    size_t curr_pos = 0;
    std::string filename;
    
    filename = output_name + ".pileup";
    mpup.open(filename, std::ios::out | std::ios::binary);
    filename = output_name + ".idx";
    midx.open(filename, std::ios::out | std::ios::binary);
    
    curr_pos += mpup.write_int32(n_sample);
    
    for(int i=0; i<n_sample; ++i)
    {
        char buf[256];
        
        pileups[i].read_sample_id(buf);
        printf("sample ID(pileup) : %s\n", buf);
        
        std::string sample_id(buf);
        
        curr_pos += mpup.write_sample_id(sample_id);
    }
    
    for(int i=0; i<n_sample; ++i)
    {
        SampleStat stat;
        pileups[i].read_sample_stat(stat);
        printf("Sample %d, AVG DP: %f, STdev: %f, AVG ISIZE: %f, STdev: %f \n", i, stat.avg_dp, stat.std_dp, stat.avg_isize, stat.std_isize);

        curr_pos += mpup.write_sample_stat(stat);
    }
    
    for(int i=0; i<n_sample; ++i)
    {
        std::vector<double> gc_factor;
        pileups[i].read_gc_factor(gc_factor, gc.num_bin);
        curr_pos += mpup.write_gc_factor(gc_factor, gc.num_bin);
    }

    // Write First Index Position
    midx.write_uint64(curr_pos);
    
    for(int c=1; c<=gc.num_chr; ++c)
    {
        int N = ceil((double)gc.chr_size[c] / 100.0) + 1;
        std::vector<uint16_t*> dp100 (n_sample, NULL);
        
        for(int i=0; i<n_sample; ++i)
        {
            dp100[i] = (uint16_t *) calloc(N, sizeof(uint16_t));
            pileups[i].read_depth(dp100[i], N);
        }
        for(int j=0; j<N; ++j)
        {
            for(int i=0; i<n_sample; ++i)
            {
                curr_pos += mpup.write_depth(&(dp100[i][j]), 1);
            }
        }
        for(int i=0;i<n_sample;++i)
        {
            delete [] dp100[i];
        }
    }
    
    for(int c=1;c<=gc.num_chr; ++c)
    {
        int N = ceil((double)gc.chr_size[c] / 10000.0) ;
        
        for(int j=1;j<=N;++j)
        {

            midx.write_uint64(curr_pos); // where each 10,000-bp interval starts;
        
            for(int i=0; i<n_sample; ++i)
            {
                uint32_t n_rp = 0;
                readpair rp;
                
                pileups[i].read_uint32(n_rp);
                curr_pos += mpup.write_uint32(n_rp);
                
                for(int k=0; k<n_rp; ++k)
                {
                    pileups[i].read_readpair(rp);
                    curr_pos += mpup.write_readpair(rp);
                }
                
                uint32_t n_sp = 0;
                splitread sp;

                pileups[i].read_uint32(n_sp);
                curr_pos += mpup.write_uint32(n_sp);
                
                for(int k=0; k<n_sp; ++k)
                {
                    pileups[i].read_splitread(sp);
                    curr_pos += mpup.write_splitread(sp);
                }
            }

        }
    }
    for(int i=0; i<n_sample; ++i)
    {
        pileups[i].close();
    }
    mpup.close();
    midx.close();
    
    std::vector<Pileup> varfiles (n_sample);
    
    // TODO: change Pileup class to file-based class, create Varfile class, inherited from the same common file class
    int n_var = 0;
    
    for(int i=0;i<n_sample; ++i)
    {
        std::string var_name = samples[i] + ".var";
        varfiles[i].open(var_name, std::ios::in | std::ios::binary);
        
        int n = 0;
        
        varfiles[i].read_int32(n); // number of samples (should be all 1)
        if (n!=1)
        {
            std::cerr << "Wrong number of samples in " << var_name << " : " << n << std::endl;
            exit(1);
        }
        varfiles[i].read_int32(n);
        if (n_var == 0)
        {
            if (n==0)
            {
                std::cerr << "Error: n_var is zero in " << var_name << std::endl;
                exit(1);
            }
            n_var = n;
        }
        else
        {
            if (n_var != n)
            {
                std::cerr << "Error: n_var is " << n << " in " << var_name << " while previous n_var is " << n_var << std::endl;
                exit(1);
            }
        }
    }
    std::cerr << "Var files are initialized" << std::endl;

    filename = output_name + ".var";
    Pileup mvar;
    mvar.open(filename, std::ios::out | std::ios::binary);
    
    mvar.write_int32(n_sample);
    mvar.write_int32(n_var);
    
    printf("n_var : %d\n", n_var);
    
    for(int i=0; i<n_sample; ++i)
    {
        char buf[256];
        varfiles[i].read_sample_id(buf);
        std::string sample_id(buf);
        mvar.write_sample_id(sample_id);
    }
    
    std::vector<uint16_t*> dp_var (n_sample, NULL);
    
    for(int i=0; i<n_sample; ++i)
    {
        dp_var[i] = (uint16_t *) calloc(n_var, sizeof(uint16_t));
        varfiles[i].read_depth(dp_var[i], n_var);
    }
    
    for(int j=0; j<n_var; ++j)
    {
        for(int i=0; i<n_sample; ++i)
        {
            mvar.write_depth(&(dp_var[i][j]), 1);
        }
    }
    for(int i=0;i<n_sample;++i)
    {
        delete [] dp_var[i];
    }
    
    for(int i=0; i<n_sample; ++i)
    {
        varfiles[i].close();
    }
    mvar.close();

    return 0;
}
