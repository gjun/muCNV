//
//  main_merge_pileup.cpp
//  muCNV
//
//  Created by Goo Jun on 11/24/18.
//  Copyright © 2018 Goo Jun. All rights reserved.
//

#include <stdio.h>

#include "muCNV.h"
#include "pileup.h"

// TCLAP headers
#include "tclap/CmdLine.h"
#include "tclap/Arg.h"


int main_merge_pileup(int argc, char** argv)
{
    string index_file;
    string output_name;
    string gc_file;
    
    // Parsing command-line arguments
    try
    {
        TCLAP::CmdLine cmd("Command description message", ' ', "0.06");
    
        TCLAP::ValueArg<string> argIndex("i","index","Text file containing list of pileup samples",true,"","string");
        TCLAP::ValueArg<string> argOut("o","output","Output base filename for merged pileup",true,"","string");
        TCLAP::ValueArg<string> argGcfile("f","gcFile","File containing GC content information",false, "GRCh38.gc", "string");
        
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
    
    std::vector<string> samples;
    read_list(index_file, samples);
    
    n_sample = (int) samples.size();
    
    GcContent gc;
    gc.initialize(gc_file);
    std::cerr << "GC content initialized" << std::endl;
    
    std::vector<std::ifstream> pileup_files;
    pileup_files.resize(n_sample);
    std::vector<Pileup> pileups (n_sample, Pileup(gc));
    std::cerr << "Pileup files are initialized" << std::endl;
    
    // TODO: change Pileup class to file-based class, create Varfile class, inherited from the same common file class
    
    for(int i=0;i<n_sample; ++i)
    {
        string pileup_name = samples[i] + ".pileup";
        
        pileup_files[i].open(pileup_name.c_str(), std::ios::in | std::ios::binary);
        
        int n = 0;
        
        pileups[i].read_number(pileup_files[i], n); // number of samples (should be all 1)
        
        if (n!=1)
        {
            std::cerr << "Wrong number of samples in " << pileup_name << " : " << n << std::endl;
            exit(1);
        }
    }
    printf("n_sample(pileup) : %d \n", n_sample);

    Pileup pup (gc);
    size_t curr_pos = 0;
    string filename;
    
    filename = output_name + ".mpileup";
    std::ofstream mpileup_file(filename.c_str(), std::ios::out | std::ios::binary);
    filename = output_name + ".midx";
    std::ofstream midx_file(filename.c_str(), std::ios::out | std::ios::binary);
    
    curr_pos += pup.write_number(mpileup_file, n_sample);
    
    for(int i=0; i<n_sample; ++i)
    {
        char buf[256];
        
        pileups[i].read_sample_id(pileup_files[i], buf);
        printf("sample ID(pileup) : %s\n", buf);
        
        std::string sample_id(buf);
        
        curr_pos += pup.write_sample_id(mpileup_file, sample_id);
    }
    
    for(int i=0; i<n_sample; ++i)
    {
        pileups[i].read_sample_stat(pileup_files[i], pileups[i].stat);
        printf("Sample %d, AVG DP: %f, STdev: %f, AVG ISIZE: %f, STdev: %f \n", i, pileups[i].stat.avg_dp, pileups[i].stat.std_dp, pileups[i].stat.avg_isize, pileups[i].stat.std_isize);

        curr_pos += pup.write_sample_stat(mpileup_file, pileups[i].stat);
    }
    
    for(int i=0; i<n_sample; ++i)
    {
        pileups[i].read_gc_factor(pileup_files[i], pileups[i].gc_factor);
        curr_pos += pup.write_gc_factor(mpileup_file, pileups[i].gc_factor);
    }

    // Write First Index Position
    midx_file.write(reinterpret_cast<char*>(&curr_pos), sizeof(size_t));
    
    for(int c=1; c<=gc.num_chr; ++c)
    {
        int N = ceil((double)gc.chr_size[c] / 100.0) + 1;
        std::vector<uint16_t*> dp100 (n_sample, NULL);
        
        for(int i=0; i<n_sample; ++i)
        {
            dp100[i] = (uint16_t *) calloc(N, sizeof(uint16_t));
            pileups[i].read_depth(pileup_files[i], dp100[i], N);
        }
        for(int j=0; j<N; ++j)
        {
            for(int i=0; i<n_sample; ++i)
            {
                curr_pos += pup.write_depth(mpileup_file, &(dp100[i][j]), 1);
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

            midx_file.write(reinterpret_cast<char*>(&curr_pos), sizeof(size_t)); // where each 10,000-bp interval starts;
        
            for(int i=0; i<n_sample; ++i)
            {
                uint32_t n_rp = 0;
                readpair rp;
                
                pileup_files[i].read(reinterpret_cast<char*>(&n_rp), sizeof(uint32_t));
                mpileup_file.write(reinterpret_cast<char*>(&n_rp), sizeof(uint32_t));
                curr_pos += sizeof(uint32_t);
                
                for(int k=0; k<n_rp; ++k)
                {
                    pileups[i].read_readpair(pileup_files[i], rp);
                    curr_pos += pup.write_readpair(mpileup_file, rp);
                }
                
                uint32_t n_sp = 0;
                splitread sp;

                pileup_files[i].read(reinterpret_cast<char*>(&n_sp), sizeof(uint32_t));
                mpileup_file.write(reinterpret_cast<char*>(&n_sp), sizeof(uint32_t));
                curr_pos += sizeof(uint32_t);
                
                for(int k=0; k<n_sp; ++k)
                {
                    pileups[i].read_splitread(pileup_files[i], sp);
                    curr_pos += pup.write_splitread(mpileup_file, sp);
                }
            }

        }
    }
    for(int i=0; i<n_sample; ++i)
    {
        pileup_files[i].close();
    }
    mpileup_file.close();
    midx_file.close();
    
    
    std::vector<std::ifstream> var_files;
    var_files.resize(n_sample);
    
    std::cerr << "Pileup files are initialized" << std::endl;
    
    // TODO: change Pileup class to file-based class, create Varfile class, inherited from the same common file class
    int n_var = 0;
    
    for(int i=0;i<n_sample; ++i)
    {
        string var_name = samples[i] + ".var";
        var_files[i].open(var_name.c_str(), std::ios::in | std::ios::binary);
        
        int n = 0;
        
        pileups[i].read_number(var_files[i], n); // number of samples (should be all 1)
        if (n!=1)
        {
            std::cerr << "Wrong number of samples in " << var_name << " : " << n << std::endl;
            exit(1);
        }
        pileups[i].read_number(var_files[i], n);
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
    filename = output_name + ".mvar";
    std::ofstream mvar_file(filename.c_str(), std::ios::out | std::ios::binary);
    
    pup.write_number(mvar_file, n_sample);
    pup.write_number(mvar_file, n_var);
    
    printf("n_var : %d\n", n_var);
    
    
    for(int i=0; i<n_sample; ++i)
    {
        char buf[256];
        pileups[i].read_sample_id(var_files[i], buf);
        std::string sample_id(buf);
        pup.write_sample_id(mvar_file, sample_id);
    }
    
    std::vector<uint16_t*> dp_var (n_sample, NULL);
    
    for(int i=0; i<n_sample; ++i)
    {
        dp_var[i] = (uint16_t *) calloc(n_var, sizeof(uint16_t));
        pileups[i].read_depth(var_files[i], dp_var[i], n_var);
    }
    
    for(int j=0; j<n_var; ++j)
    {
        for(int i=0; i<n_sample; ++i)
        {
            pup.write_depth(mvar_file, &(dp_var[i][j]), 1);
        }
    }
    for(int i=0;i<n_sample;++i)
    {
        delete [] dp_var[i];
    }
    
    for(int i=0; i<n_sample; ++i)
    {
        var_files[i].close();
    }
    mvar_file.close();

    return 0;
}
