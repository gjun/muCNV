//
//  multi_pileup.cpp
//  muCNV
//
//  Created by Goo Jun on 11/25/18.
//  Copyright © 2018 Goo Jun. All rights reserved.
//

#include "multi_pileup.h"
#include <math.h>
/*
void MultiPileup::read()
{
    string index_file;
    string vcf_file;
    string interval_file;
    string gc_file;
    string sampID;
    string region;
    
    std::vector<string> sample_ids;
    
    int n_sample = 0;
    
    std::vector<sv> vec_sv;
    std::vector<breakpoint> vec_bp;
    
    string pileup_name = sampID + ".pileup";
    string varfile_name = sampID + ".var";
    string idxfile_name = sampID + ".idx";
    
    // TODO: make this also work with VCF file
    // read out and print pileup info
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
        printf("GC-factors for sapmle %d:\n", i);
        std::vector<double> gc_factor (gc.num_bin);
        pup.read_gc_factor(gc_factor, gc.num_bin);
        for(int j=0; j<gc.num_bin; ++j)
        {
            printf("GC-bin %d: %f\n", j, gc_factor[j]);
        }
    }
    
    idx_file.read_uint64(curr_idx);
    printf("index position %d, tellg position %lu\n", (int)curr_idx, (unsigned long)pup.tellg());
    
    std::vector<uint64_t> dpsum (n_sample, 0);
    std::vector<uint64_t> n_dp (n_sample, 0);
    
    for(int c=1; c<=gc.num_chr; ++c)
    {
        int N = ceil((double)gc.chr_size[c] / 100.0) + 1;
        uint16_t dp100;
        for(int j=0;j<N;++j)
        {
            for(int i=0; i<n_sample; ++i)
            {
                pup.read_depth(&dp100, 1);
                dpsum[i] += dp100;
                n_dp[i] +=1;
            }
        }
    }
    for(int i=0; i<n_sample; ++i)
    {
        printf("Sample %d, average DP100: %d\n", i, (int)round((double)dpsum[i]/n_dp[i]/32.0));
    }
    
    
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
                for(int k=0; k<n_rp; ++k)
                {
                    readpair rp;
                    pup.read_readpair(rp);
                    printf("\t%d\t%d\t%d\t%u\t%d\n", rp.chrnum, rp.selfpos, rp.matepos, rp.matequal, rp.pairstr);
                }
                
                uint32_t n_sp = 0;
                pup.read_uint32(n_sp);
                printf("Sample %d, %d split reads\n", i, n_sp);
                for(int k=0; k<n_sp; ++k)
                {
                    splitread sp;
                    pup.read_splitread(sp);
                    printf("\t%d\t%d\t%d\t%d\t%d\n", sp.chrnum, sp.pos, sp.sapos, sp.firstclip, sp.secondclip);
                    
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
    
    for(int j=0;j<n_var;++j)
    {
        vec_sv[j].print();
        for(int i=0; i<n_sample; ++i)
        {
            uint16_t dp;
            var_file.read_depth(&dp, 1);
            printf("\t%f", (dp/32.0));
        }
        printf("\n");
    }
    var_file.close();
}

*/

/*
double MultiPileup::gcCorrected(double D, int chr, int pos)
{
    int p = pos*2 / gc.binsize;
    //    std::cerr << "pos " << pos << " p " << p << std::endl;
    int bin = gc.gc_array[chr][p];
    
    if (bin<20 && gc_factor[bin]>0.0001)
    {
        //        std::cerr << D << " at " << chr << ":" << pos << " is adjusted to " << D/gc_factor[bin] << " by gc Factor "<< gc_factor[bin] << std::endl;
        return D / gc_factor[bin];
    }
    else
    {
        // Let's not make adjustment for bins with '255' value
        return D;
    }
}
*/
