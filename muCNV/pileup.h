//
//  pileup.h
//  muCNV
//
//  Created by Goo Jun on 11/13/18.
//  Copyright © 2018 Goo Jun. All rights reserved.
//

#ifndef pileup_h
#define pileup_h

#include "sv.h"
#include "gc_content.h"
#include "debug.h"

#include <string>
#include <vector>

class readpair
{
public:
    int8_t chrnum;
    int32_t selfpos;
    int32_t matepos;
    int8_t matequal;
    int8_t pairstr;
    // pairstr: 00b, 01b, 10b, 11b
};

class splitread
{
public:
    int8_t chrnum;
    int32_t pos;
    int32_t sapos;
    int16_t firstclip; // soft clip position (+: left-side, -: right-side) in primary alignment
    int16_t secondclip; // soft clip position (+: left-side, -: right-side) in secondary alignment
};

class Pileup
{
public:
    double avg_dp;
    double std_dp;
    double avg_isize;
    double std_isize;
    double med_isize;
    double avg_rlen;
    
    GcContent& gc;
    // Get GC corrected depth for chr / pos
    double gcCorrected(double, int, int);

    std::vector<double> gc_factor;
    
    std::vector< uint16_t * > depth100; // to store depth for every 100bp interval
    std::vector<int> nbin_100;
    
    std::vector<uint64_t> gc_sum;
    std::vector<uint64_t> gc_cnt;
    
    std::vector<readpair> vec_rp;
    std::vector<splitread> vec_sp;
    
    void write(std::string &, std::vector<sv> &);
    void write_text(std::string &, std::vector<sv> &);
    
    Pileup(GcContent &x) : gc(x) {
        // Check whether gc is initialized correctly
        gc_factor.resize(gc.num_bin);
        gc_sum.resize(gc.num_bin);
        gc_cnt.resize(gc.num_bin);
        
        for(int i=0;i<gc.num_bin;++i)
        {
            gc_sum[i] = 0;
            gc_cnt[i] = 0;
        }
        
        depth100.resize(gc.num_chr + 1);
        nbin_100.resize(gc.num_chr + 1);
        nbin_100[0] = 0;
        
        for(int i=1; i<=gc.num_chr; ++i)
        {
            nbin_100[i] = ceil((double)gc.chr_size[i] / 100.0) + 1 ;
            depth100[i] = (uint16_t *) calloc(nbin_100[i], sizeof(uint16_t));
            
            DMSG("chr " << i << " bin size " << nbin_100[i]);
        }
        
        
    };
};


#endif /* pileup_h */
