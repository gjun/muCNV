//
//  bFile.h
//  muCNV
//
//  Created by Goo Jun on 11/13/18.
//  Copyright © 2018 Goo Jun. All rights reserved.
//

#ifndef bam_cram_h
#define bam_cram_h

#include "gcContent.h"

// A class for a single file BAM I/O
class BamCram
{
public:
    aux_t **data;
    hts_idx_t* idx;
    gcContent& GC;
    
    
    // Get GC corrected depth for chr / pos
    double gcCorrected(double, int, int);
    std::vector< uint16_t * > depth100; // to store depth for every 100bp interval
    std::vector<int> nbin_100;
    
    std::vector<uint64_t> gc_sum;
    std::vector<uint64_t> gc_cnt;
    
    std::vector<readpair> vec_rp;
    std::vector<splitread> vec_sp;
    
    //    void read_depth(std::vector<sv> &, std::vector<string> &);
    void read_depth_sequential(std::vector<breakpoint> &, std::vector<sv> &);
    //    void process_readpair(sv &, std::vector<int> &, string &);
    //    void get_avg_depth();
    void initialize(string &);
    void initialize_sequential(string &);
    void postprocess_depth(std::vector<sv> &);

    BamCram(gcContent &x) : GC(x) {};
    
};

#endif /* bam_cram_h */
