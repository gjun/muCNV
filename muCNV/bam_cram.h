//
//  bam_cram.h
//  muCNV
//
//  Created by Goo Jun on 11/13/18.
//  Copyright © 2018 Goo Jun. All rights reserved.
//

#ifndef bam_cram_h
#define bam_cram_h

// #include "gc_content.h"
#include "pileup.h"
#include <htslib/hts.h>
#include <htslib/sam.h>

typedef struct {     // auxiliary data structure
    samFile *fp;     // the file handle
    bam_hdr_t *hdr;  // the file header
    hts_itr_t *iter; // NULL if a region not specified
    
    uint64_t sum_isz;
    uint64_t sumsq_isz;
    uint64_t n_isz;
    
    uint32_t n_rp;
    uint32_t n_sp;
    
    std::vector<readpair> *p_vec_rp;
    std::vector<splitread> *p_vec_sp;
    std::vector<sclip> *p_vec_rclip;
    std::vector<sclip> *p_vec_lclip;
    
    int min_mapQ, min_len; // mapQ filter; length filter
} aux_t;

// A class for a single file BAM I/O
class BamCram
{
public:
    aux_t **data;
    hts_idx_t* idx;
   // GcContent& GC;
    
    std::vector<double> gc_factor;
    std::vector< std::vector<int>> gc_count;
    std::vector< uint16_t * > depth_interval; // to store depth for every 100bp interval
    
    SampleStat stat;
    std::vector<readpair> vec_rp;
    std::vector<splitread> vec_sp;
    std::vector<sclip> vec_rclip;
    std::vector<sclip> vec_lclip;
    
    //    void read_depth(std::vector<sv> &, std::vector<string> &);
    
    void read_depth_sequential(Pileup &, GcContent&, std::vector<breakpoint> &, std::vector<sv> &);
    
    int flag_softclips(std::vector<sclip> &);
    //    void process_readpair(sv &, std::vector<int> &, string &);
    //    void get_avg_depth();
    
   // void initialize(string &);
    void initialize_sequential(std::string &, GcContent&);
    void postprocess_depth(std::vector<sv> &);
    
};

#endif /* bam_cram_h */
