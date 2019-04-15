//
//  multi_pileup.hpp
//  muCNV
//
//  Created by Goo Jun on 11/25/18.
//  Copyright © 2018 Goo Jun. All rights reserved.
//

#ifndef data_reader_h
#define data_reader_h

#include <stdio.h>
#include "pileup.h"
#include "gc_content.h"

using std::string;

class ReadStat
{
public:
    std::vector<int> n_rp; // indexed by pairstr

    int n_split_inward;
    int n_split_outward;
 
    bool del_support();
    bool dup_support();
    bool inv_support();
    bool ins_support();
    
    // Vectors to store total # rp/sp in each sample per every 10 bp in +/- 500bp of the SV
    std::vector< std::vector<int> > rp_seq;
    std::vector<int> sp_seq_in;
    std::vector<int> sp_seq_out;

    
    // Vector to store total # clips in each sample per every bp in +/- 100bp of the SV
    // + : clip in forward direction (rclip)
    // - : clip in reverse direction (lclip)
    std::vector<int> lclips;
    std::vector<int> rclips;
    
    int n_lclip_start;
    int n_lclip_end;
    
    int n_rclip_start;
    int n_rclip_end;
    
    void reset()
    {
        std::fill(n_rp.begin(), n_rp.end(), 0);
        
        n_split_inward = 0;
        n_split_outward = 0;
        for(int i=0; i<4; ++i)
        {
            std::fill(rp_seq[i].begin(), rp_seq[i].end(), 0);
        }

        std::fill(sp_seq_out.begin(), sp_seq_out.end(), 0);
        std::fill(sp_seq_in.begin(), sp_seq_in.end(), 0);
        
        std::fill(lclips.begin(), lclips.end(), 0);
        std::fill(rclips.begin(), rclips.end(), 0);
        n_lclip_start = 0;
        n_rclip_start = 0;
        n_lclip_end = 0;
        n_rclip_end = 0;
    };
    
	ReadStat() 
	{
        n_rp.resize(4);
        
        rp_seq.resize(4);

        for(int i=0; i<4; ++i)
        {
            rp_seq[i].resize(200);
        }
        sp_seq_out.resize(400);
        sp_seq_in.resize(400);
        
        lclips.resize(400);
        rclips.resize(400);
    };
    
};

class DataReader
{
public:
    std::vector<std::string> sample_ids;

    // Initialize all multi-pileups, load number of samples & sapmle ids & gc factors
    int load(std::vector<string> &, std::vector<SampleStat>&, GcContent &, int );
    bool read_depth100(sv&, std::vector< std::vector<double> > &, GcContent& gc, bool);
    void read_var_depth(int, std::vector<double>&);
    void read_pair_split(sv&, std::vector<ReadStat> &, GcContent &, std::vector< std::vector<int> >&, std::vector<int> &, std::vector<int> &);
    double correct_gc(GcContent &, int, double, int, int);
    bool around_breakpoint(readpair &, sv &);
    void close();

private:
    std::vector<int32_t> n_samples; // number of samples in each multi pileup files

    int n_sample_total;
    int n_var;
    int n_pileup;
    std::vector<uint64_t> chr_idx_rp;
    std::vector< std::vector<uint64_t> > chr_bytepos_dp100;
    std::vector<Pileup> pileups;
    std::vector<Pileup> var_files;
    std::vector<BaseFile> idx_files;
    
    // Number of pileups * number of indices per sapmle (~ 300,000)
    std::vector< uint64_t* > multi_idx;
    std::vector< std::vector<double>> gc_factors;

    size_t get_dp100_offset(uint64_t, int, int);
};

#endif /* data_reader_h */
