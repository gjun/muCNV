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
	int n_pre_FR;  // DEL
	int n_pre_RF;  // DUP
	int n_post_FR; // DEL
	int n_post_RF; // DUP

	int n_pre_rp_missing; 
	int n_post_rp_missing;

	int n_pre_sp_missing; 
	int n_post_sp_missing;

    //  ---------(sclip)                  (sclip)-------- split_in, DEL
    //  (sclip)---------                  ----------(sclip) split_out, DUP
	int n_pre_split_out; // right clipped
	int n_pre_split_in; // left clipped
	int n_post_split_out; // right clipped
	int n_post_split_in; // left clipped

	ReadStat() 
	{
		n_pre_FR = 0;
		n_pre_RF = 0;
		n_post_FR = 0;
		n_post_RF = 0;
		n_pre_rp_missing = 0;
		n_post_rp_missing =0;

		n_pre_sp_missing = 0;
		n_post_sp_missing =0;
		n_pre_split_out = 0;
		n_pre_split_in = 0;
		n_post_split_out = 0;
		n_post_split_in = 0;
	}
};

class DataReader
{
public:
    std::vector<std::string> sample_ids;

    // Initialize all multi-pileups, load number of samples & sapmle ids & gc factors
    int load(std::vector<string> &, std::vector<SampleStat>&, GcContent &);
    int read_depth100(sv&, std::vector< std::vector<double> > &, std::vector< std::vector<double> >&, GcContent& gc);
    void read_var_depth(int, std::vector<double>&);
    void read_pair_split(sv&, std::vector<ReadStat> &, GcContent &);
    double correct_gc(GcContent &, int, double, int, int);

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
