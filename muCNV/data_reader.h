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

class DataReader
{
public:
    std::vector<std::string> sample_ids;

    // Initialize all multi-pileups, load number of samples & sapmle ids & gc factors
    int load(std::vector<string> &, std::vector<SampleStat>&, GcContent &);
    int read_depth100(sv&, std::vector< std::vector<double> > &, GcContent& gc);
    void read_var_depth(int, std::vector<double>&);
    void read_pair_split(sv&, std::vector< std::vector<readpair> > &, std::vector< std::vector<splitread> > &);
    double correct_gc(GcContent &, int, double, int, int);

private:
    std::vector<int32_t> n_samples; // number of samples in each multi pileup files

    int n_sample_total;
    int n_var;
    int n_pileup;
    std::vector<Pileup> pileups;
    std::vector<Pileup> var_files;
    std::vector<BaseFile> idx_files;
    
    // Number of pileups * number of indices per sapmle (~ 300,000)
    std::vector< std::vector<uint64_t> > multi_idx;
    std::vector< std::vector<double>> gc_factors;

    size_t get_dp100_offset(uint64_t, int, int);
};

#endif /* data_reader_h */
