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
#include <iostream>
#include <fstream>

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

class SampleStat
{
public:
    double avg_dp;
    double std_dp;
    double avg_isize;
    double std_isize;
    double med_isize;
    double avg_rlen;
};

class Pileup
{
public:
    SampleStat stat;
    
    GcContent& gc;
    // Get GC corrected depth for chr / pos
    double gcCorrected(double, int, int);

    std::vector<double> gc_factor;
    
    std::vector< uint16_t * > depth100; // to store depth for every 100bp interval
    std::vector<int> nbin_100;

    
    std::vector<readpair> vec_rp;
    std::vector<splitread> vec_sp;
    
    Pileup(GcContent &x) : gc(x) {};
    void initialize();
    
    int write_number(std::ofstream&, int);
    int write_sample_id(std::ofstream&, std::string &);
    int write_sample_stat(std::ofstream&, SampleStat &);
    int write_gc_factor(std::ofstream&, std::vector<double>&);
    int write_depth(std::ofstream&, uint16_t*, int);
    int write_readpair(std::ofstream&, readpair&);
    int write_splitread(std::ofstream&, splitread&);
    
    int read_number(std::ifstream&, int&);
    int read_sample_id(std::ifstream&, char *);
    int read_sample_stat(std::ifstream&, SampleStat &);
    int read_gc_factor(std::ifstream&, std::vector<double>&);
    int read_depth(std::ifstream&, uint16_t*, int);
    int read_readpair(std::ifstream&, readpair&);
    int read_splitread(std::ifstream&, splitread&);
};


#endif /* pileup_h */
