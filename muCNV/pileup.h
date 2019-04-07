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
#include "base_file.h"

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

class sclip
{
public:
    int8_t chrnum;
    int32_t pos;
    bool b_end; // to indiciate whether the clip is at the cycle end or not
    bool b_drop;
    
    bool operator < (const sclip&) const;
    bool operator == (const sclip&) const;
    bool operator <= (const sclip&) const;
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

class Pileup : public BaseFile
{
public:

    int write_sample_stat(SampleStat &);
    int write_gc_factor(std::vector<double>&, int);
    int write_depth(uint16_t*, int);
    int write_readpair(readpair&);
    int write_splitread(splitread&);
    int write_softclip(sclip&);

    int read_sample_stat(SampleStat &);
    int read_gc_factor(std::vector<double>&, int);
    int read_depth(uint16_t*, int);
    int read_readpair(readpair&);
    int read_splitread(splitread&);
    int read_softclip(sclip&);
};


#endif /* pileup_h */
