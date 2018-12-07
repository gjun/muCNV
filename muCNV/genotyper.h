//
//  genotyper.h
//  muCNV
//
//  Created by Goo Jun on 12/3/18.
//  Copyright © 2018 Goo Jun. All rights reserved.
//

#ifndef genotyper_h
#define genotyper_h

#include "sv.h"
#include "gaussian.h"
#include "data_reader.h"
#include "gaussian_mixture.h"
#include "pileup.h"
#include <math.h>

class SvGeno
{
public:
    bool b_biallelic;
    bool b_pass;
    bool dp_flag;
    bool dp2_flag;
    bool read_flag;
    
    std::string info;
    
    int n_sample;
    
    int ac;
    int ns;
    
    std::vector<int> gt; // bi-allelic genotype
    std::vector<int> cn; // copy number

    SvGeno (int);
};

class SvData
{
public:
    int n_sample;
    std::vector<ReadStat> rdstats;
    std::vector<double> var_depth;
    std::vector< std::vector<double> > dp2;
    
    SvData (int);
};

class Genotyper
{
public:
    unsigned n_sample;
    void call(sv&, SvData&, SvGeno &);
    void call_deletion(sv &, SvData &, SvGeno &);
    void call_cnv(sv &, SvData &, SvGeno &);
    void call_inversion(sv &, SvData &, SvGeno &);
    void select_model(GaussianMixture &, std::vector< std::vector<double> > &, std::vector<double> &);
    void select_model(GaussianMixture2 &, std::vector< std::vector<double> > &, std::vector<double> &, std::vector<double>&);

    std::string print(sv &, SvData&, SvGeno&);
};

bool ordered(std::vector<Gaussian> &);


#endif /* genotyper_h */
