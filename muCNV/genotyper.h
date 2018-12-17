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
    bool pd_flag;
    bool read_flag;
    
    std::string info;
    
    int n_sample;
    
    int ac;
    int ns;
	// Average depth around SVs in samples with possible het/hom genotypes
	double dp_pre_mean;
    double dp_pre_std;
    double dp_post_mean;
    double dp_post_std;
    bool b_pre;
    bool b_post;
    
	GaussianMixture gmix;
	GaussianMixture2 gmix2;
    
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
    std::vector<double> prepost_dp;
    
    SvData (int);
};

class Genotyper
{
public:
    int n_sample;
	double MAX_P_OVERLAP;
	bool b_kmeans;
	bool b_mahalanobis;

    void get_prepost_stat(SvData &, SvGeno &);

    void call( sv&,  SvData&, SvGeno &, double, bool, bool);
    void call_deletion( sv &,  SvData &, SvGeno &);
    void call_cnv( sv &,  SvData &, SvGeno &);
    void call_inversion(sv &, SvData &, SvGeno &);
    void select_model(GaussianMixture &, std::vector< std::vector<double> > &, std::vector<double> &);
    void select_model(GaussianMixture2 &, std::vector< std::vector<double> > &, std::vector<double> &, std::vector<double>&);

    std::string print(sv &, SvData&, SvGeno&);
};

bool ordered(std::vector<Gaussian> &);


#endif /* genotyper_h */
