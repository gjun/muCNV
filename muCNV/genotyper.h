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
	double MAX_P_OVERLAP;
    
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
    
    void reset();
};

class SvData
{
public:
    int n_sample;
    std::vector<ReadStat> rdstats;
    std::vector< std::vector<double> > dps;
    std::vector<double> prepost_dp;
    
    // Vectors to store total # rp/sp across all samples per every 10 bp in +/- 500bp of the SV
    std::vector<int32_t> all_rps;
    std::vector<int32_t> all_sps;
    // Vector to store total # clips cross all samples per every bp in +/- 100bp of the SV
    // + : clip in forward direction (rclip)
    // - : clip in reverse direction (lclip)
    std::vector<int32_t> all_clips;
    bool multi_dp;
    
    void reset();

    SvData (int);
};

class Genotyper
{
public:
    int n_sample;
	bool b_kmeans;
	bool b_mahalanobis;

    void get_prepost_stat(SvData &, SvGeno &);

    void call( sv&,  SvData&, SvGeno &, bool, bool, std::vector<SampleStat> &);
    void call_deletion( sv &,  SvData &, SvGeno &);
    void call_cnv( sv &,  SvData &, SvGeno &);
    void call_inversion(sv &, SvData &, SvGeno &, std::vector<SampleStat> &);
  //  void call_insertion(sv &, SvData &, SvGeno &);
    void select_model(GaussianMixture &, std::vector< std::vector<double> > &, std::vector<double> &, double);
    void select_model(GaussianMixture2 &, std::vector< std::vector<double> > &, std::vector<double> &, std::vector<double>&, double);
};

#endif /* genotyper_h */
