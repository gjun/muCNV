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
	bool rp_geno_flag;
    bool clip_flag;
	bool clip_geno_flag;
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
    
    int rp_pos;
    int rp_end;
    int clip_pos;
    int clip_end;
    
	GaussianMixture gmix;
    
	GaussianMixture2 gmix2;
    
    std::vector<int> gt; // bi-allelic genotype
    std::vector<int> rp_gt;
    std::vector<int> clip_gt;

    std::vector<int> cn; // copy number
    std::vector<int> rp_cn;
    std::vector<int> clip_cn;
    
    std::vector<int> start_clips;
    std::vector<int> end_clips;
    
    std::vector<int> start_rps;
    std::vector<int> end_rps;

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
    std::vector< std::vector<int> >all_rps;
    std::vector<int> all_sps;
    // Vector to store total # clips cross all samples per every bp in +/- 100bp of the SV

    std::vector<int> all_lclips;
    std::vector<int> all_rclips;

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
    int find_peak(std::vector<int> &, int, int);
    bool find_consensus_rp(SvData &, int, int &, int &);
    bool find_consensus_clip(SvData &, int, int &, int &);
    bool find_consensus_clip_inv(SvData &, int &, int &, int&, int&);

    void call( sv&,  SvData&, SvGeno &, bool, bool, std::vector<SampleStat> &);
    void call_deletion( sv &,  SvData &, SvGeno &);
    void call_cnv( sv &,  SvData &, SvGeno &);
    void call_inversion(sv &, SvData &, SvGeno &, std::vector<SampleStat> &);
  //  void call_insertion(sv &, SvData &, SvGeno &);
    void select_model(GaussianMixture &, std::vector< std::vector<double> > &, std::vector<double> &, double);
    void select_model(GaussianMixture2 &, std::vector< std::vector<double> > &, std::vector<double> &, std::vector<double>&, double);
    void select_model_mask(GaussianMixture &, std::vector< std::vector<double> > &, std::vector<double> &x, std::vector<bool> &mask, double MAX_P_OVERLAP);

};

#endif /* genotyper_h */
