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
#include "gaussian_mixture2.h"
#include "pileup.h"
#include <math.h>

class SvGeno
{
public:
    bool b_biallelic;
    bool b_pass;

    bool cnt_flag;
    bool dp_flag;
    bool dp2_flag;
    bool pd_flag;
    bool dpcnt_flag;
    bool readpair_flag;
    bool split_flag;
    bool clip_flag;
    double rp_sum;
    double split_sum;
    double startclip_sum;
    double endclip_sum;

    double MAX_P_OVERLAP;
    
    std::string info;
    
    int n_sample;
    int n_effect;
    
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

    std::vector<double> start_clips;
    std::vector<double> end_clips;
    
    std::vector<double> split_cnts;
    std::vector<double> rp_cnts;
    std::vector<double> all_cnts;
    
    std::vector<int> start_rps;
    std::vector<int> end_rps;
    
    std::vector<bool> nonref_mask;
    std::vector<bool> sample_mask; // true: include sample in genotyping, false: exclude sample

    SvGeno (int);
    
    void reset();
};

class BreakCluster
{
public:
    BreakCluster();
    double get_distance(std::pair<int, int> &);
    double get_distance(int, int);
    double get_distance(BreakCluster&);
    void merge(BreakCluster&);

    void add_to_cluster(int, int, int);
    void add_to_cluster(std::pair<int, int> &);

    double start_mean;
    double start_var;
    double end_mean;
    double end_var;
    int N;
};

class SvData
{
public:
    int n_sample;
    std::vector<ReadStat> rdstats;
    std::vector<double> raw_dp;
    std::vector< std::vector<double> > dps;
    std::vector<double> prepost_dp;
    
    // Vectors to store total # rp/sp across all samples per every 10 bp in +/- 500bp of the SV
    std::vector< std::vector<int> >all_rps;
    std::vector<int> all_sps;
    // Vector to store total # clips cross all samples per every bp in +/- 100bp of the SV
    std::vector<BreakCluster> vec_break_clusters;
    int clus_idx;
    
    std::vector<int> all_lclips;
    std::vector<int> all_rclips;

    bool multi_dp;
    
    void reset();

    SvData (int);
};


class SvStat
{
public:
    GaussianMixture depth_mix;
    GaussianMixture split_mix;
    GaussianMixture rp_mix;
    GaussianMixture lclip_mix;
    GaussianMixture rclip_mix;

    GaussianMixture2 dpcnt_mix;

    std::string get_summary_stat(SvGeno&, SvData&);
};


class Genotyper
{
public:
    int n_sample;

    SvStat genostat; 

    void check_biallelic(SvGeno &);
    void get_prepost_stat(SvData &, SvGeno &);
    int find_peak(std::vector<int> &, int, int);
    bool find_consensus_rp(sv &, SvData &, SvGeno &, int, int &, int &);
    
    bool find_consensus_split(sv &, SvData &, SvGeno &);
    bool is_pairsplit_oriented(sv &, PairSplit &);
    
    bool find_consensus_clip(sv &, SvData&, SvGeno &);
    bool find_consensus_clip(sv &, SvData &, SvGeno &, int, int &, int &);
    
    bool find_consensus_clip_inv(sv &, SvData &, SvGeno &, int &, int &, int&, int&);
    bool get_del_cnts(sv &, SvData &, SvGeno &);
    bool get_dup_cnts(sv &, SvData &, SvGeno &);
    bool get_inv_cnts(sv &, SvData &, SvGeno &);

    void print_genodata(sv &, SvData &, SvGeno &);

    void call( sv&,  SvData&, SvGeno &, std::vector<SampleStat> &);
    void call_deletion( sv &,  SvData &, SvGeno &, std::vector<SampleStat> &);
    void call_cnv( sv &,  SvData &, SvGeno &, std::vector<SampleStat> &);
    void call_inversion(sv &, SvData &, SvGeno &, std::vector<SampleStat> &);

    bool assign_del_genotypes(sv &, SvData &, SvGeno &, std::vector<SampleStat> &);
    bool assign_dup_genotypes(sv &, SvData &, SvGeno &, std::vector<SampleStat> &);
    bool assign_inv_genotypes(sv &, SvData &, SvGeno &, std::vector<SampleStat> &);

  //  void call_insertion(sv &, SvData &, SvGeno &);
    void select_model_1d(GaussianMixture &, std::vector< std::vector<double> > &, std::vector<double> &, double);
    void select_model_mask_1d(GaussianMixture &, std::vector< std::vector<double> > &, std::vector<double> &, std::vector<bool> &, double);
    void select_model_mask_2d(GaussianMixture2 &, std::vector< std::vector<double> > &, std::vector<double> &, std::vector<double>&, std::vector<bool> &, double);
    void select_model_dpcnt_mask(GaussianMixture2 &, std::vector< std::vector<double> > &, std::vector< std::vector<double> > &, std::vector<double> &, std::vector<double>&, std::vector<bool> &, double);

    void select_model_2d(GaussianMixture2 &, std::vector< std::vector<double> > &, std::vector<double> &, std::vector<double>&, double);
};

#endif /* genotyper_h */
