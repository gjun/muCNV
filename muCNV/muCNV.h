#ifndef __COMMON_H__
#define __COMMON_H__

#include <stdint.h>
#include <stdio.h>
#include <iostream>
#include <string>
#include <fstream>
#include <float.h>
#include <math.h>
#include <vector>
#include <algorithm>
#include <numeric>
// #include "pFile.h"

#include "hts.h"
#include "vcf.h"
#include "bgzf.h"
#include "regidx.h"
#include "sam.h"
#include "tbx.h"
#include "kseq.h"
#include "sv.h"

using std::string;

const double PI=3.14159265358979;
const double sqPI=1.77245385090552;
const double invsqrt2pi= 0.398942280401433;
const double log2pi = 0.7981798684;


typedef std::pair<uint64_t, uint64_t> interval_t;

void split(const char*, const char*, std::vector<string>&);

template <typename T>
std::vector<size_t> sort_indexes(const std::vector<T> &v) {

	// initialize original index locations
	std::vector<size_t> idx(v.size());
	iota(idx.begin(), idx.end(), 0);

	// sort indexes based on comparing values in v
	std::sort(idx.begin(), idx.end(),
			[&v](size_t i1, size_t i2) {return v[i1] > v[i2];});

	return idx;
}

bool in_centrome(sv &);
bool in_centrome(int, int);

int find_overlap_sv(sv& , std::vector<sv>&);


class svdata
{
public:
	int n;
	std::vector<double> dp;
	std::vector<double> isz;
	
	std::vector<double> cnv_pos;
	std::vector<double> cnv_neg;
	std::vector<double> inv_pos;
	std::vector<double> inv_neg;
	
	std::vector<int> n_isz;
	std::vector<int> n_cnv_pos;
	std::vector<int> n_cnv_neg;
	std::vector<int> n_inv_pos;
	std::vector<int> n_inv_neg;
	
	std::vector<double> norm_dp;
	std::vector<double> norm_readcount;
	std::vector<double> norm_cnv_pos;
	std::vector<double> norm_cnv_neg;
	std::vector<double> norm_inv_pos;
	std::vector<double> norm_inv_neg;
	
	void set_size(int num)
	{
		n = num;
		dp.resize(n);
		isz.resize(n);
		
		cnv_pos.resize(n);
		cnv_neg.resize(n);
		inv_pos.resize(n);
		inv_neg.resize(n);
		
		n_isz.resize(n);
		n_cnv_pos.resize(n);
		n_cnv_neg.resize(n);
		n_inv_pos.resize(n);
		n_inv_neg.resize(n);
		
		norm_dp.resize(n);
		norm_cnv_pos.resize(n);
		norm_cnv_neg.resize(n);
		norm_inv_pos.resize(n);
		norm_inv_neg.resize(n);
		norm_readcount.resize(n);
	};
	
	void normalize(sv &interval, std::vector<double> &avg_depth, std::vector<double> &avg_isize)
	{
		for(int i=0;i<n;++i)
		{
			norm_dp[i] = dp[i] / avg_depth[i];
			if (interval.svtype == DEL)
			{
				norm_cnv_pos[i] = (cnv_pos[i] - avg_isize[i] ) / interval.len;
				norm_cnv_neg[i] = (cnv_neg[i] - avg_isize[i] ) / interval.len;
				norm_inv_pos[i] = (inv_pos[i] - avg_isize[i] ) / interval.len;
				norm_inv_neg[i] = (inv_neg[i] - avg_isize[i] ) / interval.len;
			}
			else if (interval.svtype == DUP || interval.svtype == CNV)
			{
				norm_cnv_pos[i] = (cnv_pos[i] + avg_isize[i]) / interval.len;
				norm_cnv_neg[i] = (cnv_neg[i] + avg_isize[i]) / interval.len;
				norm_inv_pos[i] = (inv_pos[i] + avg_isize[i]) / interval.len;
				norm_inv_neg[i] = (inv_neg[i] + avg_isize[i]) / interval.len;
			}
			else if (interval.svtype == INV )
			{
				norm_cnv_pos[i] = (cnv_pos[i] ) / avg_isize[i];
				norm_cnv_neg[i] = (cnv_neg[i] ) / avg_isize[i];
				norm_inv_pos[i] = (inv_pos[i] ) / interval.len;
				norm_inv_neg[i] = (inv_neg[i] ) / interval.len;
			}
			// READLEN fixed to 150 : later!!
			norm_readcount[i] = (double)n_isz[i] * 150.0 / (interval.len + 2.0*(avg_isize[i]-75)) / avg_depth[i];
		}
	};
	
	void print(sv &interval)
	{
	}
	
};


void pick_sv_from_merged(sv &, std::vector<sv> &);

class SampleList
{
public:
    std::vector<string> vID;
    std::vector<string> vDepthFile;
    std::vector<double> AvgDepth;
    unsigned n_sample;
    
    bool readIndex(string);
    bool readAvgDepth();
};


class Gaussian
{
public:
	double Mean;
	double Stdev;
	double Alpha;
	double pdf(const double &);
	void set(const double &, const double &);
	void estimate(std::vector<double> &);
	Gaussian();
};

class Gaussian2
{
public:
	// Parameters for two-dimensional Gaussian
	double Mean[2]; // mean
	double Cov[4]; // covariance
	double Prc[4]; // precision (inverse of covariance)
	double Det; // determinant
	double Alpha;
	double pdf(const double&, const double&);
	void set(const double &, const double &, const double &, const double &, const double &, const double &);
	void estimate(std::vector<double> &, std::vector<double> &);
	double logpdf(const double&, const double&);
	void update(); // update precision matrix
	
	Gaussian2();
};

class invcfs
{
public:
	// std::vector<std::ifstream *> vfs;
	std::vector<htsFile*> vfs;
	std::vector<tbx_t*> tbs;
	std::vector<hts_itr_t *> m_itr;
	std::vector<int> num_id;
	std::vector<int> start_num;

	void parse_sv(std::vector<string> &, sv &);
	double get_second_value(string &);
	void get_value_pair(string &, int &, double &);
	int initialize(std::vector<string> &, std::vector<string> &, std::vector<double> &, std::vector<double> &, std::vector<double> &, string &);
	int read_interval_multi(sv& , svdata &, string &);	
};


class outvcf
{
public:
	FILE *fp;
	int varcnt;
	void open(string&);
	void close();
	void print(string &);
	void write_header(std::vector<string>&);
	void write_del(sv&, std::vector<int>&, std::vector<int>&, int, int, std::vector<double>&, std::vector<double>&, std::vector<Gaussian>&, double, bool);
	void write_cnv(sv&, std::vector<int>&, std::vector<int>&, int, int, std::vector<double>&, std::vector<double>&, std::vector<Gaussian>&, double, bool);
};

class svgeno 
{
public:
	bool b_biallelic;
	bool b_pass;
	bool b_dump;
	
	bool dp_flag;
	bool pos_flag;
	bool neg_flag;

	bool cnv_pos_flag;
	bool cnv_neg_flag;
	bool inv_pos_flag;
	bool inv_neg_flag;

	double p_overlap;
	string info;
	int n_sample;

	int ac;
	int ns;
	std::vector<double> bic;

	std::vector<Gaussian> Comps;

	std::vector<int> gt; // bi-allelic genotype
	std::vector<int> cn; // copy number

	void initialize(int N)
	{
		n_sample = N;
		Comps.clear();
		
		ns = 0;
		ac = 0;
		p_overlap = 1;

		gt.resize(n_sample, -1);
		cn.resize(n_sample, -1);
		bic.resize(20, 0);
		
		for(int i=0;i<n_sample; ++i)
		{
			gt[i] = -1;
			cn[i] = -1;
		}
		for(int i=0;i<20;++i)
		{
			bic[i] = 0;
		}

		b_biallelic = false;
		b_pass = false;
		b_dump = true;

		dp_flag = false;
		pos_flag = false;
		neg_flag = false;
		// For inversions
		cnv_pos_flag = false;
		cnv_neg_flag = false;
		inv_pos_flag = false;
		inv_neg_flag = false;
		info = "";

	};
	void print (sv &, svdata &, string &, std::vector<double> &);
};


class gtype
{
public:
	void call_tmp(sv &, svdata &, svgeno &, std::vector<double> &, std::vector<double> &, std::vector<double> &);
	void call_del_tmp(sv &, svdata &, svgeno &, std::vector<double> &, std::vector<double> &, std::vector<double> &);
	void call_dup_tmp(sv &, svdata &, svgeno &, std::vector<double> &, std::vector<double> &, std::vector<double> &);
	void call_del(sv &, svdata &, svgeno &, std::vector<double> &, std::vector<double> &, std::vector<double> &);
	void call_cnv(sv &, svdata &, svgeno &, std::vector<double> &, std::vector<double> &, std::vector<double> &);
	void call_inv(sv &, svdata &, svgeno &, std::vector<double> &, std::vector<double> &);

	void EM(std::vector<double>&, std::vector<Gaussian>&);
	void EM(std::vector<double>&, std::vector<double>&, std::vector<Gaussian>&);
	void fit(std::vector<double>&, std::vector<Gaussian>&);
//	void EM2(std::vector<double>&, std::vector<double> &, std::vector<Gaussian2>&);

	int assign(double, std::vector<Gaussian> &);
	void copyComps(std::vector<Gaussian> &, std::vector<Gaussian> &);
//	void call_del(sv&, std::vector<double>&, std::vector<double>&, std::vector<int>&, outvcf&, std::vector<double>&);
//	void call_cnv(sv&, std::vector<double>&, std::vector<double>&, std::vector<int>&, outvcf&, std::vector<double>&);
};


template <class T> void vprint(std::vector<T> &);

double MAX(std::vector<double>&);
double MIN(std::vector<double>&);
double normpdf(double, Gaussian&);
double mean(std::vector<double>&);
double stdev(std::vector<double>&, double);

int find_start(std::vector<sv> &, int );

void write_interval(string &, std::vector<sv> &);

void read_svs_from_vcf(string &, std::vector<breakpoint> &, std::vector<sv> &);
void read_svs_from_intfile(string &, std::vector<breakpoint> &, std::vector<sv> &);
void read_intervals_from_vcf(std::vector<string> &, std::vector<string> &, std::vector<sv> &);
int read_candidate_vcf(std::ifstream &, sv&, string& );

double RO(sv &, sv &);
void merge_svs(std::vector<sv> &, std::vector<int> &);
void cluster_svs(std::vector<sv>&, std::vector< std::vector<sv> > &, double);

double BayesError(std::vector<Gaussian>&);
double BayesError(std::vector<Gaussian2>&);

double BIC(std::vector<double>&, std::vector<Gaussian>&);
double BIC(std::vector<double>&, std::vector<Gaussian>&, std::vector<double> &);
double BIC(std::vector<double>&, std::vector<double>&, std::vector<Gaussian2>&);
double det(double*);

int median(std::vector<int> &);

bool ordered(std::vector<Gaussian> &);

void read_list(string &, std::vector<string>&);
void read_index(string, std::vector<string>&, std::vector<string>&, std::vector<string>&, std::vector<double>&);

void readIndex(string, std::vector<string>&, std::vector<string>&, std::vector<string>&, std::vector<string>&);
void readDepthOrig(std::vector<string>&, std::vector<string>&, std::vector<interval_t>&, std::vector< std::vector<double> >&, std::vector<double>&);
void readDepth(std::vector<string>&, std::vector<sv>&, std::vector< std::vector<double> >&, std::vector<double>&);

double getAvgDepth(string, string);
string getCNVsegmentFileName(string, string);
void readsv(string, int, std::vector<sv>&, std::vector<sv>&);


void printCluster(std::vector<Gaussian> &C);



#endif // __COMMON_H__
