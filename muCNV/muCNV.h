#ifndef __COMMON_H__
#define __COMMON_H__

#include <stdint.h>
#include <stdio.h>
#include <map>
#include <set>
#include <iostream>
#include <string.h>
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

const double PI=3.14159265358979;
const double sqPI=1.77245385090552;
const double invsqrt2pi= 0.398942280401433;
const double log2pi = 0.7981798684;
using namespace std;


enum svType {DEL=0, DUP=1, INV=2, CNV=3, INS=4, BND=5};
string svTypeName(svType t);

typedef pair<uint64_t, uint64_t> interval_t;

void split(const char*, const char*, std::vector<std::string>&);

template <typename T>
vector<size_t> sort_indexes(const vector<T> &v) {

	// initialize original index locations
	vector<size_t> idx(v.size());
	iota(idx.begin(), idx.end(), 0);

	// sort indexes based on comparing values in v
	sort(idx.begin(), idx.end(),
			[&v](size_t i1, size_t i2) {return v[i1] > v[i2];});

	return idx;
}

class readpair
{
public:
    int32_t selfpos;
    int32_t matepos;
    bool selfstr;
    bool matestr;
};

class splitread
{
public:
    int32_t pos;
    int32_t sapos;
    int16_t firstclip; // soft clip position (+: left-side, -: right-side) in primary alignment
    int16_t secondclip; // soft clip position (+: left-side, -: right-side) in secondary alignment
};

class sv
{
	public:
	svType svtype;
//	string source;
//	string chr;
	int chrnum;
	int pos;
	int end;
	int len;
	int supp;
//	pair<int,int> ci_pos;
//	pair<int,int> ci_end;
	void get_len()
	{
		len = end - pos + 1;
	};
	bool operator < (const sv&) const;
	bool operator == (const sv&) const;

    uint64_t dp_sum;
    int n_dp;
    uint8_t dp;
    
    vector<readpair> vec_pair;
    vector<splitread> vec_split;
	sv();
};

bool in_centrome(sv &);

class svdata
{
public:
	int n;
	vector<double> dp;
	vector<double> isz;
	
	vector<double> cnv_pos;
	vector<double> cnv_neg;
	vector<double> inv_pos;
	vector<double> inv_neg;
	
	vector<int> n_isz;
	vector<int> n_cnv_pos;
	vector<int> n_cnv_neg;
	vector<int> n_inv_pos;
	vector<int> n_inv_neg;
	
	vector<double> norm_dp;
	vector<double> norm_readcount;
	vector<double> norm_cnv_pos;
	vector<double> norm_cnv_neg;
	vector<double> norm_inv_pos;
	vector<double> norm_inv_neg;
	
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
	
	void normalize(sv &interval, vector<double> &avg_depth, vector<double> &avg_isize)
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

class gcint : public sv
{
	public:
	uint8_t gcbin;
};

class breakpoint
{
public:
    uint8_t chrnum;
    int bptype; // 0 : pos-gap, 1: pos, 2: pos+gap, 3:end-gap, 4: end, 5: end+gap
    
	int pos;
	int idx; // SV_id
	bool operator < (const breakpoint&) const;
	bool operator == (const breakpoint&) const;
    bool operator <= (const breakpoint&) const;

	breakpoint();
};

void pick_sv_from_merged(sv &, vector<sv> &);

class SampleList
{
public:
    vector<string> vID;
    vector<string> vDepthFile;
    vector<double> AvgDepth;
    unsigned n_sample;
    
    bool readIndex(string);
    bool readAvgDepth();
};

typedef struct {     // auxiliary data structure
	samFile *fp;     // the file handle
	bam_hdr_t *hdr;  // the file header
	hts_itr_t *iter; // NULL if a region not specified
    set<int> *rp_set;
    set<int> *sp_set;
    vector<sv> *vec_sv;
    uint64_t sum_isz;
    uint64_t sumsq_isz;
    uint64_t n_isz;

    //    set<int> *isz_set;
//	vector<double> *isz_sum;
//	vector<int> *isz_cnt;
//	vector<double> *isz_sum;
//	vector<int> *isz_cnt;
//	vector< vector <int> > *isz_list;
//	vector< vector <int> > *pos_list;

//	vector< vector <int> > *rev_isz_list;
//	vector< vector <int> > *rev_pos_list;
	int min_mapQ, min_len; // mapQ filter; length filter
} aux_t;

class Gaussian
{
public:
	double Mean;
	double Stdev;
	double Alpha;
	double pdf(const double &);
	void set(const double &, const double &);
	void estimate(vector<double> &);
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
	void estimate(vector<double> &, vector<double> &);
	double logpdf(const double&, const double&);
	void update(); // update precision matrix
	
	Gaussian2();
};

class invcfs
{
public:
	// vector<ifstream *> vfs;
	vector<htsFile*> vfs;
	vector<tbx_t*> tbs;
	vector<hts_itr_t *> m_itr;
	vector<int> num_id;
	vector<int> start_num;

	void parse_sv(vector<string> &, sv &);
	double get_second_value(string &);
	void get_value_pair(string &, int &, double &);
	int initialize(vector<string> &, vector<string> &, vector<double> &, vector<double> &, vector<double> &, string &);
	int read_interval_multi(sv& , svdata &, string &);	
};

class gcContent
{
public:
	void initialize(string &); // filename for GC content, populate all vectors
	uint16_t num_bin; // Number of GC bin
	uint8_t num_chr; //number of chrs
	uint16_t binsize; // Size of GC bin (bp)
	uint16_t total_bin; // Size of intervals per GC bin
    vector<size_t> chrOffset;
	vector<gcint> regions; // Double array to store list of regions for each GC bin -- non-overlapping, so let's just be it out-of-order
	vector<double> gc_dist; // Array to store proportion of Ref genome for each GC content bin
	vector<uint32_t> chrSize;
	vector<uint8_t *> gc_array; // Array to store "GC bin number" for every 400-bp (?) interval of reference genome
};

// A class for a single file BAM I/O
class bFile
{
public:
	aux_t **data;
	hts_idx_t* idx;
	gcContent& GC;
	double avg_dp;
    double std_dp;
    
	double avg_isize;
	double std_isize;
	double med_isize;
	double avg_rlen;
	
	vector<double> gc_factor;
    vector<uint64_t> gc_sum;
    vector<uint64_t> gc_cnt;
    
	// Get GC corrected depth for chr / pos
	double gcCorrected(double, int, int);
    vector< uint8_t * > depth100; // to store depth for every 100bp interval

    
//	void read_depth(vector<sv> &, vector<string> &);
    void read_depth_sequential(vector<breakpoint> &, vector<sv> &);
//	void process_readpair(sv &, vector<int> &, string &);
//	void get_avg_depth();
	void initialize(string &);
    void initialize_sequential(string &);
    void postprocess_depth(vector<sv> &);
    void write_pileup(string &, vector<sv> &);
    void write_pileup_text(string &, vector<sv> &);
    void write_interval(string &, vector<sv> &);
    void write_depth100(string &);
	bFile (gcContent &x) : GC(x) {};

};

class outvcf
{
public:
	FILE *fp;
	int varcnt;
	void open(string&);
	void close();
	void print(string &);
	void write_header(vector<string>&);
	void write_del(sv&, vector<int>&, vector<int>&, int, int, vector<double>&, vector<double>&, vector<Gaussian>&, double, bool);
	void write_cnv(sv&, vector<int>&, vector<int>&, int, int, vector<double>&, vector<double>&, vector<Gaussian>&, double, bool);
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
	vector<double> bic;

	vector<Gaussian> Comps;

	vector<int> gt; // bi-allelic genotype
	vector<int> cn; // copy number

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
	void print (sv &, svdata &, string &, vector<double> &);
};


class gtype
{
public:
	void call_tmp(sv &, svdata &, svgeno &, vector<double> &, vector<double> &, vector<double> &);
	void call_del_tmp(sv &, svdata &, svgeno &, vector<double> &, vector<double> &, vector<double> &);
	void call_dup_tmp(sv &, svdata &, svgeno &, vector<double> &, vector<double> &, vector<double> &);
	void call_del(sv &, svdata &, svgeno &, vector<double> &, vector<double> &, vector<double> &);
	void call_cnv(sv &, svdata &, svgeno &, vector<double> &, vector<double> &, vector<double> &);
	void call_inv(sv &, svdata &, svgeno &, vector<double> &, vector<double> &);

	void EM(vector<double>&, vector<Gaussian>&);
	void EM(vector<double>&, vector<double>&, vector<Gaussian>&);
	void fit(vector<double>&, vector<Gaussian>&);
//	void EM2(vector<double>&, vector<double> &, vector<Gaussian2>&);

	int assign(double, vector<Gaussian> &);
	void copyComps(vector<Gaussian> &, vector<Gaussian> &);
//	void call_del(sv&, vector<double>&, vector<double>&, vector<int>&, outvcf&, vector<double>&);
//	void call_cnv(sv&, vector<double>&, vector<double>&, vector<int>&, outvcf&, vector<double>&);
};


template <class T> void vprint(vector<T>);

double MAX(vector<double>&);
double MIN(vector<double>&);
double normpdf(double, Gaussian&);
double mean(vector<double>&);
double stdev(vector<double>&, double);

int find_start(vector<sv> &, int );

void read_svs_from_vcf(string &, vector<breakpoint> &, vector<sv> &);
void read_intervals_from_vcf(vector<string> &, vector<string> &, vector<sv> &);
int read_candidate_vcf(ifstream &, sv&, string& );

double RO(sv &, sv &);
void merge_svs(vector<sv> &, vector<int> &);
void cluster_svs(vector<sv>&, vector< vector<sv> > &);

double BayesError(vector<Gaussian>&);
double BayesError(vector<Gaussian2>&);

double BIC(vector<double>&, vector<Gaussian>&);
double BIC(vector<double>&, vector<Gaussian>&, vector<double> &);
double BIC(vector<double>&, vector<double>&, vector<Gaussian2>&);
double det(double*);

int median(vector<int> &);

bool ordered(vector<Gaussian> &);

void read_vcf_list(string &, vector<string>&);
void read_index(string, vector<string>&, vector<string>&, vector<string>&, vector<double>&);

void readIndex(string, vector<string>&, vector<string>&, vector<string>&, vector<string>&);
void readDepthOrig(vector<string>&, vector<string>&, vector<interval_t>&, vector< vector<double> >&, vector<double>&);
void readDepth(vector<string>&, vector<sv>&, vector< vector<double> >&, vector<double>&);

double getAvgDepth(string, string);
string getCNVsegmentFileName(string, string);
void readsv(string, int, vector<sv>&, vector<sv>&);


void printCluster(vector<Gaussian> &C);



#endif // __COMMON_H__
