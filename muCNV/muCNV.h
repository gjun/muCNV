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
// #include "pFile.h"

#include "hts.h"
#include "vcf.h"
#include "bgzf.h"
#include "regidx.h"
#include "sam.h"
#include "tbx.h"
#include "kseq.h"

const double PI=3.1415926535897932384626433832795028841968;
const double sqPI=1.7724538509055160272981674833411451827974;
const double log2pi = 0.7981798684;

using namespace std;

typedef pair<uint64_t, uint64_t> interval_t;

void split(const char*, const char*, std::vector<std::string>&);

class sv
{
	public:
	string svtype;
//	string source;
	string chr;
	int chrnum;
	int pos;
	int end;
	int len;
//	pair<int,int> ci_pos;
//	pair<int,int> ci_end;
	void get_len()
	{
		len = end - pos + 1;
	};
	bool operator < (const sv&) const;
	bool operator == (const sv&) const;
	
	sv();
};

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
	vector<double> norm_cnv_pos;
	vector<double> norm_cnv_neg;
	vector<double> norm_inv_pos;
	vector<double> norm_inv_neg;
	vector<double> norm_readcount;
	
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
			/*
			norm_cnv_pos[i] = (cnv_pos[i] - isz[i]) / interval.len;
			norm_cnv_neg[i] = (cnv_neg[i] - isz[i]) / interval.len;
			norm_inv_pos[i] = (inv_pos[i] - isz[i]) / interval.len;
			norm_inv_neg[i] = (inv_neg[i] - isz[i]) / interval.len;
			 */
			if (interval.svtype == "DEL")
			{
				norm_cnv_pos[i] = (cnv_pos[i] ) / interval.len;
				norm_cnv_neg[i] = (cnv_neg[i] ) / interval.len;
				norm_inv_pos[i] = (inv_pos[i] ) / interval.len;
				norm_inv_neg[i] = (inv_neg[i] ) / interval.len;
			}
			else if (interval.svtype == "DUP" || interval.svtype == "CNV")
			{
				norm_cnv_pos[i] = (cnv_pos[i] ) / avg_isize[i];
				norm_cnv_neg[i] = (cnv_neg[i] ) / avg_isize[i];
				norm_inv_pos[i] = (inv_pos[i] ) / interval.len;
				norm_inv_neg[i] = (inv_neg[i] ) / interval.len;
			}
			else if (interval.svtype == "INV" )
			{
				norm_cnv_pos[i] = (cnv_pos[i] ) / avg_isize[i];
				norm_cnv_neg[i] = (cnv_neg[i] ) / avg_isize[i];
				norm_inv_pos[i] = (inv_pos[i] ) / interval.len;
				norm_inv_neg[i] = (inv_neg[i] ) / interval.len;
			}
			// READLEN fixed to 150 : later!!
			norm_readcount[i] = (double)n_isz[i] * 150.0 / interval.len / avg_depth[i];
		}
	};
	
	void print(sv &interval)
	{
		fprintf(stderr, "Inversion, %d:%d-%d, length %d\n", interval.chrnum, interval.pos, interval.end, interval.len);
		for(int j=0;j<n;++j)
		{
			fprintf(stderr, "%f\t%d,%f\t%d,%f\t%d,%f\t%d,%f\t%d,%f\n", norm_dp[j], n_cnv_pos[j],  norm_cnv_pos[j], n_cnv_neg[j], norm_cnv_neg[j],
				   n_inv_pos[j], norm_inv_pos[j], n_inv_neg[j], norm_inv_neg[j], n_isz[j], isz[j]);
		}
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
	int pos;
	int type; // 0 : insert start, 1: startpos, 2:endpos 3: insert end
	int idx;
	bool operator < (const breakpoint&) const;
	bool operator == (const breakpoint&) const;
	
	breakpoint();
};

void pick_sv_from_merged(vector<sv> &, vector<sv> &);

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
	set<int> *isz_set;
	vector<double> *isz_sum;
	vector<int> *isz_cnt;
//	vector<double> *isz_sum;
//	vector<int> *isz_cnt;
	vector< vector <int> > *isz_list;
//	vector< vector <int> > *pos_list;

	vector< vector <int> > *rev_isz_list;
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
	int initialize(vector<string> &, vector<string> &, vector<double> &, vector<double> &, const char *);
	int read_interval_multi(sv& , svdata &, const char *);	
};

class gcContent
{
public:
	void initialize(string &); // filename for GC content, populate all vectors

	uint16_t num_bin; // Number of GC bin
	uint8_t num_chr; //number of chrs
	uint16_t binsize; // Size of GC bin (bp)
	uint16_t total_bin; // Size of intervals per GC bin

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
	double avg_isize;
	double std_isize;
	double med_isize;
	double avg_rlen;
	
	vector<double> gc_factor;

	// Get GC corrected depth for chr / pos
	double gcCorrected(double, int, int);
	
	void read_depth(vector<sv> &, vector<string> &);
	void process_readpair(sv &, vector<int> &, string &);
	void get_avg_depth();
	void initialize(string &);

	int median(vector<int> &);
	
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

class gtype
{
public:
	double min_bic;
	double p_overlap;
	bool bUseGL;
	bool dFlag;
	
//	void call_genotype(sv &, vector<double>&, vector<double>&, vector<int>&, outvcf&, vector<double>&);
	void call_del(sv &, svdata &, string& ln);
	void call_cnv(sv &, svdata &, string& ln);
	void call_inv(sv &, svdata &, string& ln);

	void eval_cluster(sv &, vector<double> &, vector<double> &, vector<Gaussian2> &, double &, double &);

	int classify_del(vector<double>&, vector<int>&, vector< vector<int> >&, vector<int>&, int&,vector<Gaussian>&, bool);
	int classify_cnv(vector<double>&, vector<int>&, vector<int>&, int&, vector<Gaussian>&);

	void EM(vector<double>&, vector<Gaussian>&, bool);
	void EM(vector<double>&, vector<Gaussian>&);
	void fit(vector<double>&, vector<Gaussian>&);
	void EM2(vector<double>&, vector<double> &, vector<Gaussian2>&);

	int assign(double, vector<Gaussian> &);
	void copyComps(vector<Gaussian> &, vector<Gaussian> &);
	void call_del(sv&, vector<double>&, vector<double>&, vector<int>&, outvcf&, vector<double>&);
	void call_cnv(sv&, vector<double>&, vector<double>&, vector<int>&, outvcf&, vector<double>&);
	void format_output(sv &, string &);
	gtype();
};


template <class T> void vprint(vector<T>);

double MAX(vector<double>&);
double MIN(vector<double>&);
double normpdf(double, Gaussian&);
double mean(vector<double>&);
double stdev(vector<double>&, double);

void read_intervals_from_vcf(vector<string> &, vector<string> &, vector<sv> &);


double RO(interval_t, interval_t);
void merge_svs(vector<sv> &, vector<int> &);
void cluster_svs(vector<sv>&, vector< vector<sv> > &);

double BayesError(vector<Gaussian>&);
double BayesError(vector<Gaussian2>&);

double BIC(vector<double>&, vector<Gaussian>&);
double BIC(vector<double>&, vector<double>&, vector<Gaussian2>&);
double det(double*);

bool ordered(vector<Gaussian>&);

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
