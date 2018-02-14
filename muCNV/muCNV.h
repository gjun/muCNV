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

#include "sam.h"

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
//	pair<int,int> ci_pos;
//	pair<int,int> ci_end;
	uint64_t len() {return (end - pos + 1);};
	bool operator < (const sv&) const;
	bool operator == (const sv&) const;
	
	sv();
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
	double logpdf(const double&, const double&);
	void update(); // update precision matrix
	
	Gaussian2();
};

class invcfs
{
public:
	vector<ifstream *> vfs;
	vector<int> id_offset;
	vector<int> n_id;

	int initialize(vector<string> &, vector<string> &, vector<double> &, vector<double> &);
	int read_interval_multi(sv& , vector<double> &, vector<double> &, vector<double> &, vector<double> &, vector<double> &);
	int read_interval(sv&, vector<double> &);
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
//	int n;
	gcContent& GC;
	double avg_dp;
	double avg_isize;
	double std_isize;
	double med_isize;
	double avg_rlen;
	
	vector<double> gc_factor;

	// Get GC corrected depth for chr / pos
	double gcCorrected(double, int, int);
	
// Handle multiple, overlapping SVs
	//void read_depth(vector<sv> &, vector<double>&, vector<double>&, vector< vector<int> >&);
	void read_depth(vector<sv> &, vector<string> &);
	//void process_readpair(sv &, vector<int> &, vector<int> &, string &);
	void process_readpair(sv &, vector<int> &, string &);
	// average depth, average gc-corrected depth, average insert size // stdev?
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
	void call_genotype(sv &, vector<double>&, vector<double>&, vector<int>&, outvcf&, vector<double>&);
	
	int classify_del(vector<double>&, vector<int>&, vector< vector<int> >&, vector<int>&, int&,vector<Gaussian>&, bool);
	int classify_cnv(vector<double>&, vector<int>&, vector<int>&, int&, vector<Gaussian>&);

	void EM(vector<double>&, vector<Gaussian>&, bool);
	void EM(vector<double>&, vector<Gaussian>&);

	void call_del(sv&, vector<double>&, vector<double>&, vector<int>&, outvcf&, vector<double>&);
	void call_cnv(sv&, vector<double>&, vector<double>&, vector<int>&, outvcf&, vector<double>&);
	
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

double BIC(vector<double>&, vector<Gaussian>&);
double BIC2(vector<double>&, vector<double>&, vector<Gaussian2>&);
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
