#ifndef __COMMON_H__
#define __COMMON_H__

#include <stdint.h>
#include <stdio.h>
#include <map>
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
	string source;
	string chr;
	int chrnum;
	int pos;
	int end;
	pair<int,int> ci_pos;
	pair<int,int> ci_end;
	uint64_t len() {return (end - pos + 1);};
	bool operator < (const sv&) const;
	bool operator == (const sv&) const;
	
	sv();
};

class breakpoint
{
public:
	int pos;
	bool type; // 0: start, 1: end
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


class bfiles
{
public:
	aux_t **data;
	vector<hts_idx_t*> idx;
	int n;

	void read_depth(sv&, vector<double>&);
	void get_avg_depth(vector<double>&);
	void initialize(vector<string>&);
};

// A class for a single file BAM I/O
class bFile
{
public:
	aux_t *data;
	hts_idx_t* idx;
//	int n;
	vector<double> gcbin;
	
// Handle multiple, overlapping SVs
	void read_depth(vector<sv> &, vector<double>&, vector<double>&);
	double read_pair(sv &);
	// average depth, average gc-corrected depth, average insert size // stdev?
	void get_avg_depth(double &, double&, double&);
	
	void initialize(string&);
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
void merge_svs(vector<sv> &, vector< vector<sv> > &);
void cluster_svs(vector<sv>&, vector< vector<sv> > &);

double BayesError(vector<Gaussian>&);

double BIC(vector<double>&, vector<Gaussian>&);
double BIC2(vector<double>&, vector<double>&, vector<Gaussian2>&);
double det(double*);

bool ordered(vector<Gaussian>&);

void read_index(string, vector<string>&, vector<string>&, vector<string>&, vector<double>&);

void readIndex(string, vector<string>&, vector<string>&, vector<string>&, vector<string>&);
void readDepthOrig(vector<string>&, vector<string>&, vector<interval_t>&, vector< vector<double> >&, vector<double>&);
void readDepth(vector<string>&, vector<sv>&, vector< vector<double> >&, vector<double>&);

double getAvgDepth(string, string);
string getCNVsegmentFileName(string, string);
void readsv(string, int, vector<sv>&, vector<sv>&);


void printCluster(vector<Gaussian> &C);



#endif // __COMMON_H__
