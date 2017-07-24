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

using namespace std;

typedef pair<uint64_t, uint64_t> interval_t;

class Gaussian
{
	public:
		double Mean;
		double Stdev;
		double Alpha;
};

static void split(const char* s, const char* delims, std::vector<std::string>& tokens)
{
	const char* p = s;
	const char* c = p;
	int ndelims = (int)strlen(delims);
	int i;
	tokens.clear();
	
	//fprintf(stderr,"s = '%s', strlen(s) = %d, delims = '%s'\n",s,(int)strlen(s), delims);
	while( *p != '\0' )
	{
		for(i=0; i < ndelims; ++i)
		{
			if ( *p == delims[i] )
				break;
		}
		if ( i != ndelims ) { // delimiter found
			if ( c < p )  { // unless delimiter is consencutive
							//std::string s1(c,p-c);
				tokens.push_back(std::string(c,p-c));
			}
			c = p+1;
		}
		++p;
	}
	if ( c < p ) {
		tokens.push_back(std::string(c,p-c));
	}
}

class sv
{
public:
	string svtype;
	string source;
	int chr;
	int pos;
	int end;
	pair<int,int> ci_pos;
	pair<int,int> ci_end;
	uint64_t len() {return (end - pos + 1);};
	bool operator < (const sv&) const;

		sv();
};

sv::sv()
{
	svtype = "NA";
	source = "NA";
	chr = -1;
	pos = -1;
	end = -1;
	ci_pos.first = 0;
	ci_pos.second = 0;
	ci_end.first = 0;
	ci_end.second = 0;
}

bool sv::operator < (const sv& s) const
{
	if (chr==s.chr)
	{
		if (pos==s.pos)
		{
			return(end<s.end);
		}
		else
		{
			return(pos<s.pos);
		}
	}
	else
	{
		return(chr<s.chr);
	}
}

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

// This function reads a BAM alignment from one BAM file.
static int read_bam(void *data, bam1_t *b) // read level filters better go here to avoid pileup
{
	aux_t *aux = (aux_t*)data; // data in fact is a pointer to an auxiliary structure
	int ret;
	while (1)
	{
		ret = aux->iter? sam_itr_next(aux->fp, aux->iter, b) : sam_read1(aux->fp, aux->hdr, b);
		if ( ret<0 ) break;
		if ( b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP) ) continue;
		if ( (int)b->core.qual < aux->min_mapQ ) continue;
		if ( aux->min_len && bam_cigar2qlen(b->core.n_cigar, bam_get_cigar(b)) < aux->min_len ) continue;
		break;
	}
	return ret;
}


class bfiles
{
public:
	aux_t **data;
	vector<hts_idx_t*> idx;
	int n;

	void read_depth(sv&, vector<double>&);

	void initialize(vector<string>&);
};


class gtype
{
public:
	double min_bic;
	double p_overlap;
	bool bUseGL;
	vector<int> geno;
	void call_genotype(sv &, vector<double>&);
	
	int classify_del(vector<double>&, vector<int>&, vector< vector<int> >&, vector<int>&, int&,vector<Gaussian>&, bool);
	int classify_cnv(vector<double>&, vector<int>&, vector<int>&, int&, vector<Gaussian>&);

	void EM(vector<double>&, vector<Gaussian>&, bool);
	void EM(vector<double>&, vector<Gaussian>&);

	void call_del(vector<double>&);
	void call_cnv(vector<double>&);
	
	void conEM(vector<double>&, vector<double>&, vector<double>&, vector<double>&);
	
	gtype();
};

gtype::gtype()
{
	min_bic = 0;
	bUseGL = false;
	p_overlap = 0;
}

template <class T> void vprint(vector<T>);

double MAX(vector<double>&);
double MIN(vector<double>&);
double normpdf(double, double, double);
double normpdf(double, Gaussian&);
double mean(vector<double>&);
double stdev(vector<double>&, double);


void read_intervals_from_vcf(vector<string> &, vector<string> &, vector<sv> &);


double RO(interval_t, interval_t);
//void clustersvs(vector<interval_t>&, vector<interval_t>&);
void clustersvs(vector<sv>&, vector<sv>&);
void cluster_svs(vector<sv>&, vector< vector<sv> > &);
//void call_deletions(vector<vector<double> > &, vector<double> &, vector<string> &, vector<interval_t>&, FILE*);
//void call_duplications(vector<vector<double> > &, vector<double> &, vector<string> &, vector<interval_t>&, FILE*);


double BayesError(vector<Gaussian>&);

double BIC(vector<double>&, vector<Gaussian>&);
bool ordered(vector<Gaussian>&);

void read_index(string, vector<string>&, vector<string>&, vector<string>&, vector<double>&);

void readIndex(string, vector<string>&, vector<string>&, vector<string>&, vector<string>&);
void readDepthOrig(vector<string>&, vector<string>&, vector<interval_t>&, vector< vector<double> >&, vector<double>&);
//void readDepth(vector<string>&, vector<interval_t>&, vector< vector<double> >&, vector<double>&);
void readDepth(vector<string>&, vector<sv>&, vector< vector<double> >&, vector<double>&);

double getAvgDepth(string, string);
string getCNVsegmentFileName(string, string);
//void readsv(string, int, vector<interval_t>&, vector<interval_t>&);
void readsv(string, int, vector<sv>&, vector<sv>&);

class Data
{
	public:
		ifstream dFile;
		void open(string, unsigned&, unsigned&, unsigned&, vector<string>&, vector<double>&);
//	 	void read(unsigned, unsigned, vector<interval_t>&,  vector< vector<double> > &);
	 	void read(unsigned, unsigned, vector<sv>&,  vector< vector<double> > &);
//		void write(string, vector<string>&, vector<interval_t>&,  vector<interval_t>&, vector< vector<double> >&, vector< vector<double> >&, vector<double>&);
		void write(string, vector<string>&, vector<sv>&,  vector<sv>&, vector< vector<double> >&, vector< vector<double> >&, vector<double>&);
};
//void write_dat(FILE* ,  vector<double>&);
//void write_vcf(FILE*, unsigned, vector<unsigned short>&, vector< vector<unsigned> >&, vector<unsigned>&, pair<uint64_t, uint64_t>&, unsigned, unsigned, vector<double>&, vector<double>&, vector<Gaussian>&, double);
void write_vcf(FILE*, unsigned, vector<unsigned short>&, vector< vector<unsigned> >&, vector<unsigned>&, sv&, unsigned, unsigned, vector<double>&, vector<double>&, vector<Gaussian>&, double, bool);
//void write_vcf_dup(FILE*, unsigned, vector<unsigned short>&, vector<unsigned>&, pair<uint64_t, uint64_t>&, unsigned, unsigned, vector<double>&, vector<double>&, vector<Gaussian>&);
void write_vcf_dup(FILE*, unsigned, vector<unsigned short>&, vector<unsigned>&, sv&, unsigned, unsigned, vector<double>&, vector<double>&, vector<Gaussian>&, bool);
void write_vcfheader(FILE*, vector<string>&);
void printCluster(vector<Gaussian> &C);



#endif // __COMMON_H__
