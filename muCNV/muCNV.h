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

static void tokenizeLine(const char* s, const char* delims, std::vector<std::string>& tokens)
{
	const char* p = s;
	const char* c = p;
	int ndelims = strlen(delims);
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

class Interval
{
	public:
		interval_t T;
		string sv_type;
		string source;
		pair<int,int> ci_pos;
		pair<int,int> ci_end;
		uint64_t len() {return (T.second - T.first + 1);};
};

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

bool compareIntervals(const Interval &, const Interval &);

template <class T> void vprint(vector<T>);

double MAX(vector<double>&);
double MIN(vector<double>&);
double normpdf(double, double, double);
double normpdf(double, Gaussian&);
double mean(vector<double>&);
double stdev(vector<double>&, double);

double RO(interval_t, interval_t);
//void clusterIntervals(vector<interval_t>&, vector<interval_t>&);
void clusterIntervals(vector<Interval>&, vector<Interval>&);

//void call_deletions(vector<vector<double> > &, vector<double> &, vector<string> &, vector<interval_t>&, FILE*);
//void call_duplications(vector<vector<double> > &, vector<double> &, vector<string> &, vector<interval_t>&, FILE*);

void call_deletions(vector<vector<double> > &, vector<double> &, vector<string> &, vector<Interval>&, FILE*);
void call_duplications(vector<vector<double> > &, vector<double> &, vector<string> &, vector<Interval>&, FILE*);

void EM(vector<double>&, vector<Gaussian>&);
void EM(vector<double>&, vector<Gaussian>&, bool);
void conEM(vector<double>&, vector<double>&, vector<double>&, vector<double>&);
double BayesError(vector<Gaussian>&);

unsigned classify(vector<double>&, vector<unsigned short>&, vector< vector<unsigned> >&, vector<unsigned>&, unsigned&,vector<Gaussian>&, bool);
unsigned classify(vector<double>&, vector<unsigned short>&, vector<unsigned>&, unsigned&,vector<Gaussian>&);
double BIC(vector<double>&, vector<Gaussian>&);
bool ordered(vector<Gaussian>&);

void readIndex(string, vector<string>&, vector<string>&, vector<string>&, vector<string>&);
void readDepthOrig(vector<string>&, vector<string>&, vector<interval_t>&, vector< vector<double> >&, vector<double>&);
//void readDepth(vector<string>&, vector<interval_t>&, vector< vector<double> >&, vector<double>&);
void readDepth(vector<string>&, vector<Interval>&, vector< vector<double> >&, vector<double>&);

double getAvgDepth(string, string);
string getCNVsegmentFileName(string, string);
//void readInterval(string, int, vector<interval_t>&, vector<interval_t>&);
void readInterval(string, int, vector<Interval>&, vector<Interval>&);

class Data
{
	public:
		ifstream dFile;
		void open(string, unsigned&, unsigned&, unsigned&, vector<string>&, vector<double>&);
//	 	void read(unsigned, unsigned, vector<interval_t>&,  vector< vector<double> > &);
	 	void read(unsigned, unsigned, vector<Interval>&,  vector< vector<double> > &);
//		void write(string, vector<string>&, vector<interval_t>&,  vector<interval_t>&, vector< vector<double> >&, vector< vector<double> >&, vector<double>&);
		void write(string, vector<string>&, vector<Interval>&,  vector<Interval>&, vector< vector<double> >&, vector< vector<double> >&, vector<double>&);
};
//void write_dat(FILE* ,  vector<double>&);
//void write_vcf(FILE*, unsigned, vector<unsigned short>&, vector< vector<unsigned> >&, vector<unsigned>&, pair<uint64_t, uint64_t>&, unsigned, unsigned, vector<double>&, vector<double>&, vector<Gaussian>&, double);
void write_vcf(FILE*, unsigned, vector<unsigned short>&, vector< vector<unsigned> >&, vector<unsigned>&, Interval&, unsigned, unsigned, vector<double>&, vector<double>&, vector<Gaussian>&, double, bool);
//void write_vcf_dup(FILE*, unsigned, vector<unsigned short>&, vector<unsigned>&, pair<uint64_t, uint64_t>&, unsigned, unsigned, vector<double>&, vector<double>&, vector<Gaussian>&);
void write_vcf_dup(FILE*, unsigned, vector<unsigned short>&, vector<unsigned>&, Interval&, unsigned, unsigned, vector<double>&, vector<double>&, vector<Gaussian>&, bool);
void write_vcfheader(FILE*, vector<string>&);
void printCluster(vector<Gaussian> &C);


#endif // __COMMON_H__
