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
#include "pFile.h"

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

bool SampleList::readIndex(string sInFile)
{
    ifstream inFile(sInFile.c_str(), ios::in);
    
    while(inFile.good())
    {
        string ln;
        getline(inFile,ln);
        if (!ln.empty())
        {
            vector<string> tokens;
            pFile::tokenizeLine(ln.c_str(), " \t\n", tokens);
            
            if (tokens[0].empty())
            {
                cerr << "Error: empty sample ID"<< endl;
                abort();
            }
            else if (tokens[1].empty())
            {
                cerr << "Error: depth file is empty" << endl;
            }
            
            vID.push_back(tokens[0]);
            vDepthFile.push_back(tokens[1]);
            
            
            // To Do: do this earlier and get vector of *pFile instead of file names
            pFile dFile;
            dFile.load(tokens[1].c_str(), NULL, true); // read the file with header
            const char *line;
            
            line = dFile.getLine();
            
            tokens.clear();
            pFile::tokenizeLine(line,"#= \t\n",tokens);
            // TMPTMP
            if (tokens[1].compare("gcCorrectedCvg") == 0)
            {
                AvgDepth.push_back(atof(tokens[2].c_str()));
            }
            else
            {
                //header is in wrong format
            }

            inFile.close();
            // NOTE:: This should be in the script
            /*
            double getAvgDepth(string smID, string asmDir)
            {
                // To do: check this file naming convention
                string sDepthFile = asmDir + "/CNV/depthOfCoverage_100000-" + smID + ".tsv";
                
                ifstream inFile(sDepthFile.c_str(), ios::in);
                int gcIdx = 4;
                double sum = 0;
                double N = 0;
                
                if (!inFile.good())
                {
                    cerr << "Cannot open  " << sDepthFile<< endl;
                    abort();
                }
                while(inFile.good())
                {
                    string ln;
                    getline(inFile,ln);
                    if (!ln.empty())
                    {
                        vector<string> tokens;
                        pFile::tokenizeLine(ln.c_str(), " \t\n", tokens);
                        if (tokens[0].at(0) == '>')
                        {
                            for(unsigned j=0; j<tokens.size() ;++j)
                            {
                                if (tokens[j].compare("gcCorrectedCvg") == 0)
                                {
                                    gcIdx = j;
                                }
                            }
                        }
                        else if (tokens[0].substr(0,3).compare("chr") == 0)
                        {
                            sum += atof(tokens[gcIdx].c_str());
                            N++;
                        }
                    }
                }
                inFile.close();
                
                return(sum/(double)N);
            }
*/
            
            
        }
    }
    inFile.close();
    return true;
}


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
