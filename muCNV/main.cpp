/*
 *   Author: Goo Jun (goo.jun@uth.tmc.edu)
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <algorithm>

// #include "pFile.h"
#include "muCNV.h"

// TCLAP headers
#include "tclap/CmdLine.h"
#include "tclap/Arg.h"


using namespace std;

bool bUseGL = false;
double P_THRESHOLD = 0.9;
double BE_THRESHOLD = 0.01;
double RO_THRESHOLD = 0.8;


int main(int argc, char** argv)
{
	cerr << "muCNV 0.5 -- Multi-sample CNV genotyper" << endl;
	cerr << "(c) 2017 Goo Jun" << endl << endl;
	cerr.setf(ios::showpoint);

//	bool bVerbose;
	string index_file;

	string out_prefix;
	string sChr;

	vector<string> sample_ids;
	vector<string> vcf_files;
	vector<string> bam_names;
	vector<double> avg_depths;
	map<string, unsigned> hIdSex;

	vector<double> AvgDepth;
	//	vector<interval_t> del_intervals;
	//	vector<interval_t> dup_intervals;
	//	vector<sv> del_intervals;
	// vector<sv> dup_intervals;
	// vector<sv> all_intervals;

	srand((unsigned int)time(NULL));

	// Parsing command-line arguments
	try 
	{
		TCLAP::CmdLine cmd("Command description message", ' ', "0.01");
		TCLAP::ValueArg<string> argIn("i","index","Input index file (sample ID, candidate VCF, BAM/CRAM)",true,"","string");
		TCLAP::ValueArg<string> argOut("o","out","Prefix for output filename",false,"muCNV","string");
		
		//		TCLAP::ValueArg<string> argsv("n","interval","File containing list of candidate intervals",false,"","string");
		TCLAP::ValueArg<double> argPos("p","posterior","(Optional) Posterior probability threshold",false,0.9,"double");
		TCLAP::ValueArg<double> argBE("b","bayes","(Optional) Bayes error threshold",false,0.01,"double");

		TCLAP::ValueArg<double> argRO("r","reciprocal","(Optional) Reciprocal overlap threshold to merge candidate intervals (default: 0.8)",false,0.8,"double");
		//		TCLAP::SwitchArg switchGL("g","gl","Use likelihood instead of posterior to call",cmd,false);
		TCLAP::SwitchArg switchVerbose("v","verbose","Turn on verbose mode",cmd,false);

		cmd.add(argIn);
		cmd.add(argOut);
		cmd.add(argPos);

		//		cmd.add(argsv);
		cmd.add(argBE);
		cmd.add(argRO);
		cmd.parse(argc, argv);

        // Index file:
		// three columns (tab separated)
		// sample ID, VCF with canddiate intervals, sequence file (BAM/CRAM)
        // SampleID, DepthFile(bgzipped, tabixed)
		
        
		index_file = argIn.getValue();
		//ssvFile = argsv.getValue();
		out_prefix = argOut.getValue();

		P_THRESHOLD = argPos.getValue();
		BE_THRESHOLD = argBE.getValue();
		RO_THRESHOLD = argRO.getValue();

      //  bVerbose = switchVerbose.getValue();

	}
	catch (TCLAP::ArgException &e)
	{
		cerr << "Error: " << e.error() << " for arg " << e.argId() << endl;
		abort();
	}
	
	// 0. Read (vcf, bam/cram) file list
	// 0.0 Calculate average sequencing depth for each BAM -- make this a pre-processing step, write a Python script, to save redundant calculations -- or just do it here?
	
	read_index(index_file, sample_ids, vcf_files, bam_names, avg_depths);
	int n = (int)sample_ids.size();

	cerr<< "index loaded" << endl;
	
	// 1. Read intervals from individual VCFs
	vector<sv> candidates;
	read_intervals_from_vcf(sample_ids, vcf_files, candidates);
	
	cerr<< candidates.size() << " intervals read" << endl;
	
	// 2. Merge intervals with RO > minRO  (store all original informattion - merged intervals will be vector of intervals)
	vector< vector<sv> > merged_candidates;
	cluster_svs(candidates, merged_candidates);

	cerr<< merged_candidates.size() << " intervals clustered " << endl;
	bfiles bf;
	// Open BAM file handles (create a class for bamfiles, method to read specific interval)
	bf.initialize(bam_names);
	//	bf.get_avg_depth(avg_depths); // maybe make this into a separate program
	
	// 3. For each interval, read depth  (+ read pair distance ?) from individual BAM/CRAM file

	string vcf_filename = out_prefix + ".vcf";
	outvcf vfile;
	vfile.open(vcf_filename);
	
	vfile.write_header(sample_ids);
	int cnt = 0;

	for(int i=0; i<(int)merged_candidates.size(); ++i)
	{
		// 4. Genotype for each variant
		vector<sv> &svlist = merged_candidates[i];
//		cerr << "this merged list contains " << svlist.size() << " items" << endl;
// TODO : pick one interval (median) from merged candidates
		
		for(int j=0; j<(int)svlist.size(); ++j)
		{
			cnt++;
			gtype g;
			vector<double> X(n, 0);
			vector<double> Y(n, 0);
			vector<int> G(n, 0);

			cerr << svlist[j].chr << "\t" << svlist[j].pos << "\t" << svlist[j].end << endl;

			bf.read_depth(svlist[j], X);
			// normalize
			for(int k=0; k<n; ++k)
			{
				X[k] = X[k] / avg_depths[k];
			}
			g.call_genotype(svlist[j], X, Y, G, vfile, avg_depths);
		}
	}
	
	vfile.close();

	return 0;
}
