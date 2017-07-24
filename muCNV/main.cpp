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
	cout << "muCNV 0.5 -- Multi-sample CNV genotyper" << endl;
	cout << "(c) 2017 Goo Jun" << endl << endl;
	cerr.setf(ios::showpoint);

	bool bVerbose;
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

        bVerbose = switchVerbose.getValue();

	}
	catch (TCLAP::ArgException &e)
	{
		cerr << "Error: " << e.error() << " for arg " << e.argId() << endl;
		abort();
	}
	
	// 0. Read (vcf, bam/cram) file list
	// 0.0 Calculate average sequencing depth for each BAM -- make this a pre-processing step, write a Python script, to save redundant calculations -- or just do it here?
	
	read_index(index_file, sample_ids, vcf_files, bam_names, avg_depths);
	int n = sample_ids.size();
	
	// 1. Read intervals from individual VCFs
	vector<sv> candidates;
	read_intervals_from_vcf(sample_ids, vcf_files, candidates);
	
	
	// 2. Merge intervals with RO > minRO  (store all original informattion - merged intervals will be vector of intervals)
	vector< vector<sv> > merged_candidates;
	cluster_svs(candidates, merged_candidates);
	
	bfiles bf;
	// Open BAM file handles (create a class for bamfiles, method to read specific interval)
	bf.initialize(bam_names);
	
	// 3. For each interval, read depth  (+ read pair distance ?) from individual BAM/CRAM file

	for(int i=0; i<merged_candidates.size(); ++i)
	{
		// 4. Genotype for each variant
		vector<sv> &svlist = merged_candidates[i];
		
		for(int j=0; j<svlist.size(); ++j)
		{
			gtype g;
			vector<double> X(n, 0);
			bf.read_depth(svlist[j], X);
			g.call_genotype(svlist[j], X);
			
			cout << svlist[j].chr << "\t" << svlist[j].pos << "\t" << svlist[j].end;
			for(int k=0; k<X.size(); ++k)
			{
				cout << "\t" << g.geno[k];
			}
			cout << endl;

		}
		

	}
	
	// 5. Output
	

	// string fname = "/Users/gjun/data/cram/NA12878.fragment.cram";
	
	//while(sam_itr_next(myFile, iter, aln)>=0)
	{
		//	cout << myHeader->target_name[aln->core.tid];
		//	cout << aln->core.pos << endl;
		//count++;
	}
	
	//	cout << "Number of reads : " << count << endl;

/*

	unsigned n_sample = 0;
	unsigned n_del = 0;
	unsigned n_dup = 0;

	vector<string> depthFiles;

	vector<sv> del_events;
	vector<sv> dup_events;

	vector<string> eventFiles;

    // Read sample index
    SampleList samples;
    samples.readIndex(sInFile);
    

    
//	n_sample = (unsigned) sampleIDs.size();
    
    if (bVerbose) cerr << samples.n_sample << " samples identified\n"<<endl;
  
	AvgDepth.resize(n_sample, 0);

	for(unsigned i=0; i<n_sample;++i)
	{
		// Read average depth per sample
		AvgDepth[i] = getAvgDepth(sampleIDs[i], sampleDirs[i]);
		if (bVerbose) cerr << "\rSample " << i+1 << " / "<< n_sample << ": " << sampleIDs[i] << ", average Depth " << AvgDepth[i];

		// Read CNV files
		string sSegmentFile = getCNVsegmentFileName(sampleIDs[i], sampleDirs[i]);
		readsv(sSegmentFile, 1, del_events, dup_events);

		// Read event files
		readsv(eventFiles[i], 2, del_events, dup_events);
	}

	if (bVerbose) 
	{
		cerr << endl << del_events.size() << " deletion and ";
		cerr << dup_events.size() << " duplication candidate intervals identified" <<  endl;
	}

	if (ssvFile != "")
	{
		readsv(ssvFile, 3, del_events, dup_events);
	}

	// Sort intervals, remove duplicates, and cluster overlapping intervals
	sort(del_events.begin(), del_events.end(), comparesvs);

	// To do: remove duplicates before clustering -- will make things faster, but results should be the same
	//		remove_duplicates(del_events);
	clustersvs(del_events, del_intervals);

	sort(dup_events.begin(), dup_events.end(), comparesvs);
	//		remove_duplicates(dup_events);
	clustersvs(dup_events, dup_intervals);

	// merge dup_intervals and del_intervals 
	all_intervals = del_intervals;
	uint32_t n_sv = (uint32_t) all_intervals.size();

//	n_del = del_intervals.size();
//	n_dup = dup_intervals.size();

	// Store average depth per interval per sample
	vector< vector<double> > X(n_del, vector<double>(n_sample,0));   
	vector< vector<double> > Y(n_dup, vector<double>(n_sample,0));   

	// Editing: read full length from one interval at a time, and do breakpoint refinement // Breakpoint refinement will be implemented in later version.
	//
//
//   Read  2 x interval (0.5 before, 0.5 after) from all samples
//   Do clustering on the original or smallest (? confident ) interval
//   If passed BIC and initial Bayes error criteria, Refine breakpoint
//   Normalize each read depth by average depth across samples ( how to adjust for allele count ? )
//   Re-do clustering, evaluate Bayes error

//	readDepth(depthFiles, del_intervals, X, AvgDepth);
//	readDepth(depthFiles, dup_intervals, Y, AvgDepth);

	uint64_t max_svlen = 0;
	for(unsigned i=0;i<n_sv;++i)
	{
		uint64_t L = all_intervals[i].len();
		// Temporary, for quick run
		if (L>1e6)
		{
			all_intervals.erase(all_intervals.begin()+i);
		}
		else if (max_svlen < L )
		{
			max_svlen = L;
		}
	}

	// Reserve memory space to read all intervals -- up to 2 * max_svlen for each sample
	vector< vector<double> > D(n_sample, vector<double> (1,0));
	vector<double>  AvgD(n_sample, 0);

	if (bVerbose)
	{
		cerr << n_sv << "intervals with max SVLEN " << max_svlen << endl;
	}

	for(unsigned i=0;i<n_sample;++i)
	{
		D[i].resize(max_svlen + 2000);
	}

	string sVcfFile = sOutPrefix + ".vcf" ;

	FILE *vcfFile;
	vcfFile = fopen(sVcfFile.c_str(), "wt");

	write_vcfheader(vcfFile, sampleIDs);

	// To Do : make a class for VCF writer, and add methods to write header and SV entries

	for(unsigned i=0;i<n_sv;++i)
	{
		// Read depth for i-th interval
	//TMPTMP	readDepth(depthFiles, all_intervals[i],  D, AvgD, AvgDepth);

		// Do clustering 
		if (all_intervals[i].sv_type == "DEL")
		{
//TMPTMP			call_deletion(all_intervals[i], AvgD, AvgDepth, vcfFile);
		}
		else
		{
//			call_duplication(all_intervals[i], vcfFile);
//		Test 
		}

		// Refine breakpoints
		
		// Re-do clustering

		// Write to VCF
	}

	call_deletions(X, AvgDepth, sampleIDs,  del_intervals, vcfFile);
	call_duplications(Y, AvgDepth, sampleIDs, dup_intervals, vcfFile);

	fclose(vcfFile);
*/
	return 0;
}
