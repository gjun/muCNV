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

bool bUseGL = true;
double P_THRESHOLD = 0.5;
double BE_THRESHOLD = 0.05;
double RO_THRESHOLD = 0.5;


int main(int argc, char** argv)
{
	cerr << "muCNV 0.6-- Multi-sample CNV genotyper" << endl;
	cerr << "(c) 2017 Goo Jun" << endl << endl;
	cerr.setf(ios::showpoint);

	bool bGenotype = false;
//	bool bVerbose;
	string index_file;
	string vcf_file;
	string out_prefix;
	string bam_file;
	string sChr;

	vector<string> sample_ids;
	vector<string> vcfs;
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
		TCLAP::CmdLine cmd("Command description message", ' ', "0.06");
//		TCLAP::ValueArg<string> argIn("i","index","Input index file (sample ID, candidate VCF, BAM/CRAM)",true,"","string");
		
		TCLAP::ValueArg<string> argBam("b","bam","Input BAM/CRAM file name",false,"","string");
		TCLAP::ValueArg<string> argOut("o","out","Prefix for output filename",false,"muCNV","string");
		TCLAP::ValueArg<string> argVcf("v","vcf","VCF file containing candidate SVs",false,"","string");
		TCLAP::ValueArg<string> argIndex("i","index","List file containing list of intermediate pileups. Required with genotype option",false,"","string");
		TCLAP::SwitchArg switchGenotype("g","genotype","Generate Genotype from intermediate pileups", cmd, false);
		
		//		TCLAP::ValueArg<string> argsv("n","interval","File containing list of candidate intervals",false,"","string");
//		TCLAP::ValueArg<double> argPos("p","posterior","(Optional) Posterior probability threshold",false,0.5,"double");
//		TCLAP::ValueArg<double> argBE("b","bayes","(Optional) Bayes error threshold",false,0.2,"double");

//		TCLAP::ValueArg<double> argRO("r","reciprocal","(Optional) Reciprocal overlap threshold to merge candidate intervals (default: 0.5)",false,0.5,"double");
		//		TCLAP::SwitchArg switchGL("g","gl","Use likelihood instead of posterior to call",cmd,false);
//		TCLAP::SwitchArg switchVerbose("v","verbose","Turn on verbose mode",cmd,false);

		cmd.add(argBam);
		cmd.add(argOut);
		cmd.add(argVcf);
		cmd.add(argIndex);
//		cmd.add(switchGenotype);

//		cmd.add(argPos);

//		cmd.add(argBE);
//		cmd.add(argRO);
		
		cmd.parse(argc, argv);
        
		bam_file = argBam.getValue();
		index_file = argIndex.getValue();
		out_prefix = argOut.getValue();
		vcf_file = argVcf.getValue();
		bGenotype = switchGenotype.getValue();
		
		if (bGenotype && index_file == "")
		{
			cerr << "Error: list file is required for genotyping" << endl;
			exit(0);
		}
		else if (bGenotype == false)
		{
			if (bam_file == "")
			{
				cerr << "Error: BAM/CRAM file is required for individual processing mode." << endl;
				exit(0);
			}
			if (vcf_file == "")
			{
				cerr << "Error: VCF file with SV events is required for individual processing mode." << endl;
				exit(0);
			}
		}
		
//		P_THRESHOLD = argPos.getValue();
//		BE_THRESHOLD = argBE.getValue();
//		RO_THRESHOLD = argRO.getValue();

      //  bVerbose = switchVerbose.getValue();

	}
	catch (TCLAP::ArgException &e)
	{
		cerr << "Error: " << e.error() << " for arg " << e.argId() << endl;
		abort();
	}
	
	// 0. Read (vcf, bam/cram) file list
	// 0.0. Calculate average sequencing depth for each BAM
	// 0.1. Calcualte GC-content statistics here

//	read_index(listrea_file, sample_ids, vcf_files, bam_names, avg_depths);
	
	int n = 0;
	if (bGenotype)
	{
		//TODO make read_index to read list files with intermediate results
		read_index(index_file, sample_ids, vcfs, bam_names, avg_depths);
		n = (int)sample_ids.size();
		cerr << "Genotyping index loaded." << endl;
		cerr << n << " samples identified." << endl;
		
		outvcf vfile;
		string out_filename = out_prefix + ".vcf";

		vfile.open(out_filename);
		
		vfile.write_header(sample_ids);
		
		vfiles V;
		V.initialize(vcfs);
		
		sv interval;
		vector<double> X (n, 0);
		vector<double> Y (n, 0);
		vector<int> G(n, 0);

		while(V.read_interval(interval, X)>=0)
		{
			gtype g;
			
			for(int k=0; k<n; ++k)
			{
				X[k] = X[k] / avg_depths[k];
			}
			g.call_genotype(interval, X, Y, G, vfile, avg_depths);
		}
		vfile.close();
	}
	else
	{
		n = 1;
		cerr << "Processing individual file." << endl;


	//	cerr<< "index loaded" << endl;
		
		// 1. Read intervals from individual VCFs
		vector<sv> candidates;
		vcfs.push_back(vcf_file);
		read_intervals_from_vcf(sample_ids, vcfs, candidates);
		
		cerr<< candidates.size() << " intervals identified from the VCF file." << endl;

		// 2. Merge intervals with RO > minRO  (store all original informattion - merged intervals will be vector of intervals)

		vector< vector<sv> > merged_candidates;
	//	cluster_svs(candidates, merged_candidates);
		// In sample-by-sample process, assume overlapping SVs are already clustered into a single event, so do simple merging without clustering
		merge_svs(candidates, merged_candidates);
		
	//	cerr<< merged_candidates.size() << " sets of intervals after clustering." << endl;
	//	bfiles bf;

		// Open BAM file handles (create a class for bamfiles, method to read specific interval)
		//bf.initialize(bam_names);

		
		bFile b;
		b.initialize(bam_file);
		
		//	bf.get_avg_depth(avg_depths); // maybe make this into a separate program // samtools flagstat works good enough for overall depth, but not for GC content
		
		// 3. For each interval, read depth  (+ read pair distance ?) from individual BAM/CRAM file

		string out_filename = out_prefix + ".vcf";
	//	outvcf vfile;
	//	vfile.open(out_filename);
	//	vfile.write_header(sample_ids);
		
	//	FILE *fp = fopen(out_filename.c_str(), "wt");
	//	fprintf(fp, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n", out_prefix.c_str());

		for(int i=0; i<(int)merged_candidates.size(); ++i)
		{
			// 4. Genotype for each variant
			vector<sv> &svlist = merged_candidates[i];

	//		pick_sv_from_merged(svlist, merged_candidates[i]); // Pick one interval (median) from merged candidates, svlist will have only one element
	//		cerr << "this merged list contains " << svlist.size() << " items" << endl;

			int m = (int)svlist.size();
			vector<double> X(m,0);
			vector<double> GX(m,0);
			vector<double> ISZ(m,0);
			
			b.read_depth(svlist, X, GX, ISZ);
			
			for(int j=0;j<m;++j)
			{
				// fprintf(fp, "%s\t%d\t.\t.\t<%s>\t.\t.\tEND=%d;SVTYPE=%s;SVLEN=%d\tDP", svlist[j].chr.c_str(), svlist[j].pos,svlist[j].svtype.c_str(), svlist[j].end,svlist[j].svtype.c_str(),svlist[j].end - svlist[j].pos);
				// fprintf(fp,"\t%f\n",X[j],GX[j],ISZ[j]);
				printf("%s\t%d\t.\t.\t<%s>\t.\t.\tEND=%d;SVTYPE=%s;SVLEN=%d\tDP", svlist[j].chr.c_str(), svlist[j].pos,svlist[j].svtype.c_str(), svlist[j].end,svlist[j].svtype.c_str(),svlist[j].end - svlist[j].pos);
				printf("\t%.1f:%.1f:%.1f\n",X[j],GX[j],ISZ[j]);
			}

		}
//		fclose(fp);

	//vfile.close();
	}
	return 0;
}
