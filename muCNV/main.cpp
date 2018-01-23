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
	string out_filename;
	string bam_file;
	string sChr;
	string gc_file;
	string sampID;

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
		TCLAP::ValueArg<string> argOut("o","out","Output filename",false,"muCNV.vcf","string");
		TCLAP::ValueArg<string> argVcf("v","vcf","VCF file containing candidate SVs",false,"","string");
		TCLAP::ValueArg<string> argSampleID("s","sample","Sample ID",false,"","string");
		TCLAP::ValueArg<string> argIndex("i","index","List file containing list of intermediate pileups. Required with genotype option",false,"","string");
		TCLAP::ValueArg<string> argGcfile("f","gcFile","File containing GC content information",false, "GRCh38.gc", "string");
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
		cmd.add(argGcfile);
		cmd.add(argIndex);
		cmd.add(argSampleID);
//		cmd.add(switchGenotype);

//		cmd.add(argPos);

//		cmd.add(argBE);
//		cmd.add(argRO);
		
		cmd.parse(argc, argv);
        
		bam_file = argBam.getValue();
		index_file = argIndex.getValue();
		out_filename = argOut.getValue();
		sampID = argSampleID.getValue();
		vcf_file = argVcf.getValue();
		gc_file = argGcfile.getValue();
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
			if (sampID == "")
			{
				cerr << "Error: Sample ID should be supplied with -s option." << endl;
				exit(0);
			}
		}
		
//		P_THRESHOLD = argPos.getValue();
//		BE_THRESHOLD = argBE.getValue();
//		RO_THRESHOLD = argRO.getValue();

	}
	catch (TCLAP::ArgException &e)
	{
		cerr << "Error: " << e.error() << " for arg " << e.argId() << endl;
		abort();
	}
	
	int n = 0;
	if (bGenotype)
	{
		//TODO make read_index to read list files with intermediate results
		read_index(index_file, sample_ids, vcfs, bam_names, avg_depths);
		n = (int)sample_ids.size();
		cerr << "Genotyping index loaded." << endl;
		cerr << n << " samples identified." << endl;
		
		outvcf vfile;

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

		gcContent GC;
		GC.initialize(gc_file);
		
		vector<sv> candidates;
		vcfs.push_back(vcf_file);
		read_intervals_from_vcf(sample_ids, vcfs, candidates);
		
		cerr<< candidates.size() << " intervals identified from the VCF file." << endl;

		vector<int> idxs;

		// In sample-by-sample process, assume overlapping SVs are already clustered into a single event, so do simple merging without clustering
		merge_svs(candidates, idxs);

		bFile b(GC);
		b.initialize(bam_file);
		cerr << "BAM file initialized" << endl;
		b.get_avg_depth();
		
		FILE *fp = fopen(out_filename.c_str(), "wt");
		fprintf(fp, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n", sampID.c_str());
		fprintf(fp, "0\t0\t.\t.\t.\t.\t.\tAVGDP;AVG_ISIZE;STD_ISIZE;MED_ISIZE\tDP:AI:SI:MI\t%.1f:%.1f:%.1f:%.1f\n", b.avg_dp, b.avg_isize,b.std_isize,b.med_isize);

		for(int i=0; i<(int)idxs.size(); ++i)
		{
			// 4. Genotype for each variant
			int last_idx;
			if (i<idxs.size()-1)
			{
				last_idx = idxs[i+1];
			}
			else
			{
				last_idx = candidates.size();
			}
			vector<sv> svlist(candidates.begin()+idxs[i], candidates.begin() + last_idx)  ;

		//	cerr << "processing SV block for " << svlist.size() << " intervals" << endl;

	//		for(int j=0;j<svlist.size(); ++j)
	//		{
	//			cerr << svlist[j].chr <<  ":" << svlist[j].pos << "-" << svlist[j].end << "\t" << svlist[j].svtype << endl;
	//		}

			int m = (int)svlist.size();
			/*
			vector<double> X(m,0);
			vector<double> GX(m,0);
			vector< vector<int> > ISZ;
			ISZ.resize(m);
			for(int j=0;j<m;++j)
			{
				ISZ[j].resize(3);
				ISZ[j][0] = 0;
				ISZ[j][1] = 0;
				ISZ[j][2] = 0;
			}
			*/
			
			vector<string> G (m, "");
		//	cerr << "reading depth... " ;
			b.read_depth(svlist, G);
		//	cerr << "done " << endl;
			
			for(int j=0;j<m;++j)
			{
				fprintf(fp, "%s\t%d\t.\t.\t.\t.\t.\tEND=%d;SVTYPE=%s\t.", svlist[j].chr.c_str(), svlist[j].pos,svlist[j].end,svlist[j].svtype.c_str());
				fprintf(fp, "\t%s\n", G[j].c_str());
//				fprintf(fp, "\t%.1f:%.1f:%d,%d,%d\n",X[j],GX[j],ISZ[j][0],ISZ[j][1],ISZ[j][2]);
			}

		}
		fclose(fp);
		cerr << "Finished without an error" << endl;
	}
	return 0;
}
