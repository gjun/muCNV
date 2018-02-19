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
	cerr << "muCNV 0.8 -- Multi-sample CNV genotyper" << endl;
	cerr << "(c) 2018 Goo Jun" << endl << endl;
	cerr.setf(ios::showpoint);

	bool bGenotype = false;
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
	map<string, unsigned> hIdSex;


	srand((unsigned int)time(NULL));

	// Parsing command-line arguments
	try 
	{
		TCLAP::CmdLine cmd("Command description message", ' ', "0.06");
		
		TCLAP::ValueArg<string> argBam("b","bam","Input BAM/CRAM file name",false,"","string");
		TCLAP::ValueArg<string> argOut("o","out","Output filename",false,"muCNV.vcf","string");
		TCLAP::ValueArg<string> argVcf("v","vcf","VCF file containing candidate SVs",false,"","string");
		TCLAP::ValueArg<string> argSampleID("s","sample","Sample ID",false,"","string");
		TCLAP::ValueArg<string> argIndex("i","index","List file containing list of intermediate pileups. Required with genotype option",false,"","string");
		TCLAP::ValueArg<string> argGcfile("f","gcFile","File containing GC content information",false, "GRCh38.gc", "string");
		TCLAP::SwitchArg switchGenotype("g","genotype","Generate Genotype from intermediate pileups", cmd, false);
		

		cmd.add(argBam);
		cmd.add(argOut);
		cmd.add(argVcf);
		cmd.add(argGcfile);
		cmd.add(argIndex);
		cmd.add(argSampleID);
		
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
			exit(1);
		}
		else if (bGenotype == false)
		{
			if (bam_file == "")
			{
				cerr << "Error: BAM/CRAM file is required for individual processing mode." << endl;
				exit(1);
			}
			if (vcf_file == "")
			{
				cerr << "Error: VCF file with SV events is required for individual processing mode." << endl;
				exit(1);
			}
			if (sampID == "")
			{
				cerr << "Error: Sample ID should be supplied with -s option." << endl;
				exit(1);
			}
		}
	}
	catch (TCLAP::ArgException &e)
	{
		cerr << "Error: " << e.error() << " for arg " << e.argId() << endl;
		abort();
	}
	
	int n = 0;
	if (bGenotype)
	{
		// Multi-sample genotyping from summary VCFs
		
		vector<double> avg_depths;
		vector<double> avg_isizes;
		
		invcfs V_list;
		read_vcf_list(index_file, vcfs);
		V_list.initialize(vcfs, sample_ids, avg_depths, avg_isizes);

		n = (int)sample_ids.size();
		cerr << "Genotyping index loaded." << endl;
		cerr << n << " samples identified." << endl;
	
		outvcf vfile;

		vfile.open(out_filename);
		
		vfile.write_header(sample_ids);
		
		sv interval;
		svdata dt;
		
		vector<int> G(n, 0);
		dt.set_size(n);
		
		while(V_list.read_interval_multi(interval, dt)>=0)
		{
			gtype g;
			dt.normalize(interval, avg_depths);
			string ln;
			
			if (interval.svtype == "DEL")
			{
				g.call_del(interval, dt, ln);
				if (ln != "")
				{
					vfile.print(ln);
				}
			}
			else if (interval.svtype == "CNV" || interval.svtype == "DUP")
			{
			//	g.call_cnv(interval, dt, g);
			}
			else if (interval.svtype == "INV")
			{
			//	g.call_inversion(interval, dt, g);
			}
			
			// Write Output
			
		}
 		vfile.close();
	}
	else
	{
		// Generate summary VCF from BAM/CRAM

		n = 1;
		cerr << "Processing individual BAM/CRAM file to genearte summary VCF." << endl;

		gcContent GC;
		GC.initialize(gc_file);
		
		vector<sv> candidates;
		vcfs.push_back(vcf_file);
		read_intervals_from_vcf(sample_ids, vcfs, candidates);
		
		cerr<< candidates.size() << " intervals identified from the VCF file." << endl;
		vector<int> idxs;

		// Here we assume overlapping SVs are already clustered into a single event, so do simple merging without clustering
		// - For large sample size, procesing merging by a separate task reduces sample-by-sample runtime
		// - For small sample size, probably genotyping every candidate event would yield better
		
		merge_svs(candidates, idxs);

		bFile b(GC);
		b.initialize(bam_file);
		cerr << "BAM/CRAM file initialized" << endl;
		b.get_avg_depth();
		
		FILE *fp = fopen(out_filename.c_str(), "wt");
		fprintf(fp, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n", sampID.c_str());
		fprintf(fp, "0\t0\t.\t.\t.\t.\t.\tAVGDP;AVG_ISIZE;STD_ISIZE;MED_ISIZE\tDP:AI:SI:MI\t%.1f:%.1f:%.1f:%.1f\n", b.avg_dp, b.avg_isize,b.std_isize,b.med_isize);

		for(int i=0; i<(int)idxs.size(); ++i)
		{
			int last_idx;
			if (i<idxs.size()-1)
			{
				last_idx = idxs[i+1];
			}
			else
			{
				last_idx = (int)candidates.size();
			}
			vector<sv> svlist(candidates.begin()+idxs[i], candidates.begin() + last_idx)  ;

			int m = (int)svlist.size();

			vector<string> G (m, "");
			b.read_depth(svlist, G);
			
			for(int j=0;j<m;++j)
			{
				fprintf(fp, "%s\t%d\t.\t.\t.\t.\t.\tEND=%d;SVTYPE=%s\t.", svlist[j].chr.c_str(), svlist[j].pos,svlist[j].end,svlist[j].svtype.c_str());
				fprintf(fp, "\t%s\n", G[j].c_str());
			}

		}
		fclose(fp);
		cerr << "Finished without an error" << endl;
	}
	return 0;
}
