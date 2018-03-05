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

bool in_centrome(sv &S)
{
	int centro[24][2] = {
		{122026460,125184587},
		{92188146,94090557},
		{90772459,93655574},
		{49708101,51743951},
		{46485901,50059807},
		{58553889,59829934},
		{58169654,60828234},
		{44033745,45877265},
		{43236168,45518558},
		{39686683,41593521},
		{51078349,54425074},
		{34769408,37185252},
		{16000001,18051248},
		{16000001,18173523},
		{17000001,19725254},
		{36311159,38280682},
		{22813680,26885980},
		{15460900,20861206},
		{24498981,27190874},
		{26436233,30038348},
		{10864561,12915808},
		{12954789,15054318},
		{58605580,62412542},
		{10316945,10544039}
	};
	int hetero[2] = {51078349, 54425074};

	if (S.pos >= centro[S.chrnum-1][0]-300000 && S.pos <= centro[S.chrnum-1][1]+300000)
//	if (S.pos >= centro[S.chrnum-1][0] && S.pos <= centro[S.chrnum-1][1])
		return true;
	if (S.end >= centro[S.chrnum-1][0]-300000 && S.end <= centro[S.chrnum-1][1]+300000)
		return true;
	if (S.chrnum == 7  && S.pos >= hetero[0] && S.pos <= hetero[1]+1000000)
		return true;
	if (S.chrnum == 7  && S.end >= hetero[0] && S.end <= hetero[1])
		return true;
	return false;
}


bool bUseGL;
double P_THRESHOLD;
double BE_THRESHOLD;
double RO_THRESHOLD;

int main(int argc, char** argv)
{
	cerr << "muCNV 0.9 -- Multi-sample CNV genotyper" << endl;
	cerr << "(c) 2018 Goo Jun" << endl << endl;
	cerr.setf(ios::showpoint);

	bool bGenotype = false;
	bool bNoHeader = false;
	bool bFail = false;
	string index_file;
	string vcf_file;
	string out_filename;
	string bam_file;
	string sChr;
	string gc_file;
	string sampID;
	string region;
	
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
//		TCLAP::ValueArg<double> argBE("e","error","Threshold for BayesError",false,0.1,"double");
//		TCLAP::ValueArg<double> argR("r","ratio","Threshold for likelihood ratio",false,5,"double");
		TCLAP::ValueArg<string> argRegion("r", "region", "Genomic region (chr:start-end)", false, "", "string" );
		TCLAP::SwitchArg switchFail("a", "all", "Print filter failed variants", cmd, false);
		TCLAP::SwitchArg switchNoHeader("n", "noheader", "Do not print header in genoptyed VCF", cmd, false);
		TCLAP::SwitchArg switchGenotype("g","genotype","Generate Genotype from intermediate pileups", cmd, false);
		

		cmd.add(argBam);
		cmd.add(argOut);
		cmd.add(argVcf);
		cmd.add(argGcfile);
		cmd.add(argIndex);
		cmd.add(argSampleID);
		//cmd.add(argBE);
		//cmd.add(argR);
		cmd.add(argRegion);
		
		cmd.parse(argc, argv);
        
		bam_file = argBam.getValue();
		index_file = argIndex.getValue();
		out_filename = argOut.getValue();
		sampID = argSampleID.getValue();
		vcf_file = argVcf.getValue();
		gc_file = argGcfile.getValue();
		bGenotype = switchGenotype.getValue();
		bNoHeader = switchNoHeader.getValue();
		bFail = switchFail.getValue();
		region = argRegion.getValue();
		
		//BE_THRESHOLD = argBE.getValue();
		//P_THRESHOLD = 1.0/ argR.getValue();
		
		BE_THRESHOLD = 0.1;
		P_THRESHOLD = 0.25;
		
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
		vector<double> std_isizes;
		
		invcfs V_list;
		read_vcf_list(index_file, vcfs);
		V_list.initialize(vcfs, sample_ids, avg_depths, avg_isizes, std_isizes, region);

		n = (int)sample_ids.size();
		cerr << "Genotyping index loaded." << endl;
		cerr << n << " samples identified." << endl;
	
		outvcf vfile;

		vfile.open(out_filename);
		
		if (!bNoHeader)
		{
			vfile.write_header(sample_ids);
		}
		
		sv S;
		svdata D;
		svgeno G;
		
		D.set_size(n);

		int val;

		while((val = V_list.read_interval_multi(S, D, region))>=0)
		{
			if (val>0 && !in_centrome(S) )
			{
				gtype T;
				G.initialize(n);
				D.normalize(S, avg_depths, avg_isizes);
				if (S.svtype == "DEL")
				{
					T.call_del(S, D, G, avg_isizes, std_isizes);
				}
				else if (S.svtype == "CNV" || S.svtype == "DUP")
				{
					T.call_cnv(S, D, G, avg_isizes, std_isizes);
				}
				else if (S.svtype == "INV")
				{
					//T.call_inv(S, D, G, avg_isizes, std_isizes);
				}
				if (S.svtype != "INV" && (G.b_pass || (bFail && !G.b_dump)) && G.ac > 0) 
				{
					string ln;
					ln.reserve(n*30);
					G.print(S, D, ln);
					vfile.print(ln);
				}

			}
			
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
