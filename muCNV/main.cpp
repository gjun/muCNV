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

int find_overlap_sv(sv &S , vector<sv>& dels)
{
	int idx = 0;
	double max_RO = 0;
	double max_idx = 0;

	if (S.pos > 5000000)
	{
		idx = find_start(dels, S.pos - 5000000);
	}
	while(dels[idx].pos < S.end && idx<(int)dels.size())
	{
		if (dels[idx].end > S.pos)
		{
			double r = RO(dels[idx], S);
			if (r>max_RO)
			{
				max_RO = r;
				max_idx = idx;
			}
		}
		idx++;
	}
	if (max_RO>0.6)
	{
		return max_idx;
	}
	else
	{
		return -1;
	}
}

bool in_centrome(sv &S)
{
	// GRCh38 Centromere coordinates , hardcoded
	int centro[24][2] = {
		{121700000, 125100000},
		{91800000, 96000000},
		{87800000, 94000000},
		{48200000, 51800000},
		{46100000, 51400000},
		{58500000, 62600000},
		{58100000, 62100000},
		{43200000, 47200000},
		{42200000, 45500000},
		{38000000, 41600000},
		{51000000, 55800000},
		{33200000, 37800000},
		{16500000, 18900000},
		{16100000, 18200000},
		{17500000, 20500000},
		{35300000, 38400000},
		{22700000, 27400000},
		{15400000, 21500000},
		{24200000, 28100000},
		{25700000, 30400000},
		{10900000, 13000000},
		{13700000, 17400000},
		{58100000, 63800000},
		{10300000, 10600000}
	};
//	int hetero[2] = {51078349, 54425074};

//	if (S.pos >= centro[S.chrnum-1][0]-300000 && S.pos <= centro[S.chrnum-1][1]+300000)
	if (S.pos >= centro[S.chrnum-1][0] && S.pos <= centro[S.chrnum-1][1])
		return true;
	if (S.end >= centro[S.chrnum-1][0] && S.end <= centro[S.chrnum-1][1])
		return true;
//	if (S.chrnum == 7  && S.pos >= hetero[0] && S.pos <= hetero[1]+1000000)
//		return true;
//	if (S.chrnum == 7  && S.end >= hetero[0] && S.end <= hetero[1])
//		return true;
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
	bool bFilter = false;
	bool bNoHeader = false;
	bool bFail = false;
	bool bMerge = false;
	string index_file;
	string vcf_file;
	string supp_file;
	string supp_id_file;
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
		TCLAP::ValueArg<string> argSupp("S","Support","Support VCF file containing suppporting info",false,"","string");
		TCLAP::ValueArg<string> argSuppID("I","IDinSupport","Sample ID list for support vectors",false,"","string");
		TCLAP::ValueArg<string> argSampleID("s","sample","Sample ID",false,"","string");
		TCLAP::ValueArg<string> argIndex("i","index","List file containing list of intermediate pileups. Required with genotype option",false,"","string");
		TCLAP::ValueArg<string> argGcfile("f","gcFile","File containing GC content information",false, "GRCh38.gc", "string");
//		TCLAP::ValueArg<double> argBE("e","error","Threshold for BayesError",false,0.1,"double");
//		TCLAP::ValueArg<double> argR("r","ratio","Threshold for likelihood ratio",false,5,"double");
		TCLAP::ValueArg<string> argRegion("r", "region", "Genomic region (chr:start-end)", false, "", "string" );
		TCLAP::SwitchArg switchFail("a", "all", "Print filter failed variants", cmd, false);
		TCLAP::SwitchArg switchMerge("m", "merge", "Merge overlapping SVs in the input VCF", cmd, false);
		TCLAP::SwitchArg switchFilter("t", "filter", "Filter candidate discovery set using supporting VCF", cmd, false);
		TCLAP::SwitchArg switchNoHeader("n", "noheader", "Do not print header in genoptyed VCF", cmd, false);
		TCLAP::SwitchArg switchGenotype("g","genotype","Generate Genotype from intermediate pileups", cmd, false);
		

		cmd.add(argBam);
		cmd.add(argOut);
		cmd.add(argVcf);
		cmd.add(argGcfile);
		cmd.add(argIndex);
		cmd.add(argSampleID);
		cmd.add(argSuppID);
		cmd.add(argSupp);
		//cmd.add(argBE);
		//cmd.add(argR);
		cmd.add(argRegion);
		
		cmd.parse(argc, argv);
        
		bam_file = argBam.getValue();
		index_file = argIndex.getValue();
		out_filename = argOut.getValue();
		sampID = argSampleID.getValue();
		vcf_file = argVcf.getValue();
		supp_file = argSupp.getValue();
		supp_id_file = argSuppID.getValue();
		gc_file = argGcfile.getValue();
		bGenotype = switchGenotype.getValue();
		bNoHeader = switchNoHeader.getValue();
		bFail = switchFail.getValue();
		bMerge = switchMerge.getValue();
		bFilter = switchFilter.getValue();
		region = argRegion.getValue();

		if (bMerge && bGenotype)
		{
			cerr << "Error: --genotype and --merge cannot be set together" << endl;
			exit(1);
		}
		
		//BE_THRESHOLD = argBE.getValue();
		//P_THRESHOLD = 1.0/ argR.getValue();
		
		BE_THRESHOLD = 0.1;
		P_THRESHOLD = 0.2;
		
		if (bGenotype && index_file == "")
		{
			cerr << "Error: list file is required for genotyping" << endl;
			exit(1);
		}
		else if (bGenotype == false && bMerge == false && bFilter == false)
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
	
	int n_sample = 0;
	if (bGenotype)
	{
		// Multi-sample genotyping from summary VCFs
		
		vector<double> avg_depths;
		vector<double> avg_isizes;
		vector<double> std_isizes;
		
		std::map<string, int> id_to_idx;


		invcfs V_list;
		read_vcf_list(index_file, vcfs);
		V_list.initialize(vcfs, sample_ids, avg_depths, avg_isizes, std_isizes, region);

		n_sample = (int)sample_ids.size();

		for(int i=0;i<n_sample; ++i)
		{
			id_to_idx[sample_ids[i]] = i;
		}

		cerr << "Genotyping index loaded." << endl;
		cerr << n_sample << " samples identified." << endl;

		ifstream sfile;
		vector<int> id_map(n_sample, 0);

		if (supp_file != "")
		{
			if (supp_id_file != "")
			{
				cerr << "Supporting vector file and ID file provided." <<endl;
			}
			else
			{
				cerr << "Error. Supporting ID file is needed for supporting vector file." << endl;
				exit(1);
			}

			ifstream idfile(supp_id_file.c_str(), ios::in);
			for(int i=0;i<n_sample;++i)
			{
				string ln;
				getline(idfile, ln);
				if (id_to_idx.find(ln)!=id_to_idx.end())
				{
					id_map[i]  = id_to_idx[ln] ;
				}
				else
				{
					cerr << "Cannot find " << ln << " in ID" << endl;
					exit(1);
				}
			}
			idfile.close();
		}
	
		sfile.open(supp_file.c_str(), ios::in);

/*
			sv S;
			string suppvec;
			int idx = -1;

			if (read_candidate_vcf(vfile, S, suppvec)>0)
				*/

		outvcf vfile;

		vfile.open(out_filename);
		
		if (!bNoHeader)
		{
			vfile.write_header(sample_ids);
		}
		
		sv S;
		svdata D;
		svgeno G;
		
		D.set_size(n_sample);

		int val = 0;

		vector<sv> supp_list;
		vector<string> suppvec_list;

		while((val = V_list.read_interval_multi(S, D, region))>=0)
		{
			
			sv suppS;
			string suppvec;
			bool bMatch = false;
			

			if (!supp_list.empty())
			{
				if (supp_list[0].pos < S.pos)
				{
					supp_list.clear();
					suppvec_list.clear();
				}
				else
				{
					for(int i=0;i<(int)supp_list.size(); ++i)
					{
						if (supp_list[i].pos == S.pos && supp_list[i].end == S.end && supp_list[i].svtype == S.svtype)
						{
							bMatch = true;
							suppS = supp_list[i];
							suppvec = suppvec_list[i];
							break;
						}
					}
				}
			}
			if (!bMatch)
			{
				while(read_candidate_vcf(sfile, suppS, suppvec)>0)
				{
					if (suppS.pos == S.pos && !(suppS.svtype == S.svtype && suppS.end == S.end) )
					{
						supp_list.push_back(suppS);
						suppvec_list.push_back(suppvec);
					}
					else if (suppS.pos == S.pos && suppS.end == S.end && suppS.svtype == S.svtype)
					{
						break;
					}
					else if (suppS.pos > S.pos) // this should never happen
					{
						cerr << "Something wrong" <<endl;
						cerr << "Curren position is " << S.pos << "-" << S.end << " while supporting vector has passed " << suppS.pos << endl;
						exit(1);
					}
				}
			}
			vector<double> wt(n_sample, 1);
			double add_wt = 0;

			if (suppS.supp > 1 )
			{
				switch(suppS.supp)
				{
					case 2:
						add_wt = 2;
						break;
					case 3:
						add_wt = 5;
						break;
					case 4:
						add_wt = 10;
						break;
					case 5:
						add_wt = 20;
						break;
				}
				for(int i=0;i<n_sample; ++i)
				{
					if (suppvec[i] == '1')
					{
						// Give additional weights to samples found in supp_vec
						wt[id_map[i]] += add_wt;
					}
				}
			}

			if (val>0  && !in_centrome(S) )
			{
				gtype T;
				G.initialize(n_sample);
				D.normalize(S, avg_depths, avg_isizes);
				G.b_pass = false;

				if (S.svtype == "DEL")
				{
					T.call_del(S, D, G, avg_isizes, std_isizes, wt);
				}
				else if (S.svtype == "CNV" || S.svtype == "DUP")
				{
					T.call_cnv(S, D, G, avg_isizes, std_isizes, wt);
				}
				/*
				else if (S.svtype == "INV")
				{
					T.call_tmp(S, D, G, avg_isizes, std_isizes, wt);
				}
				*/
//				if (S.svtype != "INV" && (G.b_pass || (bFail && !G.b_dump)) && G.ac > 0) 
				if (S.svtype != "INV"  && (bFail || G.b_pass))
//				if (G.b_pass || bFail)
				{
					string ln;
					G.info += ";SUPP=" + to_string(suppS.supp);
					ln.reserve(n_sample*30);
					G.print(S, D, ln, wt);
					vfile.print(ln);
				}

			}
			
		}
 		vfile.close();
	}
	else if (bMerge)
	{
		RO_THRESHOLD = 0.8; 

		vector<sv> dels;
		vector<sv> dups;
		vector<sv> invs;
		vector<sv> all;

		vcfs.push_back(vcf_file);

		read_intervals_from_vcf(sample_ids, vcfs, all);
		for(int i=0;i<(int)all.size();++i)
		{
			if (all[i].svtype == "DEL")
			{
				dels.push_back(all[i]);
			}
			else if (all[i].svtype == "DUP")
			{
				dups.push_back(all[i]);
			}
			else if (all[i].svtype == "INV")
			{
				invs.push_back(all[i]);
			}
		}
		all.clear();

		sort(dels.begin(), dels.end());

		vector<sv> merged_candidates;
		merged_candidates.push_back(dels[0]);
		for(int i=1; i<(int)dels.size(); ++i)
		{
			if (RO(dels[i], merged_candidates.back() ) > RO_THRESHOLD)
			{
				merged_candidates.push_back(dels[i]);
			}
			else
			{
				sv best_sv;
				pick_sv_from_merged(best_sv, merged_candidates);
				all.push_back(best_sv);

				merged_candidates.clear();
				merged_candidates.push_back(dels[i]);

			}
		}
		if (merged_candidates.size() > 0)
		{
			sv best_sv;
			pick_sv_from_merged(best_sv, merged_candidates);
			all.push_back(best_sv);
		}
		merged_candidates.clear();

		merged_candidates.push_back(dups[0]);
		for(int i=1; i<(int)dups.size(); ++i)
		{
			if (RO(dups[i], merged_candidates.back() ) > RO_THRESHOLD)
			{
				merged_candidates.push_back(dups[i]);
			}
			else
			{
				sv best_sv;
				pick_sv_from_merged(best_sv, merged_candidates);
				all.push_back(best_sv);

				merged_candidates.clear();
				merged_candidates.push_back(dups[i]);

			}
		}
		if (merged_candidates.size() > 0)
		{
			sv best_sv;
			pick_sv_from_merged(best_sv, merged_candidates);
			all.push_back(best_sv);
		}
		merged_candidates.clear();

		merged_candidates.push_back(invs[0]);
		for(int i=1; i<(int)invs.size(); ++i)
		{
			if (RO(invs[i], merged_candidates.back() ) > RO_THRESHOLD)
			{
				merged_candidates.push_back(invs[i]);
			}
			else
			{
				sv best_sv;
				pick_sv_from_merged(best_sv, merged_candidates);
				all.push_back(best_sv);

				merged_candidates.clear();
				merged_candidates.push_back(invs[i]);

			}
		}
		if (merged_candidates.size() > 0)
		{
			sv best_sv;
			pick_sv_from_merged(best_sv, merged_candidates);
			all.push_back(best_sv);
		}
		merged_candidates.clear();

		std::sort(all.begin(), all.end());

		// Write to output file
		FILE *fp = fopen(out_filename.c_str(), "w");
		fprintf(fp, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
		for(int i=0; i<(int)all.size(); ++i)
		{
			fprintf(fp, "%d\t%d\tSV%d\t.\t<%s>\t.\tPASS\tSUPP=%d;SVTYPE=%s;END=%d\n", all[i].chrnum, all[i].pos, i+1, all[i].svtype.c_str(), all[i].supp, all[i].svtype.c_str(), all[i].end);
		}
		fclose(fp);
	}
	else if (bFilter)
	{
		vector<sv> dels;
		vector<sv> dups;
		vector<sv> invs;
		vector<sv> all;

		vcfs.push_back(supp_file);
		read_intervals_from_vcf(sample_ids, vcfs, all);
		FILE *fp = fopen(out_filename.c_str(), "w");

		fprintf(fp, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");

		for(int i=0;i<(int)all.size();++i)
		{
			if (all[i].svtype == "DEL")
			{
				dels.push_back(all[i]);
			}
			else if (all[i].svtype == "DUP")
			{
				dups.push_back(all[i]);
			}
			else if (all[i].svtype == "INV")
			{
				invs.push_back(all[i]);
			}
		}
		all.clear();

		ifstream vfile(vcf_file.c_str(), ios::in);
		int cnt = 1;

		while(vfile.good())
		{
			sv S;
			string suppvec;
			int idx = -1;

			if (read_candidate_vcf(vfile, S, suppvec)>0)
			{
				sv o_S;

				if (S.svtype == "DEL")
				{
					idx = find_overlap_sv(S, dels);
					if (idx >=0)
					{
						o_S = dels[idx];
					}
				}
				else if (S.svtype == "DUP" || S.svtype=="CNV")
				{
					idx = find_overlap_sv(S, dups);
					if (idx >=0)
					{
						o_S = dups[idx];
					}
				}
				else if (S.svtype == "INV")
				{
					idx = find_overlap_sv(S, invs);
					if (idx >=0)
					{
						o_S = invs[idx];
					}
				}
				//print record
				if (idx>=0)
				{
					fprintf(fp, "%s\t%d\tSV%d\t.\t<%s>\t.\tPASS\tSVTYPE=%s;OVERLAP=%s:%d-%d;SUPP=%d;N=%d;END=%d;SUPP_VEC=%s\n", S.chr.c_str(), S.pos, cnt++, S.svtype.c_str(), S.svtype.c_str(), o_S.chr.c_str(), o_S.pos, o_S.end, o_S.supp, S.supp, S.end, suppvec.c_str());
				}
				else
				{
					fprintf(fp, "%s\t%d\tSV%d\t.\t<%s>\t.\tPASS\tSVTYPE=%s;SUPP=1;N=%d;END=%d;SUPP_VEC=%s\n", S.chr.c_str(), S.pos, cnt++, S.svtype.c_str(), S.svtype.c_str(), S.supp, S.end, suppvec.c_str());
				}
			}
		}
		fclose(fp);
		vfile.close();
	}
	else
	{
		// Generate summary VCF from BAM/CRAM
		n_sample = 1;
		cerr << "Processing individual BAM/CRAM file to genearte summary VCF." << endl;

		gcContent GC;
		GC.initialize(gc_file);
		
		vector<sv> candidates;
		vcfs.push_back(vcf_file);
		read_intervals_from_vcf(sample_ids, vcfs, candidates);
		
		cerr<< candidates.size() << " intervals identified from the VCF file." << endl;
		vector<int> idxs;

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
			if (i<(int)idxs.size()-1)
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
