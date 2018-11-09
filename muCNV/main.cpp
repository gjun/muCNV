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
#include <map>

#include "muCNV.h"

// TCLAP headers
#include "tclap/CmdLine.h"
#include "tclap/Arg.h"

int find_overlap_sv(sv &S , std::vector<sv>& dels)
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

bool bUseGL;
double P_THRESHOLD;
double BE_THRESHOLD;
double RO_THRESHOLD;

int main(int argc, char** argv)
{
	std::cerr << "muCNV 0.9 -- Multi-sample CNV genotyper" << std::endl;
	std::cerr << "(c) 2018 Goo Jun" << std::endl << std::endl;
	std::cerr.setf(std::ios::showpoint);

	bool bGenotype = false;
	bool bFilter = false;
	bool bNoHeader = false;
	bool bFail = false;
	bool bMerge = false;
    bool bWriteSV = false;
    bool bPrint = false;
	string index_file;
	string vcf_file;
    string interval_file;
	string supp_file;
	string supp_id_file;
	string out_filename;
	string bam_file;
	string sChr;
	string gc_file;
	string sampID;
	string region;
	
	std::vector<string> sample_ids;
	std::vector<string> vcfs;
	std::vector<string> bam_names;
//	std::map<string, unsigned> hIdSex;


	srand((unsigned int)time(NULL));

	// Parsing command-line arguments
	try 
	{
		TCLAP::CmdLine cmd("Command description message", ' ', "0.06");
		
		TCLAP::ValueArg<string> argBam("b","bam","Input BAM/CRAM file name",false,"","string");
		TCLAP::ValueArg<string> argOut("o","out","Output filename",false,"muCNV.vcf","string");
		TCLAP::ValueArg<string> argVcf("v","vcf","VCF file containing candidate SVs",false,"","string");
        TCLAP::ValueArg<string> argInterval("V","interVal", "Binary interval file containing candidate SVs", false, "", "string");
		TCLAP::ValueArg<string> argSupp("S","Support","Support VCF file containing suppporting info",false,"","string");
		TCLAP::ValueArg<string> argSuppID("I","IDinSupport","Sample ID list for support std::vectors",false,"","string");
		TCLAP::ValueArg<string> argSampleID("s","sample","Sample ID",false,"","string");
		TCLAP::ValueArg<string> argIndex("i","index","List file containing list of intermediate pileups. Required with genotype option",false,"","string");
		TCLAP::ValueArg<string> argGcfile("f","gcFile","File containing GC content information",false, "GRCh38.gc", "string");
//		TCLAP::ValueArg<double> argBE("e","error","Threshold for BayesError",false,0.1,"double");
//		TCLAP::ValueArg<double> argR("r","ratio","Threshold for likelihood ratio",false,5,"double");
		TCLAP::ValueArg<string> argRegion("r", "region", "Genomic region (chr:start-end)", false, "", "string" );
		TCLAP::SwitchArg switchFail("a", "all", "Report filter failed variants", cmd, false);
        TCLAP::SwitchArg switchPrint("p", "print", "Print out pileup info", cmd, false);
		TCLAP::SwitchArg switchMerge("m", "merge", "Merge overlapping SVs in the input VCF", cmd, false);
        TCLAP::SwitchArg switchWriteSV("w", "writevariant", "Write full list of SVs with pileup", cmd, false);
		TCLAP::SwitchArg switchFilter("t", "filter", "Filter candidate discovery set using supporting VCF", cmd, false);
		TCLAP::SwitchArg switchNoHeader("n", "noheader", "Do not print header in genoptyed VCF", cmd, false);
		TCLAP::SwitchArg switchGenotype("g","genotype","Generate Genotype from intermediate pileups", cmd, false);
		

		cmd.add(argBam);
		cmd.add(argOut);
		cmd.add(argVcf);
        cmd.add(argInterval);
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
        interval_file = argInterval.getValue();
		supp_file = argSupp.getValue();
		supp_id_file = argSuppID.getValue();
		gc_file = argGcfile.getValue();
		bGenotype = switchGenotype.getValue();
		bNoHeader = switchNoHeader.getValue();
		bFail = switchFail.getValue();
		bMerge = switchMerge.getValue();
		bFilter = switchFilter.getValue();
		region = argRegion.getValue();
        bWriteSV = switchWriteSV.getValue();
        bPrint = switchPrint.getValue();

		if (bMerge && bGenotype)
		{
			std::cerr << "Error: --genotype and --merge cannot be set together" << std::endl;
			exit(1);
		}
		
		//BE_THRESHOLD = argBE.getValue();
		//P_THRESHOLD = 1.0/ argR.getValue();
		
		BE_THRESHOLD = 0.1;
		P_THRESHOLD = 0.2;
		
		if (bGenotype && index_file == "")
		{
			std::cerr << "Error: list file is required for genotyping" << std::endl;
			exit(1);
		}
		else if (bGenotype == false && bMerge == false && bFilter == false)
		{
			if (bam_file == "")
			{
				std::cerr << "Error: BAM/CRAM file is required for individual processing mode." << std::endl;
				exit(1);
            }
            if (vcf_file == "" && interval_file == "")
            {
                std::cerr << "Error: VCF file or Interval file with SV events is required for individual processing mode." << std::endl;
                exit(1);
            }
            if (sampID == "")
            {
                std::cerr << "Error: Sample ID should be supplied with -s option." << std::endl;
				exit(1);
			}
		}
	}
	catch (TCLAP::ArgException &e)
	{
		std::cerr << "Error: " << e.error() << " for arg " << e.argId() << std::endl;
		abort();
	}
	
	int n_sample = 0;
	if (bGenotype)
	{
		// Multi-sample genotyping from summary VCFs
		
		std::vector<double> avg_depths;
		std::vector<double> avg_isizes;
		std::vector<double> std_isizes;
		
		std::map<string, int> id_to_idx;


		invcfs V_list;
		read_vcf_list(index_file, vcfs);
		V_list.initialize(vcfs, sample_ids, avg_depths, avg_isizes, std_isizes, region);

		n_sample = (int)sample_ids.size();

		for(int i=0;i<n_sample; ++i)
		{
			id_to_idx[sample_ids[i]] = i;
		}

		std::cerr << "Genotyping index loaded." << std::endl;
		std::cerr << n_sample << " samples identified." << std::endl;

		std::ifstream sfile;
		std::vector<int> id_map(n_sample, 0);

		if (supp_file != "")
		{
			if (supp_id_file != "")
			{
				std::cerr << "Supporting std::vector file and ID file provided." <<std::endl;
			}
			else
			{
				std::cerr << "Error. Supporting ID file is needed for supporting std::vector file." << std::endl;
				exit(1);
			}

			std::ifstream idfile(supp_id_file.c_str(), std::ios::in);
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
					std::cerr << "Cannot find " << ln << " in ID" << std::endl;
					exit(1);
				}
			}
			idfile.close();
		}
	
		sfile.open(supp_file.c_str(), std::ios::in);

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

		std::vector<sv> supp_list;
		std::vector<string> suppvec_list;

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
						std::cerr << "Something wrong" <<std::endl;
						std::cerr << "Curren position is " << S.pos << "-" << S.end << " while supporting std::vector has passed " << suppS.pos << std::endl;
						exit(1);
					}
				}
			}
			std::vector<double> wt(n_sample, 1);
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

				if (S.svtype == DEL)
				{
					T.call_del(S, D, G, avg_isizes, std_isizes, wt);
				}
                else if (S.svtype == DUP || S.svtype == CNV)
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
				if (S.svtype != INV  && (bFail || G.b_pass))
//				if (G.b_pass || bFail)
				{
					string ln;
					G.info += ";SUPP=" + std::to_string(suppS.supp);
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

		std::vector<sv> dels;
		std::vector<sv> dups;
		std::vector<sv> invs;
		std::vector<sv> all;

		vcfs.push_back(vcf_file);

		read_intervals_from_vcf(sample_ids, vcfs, all);
		for(int i=0;i<(int)all.size();++i)
		{
			if (all[i].svtype == DEL)
			{
				dels.push_back(all[i]);
			}
			else if (all[i].svtype == DUP)
			{
				dups.push_back(all[i]);
			}
			else if (all[i].svtype == INV)
			{
				invs.push_back(all[i]);
			}
		}
		all.clear();

		std::sort(dels.begin(), dels.end());

		std::vector<sv> merged_candidates;
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
			fprintf(fp, "%d\t%d\tSV%d\t.\t<%s>\t.\tPASS\tSUPP=%d;SVTYPE=%s;END=%d\n", all[i].chrnum, all[i].pos, i+1, svTypeName(all[i].svtype).c_str(), all[i].supp, svTypeName(all[i].svtype).c_str(), all[i].end);
		}
		fclose(fp);
	}
	else if (bFilter)
	{
		std::vector<sv> dels;
		std::vector<sv> dups;
		std::vector<sv> invs;
		std::vector<sv> all;

		vcfs.push_back(supp_file);
		read_intervals_from_vcf(sample_ids, vcfs, all);
		FILE *fp = fopen(out_filename.c_str(), "w");

		fprintf(fp, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");

		for(int i=0;i<(int)all.size();++i)
		{
			if (all[i].svtype == DEL)
			{
				dels.push_back(all[i]);
			}
			else if (all[i].svtype == DUP)
			{
				dups.push_back(all[i]);
			}
			else if (all[i].svtype == INV)
			{
				invs.push_back(all[i]);
			}
		}
		all.clear();

		std::ifstream vfile(vcf_file.c_str(), std::ios::in);
		int cnt = 1;

		while(vfile.good())
		{
			sv S;
			string suppvec;
			int idx = -1;

			if (read_candidate_vcf(vfile, S, suppvec)>0)
			{
				sv o_S;

				if (S.svtype == DEL)
				{
					idx = find_overlap_sv(S, dels);
					if (idx >=0)
					{
						o_S = dels[idx];
					}
				}
				else if (S.svtype == DUP || S.svtype==CNV)
				{
					idx = find_overlap_sv(S, dups);
					if (idx >=0)
					{
						o_S = dups[idx];
					}
				}
				else if (S.svtype == INV)
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
					fprintf(fp, "%d\t%d\tSV%d\t.\t<%s>\t.\tPASS\tSVTYPE=%s;OVERLAP=%d:%d-%d;SUPP=%d;N=%d;END=%d;SUPP_VEC=%s\n", S.chrnum, S.pos, cnt++, svTypeName(S.svtype).c_str(), svTypeName(S.svtype).c_str(), o_S.chrnum, o_S.pos, o_S.end, o_S.supp, S.supp, S.end, suppvec.c_str());
				}
				else
				{
					fprintf(fp, "%d\t%d\tSV%d\t.\t<%s>\t.\tPASS\tSVTYPE=%s;SUPP=1;N=%d;END=%d;SUPP_VEC=%s\n", S.chrnum, S.pos, cnt++, svTypeName(S.svtype).c_str(), svTypeName(S.svtype).c_str(), S.supp, S.end, suppvec.c_str());
				}
			}
		}
		fclose(fp);
		vfile.close();
	}
    else if (bPrint)
    {
        std::vector<sv> vec_sv;
        std::vector<breakpoint> vec_bp;
        string pileup_name = sampID + ".pileup";
        string varfile_name = sampID + ".var";
        string idxfile_name = sampID + ".idx";
        
        // read out and print pileup info
        read_svs_from_intfile(interval_file, vec_bp, vec_sv);
        
        std::ifstream pileupFile(pileup_name.c_str(), std::ios::in | std::ios::binary);
        std::ifstream varFile(varfile_name.c_str(), std::ios::in | std::ios::binary);
        std::ifstream idxFile(idxfile_name.c_str(), std::ios::in | std::ios::binary);
        
		size_t curr_idx = 0;

        char buf[256];
        
        pileupFile.read(reinterpret_cast<char*>(&n_sample), sizeof(int));
        printf("n_sample(pileup) : %d \n", n_sample);
        // TODO: read N-sample IDs
        pileupFile.read(reinterpret_cast<char *>(buf), 256);
        printf("sample ID(pileup) : %s\n", buf);

        double avg_dp, std_dp, avg_isize, std_isize;
        pileupFile.read(reinterpret_cast<char*>(&avg_dp), sizeof(double));
        pileupFile.read(reinterpret_cast<char*>(&std_dp), sizeof(double));
        pileupFile.read(reinterpret_cast<char*>(&avg_isize), sizeof(double));
        pileupFile.read(reinterpret_cast<char*>(&std_isize), sizeof(double));
        
        printf("AVG DP: %f, STdev: %f, AVG ISIZE: %f, STdev: %f \n", avg_dp, std_dp, avg_isize, std_isize);

        gcContent GC;
        GC.initialize(gc_file);
        
        std::vector<double> gc_factor (GC.num_bin);
        
        for(int i=0;i<GC.num_bin;++i)
        {
            pileupFile.read(reinterpret_cast<char*>(&(gc_factor[i])), sizeof(double));
            printf("GC-bin %d: %f\n", i, gc_factor[i]);
        }
		idxFile.read(reinterpret_cast<char*>(&curr_idx), sizeof(size_t));
		printf("index position %d, tellg position %d\n", (int)curr_idx, (int)pileupFile.tellg());
        
        uint64_t dpsum = 0;
        uint64_t n_dp = 0;
        
        for(int i=1; i<=GC.num_chr; ++i)
        {
            int N = ceil((double)GC.chrSize[i] / 100.0) ;
            uint16_t dp100;
            for(int j=0;j<N;++j)
            {
                pileupFile.read(reinterpret_cast<char*>(&dp100), sizeof(uint16_t));
                dpsum += dp100;
                n_dp+=1;
            }
        }
        printf("Average DP100: %d\n", (int)round((double)dpsum/n_dp/32.0));
        
        while(pileupFile.good())
        {
            uint32_t n_rp = 0;
            uint32_t n_sp = 0;

			idxFile.read(reinterpret_cast<char*>(&curr_idx), sizeof(size_t));
//			printf("index position %d, tellg position %d\n", (int)curr_idx, (int)pileupFile.tellg());

            pileupFile.read(reinterpret_cast<char*>(&n_rp), sizeof(uint32_t));
          	printf("%d readpairs\n", n_rp);
            for(int k=0; k<n_rp; ++k)
            {
                int8_t chrnum, pairstr;
				uint8_t matequal;
                int32_t selfpos, matepos;
                pileupFile.read(reinterpret_cast<char*>(&(chrnum)), sizeof(int8_t));
                pileupFile.read(reinterpret_cast<char*>(&(selfpos)), sizeof(int32_t));
                pileupFile.read(reinterpret_cast<char*>(&(matepos)), sizeof(int32_t));
                pileupFile.read(reinterpret_cast<char*>(&(matequal)), sizeof(uint8_t));
                pileupFile.read(reinterpret_cast<char*>(&(pairstr)), sizeof(int8_t));
                printf("\t%d\t%d\t%d\t%u\t%d\n", chrnum, selfpos, matepos, matequal, pairstr);
            }
            
            pileupFile.read(reinterpret_cast<char*>(&n_sp), sizeof(uint32_t));
			printf("%d split reads\n", n_sp);
            for(int k=0; k<n_sp; ++k)
            {
                int8_t chrnum;
                int32_t pos, sapos;
                int16_t firstclip, secondclip;
                pileupFile.read(reinterpret_cast<char*>(&(chrnum)), sizeof(int8_t));
                pileupFile.read(reinterpret_cast<char*>(&(pos)), sizeof(int32_t));
                pileupFile.read(reinterpret_cast<char*>(&(sapos)), sizeof(int32_t));
                pileupFile.read(reinterpret_cast<char*>(&(firstclip)), sizeof(int16_t));
                pileupFile.read(reinterpret_cast<char*>(&(secondclip)), sizeof(int16_t));
                printf("\t%d\t%d\t%d\t%d\t%d\n", chrnum, pos, sapos, firstclip,secondclip);
            }
        }
        pileupFile.close();
        
        int n_var = 0;
        varFile.read(reinterpret_cast<char*>(&n_sample), sizeof(int));
        varFile.read(reinterpret_cast<char*>(&n_var), sizeof(int));
		printf("n_var : %d\n", n_var);
        
        varFile.read(reinterpret_cast<char *>(buf), 256);
        printf("sample ID(var) : %s\n", buf);
        
        for(int i=0;i<n_var;++i)
        {
            uint16_t dp;
            vec_sv[i].print();
            varFile.read(reinterpret_cast<char*>(&dp), sizeof(uint16_t));
           	printf("var %d: %f\n",i, (dp/32.0));
        }
        varFile.close();

    }
	else
	{
		// Generate summary stats from BAM/CRAM
		n_sample = 1;
		std::cerr << "Processing individual BAM/CRAM file to genearte summary VCF." << std::endl;

		gcContent GC;
		GC.initialize(gc_file);
		std::cerr << "GC content initialized" << std::endl;

        std::vector<sv> vec_sv;
        std::vector<breakpoint> vec_bp;
		vcfs.push_back(vcf_file);
        
        if (vcf_file != "")
        {
            read_svs_from_vcf(vcf_file, vec_bp, vec_sv);
        }
        else if (interval_file != "")
        {
            read_svs_from_intfile(interval_file, vec_bp, vec_sv);
        }
        
        if (bWriteSV && vcf_file != "")
        {
            if (interval_file != "")
            {
                write_interval(interval_file, vec_sv);
            }
            else
            {
                std::cerr << "Error, interval file name is missing" << std::endl;
                exit(1);
            }
        }
        
		std::cerr<< vec_sv.size() << " svs and " << vec_bp.size() << " breakpoints identified from the VCF file." << std::endl;
		std::vector<int> idxs;

        std::sort(vec_bp.begin(), vec_bp.end());
        
		bFile b(GC);
        
		b.initialize_sequential(bam_file);
		std::cerr << "BAM/CRAM file initialized" << std::endl;
        
        b.read_depth_sequential(vec_bp, vec_sv);
        
        b.postprocess_depth(vec_sv);
        b.write_pileup(sampID, vec_sv);

		std::cerr << "Finished without an error" << std::endl;
	}
	return 0;
}
