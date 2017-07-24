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
#include "muCNV.h"
#include <stdlib.h>
//extern uint32_t CHR;
uint32_t CHR=1; // This is temporary, need to fix, 06/20/17

// GRCh38
int chrlen[26] = {0, 248956422,242193529,198295559,190214555,181538259,170805979,159345973,145138636,138394717,133797422,135086622,133275309,114364328,107043718,101991189,90338345,83257441,80373285,58617616,64444167,46709983,50818468,156040895,57227415, 16569};

void read_index(string index_file, vector<string> &sample_ids, vector<string> &vcf_files, vector<string> &bam_files, vector<double> &avg_depths)
{
	ifstream inFile(index_file.c_str(), ios::in);
	
	while(inFile.good())
	{
		string ln;
		getline(inFile,ln);
		if (!ln.empty())
		{
			vector<string> tokens;
			split(ln.c_str(), " \t\n", tokens);
			
			if (tokens[0].empty())
			{
				cerr << "Error loading index: empty sample ID"<< endl;
				exit(1);
			}
			else if (tokens[1].empty())
			{
				cerr << "Error loading index: cannot find vcf files" << endl;
			}
			else if (tokens[2].empty())
			{
				cerr << "Error loading index: cannot find BAM/CRAM files" << endl;
			}
			else if (tokens[3].empty())
			{
				cerr << "Error loading index: cannot find average depth info" << endl;
			}
			sample_ids.push_back(tokens[0]);
			vcf_files.push_back(tokens[1]);
			bam_files.push_back(tokens[2]);
			avg_depths.push_back(atof(tokens[3].c_str()));
		}
	}
	inFile.close();
	
}

void read_intervals_from_vcf(vector<string> &sample_ids, vector<string> &vcf_files, vector<sv> &candidates)
{
	for(int i=0;i<sample_ids.size();++i)
	{
	//	cerr << "sample ID " << sample_ids[i] << endl;
//		cerr << "opening " << vcf_files[i] << endl;
		ifstream vfile(vcf_files[i].c_str(), ios::in);
		while(vfile.good())
		{
			string ln;
			getline(vfile, ln);


			if (!ln.empty())
			{
				if (ln[0] != '#')
				{

				//	cerr << ln << endl;

					int chr;
					vector<string> tokens;
					split(ln.c_str(), " \t\n", tokens);
					// Let's add error handling later
					if (tokens[0].substr(0,3) == "chr" || tokens[0].substr(0,3) == "Chr" )
					{
						tokens[0] = tokens[0].substr(3,2);
					}
					if (tokens[0] == "X")
					{
						chr = 23;
					}
					else if (tokens[0] == "Y")
					{
						chr = 24;
					}
					else if (tokens[0] == "M" || tokens[0] == "MT")
					{
						chr = 25;
					}
					else
					{
						try
						{
							chr = atoi(tokens[0].c_str());
						}
						catch(int e)
						{
							chr = 0;
					
						}
					}
//					cerr << "chr : " << chr << endl;
					if (chr >=1 && chr <=22) // X and Y will be added later, MT will be ignored
					{
						sv new_interval;
						new_interval.chr = chr;
						new_interval.pos = atoi(tokens[1].c_str());
						new_interval.ci_pos.first = -1;
						new_interval.ci_pos.second = -1;
						new_interval.ci_end.first = -1;
						new_interval.ci_end.second = -1;

//						cerr << "pos : " << new_interval.pos << endl;

						string info = tokens[7];
						
						//cerr << "info : " << info << endl;

						vector<string> infotokens;
						split(info.c_str(), ";", infotokens);
						for(int j=0;j<infotokens.size();++j)
						{
							vector<string> infofields;
							//cerr << infotokens[j] << endl;
							split(infotokens[j].c_str(), "=", infofields);
							if (infofields.size()>1)
							{
								if (infofields[0] == "END")
								{
									new_interval.end = atoi(infofields[1].c_str());
									//cerr << "\tEND: " << new_interval.end << endl;
								}
								else if (infofields[0] == "CIPOS")
								{
									vector<string> ci;
									split(infofields[1].c_str(), ",", ci);
									new_interval.ci_pos.first = atoi(ci[0].c_str());
									new_interval.ci_pos.second = atoi(ci[1].c_str());
									//cerr << "\tCIPOS: " << new_interval.ci_pos.first << " , " << new_interval.ci_pos.second << endl;

								}
								else if (infofields[0] == "CIEND")
								{
									vector<string> ci;
									split(infofields[1].c_str(), ",", ci);
									new_interval.ci_end.first = atoi(ci[0].c_str());
									new_interval.ci_end.second = atoi(ci[1].c_str());
									//cerr << "\tCIEND: " << new_interval.ci_end.first << " , " << new_interval.ci_end.second << endl;
								}
								else if (infofields[0] == "SVTYPE")
								{
									new_interval.svtype = infofields[1];

									//cerr << "\tSVTYPE: " << new_interval.svtype << endl;
								}
							}
								
						}
						if (new_interval.pos > 0 && new_interval.end > new_interval.pos)
						{
							candidates.push_back(new_interval);
						}
					} // if chr >= 1 && chr <= 22
				}
			}
		}
		vfile.close();
		
	}
} //read_intervals_from_vcf

void bfiles::initialize(vector<string> &bnames)
{
	n = (int)bnames.size();
	
	data = (aux_t**)calloc(n, sizeof(aux_t*));
	
	for(int i=0;i<n;++i)
	{
		data[i] = (aux_t*)calloc(1, sizeof(aux_t));
		data[i]->fp = hts_open(bnames[i].c_str(), "r");

		int rf = SAM_FLAG | SAM_RNAME | SAM_POS | SAM_MAPQ | SAM_CIGAR | SAM_SEQ | SAM_QUAL;
		
		if (hts_set_opt(data[i]->fp, CRAM_OPT_REQUIRED_FIELDS, rf))
		{
			cerr << "Failed to set CRAM_OPT_REQUIRED_FIELDS value" << endl;
			exit(1);
		}
		if (hts_set_opt(data[i]->fp, CRAM_OPT_DECODE_MD, 0))
		{
			fprintf(stderr, "Failed to set CRAM_OPT_DECODE_MD value\n");
			exit(1);
		}
		data[i]->min_mapQ = 20;
		data[i]->min_len = 0;
		data[i]->hdr = sam_hdr_read(data[i]->fp);
		if (data[0]->hdr == NULL)
		{
			cerr << "Cannot open CRAM/BAM header" << endl;
			exit(1);
		}
		
		hts_idx_t* tmp_idx = sam_index_load(data[i]->fp, bnames[i].c_str());
		if (tmp_idx == NULL)
		{
			cerr << "Cannot open CRAM/BAM index" << endl;
			exit(1);
		}
		idx.push_back(tmp_idx);

	}
}

void bfiles::get_avg_depth(vector<double> &X)
{
	int n = X.size();
	vector<double> sum (n,0);
	vector<int> cnt (n,0);

	for(int c=20;c<21;++c)
	{
		int rpos=100000;
		while(rpos<chrlen[c])
		{
			char reg[100];
			sprintf(reg, "%d:%d-%d",c, rpos, rpos+10);
//			cerr << reg << endl;
			
			for(int i=0;i<n;++i)
			{

				data[i]->iter = sam_itr_querys(idx[i], data[i]->hdr, reg);
				if (data[i]->iter == NULL)
				{
					cerr << "Can't parse region" << endl;
					exit(1);
				}
			
		//		cerr << "iterator set" << endl;
			}

			bam_mplp_t mplp = bam_mplp_init(n, read_bam, (void**)data);
			int *n_plp = (int*)calloc(n, sizeof(int));
			const bam_pileup1_t **plp = (const bam_pileup1_t**)calloc(n, sizeof(bam_pileup1_t*));
			

			for(int i=0;i<n;++i)
			{
//				cerr << "reading " << i << "-th BAM file: " << endl;
		
				int tid, pos;
				
				while(bam_mplp_auto(mplp, &tid, &pos, n_plp, plp)>0)
				{
		//			cerr << "tid " << tid << " pos " << pos << endl;
					if (pos<rpos || pos >=rpos+10) continue;
					int j, m = 0;
					for(j=0; j<n_plp[i]; ++j)
					{
						const bam_pileup1_t *p = plp[i] + j;
						if (p->is_del || p->is_refskip) ++m;
						else if (bam_get_qual(p->b)[p->qpos] < 20) ++m;
					}
					sum[i] = sum[i] + n_plp[i] - m;
					cnt[i] ++;
				}
			}
			free(n_plp);
			free(plp);
			bam_mplp_destroy(mplp);
			rpos+=1000000;
		}
	}
	for(int i=0;i<n;++i)
	{
		X[i] = sum[i]/cnt[i];
	}
}


void bfiles::read_depth(sv &interval, vector<double> &X)
{
	char reg[100];
	sprintf(reg, "%d:%d-%d",interval.chr, interval.pos, interval.end);
	

	for(int i=0;i<n;++i)
	{

		data[i]->iter = sam_itr_querys(idx[i], data[i]->hdr, reg);
		if (data[i]->iter == NULL)
		{
			cerr << "Can't parse region" << endl;
			exit(1);
		}
	
//		cerr << "iterator set" << endl;
	}

	bam_mplp_t mplp = bam_mplp_init(n, read_bam, (void**)data);
	int *n_plp = (int*)calloc(n, sizeof(int));
	const bam_pileup1_t **plp = (const bam_pileup1_t**)calloc(n, sizeof(bam_pileup1_t*));
	

	for(int i=0;i<n;++i)
	{
//		cerr << "reading " << i << "-th BAM file: " << endl;
		double sum = 0;
		int cnt = 0;
	
		int tid, pos;
		
		while(bam_mplp_auto(mplp, &tid, &pos, n_plp, plp)>0)
		{
//			cerr << "tid " << tid << " pos " << pos << endl;
			if (pos<interval.pos || pos >=interval.end) continue;
			int j, m = 0;
			for(j=0; j<n_plp[i]; ++j)
			{
				const bam_pileup1_t *p = plp[i] + j;
				if (p->is_del || p->is_refskip) ++m;
				else if (bam_get_qual(p->b)[p->qpos] < 20) ++m;
			}
			sum += n_plp[i] - m;
			cnt++;
		}
		X[i] = sum/cnt;
	}
	free(n_plp);
	free(plp);
	bam_mplp_destroy(mplp);
}


void readFam(string sFamFile, map<string, unsigned> &hIdSex)
{
	ifstream inFile(sFamFile.c_str(), ios::in);

	while(inFile.good())
	{
		string ln;
		getline(inFile,ln);
		if (!ln.empty())
		{
			vector<string> tokens;
			split(ln.c_str(), " \t\n", tokens);

			if (tokens[1].empty())
			{
				cerr << "Error: empty sample ID"<< endl;
				abort();
			}
			else if (tokens[4].empty())
			{
				cerr << "Error: empty sample SEX" << endl;
			}

			if (tokens[4] == "1" || tokens[4] == "2")
			{
				hIdSex[tokens[1]] = strtoul(tokens[4].c_str(), NULL, 10);
			}
		}
	}
	inFile.close();


}



