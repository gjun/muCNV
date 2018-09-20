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

svType get_svtype(string &s)
{
    if (s == "DEL")
        return DEL;
    else if (s == "DUP")
        return DUP;
    else if (s == "INV")
        return INV;
    else if (s == "CNV")
        return CNV;
    else if (s == "INS")
        return INS;
    else
        return BND;
}
void read_vcf_list(string &index_file, vector<string> &vcf_files)
{
	ifstream inFile(index_file.c_str(), ios::in);
	
	while(inFile.good())
	{
		string ln;
		getline(inFile,ln);
		if (!ln.empty())
		{
			vcf_files.push_back(ln);
		}
	}
	cerr<< vcf_files.size() << " VCF files exist in the index file." << endl;
}
	
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
/*
int suppvcf::initialize(const char* filename, vector<string> &sample_ids, string &region)
{

	htsFile *fp = hts_open(filename, "r");
	
	if (!fp)
	{
		cerr << "Cannot open " << vcf_files[i] << endl;
		exit(1);
	}

	
	start_num.resize(n_vcf);
	num_id.resize(n_vcf);

	for(int i=0;i<n_vcf;++i)
	{
		tbx_t *tbx = tbx_index_load(vcf_files[i].c_str());
		if (!tbx)
		{
			cerr << "Cannot load tabix index " + vcf_files[i] + ".tbi" << endl;
			exit(1);
		}
		tbs.push_back(tbx);
		start_num[i] = 0;
		num_id[i] = 0;
	}

	for(int i=0;i<n_vcf;++i)
	{
		if (i>0)
		{
			start_num[i] = start_num[i-1] + num_id[i-1];
		}
		bool flag = true;
		while(flag)
		{
			kstring_t str = {0,0,0};
			
			while(hts_getline(vfs[i], KS_SEP_LINE, &str) >= 0)
			{
			
				if ( !str.l || str.s[0]!=tbs[i]->conf.meta_char )
				{
					break;
				}
				
				if (str.s[0] == '#' && str.s[1] == 'C' && str.s[2] == 'H' && str.s[3] == 'R')
				{
					// read sample ids
					vector<string> tokens;
					split(str.s, " \t\n", tokens);
					for(int j=9; j<(int)tokens.size(); ++j)
					{
						num_id[i]++;
						sample_ids.push_back(tokens[j]);
					}
				}
				else if (str.s[0] == '#' && str.s[1] == '#')
				{
					// Skip these header lines
				}
			}

			hts_itr_t *itr = tbx_itr_querys(tbs[i], "0:0");
			if (!itr)
			{
				cerr << "Cannot read averaged depth info from " << vcf_files[i] << endl;
				exit(1);
			}
			while (tbx_itr_next(vfs[i], tbs[i], itr, &str) >= 0)
			{
				// there should be only one line with chromosome 0
				// (str.s[0] == '0' && (str.s[1] == '\t' || str.s[1] == ' '))

				// Read avg depth
				flag = false;
				vector<string> tokens;
				split(str.s, " \t\n", tokens);
				for(int j=9;j<tokens.size();++j)
				{
					vector<string> fields;
					split(tokens[j].c_str(), ":", fields);

					avg_depths.push_back(atof(fields[0].c_str()));
					avg_isizes.push_back(atof(fields[1].c_str()));
					std_isizes.push_back(atof(fields[2].c_str()));
				}
			}
			tbx_itr_destroy(itr);
			
			if (region != "")
			{
				itr = tbx_itr_querys(tbs[i], region.c_str());
				if (!itr)
				{
					cerr << "Cannot parse region " << region <<  " from " << vcf_files[i] << endl;
					exit(1);
				}
				m_itr.push_back(itr);
			}
			
			free(str.s);
		}
		if (sample_ids.size() != avg_depths.size() || avg_depths.size() != avg_isizes.size())
		{
			cerr << "Error: number of fields for avegerage depth does not match" << endl;
			cerr << "sample ids has " << sample_ids.size() << " and average depths has " << avg_depths.size() << endl;
			return -1;
		}
	}
	
	cerr << "Input VCF files initialized" <<endl;
	return 0;
}

*/


int invcfs::initialize(vector<string> &vcf_files, vector<string> &sample_ids, vector<double> &avg_depths, vector<double> &avg_isizes, vector<double> &std_isizes, string &region)
{
	int n_vcf = (int)vcf_files.size();
	
	for(int i=0;i<n_vcf; ++i)
	{
		htsFile *fp = hts_open(vcf_files[i].c_str(), "r");
		
		if (!fp)
		{
			cerr << "Cannot open " << vcf_files[i] << endl;
			exit(1);
		}

		vfs.push_back(fp);
	}
	
	start_num.resize(n_vcf);
	num_id.resize(n_vcf);

	for(int i=0;i<n_vcf;++i)
	{
		tbx_t *tbx = tbx_index_load(vcf_files[i].c_str());
		if (!tbx)
		{
			cerr << "Cannot load tabix index " + vcf_files[i] + ".tbi" << endl;
			exit(1);
		}
		tbs.push_back(tbx);
		start_num[i] = 0;
		num_id[i] = 0;
	}

	for(int i=0;i<n_vcf;++i)
	{
		if (i>0)
		{
			start_num[i] = start_num[i-1] + num_id[i-1];
		}
		bool flag = true;
		while(flag)
		{
			kstring_t str = {0,0,0};
			
			while(hts_getline(vfs[i], KS_SEP_LINE, &str) >= 0)
			{
			
				if ( !str.l || str.s[0]!=tbs[i]->conf.meta_char )
				{
					break;
				}
				
				if (str.s[0] == '#' && str.s[1] == 'C' && str.s[2] == 'H' && str.s[3] == 'R')
				{
					// read sample ids
					vector<string> tokens;
					split(str.s, " \t\n", tokens);
					for(int j=9; j<(int)tokens.size(); ++j)
					{
						num_id[i]++;
						sample_ids.push_back(tokens[j]);
					}
				}
				else if (str.s[0] == '#' && str.s[1] == '#')
				{
					// Skip these header lines
				}
			}

			hts_itr_t *itr = tbx_itr_querys(tbs[i], "0:0");
			if (!itr)
			{
				cerr << "Cannot read averaged depth info from " << vcf_files[i] << endl;
				exit(1);
			}
			while (tbx_itr_next(vfs[i], tbs[i], itr, &str) >= 0)
			{
				// there should be only one line with chromosome 0
				// (str.s[0] == '0' && (str.s[1] == '\t' || str.s[1] == ' '))

				// Read avg depth
				flag = false;
				vector<string> tokens;
				split(str.s, " \t\n", tokens);
				for(int j=9;j<(int)tokens.size();++j)
				{
					vector<string> fields;
					split(tokens[j].c_str(), ":", fields);

					avg_depths.push_back(atof(fields[0].c_str()));
					avg_isizes.push_back(atof(fields[1].c_str()));
					std_isizes.push_back(atof(fields[2].c_str()));
				}
			}
			tbx_itr_destroy(itr);
			
			if (region != "")
			{
				itr = tbx_itr_querys(tbs[i], region.c_str());
				if (!itr)
				{
					cerr << "Cannot parse region " << region <<  " from " << vcf_files[i] << endl;
					exit(1);
				}
				m_itr.push_back(itr);
			}
			
			free(str.s);
		}
		if (sample_ids.size() != avg_depths.size() || avg_depths.size() != avg_isizes.size())
		{
			cerr << "Error: number of fields for avegerage depth does not match" << endl;
			cerr << "sample ids has " << sample_ids.size() << " and average depths has " << avg_depths.size() << endl;
			return -1;
		}
	}
	
	cerr << "Input VCF files initialized" <<endl;
	return 0;
}


void invcfs::parse_sv(vector<string> &tokens, sv& interval)
{
	int chr;

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
	//interval.chr = tokens[0]; // Dec 1, 2017
	interval.chrnum = chr;
	interval.pos = atoi(tokens[1].c_str());

	string info = tokens[7];

	vector<string> infotokens;

	split(info.c_str(), ";", infotokens);
	for(int j=0;j<(int)infotokens.size();++j)
	{
		vector<string> infofields;
		split(infotokens[j].c_str(), "=", infofields);
		if (infofields.size()>1)
		{
			if (infofields[0] == "END")
			{
				interval.end = atoi(infofields[1].c_str());
			}

			else if (infofields[0] == "SVTYPE")
			{
                interval.svtype = get_svtype(infofields[1]);
			}
		}
		
	} // if chr >= 1 && chr <= 22
}

double invcfs::get_second_value(string &t)
{
	if (t==".") return 0;
	
	int p=0;
	while(t[p++] != ',');
	return atof(t.substr(p,string::npos).c_str());
}

void invcfs::get_value_pair(string &t, int &n, double &x)
{
	if (t==".")
	{
		n=0;
		x=0;
		return;
	}
	
	int p=0;
	while(t[p++] != ',');
	
	n = atoi(t.substr(0, p-1).c_str());
	x = atof(t.substr(p,string::npos).c_str());
	if (x<0)
	{
		x= -1*x;
	}
}

int invcfs::read_interval_multi(sv& interval, svdata& dt, string &region)
{
	int idx = 0;
	kstring_t str = {0,0,0};
	bool bTabix = (region!="");
	int startpos = 0;
	int endpos = INT_MAX;


	if (bTabix)
	{
		vector<string> tks;
		split(region.c_str(), ":-", tks);
		if (tks.size() != 3)
		{
			cerr << "Cannot parse region " << region << endl;
			exit(1);
		}
		startpos = atoi(tks[1].c_str());
		endpos = atoi(tks[2].c_str());
	}
		
	for(int i=0;i<(int)vfs.size(); ++i)
	{
		int ret = 0;
		
		if (bTabix)
		{
			ret = tbx_itr_next(vfs[i], tbs[i], m_itr[i], &str);
		}
		else
		{
			ret = hts_getline(vfs[i], KS_SEP_LINE, &str);
		}
		if (ret>=0)
		{
			// Read per-sample depth, insert size info
			vector<string> tokens;
			split(str.s, " \t\n", tokens);
			
			// Read interval information
			if (i==0)	// Parse SV info only from the first VCF
			{
				parse_sv(tokens, interval);
				interval.get_len();
			}
			for(int j=9;j<(int)tokens.size();++j) // Parse genotype fields
			{
				vector<string> fields;
				split(tokens[j].c_str(), ":", fields);
				
				// GC corrected depth
				dt.dp[idx] = atof(fields[1].c_str());
				get_value_pair(fields[2], dt.n_cnv_pos[idx], dt.cnv_pos[idx] );
				get_value_pair(fields[3], dt.n_cnv_neg[idx], dt.cnv_neg[idx] );
				get_value_pair(fields[4], dt.n_inv_pos[idx], dt.inv_pos[idx] );
				get_value_pair(fields[5], dt.n_inv_neg[idx], dt.inv_neg[idx] );
				get_value_pair(fields[6], dt.n_isz[idx], dt.isz[idx]);
				idx++;
			}
		} // if tbx_itr_next>=0
		else
		{
			// no next record
			return -1;
		}
	}
	free(str.s);
	
	if (interval.pos < startpos || interval.pos > endpos)
	{
		return 0;
	}
	else
	{
		return 1;
	}
}

int read_candidate_vcf(ifstream &vfile, sv& new_interval, string& suppvec)
{
	string ln;
	getline(vfile, ln);

	if (!ln.empty())
	{
		if (ln[0] != '#')
		{
			int chr;
			vector<string> tokens;
			split(ln.c_str(), " \t\n", tokens);
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
			if (chr >=1 && chr <=24)  // include chrs 1-22, X, Y
			{
//				new_interval.chr = tokens[0]; // Dec 1, 2017
				new_interval.chrnum = chr;
				new_interval.pos = atoi(tokens[1].c_str());

				string info = tokens[7];
				int chr2num = new_interval.chrnum;

				vector<string> infotokens;
				split(info.c_str(), ";", infotokens);
				for(int j=0;j<(int)infotokens.size();++j)
				{
					vector<string> infofields;
					split(infotokens[j].c_str(), "=", infofields);
					if (infofields.size()>1)
					{
						if (infofields[0] == "END")
						{
							new_interval.end = atoi(infofields[1].c_str());
						}
						else if (infofields[0] == "SVTYPE")
						{
							new_interval.svtype = get_svtype(infofields[1]);

							//cerr << "\tSVTYPE: " << new_interval.svtype << endl;
						}
						else if (infofields[0] == "CHR2")
						{
                            chr2num = atoi(infofields[1].c_str()); //TODO: make chr string to int function
						}
						else if (infofields[0] == "SUPP")
						{
							new_interval.supp = atoi(infofields[1].c_str());
						}
						else if (infofields[0] == "SUPP_VEC")
						{
							suppvec = infofields[1];
						}
					}
				}
				if (chr2num == new_interval.chrnum && new_interval.pos > 0 && new_interval.end > new_interval.pos && (new_interval.end - new_interval.pos)<=10000000 )  // TEMPORARY!! 10Mb Max!
				{
					return 1;
				}
				else
				{
					return -1;
				}
			} // if chr >= 1 && chr <= 22
		}
	}
	return -1;
}

void write_interval(string &intFileName, vector<sv> &vec_sv)
{
    ofstream intFile(intFileName.c_str(), std::ios::out | std::ios::binary);
    int n_var = (int)vec_sv.size();
    
    intFile.write(reinterpret_cast<char*>(&n_var), sizeof(int));
    
    for(int i=0;i<(int)vec_sv.size();++i)
    {
        uint8_t t = (uint8_t)vec_sv[i].svtype;
        intFile.write(reinterpret_cast<char*>(&t), sizeof(uint8_t));
        intFile.write(reinterpret_cast<char*>(&i), sizeof(int));
        intFile.write(reinterpret_cast<char*>(&(vec_sv[i].chrnum)), sizeof(int));
        intFile.write(reinterpret_cast<char*>(&(vec_sv[i].pos)), sizeof(int));
        intFile.write(reinterpret_cast<char*>(&(vec_sv[i].end)), sizeof(int));
    }
    intFile.close();
}
void read_svs_from_intfile(string &intFileName, vector<breakpoint> &vec_bp, vector<sv> &vec_sv)
{
    ifstream intFile(intFileName.c_str(), ios::in | ios::binary);
    int n_var = 0;
    
    intFile.read(reinterpret_cast<char*>(&n_var), sizeof(int));
    vec_sv.resize(n_var);
    vec_bp.resize(n_var*2);
    
    for(int i=0;i<n_var;++i)
    {
        uint8_t t;
        int idx;
        intFile.read(reinterpret_cast<char*>(&t), sizeof(uint8_t));
        vec_sv[i].svtype = svTypeNum(t);
        
        intFile.read(reinterpret_cast<char*>(&idx), sizeof(int));
        // sanity check whether idx == i
        
        intFile.read(reinterpret_cast<char*>(&(vec_sv[i].chrnum)), sizeof(int));
        intFile.read(reinterpret_cast<char*>(&(vec_sv[i].pos)), sizeof(int));
        intFile.read(reinterpret_cast<char*>(&(vec_sv[i].end)), sizeof(int));
        vec_sv[i].get_len();
        vec_sv[i].dp_sum = 0;
        vec_sv[i].n_dp = 0;
        vec_sv[i].dp = 0;
        
        vec_bp[i*2].pos = vec_sv[i].pos;
        vec_bp[i*2].chrnum = vec_sv[i].chrnum;
        vec_bp[i*2].idx = i;
        vec_bp[i*2].bptype = 0;

        vec_bp[i*2+1].pos = vec_sv[i].end;
        vec_bp[i*2+1].chrnum = vec_sv[i].chrnum;
        vec_bp[i*2+1].idx = i;
        vec_bp[i*2+1].bptype = 1;
    }
    intFile.close();
}


void read_svs_from_vcf(string &vcf_file, vector<breakpoint> &v_bp, vector<sv> &v_sv)
{
    ifstream vfile(vcf_file.c_str(), ios::in);
    while(vfile.good())
    {
        string ln;
        getline(vfile, ln);
        
        if (!ln.empty() && ln[0] != '#')
        {
            int chr;
            vector<string> tokens;
            split(ln.c_str(), " \t\n", tokens);
            // Let's add error handling later (and factor this part out)
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
            if (chr >=1 && chr <=24)  // include chrs 1-22, X, Y
            {
                sv new_interval;
                new_interval.chrnum = chr;
                new_interval.pos = atoi(tokens[1].c_str());
                
                string info = tokens[7];
                
                vector<string> infotokens;
                split(info.c_str(), ";", infotokens);
                for(int j=0;j<(int)infotokens.size();++j)
                {
                    vector<string> infofields;
                    split(infotokens[j].c_str(), "=", infofields);
                    if (infofields.size()>1)
                    {
                        if (infofields[0] == "END")
                        {
                            new_interval.end = atoi(infofields[1].c_str());
                        }
                        else if (infofields[0] == "SVTYPE")
                        {
                            new_interval.svtype = get_svtype(infofields[1]);
                        }
                    }
                }
                new_interval.get_len(); // Calculate length
                
                if (new_interval.chrnum > 0 && new_interval.pos > 0 && new_interval.len > 10 && new_interval.len <= 10000000 && !in_centrome(new_interval) )  // Max SV size:  10Mb
                {
                    // 6 breakpoints:
                    // SV starting point - 500bp (or 1 if start pos < 500)
                    // SV starting point
                    // SV starting point + 500bp (or start + svlen/2 if svlen<1000)
                    // SV end point - 500bp (or end - svlen/2 if svlen < 1000)
                    // SV end point
                    // SV end point + 500bp
                    
                    v_sv.push_back(new_interval);
                    int sv_idx = (int) v_sv.size()-1;
                    breakpoint bp[2];
                    
                    bp[0].pos = new_interval.pos;
                    bp[1].pos = new_interval.end;

                    for(int k=0;k<2;++k)
                    {
                        bp[k].chrnum = new_interval.chrnum;
                        bp[k].idx = sv_idx;
                        bp[k].bptype = k;
                        v_bp.push_back(bp[k]);
                    }
                    /*
                    breakpoint bp[6];

                    bp[0].pos = (new_interval.pos > 500) ? new_interval.pos-500 : 0;
                    bp[1].pos = new_interval.pos;
                    if (new_interval.len>1001)
                    {
                        bp[2].pos = new_interval.pos + 500;
                        bp[3].pos = new_interval.end - 500;
                    }
                    else
                    {
                        bp[2].pos = new_interval.pos + new_interval.len/2-1;
                        bp[3].pos = new_interval.end - new_interval.len/2+1;
                    }
                    bp[4].pos = new_interval.end;
                    bp[5].pos = new_interval.end + 500;
                    for(int k=0;k<6;++k)
                    {
                        bp[k].chrnum = new_interval.chrnum;
                        bp[k].idx = sv_idx;
                        bp[k].bptype = k;
                        v_bp.push_back(bp[k]);
                    }
                     */
                }
            } // if chr >= 1 && chr <= 22
        }
    }
    vfile.close();
}

void read_intervals_from_vcf(vector<string> &sample_ids, vector<string> &vcf_files, vector<sv> &candidates)
{
//	for(int i=0;i<(int)sample_ids.size();++i)
	int i=0; // TEMPORARY, READ SINGLE VCF FILE
	{
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
					if (chr >=1 && chr <=24)  // include chrs 1-22, X, Y
					{
						sv new_interval;
//						new_interval.chr = tokens[0]; // Dec 1, 2017
                        new_interval.chrnum = atoi(tokens[0].c_str()); // TODO: write a function
						new_interval.pos = atoi(tokens[1].c_str());
//						new_interval.ci_pos.first = -1;
//						new_interval.ci_pos.second = -1;
//						new_interval.ci_end.first = -1;
//						new_interval.ci_end.second = -1;

//						cerr << "pos : " << new_interval.pos << endl;

						string info = tokens[7];
						int chr2num = new_interval.chrnum;
						
						//cerr << "info : " << info << endl;

						vector<string> infotokens;
						split(info.c_str(), ";", infotokens);
						for(int j=0;j<(int)infotokens.size();++j)
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
								/*
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
								*/
								else if (infofields[0] == "SVTYPE")
								{
									new_interval.svtype = get_svtype(infofields[1]);

									//cerr << "\tSVTYPE: " << new_interval.svtype << endl;
								}
								else if (infofields[0] == "CHR2")
								{
									chr2num = atoi(infofields[1].c_str()); // TODO
								}
								else if (infofields[0] == "SUPP")
								{
									new_interval.supp = atoi(infofields[1].c_str());
								}
							}
								
						}
                        if (chr2num == new_interval.chrnum && new_interval.pos > 0 && new_interval.end > new_interval.pos && (new_interval.end - new_interval.pos)<=10000000 )  // Max SV size:  10Mb
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



