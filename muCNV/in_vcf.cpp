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
#include <map>
#include <stdlib.h>
#include "in_vcf.h"
#include "common.h"
#include <iostream>
#include <fstream>

svType get_svtype(std::string &s)
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

void read_list(std::string &index_file, std::vector<std::string> &items)
{
	std::ifstream inFile(index_file.c_str(), std::ios::in);
	
	while(inFile.good())
	{
		std::string ln;
		getline(inFile,ln);
		if (!ln.empty())
		{
			items.push_back(ln);
		}
	}
	std::cerr << items.size() << " items exist in the index file " << index_file << std::endl;
}
	
void read_index(std::string index_file, std::vector<std::string> &sample_ids, std::vector<std::string> &vcf_files, std::vector<std::string> &bam_files, std::vector<double> &avg_depths)
{
	std::ifstream inFile(index_file.c_str(), std::ios::in);
	
	while(inFile.good())
	{
		std::string ln;
		getline(inFile,ln);
		if (!ln.empty())
		{
			std::vector<std::string> tokens;
			split(ln.c_str(), " \t\n", tokens);
			
			if (tokens[0].empty())
			{
				std::cerr << "Error loading index: empty sample ID"<< std::endl;
				exit(1);
			}
			else if (tokens[1].empty())
			{
				std::cerr << "Error loading index: cannot find vcf files" << std::endl;
			}
			else if (tokens[2].empty())
			{
				std::cerr << "Error loading index: cannot find BAM/CRAM files" << std::endl;
			}
			else if (tokens[3].empty())
			{
				std::cerr << "Error loading index: cannot find average depth info" << std::endl;
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
class invvcf
{
    public:
        htsFile* vcf;
        hts_itr_t* itr;
        bcf_hdr_t *hdr;
        int open(std::string, std::vector<string> &); 
        int read_next(sv&, std::vector<int> &); // read next line
};
*/

int invcf::open(std::string filename)
{
	// Open VCF file
	vcf = bcf_open(filename.c_str(), "r");
	if (vcf == NULL)
	{
		return -1;
	}
	hdr = bcf_hdr_read(vcf);
	n_sample =  bcf_hdr_nsamples(hdr);
	line = bcf_init();

	return n_sample;
}

int invcf::read_next(sv &curr_sv, std::vector<int> &genos)
{
	if (bcf_read(vcf, hdr, line) == 0)
	{
		curr_sv.pos = line->pos + 1;
		curr_sv.chrnum = line->rid + 1;
		int ngt = 0;
		//int ngt_arr = 0;

		int nend = 0;
		int ncallrate = 0;
		int32_t* end = NULL;
		float* callrate = NULL;
		bcf_get_info_int32(hdr, line, "END", &end, &nend);
		curr_sv.end = *end;

		bcf_get_info_float(hdr, line, "CALLRATE", &callrate, &ncallrate);
		curr_sv.supp = (int) ((*callrate)*100); 

		char* svt = NULL;
		int nsvt = 0;
		bcf_get_info_string(hdr, line, "SVTYPE", &svt, &nsvt);
		std::string svt_str(svt);
		curr_sv.svtype = get_svtype(svt_str);
		free(end);
		free(callrate);
		free(svt);
		// std::cerr << "First record: chrom id " << curr_sv.chrnum << " pos " << curr_sv.pos << " end " << curr_sv.end << " type " << svt_str << " callrate " << curr_sv.supp << std::endl;

		int *gt = NULL;
		int ngt_arr = 0;
		ngt = bcf_get_genotypes(hdr, line, &gt, &ngt_arr);
 		// std::cerr << ngt << " genotypes read " << "ngt_arr has " << ngt_arr << std::endl;
		for(int i=0; i<ngt; i+=2)
		{
			if (gt[i]==0 && gt[i+1]==0)
			{
				genos[i/2] = -1;
			}
			else if(gt[i]==2 && gt[i+1]==2)
			{
				genos[i/2] = 0;
			}
			else if (gt[i]==2 && gt[i+1]==4)
			{
				genos[i/2] = 1;
			}
			else if (gt[i]==4 && gt[i+1]==4)
			{
				genos[i/2] = 2;
			}
			else
			{
				genos[i/2] = -1;
			}
		}
		free(gt);
		return 0;
	}
	else
	{
		//std::cerr << "failed" <<std::endl;
		return -1;
	}
}

/*
int suppvcf::initialize(const char* filename, std::vector<std::string> &sample_ids, std::string &region)
{

	htsFile *fp = hts_open(filename, "r");
	
	if (!fp)
	{
		std::cerr << "Cannot open " << vcf_files[i] << std::endl;
		exit(1);
	}

	
	start_num.resize(n_vcf);
	num_id.resize(n_vcf);

	for(int i=0;i<n_vcf;++i)
	{
		tbx_t *tbx = tbx_index_load(vcf_files[i].c_str());
		if (!tbx)
		{
			std::cerr << "Cannot load tabix index " + vcf_files[i] + ".tbi" << std::endl;
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
					std::vector<std::string> tokens;
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
				std::cerr << "Cannot read averaged depth info from " << vcf_files[i] << std::endl;
				exit(1);
			}
			while (tbx_itr_next(vfs[i], tbs[i], itr, &str) >= 0)
			{
				// there should be only one line with chromosome 0
				// (str.s[0] == '0' && (str.s[1] == '\t' || str.s[1] == ' '))

				// Read avg depth
				flag = false;
				std::vector<std::string> tokens;
				split(str.s, " \t\n", tokens);
				for(int j=9;j<tokens.size();++j)
				{
					std::vector<std::string> fields;
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
					std::cerr << "Cannot parse region " << region <<  " from " << vcf_files[i] << std::endl;
					exit(1);
				}
				m_itr.push_back(itr);
			}
			
			free(str.s);
		}
		if (sample_ids.size() != avg_depths.size() || avg_depths.size() != avg_isizes.size())
		{
			std::cerr << "Error: number of fields for avegerage depth does not match" << std::endl;
			std::cerr << "sample ids has " << sample_ids.size() << " and average depths has " << avg_depths.size() << std::endl;
			return -1;
		}
	}
	
	std::cerr << "Input VCF files initialized" <<std::endl;
	return 0;
}

*/

/*
int invcfs::initialize(std::vector<std::string> &vcf_files, std::vector<std::string> &sample_ids, std::vector<double> &avg_depths, std::vector<double> &avg_isizes, std::vector<double> &std_isizes, std::string &region)
{
	int n_vcf = (int)vcf_files.size();
	
	for(int i=0;i<n_vcf; ++i)
	{
		htsFile *fp = hts_open(vcf_files[i].c_str(), "r");
		
		if (!fp)
		{
			std::cerr << "Cannot open " << vcf_files[i] << std::endl;
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
			std::cerr << "Cannot load tabix index " + vcf_files[i] + ".tbi" << std::endl;
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
					std::vector<std::string> tokens;
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
				std::cerr << "Cannot read averaged depth info from " << vcf_files[i] << std::endl;
				exit(1);
			}
			while (tbx_itr_next(vfs[i], tbs[i], itr, &str) >= 0)
			{
				// there should be only one line with chromosome 0
				// (str.s[0] == '0' && (str.s[1] == '\t' || str.s[1] == ' '))

				// Read avg depth
				flag = false;
				std::vector<std::string> tokens;
				split(str.s, " \t\n", tokens);
				for(int j=9;j<(int)tokens.size();++j)
				{
					std::vector<std::string> fields;
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
					std::cerr << "Cannot parse region " << region <<  " from " << vcf_files[i] << std::endl;
					exit(1);
				}
				m_itr.push_back(itr);
			}
			
			free(str.s);
		}
		if (sample_ids.size() != avg_depths.size() || avg_depths.size() != avg_isizes.size())
		{
			std::cerr << "Error: number of fields for avegerage depth does not match" << std::endl;
			std::cerr << "sample ids has " << sample_ids.size() << " and average depths has " << avg_depths.size() << std::endl;
			return -1;
		}
	}
	
	std::cerr << "Input VCF files initialized" <<std::endl;
	return 0;
}


void invcfs::parse_sv(std::vector<std::string> &tokens, sv& interval)
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

	std::string info = tokens[7];

	std::vector<std::string> infotokens;

	split(info.c_str(), ";", infotokens);
	for(int j=0;j<(int)infotokens.size();++j)
	{
		std::vector<std::string> infofields;
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

double invcfs::get_second_value(std::string &t)
{
	if (t==".") return 0;
	
	int p=0;
	while(t[p++] != ',');
	return atof(t.substr(p,std::string::npos).c_str());
}

void invcfs::get_value_pair(std::string &t, int &n, double &x)
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
	x = atof(t.substr(p,std::string::npos).c_str());
	if (x<0)
	{
		x= -1*x;
	}
}

*/
int read_candidate_vcf(std::ifstream &vfile, sv& new_interval, std::string& suppvec)
{
	std::string ln;
	getline(vfile, ln);

	if (!ln.empty())
	{
		if (ln[0] != '#')
		{
			int chr;
			std::vector<std::string> tokens;
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

				std::string info = tokens[7];
				int chr2num = new_interval.chrnum;

				std::vector<std::string> infotokens;
				split(info.c_str(), ";", infotokens);
				for(int j=0;j<(int)infotokens.size();++j)
				{
					std::vector<std::string> infofields;
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

							//std::cerr << "\tSVTYPE: " << new_interval.svtype << std::endl;
						}
						else if (infofields[0] == "CHR2")
						{
                            chr2num = atoi(infofields[1].c_str()); //TODO: make chr std::string to int function
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

void write_interval(std::string &intFileName, std::vector<sv> &vec_sv)
{
    std::ofstream intFile(intFileName.c_str(), std::ios::out | std::ios::binary);
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

void write_svs_into_vcf(std::string &vcfFileName, std::vector<sv> &vec_sv)
{
    FILE *fp = fopen(vcfFileName.c_str(), "w");

    int n_var = (int)vec_sv.size();
	fprintf(fp, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
	
    for(int i=0;i<(int)vec_sv.size();++i)
    {
		sv& S = vec_sv[i];

		const char *svtype = svTypeName(S.svtype).c_str();
		
		// TODO: convert hard-coded chr-names using real chr names from BAM/CRAM header - this should be stored in GC content file
		// current code works only with GRCh38
		if (S.chrnum < 23)
		{
			fprintf(fp, "chr%d",  S.chrnum);
		}
		else if (S.chrnum == 23)
		{
			fprintf(fp, "chrX");
		}
		else if (S.chrnum == 24)
		{
			fprintf(fp, "chrY");
		}
		else if (S.chrnum == 25)
		{
			fprintf(fp, "chrM");
		}

		fprintf(fp, "\t%d\t%s_%d:%d-%d\t.\t<%s>\t.\t", S.pos, svtype, S.chrnum, S.pos, S.end, svtype);
		fprintf(fp, "PASS\t");

		fprintf(fp, "SVTYPE=%s;END=%d;SVLEN=%d\n",  svtype, S.end, S.len);
    }
    fclose(fp);
}

void read_svs_from_intfile(std::string &intFileName, std::vector<breakpoint> &vec_bp, std::vector<sv> &vec_sv)
{
    std::ifstream intFile(intFileName.c_str(), std::ios::in | std::ios::binary);
    int n_var = 0;
    
    if (!intFile.good())
    {
        std::cerr << "Error, cannot open " << intFileName << std::endl;
        exit(1);
    }
    
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


void read_svs_from_vcf(std::string &vcf_file, std::vector<breakpoint> &v_bp, std::vector<sv> &v_sv)
{
    std::ifstream vfile(vcf_file.c_str(), std::ios::in);
    while(vfile.good())
    {
        std::string ln;
        getline(vfile, ln);
        
        if (!ln.empty() && ln[0] != '#')
        {
            int chr;
            std::vector<std::string> tokens;
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
                
                std::string info = tokens[7];
                
                std::vector<std::string> infotokens;
                split(info.c_str(), ";", infotokens);
                for(int j=0;j<(int)infotokens.size();++j)
                {
                    std::vector<std::string> infofields;
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

void read_intervals_from_vcf(std::vector<std::string> &sample_ids, std::vector<std::string> &vcf_files, std::vector<sv> &candidates)
{
//	for(int i=0;i<(int)sample_ids.size();++i)
	int i=0; // TEMPORARY, READ SINGLE VCF FILE
	{
		std::ifstream vfile(vcf_files[i].c_str(), std::ios::in);
		while(vfile.good())
		{
			std::string ln;
			getline(vfile, ln);


			if (!ln.empty())
			{
				if (ln[0] != '#')
				{

				//	std::cerr << ln << std::endl;

					int chr;
					std::vector<std::string> tokens;
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
//					std::cerr << "chr : " << chr << std::endl;
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

//						std::cerr << "pos : " << new_interval.pos << std::endl;

						std::string info = tokens[7];
						int chr2num = new_interval.chrnum;
						
						// std::cerr << "info : " << info << std::endl;

						std::vector<std::string> infotokens;
						split(info.c_str(), ";", infotokens);
						for(int j=0;j<(int)infotokens.size();++j)
						{
							std::vector<std::string> infofields;
							//std::cerr << infotokens[j] << std::endl;
							split(infotokens[j].c_str(), "=", infofields);
							if (infofields.size()>1)
							{
								if (infofields[0] == "END")
								{
									new_interval.end = atoi(infofields[1].c_str());
									//std::cerr << "\tEND: " << new_interval.end << std::endl;
								}
								/*
								else if (infofields[0] == "CIPOS")
								{
									std::vector<std::string> ci;
									split(infofields[1].c_str(), ",", ci);
									new_interval.ci_pos.first = atoi(ci[0].c_str());
									new_interval.ci_pos.second = atoi(ci[1].c_str());
									// std::cerr << "\tCIPOS: " << new_interval.ci_pos.first << " , " << new_interval.ci_pos.second << std::endl;

								}
								else if (infofields[0] == "CIEND")
								{
									std::vector<std::string> ci;
									split(infofields[1].c_str(), ",", ci);
									new_interval.ci_end.first = atoi(ci[0].c_str());
									new_interval.ci_end.second = atoi(ci[1].c_str());
									//std::cerr << "\tCIEND: " << new_interval.ci_end.first << " , " << new_interval.ci_end.second << std::endl;
								}
								*/
								else if (infofields[0] == "SVTYPE")
								{
									new_interval.svtype = get_svtype(infofields[1]);

									//std::cerr << "\tSVTYPE: " << new_interval.svtype << std::endl;
								}
								else if (infofields[0] == "CHR2")
								{
									chr2num = atoi(infofields[1].c_str()); // TODO
								}
								else if (infofields[0] == "CALLRATE")
								{
									new_interval.supp = (int)atof(infofields[1].c_str())*100;
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

void readFam(std::string sFamFile, std::map<std::string, unsigned> &hIdSex)
{
	std::ifstream inFile(sFamFile.c_str(), std::ios::in);

	while(inFile.good())
	{
		std::string ln;
		getline(inFile,ln);
		if (!ln.empty())
		{
			std::vector<std::string> tokens;
			split(ln.c_str(), " \t\n", tokens);

			if (tokens[1].empty())
			{
				std::cerr << "Error: empty sample ID"<< std::endl;
				abort();
			}
			else if (tokens[4].empty())
			{
				std::cerr << "Error: empty sample SEX" << std::endl;
			}

			if (tokens[4] == "1" || tokens[4] == "2")
			{
				hIdSex[tokens[1]] = (unsigned) strtoul(tokens[4].c_str(), NULL, 10);
			}
		}
	}
	inFile.close();


}



