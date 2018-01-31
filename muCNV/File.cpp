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


void invcfs::initialize(vector<string> &vcf_files)
{
	for(int i=0;i<(int)vcf_files.size(); ++i)
	{
		//	cerr << "sample ID " << sample_ids[i] << endl;
		//		cerr << "opening " << vcf_files[i] << endl;
		ifstream *f = new ifstream(vcf_files[i].c_str(), ios::in);
		vfs.push_back(f);
	}

	cerr << "Input VCF files initialized" <<endl;
}

int invcfs::read_interval(sv& interval, vector<double> &X)
{
	vector<string> lns (vfs.size(), "");

//	cerr << "Reading vcf " <<endl;

	for(int i=0;i<(int)vfs.size();++i)
	{
		if (!vfs[i]->good())
			return -1;
	}
	bool flag=false;
	for(int i=0;i<(int)vfs.size();++i)
	{
		getline(*vfs[i],lns[i]);

		if (lns[i].empty() || lns[i][0] == '#')
		{
//			cerr << i << "-th vcf has " << lns[i];
			flag = true;
		}
	}
	if (flag)
		return 1;

//    cerr << "Reading depth" <<endl;

	int chr;
	vector<string> tokens;
	split(lns[0].c_str(), " \t\n", tokens);
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
	interval.chr = tokens[0]; // Dec 1, 2017
	interval.chrnum = chr;
	interval.pos = atoi(tokens[1].c_str());
//	interval.ci_pos.first = -1;
//	interval.ci_pos.second = -1;
//	interval.ci_end.first = -1;
//	interval.ci_end.second = -1;
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
			/*
			else if (infofields[0] == "CIPOS")
			{
				vector<string> ci;
				split(infofields[1].c_str(), ",", ci);
				interval.ci_pos.first = atoi(ci[0].c_str());
				interval.ci_pos.second = atoi(ci[1].c_str());
			}
			else if (infofields[0] == "CIEND")
			{
				vector<string> ci;
				split(infofields[1].c_str(), ",", ci);
				interval.ci_end.first = atoi(ci[0].c_str());
				interval.ci_end.second = atoi(ci[1].c_str());
			}
			*/
			else if (infofields[0] == "SVTYPE")
			{
				interval.svtype = infofields[1];
			}
		}
		
	} // if chr >= 1 && chr <= 22

//   cerr << "SV parsed" << interval.chr << ":" << interval.pos << "-" << interval.end << " depth " << tokens[9] << endl;

	X[0] = atof(tokens[9].c_str());
	
   // cerr << "First sample depth read" <<endl;

	for(int i=1;i<(int)vfs.size();++i)
	{
		vector<string> tks;
		split(lns[i].c_str(), " \t\n", tks);

		X[i] = atof(tks[9].c_str());

    //	cerr << i << "-th sample depth read" <<endl;
	}
	return 0;
} //invcfs::read_vcf

void read_intervals_from_vcf(vector<string> &sample_ids, vector<string> &vcf_files, vector<sv> &candidates)
{
//	for(int i=0;i<(int)sample_ids.size();++i)
int i=0; // TEMPORARY, READ SINGLE VCF FILE
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
					if (chr >=1 && chr <=24)  // include chrs 1-22, X, Y
					{
						sv new_interval;
						new_interval.chr = tokens[0]; // Dec 1, 2017
						new_interval.chrnum = chr;
						new_interval.pos = atoi(tokens[1].c_str());
//						new_interval.ci_pos.first = -1;
//						new_interval.ci_pos.second = -1;
//						new_interval.ci_end.first = -1;
//						new_interval.ci_end.second = -1;

//						cerr << "pos : " << new_interval.pos << endl;

						string info = tokens[7];
						string chr2 = new_interval.chr;
						
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
									new_interval.svtype = infofields[1];

									//cerr << "\tSVTYPE: " << new_interval.svtype << endl;
								}
								else if (infofields[0] == "CHR2")
								{
									chr2 = infofields[1];
								}
							}
								
						}
						if (chr2 == new_interval.chr && new_interval.pos > 0 && new_interval.end > new_interval.pos && (new_interval.end - new_interval.pos)<=10000000 )  // TEMPORARY!! 10Mb Max!
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



