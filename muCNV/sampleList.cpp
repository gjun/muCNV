/*
 *  Copyright (C) 2010-2014  Regents of the University of Michigan
 *
 *  Author: Goo Jun
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

#include <iostream>
#include "muCNV.h"
#include "pFile.h"


bool SampleList::readIndex(string sInFile)
{
	ifstream inFile(sInFile.c_str(), ios::in);

	while(inFile.good())
	{
		string ln;
		getline(inFile,ln);
		if (!ln.empty())
		{
			vector<string> tokens;
			pFile::tokenizeLine(ln.c_str(), " \t\n", tokens);

			if (tokens[0].empty())
			{
				cerr << "Error: empty sample ID"<< endl;
				abort();
			}
			else if (tokens[1].empty())
			{
				cerr << "Error: depth file is empty" << endl;
			}

			vID.push_back(tokens[0]);
			vDepthFile.push_back(tokens[1]);


			// To Do: do this earlier and get vector of *pFile instead of file names
			pFile dFile;
			dFile.load(tokens[1].c_str(), NULL, true); // read the file with header
			const char *line;

			line = dFile.getLine();

			tokens.clear();
			pFile::tokenizeLine(line,"#= \t\n",tokens);
			// TMPTMP
			if (tokens[1].compare("gcCorrectedCvg") == 0)
			{
				AvgDepth.push_back(atof(tokens[2].c_str()));
			}
			else
			{
				//header is in wrong format
			}

			inFile.close();
			// NOTE:: This should be in the script
			/*
			   double getAvgDepth(string smID, string asmDir)
			   {
			// To do: check this file naming convention
			string sDepthFile = asmDir + "/CNV/depthOfCoverage_100000-" + smID + ".tsv";

			ifstream inFile(sDepthFile.c_str(), ios::in);
			int gcIdx = 4;
			double sum = 0;
			double N = 0;

			if (!inFile.good())
			{
			cerr << "Cannot open  " << sDepthFile<< endl;
			abort();
			}
			while(inFile.good())
			{
			string ln;
			getline(inFile,ln);
			if (!ln.empty())
			{
			vector<string> tokens;
			pFile::tokenizeLine(ln.c_str(), " \t\n", tokens);
			if (tokens[0].at(0) == '>')
			{
			for(unsigned j=0; j<tokens.size() ;++j)
			{
			if (tokens[j].compare("gcCorrectedCvg") == 0)
			{
			gcIdx = j;
			}
			}
			}
			else if (tokens[0].substr(0,3).compare("chr") == 0)
			{
			sum += atof(tokens[gcIdx].c_str());
			N++;
			}
			}
			}
			inFile.close();

			return(sum/(double)N);
			}
			 */


		}
	}
	inFile.close();
	return true;
}


