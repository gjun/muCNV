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
#include "common.h"
#include <math.h>
#include <string.h>
#include <algorithm>

void split(const char* s, const char* delims, std::vector<std::string>& tokens)
{
	const char* p = s;
	const char* c = p;
	int ndelims = (int)strlen(delims);
	int i;
	tokens.clear();
	
	//fprintf(stderr,"s = '%s', strlen(s) = %d, delims = '%s'\n",s,(int)strlen(s), delims);
	while( *p != '\0' )
	{
		for(i=0; i < ndelims; ++i)
		{
			if ( *p == delims[i] )
				break;
		}
		if ( i != ndelims ) { // delimiter found
			if ( c < p )  { // unless delimiter is consencutive
							//string s1(c,p-c);
				tokens.push_back(std::string(c,p-c));
			}
			c = p+1;
		}
		++p;
	}
	if ( c < p ) {
		tokens.push_back(std::string(c,p-c));
	}
}

double average(std::vector<double> &L)
{
    if (L.size() == 0) return 0;
    double sum = 0;

    for(int i=0; i< (int)L.size(); ++i)
    {
        sum += L[i];
    }
    return sum/L.size();
}


int median(std::vector<int> &L)
{
    // TODO: use std::n_th
	int med=0;
	if (L.size() > 0)
	{
		std::sort(L.begin(), L.end());

		if (L.size()%2)
		{
			int m = (int) (L.size()-1)/2;
			med = L[m];
		}
		else
		{
			int m = (int) L.size()/2;
			med = (int)(L[m-1] + L[m])/2;
		}
	}
	return med;
}

template <class T>
	void vprint(std::vector<T> &x)
{
	std::cerr << "(";
	for(unsigned i=0; i<x.size()-1; ++i)
	{
		std::cerr << x[i] << " ";
	}
	std::cerr << x.back() << ")" << std::endl;
}

