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

void split(const char* s, const char* delims, std::vector<string>& tokens)
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

string svTypeName(svType t)
{
    string S[6] = {"DEL", "DUP", "INV", "CNV", "INS", "BND"};
    return S[t];
}

svType svTypeNum(int t)
{
    switch(t)
    {
        case 0:
            return DEL;
            break;
        case 1:
            return DUP;
            break;
        case 2:
            return INV;
            break;
        case 3:
            return CNV;
            break;
        case 4:
            return INS;
            break;
        case 5:
            return BND;
            break;
    }
    return DEL;
}

int median(std::vector<int> &L)
{

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
	void vprint(std::vector<T> x)
{
	std::cerr << "(";
	for(unsigned i=0; i<x.size()-1; ++i)
	{
		std::cerr << x[i] << " ";
	}
	std::cerr << x.back() << ")" << std::endl;
}


double MAX(std::vector<double>& x)
{
	double ret = -1.0*DBL_MAX;
	for(unsigned i=0;i<x.size();++i)
	{
		if (ret<x[i])
		{
			ret = x[i];
		}
	}
	return ret;
}

double MIN(std::vector<double>& x)
{
	double ret = DBL_MAX;
	for(unsigned i=0;i<x.size();++i)
	{
		if (ret>x[i])
		{
			ret = x[i];
		}
	}
	return ret;
}


double det(double* M)
{
	return (M[0]*M[3]-M[1]*M[2]);
}


double mean(std::vector<double>& x)
{
	double sum = 0;
	for(unsigned j=0; j<x.size(); ++j)
	{
		sum += x[j];
	}
	return sum/x.size();
}

double stdev(std::vector<double>& x, double M)
{

	if (x.size()>1)
	{
		double sumsq = 0;
		for(unsigned j=0; j<x.size(); ++j)
		{
			sumsq += (x[j]-M) *(x[j]-M);
		}
		return sqrt(sumsq / x.size());
	}
	else
	{
		return 0.01;
	}
}
