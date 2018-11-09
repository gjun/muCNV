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

bool in_centrome(int chrnum, int pos)
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
    if (pos >= centro[chrnum-1][0] && pos <= centro[chrnum-1][1])
        return true;
	else
	    return false;
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
//  int hetero[2] = {51078349, 54425074};

//  if (S.pos >= centro[S.chrnum-1][0]-300000 && S.pos <= centro[S.chrnum-1][1]+300000)
    if (S.pos >= centro[S.chrnum-1][0] && S.pos <= centro[S.chrnum-1][1])
        return true;
    if (S.end >= centro[S.chrnum-1][0] && S.end <= centro[S.chrnum-1][1])
        return true;
//  if (S.chrnum == 7  && S.pos >= hetero[0] && S.pos <= hetero[1]+1000000)
//      return true;
//  if (S.chrnum == 7  && S.end >= hetero[0] && S.end <= hetero[1])
//      return true;
    return false;
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
	void vprint(std::vector<T> &x)
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
