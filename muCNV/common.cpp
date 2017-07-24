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

typedef pair<uint64_t, uint64_t> interval_t;


template <class T>
	void vprint(vector<T> x)
{
	cerr << "(";
	for(unsigned i=0; i<x.size()-1; ++i)
	{
		cerr << x[i] << " ";
	}
	cerr << x.back() << ")" << endl;
}


double MAX(vector<double>& x)
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

double MIN(vector<double>& x)
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



double normpdf(double x, double m, double s)
{
	double z = (x-m)/s;
	double val =  (sqPI * exp(-0.5*z*z) /s);
	if (isnan(val))
	{
		return 0;
	}
	else
	{
		return val;
	}
}

double normpdf(double x, Gaussian& C)
{
	double z = (x-C.Mean)/C.Stdev;
	double val =  (sqPI * exp(-0.5*z*z) /C.Stdev);
	if (isnan(val))
	{
		return 0;
	}
	else
	{
		return val;
	}
}

double mean(vector<double>& x)
{
	double sum = 0;
	for(unsigned j=0; j<x.size(); ++j)
	{
		sum += x[j];
	}
	return sum/x.size();
}

double stdev(vector<double>& x, double M)
{
	double sumsq = 0;
	for(unsigned j=0; j<x.size(); ++j)
	{
		sumsq += x[j]*x[j];
	}
	return sqrt(sumsq / x.size() - M*M);
}

double BayesError(vector<Gaussian>& Comps)
{
	// Returns maximum Bhattacharyya coefficient (== exp(-D) ) between components
	unsigned n_comp = Comps.size();
	double min_d = DBL_MAX;

	if (n_comp <2 )
	{
		return DBL_MAX;
	}

	// Get minimum distance between all pairs
	for(unsigned i=0; i<n_comp-1; ++i)
	{
		for(unsigned j=i+1; j<n_comp; ++j)
		{
			double s = (Comps[i].Stdev + Comps[j].Stdev)/2.0;
			double d = (Comps[i].Mean-Comps[j].Mean)*(Comps[i].Mean-Comps[j].Mean)/(8.0*s*s) + 0.5*log( s / sqrt(Comps[i].Stdev*Comps[j].Stdev));
			if (d<min_d)
			{
				min_d = d;
			}
		}
	}
	return exp(-1.0*min_d);
}

double BIC(vector<double>& x, vector<Gaussian> &C)
{
	unsigned n_sample = x.size();
	unsigned n_comp = C.size();
	double llk = 0;
	double ret = 0;

	for(unsigned j=0; j<n_sample; ++j)
	{
		double l = 0;
		for(unsigned m=0; m<n_comp; ++m)
		{
			l += C[m].Alpha * normpdf(x[j], C[m].Mean, C[m].Stdev);
		}
		if (l>0)
		{
			llk += log(l);
		}
	}

	ret = -2.0 * llk +  2*n_comp*log(n_sample);

	return ret;
}

void printCluster(vector<Gaussian> &C)
{
	cerr << C.size() << " comps: ";
	for(unsigned i=0;i<C.size();++i)
	{
		cerr << "(" << C[i].Mean << ", " << C[i].Stdev << ") " ;
	}
	cerr << endl;
}

