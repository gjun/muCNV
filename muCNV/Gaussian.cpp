//
//  Gaussian.cpp
//  muCNV
//
//  Created by Goo Jun on 7/25/17.
//  Copyright © 2017 Goo Jun. All rights reserved.
//

#include "muCNV.h"
#include <math.h>

Gaussian::Gaussian()
{
	Mean=0;
	Stdev=1;
}

double Gaussian::pdf(const double& x)
{
	double z = (x-Mean)/Stdev;
	double val =  (sqPI * exp(-0.5*z*z) /Stdev);
	if (isnan(val))
	{
		return 0;
	}
	else
	{
		return val;
	}
}

Gaussian2::Gaussian2()
{
	Mean[0]=0;
	Mean[1]=0;
	Cov[0]=1;
	Cov[1]=0;
	Cov[2]=0;
	Cov[3]=1;
	update();
}

void Gaussian2::update()
{
	Det = det(Cov);
	
	if (Det>1e-16)
	{
		Prc[0] = Cov[3]/Det;
		Prc[3] = Cov[0]/Det;
		Prc[1] = Prc[2] = -1.0*Cov[1]/Det;
	}
	else
	{
		Prc[0]=0;
		Prc[1]=0;
		Prc[2]=0;
		Prc[3]=0;
	}
}

double Gaussian2::pdf(const double& x, const double& y)
{
	double val;
	double y0 = x - Mean[0];
	double y1 = y - Mean[1];
	
	val= -0.5*(y0*(y0*Prc[0] + y1*Prc[1]) + y1*(y0*Prc[2]+y1*Prc[3]));
	val= 0.5*exp(val)/(PI*sqrt(Det));
	
	return val;
}

double Gaussian2::logpdf(const double& x, const double& y)
{
	double val;
	double y0 = x - Mean[0];
	double y1 = y - Mean[1];
	
	val= -0.5*(y0*(y0*Prc[0] + y1*Prc[1]) + y1*(y0*Prc[2]+y1*Prc[3]));
	val = val-log2pi-0.5*log(Det);
	
	return val;
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


double BayesError(vector<Gaussian>& Comps)
{
	// Returns maximum Bhattacharyya coefficient (== exp(-D) ) between components
	unsigned n_comp = (unsigned) Comps.size();
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

double BayesError(vector<Gaussian2>& Comps)
{
	int n_comp = (int) Comps.size();
	double min_d = DBL_MAX;
	
	if (n_comp <2 )
	{
		return DBL_MAX;
	}
	
	// Get minimum distance between all pairs
	for(int i=0; i<n_comp-1; ++i)
	{
		for(int j=i+1; j<n_comp; ++j)
		{
			double C[4], P[4];
			for(int k=0;k<4;++k)
			{
				C[k] = (Comps[i].Cov[k] + Comps[j].Cov[k])/2.0;
			}
			double D = det(C);
			if (D>1e-10)
			{
				P[0] = C[3]/D;
				P[3] = C[0]/D;
				P[1] = P[2] = -1.0*C[1]/D;
			}
			else continue;
			
			double m1 = Comps[i].Mean[0] - Comps[j].Mean[0];
			double m2 = Comps[i].Mean[1] - Comps[j].Mean[1];
			
			double d = 0.125*(m1 * P[0] + m2 * P[2])*m1 + (m1 * P[1] + m2 * P[3])*m2 + 0.5*log( D / sqrt(Comps[i].Det*Comps[j].Det));
			
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
	int n_sample = (int)x.size();
	int n_comp = (int)C.size();
	double llk = 0;
	double ret = 0;
	
	for(int j=0; j<n_sample; ++j)
	{
		double l = 0;
		for(int m=0; m<n_comp; ++m)
		{
			l += C[m].Alpha * normpdf(x[j], C[m]);
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

