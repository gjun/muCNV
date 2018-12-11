//
//  Gaussian.cpp
//  muCNV
//
//  Created by Goo Jun on 7/25/17.
//  Copyright © 2017 Goo Jun. All rights reserved.
//


#include "gaussian.h"
//#include <math.h>
#include <cmath>


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


Gaussian::Gaussian()
{
	Mean=0;
	Stdev=1;
}

void Gaussian::set(const double &m, const double &s)
{
	Mean = m;
	Stdev = s;
}

void Gaussian::estimate(std::vector<double> &x)
{
	Mean = mean(x);
	Stdev = stdev(x, Mean);
}


double Gaussian::pdf(const double& x)
{
	double z = (x-Mean)/Stdev;
	double val =  (invsqrt2pi * exp(-0.5*z*z) /Stdev);
	if (std::isnan(val))
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

void Gaussian2::set(const double &m0, const double & m1, const double & c0, const double & c1, const double & c2, const double & c3)
{
	Mean[0] = m0;
	Mean[1] = m1;
	Cov[0] = c0;
	Cov[1] = c1;
	Cov[2] = c2;
	Cov[3] = c3;
	
	update();
}

void Gaussian2::update()
{
	Det = det(Cov);
	
	while (Det<1e-20)
	{
		Cov[0] += 1e-10;
		Cov[3] += 1e-10;
		Det = det(Cov);
	}
	Prc[0] = Cov[3]/Det;
	Prc[3] = Cov[0]/Det;
	Prc[1] = Prc[2] = -1.0*Cov[1]/Det;
}

void Gaussian2::estimate(std::vector<double> &x, std::vector<double> &y)
{
	int n = (int)x.size();
	double sum_x=0;
	double sum_y=0;
	double sum_xx=0;
	double sum_yy=0;
	double sum_xy=0;
	
	for(int i=0;i<n;++i)
	{
		sum_x+=x[i];
		sum_y+=y[i];
		sum_xx+=x[i]*x[i];
		sum_xy+=x[i]*y[i];
		sum_yy+=y[i]*y[i];
	}
	
	Mean[0] = sum_x/(double)n;
	Mean[1] = sum_y/(double)n;
	Cov[0] = sum_xx/(double)n - (Mean[0]*Mean[0]);
	Cov[1] = sum_xy/(double)n - (Mean[0]*Mean[1]);
	Cov[2] = Cov[1];
	Cov[3] = sum_yy/(double)n - (Mean[1]*Mean[1]);
	update();
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
	double val =  (invsqrt2pi * exp(-0.5*z*z) /C.Stdev);
	if (std::isnan(val))
	{
		return 0;
	}
	else
	{
		return val;
	}
}


double lognormpdf(double x, Gaussian& C)
{
	double z = (x-C.Mean)/C.Stdev;
	if (C.Stdev>0)
	{
		return (-0.9189385332 -0.5*z*z - log(C.Stdev));
	}
	else
	{
		return 0;
	}
}




