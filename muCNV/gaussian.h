//
//  Gaussian.h
//  muCNV
//
//  Created by Goo Jun on 12/3/18.
//  Copyright © 2018 Goo Jun. All rights reserved.
//

#ifndef Gaussian_h
#define Gaussian_h

#include <stdio.h>
#include <vector>
#include <iostream>
#include <float.h>

const double PI=3.14159265358979;
const double sqPI=1.77245385090552;
const double invsqrt2pi= 0.398942280401433;
const double log2pi = 0.7981798684;

class Gaussian
{
public:
    double Mean;
    double Stdev;
    double Alpha;
    double pdf(const double &);
    void set(const double &, const double &);
    void estimate(std::vector<double> &);
    Gaussian();
};

class Gaussian2
{
public:
    // Parameters for two-dimensional Gaussian
    double Mean[2]; // mean
    double Cov[4]; // covariance
    double Prc[4]; // precision (inverse of covariance)
    double Det; // determinant
    double Alpha;
    double pdf(const double&, const double&);
    void set(const double &, const double &, const double &, const double &, const double &, const double &);
    void estimate(std::vector<double> &, std::vector<double> &);
    double logpdf(const double&, const double&);
    void update(); // update precision matrix
    
    Gaussian2();
};

double normpdf(double, Gaussian&);
double mean(std::vector<double>&);
double stdev(std::vector<double>&, double);

double det(double*);

void printCluster(std::vector<Gaussian> &C);

#endif /* Gaussian_h */
