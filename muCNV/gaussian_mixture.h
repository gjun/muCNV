//
//  mix_gaussian.hpp
//  muCNV
//
//  Created by Goo Jun on 12/3/18.
//  Copyright © 2018 Goo Jun. All rights reserved.
//

#ifndef gaussian_mixture_h
#define gaussian_mixture_h

#define MAX_P_OVERLAP 0.3

#include <stdio.h>
#include <vector>
#include "gaussian.h"

// Todo: make GaussianMixture and GaussianMixture2 to be inherited from a generic mixture of Gaussians class

class GaussianMixture
{
public:
    std::vector<Gaussian> Comps;
    int n_comp;
    double llk;
    double bic;
    double p_overlap;
    
    //1-D with weights
    void wEM(std::vector<double>&, std::vector<double>&);
    
    //1-D without weights
    void EM(std::vector<double>&);
    
    bool ordered();
    bool r_ordered();
    
    // Copy constructor
    GaussianMixture () {llk = -DBL_MAX; bic = DBL_MAX; p_overlap = 1.0; n_comp = 0;};
    GaussianMixture (const GaussianMixture &);
    GaussianMixture (std::vector<double> &, std::vector<double> &);
    GaussianMixture& operator = ( const GaussianMixture& gmix);

    int assign_copynumber(double);
};

class GaussianMixture2
{
public:
    std::vector<Gaussian2> Comps;
    int n_comp;
    double llk;
    double bic;
    double p_overlap;
    
    //2-D without weights
    void EM2(std::vector<double>&, std::vector<double>&);

    int assign_copynumber(double, double);
    bool ordered();
    bool r_ordered();
    
    GaussianMixture2 () {llk = -DBL_MAX; bic = DBL_MAX; p_overlap = 1.0; n_comp = 0;};
    GaussianMixture2 (const GaussianMixture2 &);
    GaussianMixture2 (std::vector<double> &, std::vector<double> &);
    GaussianMixture2& operator = ( const GaussianMixture2& gmix);
    
};

#endif /* gaussian_mixture_h */
