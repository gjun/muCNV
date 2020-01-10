//
//  mix_gaussian.hpp
//  muCNV
//
//  Created by Goo Jun on 12/3/18.
//  Copyright © 2018 Goo Jun. All rights reserved.
//

#ifndef gaussian_mixture_h
#define gaussian_mixture_h

#include <stdio.h>
#include <vector>
#include "gaussian.h"

// Todo: make GaussianMixture and GaussianMixture2 to be inherited from a generic mixture of Gaussians class

class GaussianMixture
{
public:
    std::vector<Gaussian> Comps;
    int n_comp;
	int zeroidx;
    double bic;
    double aic;
    double p_overlap;
    
    //1-D without weights
    void EM(std::vector<double>&);
    void EM_select(std::vector<double>&, std::vector<bool> &);

	// 1-D K-Means
    void KM(std::vector<double>&, bool);

	void print(FILE *);
    
    bool ordered();
    bool r_ordered();
    double BIC(std::vector<double>& );
    double BIC_select(std::vector<double>&, std::vector<bool>&);

    void updateAICBIC_select(std::vector<double>&, std::vector<bool> &);
    double AIC(std::vector<double>& );

    double BayesError();
    void estimate(std::vector<double> &, std::vector<int> &, int);

    
	// Constructors
    GaussianMixture () { bic = DBL_MAX; aic = DBL_MAX; p_overlap = 1.0; zeroidx = -1; n_comp = 0;};
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
    double bic;
    double p_overlap;
	int zeroidx;
    
    //2-D without weights
    void EM2(std::vector<double>&, std::vector<double>&);
    void EM2_select(std::vector<double>&, std::vector<double>&, std::vector<bool> &);

    // K-means
    void KM2(std::vector<double>&, std::vector<double>&, bool);

    int assign_copynumber(double, double);
    bool ordered();
    bool r_ordered();
    double BIC(std::vector<double>& , std::vector<double>& );
    double BIC_select(std::vector<double>& , std::vector<double>&, std::vector<bool> &);

    double BayesError();
	void print(FILE *);

    GaussianMixture2 () {bic = DBL_MAX; p_overlap = 1.0; n_comp = 0; zeroidx = -1;};
    GaussianMixture2 (const GaussianMixture2 &);
    GaussianMixture2 (std::vector<double> &, std::vector<double> &);
    GaussianMixture2& operator = ( const GaussianMixture2& gmix);
    
};

#endif /* gaussian_mixture_h */
