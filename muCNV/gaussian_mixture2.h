//
//  gaussian_mixture2.h
//  muCNV
//
//  Copyright © 2018 Goo Jun. All rights reserved.
//

#ifndef gaussian_mixture2_h
#define gaussian_mixture2_h

#include <stdio.h>
#include <vector>
#include "gaussian.h"
#include "debug.h"

class GaussianMixture2
{
public:
    std::vector<Gaussian2> Comps;
    int n_comp;
    double bic;
    double aic;
    double p_overlap;
	int zeroidx;
    
    //2-D without weights
    void EM2(std::vector<double>&, std::vector<double>&);
    void EM2_select(std::vector<double>&, std::vector<double>&, std::vector<bool> &);

    void estimate(std::vector<double> &, std::vector<double> &, std::vector<int> &, int);
    void estimate_select(std::vector<double> &, std::vector<double> &, std::vector<int> &, std::vector<bool> &, int);

    // K-means
    void KM2(std::vector<double>&, std::vector<double>&, bool);

    int assign_copynumber(double, double);
    int assign_dpcnt_copynumber(double, double);

    bool ordered();
    bool r_ordered();
    void updateAICBIC(std::vector<double>& , std::vector<double>& );
    void updateAICBIC_select(std::vector<double>& , std::vector<double>&, std::vector<bool> &);

    double BayesError();
	void print(FILE *);
    std::string print_str();

    GaussianMixture2 () {bic = DBL_MAX; p_overlap = 1.0; n_comp = 0; zeroidx = -1;};
    GaussianMixture2 (const GaussianMixture2 &);
    GaussianMixture2 (std::vector<double> &, std::vector<double> &);
    GaussianMixture2 (std::vector<double> &, std::vector<double> &, std::vector<double> &);
    GaussianMixture2& operator = ( const GaussianMixture2& gmix);
    
};

#endif /* gaussian_mixture2_h */
