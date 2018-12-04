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

class GaussianMixture
{
    std::vector<Gaussian> Comps;
    int n_comp;
    
    //1-D with weights
    void EM(std::vector<double>&, std::vector<double>&);
    
    GaussianMixture(int) {};

};

class GaussianMixture2
{
public:
    std::vector<Gaussian2> Comps;
    int n_comp;
    
    //2-D without weights
    void EM2(std::vector<double>&, std::vector<double>&);

    // TODO: keep LLK value
    
    GaussianMixture2(int) {};
};

#endif /* gaussian_mixture_h */
