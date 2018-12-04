//
//  genotyper.h
//  muCNV
//
//  Created by Goo Jun on 12/3/18.
//  Copyright © 2018 Goo Jun. All rights reserved.
//

#ifndef genotyper_h
#define genotyper_h

#include "sv.h"
#include "gaussian.h"
#include <math.h>

class Genotyper
{
public:
    void call_tmp(sv &, svdata &, svgeno &, std::vector<double> &, std::vector<double> &, std::vector<double> &);
    void call_del_tmp(sv &, svdata &, svgeno &, std::vector<double> &, std::vector<double> &, std::vector<double> &);
    void call_dup_tmp(sv &, svdata &, svgeno &, std::vector<double> &, std::vector<double> &, std::vector<double> &);
    void call_del(sv &, svdata &, svgeno &, std::vector<double> &, std::vector<double> &, std::vector<double> &);
    void call_cnv(sv &, svdata &, svgeno &, std::vector<double> &, std::vector<double> &, std::vector<double> &);
    void call_inv(sv &, svdata &, svgeno &, std::vector<double> &, std::vector<double> &);
    
    int assign(double, std::vector<Gaussian> &);
    void copyComps(std::vector<Gaussian> &, std::vector<Gaussian> &);
    //    void call_del(sv&, std::vector<double>&, std::vector<double>&, std::vector<int>&, outvcf&, std::vector<double>&);
    //    void call_cnv(sv&, std::vector<double>&, std::vector<double>&, std::vector<int>&, outvcf&, std::vector<double>&);
};

bool ordered(std::vector<Gaussian> &);


#endif /* genotyper_h */
