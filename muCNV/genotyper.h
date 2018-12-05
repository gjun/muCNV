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

class svgeno
{
public:
    bool b_biallelic;
    bool b_pass;
    bool b_dump;
    
    bool dp_flag;
    bool pos_flag;
    bool neg_flag;
    
    bool cnv_pos_flag;
    bool cnv_neg_flag;
    bool inv_pos_flag;
    bool inv_neg_flag;
    
    double p_overlap;
    std::string info;
    int n_sample;
    
    int ac;
    int ns;
    std::vector<double> bic;
    
    std::vector<Gaussian> Comps;
    
    std::vector<int> gt; // bi-allelic genotype
    std::vector<int> cn; // copy number
    
    void initialize(int N)
    {
        n_sample = N;
        Comps.clear();
        
        ns = 0;
        ac = 0;
        p_overlap = 1;
        
        gt.resize(n_sample, -1);
        cn.resize(n_sample, -1);
        bic.resize(20, 0);
        
        for(int i=0;i<n_sample; ++i)
        {
            gt[i] = -1;
            cn[i] = -1;
        }
        for(int i=0;i<20;++i)
        {
            bic[i] = 0;
        }
        
        b_biallelic = false;
        b_pass = false;
        b_dump = true;
        
        dp_flag = false;
        pos_flag = false;
        neg_flag = false;
        // For inversions
        cnv_pos_flag = false;
        cnv_neg_flag = false;
        inv_pos_flag = false;
        inv_neg_flag = false;
        info = "";
    };
    void print (sv &, svdata &, std::string &, std::vector<double> &);
};


class Genotyper
{
public:
    void call_del(sv &, svdata &, svgeno &, std::vector<double> &, std::vector<double> &, std::vector<double> &);
    void call_cnv(sv &, svdata &, svgeno &, std::vector<double> &, std::vector<double> &, std::vector<double> &);
    
    int assign(double, std::vector<Gaussian> &);
    void copyComps(std::vector<Gaussian> &, std::vector<Gaussian> &);
};

bool ordered(std::vector<Gaussian> &);


#endif /* genotyper_h */
