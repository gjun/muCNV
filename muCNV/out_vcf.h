//
//  out_vcf.h
//  muCNV
//
//  Created by Goo Jun on 12/3/18.
//  Copyright © 2018 Goo Jun. All rights reserved.
//

#ifndef out_vcf_h
#define out_vcf_h

#include <stdio.h>
#include <vector>
#include "gaussian.h"
#include "sv.h"

class outvcf
{
public:
    FILE *fp;
    int varcnt;
    void open(std::string&);
    void close();
    void print(std::string &);
    void write_header(std::vector<std::string>&);
    void write_del(sv&, std::vector<int>&, std::vector<int>&, int, int, std::vector<double>&, std::vector<double>&, std::vector<Gaussian>&, double, bool);
    void write_cnv(sv&, std::vector<int>&, std::vector<int>&, int, int, std::vector<double>&, std::vector<double>&, std::vector<Gaussian>&, double, bool);
};


#endif /* out_vcf_h */
