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
#include "genotyper.h"

class OutVcf
{
public:
    FILE *fp;
    int varcnt;
    void open(std::string&);
    void close();
    void print(std::string &);
    void write_header(std::vector<std::string>&, std::vector<bool>&, GcContent &);

	void write_sv(sv &, SvData &, SvGeno &, std::vector<int> &);
};


#endif /* out_vcf_h */
