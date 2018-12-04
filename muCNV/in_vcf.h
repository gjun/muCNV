//
//  in_vcf.h
//  muCNV
//
//  Created by Goo Jun on 12/3/18.
//  Copyright © 2018 Goo Jun. All rights reserved.
//

#ifndef in_vcf_h
#define in_vcf_h

#include <stdio.h>
#include <vector>
#include <string>

#include "hts.h"
#include "tbx.h"
#include "kseq.h"

#include "sv.h"

class invcfs
{
public:
    // std::vector<std::ifstream *> vfs;
    std::vector<htsFile*> vfs;
    std::vector<tbx_t*> tbs;
    std::vector<hts_itr_t *> m_itr;
    std::vector<int> num_id;
    std::vector<int> start_num;
    
    void parse_sv(std::vector<std::string> &, sv &);
    double get_second_value(std::string &);
    void get_value_pair(std::string &, int &, double &);
    int initialize(std::vector<std::string> &, std::vector<std::string> &, std::vector<double> &, std::vector<double> &, std::vector<double> &, std::string &);
    int read_interval_multi(sv& , svdata &, std::string &);
};


void read_svs_from_vcf(std::string &, std::vector<breakpoint> &, std::vector<sv> &);
void read_svs_from_intfile(std::string &, std::vector<breakpoint> &, std::vector<sv> &);
void read_intervals_from_vcf(std::vector<std::string> &, std::vector<std::string> &, std::vector<sv> &);
int read_candidate_vcf(std::ifstream &, sv&, std::string& );

void read_list(std::string &, std::vector<std::string>&);
void read_index(std::string, std::vector<std::string>&, std::vector<std::string>&, std::vector<std::string>&, std::vector<double>&);

void readIndex(std::string, std::vector<std::string>&, std::vector<std::string>&, std::vector<std::string>&, std::vector<std::string>&);

void readDepth(std::vector<std::string>&, std::vector<sv>&, std::vector< std::vector<double> >&, std::vector<double>&);

void write_interval(std::string &, std::vector<sv> &);


#endif /* in_vcf_h */
