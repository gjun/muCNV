//
//  sampleStat.hpp
//  muCNV
//
//  Created by Goo Jun on 11/13/18.
//  Copyright © 2018 Goo Jun. All rights reserved.
//

#ifndef sample_stat_h
#define sample_stat_h

#include <vector>

class SampleStat
{
public:
    double avg_dp;
    double std_dp;
    double avg_isize;
    double std_isize;
    double med_isize;
    double avg_rlen;
    
    std::vector<double> gc_factor;
};
#endif /* sample_stat_h */
