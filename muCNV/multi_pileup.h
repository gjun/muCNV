//
//  multi_pileup.hpp
//  muCNV
//
//  Created by Goo Jun on 11/25/18.
//  Copyright © 2018 Goo Jun. All rights reserved.
//

#ifndef multi_pileup_hpp
#define multi_pileup_hpp

#include <stdio.h>
#include "pileup.h"

using std::string;

class MultiPileup : public Pileup
{
/*    void read();
    void read_depth100();
    void read_var_depth();
    void read_pair_split();
  */
    // Get GC corrected depth for chr / pos
  //  double gcCorrected(double, int, int);
    
};

#endif /* multi_pileup_hpp */
