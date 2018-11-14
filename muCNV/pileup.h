//
//  pileup.h
//  muCNV
//
//  Created by Goo Jun on 11/13/18.
//  Copyright © 2018 Goo Jun. All rights reserved.
//

#ifndef pileup_h
#define pileup_h

#include "sv.h"
#include <string>
#include <vector>

class pileup
{
public:
    void write(std::string &, std::vector<sv> &);
    void write_text(std::string &, std::vector<sv> &);
};


#endif /* pileup_h */
