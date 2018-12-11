//
//  cluster.h
//  muCNV
//
//  Created by Goo Jun on 12/3/18.
//  Copyright © 2018 Goo Jun. All rights reserved.
//

#ifndef cluster_h
#define cluster_h

#include <stdio.h>
#include <vector>
#include "sv.h"

void merge_svs(std::vector<sv> &, std::vector<int> &);
void cluster_svs(std::vector<sv>&, std::vector< std::vector<sv> > &, double);


#endif /* cluster_h */
