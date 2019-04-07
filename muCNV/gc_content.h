//
//  gcContent.h
//
//  muCNV
//
//  Created by Goo Jun on 11/13/18.
//  Copyright © 2018 Goo Jun. All rights reserved.
//

#ifndef gcContent_h
#define gcContent_h

#include <string>
#include <vector>
#include "sv.h"

class gcint : public sv
{
public:
    uint8_t gcbin;
};


class GcContent
{
public:
    void initialize(std::string &); // filename for GC content, populate all std::vectors
    double get_gc_content(int, int, int);
    uint16_t num_bin; // Number of GC bin
    uint8_t num_chr; //number of chrs
    uint16_t interval_width; // Size of GC averaged over region (bp)
    uint16_t interval_dist; // distance between GC content measured points
    std::vector<size_t> chr_offset;
//    std::vector<gcint> regions; // Double array to store list of regions for each GC bin -- non-overlapping
    std::vector<double> bin_fraction; // Array to store fraction of Ref genome for each GC content bin
    std::vector<uint32_t> chr_size; // size of chromosomes (bp)
    std::vector<uint32_t> n_interval;
    std::vector<uint8_t *> gc_array; // Array to store "GC bin number" for every 400-bp (?) interval of reference genome
};


#endif /* gcContent_h */
