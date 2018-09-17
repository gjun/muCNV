//
//  gcContent.cpp
//  muCNV
//
//  Created by Goo Jun on 9/17/18.
//  Copyright © 2018 Goo Jun. All rights reserved.
//

#include <stdio.h>
#include "muCNV.h"



void readmagic(ifstream &F)
{
    char buf[100];
    F.read(reinterpret_cast<char *>(buf), 7);
    buf[7] = '\0';
    if (strcmp(buf, "mCNVMGC"))
    {
        cerr << "Error: GC content file is corrupted." << endl;
        exit(1);
    }
}

void gcContent::initialize(string &gcFile)
{// filename for GC content, populate all vectors
    ifstream inFile(gcFile.c_str(), std::ios::in | std::ios::binary);
    if (!inFile.good())
    {
        cerr << "Error: cannot open GC file."<< endl;
        exit(1);
    }
    
    readmagic(inFile);
    // read number of chrs
    inFile.read(reinterpret_cast <char *> (&num_chr), sizeof(uint8_t));
    
    //    cerr << (int)num_chr << " chromosomes in GC content file." << endl;
    
    // read size of chrs
    chrSize.resize(num_chr+1);
    chrSize[0] = 0;
    chrOffset.resize(num_chr+1);
    chrOffset[0] = 0;
    
    // gc_array[0] = NULL;
    
    for(int i=1;i<=num_chr;++i)
    {
        inFile.read(reinterpret_cast <char *> (&chrSize[i]), sizeof(uint32_t));
        chrOffset[i] = chrOffset[i-1] + chrSize[i-1];
        //cerr << "Chr " << i << " size: " << chrSize[i] <<  "offset: " << chrOffset[i] << endl;
    }
    
    // read size of GC-interval bin
    inFile.read(reinterpret_cast <char *> (&binsize), sizeof(uint16_t));
    cerr << "Bin size: " << (int) binsize << endl;
    
    // read number of GC bins
    inFile.read(reinterpret_cast <char *> (&num_bin), sizeof(uint16_t));
    //    cerr << "Num_bin : " << num_bin << endl;
    
    // read number of total intervals
    inFile.read(reinterpret_cast <char *> (&total_bin), sizeof(uint16_t));
    //    cerr << "Total bin : " << total_bin << endl;
    
    regions.resize(total_bin);
    
    readmagic(inFile);
    
    gc_array.resize(num_chr+1);
    
    // read 'sampled' intervals (chr, start, end)
    for(int i=0; i<total_bin; ++i)
    {
        uint8_t c;
        uint32_t pos1, pos2;
        uint8_t gc;
        inFile.read(reinterpret_cast<char *>(&c), sizeof(uint8_t));
        inFile.read(reinterpret_cast<char *>(&pos1), sizeof(uint32_t));
        inFile.read(reinterpret_cast<char *>(&pos2), sizeof(uint32_t));
        inFile.read(reinterpret_cast<char *>(&gc), sizeof(uint8_t));
        
        regions[i].chrnum = c;
        //        regions[i].chr = "chr" + to_string(c);
        regions[i].pos = pos1;
        regions[i].end = pos2;
        regions[i].gcbin = gc;
    }
    
    readmagic(inFile);
    // read GC content for each (bin size)-bp interval for each chromosome
    // Current : 400-bp with 200bp overlap
    // cerr << "Currnet position : " << inFile.tellg() << endl;
    for(int i=1; i<=num_chr; ++i)
    {
        int N = ceil((chrSize[i] / (double)binsize)*2.0) ;
        gc_array[i] = (uint8_t *) calloc(N, sizeof(uint8_t));
        inFile.read(reinterpret_cast<char *>(gc_array[i]), sizeof(uint8_t)*N);
        readmagic(inFile);
        
        if (!inFile.good())
        {
            cerr << "Cannot finish reading GC content file." <<endl;
            exit(1);
        }
        //        cerr << "Chr " << i << " GC content array loaded for "<<  N << " segments. " <<  endl;
    }
    //    cerr << "Currnet position : " << inFile.tellg() << endl;
    
    for(int i=0;i<num_bin;++i)
    {
        if (!inFile.good())
        {
            cerr << "Cannot finish reading GC content file." <<endl;
            exit(1);
        }
        double v;
        inFile.read(reinterpret_cast<char *>(&v), sizeof(double));
        //        cerr << "Bin " << i << " GC content proportion: "<< v << endl;
        gc_dist.push_back(v);
    }
    readmagic(inFile);
    inFile.close();
}
