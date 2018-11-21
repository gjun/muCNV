//
//  gcContent.cpp
//  muCNV
//
//  Created by Goo Jun on 9/17/18.
//  Copyright © 2018 Goo Jun. All rights reserved.
//

#include <stdio.h>
#include <iostream>
#include <fstream>
#include "gc_content.h"
#include "debug.h"
#include <math.h>

void readmagic(std::ifstream &F)
{
    char buf[100];
    F.read(reinterpret_cast<char *>(buf), 7);
    buf[7] = '\0';
    if (strcmp(buf, "mCNVMGC"))
    {
        std::cerr << "Error: GC content file is corrupt." << std::endl;
        exit(1);
    }
}

void GcContent::initialize(std::string &gcFile)
{
	// gcFile: filename for GC content file

    std::ifstream inFile(gcFile.c_str(), std::ios::in | std::ios::binary);
    if (!inFile.good())
    {
        std::cerr << "Error: cannot open GC file."<< std::endl;
        exit(1);
    }
    
    readmagic(inFile);

    // read number of chrs
    inFile.read(reinterpret_cast <char *> (&num_chr), sizeof(uint8_t));
    
    DMSG((int)num_chr << " chromosomes in GC content file.");
    
    // read size of chrs 
    chr_size.resize(num_chr+1);
    chr_size[0] = 0;
    chr_offset.resize(num_chr+1);
    chr_offset[0] = 0;
    
    
    for(int i=1;i<=num_chr;++i)
    {
        inFile.read(reinterpret_cast <char *> (&chr_size[i]), sizeof(uint32_t));
        chr_offset[i] = chr_offset[i-1] + chr_size[i-1];
        DMSG("Chr " << i << " size: " << chr_size[i] <<  " offset: " << chr_offset[i]);
    }
    
    // read size of GC-interval bin
    inFile.read(reinterpret_cast <char *> (&binsize), sizeof(uint16_t));
    DMSG( "Bin size: " << (int) binsize);
    
    // read number of GC bins
    inFile.read(reinterpret_cast <char *> (&num_bin), sizeof(uint16_t));
    DMSG( "Num_bin : " << num_bin );
    
    // read number of total sampled intervals
    inFile.read(reinterpret_cast <char *> (&num_interval), sizeof(uint16_t));
    DMSG( "Num_interval : " << num_interval);
    
    regions.resize(num_interval);
    
    readmagic(inFile);
    
    // read 'sampled' intervals (chr, start, end)
    for(int i=0; i<num_interval; ++i)
    {
        uint8_t c;
        uint32_t pos1, pos2;
        uint8_t gc;
        inFile.read(reinterpret_cast<char *>(&c), sizeof(uint8_t));
        inFile.read(reinterpret_cast<char *>(&pos1), sizeof(uint32_t));
        inFile.read(reinterpret_cast<char *>(&pos2), sizeof(uint32_t));
        inFile.read(reinterpret_cast<char *>(&gc), sizeof(uint8_t));
        
        regions[i].chrnum = c;
        regions[i].pos = pos1;
        regions[i].end = pos2;
        regions[i].gcbin = gc;
    }
    
	DMSG("Sampled intervals read"); 
    readmagic(inFile);

    // read GC content for each (bin size)-bp interval for each chromosome
    // Current : 400-bp with 200bp overlap
    // std::cerr << "Currnet position : " << inFile.tellg() << std::endl;

    gc_array.resize(num_chr+1);
    gc_array[0] = NULL;

    for(int i=1; i<=num_chr; ++i)
    {
		// Fixed binsize
        int num_bins_in_chr = ceil(chr_size[i] / 200.0);

        gc_array[i] = (uint8_t *) calloc(num_bins_in_chr + 1, sizeof(uint8_t));
		DMSG("GC array " << i << " is allocated");

        inFile.read(reinterpret_cast<char *>(gc_array[i]), sizeof(uint8_t)*num_bins_in_chr);
        // TEMPORARY, TODO: update GC-content file to include end-of-chr bin
		gc_array[i][num_bins_in_chr] = 0;
        readmagic(inFile);
        
        if (!inFile.good())
        {
            std::cerr << "Cannot finish reading GC content file." <<std::endl;
            exit(1);
        }
        DMSG("Chr " << i << " GC content array loaded for "<<  num_bins_in_chr << " segments. ");
    }
    
    //    std::cerr << "Currnet position : " << inFile.tellg() << std::endl;
    
    for(int i=0;i<num_bin;++i)
    {
        if (!inFile.good())
        {
            std::cerr << "Cannot finish reading GC content file." <<std::endl;
            exit(1);
        }
        double v;
        inFile.read(reinterpret_cast<char *>(&v), sizeof(double));
        DMSG("Bin " << i << " GC content proportion: "<< v);
        gc_dist.push_back(v);
    }
    readmagic(inFile);
    inFile.close();
}
