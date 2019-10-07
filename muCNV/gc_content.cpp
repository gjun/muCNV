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
#include <string.h>

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

// GC_Content File Structure
// ----------
// MAGIC
// uint8_t : Number of Chr
// uint32_t * (Number of Chr) : Size of each Chr.
// uint16_t : bin width (how much bp used to average GC content), 400bp
// uint16_t : bin_dist (how much distance in bp between recorded GC contents, 100bp
// uint16_t : number of bins in GC content curve, default: 100 (0 means GC content from 0 to 1%)
// MAGIC
// (# Chr)
//  | (Chr Size) * uint8_t : GC content for each genomic position (every bin_dist-th bp)
//  | MAGIC
// double * (number of bins) : fraction of genomic bin_dist intervals in each GC-bin
// MAGIC

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
    n_interval.resize(num_chr+1);
    n_interval[0] = 0;
    
    
    for(int i=1;i<=num_chr;++i)
    {
        inFile.read(reinterpret_cast <char *> (&chr_size[i]), sizeof(uint32_t));
        chr_offset[i] = chr_offset[i-1] + chr_size[i-1];
        DMSG("Chr " << i << " size: " << chr_size[i] <<  " offset: " << chr_offset[i]);
        inFile.read(reinterpret_cast <char *> (&n_interval[i]), sizeof(uint32_t));
    }
    
    // read size of GC-interval bin
    inFile.read(reinterpret_cast <char *> (&interval_width), sizeof(uint16_t));
    DMSG( "Length of genomic regions that GC content is averaged over: " << (int) interval_width);
    
    // read distance between GC sampling points
    inFile.read(reinterpret_cast <char *> (&interval_dist), sizeof(uint16_t));
    DMSG( "Distance between GC content measuring points: " << interval_dist);

    // read number of total sampled intervals
    inFile.read(reinterpret_cast <char *> (&num_bin), sizeof(uint16_t));
    DMSG( "Num_bin : " << num_bin );
    
    readmagic(inFile);

    /*
    regions.resize(num_interval);
    
    
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
*/
    // read GC content for each (bin size)-bp interval for each chromosome
    // Current : 400-bp with 200bp overlap
    // std::cerr << "Currnet position : " << inFile.tellg() << std::endl;

    gc_array.resize(num_chr+1);
    gc_array[0] = NULL;

    for(int i=1; i<=num_chr; ++i)
    {
		// Fixed binsize
        int num_bins_in_chr = ceil((chr_size[i] + 1.0)/ (double)interval_dist); // Number of intervals in a chromosome
        gc_array[i] = (uint8_t *) calloc(num_bins_in_chr, sizeof(uint8_t));
        
		DMSG("GC array " << i << " is allocated");

        inFile.read(reinterpret_cast<char *>(gc_array[i]), sizeof(uint8_t)*num_bins_in_chr);
        readmagic(inFile);

        // TEMPORARY, TODO: update GC-content file to include end-of-chr bin
		// gc_array[i][num_bins_in_chr] = 0;
        
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
        DMSG("Bin " << i << " Fraction : "<< v);
        bin_fraction.push_back(v);
    }
    readmagic(inFile);
    inFile.close();
}

void GcContent::initialize_quick(std::string &gcFile)
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
    n_interval.resize(num_chr+1);
    n_interval[0] = 0;
    
    
    for(int i=1;i<=num_chr;++i)
    {
        inFile.read(reinterpret_cast <char *> (&chr_size[i]), sizeof(uint32_t));
        chr_offset[i] = chr_offset[i-1] + chr_size[i-1];
        DMSG("Chr " << i << " size: " << chr_size[i] <<  " offset: " << chr_offset[i]);
        inFile.read(reinterpret_cast <char *> (&n_interval[i]), sizeof(uint32_t));
    }
    
    // read size of GC-interval bin
    inFile.read(reinterpret_cast <char *> (&interval_width), sizeof(uint16_t));
    DMSG( "Length of genomic regions that GC content is averaged over: " << (int) interval_width);
    
    // read distance between GC sampling points
    inFile.read(reinterpret_cast <char *> (&interval_dist), sizeof(uint16_t));
    DMSG( "Distance between GC content measuring points: " << interval_dist);

    // read number of total sampled intervals
    inFile.read(reinterpret_cast <char *> (&num_bin), sizeof(uint16_t));
    DMSG( "Num_bin : " << num_bin );
    
    readmagic(inFile);

    inFile.close();
}

double GcContent::get_gc_content(int c, int startpos, int endpos)
{
    if (c < 1 || c > num_chr)
    {
        return 0.0;
    }
    
    // TODO: check this logic
    
    if (round(startpos/(double)interval_dist) == round(endpos/(double)interval_dist))
    {
        return ((double)gc_array[c][(int)round(startpos/(double)interval_dist)-1]/ (double)num_bin + 0.005);
    }

    double gc_sum = 0 ;
    int gc_cnt = 0;

    for(int i=round(startpos/(double)interval_dist)-1; i <= round(endpos/(double)interval_dist)-1; ++i)
    {
        if (gc_array[c][i] >=0 && gc_array[c][i]<num_bin)
        {
            gc_sum += (double)((gc_array[c][i]) / (double)num_bin) + 0.005 ;
            gc_cnt ++;
        }
    }
    if (gc_cnt>0)
        return ((double)gc_sum/gc_cnt);
    else
        return 0.0;
}
