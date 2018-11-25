//
//  pileup.cpp
//  muCNV
//
//  Created by Goo Jun on 11/13/18.
//  Copyright © 2018 Goo Jun. All rights reserved.
//

#include <stdio.h>
#include <math.h>

#include <iostream>
#include <fstream>
#include "pileup.h"

void Pileup::initialize()
{
    // TODO: check whether gc is initialized correctly
    gc_factor.resize(gc.num_bin);
    
    depth100.resize(gc.num_chr + 1);
    nbin_100.resize(gc.num_chr + 1);
    nbin_100[0] = 0;
    
    for(int i=1; i<=gc.num_chr; ++i)
    {
        nbin_100[i] = ceil((double)gc.chr_size[i] / 100.0) + 1 ;
        depth100[i] = (uint16_t *) calloc(nbin_100[i], sizeof(uint16_t));
        
        DMSG("chr " << i << " bin size " << nbin_100[i]);
    }
}


int Pileup::write_number(std::ofstream &pileup_file, int number)
{
    pileup_file.write(reinterpret_cast<char*>(&number), sizeof(int));
    return(sizeof(int));
}

int Pileup::write_sample_id(std::ofstream &pileup_file, std::string &sample_id)
{
    char pad[256] = {0};

    // Write sample ID (each with 256 bytes)
    if (sample_id.length() > 255)
    {
        std::cerr << "Error, sample ID " << sample_id << " is too long." << std::endl;
        exit(1);
    }
    
    pileup_file.write(sample_id.c_str(), sample_id.length());
    pileup_file.write(pad, 256-sample_id.length());
    
    return 256;
}

int Pileup::write_sample_stat(std::ofstream &pileup_file, SampleStat &s)
{
    // Write depth and isize stats
    pileup_file.write(reinterpret_cast<char*>(&s.avg_dp), sizeof(double));
    pileup_file.write(reinterpret_cast<char*>(&s.std_dp), sizeof(double));
    pileup_file.write(reinterpret_cast<char*>(&s.avg_isize), sizeof(double));
    pileup_file.write(reinterpret_cast<char*>(&s.std_isize), sizeof(double));
    return sizeof(double) * 4;
}

int Pileup::write_gc_factor(std::ofstream &pileup_file, std::vector<double>& gc_f)
{
    int ret=0;
    // Write GC curve
    for(int i=0;i<gc.num_bin;++i)
    {
        pileup_file.write(reinterpret_cast<char*>(&(gc_f[i])), sizeof(double));
        ret += sizeof(double);
    }
    return ret;
}

int Pileup::write_depth(std::ofstream &pileup_file, uint16_t *dp, int N)
{
    pileup_file.write(reinterpret_cast<char*>(dp), sizeof(uint16_t)*N);
    return (sizeof(uint16_t)*N);
}

int Pileup::write_readpair(std::ofstream &pileup_file, readpair &rp)
{
    int ret = 0;
    pileup_file.write(reinterpret_cast<char*>(&(rp.chrnum)), sizeof(int8_t));
    pileup_file.write(reinterpret_cast<char*>(&(rp.selfpos)), sizeof(int32_t));
    pileup_file.write(reinterpret_cast<char*>(&(rp.matepos)), sizeof(int32_t));
    pileup_file.write(reinterpret_cast<char*>(&(rp.matequal)), sizeof(int8_t));
    pileup_file.write(reinterpret_cast<char*>(&(rp.pairstr)), sizeof(int8_t));
    ret += sizeof(int8_t) + sizeof(int32_t) + sizeof(int32_t) + sizeof(int8_t) + sizeof(int8_t);
    return ret;
}

int Pileup::write_splitread(std::ofstream &pileup_file, splitread& sp)
{
    int ret = 0;
    
    pileup_file.write(reinterpret_cast<char*>(&(sp.chrnum)), sizeof(int8_t));
    pileup_file.write(reinterpret_cast<char*>(&(sp.pos)), sizeof(int32_t));
    pileup_file.write(reinterpret_cast<char*>(&(sp.sapos)), sizeof(int32_t));
    pileup_file.write(reinterpret_cast<char*>(&(sp.firstclip)), sizeof(int16_t));
    pileup_file.write(reinterpret_cast<char*>(&(sp.secondclip)), sizeof(int16_t));
    ret += sizeof(int8_t) + sizeof(int32_t) + sizeof(int32_t) + sizeof(int16_t) + sizeof(int16_t);
    
    return ret;
}

int Pileup::read_number(std::ifstream &pileup_file, int number)
{
    pileup_file.read(reinterpret_cast<char*>(&number), sizeof(int));
    return(sizeof(int));
}

int Pileup::read_sample_id(std::ifstream &pileup_file, char* buf)
{
    // read sample ID (each with 256 bytes)
    pileup_file.read(buf, 256);
    return 256;
}


int Pileup::read_sample_stat(std::ifstream &pileup_file, SampleStat &s)
{
    // read depth and isize stats
    pileup_file.read(reinterpret_cast<char*>(&s.avg_dp), sizeof(double));
    pileup_file.read(reinterpret_cast<char*>(&s.std_dp), sizeof(double));
    pileup_file.read(reinterpret_cast<char*>(&s.avg_isize), sizeof(double));
    pileup_file.read(reinterpret_cast<char*>(&s.std_isize), sizeof(double));
    return sizeof(double) * 4;
}

int Pileup::read_gc_factor(std::ifstream &pileup_file, std::vector<double>& gc_f)
{
    int ret=0;
    // read GC curve
    gc_f.resize(gc.num_bin);
    
    for(int i=0;i<gc.num_bin;++i)
    {
        pileup_file.read(reinterpret_cast<char*>(&(gc_f[i])), sizeof(double));
        ret += sizeof(double);
    }
    return ret;
}

int Pileup::read_depth(std::ifstream &pileup_file, uint16_t *dp, int N)
{
    pileup_file.read(reinterpret_cast<char*>(dp), sizeof(uint16_t)*N);
    return (sizeof(uint16_t)*N);
}

int Pileup::read_readpair(std::ifstream &pileup_file, readpair &rp)
{
    int ret = 0;
    pileup_file.read(reinterpret_cast<char*>(&(rp.chrnum)), sizeof(int8_t));
    pileup_file.read(reinterpret_cast<char*>(&(rp.selfpos)), sizeof(int32_t));
    pileup_file.read(reinterpret_cast<char*>(&(rp.matepos)), sizeof(int32_t));
    pileup_file.read(reinterpret_cast<char*>(&(rp.matequal)), sizeof(int8_t));
    pileup_file.read(reinterpret_cast<char*>(&(rp.pairstr)), sizeof(int8_t));
    ret += sizeof(int8_t) + sizeof(int32_t) + sizeof(int32_t) + sizeof(int8_t) + sizeof(int8_t);
    return ret;
}

int Pileup::read_splitread(std::ifstream &pileup_file, splitread& sp)
{
    int ret = 0;
    
    pileup_file.read(reinterpret_cast<char*>(&(sp.chrnum)), sizeof(int8_t));
    pileup_file.read(reinterpret_cast<char*>(&(sp.pos)), sizeof(int32_t));
    pileup_file.read(reinterpret_cast<char*>(&(sp.sapos)), sizeof(int32_t));
    pileup_file.read(reinterpret_cast<char*>(&(sp.firstclip)), sizeof(int16_t));
    pileup_file.read(reinterpret_cast<char*>(&(sp.secondclip)), sizeof(int16_t));
    ret += sizeof(int8_t) + sizeof(int32_t) + sizeof(int32_t) + sizeof(int16_t) + sizeof(int16_t);
    
    return ret;
}


double Pileup::gcCorrected(double D, int chr, int pos)
{
    int p = pos*2 / gc.binsize;
    //    std::cerr << "pos " << pos << " p " << p << std::endl;
    int bin = gc.gc_array[chr][p];
    
    if (bin<20 && gc_factor[bin]>0.0001)
    {
        //        std::cerr << D << " at " << chr << ":" << pos << " is adjusted to " << D/gc_factor[bin] << " by gc Factor "<< gc_factor[bin] << std::endl;
        return D / gc_factor[bin];
    }
    else
    {
        // Let's not make adjustment for bins with '255' value
        return D;
    }
}

