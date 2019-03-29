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
 
int Pileup::write_sample_stat(SampleStat &s)
{
    // Write depth and isize stats
    fs.write(reinterpret_cast<char*>(&s.avg_dp), sizeof(double));
    fs.write(reinterpret_cast<char*>(&s.std_dp), sizeof(double));
    fs.write(reinterpret_cast<char*>(&s.avg_isize), sizeof(double));
    fs.write(reinterpret_cast<char*>(&s.std_isize), sizeof(double));
    return sizeof(double) * 4;
}

int Pileup::write_gc_factor(std::vector<double>& gc_f, int num_bin)
{
    int ret=0;
    // Write GC curve
    for(int i=0;i<num_bin;++i)
    {
        fs.write(reinterpret_cast<char*>(&(gc_f[i])), sizeof(double));
        ret += sizeof(double);
    }
    return ret;
}

int Pileup::write_depth(uint16_t *dp, int N)
{
    fs.write(reinterpret_cast<char*>(dp), sizeof(uint16_t)*N);
    return (sizeof(uint16_t)*N);
}

int Pileup::write_readpair(readpair &rp)
{
    int ret = 0;
    fs.write(reinterpret_cast<char*>(&(rp.chrnum)), sizeof(int8_t));
    fs.write(reinterpret_cast<char*>(&(rp.selfpos)), sizeof(int32_t));
    fs.write(reinterpret_cast<char*>(&(rp.matepos)), sizeof(int32_t));
    fs.write(reinterpret_cast<char*>(&(rp.matequal)), sizeof(int8_t));
    fs.write(reinterpret_cast<char*>(&(rp.pairstr)), sizeof(int8_t));
    ret += sizeof(int8_t) + sizeof(int32_t) + sizeof(int32_t) + sizeof(int8_t) + sizeof(int8_t);
    return ret;
}

int Pileup::write_splitread(splitread& sp)
{
    int ret = 0;
    
    fs.write(reinterpret_cast<char*>(&(sp.chrnum)), sizeof(int8_t));
    fs.write(reinterpret_cast<char*>(&(sp.pos)), sizeof(int32_t));
    fs.write(reinterpret_cast<char*>(&(sp.sapos)), sizeof(int32_t));
    
    // softclips are separate entries now, 03/27/19
    // fs.write(reinterpret_cast<char*>(&(sp.firstclip)), sizeof(int16_t));
    // fs.write(reinterpret_cast<char*>(&(sp.secondclip)), sizeof(int16_t));
    // ret += sizeof(int8_t) + sizeof(int32_t) + sizeof(int32_t) + sizeof(int16_t) + sizeof(int16_t);
    ret += sizeof(int8_t) + sizeof(int32_t) + sizeof(int32_t);
    
    return ret;
}

int Pileup::write_softclip(sclip &sc)
{
    int ret = 0;
    
    fs.write(reinterpret_cast<char*>(&(sc.chrnum)), sizeof(int8_t));
    fs.write(reinterpret_cast<char*>(&(sc.pos)), sizeof(int32_t)); // The position where M and S divides left/right will be recoreded in separate lists
    ret += sizeof(int8_t) + sizeof(int32_t);
    
    return ret;
}

int Pileup::read_sample_stat(SampleStat &s)
{
    // read depth and isize stats
    fs.read(reinterpret_cast<char*>(&s.avg_dp), sizeof(double));
    fs.read(reinterpret_cast<char*>(&s.std_dp), sizeof(double));
    fs.read(reinterpret_cast<char*>(&s.avg_isize), sizeof(double));
    fs.read(reinterpret_cast<char*>(&s.std_isize), sizeof(double));
    return sizeof(double) * 4;
}

int Pileup::read_gc_factor(std::vector<double>& gc_f, int num_bin)
{
    int ret=0;
    // read GC curve
    gc_f.resize(num_bin);
    
    for(int i=0;i<num_bin;++i)
    {
        fs.read(reinterpret_cast<char*>(&(gc_f[i])), sizeof(double));
        ret += sizeof(double);
    }
    return ret;
}

int Pileup::read_depth(uint16_t *dp, int N)
{
    fs.read(reinterpret_cast<char*>(dp), sizeof(uint16_t)*N);
    return (sizeof(uint16_t)*N);
}

int Pileup::read_readpair(readpair &rp)
{
    int ret = 0;
    fs.read(reinterpret_cast<char*>(&(rp.chrnum)), sizeof(int8_t));
    fs.read(reinterpret_cast<char*>(&(rp.selfpos)), sizeof(int32_t));
    fs.read(reinterpret_cast<char*>(&(rp.matepos)), sizeof(int32_t));
    fs.read(reinterpret_cast<char*>(&(rp.matequal)), sizeof(int8_t));
    fs.read(reinterpret_cast<char*>(&(rp.pairstr)), sizeof(int8_t));
    ret += sizeof(int8_t) + sizeof(int32_t) + sizeof(int32_t) + sizeof(int8_t) + sizeof(int8_t);
    return ret;
}

int Pileup::read_splitread(splitread& sp)
{
    int ret = 0;
    
    fs.read(reinterpret_cast<char*>(&(sp.chrnum)), sizeof(int8_t));
    fs.read(reinterpret_cast<char*>(&(sp.pos)), sizeof(int32_t));
    fs.read(reinterpret_cast<char*>(&(sp.sapos)), sizeof(int32_t));
    fs.read(reinterpret_cast<char*>(&(sp.firstclip)), sizeof(int16_t));
    fs.read(reinterpret_cast<char*>(&(sp.secondclip)), sizeof(int16_t));
    ret += sizeof(int8_t) + sizeof(int32_t) + sizeof(int32_t) + sizeof(int16_t) + sizeof(int16_t);
    
    return ret;
}

int Pileup::read_softclip(sclip &sc)
{
    int ret = 0;
    
    fs.read(reinterpret_cast<char*>(&(sc.chrnum)), sizeof(int8_t));
    fs.read(reinterpret_cast<char*>(&(sc.pos)), sizeof(int32_t)); // The position where M and S divides, and also direction by the sign
    ret += sizeof(int8_t) + sizeof(int32_t);
    
    return ret;
}
