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

bool sclip::operator < (const sclip& b) const
{
    return (chrnum<b.chrnum || (chrnum == b.chrnum && pos<b.pos));
}

bool sclip::operator == (const sclip& b) const
{
    return (chrnum == b.chrnum && pos == b.pos);
}

bool sclip::operator <= (const sclip& b) const
{
    return (chrnum<b.chrnum || (chrnum==b.chrnum && pos<=b.pos));
    //    return (*this <b || *this ==b);
}

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
    //fs.write(reinterpret_cast<char*>(&(rp.chrnum)), sizeof(int8_t));
    int16_t selfpos = (int16_t) (rp.selfpos % 10000);
    int16_t matepos = (int16_t) (rp.matepos % 10000);
    
    fs.write(reinterpret_cast<char*>(&(selfpos)), sizeof(int16_t));
    fs.write(reinterpret_cast<char*>(&(matepos)), sizeof(int16_t));
    fs.write(reinterpret_cast<char*>(&(rp.matequal)), sizeof(int8_t));
    fs.write(reinterpret_cast<char*>(&(rp.pairstr)), sizeof(int8_t));
//    ret += sizeof(int8_t) + sizeof(int32_t) + sizeof(int32_t) + sizeof(int8_t) + sizeof(int8_t);
 //   ret += sizeof(int32_t) + sizeof(int32_t) + sizeof(int8_t) + sizeof(int8_t);
    ret += sizeof(int16_t) + sizeof(int16_t) + sizeof(int8_t) + sizeof(int8_t);

    return ret;
}

int Pileup::write_splitread(splitread& sp)
{
    int ret = 0;
    int16_t pos = (int16_t) (sp.pos % 10000);
    int16_t sapos = (int16_t) (sp.sapos % 10000);
    
 //   fs.write(reinterpret_cast<char*>(&(sp.chrnum)), sizeof(int8_t));
    fs.write(reinterpret_cast<char*>(&(pos)), sizeof(int16_t));
    fs.write(reinterpret_cast<char*>(&(sapos)), sizeof(int16_t));
    
    fs.write(reinterpret_cast<char*>(&(sp.firstclip)), sizeof(int16_t));
    fs.write(reinterpret_cast<char*>(&(sp.secondclip)), sizeof(int16_t));
//    ret += sizeof(int8_t) + sizeof(int32_t) + sizeof(int32_t) + sizeof(int16_t) + sizeof(int16_t);
    ret += sizeof(int16_t) + sizeof(int16_t) + sizeof(int16_t) + sizeof(int16_t);

    return ret;
}

int Pileup::write_softclip(sclip &sc)
{
    int ret = 0;
    int16_t pos16 = (int16_t)(sc.pos % 10000); // from the offset

 //   fs.write(reinterpret_cast<char*>(&(sc.chrnum)), sizeof(int8_t));
 //   fs.write(reinterpret_cast<char*>(&(sc.pos)), sizeof(int32_t)); // The position where M and S divides left/right will be recoreded in separate lists
    fs.write(reinterpret_cast<char*>(&pos16), sizeof(int16_t));
 //   ret += sizeof(int8_t) + sizeof(int32_t);
 //   ret += sizeof(int32_t);
    ret += sizeof(int16_t);

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
    int16_t selfpos=0;
    int16_t matepos=0;
    
//    fs.read(reinterpret_cast<char*>(&(rp.chrnum)), sizeof(int8_t));
    fs.read(reinterpret_cast<char*>(&(selfpos)), sizeof(int16_t));
    fs.read(reinterpret_cast<char*>(&(matepos)), sizeof(int16_t));
    fs.read(reinterpret_cast<char*>(&(rp.matequal)), sizeof(int8_t));
    fs.read(reinterpret_cast<char*>(&(rp.pairstr)), sizeof(int8_t));
//    ret += sizeof(int8_t) + sizeof(int32_t) + sizeof(int32_t) + sizeof(int8_t) + sizeof(int8_t);
    ret += sizeof(int16_t) + sizeof(int16_t) + sizeof(int8_t) + sizeof(int8_t);

    rp.selfpos = selfpos;
    rp.matepos = matepos;
    
    return ret;
}

int Pileup::read_splitread(splitread& sp)
{
    int ret = 0;
    int16_t pos;
    int16_t sapos;
    
//    fs.read(reinterpret_cast<char*>(&(sp.chrnum)), sizeof(int8_t));
    fs.read(reinterpret_cast<char*>(&(pos)), sizeof(int16_t));
    fs.read(reinterpret_cast<char*>(&(sapos)), sizeof(int16_t));
    fs.read(reinterpret_cast<char*>(&(sp.firstclip)), sizeof(int16_t));
    fs.read(reinterpret_cast<char*>(&(sp.secondclip)), sizeof(int16_t));
//    ret += sizeof(int8_t) + sizeof(int32_t) + sizeof(int32_t) + sizeof(int16_t) + sizeof(int16_t);
    ret += sizeof(int16_t) + sizeof(int16_t) + sizeof(int16_t) + sizeof(int16_t);
    
    sp.pos = pos;
    sp.sapos = sapos;

    return ret;
}

int Pileup::read_softclip(sclip &sc)
{
    int ret = 0;
    int16_t pos16;
    
//    fs.read(reinterpret_cast<char*>(&(sc.chrnum)), sizeof(int8_t));
//    fs.read(reinterpret_cast<char*>(&(sc.pos)), sizeof(int32_t)); // The position where M and S divides, and also direction by the sign
    fs.read(reinterpret_cast<char*>(&(pos16)), sizeof(int16_t)); // The position where M and S divides, and also direction by the sign

 //   ret += sizeof(int8_t) + sizeof(int32_t);
  //  ret += sizeof(int32_t);
    ret += sizeof(int16_t);
    sc.pos = (int32_t) pos16;
    
    return ret;
}

int Pileup::fix_offset_pos(int32_t j, int32_t r)
{
    int ans = 0;
    
    if (j<30000)
    {
        ans = r;
    }
    else
    {
        if ( ((j+10000) + 32768) / 65536 > (j + 32768) / 65536)
        {
            int32_t kk = (int32_t)((j+10000+32768)/65536)*65536 - 32768 - j;
            
            if (r>=0 && r<2768)
            {
                ans = kk - (2768-r);
            }
            else if (r>2767)
            {
                ans = kk - (2768-(r-10000));
            }
            else if (r<0 && r>= -2768)
            {
                ans = kk + (r+2768);
            }
            else
            {
                ans = kk + 10000 + r + 2768;
            }
        }
        else
        {
            int16_t mm = (int16_t) j;
            int16_t rr = mm % 10000;
            if (r >= rr)
            {
                ans = r-rr;
            }
            else
            {
                ans = r-rr+10000;
            }
        }
    }
    return ans;
}
