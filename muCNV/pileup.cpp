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

void Pileup::write(std::string &sampID, std::vector<sv> &vec_sv)
{
    std::string pileup_name = sampID + ".pileup";
    std::string varfile_name = sampID + ".var";
    std::string idxfile_name = sampID + ".idx";
    
    size_t curr_pos = 0;
    
    std::ofstream pileupFile(pileup_name.c_str(), std::ios::out | std::ios::binary);
    std::ofstream idxFile(idxfile_name.c_str(), std::ios::out | std::ios::binary);
    
    int n_sample = 1;
    char pad[256] = {0};
    
    pileupFile.write(reinterpret_cast<char*>(&n_sample), sizeof(int));
    curr_pos += sizeof(int);
    
    // Sample ID (each with 256 bytes)
    if (sampID.length() > 255)
    {
        std::cerr << "Error, sample ID " << sampID << " is too long." << std::endl;
        exit(1);
    }
    pileupFile.write(sampID.c_str(), sampID.length());
    pileupFile.write(pad, 256-sampID.length());
    curr_pos += 256;
    
    // Write depth and isize stats
    pileupFile.write(reinterpret_cast<char*>(&avg_dp), sizeof(double));
    pileupFile.write(reinterpret_cast<char*>(&std_dp), sizeof(double));
    pileupFile.write(reinterpret_cast<char*>(&avg_isize), sizeof(double));
    pileupFile.write(reinterpret_cast<char*>(&std_isize), sizeof(double));
    
    curr_pos += sizeof(double) * 4;
    
    // Write GC curve
    for(int i=0;i<gc.num_bin;++i)
    {
        pileupFile.write(reinterpret_cast<char*>(&(gc_factor[i])), sizeof(double));
        curr_pos += sizeof(double);
    }
    
    // Write Index of var files (every chr offset, 1000-th variants)
    //    std::cerr << "Sample " << sampID << ", header length " << curr_pos << std::endl;
    
    idxFile.write(reinterpret_cast<char*>(&curr_pos), sizeof(size_t)); // where SV DP starts
    
    // Write DP100
    for(int i=1; i<=gc.num_chr; ++i)
    {
        pileupFile.write(reinterpret_cast<char*>(depth100[i]), sizeof(uint16_t)*(nbin_100[i]));
        curr_pos += sizeof(uint16_t)*nbin_100[i];
    }
    
    //   std::cerr << "After DP100 written, curr_pos is at " << curr_pos << std::endl;
    
    int sp_idx = 0;
    int rp_idx = 0;
    int prev_sp = 0;
    int prev_rp = 0;
    
    int cnt_rp = 0;
    int cnt_sp = 0;
    
    for(int i=1;i<=gc.num_chr; ++i)
    {
        int N = ceil((double)gc.chr_size[i] / 10000.0) ;
        
        for(int j=1;j<=N;++j)
        {
            idxFile.write(reinterpret_cast<char*>(&curr_pos), sizeof(size_t)); // where each 10,000-bp interval starts;
            // RP
            while(rp_idx < (int)vec_rp.size() && vec_rp[rp_idx].chrnum == i && vec_rp[rp_idx].selfpos <= j*10000)
            {
                rp_idx ++;
            }
            uint32_t  n_rp = 0;
            if (rp_idx - prev_rp < 0)
            {
                std::cerr << "Wrong read pair index while writing... " << std::endl;
                exit(1);
            }
            else
            {
                n_rp = (uint32_t) rp_idx - prev_rp;
            }
            pileupFile.write(reinterpret_cast<char*>(&n_rp), sizeof(uint32_t));
            curr_pos += sizeof(uint32_t);
            
            for(int k=prev_rp; k<rp_idx; ++k)
            {
                pileupFile.write(reinterpret_cast<char*>(&(vec_rp[k].chrnum)), sizeof(int8_t));
                pileupFile.write(reinterpret_cast<char*>(&(vec_rp[k].selfpos)), sizeof(int32_t));
                pileupFile.write(reinterpret_cast<char*>(&(vec_rp[k].matepos)), sizeof(int32_t));
                pileupFile.write(reinterpret_cast<char*>(&(vec_rp[k].matequal)), sizeof(int8_t));
                pileupFile.write(reinterpret_cast<char*>(&(vec_rp[k].pairstr)), sizeof(int8_t));
                curr_pos += sizeof(int8_t) + sizeof(int32_t) + sizeof(int32_t) + sizeof(int8_t) + sizeof(int8_t);
                cnt_rp ++;
            }
            prev_rp = rp_idx;
            
            
            // SP
            while(sp_idx < (int)vec_sp.size() && vec_sp[sp_idx].chrnum == i && vec_sp[sp_idx].pos <= j*10000)
            {
                sp_idx ++;
            }
            uint32_t n_sp = 0;
            if (sp_idx - prev_sp < 0)
            {
                std::cerr << "Wrong split read index while writing... " << std::endl;
                exit(1);
            }
            else
            {
                n_sp = (uint32_t) sp_idx - prev_sp;
            }
            pileupFile.write(reinterpret_cast<char*>(&n_sp), sizeof(uint32_t));
            curr_pos += sizeof(uint32_t);
            
            for(int k=prev_sp; k<sp_idx; ++k)
            {
                pileupFile.write(reinterpret_cast<char*>(&(vec_sp[k].chrnum)), sizeof(int8_t));
                pileupFile.write(reinterpret_cast<char*>(&(vec_sp[k].pos)), sizeof(int32_t));
                pileupFile.write(reinterpret_cast<char*>(&(vec_sp[k].sapos)), sizeof(int32_t));
                pileupFile.write(reinterpret_cast<char*>(&(vec_sp[k].firstclip)), sizeof(int16_t));
                pileupFile.write(reinterpret_cast<char*>(&(vec_sp[k].secondclip)), sizeof(int16_t));
                curr_pos += sizeof(int8_t) + sizeof(int32_t) + sizeof(int32_t) + sizeof(int16_t) + sizeof(int16_t);
                cnt_sp ++;
            }
            prev_sp = sp_idx;
        }
        // fprintf(stderr, "\rChr %d, Index %d, cnt_rp %d, cnt_sp %d\n", i, (int)curr_pos, cnt_rp, cnt_sp);
    }
    pileupFile.close();
    idxFile.close();
    
    std::ofstream varFile(varfile_name.c_str(), std::ios::out | std::ios::binary);
    
    // Number of samples in this varFile (can include multiple samples)
    varFile.write(reinterpret_cast<char*>(&n_sample), sizeof(int));
    
    // Number of SVs in this varFile (can include multiple samples)
    int n_var = (int)vec_sv.size();
    varFile.write(reinterpret_cast<char*>(&n_var), sizeof(int));
    
    // Sample ID (each with 256 bytes)
    if (sampID.length() > 255)
    {
        std::cerr << "Error, sample ID " << sampID << " is too long." << std::endl;
        exit(1);
    }
    varFile.write(sampID.c_str(), sampID.length());
    varFile.write(pad, 256-sampID.length());
    
    for(int i=0;i<(int)vec_sv.size();++i)
    {
        varFile.write(reinterpret_cast<char*>(&(vec_sv[i].dp)), sizeof(uint16_t));
    }
    varFile.close();
}



void Pileup::write_text(std::string &sampID, std::vector<sv> &vec_sv)
{
    // TODO: add write GC curve
    // TODO: write average 100-bp depth
    std::string fname = sampID + ".pileup.txt";
    FILE *fp = fopen(fname.c_str(), "w");
    
    fclose(fp);
    
    fname = sampID + ".var.txt";
    fp = fopen(fname.c_str(), "w");
    for(int i=0;i<(int)vec_sv.size();++i)
    {
        /*
         uint16_t n_rp = (uint16_t) vec_sv[i].vec_pair.size();
         uint16_t n_sp = (uint16_t) vec_sv[i].vec_split.size();
         */
        fprintf(fp, "%d", vec_sv[i].dp);
        /*
         fprintf(fp, "\t%d\t%d\t", n_rp, n_sp);
         
         for(int j=0;j<vec_sv[i].vec_pair.size();++j)
         {
         fprintf(fp, "%d,%d,%d,%d;", vec_sv[i].vec_pair[j].selfpos, vec_sv[i].vec_pair[j].matepos, vec_sv[i].vec_pair[j].selfstr, vec_sv[i].vec_pair[j].matestr);
         }
         fprintf(fp,"\t");
         for(int j=0;j<vec_sv[i].vec_split.size();++j)
         {
         fprintf(fp, "%d,%d,%d,%d;", vec_sv[i].vec_split[j].pos, vec_sv[i].vec_split[j].sapos, vec_sv[i].vec_split[j].firstclip, vec_sv[i].vec_split[j].secondclip);
         }
         */
        fprintf(fp, "\n");
        
    }
    fclose(fp);
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

