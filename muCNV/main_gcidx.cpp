//
//  main_gcidx.cpp
//  muCNV
//
//  Created by Goo Jun on 11/24/18.
//  Copyright © 2018 Goo Jun. All rights reserved.
//

#include <stdio.h>

// TCLAP headers
#include "tclap/CmdLine.h"
#include "tclap/Arg.h"


#include <iostream>
#include <string>
#include <cstdlib>
#include <sstream>
#include <vector>
#include <math.h>
#include "fasta.h"

double gc(std::string &S)
{
    int AT = 0;
    int GC = 0;
    
    for(int j=0;j<(int)S.length();++j)
    {
        if (S[j] == 'A'  || S[j] == 'T')
        {
            AT++;
        }
        else if (S[j] == 'G' || S[j] == 'C')
        {
            GC++;
        }
    }
    if (AT+GC>(int)(S.length()/2))
    {
        return (GC / double(AT+GC));
    }
    else
    {
        return -1;
    }
}

void writemagic(std::ofstream &F)
{
    char magic [] = "mCNVMGC";
    
    F.write(reinterpret_cast<char *>(magic), 7);
}

int main_gcidx(int argc, char** argv)
{
    fasta F;
    std::string fastaFileName;
    std::string faidxFileName;
    std::string gcOutFileName = "GRCh38.gc";
    
    //    fastaFileName = getenv("HOME");
    fastaFileName = "/data/ref/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa";
    faidxFileName = fastaFileName + ".fai";
    
    F.load(fastaFileName.c_str(), faidxFileName.c_str());
    //    F.printIndex();
    
    std::ofstream outFile(gcOutFileName.c_str(), std::ios::out | std::ios::binary);
    if (!outFile.is_open())
    {
        std::cerr << "Error: Cannot open " << gcOutFileName << " for output."  << std::endl;
        exit(0);
    }
    
    writemagic(outFile);
    
    uint8_t n_chr = 24;
    outFile.write(reinterpret_cast <char*> (&n_chr), sizeof(uint8_t));
    for(int i=1;i<25;++i)
    {
        std::ostringstream sStr;
        uint32_t L = 0;
        sStr << "chr" << i;
        if (i<23)
        {
            L = (uint32_t)F.chrlen(sStr.str());
        }
        else if (i==23)
        {
            std::string s = "chrX";
            L = (uint32_t)F.chrlen(s);
        }
        else if (i==24)
        {
            std::string s = "chrY";
            L = (uint32_t)F.chrlen(s);
        }
        outFile.write(reinterpret_cast <char*> (&L), sizeof(uint32_t));
    }
    
    uint16_t binsize = 400;
    // size of GC bins
    outFile.write(reinterpret_cast <char *> (&binsize), sizeof(uint16_t));
    
    uint16_t num_bin = 20;
    // number of GC bins
    outFile.write(reinterpret_cast <char *> (&num_bin), sizeof(uint16_t));
    
    uint16_t total_bin = 1800;
    // number of intervals per GC bin
    outFile.write(reinterpret_cast <char *> (&total_bin), sizeof(uint16_t));
    
    writemagic(outFile);
    int cnt = 0;
    std::ifstream inFile("sorted.selected.txt", std::ios::in);
    while(inFile.good())
    {
        int chr = 0;
        uint32_t pos = 0;
        double val = 0;
        
        if (inFile>>chr)
        {
            //        cerr << chr << "\t";
            uint8_t c = (uint8_t)chr;
            outFile.write(reinterpret_cast <char *> (&c), sizeof(uint8_t));
            
            inFile >> pos;
            
            pos -= 200;
            outFile.write(reinterpret_cast <char *> (&pos), sizeof(uint32_t));
            //        cerr << pos << "-";
            
            pos += 399;
            outFile.write(reinterpret_cast <char *> (&pos), sizeof(uint32_t));
            //        cerr << pos << "\t";
            
            inFile >> val;
            uint8_t gcbin = (uint8_t)floor(val*20);
            outFile.write(reinterpret_cast <char *> (&gcbin), sizeof(uint8_t));
            //        cerr <<val << endl;
            cnt++;
        }
    }
    inFile.close();
    
    //    cerr << "Current position: " << outFile.tellp() << endl;
    
    writemagic(outFile);
    
    double gc_dist[20];
    for(int j=0;j<20;++j)
    {
        gc_dist[j] = 0;
    }
    
    // 8 byte per each chr len ? or 4 byte ?
    for(int i=1;i<25;++i)
    {
        std::string chr = "chr" + std::to_string(i);
        if (i==23)
            chr = "chrX";
        else if (i==24)
            chr = "chrY";
        
        F.seek(chr, 1);
        double prev_g = -1;
        
        int N = ceil(F.chrlen(chr) / 200.0);
        uint8_t* gc_array = (uint8_t *) calloc(N, sizeof(uint8_t));
        std::cerr << "Current position: " << outFile.tellp() << std::endl;
        
        //        for(size_t pos=1; pos<F.chrlen(chr); pos+=200)
        for(int j=0, pos=1;j<N;++j, pos+=200)
        {
            
            std::string S;
            F.read(200, S);
            double g = gc(S);
            //            if (g>=0 && prev_g>=0)
            {
                //                printf("%d\t%d\t%f\n", i, pos, (g+prev_g)/2.0);
            }
            
            uint8_t gcbin = (uint8_t)floor((g+prev_g) * 10.0);
            if (gcbin == 20) gcbin = 19;
            if (g<=0 || prev_g<=0)
            {
                gcbin=255;
            }
            else
            {
                gc_dist[gcbin] += 1;
            }
            gc_array[j] = gcbin;
            //outFile.write(reinterpret_cast <char *>(&gcbin), sizeof(uint8_t));
            prev_g = g;
        }
        outFile.write(reinterpret_cast <char*>(gc_array), sizeof(uint8_t) * N);
        writemagic(outFile);
        std::cerr << N << " written." << std::endl;
        
        free(gc_array);
    }
    
    std::cerr << "Current position: " << outFile.tellp() << std::endl;
    
    double sum = 0;
    for (int i=0;i<20;++i)
    {
        sum += gc_dist[i];
    }
    for (int i=0;i<20;++i)
    {
        gc_dist[i] /= sum;
        std::cerr << "bin " << i << " gc dist " << gc_dist[i] << std::endl;
    }
    std::cerr << "Current position: " << outFile.tellp() << std::endl;
    
    outFile.write(reinterpret_cast <char *>(gc_dist), sizeof(double)*20);
    writemagic(outFile);
    
    outFile.close();
    return 0;
}
