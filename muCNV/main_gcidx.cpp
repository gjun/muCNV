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

// GC_Content File Structure
// ----------
// MAGIC
// uint8_t : Number of Chr
// (# Chr)
//  | uint32_t : Size of each Chr.
//  | uint32_t : Number of bins in each chr, ceil( (chrssize + 1.0) / bin_dist )
// uint16_t : bin width (how much bp used to average GC content), 400bp
// uint16_t : bin_dist (how much distance in bp between recorded GC contents, 100bp
// uint16_t : number of bins in GC content curve, default: 100 (0 means GC content from 0 to 1%)
// MAGIC
// (# Chr)
//  | (Chr Size) * uint8_t : GC content for each genomic position (every bin_dist-th bp)
//  | MAGIC
// double * (number of bins) : fraction of genomic bin_dist intervals in each GC-bin
// MAGIC

int main_gcidx(int argc, char** argv)
{
    fasta F;
    std::string fasta_file;
    std::string fai_file;
    std::string GC_file;
    
    try
    {
        TCLAP::CmdLine cmd("Command description message", ' ', "0.06");
        
        TCLAP::ValueArg<std::string> argFasta("f","fasta","Genome reference FASTA file",false,"/data/ref/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa","string");
        TCLAP::ValueArg<std::string> argGCout("o","out", "Output GC-content filename", false, "GRCh38.gc", "string");
        TCLAP::SwitchArg switchPrint("p","print", "Print out SV variants", cmd, false);
        
        cmd.add(argFasta);
        cmd.add(argGCout);
        
        cmd.parse(argc, argv);
        
        fasta_file = argFasta.getValue();
        GC_file = argGCout.getValue();
    }
    catch (TCLAP::ArgException &e)
    {
        std::cerr << "Error: " << e.error() << " for arg " << e.argId() << std::endl;
        abort();
    }
    
    //    fastaFileName = getenv("HOME");
    fai_file= fasta_file + ".fai";
    
    F.load(fasta_file.c_str(), fai_file.c_str());
    //    F.printIndex();
    
    std::ofstream outFile(GC_file.c_str(), std::ios::out | std::ios::binary);
    if (!outFile.is_open())
    {
        std::cerr << "Error: Cannot open " << GC_file << " for output."  << std::endl;
        exit(0);
    }
    writemagic(outFile);
    
    uint8_t n_chr = 24;
    
    uint16_t interval_width = 400; // Averaged over 400bp
    uint16_t interval_dist = 100; // One bin for every 100-bp point, approx. 30M bins
    
    std::vector<uint32_t> n_interval;
    
    n_interval.resize(n_chr+1);
    n_interval[0] = 0;
    
    outFile.write(reinterpret_cast <char*> (&n_chr), sizeof(uint8_t));
    for(int i=1;i<=n_chr;++i)
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
        n_interval[i] = ceil(((double)L+1.0) / (double)interval_dist); // Number of intervals in a chromosome
        outFile.write(reinterpret_cast <char*> (&n_interval[i]), sizeof(uint32_t));
    }

    // size of GC bins
    outFile.write(reinterpret_cast <char *> (&interval_width), sizeof(uint16_t));
    outFile.write(reinterpret_cast <char *> (&interval_dist), sizeof(uint16_t));

    uint16_t num_bin = 100;
    // number of GC bins, changed to 100 from 20 on Mar 20, 2019, for GC contents 0-1% bin to 99-100% bin
    outFile.write(reinterpret_cast <char *> (&num_bin), sizeof(uint16_t));
    
    // This part has been commented out on Mar 20, 2019
//
//    uint16_t total_bin = 1800;
//    // number of intervals per GC bin, for sampling
//    outFile.write(reinterpret_cast <char *> (&total_bin), sizeof(uint16_t));
//
//    writemagic(outFile);
//    int cnt = 0;
//    std::ifstream inFile("sorted.selected.txt", std::ios::in);
//    while(inFile.good())
//    {
//        int chr = 0;
//        uint32_t pos = 0;
//        double val = 0;
//
//        if (inFile>>chr)
//        {
//            //        cerr << chr << "\t";
//            uint8_t c = (uint8_t)chr;
//            outFile.write(reinterpret_cast <char *> (&c), sizeof(uint8_t));
//
//            inFile >> pos;
//
//            pos -= 200;
//            outFile.write(reinterpret_cast <char *> (&pos), sizeof(uint32_t));
//            //        cerr << pos << "-";
//
//            pos += 399;
//            outFile.write(reinterpret_cast <char *> (&pos), sizeof(uint32_t));
//            //        cerr << pos << "\t";
//
//            inFile >> val;
//            uint8_t gcbin = (uint8_t)floor(val*20);
//            outFile.write(reinterpret_cast <char *> (&gcbin), sizeof(uint8_t));
//            //        cerr <<val << endl;
//            cnt++;
//        }
//    }
//    inFile.close();
    
    //    cerr << "Current position: " << outFile.tellp() << endl;
    
    writemagic(outFile);
    
    double gc_dist[num_bin];
    for(int j=0;j<num_bin;++j)
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
        double g_buf[4] = {-1.0};
        int cnt = 0;

        uint8_t* gc_array = (uint8_t *) calloc(n_interval[i], sizeof(uint8_t));
        
        std::cerr << "Current position: " << outFile.tellp() << std::endl;
        
        //        for(size_t pos=1; pos<F.chrlen(chr); pos+=200)
        std::vector<double> gc_array_raw (n_interval[i], 0);
        
        for(int j=0, pos=1; j<n_interval[i]; ++j, pos+=interval_dist)
        {
            std::string S;
            F.read(interval_dist, S);
            
            gc_array_raw[j] = gc(S);
        }
        
        for(int j=2; j<n_interval[i]-2; ++j)
        {
            gc_array[j] =0;
            double sum_gc = 0;
            int n_gc = 0;
            for(int k=j-2; k<j+2; ++k)
            {
                if (gc_array_raw[k]>=0)
                {
                    sum_gc += gc_array_raw[k];
                    n_gc ++;
                }
            }
            if (n_gc>0)
            {
                gc_array[j] = round(sum_gc/(double)n_gc);
            }
            else
            {
                gc_array[j] = 255;
            }
            if (gc_array[j] != 255)
                gc_dist[gc_array[j]]++;
        }
        gc_array[0] = gc_array[1] = gc_array[2];
        gc_array[n_interval[i]-1] = gc_array[n_interval[i]-2] = gc_array[n_interval[i]-3];
        
        outFile.write(reinterpret_cast <char*>(gc_array), sizeof(uint8_t) * n_interval[i]);
        writemagic(outFile);
        std::cerr << n_interval[i] << " written." << std::endl;
        
        free(gc_array);
    }
    
    std::cerr << "Current position: " << outFile.tellp() << std::endl;
    
    double sum = 0;
    for (int i=0;i<num_bin;++i)
    {
        sum += gc_dist[i];
    }
    for (int i=0;i<num_bin;++i)
    {
        gc_dist[i] /= sum;
        std::cerr << "bin " << i << " gc dist " << gc_dist[i] << std::endl;
    }
    std::cerr << "Current position: " << outFile.tellp() << std::endl;
    
    outFile.write(reinterpret_cast <char *>(gc_dist), sizeof(double)*num_bin);
    writemagic(outFile);
    
    outFile.close();
    return 0;
}
