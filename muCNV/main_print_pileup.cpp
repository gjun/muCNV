//
//  main_print_pileup.cpp
//  muCNV
//
//  Created by Goo Jun on 11/23/18.
//  Copyright © 2018 Goo Jun. All rights reserved.
//

#include <stdio.h>
// TCLAP headers
#include "tclap/CmdLine.h"
#include "tclap/Arg.h"

#include "muCNV.h"
#include "gc_content.h"

int main_print_pileup(int argc, char** argv)
{
    string index_file;
    string vcf_file;
    string interval_file;
    string gc_file;
    string sampID;
    string region;
    
    std::vector<string> sample_ids;
    
    // Parsing command-line arguments
    try
    {
        TCLAP::CmdLine cmd("Command description message", ' ', "0.06");
        

        TCLAP::ValueArg<string> argSampleID("s","sample","Sample ID",false,"","string");
        TCLAP::ValueArg<string> argVcf("v","vcf","VCF file containing candidate SVs",false,"","string");
        TCLAP::ValueArg<string> argInterval("V","interVal", "Binary interval file containing candidate SVs", false, "", "string");
        TCLAP::ValueArg<string> argGcfile("f","gcFile","File containing GC content information",false, "GRCh38.gc", "string");
        TCLAP::ValueArg<string> argRegion("r", "region", "Genomic region (chr:start-end)", false, "", "string" );
        
        cmd.add(argVcf);
        cmd.add(argInterval);
        cmd.add(argGcfile);
        cmd.add(argSampleID);
        cmd.add(argRegion);
        
        cmd.parse(argc, argv);
        
        sampID = argSampleID.getValue();
        vcf_file = argVcf.getValue();
        interval_file = argInterval.getValue();
        gc_file = argGcfile.getValue();
        region = argRegion.getValue();
       
    }
    catch (TCLAP::ArgException &e)
    {
        std::cerr << "Error: " << e.error() << " for arg " << e.argId() << std::endl;
        abort();
    }
    
    int n_sample = 0;
    
    std::vector<sv> vec_sv;
    std::vector<breakpoint> vec_bp;
    
    string pileup_name = sampID + ".pileup";
    string varfile_name = sampID + ".var";
    string idxfile_name = sampID + ".idx";
    
    // TODO: make this also work with VCF file
    // read out and print pileup info
    read_svs_from_intfile(interval_file, vec_bp, vec_sv);
    
    std::ifstream pileupFile(pileup_name.c_str(), std::ios::in | std::ios::binary);
    std::ifstream varFile(varfile_name.c_str(), std::ios::in | std::ios::binary);
    std::ifstream idxFile(idxfile_name.c_str(), std::ios::in | std::ios::binary);
    
    size_t curr_idx = 0;
    
    char buf[256];
    
    pileupFile.read(reinterpret_cast<char*>(&n_sample), sizeof(int));
    printf("n_sample(pileup) : %d \n", n_sample);
    // TODO: read N-sample IDs
    pileupFile.read(reinterpret_cast<char *>(buf), 256);
    printf("sample ID(pileup) : %s\n", buf);
    
    double avg_dp, std_dp, avg_isize, std_isize;
    pileupFile.read(reinterpret_cast<char*>(&avg_dp), sizeof(double));
    pileupFile.read(reinterpret_cast<char*>(&std_dp), sizeof(double));
    pileupFile.read(reinterpret_cast<char*>(&avg_isize), sizeof(double));
    pileupFile.read(reinterpret_cast<char*>(&std_isize), sizeof(double));
    
    printf("AVG DP: %f, STdev: %f, AVG ISIZE: %f, STdev: %f \n", avg_dp, std_dp, avg_isize, std_isize);
    
    GcContent gc;
    gc.initialize(gc_file);
    
    // TODO: Use Pileup::gc_factor
    
    std::vector<double> gc_factor (gc.num_bin);
    
    for(int i=0;i<gc.num_bin;++i)
    {
        pileupFile.read(reinterpret_cast<char*>(&(gc_factor[i])), sizeof(double));
        printf("GC-bin %d: %f\n", i, gc_factor[i]);
    }
    idxFile.read(reinterpret_cast<char*>(&curr_idx), sizeof(size_t));
    printf("index position %d, tellg position %d\n", (int)curr_idx, (int)pileupFile.tellg());
    
    uint64_t dpsum = 0;
    uint64_t n_dp = 0;
    
    for(int i=1; i<=gc.num_chr; ++i)
    {
        int N = ceil((double)gc.chr_size[i] / 100.0) + 1;
        uint16_t dp100;
        for(int j=0;j<N;++j)
        {
            pileupFile.read(reinterpret_cast<char*>(&dp100), sizeof(uint16_t));
            dpsum += dp100;
            n_dp+=1;
        }
    }
    printf("Average DP100: %d\n", (int)round((double)dpsum/n_dp/32.0));
    
    while(pileupFile.good())
    {
        uint32_t n_rp = 0;
        uint32_t n_sp = 0;
        
        idxFile.read(reinterpret_cast<char*>(&curr_idx), sizeof(size_t));
        printf("index position %d, tellg position %d\n", (int)curr_idx, (int)pileupFile.tellg());
        
        pileupFile.read(reinterpret_cast<char*>(&n_rp), sizeof(uint32_t));
        printf("%d readpairs\n", n_rp);
        for(int k=0; k<n_rp; ++k)
        {
            int8_t chrnum, pairstr;
            uint8_t matequal;
            int32_t selfpos, matepos;
            pileupFile.read(reinterpret_cast<char*>(&(chrnum)), sizeof(int8_t));
            pileupFile.read(reinterpret_cast<char*>(&(selfpos)), sizeof(int32_t));
            pileupFile.read(reinterpret_cast<char*>(&(matepos)), sizeof(int32_t));
            pileupFile.read(reinterpret_cast<char*>(&(matequal)), sizeof(uint8_t));
            pileupFile.read(reinterpret_cast<char*>(&(pairstr)), sizeof(int8_t));
            printf("\t%d\t%d\t%d\t%u\t%d\n", chrnum, selfpos, matepos, matequal, pairstr);
        }
        
        pileupFile.read(reinterpret_cast<char*>(&n_sp), sizeof(uint32_t));
        printf("%d split reads\n", n_sp);
        for(int k=0; k<n_sp; ++k)
        {
            int8_t chrnum;
            int32_t pos, sapos;
            int16_t firstclip, secondclip;
            pileupFile.read(reinterpret_cast<char*>(&(chrnum)), sizeof(int8_t));
            pileupFile.read(reinterpret_cast<char*>(&(pos)), sizeof(int32_t));
            pileupFile.read(reinterpret_cast<char*>(&(sapos)), sizeof(int32_t));
            pileupFile.read(reinterpret_cast<char*>(&(firstclip)), sizeof(int16_t));
            pileupFile.read(reinterpret_cast<char*>(&(secondclip)), sizeof(int16_t));
            printf("\t%d\t%d\t%d\t%d\t%d\n", chrnum, pos, sapos, firstclip,secondclip);
        }
    }
    pileupFile.close();
    
    int n_var = 0;
    varFile.read(reinterpret_cast<char*>(&n_sample), sizeof(int));
    varFile.read(reinterpret_cast<char*>(&n_var), sizeof(int));
    printf("n_var : %d\n", n_var);
    
    varFile.read(reinterpret_cast<char *>(buf), 256);
    printf("sample ID(var) : %s\n", buf);
    
    for(int i=0;i<n_var;++i)
    {
        uint16_t dp;
        vec_sv[i].print();
        varFile.read(reinterpret_cast<char*>(&dp), sizeof(uint16_t));
        printf("var %d: %f\n",i, (dp/32.0));
    }
    varFile.close();
    
    return 0;
}
