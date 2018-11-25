//
//  pileup_main.cpp
//  muCNV
//
//  Created by Goo Jun on 11/22/18.
//  Copyright © 2018 Goo Jun. All rights reserved.
//

#include <stdio.h>
// TCLAP headers
#include "tclap/CmdLine.h"
#include "tclap/Arg.h"

#include "muCNV.h"

#include "gc_content.h"
#include "pileup.h"
#include "bam_cram.h"
#include "sv.h"

int main_pileup(int argc, char** argv)
{
    // Parsing command-line arguments
    std::string bam_file;
    std::string out_filename;
    std::string vcf_file;
    std::string interval_file;
    std::string gc_file;
    std::string sample_id;

    int n_sample = 0;

    std::vector<std::string> vcfs;
    
    try
    {
        TCLAP::CmdLine cmd("Command description message", ' ', "0.06");
        
        TCLAP::ValueArg<std::string> argBam("b","bam","Input BAM/CRAM file name",false,"","std::string");
        TCLAP::ValueArg<std::string> argOut("o","out","Output filename",false,"muCNV.vcf","std::string");
        TCLAP::ValueArg<std::string> argVcf("v","vcf","VCF file containing candidate SVs",false,"","std::string");
        TCLAP::ValueArg<std::string> argInterval("V","interVal", "Binary interval file containing candidate SVs", false, "", "std::string");
        TCLAP::ValueArg<std::string> argSampleID("s","sample","Sample ID",false,"","std::string");
        TCLAP::ValueArg<std::string> argGcfile("f","gcFile","File containing GC content information",false, "GRCh38.gc", "std::string");
        
        cmd.add(argBam);
        cmd.add(argOut);
        cmd.add(argVcf);
        cmd.add(argInterval);
        cmd.add(argGcfile);
        cmd.add(argSampleID);
        cmd.parse(argc, argv);
        
        bam_file = argBam.getValue();
        out_filename = argOut.getValue();
        sample_id = argSampleID.getValue();
        vcf_file = argVcf.getValue();
        interval_file = argInterval.getValue();
        gc_file = argGcfile.getValue();
    }
    catch (TCLAP::ArgException &e)
    {
        std::cerr << "Error: " << e.error() << " for arg " << e.argId() << std::endl;
        abort();
    }
    
    // Generate summary stats from BAM/CRAM
    n_sample = 1;
    std::cerr << "Processing individual BAM/CRAM file to genearte pileup information" << std::endl;
    
    GcContent gc;
    gc.initialize(gc_file);
    std::cerr << "GC content initialized" << std::endl;
    
    Pileup pup (gc);
    std::cerr << "Pileup file initialized" << std::endl;
    
    std::vector<sv> vec_sv;
    std::vector<breakpoint> vec_bp;
    vcfs.push_back(vcf_file);
    
    if (vcf_file != "")
    {
        read_svs_from_vcf(vcf_file, vec_bp, vec_sv);
    }
    else if (interval_file != "")
    {
        read_svs_from_intfile(interval_file, vec_bp, vec_sv);
    }
    else
    {
        std::cerr << "Error, VCF or Interval file is required." << std::endl;
    }

    std::cerr<< vec_sv.size() << " svs and " << vec_bp.size() << " breakpoints identified from the VCF/Interval file." << std::endl;
    std::vector<int> idxs;
    
    std::sort(vec_bp.begin(), vec_bp.end());
    
    BamCram b(pup);
    
    b.initialize_sequential(bam_file);
    std::cerr << "BAM/CRAM file initialized" << std::endl;
    
    b.read_depth_sequential(vec_bp, vec_sv);
   
    b.postprocess_depth(vec_sv);
    
    std::string pileup_name = sample_id + ".pileup";
    std::string var_file_name = sample_id + ".var";
    std::string idxfile_name = sample_id + ".idx";
    
    size_t curr_pos = 0;
    
    std::ofstream pileup_file(pileup_name.c_str(), std::ios::out | std::ios::binary);
    std::ofstream idxFile(idxfile_name.c_str(), std::ios::out | std::ios::binary);
    
    curr_pos += pup.write_number(pileup_file, 1);
    curr_pos += pup.write_sample_id(pileup_file, sample_id);
    curr_pos += pup.write_sample_stat(pileup_file, pup.stat);
    curr_pos += pup.write_gc_factor(pileup_file, pup.gc_factor);
    
    idxFile.write(reinterpret_cast<char*>(&curr_pos), sizeof(size_t)); // 1-st index: where DP100 starts, should be the same for all samples
    
    // Write DP100
    for(int i=1; i<=gc.num_chr; ++i)
    {
        curr_pos += pup.write_depth(pileup_file, pup.depth100[i], pup.nbin_100[i]);
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
            while(rp_idx < (int)pup.vec_rp.size() && pup.vec_rp[rp_idx].chrnum == i && pup.vec_rp[rp_idx].selfpos <= j*10000)
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
            pileup_file.write(reinterpret_cast<char*>(&n_rp), sizeof(uint32_t));
            curr_pos += sizeof(uint32_t);
            
            for(int k=prev_rp; k<rp_idx; ++k)
            {
                curr_pos += pup.write_readpair(pileup_file, pup.vec_rp[k]);
                cnt_rp ++;
            }
            prev_rp = rp_idx;
            
            // SP
            while(sp_idx < (int)pup.vec_sp.size() && pup.vec_sp[sp_idx].chrnum == i && pup.vec_sp[sp_idx].pos <= j*10000)
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
            pileup_file.write(reinterpret_cast<char*>(&n_sp), sizeof(uint32_t));
            curr_pos += sizeof(uint32_t);
            
            for(int k=prev_sp; k<sp_idx; ++k)
            {
                curr_pos += pup.write_splitread(pileup_file, pup.vec_sp[k]);
                cnt_sp ++;
            }
            prev_sp = sp_idx;
        }
        // fprintf(stderr, "\rChr %d, Index %d, cnt_rp %d, cnt_sp %d\n", i, (int)curr_pos, cnt_rp, cnt_sp);
    }
    pileup_file.close();
    idxFile.close();
    
    std::ofstream var_file(var_file_name.c_str(), std::ios::out | std::ios::binary);
    // Number of samples in this var_file (can include multiple samples)
    pup.write_number(var_file, 1);
    
    // Number of SVs in this var_file (can include multiple samples)
    int n_var = (int)vec_sv.size();
    pup.write_number(var_file, n_var);
    pup.write_sample_id(var_file, sample_id);
    
    for(int i=0;i<(int)vec_sv.size();++i)
    {
        pup.write_depth(var_file, &(vec_sv[i].dp), 1);        
    }
    var_file.close();
    
    
    std::cerr << "Finished without an error" << std::endl;
    return 0;
}
