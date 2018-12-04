//
//  pileup_main.cpp
//  muCNV
//
//  Created by Goo Jun on 11/22/18.
//  Copyright © 2018 Goo Jun. All rights reserved.
//

#include <stdio.h>
#include <cmath>
// TCLAP headers
#include "tclap/CmdLine.h"
#include "tclap/Arg.h"

#include "in_vcf.h"

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
        
        TCLAP::ValueArg<std::string> argBam("b","bam","Input BAM/CRAM file name",true,"","std::string");
        TCLAP::ValueArg<std::string> argVcf("v","vcf","VCF file containing candidate SVs",false,"","std::string");
        TCLAP::ValueArg<std::string> argInterval("V","interVal", "Binary interval file containing candidate SVs", false, "", "std::string");
        TCLAP::ValueArg<std::string> argSampleID("s","sample","Sample ID for output filename base",false,"","std::string");
        TCLAP::ValueArg<std::string> argGcfile("f","gcFile","File containing GC content information",false, "GRCh38.gc", "std::string");
        
        cmd.add(argBam);
        cmd.add(argVcf);
        cmd.add(argInterval);
        cmd.add(argGcfile);
        cmd.add(argSampleID);
        cmd.parse(argc, argv);
        
        bam_file = argBam.getValue();
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
    
    // TODO: check whether gc is initialized correctly

    
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

    // Read discovery variant list
    std::cerr<< vec_sv.size() << " SVs and " << vec_bp.size() << " breakpoints identified from the VCF/Interval file." << std::endl;
    std::vector<int> idxs;
    std::sort(vec_bp.begin(), vec_bp.end());
    
    std::string pileup_name = sample_id + ".pileup";
    std::string varfile_name = sample_id + ".var";
    std::string idxfile_name = sample_id + ".idx";
    
    // Initialize Pileup
    Pileup pup;
    pup.open(pileup_name, std::ios::out | std::ios::binary);
    std::cerr << "Pileup file initialized" << std::endl;

    //Initialize BAM/CRAM
    BamCram b;
    b.initialize_sequential(bam_file, gc);
    std::cerr << "BAM/CRAM file initialized" << std::endl;
    // Read BAM/CRAM
    b.read_depth_sequential(pup, gc, vec_bp, vec_sv);
    b.postprocess_depth(vec_sv);

    // Write Pileup
    size_t curr_pos = 0;

    BaseFile idx_file;
    idx_file.open(idxfile_name, std::ios::out | std::ios::binary);
    
    curr_pos += pup.write_int32(1);
    curr_pos += pup.write_sample_id(sample_id);
    curr_pos += pup.write_sample_stat(b.stat);
    curr_pos += pup.write_gc_factor(b.gc_factor, gc.num_bin);
    
    idx_file.write_uint64(curr_pos);

    // Write DP100
    for(int i=1; i<=gc.num_chr; ++i)
    {
        curr_pos += pup.write_depth(b.depth100[i], b.nbin_100[i]);
    }
    
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
            idx_file.write_uint64(curr_pos); // where each 10,000-bp interval starts;
            // RP
            while(rp_idx < (int)b.vec_rp.size() && b.vec_rp[rp_idx].chrnum == i && b.vec_rp[rp_idx].selfpos <= j*10000)
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
            curr_pos += pup.write_uint32(n_rp);
            
            for(int k=prev_rp; k<rp_idx; ++k)
            {
                curr_pos += pup.write_readpair(b.vec_rp[k]);
                cnt_rp ++;
            }
            prev_rp = rp_idx;
            
            // SP
            while(sp_idx < (int)b.vec_sp.size() && b.vec_sp[sp_idx].chrnum == i && b.vec_sp[sp_idx].pos <= j*10000)
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
            curr_pos += pup.write_uint32(n_sp);
            
            for(int k=prev_sp; k<sp_idx; ++k)
            {
                curr_pos += pup.write_splitread(b.vec_sp[k]);
                cnt_sp ++;
            }
            prev_sp = sp_idx;
        }
        // fprintf(stderr, "\rChr %d, Index %d, cnt_rp %d, cnt_sp %d\n", i, (int)curr_pos, cnt_rp, cnt_sp);
    }
    pup.close();
    idx_file.close();
    
    Pileup var_file;
    var_file.open(varfile_name, std::ios::out | std::ios::binary);
    
    // Number of samples in this var_file (can include multiple samples)
    var_file.write_int32(1);
    
    // Number of SVs in this var_file (can include multiple samples)
    int n_var = (int)vec_sv.size();
    var_file.write_int32(n_var);
    var_file.write_sample_id(sample_id);
    
    for(int i=0;i<(int)vec_sv.size();++i)
    {
        var_file.write_depth(&(vec_sv[i].dp), 1);
    }
    var_file.close();
    
    std::cerr << "Finished without an error" << std::endl;
    return 0;
}
