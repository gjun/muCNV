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
    std::string out_dir;
    std::string gc_file;
    std::string sample_id;

    std::vector<std::string> vcfs;
    
    try
    {
        TCLAP::CmdLine cmd("Command description message", ' ', "0.06");
        
        TCLAP::ValueArg<std::string> argBam("b","bam","Input BAM/CRAM file name",true,"","string");
        TCLAP::ValueArg<std::string> argVcf("v","vcf","VCF file containing candidate SVs",false,"","string");
        TCLAP::ValueArg<std::string> argInterval("V","interVal", "Binary interval file containing candidate SVs", false, "", "string");
        TCLAP::ValueArg<std::string> argSampleID("s","sample","Sample ID, also used for output filenames (.pileup, .var, .idx)",true,"","string");
        TCLAP::ValueArg<std::string> argOutDir("o","outdir","Output directory, current directory if omitted",false,".","string");
        TCLAP::ValueArg<std::string> argGcfile("f","gcFile","File containing GC content information (default: GRCh38.gc)",false, "GRCh38.gc", "string");
        
        cmd.add(argBam);
        cmd.add(argVcf);
        cmd.add(argInterval);
        cmd.add(argGcfile);
        cmd.add(argSampleID);
        cmd.add(argOutDir);
        cmd.parse(argc, argv);
        
        bam_file = argBam.getValue();
        sample_id = argSampleID.getValue();
        vcf_file = argVcf.getValue();
        out_dir = argOutDir.getValue();
        interval_file = argInterval.getValue();
        gc_file = argGcfile.getValue();
    }
    catch (TCLAP::ArgException &e)
    {
        std::cerr << "Error: " << e.error() << " for arg " << e.argId() << std::endl;
        abort();
    }
    
    // Generate summary stats from BAM/CRAM
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
    
    if (out_dir.back() != '/')
        out_dir += '/';
    

    std::string pileup_name = out_dir + sample_id + ".pileup";
    std::string varfile_name = out_dir + sample_id + ".var";
    std::string idxfile_name = out_dir + sample_id + ".idx";
    
    std::cerr << "Writing to " << pileup_name << std::endl;
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
        curr_pos += pup.write_depth(b.depth_interval[i], gc.n_interval[i]);
    }
    
    int sp_idx = 0;
    int rp_idx = 0;
    int lclip_idx = 0;
    int rclip_idx = 0;
    
    int prev_sp = 0;
    int prev_rp = 0;
    int prev_lclip = 0;
    int prev_rclip = 0;
    
    int cnt_rp = 0;
    int cnt_sp = 0;
    int cnt_lclip = 0;
    int cnt_rclip = 0;
    
    for(int i=1;i<=gc.num_chr; ++i)
    {
        int N = ceil((double)gc.chr_size[i] / 10000.0) ;
        
        for(int j=1;j<=N;++j)
        {
            // TODO: Write chr:pos at the beginning of each interval, both in the pileup and index file?
            // 5 byte * 3G / 10K = 1.5 MB, negligible increase in file size
            
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
            
            //TODO: Templatize these...
            
            // L-Clip
            uint32_t n_lclip = 0;
            uint32_t n_ldrop = 0;
            while(lclip_idx < (int)b.vec_lclip.size() && b.vec_lclip[lclip_idx].chrnum == i && b.vec_lclip[lclip_idx].pos <= j*10000)
            {
                if (b.vec_lclip[lclip_idx].b_drop == true)
                    n_ldrop++;
                lclip_idx ++;
            }
            if (lclip_idx - prev_lclip - n_ldrop < 0)
            {
                std::cerr << "Wrong lclip index while writing... " << std::endl;
                exit(1);
            }
            else
            {
                n_lclip = (uint32_t) lclip_idx - prev_lclip - n_ldrop;
            }
            curr_pos += pup.write_uint32(n_lclip);
            for(int k=prev_lclip; k<lclip_idx; ++k)
            {
                if (b.vec_lclip[k].b_drop == false)
                {
                    curr_pos += pup.write_softclip(b.vec_lclip[k]);
                    cnt_lclip ++;
                }
            }
            prev_lclip = lclip_idx;
            
            // R-Clip
            uint32_t n_rclip = 0;
            uint32_t n_rdrop = 0;
            while(rclip_idx < (int)b.vec_rclip.size() && b.vec_rclip[rclip_idx].chrnum == i && b.vec_rclip[rclip_idx].pos <= j*10000)
            {
                if (b.vec_rclip[rclip_idx].b_drop == true)
                    n_rdrop++;
                rclip_idx ++;
            }
            if (rclip_idx - prev_rclip - n_rdrop < 0)
            {
                std::cerr << "Wrong rclip index while writing... " << std::endl;
                exit(1);
            }
            else
            {
                n_rclip = (uint32_t) rclip_idx - prev_rclip - n_rdrop;
            }
            curr_pos += pup.write_uint32(n_rclip);
            for(int k=prev_rclip; k<rclip_idx; ++k)
            {
                if (b.vec_rclip[k].b_drop == false)
                {
                    curr_pos += pup.write_softclip(b.vec_rclip[k]);
                    cnt_rclip ++;
                }
            }
            prev_rclip = rclip_idx;
            
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
