//
//  main_varpileup
//  muCNV
//
//  Created by Goo Jun on Jan 2020
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

int main_varpileup(int argc, char** argv)
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
    int chr;
    int startpos;
    int endpos;
    
    try
    {
        TCLAP::CmdLine cmd("Command description message", ' ', "0.06");
        
        TCLAP::ValueArg<std::string> argBam("b","bam","Input BAM/CRAM file name",true,"","string");
        TCLAP::ValueArg<std::string> argVcf("v","vcf","VCF file containing candidate SVs",false,"","string");
        TCLAP::ValueArg<int> argChr("c","chr","Chromosome number (1-24) if pileup contains only a single chromosome",false,0,"integer (1-24)");
        TCLAP::ValueArg<std::string> argInterval("V","interVal", "Binary interval file containing candidate SVs", false, "", "string");
        TCLAP::ValueArg<std::string> argSampleID("s","sample","Sample ID, also used for output filenames (.pileup, .var, .idx)",true,"","string");
        TCLAP::ValueArg<std::string> argOutDir("o","outdir","Output directory, current directory if omitted",false,".","string");
        TCLAP::ValueArg<std::string> argGcfile("f","gcFile","File containing GC content information (default: GRCh38.gc)",false, "GRCh38.gc", "string");
        
        cmd.add(argBam);
        cmd.add(argVcf);
        cmd.add(argChr);
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
        chr = argChr.getValue();
    }
    catch (TCLAP::ArgException &e)
    {
        std::cerr << "Error: " << e.error() << " for arg " << e.argId() << std::endl;
        abort();
    }
    
    // Generate summary stats from BAM/CRAM
    std::cerr << "Processing individual BAM/CRAM file to genearte var pileup on a region" << std::endl;
    
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
    
    
    std::string varfile_name = out_dir + sample_id + ".var";
    if (chr > 0)
        varfile_name = out_dir + sample_id + ".chr20.var";


    std::cerr << "Writing to " << varfile_name << std::endl;

    //Initialize BAM/CRAM
    BamCram b;
    b.initialize_idx(bam_file);
    std::cerr << "BAM/CRAM file initialized" << std::endl;
    // Read BAM/CRAM
    startpos = 1;
    endpos = gc.chr_size[chr];
    b.read_vardepth(gc, vec_bp, vec_sv, chr, startpos, endpos);
    b.postprocess_depth(vec_sv);

    // Write Pileup 
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
        //printf("%d\t%d-%d\t%2f\n", vec_sv[i].chrnum, vec_sv[i].pos, vec_sv[i].end, (double)vec_sv[i].dp/32.0);
    }
    var_file.close();
    
    std::cerr << "Finished without an error" << std::endl;
    return 0;
}
