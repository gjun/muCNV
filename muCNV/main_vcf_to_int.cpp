//
//  vcf_to_idx.cpp
//  muCNV
//
//  Created by Goo Jun on 11/22/18.
//  Copyright © 2018 Goo Jun. All rights reserved.
//

#include <stdio.h>
// TCLAP headers
#include "tclap/CmdLine.h"
#include "tclap/Arg.h"

#include "in_vcf.h"

int main_vcf_to_int(int argc, char** argv)
{
    std::string vcf_file;
    std::string interval_file;
    bool bPrint;
    bool bReverse;
    std::vector<sv> vec_sv;
    std::vector<breakpoint> vec_bp;
    
    try
    {
        TCLAP::CmdLine cmd("Command description message", ' ', "0.06");
        
        TCLAP::ValueArg<std::string> argVcf("v","vcf","VCF file containing candidate SVs",true,"","string");
        TCLAP::ValueArg<std::string> argInterval("i","interVal", "Binary interval file containing candidate SVs", false, "", "string");
        TCLAP::SwitchArg switchPrint("p","print", "Print out SV variants", cmd, false);
        TCLAP::SwitchArg switchReverse("r","reverse", "Generate VCF from interval", cmd, false);


        cmd.add(argVcf);
        cmd.add(argInterval);
        
        cmd.parse(argc, argv);

        vcf_file = argVcf.getValue();
        interval_file = argInterval.getValue();
        bPrint = switchPrint.getValue();
        bReverse = switchReverse.getValue();
    }
    catch (TCLAP::ArgException &e)
    {
        std::cerr << "Error: " << e.error() << " for arg " << e.argId() << std::endl;
        abort();
    }
    // TODO: error checking
    if (!bReverse)
    {
        read_svs_from_vcf(vcf_file, vec_bp, vec_sv);
        if (interval_file != "")
        {
            write_interval(interval_file, vec_sv);
        }
        if (bPrint)
        {
            for(int i=0; i<(int)vec_sv.size(); ++i)
            {
                printf("%d\t", i);
                vec_sv[i].print(stdout);
                printf("\n");
            }

        }
    }
    else
    {
        read_svs_from_intfile(interval_file, vec_bp, vec_sv);
        write_svs_into_vcf(vcf_file, vec_sv);
    }
    return 0;
}
