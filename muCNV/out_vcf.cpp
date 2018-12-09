//
//  vcf.cpp
//  muCNV
//
//  Created by Goo Jun on 7/24/17.
//  Copyright © 2017 Goo Jun. All rights reserved.
//

#include "out_vcf.h"
#include <math.h>

void OutVcf::open(std::string &fname)
{
	fp = fopen(fname.c_str(), "w");
	varcnt = 0;
}

void OutVcf::close()
{
	fclose(fp);
}

void OutVcf::write_header(std::vector<std::string> &sampleIDs)
{
	fprintf(fp,"##fileformat=VCFv4.1\n");
	fprintf(fp,"##source=muCNV_pipeline_v0.9.1\n");
	fprintf(fp,"##INFO=<ID=AC,Number=1,Type=Integer,Description=\"Number of alternative allele\">\n");
	fprintf(fp,"##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">\n");
	fprintf(fp,"##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n");
	fprintf(fp,"##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">\n");
	fprintf(fp,"##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n");
	fprintf(fp,"##ALT=<ID=DEL,Description=\"Deletion\">\n");
	fprintf(fp,"##ALT=<ID=DUP,Description=\"Duplication\">\n");
	fprintf(fp,"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
	fprintf(fp,"##FORMAT=<ID=CN,Number=1,Type=Integer,Description=\"Copy Number\">\n");
	fprintf(fp,"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
	for(unsigned j=0; j<sampleIDs.size();++j)
	{
		fprintf(fp,"\t%s",sampleIDs[j].c_str());
	}
	fprintf(fp,"\n");
	fflush(fp);
	
}


void OutVcf::write_sv(sv &S, SvData &D, SvGeno &G)
{
	const char *svtype = svTypeName(S.svtype).c_str();

	if (S.chrnum < 23)
	{
    	fprintf(fp, "%d",  S.chrnum);
	}
	else if (S.chrnum == 23)
	{
		fprintf(fp, "X");
	}
	else if (S.chrnum == 24)
	{
		fprintf(fp, "Y");
	}

	fprintf(fp, "\t%d\t%s_%d:%d-%d\t.\t<%s>\t.\t", S.pos, svtype, S.chrnum, S.pos, S.end, svtype);

    if (G.b_pass)
    {
        fprintf(fp, "PASS\t");
    }
    else
    {
        fprintf(fp, "FAIL\t");
    }

    double af = 0;
    if (G.ns>0)
        af = (double)G.ac/(2.0*G.ns);
    
	fprintf(fp, "SVTYPE=%s;END=%d;SVLEN=%d;AC=%d;NS=%d;AF=%f",  svtype, S.end, S.len, G.ac, G.ns, af);
	fprintf(fp, "%s", G.info.c_str());

    if (G.dp_flag)
        fprintf(fp, ";DP");
    if (G.dp2_flag)
        fprintf(fp, ";DP2");
    if (G.read_flag)
        fprintf(fp, ";READ");

    fprintf(fp, "\tGT:CN");

    for (int i=0; i<G.gt.size(); ++i)
    {
        switch(G.gt[i])
        {
            case 0:
                fprintf(fp, "\t0/0");
                break;
            case 1:
                fprintf(fp, "\t0/1");
                break;
            case 2:
                fprintf(fp, "\t1/1");
                break;
            default:
                fprintf(fp, "\t.");
                break;
        }
        if (G.cn[i]<0)
            fprintf(fp, ":.");
        else
            fprintf(fp, ":%d",G.cn[i]);

    }
	fprintf(fp, "\n");
}

void OutVcf::print(std::string &ln)
{
	fprintf(fp, "%s\n", ln.c_str());
}
