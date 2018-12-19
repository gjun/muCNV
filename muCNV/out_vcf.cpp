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
	fprintf(fp,"##source=muCNV_pipeline_v0.9.2\n");
	fprintf(fp,"##INFO=<ID=AC,Number=1,Type=Integer,Description=\"Number of alternative allele\">\n");
	fprintf(fp,"##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">\n");
	fprintf(fp,"##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n");
	fprintf(fp,"##INFO=<ID=CALLRATE,Number=1,Type=Float,Description=\"Call rate\">\n");
	fprintf(fp,"##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">\n");
	fprintf(fp,"##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n");
	fprintf(fp,"##INFO=<ID=DP,Number=1,Type=String,Description=\"1-D Depth clustering\">\n");
	fprintf(fp,"##INFO=<ID=DP2,Number=1,Type=String,Description=\"2-D Depth clustering\">\n");
	fprintf(fp,"##INFO=<ID=READ,Number=0,Type=String,Description=\"Read pair based genotyping\">\n");
	fprintf(fp,"##INFO=<ID=Biallelic,Number=0,Type=String,Description=\"Biallelic variant\">\n");
	fprintf(fp,"##ALT=<ID=DEL,Description=\"Deletion\">\n");
	fprintf(fp,"##ALT=<ID=DUP,Description=\"Duplication\">\n");
	fprintf(fp,"##ALT=<ID=INV,Description=\"Inversion\">\n");
	fprintf(fp,"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
	fprintf(fp,"##FORMAT=<ID=CN,Number=1,Type=Integer,Description=\"Copy Number\">\n");
	fprintf(fp,"##FORMAT=<ID=DP,Number=1,Type=Float,Description=\"Normalized depth\">\n");
	fprintf(fp,"##FORMAT=<ID=DD,Number=A,Type=Float,Description=\"Normalized depth in segments\">\n");
	fprintf(fp,"##FORMAT=<ID=RP,Number=1,Type=Integer,Description=\"Number of supporting read pairs\">\n");
	fprintf(fp,"##FORMAT=<ID=SP,Number=1,Type=Integer,Description=\"Number of supporting split reads\">\n");
	fprintf(fp,"##FORMAT=<ID=SC,Number=1,Type=Integer,Description=\"Number of supporting soft clips\">\n");
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
	{
        af = (double)G.ac/(2.0*G.ns);
	}
    
	fprintf(fp, "SVTYPE=%s;END=%d;SVLEN=%d;AC=%d;NS=%d;CALLRATE=%.2f:AF=%f",  svtype, S.end, S.len, G.ac, G.ns, G.ns / (float)G.gt.size(), af);
	fprintf(fp, ";%s", G.info.c_str());

    if (G.pd_flag)
    {
        fprintf(fp, ";PD");
    }
    fprintf(fp, ";PRE=(%.2f,%2f)", G.dp_pre_mean, G.dp_pre_std);
    fprintf(fp, ";POST=(%.2f,%2f)", G.dp_post_mean, G.dp_post_std);
    if (G.dp_flag)
	{
        fprintf(fp, ";DP=(%.2f:%.2f:%.2f", G.gmix.Comps[0].Mean, G.gmix.Comps[0].Stdev, G.gmix.Comps[0].Alpha);
		for(int j=1;j<G.gmix.n_comp;++j)
		{
			fprintf(fp,"::%.2f:%.2f:%.2f", G.gmix.Comps[j].Mean, G.gmix.Comps[j].Stdev, G.gmix.Comps[j].Alpha);
		}
		fprintf(fp, ");DPoverlap=%.2f", G.gmix.p_overlap);
	}
    if (G.dp2_flag)
	{
        fprintf(fp, ";DP2=(%.2f,%.2f:%.2f,%.2f:%.2f", G.gmix2.Comps[0].Mean[0], G.gmix2.Comps[0].Mean[1], G.gmix2.Comps[0].Cov[0], G.gmix2.Comps[0].Cov[3], G.gmix2.Comps[0].Alpha);
		for(int j=1;j<G.gmix2.n_comp;++j)
		{
			fprintf(fp, "::%.2f,%.2f:%.2f,%.2f:%.2f", G.gmix2.Comps[j].Mean[0], G.gmix2.Comps[j].Mean[1], G.gmix2.Comps[j].Cov[0], G.gmix2.Comps[j].Cov[3], G.gmix2.Comps[j].Alpha);
		}
		fprintf(fp, "):DP2overlap=%2f", G.gmix2.p_overlap);
	}
    if (G.read_flag)
	{
        fprintf(fp, ";READ");
	}
	if (G.b_biallelic)
	{
		fprintf(fp, ";Biallelic");
	}

    fprintf(fp, "\tGT:CN:DP:DD:RP:SP:SC");

    for (int i=0; i<(int)G.gt.size(); ++i)
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
		fprintf(fp, ":%.2f",D.dps[2][i]);

		if (S.svtype == DEL || S.svtype == CNV || S.svtype == DUP )
		{
			fprintf(fp, ":%.2f", D.dps[0][i]);
			if (D.dps.size() == 5)
			{
				fprintf(fp, ",%.2f,%.2f", D.dps[3][i], D.dps[4][i]);
			}
			else
			{
				fprintf(fp, ",.,.");
			}
			fprintf(fp, ",%.2f", D.dps[1][i]);
		}
		else
		{
			fprintf(fp, ":,");
		}

		if (S.svtype == DEL)
		{
			fprintf(fp, ":%d:%d:%d", D.rdstats[i].n_pre_FR + D.rdstats[i].n_post_FR, D.rdstats[i].n_pre_clip_in + D.rdstats[i].n_post_clip_in,  D.rdstats[i].n_pre_split_in + D.rdstats[i].n_post_split_in);
		}
		else if (S.svtype == DUP || S.svtype==CNV)
		{
			fprintf(fp, ":%d:%d:%d", D.rdstats[i].n_pre_RF + D.rdstats[i].n_post_RF, D.rdstats[i].n_pre_clip_out + D.rdstats[i].n_post_clip_out, D.rdstats[i].n_pre_split_out + D.rdstats[i].n_post_split_out);
		}
		else if (S.svtype == INV)
		{
			fprintf(fp, ":%d:%d:%d", D.rdstats[i].n_pre_FF + D.rdstats[i].n_post_RR, D.rdstats[i].n_pre_clip_in + D.rdstats[i].n_pre_clip_out + D.rdstats[i].n_post_clip_in + D.rdstats[i].n_post_clip_out, D.rdstats[i].n_pre_split_out + D.rdstats[i].n_post_split_out + D.rdstats[i].n_pre_split_in + D.rdstats[i].n_pre_split_out);
		}
		else if (S.svtype == INS)
		{
			fprintf(fp, ":%d:%d", D.rdstats[i].n_pre_INS, D.rdstats[i].n_pre_clip_in + D.rdstats[i].n_post_clip_in);
		}

    }
	fprintf(fp, "\n");
}

void OutVcf::print(std::string &ln)
{
	fprintf(fp, "%s\n", ln.c_str());
}
