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
	fprintf(fp,"##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of the structural variant\">\n");
	fprintf(fp,"##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n");
    fprintf(fp,"##INFO=<ID=VarID,Number=1,Type=Float,Description=\"Variant ID\">\n");
	fprintf(fp,"##INFO=<ID=DP,Number=1,Type=String,Description=\"1-D Depth clustering\">\n");
    fprintf(fp,"##INFO=<ID=DPOverlap,Number=1,Type=Float,Description=\"Overlap between 1-D depth clusters\">\n");
	fprintf(fp,"##INFO=<ID=DP2,Number=1,Type=String,Description=\"2-D Depth clustering\">\n");
    fprintf(fp,"##INFO=<ID=DP2Overlap,Number=1,Type=Float,Description=\"Overlap between 2-D depth clusters\">\n");
	fprintf(fp,"##INFO=<ID=SPLIT,Number=1,Type=String,Description=\"Breakpoints estimated by split reads\">\n");
    fprintf(fp,"##INFO=<ID=CLIP,Number=1,Type=String,Description=\"Breakpoints estimated by soft clips\">\n");
	fprintf(fp,"##INFO=<ID=PRE, Number=1,Type=String,Description=\"Read depth statistic before SV\">\n");
	fprintf(fp,"##INFO=<ID=POST, Number=1,Type=String,Description=\"Read depth statistic after SV\">\n");
	fprintf(fp,"##INFO=<ID=Biallelic,Number=0,Type=String,Description=\"Biallelic variant\">\n");
	fprintf(fp,"##ALT=<ID=DEL,Description=\"Deletion\">\n");
	fprintf(fp,"##ALT=<ID=DUP,Description=\"Duplication\">\n");
	fprintf(fp,"##ALT=<ID=INV,Description=\"Inversion\">\n");
	fprintf(fp,"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
	fprintf(fp,"##FORMAT=<ID=CN,Number=1,Type=Integer,Description=\"Copy Number\">\n");
	fprintf(fp,"##FORMAT=<ID=DP,Number=1,Type=Float,Description=\"Normalized depth\">\n");
	fprintf(fp,"##FORMAT=<ID=DD,Number=A,Type=Float,Description=\"Normalized depth in before, inside first, inside second, and after SV\">\n");
	fprintf(fp,"##FORMAT=<ID=RP,Number=1,Type=Integer,Description=\"Number of supporting read pairs\">\n");
	fprintf(fp,"##FORMAT=<ID=SP,Number=1,Type=Integer,Description=\"Number of supporting split reads\">\n");
	fprintf(fp,"##FORMAT=<ID=SC,Number=1,Type=Integer,Description=\"Number of supporting soft clips around starting position of SV\">\n");
    fprintf(fp,"##FORMAT=<ID=EC,Number=1,Type=Integer,Description=\"Number of supporting soft clips around end position of SV\">\n");
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
    
	fprintf(fp, "SVTYPE=%s;END=%d;SVLEN=%d;AC=%d;NS=%d;CALLRATE=%.2f;AF=%f",  svtype, S.end, S.len, G.ac, G.ns, G.ns / (float)G.gt.size(), af);
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

	if (G.b_biallelic)
	{
		fprintf(fp, ";Biallelic");
	}
    
    if (G.split_flag)
    {
        fprintf(fp, ";SPLIT=(");
        fprintf(fp, "%d,%d", (int)(D.vec_break_clusters[0].start_mean + 0.5), (int)(D.vec_break_clusters[0].end_mean+0.5));
        for(int j=1; j<D.vec_break_clusters.size(); ++j)
        {
            fprintf(fp, ";%d,%d", (int)(D.vec_break_clusters[j].start_mean + 0.5), (int)(D.vec_break_clusters[j].end_mean+0.5));
        }
        fprintf(fp, ")");
    }
    else if (G.clip_flag)
    {
        fprintf(fp, ";CLIP=(");
        fprintf(fp, "%d,%d", (int)(D.vec_break_clusters[0].start_mean + 0.5), (int)(D.vec_break_clusters[0].end_mean+0.5));
        for(int j=1; j<D.vec_break_clusters.size(); ++j)
        {
            fprintf(fp, ";%d,%d", (int)(D.vec_break_clusters[j].start_mean + 0.5), (int)(D.vec_break_clusters[j].end_mean+0.5));
        }
        fprintf(fp, ")");
    }
    
    
    fprintf(fp, "\tGT:CN:DP:DD:SP:RP:SC:EC");
    
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

		if (S.svtype == DEL || S.svtype == DUP || S.svtype == CNV)
		{
            fprintf(fp, ":%d:%d:%d:%d", (int)G.split_cnts[i], (int)G.rp_cnts[i], (int)G.start_clips[i], (int)G.end_clips[i]);
		}
		else if (S.svtype == INV)
		{
            fprintf(fp, ":%d:%d:%d:%d", (int)G.split_cnts[i], (int)G.rp_cnts[i], (int)G.start_clips[i], (int)G.end_clips[i]);
		}
    }
	fprintf(fp, "\n");
}

void OutVcf::print(std::string &ln)
{
	fprintf(fp, "%s\n", ln.c_str());
}
