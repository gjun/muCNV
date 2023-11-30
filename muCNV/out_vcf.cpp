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

void OutVcf::write_header(std::vector<std::string> &sampleIDs, std::vector<bool> &mask, GcContent &gc)
{
	fprintf(fp,"##fileformat=VCFv4.1\n");
	fprintf(fp,"##source=muCNV_pipeline_v0.9.6\n");
	for(int i=1; i<=22; ++i)
	{
		fprintf(fp, "##contig=<ID=chr%d,length=%d>\n", i, gc.chr_size[i]);
	}
	fprintf(fp,"##contig=<ID=chrX,length=%d>\n", gc.chr_size[23]);
	fprintf(fp,"##contig=<ID=chrY,length=%d>\n", gc.chr_size[24]);
	fprintf(fp,"##INFO=<ID=AC,Number=1,Type=Integer,Description=\"Number of alternative allele\">\n");
	fprintf(fp,"##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">\n");
	fprintf(fp,"##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n");
	fprintf(fp,"##INFO=<ID=CALLRATE,Number=1,Type=Float,Description=\"Call rate\">\n");
	fprintf(fp,"##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">\n");
	fprintf(fp,"##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of the structural variant\">\n");
	fprintf(fp,"##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n");
 //   fprintf(fp,"##INFO=<ID=VarID,Number=1,Type=Float,Description=\"Variant ID\">\n");
	fprintf(fp,"##INFO=<ID=DP,Number=1,Type=String,Description=\"1-D Depth clustering\">\n");
    fprintf(fp,"##INFO=<ID=DPOverlap,Number=1,Type=Float,Description=\"Overlap between 1-D depth clusters\">\n");
	fprintf(fp,"##INFO=<ID=DP2,Number=1,Type=String,Description=\"2-D Depth clustering\">\n");
	fprintf(fp,"##INFO=<ID=DP2Overlap,Number=1,Type=Float,Description=\"Overlap between 2-D depth clusters\">\n");
	fprintf(fp,"##INFO=<ID=DPCNT,Number=1,Type=String,Description=\"2-D Depth and Read Count clustering\">\n");
    fprintf(fp,"##INFO=<ID=DP2Overlap,Number=1,Type=Float,Description=\"Overlap between 2-D Depth-ReadCount clusters\">\n");
	fprintf(fp,"##INFO=<ID=SPLIT,Number=1,Type=String,Description=\"Breakpoints estimated by split reads\">\n");
    fprintf(fp,"##INFO=<ID=CLIP,Number=1,Type=String,Description=\"Breakpoints estimated by soft clips\">\n");
    fprintf(fp,"##INFO=<ID=PRE_FAIL,Number=0,Type=String,Description=\"Pre-depth distribution is out of bound\">\n");
    fprintf(fp,"##INFO=<ID=POST_FAIL,Number=0,Type=String,Description=\"Post-depth distribution is out of bound\">\n");
    fprintf(fp,"##INFO=<ID=READPAIR,Number=0,Type=String,Description=\"Genotyped by reported breakpoints and read pairs\">\n");
    fprintf(fp,"##INFO=<ID=CNT,Number=1,Type=String,Description=\"Clustering stats by split read or read pair counts\">\n");
    fprintf(fp,"##INFO=<ID=CNTOverlap,Number=1,Type=Float,Description=\"Overlap between 1-D depth clusters using split/readpair counts\">\n");
	fprintf(fp,"##INFO=<ID=PRE, Number=1,Type=String,Description=\"Read depth statistic before SV\">\n");
	fprintf(fp,"##INFO=<ID=POST, Number=1,Type=String,Description=\"Read depth statistic after SV\">\n");
	fprintf(fp,"##INFO=<ID=Biallelic,Number=0,Type=String,Description=\"Biallelic variant\">\n");
	fprintf(fp,"##INFO=<ID=Biallelic,Number=0,Type=String,Description=\"Biallelic variant\">\n");
	fprintf(fp,"##INFO=<ID=N_MISS,Number=1,Type=Integer,Descrption=\"Number of missing genotypes (in females for chrX)\">\n");
	fprintf(fp,"##INFO=<ID=N_HOMREF,Number=1,Type=Integer,Descrption=\"Number of REF/REF genotypes (in females for chrX)\">\n");
	fprintf(fp,"##INFO=<ID=N_HET,Number=1,Type=Integer,Descrption=\"Number of REF/ALT genotypes (in females for chrX)\">\n");
	fprintf(fp,"##INFO=<ID=N_HOMALT,Number=1,Type=Integer,Descrption=\"Number of ALT/ALT genotypes (in females for chrX)\">\n");
	fprintf(fp,"##INFO=<ID=N_MALE_MISS,Number=1,Type=Integer,Descrption=\"Number of missing genotypes in males for chrX\">\n");
	fprintf(fp,"##INFO=<ID=N_MALE_REF,Number=1,Type=Integer,Descrption=\"Number of REF genotypes in males for chrX\">\n");
	fprintf(fp,"##INFO=<ID=N_MALE_ALT,Number=1,Type=Integer,Descrption=\"Number of ALT genotypes in males for chrX\">\n");
	fprintf(fp,"##INFO=<ID=HWECHISQ,Number=1,Type=Float,Descrption=\"Chi-square value of HWE\">\n");
	fprintf(fp,"##ALT=<ID=DEL,Description=\"Deletion\">\n");
	fprintf(fp,"##ALT=<ID=DUP,Description=\"Duplication\">\n");
	fprintf(fp,"##ALT=<ID=INV,Description=\"Inversion\">\n");
	fprintf(fp,"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
	fprintf(fp,"##FORMAT=<ID=CN,Number=1,Type=Integer,Description=\"Copy Number\">\n");
	fprintf(fp,"##FORMAT=<ID=DP,Number=1,Type=Float,Description=\"Raw depth\">\n");
	fprintf(fp,"##FORMAT=<ID=ND,Number=1,Type=Float,Description=\"Normalized depth\">\n");
	fprintf(fp,"##FORMAT=<ID=DD,Number=A,Type=Float,Description=\"Normalized depth in before, inside first, inside second, and after SV\">\n");
	fprintf(fp,"##FORMAT=<ID=RP,Number=1,Type=Integer,Description=\"Number of supporting read pairs\">\n");
	fprintf(fp,"##FORMAT=<ID=SP,Number=1,Type=Integer,Description=\"Number of supporting split reads\">\n");
	fprintf(fp,"##FORMAT=<ID=SC,Number=1,Type=Integer,Description=\"Number of supporting soft clips around starting position of SV\">\n");
    fprintf(fp,"##FORMAT=<ID=EC,Number=1,Type=Integer,Description=\"Number of supporting soft clips around end position of SV\">\n");
	fprintf(fp,"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
	for(unsigned j=0; j<sampleIDs.size();++j)
	{
		if (mask[j])
			fprintf(fp,"\t%s",sampleIDs[j].c_str());
	}
	fprintf(fp,"\n");
	fflush(fp);
	
}


void OutVcf::write_sv(sv &S, SvData &D, SvGeno &G)
{
	const char *svtype = svTypeName(S.svtype).c_str();
    
    // TODO: convert hard-coded chr-names using real chr names from BAM/CRAM header - this should be stored in GC content file
    // current code works only with GRCh38
    
	if (S.chrnum < 23)
	{
    	fprintf(fp, "chr%d",  S.chrnum);
	}
	else if (S.chrnum == 23)
	{
		fprintf(fp, "chrX");
	}
	else if (S.chrnum == 24)
	{
		fprintf(fp, "chrY");
	}
    else if (S.chrnum == 25)
    {
        fprintf(fp, "chrM");
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

	if (S.svtype == DEL)
		S.len = -S.len;	
	fprintf(fp, "SVTYPE=%s;END=%d;SVLEN=%d;AC=%d;NS=%d;CALLRATE=%.2f;AF=%f",  svtype, S.end, S.len, G.ac, G.ns, G.ns / (float)G.n_effect, G.af);
    fprintf(fp, ";N_MISS=%d;N_HOMREF=%d;N_HET=%d;N_HOMALT=%d;HWECHISQ=%.2f", G.gts[0], G.gts[1], G.gts[2], G.gts[3], G.chisq);
    if (S.chrnum == 23)
    {
        fprintf(fp, ";N_MALE_MISS=%d;N_MALE_REF=%d;N_MALE_ALT=%d", G.m_gts[0], G.m_gts[1], G.m_gts[2]+G.m_gts[3]);
    }

    fprintf(fp, ";%s", G.info.c_str());
    
    if (G.pd_flag)
    {
        fprintf(fp, ";PD");
    }
    fprintf(fp, ";PRE=(%.2f,%2f)", G.dp_pre_mean, G.dp_pre_std);
    fprintf(fp, ";POST=(%.2f,%2f)", G.dp_post_mean, G.dp_post_std);
    if (G.cnt_flag)
    {
        fprintf(fp, ";CNT=(%.3f:%.3f:%.3f", G.gmix.Comps[0].Mean, G.gmix.Comps[0].Stdev, G.gmix.Comps[0].Alpha);
        for(int j=1;j<G.gmix.n_comp;++j)
        {
            fprintf(fp,"::%.3f:%.3f:%.3f", G.gmix.Comps[j].Mean, G.gmix.Comps[j].Stdev, G.gmix.Comps[j].Alpha);
        }
        fprintf(fp, ");CNTOverlap=%.3f", G.gmix.p_overlap);
    }
    if (G.dp_flag)
	{
        fprintf(fp, ";DP=(%.3f:%.3f:%.3f", G.gmix.Comps[0].Mean, G.gmix.Comps[0].Stdev, G.gmix.Comps[0].Alpha);
		for(int j=1;j<G.gmix.n_comp;++j)
		{
			fprintf(fp,"::%.3f:%.3f:%.3f", G.gmix.Comps[j].Mean, G.gmix.Comps[j].Stdev, G.gmix.Comps[j].Alpha);
		}
		fprintf(fp, ");DPOverlap=%.3f", G.gmix.p_overlap);
	}
    if (G.dp2_flag)
	{
        fprintf(fp, ";DP2=(%.3f,%.3f:%.3f,%.3f,%.3f:%.3f", G.gmix2.Comps[0].Mean[0], G.gmix2.Comps[0].Mean[1], G.gmix2.Comps[0].Cov[0], G.gmix2.Comps[0].Cov[1], G.gmix2.Comps[0].Cov[3], G.gmix2.Comps[0].Alpha);
		for(int j=1;j<G.gmix2.n_comp;++j)
		{
			fprintf(fp, "::%.3f,%.3f:%.3f,%.3f,%.3f:%.3f", G.gmix2.Comps[j].Mean[0], G.gmix2.Comps[j].Mean[1], G.gmix2.Comps[j].Cov[0], G.gmix2.Comps[j].Cov[1], G.gmix2.Comps[j].Cov[3], G.gmix2.Comps[j].Alpha);
		}
		fprintf(fp, "):DP2overlap=%.3f", G.gmix2.p_overlap);
	}

    if (G.dpcnt_flag)
	{
        fprintf(fp, ";DPCNT=(%.3f,%.3f:%.3f,%.3f,%.3f:%.3f", G.gmix2.Comps[0].Mean[0], G.gmix2.Comps[0].Mean[1], G.gmix2.Comps[0].Cov[0], G.gmix2.Comps[0].Cov[1], G.gmix2.Comps[0].Cov[3], G.gmix2.Comps[0].Alpha);
		for(int j=1;j<G.gmix2.n_comp;++j)
		{
			fprintf(fp, "::%.3f,%.3f:%.3f,%.3f,%.3f:%.3f", G.gmix2.Comps[j].Mean[0], G.gmix2.Comps[j].Mean[1], G.gmix2.Comps[j].Cov[0], G.gmix2.Comps[j].Cov[1], G.gmix2.Comps[j].Cov[3], G.gmix2.Comps[j].Alpha);
		}
		fprintf(fp, "):DPCNToverlap=%.3f", G.gmix2.p_overlap);
	}

	if (G.b_biallelic)
	{
		fprintf(fp, ";Biallelic");
	}
    
    if (G.split_flag)
    {
        fprintf(fp, ";SPLIT=(");
        fprintf(fp, "%d,%d", (int)(D.vec_break_clusters[0].start_mean + 0.5), (int)(D.vec_break_clusters[0].end_mean+0.5));
        for(int j=1; j<(int)D.vec_break_clusters.size(); ++j)
        {
            fprintf(fp, ":%d,%d", (int)(D.vec_break_clusters[j].start_mean + 0.5), (int)(D.vec_break_clusters[j].end_mean+0.5));
        }
        fprintf(fp, ")");
    }
    else if (G.clip_flag)
    {
        fprintf(fp, ";CLIP=(");
        fprintf(fp, "%d,%d", (int)(D.vec_break_clusters[0].start_mean + 0.5), (int)(D.vec_break_clusters[0].end_mean+0.5));
        for(int j=1; j<(int)D.vec_break_clusters.size(); ++j)
        {
            fprintf(fp, ":%d,%d", (int)(D.vec_break_clusters[j].start_mean + 0.5), (int)(D.vec_break_clusters[j].end_mean+0.5));
        }
        fprintf(fp, ")");
    }
    
    
    fprintf(fp, "\tGT:CN:ND:DP:DD:SP:RP:SC:EC");
    
    for (int i=0; i<(int)G.gt.size(); ++i)
    {
		if (!G.sample_mask[i])
			continue;

		
		if (S.chrnum < 23)
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
					fprintf(fp, "\t./.");
					break;
			}
		}
		else if (S.chrnum == 23)
		{
			if (G.sex[i] == 1)
			{
				switch(G.gt[i])
				{
					case 0:
						fprintf(fp, "\t0");
						break;
					case 1:
						fprintf(fp, "\t1");
						break;
					case 2:
						fprintf(fp, "\t1");
						break;
					default:
						fprintf(fp, "\t.");
						break;
				}
			}
			else if (G.sex[i] == 2)
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
						fprintf(fp, "\t./.");
						break;
				}
			}
			else
			{
				fprintf(fp, "\t./.");
			}
		}

        if (G.cn[i]<0)
            fprintf(fp, ":.");
        else
            fprintf(fp, ":%d",G.cn[i]);	
				
		fprintf(fp, ":%.2f",D.dps[2][i]);
		
		fprintf(fp, ":%.2f",D.raw_dp[i]);

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
