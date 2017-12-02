//
//  vcf.cpp
//  muCNV
//
//  Created by Goo Jun on 7/24/17.
//  Copyright © 2017 Goo Jun. All rights reserved.
//

#include "muCNV.h"
#include <stdio.h>
#include <math.h>

void outvcf::open(string &fname)
{
	fp = fopen(fname.c_str(), "w");
	varcnt = 0;
}

void outvcf::close()
{
	fclose(fp);
}

void outvcf::write_header(vector<string> &sampleIDs)
{
	fprintf(fp,"##fileformat=VCFv4.1\n");
	fprintf(fp,"##source=UM_CGI_CNV_pipeline_v0.1\n");
	fprintf(fp,"##INFO=<ID=AC,Number=1,Type=Integer,Description=\"Number of alternative allele\">\n");
	fprintf(fp,"##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">\n");
	fprintf(fp,"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n");
	fprintf(fp,"##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n");
	fprintf(fp,"##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">\n");
	fprintf(fp,"##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">\n");
	fprintf(fp,"##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n");
	//	fprintf(fp,"##INFO=<ID=CLUS,Number=1,Type=IntegerDescription=\"Number of depth clusters\">\n");
	//	fprintf(fp,"##INFO=<ID=MEAN,Number=.,Type=,Description=\"Means of each depth cluster\">\n");
	//	fprintf(fp,"##INFO=<ID=STDEV,Number=.,Type=,Description=\"Std. Dev. of each depth cluster\">\n");
	//	fprintf(fp,"##INFO=<ID=PR,Number=.,Type=,Description=\"Probability of each depth cluster\">\n");
	//	fprintf(fp,"##ALT=<ID=CNV,Description=\"Copy number variable region\">\n");
	fprintf(fp,"##ALT=<ID=DEL,Description=\"Deletion\">\n");
	fprintf(fp,"##ALT=<ID=DUP,Description=\"Duplication\">\n");
	fprintf(fp,"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
	fprintf(fp,"##FORMAT=<ID=CN,Number=1,Type=Integer,Description=\"Copy Number\">\n");
	fprintf(fp,"##FORMAT=<ID=CNQ,Number=G,Type=Integer,Description=\"Copy Number Quality\">\n");
	fprintf(fp,"##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n");
	fprintf(fp,"##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n");
	fprintf(fp,"##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Genotype Likelihood\">\n");
	fprintf(fp,"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
	for(unsigned j=0; j<sampleIDs.size();++j)
	{
		fprintf(fp,"\t%s",sampleIDs[j].c_str());
	}
	fprintf(fp,"\n");
	fflush(fp);
	
}


void outvcf::write_del(sv& interval, vector<int>& gt, vector<int>& GQ, int ac, int ns, vector<double>& X, vector<double>& AvgDepth, vector<Gaussian>& C, double be, bool bFilter)
{
	// For Deletions
	int n_comp=(int)C.size();
	int chrnum = interval.chrnum;
	string chr = interval.chr;
	int pos = interval.pos;
	int svend = interval.end;
	
	fprintf(fp,"%s\t", chr.c_str());
	fprintf(fp,"%d\t", pos);
	if (interval.svtype == "DEL")
	{
		fprintf(fp,"muCNV%d\t.\t<DEL>\t.\t", ++varcnt);
	}
	else
	{
		fprintf(fp,"muCNV%d\t.\t<INV>\t.\t", ++varcnt);
	}

	if (bFilter)
	{
		fprintf(fp, "PASS");
	}
	else
	{
		fprintf(fp, "FAIL");
	}
	

	double af=0;
	if (ns>0)
	{
		af = 0.5*ac/ns;
	}

	if (interval.svtype == "DEL")
	{
		fprintf(fp,"\tIMPRECISE;CIPOS=%d,%d;CIEND=%d,%d;VT=SV;END=%d;SVLEN=-%d;AC=%d;AF=%1.4f;AN=%d;NS=%d;P_OVERLAP=%1.8f;SVTYPE=DEL", interval.ci_pos.first, interval.ci_pos.second, interval.ci_end.first, interval.ci_end.second, svend, svend-pos, ac, af, ns*2, ns, be);
	}
	else
	{
		fprintf(fp,"\tIMPRECISE;CIPOS=%d,%d;CIEND=%d,%d;VT=SV;END=%d;SVLEN=-%d;AC=%d;AF=%1.4f;AN=%d;NS=%d;P_OVERLAP=%1.8f;SVTYPE=INV", interval.ci_pos.first, interval.ci_pos.second, interval.ci_end.first, interval.ci_end.second, svend, svend-pos, ac, af, ns*2, ns, be);
	}
	fprintf(fp,";CLUS=%d", n_comp);
	
	fprintf(fp,";MEAN=%1.4f",C[0].Mean);
	for(int i=1;i<n_comp;++i)
	{
		fprintf(fp,",%1.4f",C[i].Mean);
	}
	
	fprintf(fp,";STDEV=%1.4f",C[0].Stdev);
	for(int i=1;i<n_comp;++i)
	{
		fprintf(fp,",%1.4f",C[i].Stdev);
	}
	
	fprintf(fp,";PR=%1.4f",C[0].Alpha);
	for(int i=1;i<n_comp;++i)
	{
		fprintf(fp,",%1.4f",C[i].Alpha);
	}
	
	fprintf(fp,"\tGT:DP:GQ" );
	for(unsigned j=0; j<gt.size(); ++j)
	{
		switch(gt[j])
		{
			case 0:
				fprintf(fp, "\t./.:");
				break;
			case 1:
				fprintf(fp, "\t0/0:");
				break;
			case 2:
				fprintf(fp, "\t0/1:");
				break;
			case 3:
				fprintf(fp, "\t1/1:");
				break;
		}
		fprintf(fp, "%u:", (unsigned)round(X[j]*AvgDepth[j]));
		fprintf(fp, "%u", GQ[j]);
	}
	fprintf(fp, "\n");
	fflush(fp);
}


void outvcf::write_cnv(sv& interval, vector<int>& gt, vector<int>& GQ, int ac, int ns, vector<double>& X, vector<double>& AvgDepth, vector<Gaussian>& C, double be, bool bFilter)
{
	int n_comp=(int)C.size();
	string chr = interval.chr;
	int pos = interval.pos;
	int svend = interval.end;
	int alts[11] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};;
	
	for(int j=0;j<(int)gt.size(); ++j)
	{
		if (gt[j] >0 && gt[j]<=10)
		{
			alts[gt[j]] = 1;
		}
	}

	fprintf(fp,"%s\t", chr.c_str());
	fprintf(fp,"%d\t", pos);
	
	fprintf(fp,"muCNV%d\t.\t",++varcnt);

	bool flag=false;
	for(int j=1; j<11; ++j)
	{
		if (alts[j])
		{
			if (flag)
			{
				fprintf(fp,",<CN%d>",j);
			}
			else
			{
				flag=true;
				fprintf(fp,"<CN%d>",j);
			}
		}
	}
	if (!flag)
	{
		fprintf(fp,"<CNV>");
	}
	fprintf(fp, "\t.\t");

	if (bFilter)
	{
		fprintf(fp, "PASS");
	}
	else
	{
		fprintf(fp, "FAIL");
	}
	double af=0;
	if (ns>0)
	{
		af = 0.5*ac/ns;
	}
	fprintf(fp,"\tIMPRECISE;CIPOS=%d,%d;CIEND=%d,%d;VT=SV;END=%d;SVLEN=%d;AC=%d;AF=%1.4f;AN=%d;NS=%d;SVTYPE=CNV", interval.ci_pos.first, interval.ci_pos.second, interval.ci_end.first, interval.ci_end.second, svend, svend-pos, ac, af, ns*2, ns);
	
	fprintf(fp,";CLUS=%d", n_comp);
	
	fprintf(fp,";MEAN=%1.4f",C[0].Mean);
	for(int i=1;i<n_comp;++i)
	{
		fprintf(fp,",%1.4f",C[i].Mean);
	}
	
	fprintf(fp,";STDEV=%1.4f",C[0].Stdev);
	for(int i=1;i<n_comp;++i)
	{
		fprintf(fp,",%1.4f",C[i].Stdev);
	}
	
	fprintf(fp,";PR=%1.4f",C[0].Alpha);
	for(int i=1;i<n_comp;++i)
	{
		fprintf(fp,",%1.4f",C[i].Alpha);
	}
	
	fprintf(fp,"\tGT:CN:DP:GQ" );


	for(int j=0; j<(int)gt.size(); ++j)
	{
		switch(gt[j])
		{
			case 0:
				fprintf(fp,"\t./.:");
				break;
			case 2:
				fprintf(fp,"\t0/0:");
				break;
			case 3:
				fprintf(fp,"\t0/1:");
				break;
			case 4:
				fprintf(fp,"\t1/1:");
				break;
			default:
				fprintf(fp,"\t./.:");
				break;
		}
		if (gt[j] == 0)
		{
			fprintf(fp, ".:");
		}
		else
		{
			if (gt[j]<=4)
				fprintf(fp, "%d:",gt[j]);
			else
				fprintf(fp, ".:");
		}
		
		fprintf(fp, "%u:", (unsigned)round(X[j]*AvgDepth[j]));
		fprintf(fp, "%u", GQ[j]);
	}
	fprintf(fp, "\n");
	fflush(fp);
}


