/*
 *   Author: Goo Jun (goo.jun@uth.tmc.edu)
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "genotyper.h"

SvGeno::SvGeno(int n)
{
    n_sample = n;
    ns = 0;
    ac = 0;
    gt.resize(n_sample, -1);
    cn.resize(n_sample, -1);
    b_biallelic = false;
    b_pass = false;
    dp_flag = false;
    dp2_flag = false;
    read_flag = false;
    info = "";
}

SvData::SvData(int n)
{
    n_sample = n;
    rdstats.resize(n);
    var_depth.resize(n);
    dp2.resize(n);
}


void copyComps(std::vector<Gaussian> &C, std::vector<Gaussian> &C0)
{
	C.clear();
	C.resize(C0.size());
	for(unsigned j=0;j<C.size(); ++j)
	{
		C[j].set(C0[j].Mean, C0[j].Stdev);
		C[j].Alpha = C0[j].Alpha;
	}
}

void Genotyper::select_model(GaussianMixture &ret_gmix, std::vector< std::vector<double> > &means, std::vector<double> &x)
{
    double best_bic = DBL_MAX;
    
    // number of models
    for(int m=0; m<means.size(); ++m)
    {
        std::vector<double> s (means[m].size(), 0.05);
        GaussianMixture gmix(means[m], s);
        gmix.EM(x); // fit mixture model
        if (gmix.bic < best_bic && gmix.p_overlap < MAX_P_OVERLAP)
        {
            best_bic = gmix.bic;
            ret_gmix = gmix; // assignment
        }
    }
    return;
}

void Genotyper::select_model(GaussianMixture2 &ret_gmix2, std::vector< std::vector<double> > &means, std::vector<double> &x, std::vector<double> &y)
{
    double best_bic = DBL_MAX;
    
    // number of models
    for(int m=0; m<means.size(); ++m)
    {
        std::vector<double> s (means[m].size(), 0.05);
        GaussianMixture2 gmix2(means[m], s);
        gmix2.EM2(x, y); // fit mixture model
        if (gmix2.bic < best_bic && gmix2.p_overlap < MAX_P_OVERLAP )
        {
            best_bic = gmix2.bic;
            ret_gmix2 = gmix2; // assignment (copy operation)
        }
    }
    return;
}

void Genotyper::call(sv &S, SvData &D, SvGeno &G)
{
    if (S.svtype == DEL)
        call_deletion(S, D, G);
    else if (S.svtype == DUP || S.svtype == CNV)
        call_cnv(S, D, G);
    // TODO: inversion and insertion
}

void Genotyper::call_deletion(sv &S, SvData &D, SvGeno &G)
{
    int n_sample = D.n_sample;
    
	// Fit Gaussian mixture models with 1, 2, and 3 components, compare BIC
    GaussianMixture gmix;
    std::vector< std::vector<double> > means = { {1.0}, {1.0, 0.5}, {0.5, 0.0}, {1.0, 0.5, 0.0}};
    select_model(gmix, means, D.var_depth);
    
    G.b_biallelic = true;
    if (gmix.n_comp > 1)
    {
        // success
        G.dp_flag = true;
        // assign dp-genotypes
        for(int i=0; i<n_sample; ++i)
        {
            G.cn[i] = gmix.assign_copynumber(D.var_depth[i]);
            if (G.cn[i] == 2)
                G.gt[i] = 0; // 0/0
            else if (G.cn[i] == 1)
                G.gt[i] = 1; // 0/1
            else if (G.cn[i] == 0)
                G.gt[i] = 0; // 1/1
        }
    }

    if (D.dp2.size()>1) //dp2 has more than 2 vectors
    {
        // DP100 genotyping
        GaussianMixture2 gmix2;

        G.dp2_flag = true;
        int dp2_idx = 0;
        if (D.dp2.size()==4)
            dp2_idx = 1;
        
        select_model(gmix2, means, D.dp2[dp2_idx], D.dp2[dp2_idx+1]);

        // 2-D genotyping
        if (gmix2.n_comp>1)
        {
            // success
            //assign dp2 genotypes
            for(int i=0; i<n_sample; ++i)
            {
                G.cn[i] = gmix2.assign_copynumber(D.dp2[dp2_idx][i], D.dp2[dp2_idx+1][i]);
                if (G.cn[i] == 2)
                    G.gt[i] = 0; // 0/0
                else if (G.cn[i] == 1)
                    G.gt[i] = 1; // 0/1
                else if (G.cn[i] == 0)
                    G.gt[i] = 0; // 1/1
            }
        }
    }
    
    // Readpair genotyping
    for(int i=0; i<n_sample; ++i)
    {
        if (G.dp_flag || G.dp2_flag)
        {
            if (D.rdstats[i].sv_support() == DEL)
            {
                G.read_flag = true;
                if (D.var_depth[i] > 0.25 && D.var_depth[i] < 0.65 && G.gt[i]== -1) // TEMPORARY
                {
                    G.cn[i] = 1;
                    G.gt[i] = 1;
                }
                else if (D.var_depth[i]<=0.25 && G.gt[i] == -1)
                {
                    G.cn[i] = 0;
                    G.gt[i] = 2;
                }
            }
        }
        else
        {
            if (D.rdstats[i].sv_support() != DEL &&  D.var_depth[i] > 0.8 && D.var_depth[i] < 1.2)
            {
                G.cn[i] = 2;
                G.gt[i] = 0;
            }
            else if (D.rdstats[i].sv_support() == DEL)
            {
                G.read_flag = true;
                if (D.var_depth[i] > 0.25 && D.var_depth[i] < 0.65) // TEMPORARY
                {
                    G.cn[i] = 1;
                    G.gt[i] = 1;
                }
                else if (D.var_depth[i]<0.25)
                {
                    G.cn[i] = 0;
                    G.gt[i] = 2;
                }
            }
        }
        if (G.gt[i] >=0)
        {
            G.ns += 1;
            G.ac += G.gt[i];
        }
    }
    if (G.dp_flag || G.dp2_flag || G.read_flag)
        G.b_pass = true;
    
}

void Genotyper::call_cnv(sv &S, SvData& D, SvGeno &G)
{
    // Fit Gaussian mixture models with 1, 2, and 3 components, compare BIC
    GaussianMixture gmix;
    
    std::vector< std::vector<double> > means = { {1.0}, {1.0, 1.5}, {1.0, 1.5, 2.0}, {1.0, 1.5, 2.0, 2.5}, {1.0, 1.5, 2.0, 2.5, 3.0}, {1.0, 1.5, 2.0, 2.5, 3.0, 3.5}, {1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0}, {1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5}, {1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0}};
    
    select_model(gmix, means, D.var_depth);
    
    if (gmix.n_comp > 1)
    {
        // success
        G.dp_flag = true;
        // assign genotypes

        for(int i=0; i<n_sample; ++i)
            G.cn[i] = gmix.assign_copynumber(D.var_depth[i]);
  
    }
    
    if (D.dp2.size()>1)
    {
        // DP100 genotyping
        GaussianMixture2 gmix2;
        int dp2_idx = 0;
        if (D.dp2.size()==4)
            dp2_idx = 1;
        
        select_model(gmix2, means, D.dp2[dp2_idx], D.dp2[dp2_idx+1]);
        
        // 2-D genotyping
        if (gmix2.n_comp>1)
        {
            // success
            G.dp2_flag = true;
            //assign dp2 genotypes
            for(int i=0; i<n_sample; ++i)
                G.cn[i] = gmix2.assign_copynumber(D.dp2[dp2_idx][i], D.dp2[dp2_idx+1][i]);
        }
    }
    
    // Readpair genotyping
    for(int i=0; i<n_sample; ++i)
    {
        if (G.dp_flag || G.dp2_flag)
        {
            if (D.rdstats[i].sv_support() == DUP)
            {
                G.read_flag = true;
                if (D.var_depth[i] > 1.35 && D.var_depth[i] < 1.75 && G.gt[i]== -1) // TEMPORARY
                {
                    G.cn[i] = 3;
                }
                else if (D.var_depth[i]>1.85 && G.gt[i] == -1)
                {
                    G.cn[i] = round(D.var_depth[i] * 2.0);
                }
            }
        }
        else
        {
            if (D.rdstats[i].sv_support() != DUP &&  D.var_depth[i] > 0.8 && D.var_depth[i] < 1.2)
            {
                G.cn[i] = 2;
            }
            else if (D.rdstats[i].sv_support() == DUP)
            {
                G.read_flag = true;
                if (D.var_depth[i] > 1.35)
                    G.cn[i] = round(D.var_depth[i] * 2.0);
            }
        }
        if (G.gt[i] >=0)
        {
            G.ns += 1;
            G.ac += G.gt[i];
        }
    }
    
    int max_cn = 2;
    int min_cn = 2;
    for(int i=0; i<n_sample; ++i)
    {
        if (G.cn[i] > max_cn)
        {
            max_cn = G.cn[i];
        }
        if (G.cn[i] < min_cn)
        {
            min_cn = G.cn[i];
        }
    }
    
    if (min_cn >=2 && max_cn <= 4)
        G.b_biallelic = true;
    else
        G.b_biallelic = false;
    
    if (G.b_biallelic)
    {
        for(int i=0; i<n_sample; ++i)
        {
            if (G.cn[i] == 2)
                G.gt[i] = 0; // 0/0
            else if (G.cn[i] == 3)
                G.gt[i] = 1; // 0/1
            else if (G.cn[i] == 4)
                G.gt[i] = 0; // 1/1
        }
    }
    if (G.dp_flag || G.dp2_flag || G.read_flag)
        G.b_pass = true;
}

std::string Genotyper::print(sv &S, SvData &D, SvGeno &G)
{
    std::string ln = std::to_string(S.chrnum);
    ln += "\t" + std::to_string(S.pos) + "\t" + svTypeName(S.svtype) + "_" + std::to_string(S.chrnum) + ":" + std::to_string(S.pos) + "-" + std::to_string(S.end) + "\t.\t<" + svTypeName(S.svtype) + ">\t.\t";
    
    if (G.b_pass)
    {
        ln += "PASS\t";
    }
    else
    {
        ln += "FAIL\t";
    }
    
    ln += "SVTYPE=" + std::string(svTypeName(S.svtype)) + ";END=" + std::to_string(S.end) + ";SVLEN=" + std::to_string(S.len) + ";AC=" + std::to_string(G.ac) + ";NS=" + std::to_string(G.ns) + ";AF=";
    if (G.ns>0)
        ln+=std::to_string((double)G.ac/(double)(2.0*G.ns));
    else
        ln+="0";
    
    ln+= G.info;

    
    if (G.dp_flag)
        ln += ";DP";
    if (G.dp2_flag)
        ln += ";DP2";
    if (G.read_flag)
        ln += ";READ";
    
    ln += "\tGT:CN";
    
    for (int i=0; i<n_sample; ++i)
    {
        switch(G.gt[i])
        {
            case 0:
                ln += "\t0/0";
                break;
            case 1:
                ln += "\t0/1";
                break;
            case 2:
                ln += "\t1/1";
                break;
            default:
                ln += "\t.";
                break;
        }
        if (G.cn[i]<0)
            ln += ":.";
        else
            ln += ":" + std::to_string(G.cn[i]);

    }
    ln += "\n";
    return ln;
}
