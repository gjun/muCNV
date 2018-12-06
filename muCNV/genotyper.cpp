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
    rp_flag = false;
    sp_flag = false;
    info = "";
}

SvData::SvData(int n)
{
    n_sample = n;
    rdstats.resize(n);
    sampstats.resize(n);
    var_depth.resize(n);
    dp_x.resize(n);
    dp_y.resize(n);
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
            ret_gmix2 = gmix2; // assignment
        }
    }
    return;
}

void Genotyper::call_deletion(sv &S, SvData &D, SvGeno &G)
{
    int n_sample = D.n_sample;
    
	// Fit Gaussian mixture models with 1, 2, and 3 components, compare BIC
    GaussianMixture gmix;
    std::vector< std::vector<double> > means = { {1.0}, {1.0, 0.5}, {0.5, 0.0}, {1.0, 0.5, 0.0}};
    select_model(gmix, means, D.var_depth);
    
    if (gmix.n_comp > 1)
    {
        // success
        // assign dp-genotypes
    }

    if (S.len > 300) // this should be dp_x and dp_y valid or not
    {
        // DP100 genotyping
        GaussianMixture2 gmix2;
        select_model(gmix2, means, D.dp_x, D.dp_y);
        
        // 2-D genotyping
        if (gmix2.n_comp>1)
        {
            // success
            //assign dp2 genotypes
        }
    }
    
    // Readpair genotyping
    
	std::vector<int> dp_gt(n_sample, -1);
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
        // assign genotypes
    }
    
    if (S.len > 300)
    {
        // DP100 genotyping
        GaussianMixture2 gmix2;
        select_model(gmix2, means, D.dp_x, D.dp_y);
        
        // 2-D genotyping
        if (gmix2.n_comp>1)
        {
            // success
            //assign dp2 genotypes
        }
    }
    
    // Readpair genotyping

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
    if (G.rp_flag)
        ln += ";RP";
    if (G.sp_flag)
        ln += ";SP";
    
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
        ln += ":" + std::to_string(G.cn[i]);
    }
    ln += "\n";
    return ln;
}
