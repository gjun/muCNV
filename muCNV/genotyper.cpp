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
    pd_flag = false;
    read_flag = false;
    dp_pre_mean = 1.0;
    dp_pre_std = 0.1;
    dp_post_mean = 1.0;
    dp_post_std = 0.1;
    b_pre = false;
    b_post = false;
    info = "";
}

SvData::SvData(int n)
{
    n_sample = n;
    rdstats.resize(n);
    prepost_dp.resize(n_sample, 1);
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

void Genotyper::get_prepost_stat(SvData &D, SvGeno &G)
{
    // depth context genotyping
    double dp_pre_sum = 0;
    double dp_post_sum = 0;
    double dp_pre_sumsq = 0;
    double dp_post_sumsq = 0;
    int dp_cnt = 0;

    for(int i=0; i<n_sample; ++i)
    {
        dp_pre_sum += D.dps[0][i];
        dp_post_sum += D.dps[1][i];
        dp_pre_sumsq += D.dps[0][i] * D.dps[0][i];
        dp_post_sumsq += D.dps[1][i] * D.dps[1][i];
        dp_cnt ++;
    }

    G.dp_pre_mean = dp_pre_sum/dp_cnt;
    G.dp_post_mean = dp_post_sum/dp_cnt;
    G.dp_pre_std = sqrt(dp_pre_sumsq / dp_cnt - G.dp_pre_mean * G.dp_pre_mean);
    G.dp_post_std = sqrt(dp_post_sumsq / dp_cnt - G.dp_post_mean * G.dp_post_mean);

    if (G.dp_pre_mean > 0.8 && G.dp_pre_mean < 1.2 && G.dp_pre_std < 0.2 )
    {
        G.b_pre = true;
    }
    if (G.dp_post_mean > 0.8 && G.dp_post_mean < 1.2 && G.dp_post_std < 0.2 )
    {
        G.b_post = true;
    }

    for(int i=0; i<G.n_sample; ++i)
    {
        double sum = 0;
        int cnt = 0;
        if (G.b_pre)
        {
            sum += D.dps[0][i];
            cnt ++;
        }
        if (G.b_post)
        {
            sum += D.dps[1][i];
            cnt ++;
        }

        D.prepost_dp[i] = (double) sum / cnt;
        if (G.b_pre && abs(D.prepost_dp[i] - G.dp_pre_mean)  > G.dp_pre_std * 2.0 )
            D.prepost_dp[i] = -1;
        if (G.b_post && abs(D.prepost_dp[i] - G.dp_post_mean)  > G.dp_post_std * 2.0 )
            D.prepost_dp[i] = -1;
    }
}

void Genotyper::select_model(GaussianMixture &ret_gmix, std::vector< std::vector<double> > &means, std::vector<double> &x)
{
    double best_bic = DBL_MAX;

    // number of models
    for(int m=0; m<(int)means.size(); ++m)
    {
        std::vector<double> s (means[m].size(), 0.1);
        GaussianMixture gmix(means[m], s);
        if (b_kmeans)
        {
            gmix.KM(x, b_mahalanobis);
        }
        else
        {
            gmix.EM(x); // fit mixture model
        }
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
    for(int m=0; m<(int)means.size(); ++m)
    {
        std::vector<double> s (means[m].size(), 0.01);
        GaussianMixture2 gmix2(means[m], s);
        if (b_kmeans)
        {
            gmix2.KM2(x, y, b_mahalanobis);
        }
        else
        {
            gmix2.EM2(x, y); // fit mixture model
        }
        if (gmix2.bic < best_bic && gmix2.p_overlap < MAX_P_OVERLAP )
        {
            best_bic = gmix2.bic;
            ret_gmix2 = gmix2; // assignment (copy operation)
        }
    }
    return;
}

void Genotyper::call(sv &S, SvData &D, SvGeno &G, double p, bool bk, bool bm)
{
    n_sample = D.n_sample;
    MAX_P_OVERLAP = p;

    b_kmeans = bk;
    b_mahalanobis = bm;

    get_prepost_stat(D, G);

    if (S.svtype == DEL)
    {
        call_deletion(S, D, G);
    }
    else if (S.svtype == DUP || S.svtype == CNV)
    {
        // TODO: Temporary. P_overlap does not work well for dups
        MAX_P_OVERLAP = 2.5*p;
        call_cnv(S, D, G);
    }
    else if (S.svtype == INV)
    {
        call_inversion(S, D, G);
    }
    /*
    else if (S.svtype == INS)
    {
        call_insertion(S, D, G);
    }
    */  // TODO: insertion calling
}

void Genotyper::call_inversion(sv &S, SvData &D, SvGeno &G)
{
    G.b_biallelic = true;
    for(int i=0; i<n_sample; ++i)
    {
        if (D.rdstats[i].inv_support())
        {
            G.read_flag = true;
            if (D.rdstats[i].n_pre_FF + D.rdstats[i].n_post_RR > 20) // todo: arbitrary, maybe clustering?
                G.gt[i] = 2;
            else
                G.gt[i] = 1;
            G.ac ++;
            G.ns ++;
        }
        else if (D.rdstats[i].n_pre_FF == 0 && D.rdstats[i].n_post_RR == 0)
        {
            G.gt[i] = 0;
            G.ns ++;
        }
    }

    double callrate = (double)G.ns / n_sample;

    if (callrate>0.5 && G.ac > 0 && G.ac < (G.ns*2))
        G.b_pass = true;
}

void Genotyper::call_insertion(sv &S, SvData &D, SvGeno &G)
{
    G.b_biallelic = true;
    for(int i=0; i<n_sample; ++i)
    {
        if (D.rdstats[i].ins_support())
        {
            G.read_flag = true;
            if (D.rdstats[i].n_pre_INS > 5 && D.rdstats[i].n_pre_INS + D.rdstats[i].n_post_clip_in + D.rdstats[i].n_pre_clip_in > 20) // todo: arbitrary, maybe clustering?
                G.gt[i] = 2;
            else
                G.gt[i] = 1;
            G.ac ++;
            G.ns ++;
        }
        else if (D.rdstats[i].n_pre_INS == 0)
        {
            G.gt[i] = 0;
            G.ns ++;
        }
    }

    double callrate = (double)G.ns / n_sample;

    if (callrate>0.5 && G.ac > 0 && G.ac < (G.ns*2))
        G.b_pass = true;
}

void Genotyper::call_deletion(sv &S, SvData &D, SvGeno &G)
{
    // fit gaussian mixture models with 1, 2, and 3 components, compare bic
    std::vector< std::vector<double> > means = { {1.0}, {1.0, 0.5}, {1.0, 0.5, 0.0}};

    G.b_biallelic = true;

    double best_BIC = DBL_MAX;
    double best_dp_idx = 2;

    // Genotyping based on 'variation' of depth: --__--

/*
    if (D.dps.size()>3)
    {
        for(int i=2; i<5; ++i)
        {
            GaussianMixture gm;
            select_model(gm, means, D.dps[i]);
            if (gm.bic < best_BIC && gm.p_overlap < MAX_P_OVERLAP)
            {
                best_BIC = gm.bic;
                G.gmix = gm;
                best_dp_idx = i;
            }
        }
    }

*/
    select_model(G.gmix, means, D.dps[best_dp_idx]);

    std::vector<double> &var_depth = D.dps[best_dp_idx];

    // depth clustering
    if (G.gmix.n_comp > 1 && G.gmix.ordered())
    {
        G.dp_flag = true;
        // success
        // assign dp-genotypes
        for(int i=0; i<n_sample; ++i)
        {
            int cn = G.gmix.assign_copynumber(var_depth[i]);

            if (cn == 2)
            {
                G.cn[i] = 2;
                G.gt[i] = 0; // 0/0
            }
            else if (cn == 1)
            {
                G.cn[i] = 1;
                G.gt[i] = 1; // 0/1

            }
            else if (cn == 0)
            {
                G.cn[i] = 0;
                G.gt[i] = 2; // 1/1
            }
        }
    }

    int dp2_idx = 3;
    if (D.dps.size()>3) //dp2 has more than 2 vectors
    {
        // dp100 genotyping
        select_model(G.gmix2, means, D.dps[dp2_idx], D.dps[dp2_idx+1]);

        // 2-d genotyping
        if (G.gmix2.n_comp>1 && G.gmix2.ordered())
        {
            G.dp2_flag = true;
            //assign dp2 genotypes

            for(int i=0; i<n_sample; ++i)
            {
                int cn = G.gmix2.assign_copynumber(D.dps[dp2_idx][i], D.dps[dp2_idx+1][i]);
                if (cn == 2)
                {
                    G.cn[i] = 2;
                    G.gt[i] = 0; // 0/0
                }
                else if (cn == 1)
                {
                    G.cn[i] = 1;
                    G.gt[i] = 1; // 0/1
                }
                else if (cn == 0)
                {
                    G.cn[i] = 0;
                    G.gt[i] = 2; // 1/1
                }
            }
        }
    }

    // Let's use peripheral depth to identify 'false positve' only

    if ((G.b_pre && G.b_post) && (G.dp_flag || G.dp2_flag))
    {
        // If clustered, filter false positive variants
        for (int i=0;i<n_sample; ++i)
        {
            if (G.gt[i] != 0) // missing or variant
            {
                if (D.prepost_dp[i] - var_depth[i] < 0.2) // if there's no 'dip', 0/0
                {
                    if (var_depth[i] > 0.8) // TODO: arbitrary
                    {
                        G.cn[i] = 2;
                        G.gt[i] = 0;
                    }
                    else
                    {
                        G.cn[i] = -1;
                        G.gt[i] = -1;
                    }
                }
                else if (D.prepost_dp[i] > 0 && D.prepost_dp[i] - var_depth[i] > 0.3)
                {
                    if (var_depth[i] >=0.15 && var_depth[i] < 0.6) // if there's a dip and depth is low
                    {
                        G.pd_flag = true;
                        G.cn[i] = 1;
                        G.gt[i] = 1;
                    }
                    else if (var_depth[i] < 0.15) // if there's a dip and depth is very low
                    {
                        G.pd_flag = true;
                        G.cn[i] = 0;
                        G.gt[i] = 2;
                    }
                    else // otherwise, missing
                    {
                        G.cn[i] = -1;
                        G.gt[i] = -1;
                    }
                }
                else // otherwise, missing
                {
                    G.cn[i] = -1;
                    G.gt[i] = -1;
                }
            }
        }
    }
    

   // readpair genotyping
    if (S.len >= 50)
    {
        for(int i=0; i<n_sample; ++i)
        {
            if (G.dp_flag || G.dp2_flag || G.pd_flag)
            {
                if (D.rdstats[i].del_support())
                {
                    if (var_depth[i] > 0.15 && var_depth[i] < 0.7) // TODO: This is arbitrary - cluster read pair numbers to get genotypes
                    {
                        G.read_flag = true;
                        G.cn[i] = 1;
                        G.gt[i] = 1;
                    }
                    else if (var_depth[i]<=0.15)
                    {
                        G.read_flag = true;
                        G.cn[i] = 0;
                        G.gt[i] = 2;
                    }
                    else if (var_depth[i]>0.8)
                    {
                        G.cn[i] = -1;
                        G.gt[i] = -1;
                    }
                }
                else if (var_depth[i] > 0.8 && var_depth[i] < 1.2)
                {
                    G.cn[i] = 2;
                    G.gt[i] = 0;
                }
            }
            else
            {
                if (D.rdstats[i].del_support())
                {
                    if (var_depth[i] > 0.2 && var_depth[i] < 0.7) // TEMPORARY
                    {
                        G.read_flag = true;
                        G.cn[i] = 1;
                        G.gt[i] = 1;
                    }
                    else if (var_depth[i]<=0.2)
                    {
                        G.read_flag = true;
                        G.cn[i] = 0;
                        G.gt[i] = 2;
                    }
                }
                else if (var_depth[i] > 0.8 && var_depth[i] < 1.2)
                {
                    G.cn[i] = 2;
                    G.gt[i] = 0;
                }
            }

        }
    }
 

    int n_het = 0;
    int n_hom = 0;
    int n_ref = 0;

    double het_dp = 0;
    double hom_dp = 0;
    double ref_dp = 0;


    for(int i=0; i<n_sample; ++i)
    {
        if (G.gt[i] >=0)
        {
            G.ns += 1;
            G.ac += G.gt[i];

            if (G.gt[i] == 0)
            {
                n_ref ++;
                ref_dp += var_depth[i];

            }
            else if (G.gt[i] == 1)
            {
                n_het ++;
                het_dp += var_depth[i];
            }
            else if (G.gt[i] == 2)
            {
                n_hom ++;
                hom_dp += var_depth[i];
            }
        }
    }


    if (n_ref == 0)
    {
        return;
    }

    if (n_hom >0 )
        hom_dp /= (double)n_hom;

    if (n_het >0 )
        het_dp /= (double)n_het;

    if (n_ref >0 )
        ref_dp /= (double)n_ref;

    if (het_dp < 0.3 || hom_dp > 0.2)
    {
        return;
    }

    double callrate = (double)G.ns / n_sample;

    if ( ( ((G.dp_flag || G.dp2_flag) && (G.pd_flag || G.read_flag)) || G.read_flag ) &&  callrate>0.5 && G.ac > 0 && G.ac < G.ns*2)
        G.b_pass = true;

    // Excessive heterozygosity (basically all-het case)
    if ( n_het > ((double)G.ns * 0.75) && n_hom < ((double)0.05 * G.ns) )
        G.b_pass = false;
}

void Genotyper::call_cnv(sv &S, SvData& D, SvGeno &G)
{
    // Fit Gaussian mixture models with 1, 2, and 3 components, compare BIC

    std::vector< std::vector<double> > means = { {1.0}, {1.0, 1.5}, {1.0, 1.5, 2.0}, {1.0, 1.5, 2.0, 2.5}, {1.0, 1.5, 2.0, 2.5, 3.0}, {1.0, 1.5, 2.0, 2.5, 3.0, 3.5}, {1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0} };

    std::vector<int> dp_cn (n_sample, -1);
    std::vector<int> dp2_cn (n_sample, -1);

    int dp_ns = 0;
    int dp2_ns = 0;


    double best_BIC = DBL_MAX;
    double best_dp_idx = 2;

/*
    if (D.dps.size() > 3)
    {
        for(int i=3; i<5; ++i)
        {
            GaussianMixture gm;
            select_model(gm, means, D.dps[i]);
            if (gm.bic < best_BIC && gm.p_overlap < MAX_P_OVERLAP)
            {
                best_BIC = gm.bic;
                G.gmix = gm;
                best_dp_idx = i;
            }
        }
    }
    */
    select_model(G.gmix, means, D.dps[best_dp_idx]);

    std::vector<double> &var_depth = D.dps[best_dp_idx];

    if (G.gmix.n_comp > 1 && G.gmix.r_ordered() )
    {
        // success

        for(int i=0; i<(int)n_sample; ++i)
        {
            dp_cn[i] = G.gmix.assign_copynumber(var_depth[i]);
            if (dp_cn[i] >= 2)
            {
                dp_ns++;
            }

            // TODO: arbitrary
            /*
            if (var_depth[i] > 0.7 && var_depth[i] < 1.25)
            {
                dp_cn[i] = 2;
                dp_ns++;
            }
            else if (var_depth[i] > 1.4 &&  var_depth[i] < 1.8)
            {
                G.dp_flag = true;
                dp_cn[i] = 3;
                dp_ns ++;
            }
            else if (var_depth[i] >= 1.85)
            {
                G.dp_flag = true;
                dp_cn[i] = round(var_depth[i] * 2.0);
                dp_ns++;
            }
            */
        }
    }

/*
    if (D.dps.size()>3)
    {
        // DP100 genotyping
        int dp2_idx = 2;

        select_model(G.gmix2, means, D.dps[dp2_idx], D.dps[dp2_idx+1]);

        // 2-D genotyping
        if (G.gmix2.n_comp>1 && G.gmix2.r_ordered() )
        {
            //TODO : dp2 geotyping seems unstable, maybe make this to post-processing for variants with pd_flag && rd_flag
            for(int i=0; i<(int)n_sample; ++i)
            {
//                dp2_cn[i] = G.gmix2.assign_copynumber(D.dps[dp2_idx][i], D.dps[dp2_idx+1][i]);

                double dp = (D.dps[dp2_idx][i] + D.dps[dp2_idx+1][i]) / 2.0;
                if (dp > 0.7 && dp < 1.25)
                {
                    dp2_cn[i] = 2;
                    dp2_ns++;
                }
                else if (dp > 1.4 &&  dp < 1.8)
                {
                     G.dp2_flag = true;
                    dp2_cn[i] = 3;
                    dp2_ns ++;
                }
                else if (dp >= 1.85)
                {
                    G.dp2_flag = true;
                    dp2_cn[i] = round(var_depth[i] * 2.0);
                    dp2_ns++;
                }
            }
        }
    }
    */

    if (dp2_ns > dp_ns)
    {
        for(int i=0; i<n_sample; ++i)
            G.cn[i] = dp2_cn[i];
    }
    else
    {
        for(int i=0; i<n_sample; ++i)
            G.cn[i] = dp_cn[i];
    }

    get_prepost_stat(D, G);

    if (G.b_pre && G.b_post  && (G.dp_flag || G.dp2_flag))
    {
        // If clustered, filter false positive variants
        for (int i=0;i<n_sample; ++i)
        {
            if (G.gt[i] != 0) // missing or variant
            {
                if (var_depth[i] - D.prepost_dp[i] < 0.2)
                {
                    // TODO: arbitrary
                    if (var_depth[i] < 1.2 && var_depth[i] > 0.8 )
                    {
                        G.cn[i] = 2;
                    }
                    else
                    {
                        G.cn[i] = -1;
                    }
                }
                else if (D.prepost_dp[i] > 0 && var_depth[i] - D.prepost_dp[i] > 0.3)
                {
                    // TODO: arbitrary
                    if (var_depth[i] >= 1.4 && var_depth[i] <= 1.85)
                    {
                        G.pd_flag = true;
                        G.cn[i] = 3;
                    }
                    else if (var_depth[i] > 1.85)
                    {
                        G.pd_flag = true;
                        G.cn[i] = round(var_depth[i] * 2.0);
                    }
                    else
                    {
                        G.cn[i] = -1;
                    }
                }
                else
                {
                    G.cn[i] = -1;
                }
            }
        }
    }
    // Readpair genotyping
    for(int i=0; i<n_sample; ++i)
    {
        if (G.dp_flag || G.dp2_flag || G.pd_flag)
        {
            if (D.rdstats[i].dup_support())
            {
                if (var_depth[i] > 1.3 && var_depth[i] <= 1.85) // TEMPORARY
                {
                    G.read_flag = true;
                    G.cn[i] = 3;
                }
                else if (var_depth[i]>1.85)
                {
                    G.read_flag = true;
                    G.cn[i] = round(var_depth[i] * 2.0);
                }
            }
            else if (var_depth[i] > 0.8 && var_depth[i] <= 1.3 && G.cn[i] < 0)
            {
                G.cn[i] = 2.0;
            }
        }
        else
        {
            if (!D.rdstats[i].dup_support() &&  var_depth[i] > 0.75 && var_depth[i] <= 1.3)
            {
                G.cn[i] = 2;
            }
            else if (D.rdstats[i].dup_support())
            {
                if (var_depth[i] >= 1.35)
                {
                    G.read_flag = true;
                    G.cn[i] = round(var_depth[i] * 2.0);
                }
            }
        }
        if (G.cn[i] >=0)
        {
            G.ns += 1;
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
    {
        G.b_biallelic = true;
    }
    else
    {
        G.b_biallelic = false;
    }

    for(int i=0; i<n_sample; ++i)
    {
        if (G.cn[i] == 2)
        {
            G.gt[i] = 0; // 0/0
        }
        else if (G.cn[i] == 3)
        {
            G.gt[i] = 1; // 0/1
            G.ac++;
        }
        else if (G.cn[i] >3)
        {
            G.gt[i] = 2; // 1/1
            G.ac+=2;
        }
    }

    double callrate = (double)G.ns / n_sample;

    if ((G.dp_flag || G.dp2_flag || G.read_flag ) && callrate>0.5 && G.ac>0)
        G.b_pass = true;
}
