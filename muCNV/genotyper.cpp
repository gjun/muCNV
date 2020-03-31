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
#include <unordered_map>

BreakCluster::BreakCluster()
{
    N = 0;
    start_mean = 0;
    start_var = 0;
    end_mean = 0;
    end_var = 0;
}

double BreakCluster::get_distance(std::pair<int, int> &span)
{
    double dx = (span.first - start_mean);
    double dy = (span.second - end_mean);
    
    return abs(dx) + abs(dy);
}

double BreakCluster::get_distance(int first, int second)
{
    double dx = (first - start_mean);
    double dy = (second - end_mean);
    
    return abs(dx) + abs(dy);
}

double BreakCluster::get_distance(BreakCluster& br)
{
    double dx = (br.start_mean - start_mean);
    double dy = (br.end_mean - end_mean);
    
    return abs(dx) + abs(dy);
}


void BreakCluster::add_to_cluster(std::pair<int, int> &span)
{
    if (N == 0)
    {
        start_mean = span.first;
        end_mean = span.second;
        start_var = 0;
        end_var = 0;
        N = 1;
    }
    else
    {
        double dx = span.first - start_mean;
        double dy = span.second - end_mean;
        start_var += dx*dx /(double)(N+1) - start_var/(double)N;
        start_mean += dx / (double)(N+1);
        
        end_var += dy*dy/(double)(N+1) - end_var/(double)N;
        end_mean += dy / (double)(N+1);
        N = N + 1;
    }
}

void BreakCluster::add_to_cluster(int first, int second, int k)
{
    if (N == 0)
    {
        start_mean = first;
        end_mean = second;
        start_var = 0;
        end_var = 0;
        N = k;
    }
    else
    {
        double dx = first - start_mean;
        double dy = second - end_mean;
        start_var += dx*dx /(double)(N+1) - start_var/(double)N; // TODO: Fix this later, variance is not correct now but it's not being used.
        start_mean += (k * dx) / (double)(N+k);
        
        end_var += dy*dy/(double)(N+1) - end_var/(double)N;
        end_mean += ( k * dy ) / (double)(N+k);
        N = N + k;
    }
}


void BreakCluster::merge(BreakCluster& br)
{
    start_mean = (start_mean * N + br.start_mean * br.N) / (double)(N + br.N);
    end_mean = (end_mean * N + br.end_mean * br.N) / (double)(N + br.N);
    
    start_var = ((N - 1) * start_var + (br.N -1) * br.start_var + (N * br.N) / (N + br.N) * (start_mean - br.start_mean) * (start_mean - br.start_mean)) / (double)(N + br.N - 1.0);
    end_var = ((N - 1) * end_var + (br.N -1) * br.end_var + (N * br.N) / (N + br.N) * (end_mean - br.end_mean) * (end_mean - br.end_mean)) / (double)(N + br.N - 1.0);
    N += br.N;
}


SvGeno::SvGeno(int n)
{
    n_sample = n;
    n_effect = n_sample;
    
    MAX_P_OVERLAP = 1.0;

    gt.resize(n_sample, -1);
    cn.resize(n_sample, -1);

    start_clips.resize(n_sample, 0);
    end_clips.resize(n_sample, 0);
    
    split_cnts.resize(n_sample, 0);
    rp_cnts.resize(n_sample, 0);
    
    all_cnts.resize(n_sample, 0);

    start_rps.resize(n_sample, 0);
    end_rps.resize(n_sample, 0);
    
    nonref_mask.resize(n_sample, false);
    sample_mask.resize(n_sample, true);

    ns = 0;
    ac = 0;

    b_biallelic = false;
    b_pass = false;
    
    cnt_flag = false;
    dp_flag = false;
    dp2_flag = false;
    dpcnt_flag = false;
    pd_flag = false;
    readpair_flag = false;
    clip_flag = false;
    split_flag = false;
    
    dp_pre_mean = 1.0;
    dp_pre_std = 0.1;
    dp_post_mean = 1.0;
    dp_post_std = 0.1;
    b_pre = false;
    b_post = false;
    info = "";
}

void SvGeno::reset()
{
    std::fill(gt.begin(), gt.end(), -1);
    std::fill(cn.begin(), cn.end(), -1);
    
    std::fill(split_cnts.begin(), split_cnts.end(), 0);
    std::fill(rp_cnts.begin(), rp_cnts.end(), 0);    
    std::fill(nonref_mask.begin(), nonref_mask.end(), false);
    
    std::fill(all_cnts.begin(), all_cnts.end(), 0);

    std::fill(start_clips.begin(), start_clips.end(), 0);
    std::fill(end_clips.begin(), end_clips.end(), 0);
    std::fill(start_rps.begin(), start_rps.end(), 0);
    std::fill(end_rps.begin(), end_rps.end(), 0);
    
    ns = 0;
    ac = 0;

    b_biallelic = false;
    b_pass = false;
    cnt_flag = false;
    dp_flag = false;
    dp2_flag = false;
    dpcnt_flag = false;
    pd_flag = false;
    readpair_flag = false;
    clip_flag = false;
    split_flag = false;
    
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
    raw_dp.resize(n_sample);

    dps.resize(5);
    
    
    multi_dp = false;
    
    all_rps.resize(4);
    
    for(int i=0; i<4; ++i)
    {
        all_rps[i].resize(200, 0);
    }
    
    all_sps.resize(200, 0);
    
    all_lclips.resize(400, 0);
    all_rclips.resize(400, 0);
    
    clus_idx = -1;
    
    for(int j=0;j<5;++j)
    {
        dps[j].resize(n_sample, 0);
        // 0 : pre-depth
        // 1 : post-depth
        // 2 : 1-D depth (var_depth)
        // 3 : avg depth in some place of first half
        // 4 : avg depth in some place of second half
    }
}

void SvData::reset()
{
    for(int i=0;i<n_sample;++i)
        rdstats[i].reset();
    multi_dp = false;
    
    for(int i=0; i<4; ++i)
        std::fill(all_rps[i].begin(), all_rps[i].end(), 0);
    
    std::fill(all_sps.begin(), all_sps.end(), 0);
    std::fill(all_lclips.begin(), all_lclips.end(), 0);
    std::fill(all_rclips.begin(), all_rclips.end(), 0);

    for(int j=0;j<5;++j)
    {
        std::fill(dps[j].begin(), dps[j].end(), 0);
    }
    
    clus_idx = -1;
    
    std::fill(prepost_dp.begin(), prepost_dp.end(), 1);
    std::fill(raw_dp.begin(), raw_dp.end(), 0);
    vec_break_clusters.clear();
}

std::string SvStat::get_summary_stat(SvGeno &G, SvData &D)
{
    std::string outstr = "";
    depth_mix.estimate(D.dps[2], G.gt, 3);
    outstr += depth_mix.print_str();

    split_mix.estimate(G.split_cnts, G.gt, 3);
    outstr += ";SPSTAT=" + split_mix.print_str();

    rp_mix.estimate(G.rp_cnts, G.gt, 3);    
    outstr += ";RPSTAT=" + rp_mix.print_str();

    lclip_mix.estimate(G.start_clips, G.gt, 3);
    outstr += ";LCSTAT=" + lclip_mix.print_str();

    rclip_mix.estimate(G.end_clips, G.gt, 3);
    outstr += ";RCSTAT=" + rclip_mix.print_str();

    dpcnt_mix.estimate(D.dps[2], G.all_cnts, G.gt, 3);
    outstr += ";DPCNTSTAT=" + dpcnt_mix.print_str();

    return outstr;
    
};


void copyComps(std::vector<Gaussian> &C, std::vector<Gaussian> &C0)
{
    C.clear();
    C.resize(C0.size());
    for(unsigned j=0;j<(unsigned)C.size(); ++j)
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
        if (!G.sample_mask[i]) continue;

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
        if (!G.sample_mask[i]) continue;

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

void Genotyper::select_model_1d(GaussianMixture &ret_gmix, std::vector< std::vector<double> > &means, std::vector<double> &x, double MAX_P_OVERLAP)
{
    double best_bic = DBL_MAX;

    // number of models
    for(int m=0; m<(int)means.size(); ++m)
    {
        std::vector<double> s (means[m].size(), 0.1);
        GaussianMixture gmix(means[m], s);

        gmix.EM(x); // fit mixture model
        // gmix.print(stderr);
        if (gmix.bic < best_bic && gmix.p_overlap < MAX_P_OVERLAP)
        {
            best_bic = gmix.bic;
            ret_gmix = gmix; // assignment
        }
    }
    return;
}

void Genotyper::select_model_mask_1d(GaussianMixture &ret_gmix, std::vector< std::vector<double> > &means, std::vector<double> &x, std::vector<bool> &mask, double MAX_P_OVERLAP)
{
    double best_bic = DBL_MAX;

    // number of models
    for(int m=0; m<(int)means.size(); ++m)
    {
        std::vector<double> s (means[m].size(), 0.1);
        GaussianMixture gmix(means[m], s);
        gmix.EM_select(x, mask); // fit mixture model
        // gmix.print(stderr);

        if (gmix.bic < best_bic && gmix.p_overlap < MAX_P_OVERLAP)
        {
            best_bic = gmix.bic;
            ret_gmix = gmix; // assignment
        }
    }
    return;
}

void Genotyper::select_model_mask_2d(GaussianMixture2 &ret_gmix2, std::vector< std::vector<double> > &means, std::vector<double> &x, std::vector<double> &y, std::vector<bool> &mask, double MAX_P_OVERLAP)
{
    double best_bic = DBL_MAX;

    // number of models
    for(int m=0; m<(int)means.size(); ++m)
    {
        std::vector<double> s (means[m].size(), 0.1);
        GaussianMixture2 gmix2(means[m], s);
        gmix2.EM2_select(x, y, mask); // fit mixture model
        
        if (gmix2.bic < best_bic && gmix2.p_overlap < MAX_P_OVERLAP )
        {
            best_bic = gmix2.bic;
            ret_gmix2 = gmix2; // assignment (copy operation)
        }
    }
    return;
}

void Genotyper::select_model_2d(GaussianMixture2 &ret_gmix2, std::vector< std::vector<double> > &means, std::vector<double> &x, std::vector<double> &y, double MAX_P_OVERLAP)
{
    double best_bic = DBL_MAX;

    // number of models
    for(int m=0; m<(int)means.size(); ++m)
    {
        std::vector<double> s (means[m].size(), 0.1);
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

void Genotyper::select_model_dpcnt_mask(GaussianMixture2 &ret_gmix2, std::vector<double> &dps, std::vector<double> &cnts, std::vector<int> &lbl, int max_n_comp,
            std::vector<bool> &mask, double MAX_P_OVERLAP)
{
    std::vector<double> m = {1.0};
    std::vector<double> s = {0.1};
    GaussianMixture2 tmp_mix(m, s);

    DDMSG("DPCNT model with 1 component");
    tmp_mix.EM2_select(dps, cnts, mask);
    ret_gmix2 = tmp_mix;

    for(int n_comp=2; n_comp<=max_n_comp; ++n_comp)
    {
        DDPRINT("DPCNT model with %d components\n", n_comp);
        tmp_mix.estimate_select(dps, cnts, lbl, mask, n_comp);
    
        if (tmp_mix.bic < ret_gmix2.bic && tmp_mix.p_overlap < MAX_P_OVERLAP)
        {
            ret_gmix2 = tmp_mix;
        }
    }
    return;
}

void Genotyper::call(sv &S, SvData &D, SvGeno &G, std::vector<SampleStat> &stats)
{
    n_sample = G.n_sample;

    get_prepost_stat(D, G);

	// if (!G.b_pre || !G.b_post)
	//	return;
#ifdef DDEBUG
    S.print(stderr);
#endif
    DDMSG("");
    if (S.svtype == DEL)
	{
       call_deletion(S, D, G, stats);
        // HWE: chisq = G.ns *  (( 4*n_AA*n_aa - n_Aa^2) / ((2*n_AA + n_Aa) * (2*n_aa + n_Aa)) )^2 

    }
    else if (S.svtype == DUP || S.svtype == CNV)
    {
       call_cnv(S, D, G, stats);
       // How to check HWE for CNVs?
       // Should we try to identify CN alleles by taking GCD of CNs 
       // Fraction_biallelic ? (Most samples are biallelic except for a few?)
    }
    else if (S.svtype == INV)
    {
        call_inversion(S, D, G, stats);
         // HWE: chisq = G.ns *  (( 4*n_AA - n_Aa^2) / ((2*n_AA + n_Aa) * (2n_aa + n_Aa)) )^2 

    }
}

bool Genotyper::assign_inv_genotypes(sv &S, SvData &D, SvGeno &G, std::vector<SampleStat> &stats)
{
    std::vector<double> norm_cnts (n_sample, 0);
    GaussianMixture tmp_mix;

    double tmp_sum = 0;
    double tmp_cnt = 0;
    double cnt_max = -1;
    double cnt_mean = 0;
    double best_bic = DBL_MAX;

    DDMSG("assigning INV genotypes");

    for(int i=0; i<n_sample; ++i)
    {
        if (G.sample_mask[i])
        {
            norm_cnts[i] = G.all_cnts[i] / stats[i].avg_dp;

            if (G.gt[i] > 0)
            {
                tmp_sum += norm_cnts[i];
                tmp_cnt += 1;
            } 
        }
    }
    if (tmp_cnt > 0)
    {
        cnt_mean = tmp_sum / tmp_cnt;
    }
    else
    {
        return false;
    }

    G.gmix.Comps.resize(1);
    G.gmix.n_comp = 1;
    G.gmix.Comps[0].estimate_select(norm_cnts, G.sample_mask);
    G.gmix.Comps[0].Alpha = 1.0;
    G.gmix.updateAICBIC_select(norm_cnts, G.sample_mask);
    best_bic = G.gmix.bic;
    double MAX_P_OVERLAP = 0.3; // NOTE: ARBITRARY

    tmp_mix.estimate(norm_cnts, G.gt, 2);
    if (tmp_mix.bic < best_bic && tmp_mix.p_overlap < MAX_P_OVERLAP && tmp_mix.Comps[0].Mean > 0.15)
    {
        G.gmix = tmp_mix;
        best_bic = tmp_mix.bic;
    }

    tmp_sum = 0;
    tmp_cnt = 0;
    for(int i=0; i<n_sample; ++i)
    {
        if (G.sample_mask[i])
        {
            if (G.gt[i] > 0 && norm_cnts[i] > cnt_mean * 1.5)
            {
                tmp_sum += norm_cnts[i];
                tmp_cnt += 1;
                G.gt[i] = 2;
            }
        }
    }
    if (tmp_cnt>0)
    {
        cnt_max = tmp_sum / tmp_cnt;

        tmp_mix.estimate(norm_cnts, G.gt, 3);
        if (tmp_mix.bic < best_bic && tmp_mix.p_overlap < MAX_P_OVERLAP)
        {
            G.gmix = tmp_mix;
            best_bic = tmp_mix.bic;
        }
    }

    DDPRINT("cnt_mean %f, cnt_max %f \n", cnt_mean, cnt_max);

    G.dp_flag = false;
    
    if (G.gmix.n_comp > 1)
    {     
        G.ns = 0;
        G.ac = 0;
        
        DDPRINT("%d components\n", genostat.dpcnt_mix.n_comp);
        for(int i=0; i<n_sample; ++i)
        {
            G.cn[i] = -1;
            G.gt[i] = -1;

            if (!G.sample_mask[i]) continue;

            int gt = G.gmix.assign_cluster(norm_cnts[i]);

            if (gt == 0)
            {
                if (norm_cnts[i] > 0.1)
                {
                    G.cn[i] = -1;
                    G.gt[i] = -1;
                }
                else
                {
                    G.cn[i] = 2;
                    G.gt[i] = 0; // 0/0
                    G.ns ++;
                }
            }
            else
            {
                if (gt == -1 && norm_cnts[i] < 0.05)
                {
                    G.gt[i] = 0;
                    G.ns ++;
                }
                else if (norm_cnts[i] < 0.05)
                {
                    G.cn[i] = -1;
                    G.gt[i] = -1;
                }
                else if (gt == 1)
                { 
                    G.gt[i] = 1; // 0/1
                    G.ac ++;
                    G.ns ++;
                }
                else if (gt == 2)
                { 
                    G.gt[i] = 2; // 1/1
                    G.ns ++;
                    G.ac += 2;
                }
            }
        
        }

        double callrate = (double)G.ns / G.n_effect;
        DDPRINT("NS %d, AC %d, Call rate %f\n", G.ns, G.ac, callrate);
        
        if (callrate>0.5 && G.ac > 0 && G.ac < G.ns*2) // if successful, return genotypes
        { 
            G.cnt_flag = true; 
            G.b_pass = true;
            DDMSG("INV CNT clustering successful");
            return true;
        }
        else
        {
            G.dp_flag = false;
            G.cnt_flag = false;
            return false;
        }
    }
    return false;    
}

// bool Genotyper::assign_inv_genotypes(sv &S, SvData &D, SvGeno &G, std::vector<SampleStat> &stats)
// {
//     GaussianMixture all_mix;
    
//     std::vector<double> norm_cnts (n_sample, 0);
//     for(int i=0; i<n_sample; ++i)
//     {
//         if (G.sample_mask[i])
//             norm_cnts[i] = G.all_cnts[i] / stats[i].avg_dp;
//     }
    
//     all_mix.estimate(norm_cnts, G.gt, 2);

//     double err_bound = all_mix.Comps[0].Stdev;
//     if (err_bound < 0.15) err_bound = 0.15;
    
//     // TODO: Arbitrary cutoffs
//     if (all_mix.Comps[0].Mean>0.1 || all_mix.Comps[1].Mean - err_bound < all_mix.Comps[0].Mean || all_mix.p_overlap > 0.5)
//         return false;

//    std::vector<std::vector<double> > alt_means = { {all_mix.Comps[1].Mean}, {all_mix.Comps[1].Mean, all_mix.Comps[1].Mean * 2.0} };
//    select_model_mask_1d(G.gmix, alt_means, norm_cnts, G.nonref_mask, 0.3);
//     //   G.gmix.print(stdout);
    
//     G.ns = 0;
//     G.ac = 0;
    
//     for(int i=0; i<n_sample; ++i)
//     {
//         if (!G.sample_mask[i])
//             continue;

//         double d0 = norm_cnts[i] - all_mix.Comps[0].Mean ;
//         double d1 = G.gmix.Comps[0].Mean - norm_cnts[i];
        
//         G.gt[i] = -1;
        
//         //if (d0>0 &&  G.all_cnts[i]>0 && d1 < d0 * ((G.split_cnts[i] + G.rp_cnts[i] + G.all_cnts[i] - 1) / 2.0))
//         if (d0>0 &&  G.all_cnts[i]>=5 && d1 < d0 && (G.split_cnts[i] + G.rp_cnts[i] > 1))
//         {
//             if (G.gmix.n_comp > 1  && (G.gmix.Comps[1].Mean - norm_cnts[i] < (norm_cnts[i] - G.gmix.Comps[0].Mean) * 0.5))
//             {
//                 G.gt[i] = 2;
//                 G.ns++;
//                 G.ac += 2;
//             }
//             else
//             {
//                 G.gt[i] = 1;
//                 G.ns++;
//                 G.ac++;
//             }
//         }
//         else
//         {
//             if (norm_cnts[i] > all_mix.Comps[0].Mean + (err_bound / (G.all_cnts[i]+1.0) ) )  // missing
//                 G.gt[i] = -1;
//             else if (d0 < d1)
//             {
//                 G.gt[i] = 0;
//                 G.ns++;
//             }
//             else
//             {
//                 G.gt[i] = -1;
//             }
//         }
//     }
// /*
//     for(int i=0; i<n_sample;++i)
//     {
//         printf("%d\t%d\t%d\t%d\t%d\t%f\t%f\t%f\t%d\n", (int)G.split_cnts[i], (int)G.rp_cnts[i], (int)G.start_clips[i], (int)G.end_clips[i], (int)G.all_cnts[i], G.rp_cnts[i] / stats[i].avg_dp, G.start_clips[i] / stats[i].avg_dp, G.end_clips[i] / stats[i].avg_dp , G.gt[i]);
//     }
//   */
//     double callrate = (double)G.ns / G.n_effect;;

//     if (callrate>0.5 && G.ac > 0 && G.ac < G.ns*2) // if successful, return genotypes
//     {
//         G.b_pass = true;
//         G.cnt_flag = true;
//         G.gmix.n_comp = G.gmix.n_comp + 1;
//         G.gmix.Comps.resize(G.gmix.n_comp);
//         G.gmix.Comps[G.gmix.n_comp-1].Mean = all_mix.Comps[0].Mean;
//         G.gmix.Comps[G.gmix.n_comp-1].Stdev = all_mix.Comps[0].Stdev;
//         G.gmix.Comps[G.gmix.n_comp-1].Alpha = all_mix.Comps[0].Alpha;
        
//         return true;
//     }
//     return false;
// }
        
bool Genotyper::get_inv_cnts(sv &S, SvData &D, SvGeno &G)
{
    G.startclip_sum = 0;
    G.endclip_sum = 0;
    G.rp_sum = 0;
    G.split_sum = 0;
    G.ns = 0;
    G.ac = 0;
    
    int start_pos  =(int)(D.vec_break_clusters[D.clus_idx].start_mean+0.5);
    int end_pos = (int)(D.vec_break_clusters[D.clus_idx].end_mean+0.5);
    
    int n_clus = (int)D. vec_break_clusters.size();
    
    std::vector<int> start_idx (n_clus, 0);
    std::vector<int> end_idx (n_clus, 0);
    
    for(int j=0; j<n_clus; ++j)
    {
        start_idx[j] = ((int)(D.vec_break_clusters[j].start_mean+0.5) - S.pos + 100);
        end_idx[j] = ((int)(D.vec_break_clusters[j].end_mean+0.5) - S.end + 300);
        
        if (start_idx[j] < 3) start_idx[j] = 3;
        if (start_idx[j] > 196) start_idx[j] = 196;
        if (end_idx[j] < 203) end_idx[j] = 203;
        if (end_idx[j] > 396) end_idx[j] = 396;
    }
    
    for(int i=0; i< n_sample; ++i)
    {
        G.nonref_mask[i] = false;
        G.gt[i] = -1;
        G.cn[i] = -1;
        if (!G.sample_mask[i])
            continue;

        G.split_cnts[i] = 0;
        // Count the number of split reads
        for(auto &split : D.rdstats[i].splits)
        {
            for (int j=0; j<n_clus; ++j)
            {
                if (D.vec_break_clusters[j].get_distance(split.positions) <=6)
                {
                    G.split_cnts[i] ++;
                }
            }
        }
        
        for(int j=0; j<n_clus; ++j)
        {
            // Count the number of soft clips around breakpoints
            G.start_clips[i] = 0;
            G.end_clips[i] = 0;
            for (int k=-3; k<=3; k++)
            {
                G.start_clips[i] += D.rdstats[i].rclips[start_idx[j] - k] + D.rdstats[i].lclips[start_idx[j]-k];
                G.end_clips[i] += D.rdstats[i].lclips[end_idx[j]-k] + D.rdstats[i].rclips[end_idx[j]-k];
            }
        }
        
        G.rp_cnts[i] = 0;
        for(auto &rp : D.rdstats[i].readpairs)
        {
            // Because read pair counting has enough buffer size, we do not iterate through breakpoints to avoid double-counting
            if (rp.directions.first && rp.directions.second) // FF
            {
        //        printf("FF RP: %d, %d\n", rp.positions.first, rp.positions.second);
                int dx = (int) (start_pos - rp.positions.first);
                int dy = (int) (end_pos - rp.positions.second);
                if (dx > 50 && dx < 400 && dy > 50 && dy < 400)
                {
                    G.rp_cnts[i] ++;
                }
            }
            else if (!rp.directions.first && !rp.directions.second) // RR
            {
        //        printf("RR RP: %d, %d\n", rp.positions.first, rp.positions.second);

                int dx = (int) (rp.positions.first - start_pos);
                int dy = (int) (rp.positions.second - end_pos);
                if (dx > -50 && dx < 375 && dy > -50 && dy < 375)
                {
                    G.rp_cnts[i] ++;
                }
            }
        }
        
        G.all_cnts[i] = G.split_cnts[i] + G.rp_cnts[i] + (G.start_clips[i] + G.end_clips[i])/2.0;
        
        // TODO: cut-off values are arbitrary
        if (G.rp_cnts[i] > 0 && G.all_cnts[i] > 3)
        {
            G.gt[i] = 1;
            G.nonref_mask[i] = true;
            G.ac++;
            G.ns++;
            G.startclip_sum += G.start_clips[i];
            G.endclip_sum += G.end_clips[i];
            G.split_sum += G.split_cnts[i];
            G.rp_sum += G.rp_cnts[i];
        }
        else if (G.all_cnts[i] < 1)
        {
            G.gt[i] = 0;
            G.ns++;
        }
    }
    
    if (G.ac>0 && G.ns > G.ac)
        return true;
    else
        return false;
}



void Genotyper::call_inversion(sv &S, SvData &D, SvGeno &G, std::vector<SampleStat> &stats)
{
    G.b_biallelic = true;
    
    D.clus_idx = -1;
    if (find_consensus_split(S, D, G))
    {
        G.split_flag = true;
        if (get_inv_cnts(S, D, G) && assign_inv_genotypes(S, D, G, stats))
        {
            return;
        }
        G.split_flag = false;
    }
    if (find_consensus_clip(S, D, G))
    {
        if (get_inv_cnts(S, D, G) && assign_inv_genotypes(S, D, G, stats))
        {
            G.clip_flag = true;
            return;
        }
    }
    
    // If no breakpoints identified by split reads nor by soft clips, just use the SV breakpoints as given
    D.vec_break_clusters.clear();
    D.clus_idx = -1;
    BreakCluster br;
    br.add_to_cluster(S.pos, S.end, 1);
    D.vec_break_clusters.push_back(br);

    if (get_inv_cnts(S, D, G) && assign_inv_genotypes(S, D, G, stats))
    {
        G.readpair_flag = true;
        return;
    }
    
    /*
    G.b_biallelic = true;
    
    int start_peak_ff = -1;
    int end_peak_ff = -1;
    int start_peak_rr = -1;
    int end_peak_rr = -1;
    int pairstr = 0;
	bool b_ff = find_consensus_rp(S, D, pairstr, start_peak_ff, end_peak_ff);
	pairstr = 3;
	bool b_rr = find_consensus_rp(S, D, pairstr, start_peak_rr, end_peak_rr);

	if (b_ff && b_rr && start_peak_ff < 50 && end_peak_ff < 150 && start_peak_rr > 40 && end_peak_rr > 140  && start_peak_ff < start_peak_rr && end_peak_rr > end_peak_ff)
    {
		int b_cnt = 1;
		if (b_ff && b_rr)
		{
			G.rp_pos = S.pos + (start_peak_ff+start_peak_rr)*5 - 500;
			G.rp_end = S.end + (end_peak_ff+end_peak_rr)*5 - 1500;
			b_cnt = 2;
		}
		else if (b_ff)
		{
			G.rp_pos = S.pos + start_peak_ff*10 - 500;
			G.rp_end = S.end + end_peak_ff*10 - 1500;
		}
		else if (b_rr)
		{
			G.rp_pos = S.pos + start_peak_rr*10 - 500;
			G.rp_end = S.end + end_peak_rr*10 - 1500;
		}
        
		if (G.rp_pos < G.rp_end)
		{
			for(int i=0; i<n_sample; ++i)
			{
				if (b_ff)
				{
					pairstr = 0;
					G.start_rps[i] += D.rdstats[i].rp_seq[pairstr][start_peak_ff] + D.rdstats[i].rp_seq[pairstr][start_peak_ff-1] + D.rdstats[i].rp_seq[pairstr][start_peak_ff+1];
					G.end_rps[i] += D.rdstats[i].rp_seq[pairstr][end_peak_ff] + D.rdstats[i].rp_seq[pairstr][end_peak_ff-1] + D.rdstats[i].rp_seq[pairstr][end_peak_ff+1];
				}
				if (b_rr)
				{
					pairstr = 3;
					G.start_rps[i] += D.rdstats[i].rp_seq[pairstr][start_peak_rr] + D.rdstats[i].rp_seq[pairstr][start_peak_rr-1] + D.rdstats[i].rp_seq[pairstr][start_peak_rr+1];
					G.end_rps[i] += D.rdstats[i].rp_seq[pairstr][end_peak_rr] + D.rdstats[i].rp_seq[pairstr][end_peak_rr-1] + D.rdstats[i].rp_seq[pairstr][end_peak_rr+1];
				}
				
				if (G.start_rps[i] > 15*b_cnt && G.end_rps[i] > 15*b_cnt)
				{
					G.rp_gt[i] = 2;
				}
				else if (G.start_rps[i] > 5 && G.end_rps[i] > 5 && (G.start_rps[i] + G.end_rps[i] > b_cnt * 5))
				{
					G.rp_gt[i] = 1;
				}
				else if( G.start_rps[i] < 2 && G.end_rps[i] < 2)
				{
					G.rp_gt[i] = 0;
				}
			}
			G.read_flag = true;
		}
    }
    
    int l_start = -1;
    int l_end = -1;
    int r_start = -1;
    int r_end = -1;
    
    if (find_consensus_clip_inv(S, D, l_start, l_end, r_start, r_end))
    {
        int start_clip = -1;
        int end_clip = -1;

        if (r_start >= 0)
        {
            start_clip = r_start;
        }
        else if (l_start>=0)
        {
            start_clip = l_start;
        }
        if (l_end >= 0)
        {
            end_clip = l_end;
        }
        else if (r_end >= 0)
        {
            end_clip = r_end;
        }
        
		G.clip_pos = S.pos + start_clip - 100;
		G.clip_end = S.end + end_clip - 300;

        if (G.clip_pos < G.clip_end && start_clip > 0 && start_clip<199 && end_clip>200 && end_clip<399)
        {
            for(int i=0; i<n_sample; ++i)
            {
                int n_start = 0;
                int n_end = 0;
                
                if (r_start>0)
                {
                    n_start += D.rdstats[i].rclips[r_start] + D.rdstats[i].rclips[r_start-1] + D.rdstats[i].rclips[r_start+1];
                }
                if (l_start>0)
                {
                    n_start += D.rdstats[i].lclips[l_start] + D.rdstats[i].lclips[l_start-1] + D.rdstats[i].lclips[l_start+1];
                }
                if (r_end >0)
                {
                    n_end += D.rdstats[i].rclips[r_end] + D.rdstats[i].rclips[r_end-1] + D.rdstats[i].rclips[r_end+1];
                }
                if (l_end >0)
                {
                    n_end += D.rdstats[i].lclips[l_end] + D.rdstats[i].lclips[l_end-1] + D.rdstats[i].lclips[l_end+1];
                }
                
                G.start_clips[i] = n_start;
                G.end_clips[i] = n_end;
                if (n_start > 15 && n_end > 15 && G.start_rps[i] > 0 && G.end_rps[i] > 0)
                {
                    G.clip_gt[i] = 2;
                }
                else if (n_start >= 5 && n_end >= 5 && G.start_rps[i] > 0 && G.end_rps[i] > 0)
                {
                    G.clip_gt[i] = 1;
                }
                else if (n_start < 2 && n_end < 2 && G.start_rps[i] == 0 && G.end_rps[i] == 0)
                {
                    G.clip_gt[i] = 0;
                }
            }
            G.clip_flag = true;
        }
    }
    
    G.ns = 0;
    G.ac = 0;

    if (G.read_flag || (G.read_flag && G.clip_flag))
    {
        for(int i=0; i<n_sample; ++i)
        {
            if (G.rp_gt[i] < 0)
            {
                G.gt[i] = G.clip_gt[i];
            }
            else if (G.clip_gt[i] < 0)
            {
                G.gt[i] = G.rp_gt[i];
            }
            else
            {
                if (G.clip_gt[i] == 0 && G.rp_gt[i] == 0)
                    G.gt[i] = 0;
                else
                    G.gt[i] = (G.clip_gt[i] > G.rp_gt[i]) ? G.clip_gt[i] : G.rp_gt[i];
            }
            
            if (G.gt[i] >= 0)
            {
                G.ns++;
                G.ac += G.gt[i];
            }
        }

        double callrate = (double)G.ns / n_sample;

        if (callrate>0.3 && G.ac > 0 && G.ac < (G.ns*2))
            G.b_pass = true;
    }
     */
}

int Genotyper::find_peak(std::vector<int> &seq, int start, int end)
{
    std::vector<double> peak_vec (seq.size(), 0);
    double max_peak;
    int peak_idx  = -1;
    
    int sumsq = seq[start]*seq[start] + seq[end-1]*seq[end-1];
    
    for(int i=start+1; i<end-1; ++i)
    {
        sumsq += seq[i] * seq[i];
        peak_vec[i] = (seq[i] + (double)(seq[i-1]) * 0.5 + (double)(seq[i+1]) * 0.5) / 2.0;
    }
    double stdev =sqrt( (double)(sumsq)/ (double)(end-start));
//	fprintf(stderr, "peak stdev: %f\n", stdev);
    
    max_peak = 0;
    for(int i=start+1; i<end-1; ++i)
    {
        if (peak_vec[i] >= peak_vec[i-1] && peak_vec[i] >= peak_vec[i+1] && peak_vec[i] >= stdev)
        {
		//	fprintf(stderr, "peak at %d, %d, %f\n", i, seq[i], peak_vec[i]);
			if (peak_vec[i] > max_peak)
			{
				max_peak = peak_vec[i];
				peak_idx = i;
			}
        }
    }
    return peak_idx;
}

bool Genotyper::is_pairsplit_oriented(sv &S, PairSplit &split)
{
    if (S.svtype == DEL && split.directions.first  && !split.directions.second)
        return true;
    else if (S.svtype == DUP && !split.directions.first && split.directions.second)
        return true;
    else if (S.svtype == INV && ((!split.directions.first && !split.directions.second) || (split.directions.first && split.directions.second)))
        return true;
    
    return false;
}

bool Genotyper::find_consensus_split(sv &S, SvData &D, SvGeno &G)
{
    D.vec_break_clusters.clear();

    for(int i=0; i<n_sample; ++i)
    {
        if (!G.sample_mask[i])
            continue;

        std::vector<BreakCluster> vec_br;
       // std::cerr << "Sample " << i << " has " << D.rdstats[i].splits.size() << " splits\n";

        for(auto &split : D.rdstats[i].splits)
        {
         //   std::cerr << "Sample " << i  << ": " << split.positions.first << "\t" << split.positions.second << std::endl;
            double min_dist = 1000;
            int min_idx = -1;
            
            for(int j=0; j<(int)vec_br.size(); ++j)
            {
                double dist = vec_br[j].get_distance(split.positions);
                if (dist < min_dist)
                {
                    min_dist = dist;
                    min_idx = j;
                }
            }
            if (min_idx >= 0 && min_dist <= 6) // TODO: arbitrary cutoff. Check smaller SVs how this affects.
            {
                vec_br[min_idx].add_to_cluster(split.positions);
            }
            else
            {
                BreakCluster br;
                br.add_to_cluster(split.positions);
                vec_br.push_back(br);
           //     std::cerr << "Adding to a new cluster: " << split.positions.first << "\t" << split.positions.second << std::endl;
            }
        }
        
        int max_cnt = 0;
        int max_idx = 0;
        
        for(int j=0; j<(int)vec_br.size(); ++j)
        {
            int cnt = vec_br[j].N;
            if (cnt  > max_cnt)
            {
                max_cnt = cnt;
                max_idx = j;
            }
        }
        
        if (max_cnt >= 3) // TODO: arbitrary cutoff
        {
            // add to the overall cluster
            double min_dist = 1000;
            int min_idx = -1;
            
            for(int j=0; j<(int)D.vec_break_clusters.size(); ++j)
            {
                double dist = D.vec_break_clusters[j].get_distance(vec_br[max_idx]);
                
                if (dist < min_dist)
                {
                    min_dist = dist;
                    min_idx = j;
                }
            }
            if (min_idx >= 0 && min_dist <= 6) // TODO: arbitrary cutoff. Check smaller SVs how this affects.
            {
                D.vec_break_clusters[min_idx].merge(vec_br[max_idx]);
            }
            else
            {
                D.vec_break_clusters.push_back(vec_br[max_idx]);
                //std::cerr << "Adding to a new cluster: " << split.positions.first << "\t" << split.positions.second << std::endl;
            }
        }
        //std::cout << "SV " << S.pos << "-" << S.end << ", sample " << i<< " startsplit " << vec_br[maxidx].start_mean << ", endsplit " << vec_br[maxidx].end_mean << " cnt " << maxcnt << std::endl;
        vec_br.clear();
    }
    if (D.vec_break_clusters.size() > 0)
    {
        //std::cout << "SV " << S.pos << "-" << S.end << ", " << D.vec_break_clusters.size() << " split clusters found\n";

        int max_idx = -1;
        int max_cnt = 0;
        for(int j=0; j<(int)D.vec_break_clusters.size(); ++j)
        {
           // std::cout << "Cluster " << j << " " << (int)D.vec_break_clusters[j].start_mean << "(" << D.vec_break_clusters[j].start_var << ")" << "-" << (int)D.vec_break_clusters[j].end_mean <<  "(" << D.vec_break_clusters[j].end_var << ")" << ", " << D.vec_break_clusters[j].N << " elements\n";
            if (D.vec_break_clusters[j].N > max_cnt)
            {
                max_cnt = D.vec_break_clusters[j].N;
                max_idx = j;
            }
        }
        D.clus_idx = max_idx;

        return true;
    }
    return false;
}

bool Genotyper::find_consensus_rp(sv &S, SvData &D, SvGeno &G, int pairstr, int &start_peak, int &end_peak)
{
    std::vector<int> &seq = D.all_rps[pairstr];
    int N = (int) seq.size();
	int N_first = (int)N/2;
	int N_second = (int)N/2;

	if (S.len < 1000)
	{
		N_first = 50 + (int)(S.len/20);
		N_second = 150 - (int)(S.len/20);
	}

	if (pairstr == 0)
	{
        start_peak = find_peak(D.all_rps[pairstr], 0, 50); // FF
        end_peak = find_peak(D.all_rps[pairstr], 100, 150); // FF
	}
    else if (pairstr == 3)
    {
        start_peak = find_peak(D.all_rps[pairstr], 40, 100); // RR
        end_peak = find_peak(D.all_rps[pairstr], 140, 200); // RR
    }
    else
    {
        start_peak = find_peak(seq, 0, N_first);
        end_peak = find_peak(seq, N_second, N);
    }
    
    if (start_peak>=0 && seq[start_peak]>=5 && end_peak>=0 && seq[end_peak]>=5)
    {
        return true;
    }
    else
    {
        int start_sum = 0;
        int start_sumsq = 0;
        
        int end_sum = 0;
        int end_sumsq = 0;
        
        int cnt = 0;
        
        // From each sample
        for(int i=0; i<n_sample; ++i)
        {
            if (!G.sample_mask[i])
                continue;

            if (D.rdstats[i].n_rp[pairstr] > 4)
            {
                // Supporting read pairs in pairstr in first 100
                int peak1 = -1;
				int peak2 = -1;
                
                if (pairstr == 0)
				{
					peak1 = find_peak(D.rdstats[i].rp_seq[pairstr], 0, 50);
				}
                else if (pairstr == 3)
				{
                    peak1 = find_peak(D.rdstats[i].rp_seq[pairstr], 40, 100);
				}
                else
				{
                    peak1 = find_peak(D.rdstats[i].rp_seq[pairstr], 0, N_first);
				}

                // Supporting read pairs in pairstr, in last 100
				if (pairstr == 0)
				{
					peak2 = find_peak(D.rdstats[i].rp_seq[pairstr], 100,150);
				}
				if (pairstr == 3)
				{
					peak2 = find_peak(D.rdstats[i].rp_seq[pairstr], 140, 200);
				}
				else
				{
                	peak2 = find_peak(D.rdstats[i].rp_seq[pairstr], N_second, 200);
				}
                
                if (peak1>=0 && D.rdstats[i].rp_seq[pairstr][peak1]>2 && peak2>=0 && D.rdstats[i].rp_seq[pairstr][peak2]>2)
                {
                    start_sum += peak1;
                    start_sumsq += peak1*peak1;
                    end_sum += peak2;
                    end_sumsq += peak2*peak2;
                    cnt ++;
                }
            }
        }
        
        if (cnt>0)
        {
            double mean_start = start_sum/(double)cnt;
            double std_start = sqrt(start_sumsq/(double)cnt - mean_start*mean_start);
            
            double mean_end = end_sum/(double)cnt;
            double std_end = sqrt(end_sumsq/(double)cnt - mean_end*mean_end);
            
            if (std_start<10 && std_end<10) // maybe too loose...
            {
                start_peak= round(mean_start);
                end_peak = round(mean_end);
                return true;
            }
            else
            {
                start_peak = -1;
                end_peak = -1;
            }
        }
    }
    return false;
}

bool Genotyper::find_consensus_clip(sv &S, SvData &D, SvGeno &G)
{
    D.vec_break_clusters.clear();

    if (S.svtype == DEL)
    {
        for(int i=0;i<n_sample; ++i)
        {
            if (!G.sample_mask[i])
                continue;
                
            int max_lclip = 0;
            int max_l_idx = 0;
            int max_rclip = 0;
            int max_r_idx = 0;
            for(int j=0;j<200; ++j)
            {
                if (D.rdstats[i].rclips[j] > max_rclip)
                {
                    max_r_idx = j;
                    max_rclip = D.rdstats[i].rclips[j];
                }
            }
            for(int j=200; j<400; ++j)
            {
                if (D.rdstats[i].lclips[j] > max_lclip)
                {
                    max_l_idx = j;
                    max_lclip = D.rdstats[i].lclips[j];
                }
            }
            int start_pos = S.pos + max_r_idx - 100;
            int end_pos = S.end + max_l_idx - 300;
            
            int cnt = (max_lclip > max_rclip) ? max_rclip : max_lclip; // minimum of both values

            if (cnt >= 3 && start_pos < end_pos)
            {
                double min_dist = 1000;
                int min_idx = -1;
                
                for(int j=0; j<(int)D.vec_break_clusters.size(); ++j)
                {
                    double dist = D.vec_break_clusters[j].get_distance(start_pos, end_pos);
                    
                    if (dist < min_dist)
                    {
                        min_dist = dist;
                        min_idx = j;
                    }
                }
                if (min_idx >= 0 && min_dist <= 6) // TODO: arbitrary cutoff
                {
                    D.vec_break_clusters[min_idx].add_to_cluster(start_pos, end_pos, cnt);
                }
                else
                {
                    BreakCluster br;
                    br.add_to_cluster(start_pos, end_pos, cnt);
                    D.vec_break_clusters.push_back(br);
                }
                
            }
        }
    }
    else if (S.svtype == DUP)
    {
        for(int i=0;i<n_sample; ++i)
        {
            if (!G.sample_mask[i])
                continue;
            int max_lclip = 0;
            int max_l_idx = 0;
            int max_rclip = 0;
            int max_r_idx = 0;
            for(int j=200;j<400; ++j)
            {
                if (D.rdstats[i].rclips[j] > max_rclip)
                {
                    max_r_idx = j;
                    max_rclip = D.rdstats[i].rclips[j];
                }
            }
            for(int j=0; j<200; ++j)
            {
                if (D.rdstats[i].lclips[j] > max_lclip)
                {
                    max_l_idx = j;
                    max_lclip = D.rdstats[i].lclips[j];
                }
            }
            int start_pos = S.pos + max_l_idx - 100;
            int end_pos = S.end + max_r_idx - 300;
            
            int cnt = (max_lclip > max_rclip) ? max_rclip : max_lclip; // minimum of both values

            if (cnt >= 3 && start_pos < end_pos)
            {
                double min_dist = 1000;
                int min_idx = -1;
                
                for(int j=0; j<(int)D.vec_break_clusters.size(); ++j)
                {
                    double dist = D.vec_break_clusters[j].get_distance(start_pos, end_pos);
                    
                    if (dist < min_dist)
                    {
                        min_dist = dist;
                        min_idx = j;
                    }
                }
                if (min_idx >= 0 && min_dist <=6) // TODO: arbitrary cutoff
                {
                    D.vec_break_clusters[min_idx].add_to_cluster(start_pos, end_pos, cnt);
                }
                else
                {
                    BreakCluster br;
                    br.add_to_cluster(start_pos, end_pos, cnt);
                    D.vec_break_clusters.push_back(br);
                }
                
            }
        }
    }
    else if (S.svtype == INV)
    {
        for(int i=0;i<n_sample; ++i)
        {
            if (!G.sample_mask[i])
                continue;
            int max_lclip = 0;
            int max_l_idx = 0;
            int max_rclip = 0;
            int max_r_idx = 0;
            
            for(int j=0;j<200; ++j)
            {
                if (D.rdstats[i].rclips[j] > max_rclip)
                {
                    max_r_idx = j;
                    max_rclip = D.rdstats[i].rclips[j];
                }
            }
            for(int j=0; j<200; ++j)
            {
                if (D.rdstats[i].lclips[j] > max_lclip)
                {
                    max_l_idx = j;
                    max_lclip = D.rdstats[i].lclips[j];
                }
            }
            int start_pos = S.pos;
            int start_cnt = 0;
            
            if (max_lclip >= 1 && max_rclip >= 1 && abs(max_l_idx - max_r_idx) < 10 && max_lclip + max_rclip >= 3)
            {
                start_pos = S.pos + (int)((max_l_idx + max_r_idx)/2.0) - 100;
                start_cnt = max_lclip + max_rclip;
            }
            else if (max_lclip >= 3 && max_rclip <= 1)
            {
                start_pos = S.pos + max_l_idx - 100;
                start_cnt = max_lclip;
            }
            else if (max_rclip >= 3 && max_lclip <= 1)
            {
                start_pos = S.pos + max_r_idx - 100;
                start_cnt = max_rclip;
            }
            
            max_lclip = 0;
            max_l_idx = 0;
            max_rclip = 0;
            max_r_idx = 0;
            
            for(int j=200;j<400; ++j)
            {
                if (D.rdstats[i].rclips[j] > max_rclip)
                {
                    max_r_idx = j;
                    max_rclip = D.rdstats[i].rclips[j];
                }
            }

            for(int j=200; j<400; ++j)
            {
                if (D.rdstats[i].lclips[j] > max_lclip)
                {
                    max_l_idx = j;
                    max_lclip = D.rdstats[i].lclips[j];
                }
            }

            int end_pos = S.end;
            int end_cnt = 0;
            
            if (max_lclip >= 1 && max_rclip >= 1 && abs(max_l_idx - max_r_idx) < 10 && max_lclip + max_rclip >= 3)
            {
                end_pos = S.end + (int)((max_l_idx + max_r_idx)/2.0) - 300;
                end_cnt = max_lclip + max_rclip;
            }
            else if (max_lclip >= 3 && max_rclip <= 1)
            {
                end_pos = S.end + max_l_idx - 300;
                end_cnt = max_lclip;
            }
            else if (max_rclip >= 3 && max_lclip <= 1)
            {
                end_pos = S.end + max_r_idx - 300;
                end_cnt = max_rclip;
            }
            
            int cnt = start_cnt > end_cnt ? end_cnt : start_cnt;
                        
            if (cnt >=3 && start_pos < end_pos)
            {
                double min_dist = 1000;
                int min_idx = -1;
                
                for(int j=0; j<(int)D.vec_break_clusters.size(); ++j)
                {
                    double dist = D.vec_break_clusters[j].get_distance(start_pos, end_pos);
                    
                    if (dist < min_dist)
                    {
                        min_dist = dist;
                        min_idx = j;
                    }
                }
                if (min_idx >= 0 && min_dist <=6) // TODO: arbitrary cutoff
                {
                    D.vec_break_clusters[min_idx].add_to_cluster(start_pos, end_pos, cnt);
                }
                else
                {
                    BreakCluster br;
                    br.add_to_cluster(start_pos, end_pos, cnt);
                    D.vec_break_clusters.push_back(br);
                }
            }
        }
    }
    if (D.vec_break_clusters.size() > 0)
    {
//        std::cout << "SV " << S.pos << "-" << S.end << ", " << D.vec_break_clusters.size() << " softclip clusters found\n";

        int max_idx = -1;
        int max_cnt = 0;
        for(int j=0; j<(int)D.vec_break_clusters.size(); ++j)
        {
  //          std::cout << "Cluster " << j << " " << (int)D.vec_break_clusters[j].start_mean << "(" << D.vec_break_clusters[j].start_var << ")" << "-" << (int)D.vec_break_clusters[j].end_mean <<  "(" << D.vec_break_clusters[j].end_var << ")" << ", " << D.vec_break_clusters[j].N << " elements\n";
            if (D.vec_break_clusters[j].N > max_cnt)
            {
                max_cnt = D.vec_break_clusters[j].N;
                max_idx = j;
            }
        }
        D.clus_idx = max_idx;

        return true;
    }
    return false;
}

bool Genotyper::find_consensus_clip(sv &S, SvData &D, SvGeno &G, int pairstr, int &l_peak, int &r_peak)
{
	int N_buf = 100;
	
	if (S.len < 200) 
	{
		N_buf = (int)S.len/2;
	}

    switch(pairstr)
    {
        case 1:
            r_peak = find_peak(D.all_rclips, 100 - N_buf, 100 + N_buf);
            l_peak = find_peak(D.all_lclips, 300 - N_buf, 300 + N_buf);
            break;
        case 2:
            l_peak = find_peak(D.all_lclips, 100 - N_buf, 100 + N_buf);
            r_peak = find_peak(D.all_rclips, 300 - N_buf, 300 + N_buf);
            break;
    }
    
    if (l_peak>=0 && D.all_lclips[l_peak]>=5 && r_peak>=0 && D.all_rclips[r_peak]>=5)
    {
        return true;
    }
    else
    {
        int l_sum = 0;
        int l_sumsq = 0;
        
        int r_sum = 0;
        int r_sumsq = 0;
        
        int cnt = 0;
        
        // From each sample
        for(int i=0; i<n_sample; ++i)
        {
            if (!G.sample_mask[i])
                continue;
            int n_start = 0;
            int n_end = 0;
            
            if (pairstr == 1)
                n_start = D.rdstats[i].n_rclip_start;
            else
                n_start = D.rdstats[i].n_lclip_start;
            
            if (pairstr == 1)
                n_end = D.rdstats[i].n_lclip_end;
            else
                n_end = D.rdstats[i].n_rclip_end;
            
            if (n_start >= 2 && n_end >= 2)
            {
                l_peak = -1;
                r_peak = -1;
                switch(pairstr)
                {
                    case 1:
                        r_peak = find_peak(D.rdstats[i].rclips, 100 - N_buf, 100 + N_buf);
                        l_peak = find_peak(D.rdstats[i].lclips, 300 - N_buf, 300 + N_buf);
                        break;
                    case 2:
                        l_peak = find_peak(D.rdstats[i].lclips, 100 - N_buf, 100 + N_buf);
                        r_peak = find_peak(D.rdstats[i].rclips, 300 - N_buf, 300 + N_buf);
                        break;
                    case 0:
                    case 3:
                        // inversion
                        break;
                }
                
                if (l_peak>=0 && D.rdstats[i].lclips[l_peak]>=2 && r_peak>=0 && D.rdstats[i].rclips[r_peak]>=2)
                {
                    l_sum += l_peak;
                    l_sumsq += l_peak * l_peak;
                    r_sum += r_peak;
                    r_sumsq += r_peak * r_peak;
                    cnt ++;
                }
            }
        }
        
        if (cnt>0)
        {
            double mean_l = l_sum/(double)cnt;
            double std_l = sqrt(l_sumsq/(double)cnt - mean_l * mean_l);
            
            double mean_r = r_sum/(double)cnt;
            double std_r = sqrt(r_sumsq/(double)cnt - mean_r * mean_r);
            
            if (std_l<10 && std_r<10)
            {
                l_peak = round(mean_l);
                r_peak = round(mean_r);
                return true;
            }
        }
    }
    l_peak = -1;
    r_peak = -1;
    return false;
}


bool Genotyper::find_consensus_clip_inv(sv &S, SvData &D, SvGeno &G, int &l_start, int &l_end, int &r_start, int &r_end)
{

	int N_buf = 100;

	if (S.len < 200)
	{
		N_buf = (int) S.len / 2;
	}

    l_start = find_peak(D.all_lclips, 100 - N_buf, 100 + N_buf);
    l_end = find_peak(D.all_lclips, 300 - N_buf, 300 + N_buf);
    r_start = find_peak(D.all_rclips, 100 - N_buf, 100 + N_buf);
    r_end = find_peak(D.all_rclips, 300 - N_buf, 300 + N_buf);
    
    if ( ( (l_start>=0 && D.all_lclips[l_start]>=5) || (r_start >=0 && D.all_rclips[r_start] >= 5) ||
           (r_start >= 0 && l_start >= 0 && (D.all_lclips[l_start] + D.all_rclips[r_start] >=5 ) ) )  &&
         ( (l_end  >=0 && D.all_lclips[ l_end ]>=5) || (r_end >= 0  && D.all_rclips[ r_end ] >= 5) ||
           (l_end >= 0 && r_end >= 0  && (D.all_lclips[l_end] + D.all_rclips[r_end] >= 5 ) ) ) )
    {
        return true;
    }
    else
    {
        int lstart_sum = 0;
        int lstart_sumsq = 0;
        
        int rstart_sum = 0;
        int rstart_sumsq = 0;
        
        int lend_sum = 0;
        int lend_sumsq = 0;
        
        int rend_sum = 0;
        int rend_sumsq = 0;
        
        int lstart_cnt = 0;
        int lend_cnt = 0;
        
        int rstart_cnt = 0;
        int rend_cnt = 0;
        
        // From each sample
        for(int i=0; i<n_sample; ++i)
        {
            if (!G.sample_mask[i])
                continue;
            if (D.rdstats[i].n_lclip_start + D.rdstats[i].n_rclip_start >= 3 && D.rdstats[i].n_rclip_end + D.rdstats[i].n_lclip_end >= 2)
            {
                l_start = find_peak(D.rdstats[i].lclips, 100 - N_buf, 100 + N_buf); // lstart
                l_end = find_peak(D.rdstats[i].lclips, 300 - N_buf, 300 + N_buf); // lend
                r_start = find_peak(D.rdstats[i].rclips, 100 - N_buf, 100 + N_buf); //rstart
                r_end  = find_peak(D.rdstats[i].rclips, 300 - N_buf, 300 + N_buf); //rend

                if (l_start >=0 && D.rdstats[i].lclips[l_start] >=3)
                {
                    lstart_sum += l_start;
                    lstart_sumsq += l_start * l_start;
                    lstart_cnt ++;
                }
                if (l_end >=0 && D.rdstats[i].lclips[l_end] >=3)
                {
                    lend_sum += l_end;
                    lend_sumsq += l_end * l_end;
                    lend_cnt ++;
                }
                if (r_start >=0 && D.rdstats[i].rclips[r_start] >=3)
                {
                    rstart_sum += r_start;
                    rstart_sumsq += r_start * r_start;
                    rstart_cnt ++;
                }
                if (r_end >=0 && D.rdstats[i].rclips[r_end] >=3)
                {
                    rend_sum += r_end;
                    rend_sumsq += r_end * r_end;
                    rend_cnt ++;
                }
            }
        }
        
        if (lstart_cnt > 0)
        {
            double mean_l = lstart_sum/(double)lstart_cnt;
            double std_l = sqrt(lstart_sumsq/(double)lstart_cnt - mean_l * mean_l);
            
            if (std_l < 10)
                l_start = round(mean_l);
            else
                l_start = -1;
        }
        
        if (lend_cnt > 0)
        {
            double mean_l = lend_sum/(double)lend_cnt;
            double std_l = sqrt(lend_sumsq/(double)lend_cnt - mean_l * mean_l);
            
            if (std_l < 10)
                l_end = round(mean_l);
            else
                l_end = -1;
        }
        
        if (rstart_cnt > 0)
        {
            double mean_r = rstart_sum/(double)rstart_cnt;
            double std_r = sqrt(rstart_sumsq/(double)rstart_cnt - mean_r * mean_r);
            
            if (std_r < 10)
                r_start = round(mean_r);
            else
                r_start = -1;
        }
        if (rend_cnt > 0)
        {
            double mean_r = rend_sum/(double)rend_cnt;
            double std_r = sqrt(rend_sumsq/(double)rend_cnt - mean_r * mean_r);
            
            if (std_r < 10)
                r_end = round(mean_r);
            else
                r_end = -1;
        }
        if ((l_start >=0 || r_start >=0) && (l_end >=0 || r_end>=0))
            return true;
    }
    l_start = -1;
    l_end = -1;
    r_start = -1;
    r_end = -1;
    return false;
}

bool Genotyper::assign_del_genotypes(sv &S, SvData &D, SvGeno &G, std::vector<SampleStat> &stats)
{
    // GaussianMixture dp_mix, all_mix;

    std::vector<double> norm_cnts (n_sample, 0);
    std::vector<int> lbl (n_sample, -1);

    double tmp_sum = 0;
    double tmp_cnt = 0;
    double cnt_max = -1;
    double cnt_mean = 0;

    DDMSG("assigning DEL genotypes");

    int max_ncomp = 1;

    for(int i=0; i<n_sample; ++i)
    {
        if (G.sample_mask[i])
        {
            norm_cnts[i] = G.all_cnts[i] / stats[i].avg_dp;
            if (D.dps[2][i] <= 0.65 && norm_cnts[i] >= 0.1)
            {
                tmp_sum += norm_cnts[i];
                tmp_cnt += 1;
                lbl[i] = 1;

                if (D.dps[2][i] < 0.2)
                {
                    if (norm_cnts[i] > cnt_max)
                    {
                        cnt_max = norm_cnts[i];
                    }
                    lbl[i] = 2;
                }
            }
            else if (norm_cnts[i] < 0.05 && D.dps[2][i]>0.85)
            {
                lbl[i] = 0;
            }
        }
    }
    if (tmp_cnt > 0)
    {
        cnt_mean = tmp_sum / tmp_cnt;
        max_ncomp = 2;
    }
    else
    {
        return false;
    }

    if (cnt_max > 0 )
    {
        max_ncomp = 3;
    }

    DDPRINT("cnt_mean %f, cnt_max %f \n", cnt_mean, cnt_max);

    select_model_dpcnt_mask(genostat.dpcnt_mix, D.dps[2], norm_cnts, lbl, max_ncomp, G.sample_mask, G.MAX_P_OVERLAP);
    G.dp_flag = false;
    
    if (genostat.dpcnt_mix.n_comp > 1)
    {     
        G.ns = 0;
        G.ac = 0;
        
        DDPRINT("%d components\n", genostat.dpcnt_mix.n_comp);
        for(int i=0; i<n_sample; ++i)
        {
            G.cn[i] = -1;
            G.gt[i] = -1;
            if (!G.sample_mask[i]) continue;

            int cn = genostat.dpcnt_mix.assign_dpcnt_copynumber(D.dps[2][i], norm_cnts[i]);

            if (cn == 2)
            {
                if (D.dps[2][i] < 0.8 || norm_cnts[i] > 0.1)
                {
                    G.cn[i] = -1;
                    G.gt[i] = -1;
                }
                else
                {
                    G.cn[i] = 2;
                    G.gt[i] = 0; // 0/0
                    G.ns ++;
                }
            }
            else
            {
                if (cn == -1 && D.dps[2][i] > 0.85 && norm_cnts[i] < 0.05)
                {
                    G.cn[i] = 2;
                    G.gt[i] = 0;
                    G.ns ++;
                }
                else if (D.dps[2][i] > 0.7 || norm_cnts[i] < 0.05)
                {
                    G.cn[i] = -1;
                    G.gt[i] = -1;
                }
                else if (cn == 1)
                {
                    G.cn[i] = 1;
                    G.gt[i] = 1; // 0/1
                    G.ac ++;
                    G.ns ++;
                }
                else if (cn == 0)
                {
                    G.cn[i] = 0;
                    G.gt[i] = 2; // 1/1
                    G.ns ++;
                    G.ac += 2;
                }
            }
        
        }

        double callrate = (double)G.ns / G.n_effect;
        DDPRINT("NS %d, AC %d, Call rate %f\n", G.ns, G.ac, callrate);
        
        if (callrate>0.5 && G.ac > 0 && G.ac < G.ns*2) // if successful, return genotypes
        {
            G.dpcnt_flag = true;
            G.cnt_flag = false;
//            G.gmix.n_comp = G.gmix.n_comp + 1;
//           G.gmix.Comps.resize(G.gmix.n_comp);
//            G.gmix.Comps[G.gmix.n_comp-1].Mean = dpcnt_mix.Comps[0].Mean[0];
//            G.gmix.Comps[G.gmix.n_comp-1].Stdev = dpcnt_mix.Comps[0].Cov[0];
//           G.gmix.Comps[G.gmix.n_comp-1].Alpha = dpcnt_mix.Comps[0].Alpha;
            G.gmix2 = genostat.dpcnt_mix;
            G.b_pass = true;
            DDMSG("DPCNT clustering successful");
            return true;
        }
        else
        {
            G.dp_flag = false;
            G.cnt_flag = false;
            return false;
        }
    }
    return false;    
}
// bool Genotyper::assign_del_genotypes(sv &S, SvData &D, SvGeno &G, std::vector<SampleStat> &stats)
// {
//     // GaussianMixture dp_mix, all_mix;
//     GaussianMixture2 dpcnt_mix;
//     std::vector<double> norm_cnts (n_sample, 0);
//     for(int i=0; i<n_sample; ++i)
//     {
//         if (G.sample_mask[i])
//             norm_cnts[i] = (G.split_cnts[i]  + G.rp_cnts[i]) / stats[i].avg_dp;
//     }
        
//     //if (G.split_sum / G.ac < 1 && G.rp_sum / G.ac < 1 && ((G.startclip_sum/(G.endclip_sum+0.001) > 5) || (G.endclip_sum/(G.startclip_sum+0.001) > 5)))
//     //    return false;
        
//     // Check whether the 'variant' depth is separated from 'non-variant' depth
//     // dp_mix.estimate(D.dps[2], G.gt, 2);

// //    dpcnt_mix.estimate(D.dps[2], norm_cnts, G.gt, 2);

// //    fprintf(stderr, "dpmix cluster\n");
//  //   dp_mix.print(stderr);
    
// //    double err_bound = dp_mix.Comps[0].Stdev;
//     double err_bound = dpcnt_mix.Comps[0].Cov[0];
//     if (err_bound < 0.2) err_bound = 0.2;
    
//     // BIC? 

//     // double P_cutoff = 0.5;
//     // if (G.split_sum / G.ac > 2)
//     // {
//     //     P_cutoff += 0.2;
//     // }
//     // if (G.rp_sum / G.ac > 2)
//     // {
//     //     P_cutoff += 0.2;
//     // }
    
// //    printf("Sum/AC: %f, %f, %f, %f\n", G.split_sum/G.ac, G.rp_sum/G.ac, G.startclip_sum/G.ac, G.endclip_sum/G.ac);

//     G.dp_flag = false;
// //    if (dp_mix.Comps[0].Mean>=0.8 && dp_mix.Comps[0].Mean - err_bound > dp_mix.Comps[1].Mean && dp_mix.p_overlap <= P_cutoff)
// //    if (dp_mix.Comps[0].Mean>=0.8 && dp_mix.Comps[0].Mean - err_bound > dp_mix.Comps[1].Mean)
//     if (1)
//     {
//         // if yes, then try to EM with masks, variants only
//         std::vector<std::vector<double> > alt_means = { {0.5}, {0.5, 0.0} };
//         select_model(G.gmix, alt_means, D.dps[2], G.nonref_mask, 0.3);
        
//         G.ns = 0;
//         G.ac = 0;
        
//         for(int i=0; i<n_sample; ++i)
//         {
//             if (!G.sample_mask[i]) continue;

//             double d0 = dp_mix.Comps[0].Mean - D.dps[2][i];
//             double d1 = D.dps[2][i] - G.gmix.Comps[0].Mean;
            
//             G.gt[i] = -1;
//             G.cn[i] = -1;
            
//             if (d0>0 &&  G.all_cnts[i]>0 && d1 < d0 * ((G.split_cnts[i] + G.all_cnts[i] - 1) / 2.0))
//             {
//                 if (D.dps[2][i] < 0.1 || (G.gmix.Comps.size()>1 && D.dps[2][i] < G.gmix.Comps[1].Stdev))
//                 {
//                     G.gt[i] = 2;
//                     G.cn[i] = 0;
//                     G.ns++;
//                     G.ac += 2;
//                 }
//                 else // if there are few ambiguous samples, do not code alternative allele solely based on depth
//                 {
//                     G.gt[i] = 1;
//                     G.cn[i] = 1;
//                     G.ns++;
//                     G.ac++;
//                 }
//             }
//             else
//             {
//                 if (D.dps[2][i] < dp_mix.Comps[0].Mean - (err_bound / (G.all_cnts[i]+1.0) ) )  // missing
//                     G.gt[i] = -1;
//                 else if (d0 < d1)
//                 {
//                     G.gt[i] = 0;
//                     G.cn[i] = 2;
//                     G.ns++;
//                 }
//                 else
//                 {
//                     G.gt[i] = -1;
//                 }
//             }
//         }
        
//         double callrate = (double)G.ns / G.n_effect;
        
//         if (callrate>0.5 && G.ac > 0 && G.ac < G.ns*2) // if successful, return genotypes
//         {
//             G.dp_flag = true;
//             G.cnt_flag = false;
//             G.gmix.n_comp = G.gmix.n_comp + 1;
//             G.gmix.Comps.resize(G.gmix.n_comp);
//             G.gmix.Comps[G.gmix.n_comp-1].Mean = dp_mix.Comps[0].Mean;
//             G.gmix.Comps[G.gmix.n_comp-1].Stdev = dp_mix.Comps[0].Stdev;
//             G.gmix.Comps[G.gmix.n_comp-1].Alpha = dp_mix.Comps[0].Alpha;
//             G.b_pass = true;
//             return true;
//         }
//         else
//         {
//             G.dp_flag = false;
//             G.cnt_flag = false;
//             return false;
//         }
//     }
    
//     if (dp_mix.Comps[0].Mean - 0.1 > dp_mix.Comps[1].Mean)
//     {
//         // genotype using read counts only, not using depth at all

//         all_mix.estimate(norm_cnts, G.gt, 2);
//      //   fprintf(stderr, "allmix estimates:");
//      //   all_mix.print(stderr);
        
//         std::vector<std::vector<double> > alt_means = { {all_mix.Comps[1].Mean}, {all_mix.Comps[1].Mean, all_mix.Comps[1].Mean * 2.0} };

//         select_model_mask_1d(G.gmix, alt_means, norm_cnts, G.nonref_mask, 0.3);
        
//         double err_bound = all_mix.Comps[0].Stdev;
        
//         if (err_bound < 0.1) err_bound = 0.1;

//         // TODO: Arbitrary cutoffs
//         if (all_mix.Comps[0].Mean>0.1 || all_mix.Comps[1].Mean - err_bound < all_mix.Comps[0].Mean || all_mix.p_overlap > 0.5)
//             return false;

//         G.ns = 0;
//         G.ac = 0;
        
//         for(int i=0; i<n_sample; ++i)
//         {
//             if (!G.sample_mask[i])
//                 continue;
//             double d0 = norm_cnts[i] - all_mix.Comps[0].Mean ;
//             double d1 = G.gmix.Comps[0].Mean - norm_cnts[i];
            
//             G.gt[i] = -1;
//             G.cn[i] = -1;
//             // TODO: test genotypes
//             if (d0>0 &&  G.all_cnts[i]>0 && d1 < d0 * ((G.split_cnts[i] + G.rp_cnts[i] + G.all_cnts[i] - 1) / 2.0))
//             {
//                 if (G.gmix.Comps.size()>1  && ((G.gmix.Comps[1].Mean - norm_cnts[i] < (norm_cnts[i] - G.gmix.Comps[0].Mean) * 0.5) || D.dps[2][i] < 0.1))
//                 {
//                     G.gt[i] = 2;
//                     G.cn[i] = 0;
//                     G.ns++;
//                     G.ac += 2;
//                 }
//                 else // if there are few ambiguous samples, do not code alternative allele solely based on depth
//                 {
//                     G.gt[i] = 1;
//                     G.cn[i] = 1;
//                     G.ns++;
//                     G.ac++;
//                 }
//             }
//             else
//             {
//                 if (norm_cnts[i] > all_mix.Comps[0].Mean + (err_bound / (G.all_cnts[i]+1.0) ) )  // missing
//                 {
//                     G.gt[i] = -1;
//                     G.cn[i] = -1;
//                 }
//                 else if (d0 < d1)
//                 {
//                     G.cn[i] = 2;
//                     G.gt[i] = 0;
//                     G.ns++;
//                 }
//                 else
//                 {
//                     G.cn[i] = -1;
//                     G.gt[i] = -1;
//                 }
//             }
//         }
//         double callrate = (double)G.ns / G.n_effect;

//         if (callrate>0.5 && G.ac > 0 && G.ac < G.ns*2) // if successful, return genotypes
//         {
//             G.dp_flag = false;
//             G.cnt_flag = true;
//             G.gmix.n_comp = G.gmix.n_comp + 1;
//             G.gmix.Comps.resize(G.gmix.n_comp);
//             G.gmix.Comps[G.gmix.n_comp-1].Mean = dp_mix.Comps[0].Mean;
//             G.gmix.Comps[G.gmix.n_comp-1].Stdev = dp_mix.Comps[0].Stdev;
//             G.gmix.Comps[G.gmix.n_comp-1].Alpha = dp_mix.Comps[0].Alpha;
//             G.b_pass = true;
//             return true;
//         }
//         else
//         {
//             G.dp_flag = false;
//             G.cnt_flag = false;
//             return false;
//         }
//     }

//     return false;
// }
        
bool Genotyper::get_del_cnts(sv &S, SvData &D, SvGeno &G)
{
    G.startclip_sum = 0;
    G.endclip_sum = 0;
    G.rp_sum = 0;
    G.split_sum = 0;
    G.ns = 0;
    G.ac = 0;
    
    int start_pos  =(int)(D.vec_break_clusters[D.clus_idx].start_mean+0.5);
    int end_pos = (int)(D.vec_break_clusters[D.clus_idx].end_mean+0.5);
    
    int n_clus = (int)D. vec_break_clusters.size();
    
    std::vector<int> r_clip_idx (n_clus, 0);
    std::vector<int> l_clip_idx (n_clus, 0);
    
    for(int j=0; j<n_clus; ++j)
    {
        r_clip_idx[j] = ((int)(D.vec_break_clusters[j].start_mean+0.5) - S.pos + 100);
        l_clip_idx[j] = ((int)(D.vec_break_clusters[j].end_mean+0.5) - S.end + 300);
        
        if (r_clip_idx[j] < 3) r_clip_idx[j] = 3;
        if (r_clip_idx[j] > 196) r_clip_idx[j] = 196;
        if (l_clip_idx[j] < 203) l_clip_idx[j] = 203;
        if (l_clip_idx[j] > 396) l_clip_idx[j] = 396;
    }
    
    for(int i=0; i< n_sample; ++i)
    {
        G.nonref_mask[i] = false;
        G.gt[i] = -1;
        G.cn[i] = -1;
        if (!G.sample_mask[i])
            continue;

        G.split_cnts[i] = 0;
        // Count the number of split reads
        for(auto &split : D.rdstats[i].splits)
        {
            for (int j=0; j<n_clus; ++j)
            {
                if (D.vec_break_clusters[j].get_distance(split.positions) <= 6)
                {
                    G.split_cnts[i] ++;
                }
            }
        }
        
        for(int j=0; j<n_clus; ++j)
        {
            // Count the number of soft clips around breakpoints
            for(int k=r_clip_idx[j]-3; k<r_clip_idx[j] + 4; ++k)
            {
                G.start_clips[i] += D.rdstats[i].rclips[k];
            }
            for(int k=l_clip_idx[j]-3; k<l_clip_idx[j] + 4; ++k)
            {
                G.end_clips[i] += D.rdstats[i].lclips[k];
            }
        }
        
        G.rp_cnts[i] = 0;
        for(auto &rp : D.rdstats[i].readpairs)
        {
            // Because read pair counting has enough buffer size, we do not iterate through breakpoints to avoid double-counting
            int dx = (int) (start_pos - rp.positions.first);
            int dy = (int) (rp.positions.second - end_pos);
            if (dx > 75 && dx < 375 && dy > -75 && dy < 300)
            {
                G.rp_cnts[i] ++;
            }
        }
        
        G.all_cnts[i] = G.split_cnts[i] + G.rp_cnts[i] + (G.start_clips[i] + G.end_clips[i])/2.0;
        
        // TODO: cut-off values are arbitrary
        if (G.split_cnts[i] >= 3 || G.split_cnts[i]+G.rp_cnts[i] >= 5 || (G.split_cnts[i]+G.rp_cnts[i] > 0 && G.all_cnts[i] >= 5))
        {
            G.gt[i] = 1;
            G.nonref_mask[i] = true;
            G.ac++;
            G.ns++;
            G.startclip_sum += G.start_clips[i];
            G.endclip_sum += G.end_clips[i];
            G.split_sum += G.split_cnts[i];
            G.rp_sum += G.rp_cnts[i];
        }
        else if (G.all_cnts[i] < 1)
        {
            G.gt[i] = 0;
            G.ns++;
        }
    }
    
    if (G.ac>0 && G.ns > G.ac)
        return true;
    else
        return false;
}

void Genotyper::print_genodata(sv &S, SvData &D, SvGeno &G)
{

    for(int i=0; i<n_sample;++i)
    {
        if (G.sample_mask[i])
        {
            printf("TRUE\t");
        }
        else
        {
            printf("FALSE\t");
        }
        

        printf("%d\t%d\t%d\t%d\t%d\t%f\t%f\t%f\t%d\n", (int)G.split_cnts[i], (int)G.rp_cnts[i], (int)G.start_clips[i], (int)G.end_clips[i], (int)G.all_cnts[i], D.dps[0][i], D.dps[1][i], D.dps[2][i], G.gt[i]);
    }
}

void Genotyper::call_deletion(sv &S, SvData &D, SvGeno &G, std::vector<SampleStat> &stats)
{
    G.b_biallelic = true;
    
    bool split_flag = false;
    
    // First, try split read - best evidence
    if (find_consensus_split(S, D, G))
    {
        if (get_del_cnts(S, D, G) && assign_del_genotypes(S, D, G, stats))
        {
            // print_genodata(S, D, G);

            G.split_flag = true;
            return;
        }
        else
        {
            G.split_flag = false;
        }
        split_flag = true;
        
#ifdef DDEBUG
        for(int i=0; i<(int)D.vec_break_clusters.size(); ++i)
            fprintf(stderr, "SPLIT found: %f-%f, %d\n", D.vec_break_clusters[i].start_mean, D.vec_break_clusters[i].end_mean, D.vec_break_clusters[i].N);  
#endif
    }
    
    if (!split_flag && find_consensus_clip(S, D, G))
    {
        D.clus_idx = -1;
        G.clip_flag = true;
        if (get_del_cnts(S, D, G) && assign_del_genotypes(S, D, G, stats))
        {
            return;
        }
        /*
        for(int i=0; i<(int)D.vec_break_clusters.size(); ++i)
            fprintf(stderr, "CLIP found: %f-%f, %d\n", D.vec_break_clusters[i].start_mean, D.vec_break_clusters[i].end_mean, D.vec_break_clusters[i].N); */
    }
    
    if (!split_flag && !G.clip_flag)
    {
        // If no breakpoints identified by split reads nor by soft clips, just use the SV breakpoints as given
        D.vec_break_clusters.clear();
        D.clus_idx = -1;
        BreakCluster br;
        br.add_to_cluster(S.pos, S.end, 1);
        D.vec_break_clusters.push_back(br);

        if (get_del_cnts(S, D, G) && assign_del_genotypes(S, D, G, stats))
        {
            G.readpair_flag = true;
            return;
        }
    }

    // use pre_post filtering for depth-only
    if (!G.b_pre || !G.b_post) return;
    
    // If nothing works, use Depth-only clustering with stringent cutoff
    std::vector< std::vector<double> > means = { {1.0}, {1.0, 0.5}, {1.0, 0.5, 0.0}};
    G.MAX_P_OVERLAP = 0.2;
	double best_dp_idx = 2;
	std::vector<double> &var_depth = D.dps[best_dp_idx];

    G.ns = 0;
    G.ac = 0;
    
    // Depth-based clustering genotyping
    select_model_mask_1d(G.gmix, means, D.dps[best_dp_idx], G.sample_mask, G.MAX_P_OVERLAP);
    // select_model(G.gmix, means, D.dps[best_dp_idx], G.MAX_P_OVERLAP);
    // depth clustering
    if (G.gmix.n_comp > 1 && G.gmix.ordered() && G.gmix.Comps[0].Alpha > 0.5)
    {
        G.dp_flag = true;
        // success
        // assign dp-genotypes
        for(int i=0; i<n_sample; ++i)
        {
            if (!G.sample_mask[i]) continue;

            int cn = G.gmix.assign_copynumber(var_depth[i]);
            
            if (cn == 2)
            {
                G.cn[i] = 2;
                G.gt[i] = 0; // 0/0
            }
            else if (var_depth[i] > 0.6)
            {
                G.cn[i] = -1;
                G.gt[i] = -1;
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
    if (D.multi_dp) //dp2 has more than 2 vectors
    {
        // dp100 genotyping
        select_model_mask_2d(G.gmix2, means, D.dps[dp2_idx], D.dps[dp2_idx+1], G.sample_mask, G.MAX_P_OVERLAP);

        // 2-d genotyping
        if (G.gmix2.n_comp>1 && G.gmix2.ordered())
        {
            G.dp2_flag = true;
            //assign dp2 genotypes

            for(int i=0; i<n_sample; ++i)
            {
                if (!G.sample_mask[i]) continue;

                int cn = G.gmix2.assign_copynumber(D.dps[dp2_idx][i], D.dps[dp2_idx+1][i]);
                if (cn == 2)
                {
                    G.cn[i] = 2;
                    G.gt[i] = 0; // 0/0
                }
                else if (var_depth[i] > 0.6)
                {
                    G.cn[i] = -1;
                    G.gt[i] = -1;
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

    // If clustered, filter false positive variants
    for (int i=0;i<n_sample; ++i)
    {        
        if (!G.sample_mask[i]) continue;
        if (abs(D.dps[0][i] - G.dp_pre_mean) > 2.0*G.dp_pre_std || abs(D.dps[1][i] - G.dp_post_mean) > 2.0* G.dp_post_std || (G.gt[i] > 0 && abs(var_depth[i] - D.dps[1][i]) < 2.0*G.dp_post_std) || (G.gt[i] > 0 && abs(var_depth[i] - D.dps[0][i]) < 2.0*G.dp_pre_std))
        {
            G.cn[i] = -1;
            G.gt[i] = -1;
        }
    }
    
    int n_het = 0;
    int n_hom = 0;
    int n_ref = 0;

    for(int i=0; i<n_sample; ++i)
    {
        if (!G.sample_mask[i]) continue;

        if (G.gt[i] >=0)
        {
            G.ns += 1;
            G.ac += G.gt[i];

            if (G.gt[i] == 0)
            {
                n_ref ++;
            }
            else if (G.gt[i] == 1)
            {
                n_het ++;
            }
            else if (G.gt[i] == 2)
            {
                n_hom ++;
            }
        }
    }

    double callrate = (double)G.ns / G.n_effect;

    // Excessive heterozygosity (basically all-het case)
    if ( n_het > ((double)G.ns * 0.75) && n_hom < ((double)0.05 * G.ns) )
    {
        G.b_pass = false;
    }
    else if ((G.dp_flag  || (G.dp_flag && G.dp2_flag)) &&  callrate>0.5 && G.ac > 0 && G.ac < G.ns*2)
    {
        G.b_pass = true;

        //G.gmix.print(stdout);
        //G.gmix2.print(stdout);

    }
}

bool Genotyper::assign_dup_genotypes(sv &S, SvData &D, SvGeno &G, std::vector<SampleStat> &stats)
{
  // GaussianMixture dp_mix, all_mix;

    std::vector<double> norm_cnts (n_sample, 0);
    std::vector<int> lbl (n_sample, -1);

    double tmp_sum = 0;
    double tmp_cnt = 0;
    double cnt_max = -1;
    double cnt_mean = 0;

    DDMSG("assigning DUP genotypes");

    int max_ncomp = 1;

    for(int i=0; i<n_sample; ++i)
    {
        if (G.sample_mask[i])
        {
            norm_cnts[i] = G.all_cnts[i] / stats[i].avg_dp;
            if (D.dps[2][i] >= 1.35 && norm_cnts[i] >= 0.1)
            {
                tmp_sum += norm_cnts[i];
                tmp_cnt += 1;
                lbl[i] = 1;

                if (D.dps[2][i] > 1.85)
                {
                    if (norm_cnts[i] > cnt_max)
                    {
                        cnt_max = norm_cnts[i];
                    }
                    
                    lbl[i] = round(D.dps[2][i]*2.0) - 2;
                }
            }
            else if (norm_cnts[i] < 0.05 && D.dps[2][i]<1.2 && D.dps[2][i] > 0.8)
            {
                lbl[i] = 0;
            }
            if (lbl[i] > max_ncomp-1)
            {
                max_ncomp = lbl[i] + 1;
            }
        }
    }
    if (tmp_cnt > 0)
    {
        cnt_mean = tmp_sum / tmp_cnt;
    }
    else
    {
        return false;
    }

    DDPRINT("cnt_mean %f, cnt_max %f \n", cnt_mean, cnt_max);

    select_model_dpcnt_mask(genostat.dpcnt_mix, D.dps[2], norm_cnts, lbl, max_ncomp, G.sample_mask, G.MAX_P_OVERLAP);
    G.dp_flag = false;
    
    if (genostat.dpcnt_mix.n_comp > 1)
    {     
        G.ns = 0;
        G.ac = 0;
        
        DDPRINT("%d components\n", genostat.dpcnt_mix.n_comp);
        for(int i=0; i<n_sample; ++i)
        {
            G.cn[i] = -1;
            G.gt[i] = -1;
            if (!G.sample_mask[i]) continue;

            int cn = genostat.dpcnt_mix.assign_dpcnt_copynumber(D.dps[2][i], norm_cnts[i]);

            if (cn == 2)
            {
                if (D.dps[2][i] < 0.8 || D.dps[2][i] > 1.2 || norm_cnts[i] > 0.1)
                {
                    G.cn[i] = -1;
                    G.gt[i] = -1;
                }
                else
                {
                    G.cn[i] = 2;
                    G.gt[i] = 0; // 0/0
                    G.ns ++;
                }
            }
            else
            {
                if (cn == -1 && D.dps[2][i] < 1.2 && D.dps[2][i] > 0.8 && norm_cnts[i] < 0.05)
                {
                    G.cn[i] = 2;
                    G.gt[i] = 0;
                    G.ns ++;
                }
                else if (D.dps[2][i] < 1.3 || norm_cnts[i] < 0.05)
                {
                    G.cn[i] = -1;
                    G.gt[i] = -1;
                }
                else if (cn == 3)
                {
                    G.cn[i] = 3;
                    G.gt[i] = 1; // 0/1
                    G.ac ++;
                    G.ns ++;
                }
                else if (cn > 3)
                {
                    G.cn[i] = cn;
                    G.gt[i] = 2; // 1/1
                    G.ns ++;
                    G.ac += 2;
                }
                else if (cn >= 0 && cn < 2)
                {
                    G.cn[i] = cn;
                    G.gt[i] = -1;
                    G.ns ++;
                }
            }
        
        }

        double callrate = (double)G.ns / G.n_effect;
        DDPRINT("NS %d, AC %d, Call rate %f\n", G.ns, G.ac, callrate);
        
        if (callrate>0.5 && G.ac > 0 && G.ac < G.ns*2) // if successful, return genotypes
        {
            G.dpcnt_flag = true;
            G.cnt_flag = false;
//            G.gmix.n_comp = G.gmix.n_comp + 1;
//           G.gmix.Comps.resize(G.gmix.n_comp);
//            G.gmix.Comps[G.gmix.n_comp-1].Mean = dpcnt_mix.Comps[0].Mean[0];
//            G.gmix.Comps[G.gmix.n_comp-1].Stdev = dpcnt_mix.Comps[0].Cov[0];
//           G.gmix.Comps[G.gmix.n_comp-1].Alpha = dpcnt_mix.Comps[0].Alpha;
            G.gmix2 = genostat.dpcnt_mix;
            G.b_pass = true;
            check_biallelic(G);
            DDMSG("DPCNT clustering successful");
            return true;
        }
        else
        {
            G.dp_flag = false;
            G.cnt_flag = false;
            return false;
        }
    }
    return false;    
}

// bool Genotyper::assign_dup_genotypes(sv &S, SvData &D, SvGeno &G, std::vector<SampleStat> &stats)
// {
//     GaussianMixture dp_mix, all_mix;

//     if (G.split_sum / G.ac < 1 && G.rp_sum / G.ac < 1 && ((G.startclip_sum/(G.endclip_sum+0.001) > 5) || (G.endclip_sum/(G.startclip_sum+0.001) > 5)))
//         return false;

//     // Check whether the 'variant' depth is separated from 'non-variant' depth
//     dp_mix.estimate(D.dps[2], G.gt, 2);
    
//     double err_bound = dp_mix.Comps[0].Stdev;
    
//     // TODO: Arbitrary cutoffs
//     if (err_bound < 0.2) err_bound = 0.2;
    
//     if (dp_mix.Comps[0].Mean>=0.8 && dp_mix.Comps[0].Mean <= 1.2 && dp_mix.Comps[1].Mean - err_bound >= dp_mix.Comps[0].Mean)
//     {
//         // if yes, then try to EM with masks, variants only
//         std::vector<std::vector<double> > alt_means =  { {1.5}, {1.5, 2.0}, {1.5, 2.0, 2.5}, {1.5, 2.0, 2.5, 3.0}, {1.5, 2.0, 2.5, 3.0, 3.5}, {1.5, 2.0, 2.5, 3.0, 3.5, 4.0} };
//         select_model_mask_1d(G.gmix, alt_means, D.dps[2], G.nonref_mask, 0.3);
            
//         G.ns = 0;
//         G.ac = 0;
        
//         for(int i=0; i<n_sample; ++i)
//         {
//             if (!G.sample_mask[i]) continue;

//             double d0 = D.dps[2][i] - dp_mix.Comps[0].Mean;
//             double d1 = G.gmix.Comps[0].Mean - D.dps[2][i];
            
//             G.gt[i] = -1;
//             G.cn[i] = -1;
            
//             if (d0>0 && G.all_cnts[i]>0 && d1 < d0 * ((G.split_cnts[i] + G.all_cnts[i] - 1) / 2.0))
//             {
//                 if (G.gmix.Comps.size()>1 && D.dps[2][i] > G.gmix.Comps[1].Mean - G.gmix.Comps[1].Stdev)
//                 {
//                     G.gt[i] = 2;
//                     G.cn[i] = (int)( D.dps[2][i] * 2.0 );
//                     G.ns++;
//                     G.ac += 2;
//                 }
//                 else
//                 {
//                     G.gt[i] = 1;
//                     G.cn[i] = 1;
//                     G.ns++;
//                     G.ac++;
//                 }
//             }
//             else
//             {
//                 if (D.dps[2][i] > dp_mix.Comps[0].Mean + (err_bound / (G.all_cnts[i]+1.0) ) || G.split_cnts[i]>3 || (G.split_cnts[i]>0 && G.all_cnts[i]>5) || G.all_cnts[i] > 10)  // missing
//                 {
//                     G.gt[i] = -1;
//                 }
//                 else if (d0 < d1)
//                 {
//                     G.gt[i] = 0;
//                     G.cn[i] = 2;
//                     G.ns++;
//                 }
//                 else
//                 {
//                     G.gt[i] = -1;
//                 }
//             }
//         }

//         fflush(stdout);
//         double callrate = (double)G.ns / G.n_effect;

//         if (callrate>0.5 && G.ac > 0 && G.ac < G.ns*2) // if successful, return genotypes
//         {
//             G.b_pass = true;
//             check_biallelic(G);
//             return true;
//         }
//     }
//     if (dp_mix.Comps[1].Mean - 0.1 > dp_mix.Comps[0].Mean)
//     {
//         // genotype using read counts only, not using depth at all
//         std::vector<double> norm_cnts (n_sample, 0);
//         for(int i=0; i<n_sample; ++i)
//         {
//             if (G.sample_mask[i])
//                norm_cnts[i] = (G.split_cnts[i]  + G.rp_cnts[i]) / stats[i].avg_dp ;
//         }
        
//         all_mix.estimate(norm_cnts, G.gt, 2);
        
//         std::vector<std::vector<double> > alt_means = { {all_mix.Comps[1].Mean}, {all_mix.Comps[1].Mean, all_mix.Comps[1].Mean * 2.0} };

//         select_model_mask_1d(G.gmix, alt_means, norm_cnts, G.nonref_mask, 0.3);
        
//         double err_bound = all_mix.Comps[0].Stdev;
        
//         if (err_bound < 0.1) err_bound = 0.1;

//         // TODO: Arbitrary cutoffs
//         if (all_mix.Comps[0].Mean>0.1 || all_mix.Comps[1].Mean - err_bound < all_mix.Comps[0].Mean || all_mix.p_overlap > 0.5)
//             return false;

//         G.ns = 0;
//         G.ac = 0;
        
//         for(int i=0; i<n_sample; ++i)
//         {
//             if (!G.sample_mask[i]) continue;

//             double d0 = norm_cnts[i] - all_mix.Comps[0].Mean ;
//             double d1 = G.gmix.Comps[0].Mean - norm_cnts[i];
            
//             G.gt[i] = -1;
//             G.cn[i] = -1;
//             // TODO: test genotypes
//             if (d0>0 &&  G.all_cnts[i]>0 && d1 < d0 * ((G.split_cnts[i] + G.rp_cnts[i] + G.all_cnts[i] - 1) / 2.0))
//             {
//                 if (G.gmix.Comps.size()>1  && ((G.gmix.Comps[1].Mean - norm_cnts[i] < (norm_cnts[i] - G.gmix.Comps[0].Mean) * 0.5) || D.dps[2][i] > 2.0))
//                 {
//                     G.gt[i] = 2;
//                     G.cn[i] = (int)( D.dps[2][i] * 2.0 );
//                     G.ns++;
//                     G.ac += 2;
//                 }
//                 else // if there are few ambiguous samples, do not code alternative allele solely based on depth
//                 {
//                     G.gt[i] = 1;
//                     G.cn[i] = 3;
//                     G.ns++;
//                     G.ac++;
//                 }
//             }
//             else
//             {
//                 if (norm_cnts[i] > all_mix.Comps[0].Mean + (err_bound / (G.all_cnts[i]+1.0) ) )  // missing
//                 {
//                     G.gt[i] = -1;
//                     G.cn[i] = -1;
//                 }
//                 else if (d0 < d1)
//                 {
//                     G.cn[i] = 2;
//                     G.gt[i] = 0;
//                     G.ns++;
//                 }
//                 else
//                 {
//                     G.cn[i] = -1;
//                     G.gt[i] = -1;
//                 }
//             }
//         }
//         double callrate = (double)G.ns / G.n_effect;


//         if (callrate>0.5 && G.ac > 0 && G.ac < G.ns*2) // if successful, return genotypes
//         {
//             G.dp_flag = false;
//             G.cnt_flag = true;
//             G.gmix.n_comp = G.gmix.n_comp + 1;
//             G.gmix.Comps.resize(G.gmix.n_comp);
//             G.gmix.Comps[G.gmix.n_comp-1].Mean = dp_mix.Comps[0].Mean;
//             G.gmix.Comps[G.gmix.n_comp-1].Stdev = dp_mix.Comps[0].Stdev;
//             G.gmix.Comps[G.gmix.n_comp-1].Alpha = dp_mix.Comps[0].Alpha;
//             G.b_pass = true;
//             return true;
//         }
//         else
//         {
//             G.dp_flag = false;
//             G.cnt_flag = false;
//             return false;
//         }
//     }
//     return false;
// }
        
bool Genotyper::get_dup_cnts(sv &S, SvData &D, SvGeno &G)
{
    G.startclip_sum = 0;
    G.endclip_sum = 0;
    G.rp_sum = 0;
    G.split_sum = 0;
    G.ns = 0;
    G.ac = 0;
    
    int start_pos  =(int)(D.vec_break_clusters[D.clus_idx].start_mean+0.5);
    int end_pos = (int)(D.vec_break_clusters[D.clus_idx].end_mean+0.5);
    
    int n_clus = (int)D. vec_break_clusters.size();
    
    std::vector<int> r_clip_idx (n_clus, 0);
    std::vector<int> l_clip_idx (n_clus, 0);
    
    for(int j=0; j<n_clus; ++j)
    {
        l_clip_idx[j] = ((int)(D.vec_break_clusters[j].start_mean+0.5) - S.pos + 100);
        r_clip_idx[j] = ((int)(D.vec_break_clusters[j].end_mean+0.5) - S.end + 300);
        
        if (l_clip_idx[j] < 3) l_clip_idx[j] = 3;
        if (l_clip_idx[j] > 196) l_clip_idx[j] = 196;
        if (r_clip_idx[j] < 203) r_clip_idx[j] = 203;
        if (r_clip_idx[j] > 396) r_clip_idx[j] = 396;
    }
    
    for(int i=0; i< n_sample; ++i)
    {
        G.nonref_mask[i] = false;
        G.gt[i] = -1;
        G.cn[i] = -1;
        if (!G.sample_mask[i])
            continue;

        G.split_cnts[i] = 0;
        // Count the number of split reads
        for(auto &split : D.rdstats[i].splits)
        {

            for (int j=0; j<n_clus; ++j)
            {
                if (D.vec_break_clusters[j].get_distance(split.positions) <=6)
                {
                    G.split_cnts[i] ++;
                }
            }
        }
        
        for(int j=0; j<n_clus; ++j)
        {
            // Count the number of soft clips around breakpoints
            for(int j=0; j<n_clus; ++j)
            {
                // Count the number of soft clips around breakpoints
                for(int k=r_clip_idx[j]-3; k<r_clip_idx[j] + 3; ++k)
                {
                    G.end_clips[i] += D.rdstats[i].rclips[k];
                }
                for(int k=l_clip_idx[j]-3; k<l_clip_idx[j] + 3; ++k)
                {
                    G.start_clips[i] += D.rdstats[i].lclips[k];
                }
            }
        }
        G.rp_cnts[i] = 0;
        for(auto &rp : D.rdstats[i].readpairs)
        {
            // Because read pair counting has enough buffer size, we do not iterate through breakpoints to avoid double-counting
            int dx = (int) (start_pos - rp.positions.first);
            int dy = (int) (rp.positions.second - end_pos);
            if (dx > 75 && dx < 375 && dy > -75 && dy < 300)
            {
                G.rp_cnts[i] ++;
            }
        }
        
        G.all_cnts[i] = G.split_cnts[i] + G.rp_cnts[i] + (G.start_clips[i] + G.end_clips[i])/2.0;
        
        // TODO: cut-off values are arbitrary
        if (G.split_cnts[i] > 3 || (G.split_flag && G.split_cnts[i] > 0 && G.all_cnts[i] > 5) || (!G.split_flag && G.all_cnts[i] > 5))
        {
            G.gt[i] = 1;
            G.nonref_mask[i] = true;
            G.ac++;
            G.ns++;
            G.startclip_sum += G.start_clips[i];
            G.endclip_sum += G.end_clips[i];
            G.split_sum += G.split_cnts[i];
            G.rp_sum += G.rp_cnts[i];
        }
        else if (G.all_cnts[i] < 1)
        {
            G.gt[i] = 0;
            G.ns++;
        }
    }
    
    if (G.ac>0 && G.ns > G.ac)
        return true;
    else
        return false;
}

void Genotyper::call_cnv(sv &S, SvData& D, SvGeno &G, std::vector<SampleStat> &stats)
{
    bool split_flag = false;
    // First, try split read - best evidence
    if (find_consensus_split(S, D, G))
    {
        if (get_dup_cnts(S, D, G) && assign_dup_genotypes(S, D, G, stats))
        {
            G.split_flag = true;
            /*
            for(int i=0; i<n_sample;++i)
            {
                printf("%d\t%d\t%d\t%d\t%d\t%f\t%f\t%f\t%d\n", (int)G.split_cnts[i], (int)G.rp_cnts[i], (int)G.start_clips[i], (int)G.end_clips[i], (int)G.all_cnts[i], D.dps[0][i], D.dps[1][i], D.dps[2][i], G.gt[i]);
                
            }*/
            return;
        }
        else
        {
            G.split_flag = false;
        }
        split_flag = true;
    }

    D.clus_idx = -1;
    if (!split_flag && find_consensus_clip(S, D, G))
    {
        G.clip_flag = true;

        if (get_dup_cnts(S, D, G) && assign_dup_genotypes(S, D, G, stats))
        {
            return;
        }
    }
    
    // If no breakpoints identified by split reads nor by soft clips, just use the SV breakpoints as given
    if (!split_flag && !G.clip_flag)
    {
        D.vec_break_clusters.clear();
        D.clus_idx = -1;
        BreakCluster br;
        br.add_to_cluster(S.pos, S.end, 1);
        D.vec_break_clusters.push_back(br);

        if (get_dup_cnts(S, D, G) && assign_dup_genotypes(S, D, G, stats))
        {
            G.readpair_flag = true;
            /*
            for(int i=0; i<n_sample;++i)
            {
                printf("%d\t%d\t%d\t%d\t%d\t%f\t%f\t%f\t%d\n", (int)G.split_cnts[i], (int)G.rp_cnts[i], (int)G.start_clips[i], (int)G.end_clips[i], (int)G.all_cnts[i], D.dps[0][i], D.dps[1][i], D.dps[2][i], G.gt[i]);
                
            }*/
            return;
        }
    }
    // use pre_post filtering for depth-only
    if (!G.b_pre || !G.b_post) return;
    
    // Fit Gaussian mixture models with 1, 2, and 3 components, compare BIc
    std::vector< std::vector<double> > means = { {1.0}, {1.0, 1.5}, {1.0, 1.5, 2.0}, {1.0, 1.5, 2.0, 2.5}, {1.0, 1.5, 2.0, 2.5, 3.0}, {1.0, 1.5, 2.0, 2.5, 3.0, 3.5}, {1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0} };

    std::vector<int> dp_cn (n_sample, -1);
    std::vector<int> dp2_cn (n_sample, -1);

    int dp_ns = 0;
    int dp2_ns = 0;

	double best_dp_idx = 2;
	std::vector<double> &var_depth = D.dps[best_dp_idx];

    select_model_mask_1d(G.gmix, means, D.dps[best_dp_idx], G.sample_mask, G.MAX_P_OVERLAP);
    
    if (G.gmix.n_comp > 1 && G.gmix.r_ordered() )
    {
        // success
        G.dp_flag = true;

        for(int i=0; i<(int)n_sample; ++i)
        {
            if (!G.sample_mask[i]) continue;
            // TODO: arbitrary
            if (var_depth[i] > 0.7 && var_depth[i] < 1.25)
            {
                dp_cn[i] = 2;
                dp_ns++;
            }
            else if (var_depth[i] > 1.4 &&  var_depth[i] < 1.9)
            {
                dp_cn[i] = 3;
                dp_ns ++;
            }
            else if (var_depth[i] >= 1.9)
            {
                dp_cn[i] = round(var_depth[i] * 2.0);
                dp_ns++;
            }
        }
    }
    if (D.multi_dp)
    {
        // DP100 genotyping
        int dp2_idx = 3;

        select_model_mask_2d(G.gmix2, means, D.dps[dp2_idx], D.dps[dp2_idx+1], G.sample_mask, G.MAX_P_OVERLAP);

        // 2-D genotyping
        if (G.gmix2.n_comp>1 && G.gmix2.r_ordered() )
        {
            G.dp2_flag = true;
            for(int i=0; i<(int)n_sample; ++i)
            {
                if (!G.sample_mask[i]) continue;
//                dp2_cn[i] = G.gmix2.assign_copynumber(D.dps[dp2_idx][i], D.dps[dp2_idx+1][i]);

                double dp1 = D.dps[dp2_idx][i];
                double dp2 = D.dps[dp2_idx+1][i];
                
                if (dp1 > 0.8 && dp1 < 1.25 && dp2>0.8 && dp2<1.25)
                {
                    dp2_cn[i] = 2;
                    dp2_ns++;
                }
                else if (dp1>1.4 && dp2 > 1.4 &&  dp2 < 1.9 && dp2<1.9)
                {
                    dp2_cn[i] = 3;
                    dp2_ns ++;
                }
                else if (dp1 >= 1.9 && dp2 >= 1.9)
                {
                    dp2_cn[i] = round(var_depth[i] * 2.0);
                    dp2_ns++;
                }
            }
        }
    }

    if (G.dp_flag && G.dp2_flag)
    {
        for(int i=0; i<n_sample; ++i)
        {
            if (!G.sample_mask[i]) continue;
            if (dp_cn[i] <0 && dp2_cn[i]>=0)
                G.cn[i] = dp2_cn[i];
            else if (dp_cn[i] >= 0 && dp2_cn[i] < 0)
                G.cn[i] = dp_cn[i];
            else
                G.cn[i] = round((dp_cn[i] + dp2_cn[i])/2.0);
        }
    }
    else if (G.dp_flag)
    {
        for(int i=0; i<n_sample; ++i)
        {
            if (!G.sample_mask[i]) continue;
            G.cn[i] = dp_cn[i];
        }
    }
    else if (G.dp2_flag)
    {
        for(int i=0; i<n_sample; ++i)
        {
            if (!G.sample_mask[i]) continue;
            G.cn[i] = dp2_cn[i];
        }
    }
    // If clustered, filter false positive variants
    for (int i=0;i<n_sample; ++i)
    {
        if (!G.sample_mask[i]) continue;

        if (abs(D.dps[0][i] - G.dp_pre_mean) > 2.0*G.dp_pre_std || abs(D.dps[1][i] - G.dp_post_mean) > 2.0* G.dp_post_std || (G.cn[i] > 0 && abs(var_depth[i] - D.dps[1][i]) < 2.0*G.dp_post_std) || (G.cn[i] > 0 && abs(var_depth[i] - D.dps[0][i]) < 2.0*G.dp_pre_std))
        {
            G.cn[i] = -1;
        }
    }

    check_biallelic(G);
  
    G.ns = 0;
    G.ac = 0;
    for(int i=0; i<n_sample; ++i)
    {
        if (!G.sample_mask[i]) continue;
        if (G.cn[i] == 2)
        {
            G.gt[i] = 0; // 0/0
            G.ns++;
        }
        else if (G.cn[i] == 3)
        {
            G.gt[i] = 1; // 0/1
            G.ns++;
            G.ac++;
        }
        else if (G.cn[i] >3)
        {
            G.gt[i] = 2; // 1/1
            G.ns++;
            G.ac+=2;
        }
    }

    double callrate = (double)G.ns / G.n_effect;
    
	if ( (G.dp_flag || G.dp2_flag) && callrate>0.5 && G.ac > 0 && G.ac < G.ns*2)
    {
        /*
        for(int i=0; i<n_sample;++i)
        {
            printf("%d\t%d\t%d\t%d\t%d\t%f\t%f\t%f\t%d\n", (int)G.split_cnts[i], (int)G.rp_cnts[i], (int)G.start_clips[i], (int)G.end_clips[i], (int)G.all_cnts[i], D.dps[0][i], D.dps[1][i], D.dps[2][i], G.gt[i]);
            
        }*/
        G.b_pass = true;
    }
}

void Genotyper::check_biallelic(SvGeno &G)
{
    int max_cn = 2;
    int min_cn = 2;
    for(int i=0; i<n_sample; ++i)
    {
        if (!G.sample_mask[i]) continue;
        if (G.cn[i] > max_cn)
        {
            max_cn = G.cn[i];
        }
        if (G.cn[i] >=0 && G.cn[i] < min_cn)
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

}
