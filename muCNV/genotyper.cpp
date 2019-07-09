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
    
    MAX_P_OVERLAP = 1.0;

    gt.resize(n_sample, -1);
    cn.resize(n_sample, -1);
    
    rp_gt.resize(n_sample, -1);
    rp_cn.resize(n_sample, -1);
    
    clip_gt.resize(n_sample, -1);
    clip_cn.resize(n_sample, -1);
    
    start_clips.resize(n_sample, 0);
    end_clips.resize(n_sample, 0);
    
    start_rps.resize(n_sample, 0);
    end_rps.resize(n_sample, 0);

}

void SvGeno::reset()
{
    std::fill(gt.begin(), gt.end(), -1);
    std::fill(cn.begin(), cn.end(), -1);
    
    std::fill(rp_gt.begin(), rp_gt.end(), -1);
    std::fill(rp_cn.begin(), rp_cn.end(), -1);
    
    std::fill(clip_gt.begin(), clip_gt.end(), -1);
    std::fill(clip_cn.begin(), clip_cn.end(), -1);
    
    std::fill(start_clips.begin(), start_clips.end(), 0);
    std::fill(end_clips.begin(), end_clips.end(), 0);
    std::fill(start_rps.begin(), start_rps.end(), 0);
    std::fill(end_rps.begin(), end_rps.end(), 0);
    
    ns = 0;
    ac = 0;

    b_biallelic = false;
    b_pass = false;
    dp_flag = false;
    dp2_flag = false;
    pd_flag = false;
    read_flag = false;
	rp_geno_flag = false;
    clip_flag = false;
	clip_geno_flag = false;
	
    
    dp_pre_mean = 1.0;
    dp_pre_std = 0.1;
    dp_post_mean = 1.0;
    dp_post_std = 0.1;
    b_pre = false;
    b_post = false;
    info = "";
    rp_pos = -1;
    rp_end = -1;
    clip_pos = -1;
    clip_end = -1;
}

SvData::SvData(int n)
{
    n_sample = n;
    rdstats.resize(n);
    prepost_dp.resize(n_sample, 1);
    
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
    
    // ??
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
    
    std::fill(prepost_dp.begin(), prepost_dp.end(), 1);
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

void Genotyper::select_model(GaussianMixture &ret_gmix, std::vector< std::vector<double> > &means, std::vector<double> &x, double MAX_P_OVERLAP)
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
    //        gmix.print();
        }
       // if (gmix.bic < best_bic)
        if (gmix.bic < best_bic && gmix.p_overlap < MAX_P_OVERLAP)
        {
            best_bic = gmix.bic;
            ret_gmix = gmix; // assignment
        }
    }
    return;
}

void Genotyper::select_model_mask(GaussianMixture &ret_gmix, std::vector< std::vector<double> > &means, std::vector<double> &x, std::vector<bool> &mask, double MAX_P_OVERLAP)
{
    double best_bic = DBL_MAX;
    
    // number of models
    for(int m=0; m<(int)means.size(); ++m)
    {
        std::vector<double> s (means[m].size(), 0.1);
        GaussianMixture gmix(means[m], s);
        if (b_kmeans)  //should not be true for this version (04/13/2019)
        {
            gmix.KM(x, b_mahalanobis);
        }
        else
        {
            gmix.EM_select(x, mask); // fit mixture model
                        //        gmix.print();
        }
        // if (gmix.bic < best_bic)
        if (gmix.bic < best_bic && gmix.p_overlap < MAX_P_OVERLAP)
        {
            best_bic = gmix.bic;
            ret_gmix = gmix; // assignment
        }
    }
    return;
}


void Genotyper::select_model(GaussianMixture2 &ret_gmix2, std::vector< std::vector<double> > &means, std::vector<double> &x, std::vector<double> &y, double MAX_P_OVERLAP)
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
        //if (gmix2.bic < best_bic )
        if (gmix2.bic < best_bic && gmix2.p_overlap < MAX_P_OVERLAP )
        {
            best_bic = gmix2.bic;
            ret_gmix2 = gmix2; // assignment (copy operation)
        }
    }
    return;
}

void Genotyper::call(sv &S, SvData &D, SvGeno &G, bool bk, bool bm, std::vector<SampleStat> &stats)
{
    n_sample = D.n_sample;

    b_kmeans = bk;
    b_mahalanobis = bm;

    get_prepost_stat(D, G);

	if (!G.b_pre || !G.b_post)
		return;

    if (S.svtype == DEL)
	{
        call_deletion(S, D, G);
    }
    else if (S.svtype == DUP || S.svtype == CNV)
    {
        call_cnv(S, D, G);
    }
    else if (S.svtype == INV)
    {
        call_inversion(S, D, G, stats);
    }
    /*
    else if (S.svtype == INS)
    {
        call_insertion(S, D, G);
    }
    */  // TODO: insertion calling
}

void Genotyper::call_inversion(sv &S, SvData &D, SvGeno &G, std::vector<SampleStat> &stats)
{
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

bool Genotyper::find_consensus_rp(sv &S, SvData &D, int pairstr, int &start_peak, int &end_peak)
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

bool Genotyper::find_consensus_clip(sv &S, SvData &D, int pairstr, int &l_peak, int &r_peak)
{
	int N_buf = 100;
	
	if (S.len < 200) 
	{
		N_buf = (int)S.len/2;
	//	fprintf(stderr, "first %d\tsecond %d\n", N_first, N_second);
	}

/*
	fprintf(stderr, "RCLIP\tLCLIP\n");
	for(int i=0;i<D.all_rclips.size(); ++i)
	{
		fprintf(stderr, "%d: %d\t%d\n", i, D.all_rclips[i], D.all_lclips[i]);
	}
	*/

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


bool Genotyper::find_consensus_clip_inv(sv &S, SvData &D, int &l_start, int &l_end, int &r_start, int &r_end)
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

void Genotyper::call_deletion(sv &S, SvData &D, SvGeno &G)
{
    // fit gaussian mixture models with 1, 2, and 3 components, compare bic
    std::vector< std::vector<double> > means = { {1.0}, {1.0, 0.5}, {1.0, 0.5, 0.0}};
    G.b_biallelic = true;

    int start_peak = -1;
    int end_peak = -1;
    int pairstr = 1;
    
   if (find_consensus_rp(S, D, pairstr, start_peak, end_peak) )
    {
        G.rp_pos = S.pos + start_peak*10 - 500;
        G.rp_end = S.end + end_peak*10 - 1500;

		if (G.rp_pos < G.rp_end)
		{
			G.read_flag = true;
			for(int i=0; i<n_sample; ++i)
			{
				G.start_rps[i]= D.rdstats[i].rp_seq[pairstr][start_peak] + D.rdstats[i].rp_seq[pairstr][start_peak-1] + D.rdstats[i].rp_seq[pairstr][start_peak+1];
				
				G.end_rps[i] = D.rdstats[i].rp_seq[pairstr][end_peak] + D.rdstats[i].rp_seq[pairstr][end_peak-1] + D.rdstats[i].rp_seq[pairstr][end_peak+1];
			
				if ( G.start_rps[i] > 5 && G.end_rps[i] > 5 && D.dps[2][i] < 0.15)
				{
					G.rp_gt[i] = 2;
					G.rp_cn[i] = 0;
					G.rp_geno_flag = true;
				}
				else if ( G.start_rps[i] > 3 && G.end_rps[i] > 3 && D.dps[2][i] < 0.75)
				{
					G.rp_gt[i] = 1;
					G.rp_cn[i] = 1;
					G.rp_geno_flag = true;
				}
				else
				{
					if (D.dps[2][i]>0.85 && D.dps[2][i] <= 1.5)
					{
						G.rp_gt[i] = 0;
						G.rp_cn[i] = 2;
					}
				}
			}
		}
    }

    int l_clip = -1;
    int r_clip = -1;
    
    if (find_consensus_clip(S, D, pairstr, l_clip, r_clip))
    {
        int start_clip = r_clip;
        int end_clip = l_clip;
        
        if ((S.pos + start_clip - 100) < (S.end + end_clip - 300) && start_clip > 0 && start_clip<199 && end_clip>200 && end_clip<399)
        {
            G.clip_pos = S.pos + start_clip - 100;
            G.clip_end = S.end + end_clip - 300;
			G.clip_flag = true;

            for(int i=0; i<n_sample; ++i)
            {
                G.start_clips[i] = D.rdstats[i].rclips[r_clip] + D.rdstats[i].rclips[r_clip-1] + D.rdstats[i].rclips[r_clip+1];
                G.end_clips[i] = D.rdstats[i].lclips[l_clip] + D.rdstats[i].lclips[l_clip-1] + D.rdstats[i].lclips[l_clip+1];
                
                if ( G.start_clips[i] >= 5 && G.end_clips[i] >=5 && D.dps[2][i] < 0.15 )
                {
                    G.clip_gt[i] = 2;
                    G.clip_cn[i] = 0;
                    G.clip_geno_flag = true;}
                else if ( G.start_clips[i] >= 2 && G.end_clips[i] >=2 && D.dps[2][i] < 0.75 )
                {
                    G.clip_gt[i] = 1;
                    G.clip_cn[i] = 1;
                    G.clip_geno_flag = true;
                }
                else
                {
                    if ( D.dps[2][i] >0.8 && D.dps[2][i] <= 1.5)
                    {
                        G.clip_gt[i] = 0;
                        G.clip_cn[i] = 2;
                    }
                }
            }
        }
    }
    
	double best_dp_idx = 2;
	std::vector<double> &var_depth = D.dps[best_dp_idx];

	if (G.b_pre && G.b_post)
	{
		// Depth-based clustering genotyping
		select_model(G.gmix, means, D.dps[best_dp_idx], G.MAX_P_OVERLAP);


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
		if (D.multi_dp) //dp2 has more than 2 vectors
		{
			// dp100 genotyping
			select_model(G.gmix2, means, D.dps[dp2_idx], D.dps[dp2_idx+1], G.MAX_P_OVERLAP);

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
	} // G.b_pre && G.b_post

    // Let's use peripheral depth to identify 'false positve' only
    if (G.b_pre && G.b_post)
    {
        for (int i=0;i<n_sample; ++i)
        {
          	if (abs(D.dps[0][i] - G.dp_pre_mean) > 2.0*G.dp_pre_std || abs(D.dps[1][i] - G.dp_post_mean) > 2.0* G.dp_post_std || abs(var_depth[i] - D.dps[1][i]) < 0.3 || abs(var_depth[i] - D.dps[0][i]) < 0.3) 
            {
				G.gt[i] = -1;
                G.cn[i] = -1;
                G.rp_cn[i] = -1;
                G.clip_cn[i] = -1;
            }

			// fill in missing values
			if (G.cn[i] < 0 && D.prepost_dp[i] > 0 && abs(var_depth[i] - D.prepost_dp[i]) < 0.15)
			{
				G.cn[i] = 2;
				G.gt[i] = 0;
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
        int sum_cn = 0;
        int cnt_cn = 0;
        // get consensus

        if (G.rp_geno_flag && G.rp_cn[i] >= 0)
        {
            sum_cn += G.rp_cn[i];
            cnt_cn ++;
        }
        if (G.clip_geno_flag && G.clip_cn[i] >= 0)
        {
            sum_cn += G.clip_cn[i];
            cnt_cn ++;
        }
        if (cnt_cn>0)
        {
            G.cn[i] = round(sum_cn / (double)cnt_cn);
            G.gt[i] = 2-G.cn[i];
        }
        
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

/*
    if (n_ref == 0)
    {
        G.b_pass = false;
        return;
    }

    if (n_hom >0 )
    {
        hom_dp /= (double)n_hom;
    }

    if (n_het >0 )
    {
        het_dp /= (double)n_het;
    }

    if (n_ref >0 )
    {
        ref_dp /= (double)n_ref;
    }

    if (het_dp < 0.3 || hom_dp > 0.2)
    {
        G.b_pass = false;
        return;
    }
	*/

    double callrate = (double)G.ns / n_sample;

//	fprintf(stderr, " callrate %f , ns %d, ac %d\n", callrate, G.ns, G.ac); 
    if ( (((G.dp_flag  || G.dp2_flag) && G.clip_flag && G.read_flag) || G.clip_geno_flag || G.rp_geno_flag ) &&  callrate>0.3 && G.ac > 0 && G.ac < G.ns*2)
        G.b_pass = true;

    // Excessive heterozygosity (basically all-het case)
    if ( n_het > ((double)G.ns * 0.75) && n_hom < ((double)0.05 * G.ns) )
        G.b_pass = false;
}

void Genotyper::call_cnv(sv &S, SvData& D, SvGeno &G)
{
    int start_peak = -1;
    int end_peak = -1;
    int pairstr = 2;
    
    if (find_consensus_rp(S, D, pairstr, start_peak, end_peak)  )
    {
        G.rp_pos = S.pos + start_peak*10 - 500;
        G.rp_end = S.end + end_peak*10 - 1500;
        
		if (G.rp_pos < G.rp_end)
		{
			G.read_flag = true;
			for(int i=0; i<n_sample; ++i)
			{
				G.start_rps[i] = D.rdstats[i].rp_seq[pairstr][start_peak] + D.rdstats[i].rp_seq[pairstr][start_peak-1] + D.rdstats[i].rp_seq[pairstr][start_peak+1];
				G.end_rps[i] = D.rdstats[i].rp_seq[pairstr][end_peak] + D.rdstats[i].rp_seq[pairstr][end_peak-1] + D.rdstats[i].rp_seq[pairstr][end_peak+1];
				
				if (G.start_rps[i] > 5 && G.end_rps[i] > 5 && D.dps[2][i] > 1.85)
				{
					G.rp_cn[i] = round(D.dps[2][i]*2.0);
					G.rp_geno_flag = true;
				}
				else if (G.start_rps[i] > 3 && G.end_rps[i] > 3 && D.dps[2][i] > 1.3 && D.dps[2][i]<1.8 )
				{
					G.rp_cn[i] = 3;
					G.rp_geno_flag = true;
				}
				else if (G.start_rps[i] < 3 && G.end_rps[i] < 3 && D.dps[2][i]>0.85 && D.dps[2][i] <= 1.3)
				{
					G.rp_gt[i] = 0;
					G.rp_cn[i] = 2;
				}
			}
		}
    }
    
    int l_clip = -1;
    int r_clip = -1;
    
    if (find_consensus_clip(S, D, pairstr, l_clip, r_clip))
    {
        int start_clip = l_clip;
        int end_clip = r_clip;
        
        if ((S.pos + start_clip - 100) < (S.end + end_clip - 300) && start_clip > 0 && start_clip<199 && end_clip>200 && end_clip<399)
        {
            G.clip_pos = S.pos + start_clip - 100;
            G.clip_end = S.end + end_clip - 300;
			G.clip_flag = true;

            for(int i=0; i<n_sample; ++i)
            {
                G.start_clips[i] = D.rdstats[i].lclips[l_clip] + D.rdstats[i].lclips[l_clip-1] + D.rdstats[i].lclips[l_clip+1];
                G.end_clips[i] = D.rdstats[i].rclips[r_clip] + D.rdstats[i].rclips[r_clip-1] + D.rdstats[i].rclips[r_clip+1];
                if (G.start_clips[i] >=5 && G.end_clips[i] >= 5 && D.dps[2][i] > 1.85 )
                {
                    G.clip_gt[i] = 2;
                    G.clip_cn[i] = round(D.dps[2][i] * 2.0);
                    G.clip_geno_flag = true;
                }
                else if (G.start_clips[i] >=2 && G.end_clips[i] >=2 && D.dps[2][i] < 1.8 && D.dps[2][i] >1.3 )
                {
                    G.clip_gt[i] = 1;
                    G.clip_cn[i] = 3;
                    G.clip_geno_flag = true;
                }
                else if (G.start_clips[i] < 2 && G.end_clips[i] < 2 && D.dps[2][i]>0.8 && D.dps[2][i] <= 1.3)
				{
					G.clip_gt[i] = 0;
					G.clip_cn[i] = 2;
				}
            }
        }
    }
    
    // Fit Gaussian mixture models with 1, 2, and 3 components, compare BIC
    std::vector< std::vector<double> > means = { {1.0}, {1.0, 1.5}, {1.0, 1.5, 2.0}, {1.0, 1.5, 2.0, 2.5}, {1.0, 1.5, 2.0, 2.5, 3.0}, {1.0, 1.5, 2.0, 2.5, 3.0, 3.5}, {1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0} };

    std::vector<int> dp_cn (n_sample, -1);
    std::vector<int> dp2_cn (n_sample, -1);

    int dp_ns = 0;
    int dp2_ns = 0;

	double best_dp_idx = 2;
	std::vector<double> &var_depth = D.dps[best_dp_idx];
	if (G.b_pre && G.b_post)
	{

		select_model(G.gmix, means, D.dps[best_dp_idx], G.MAX_P_OVERLAP);


		if (G.gmix.n_comp > 1 && G.gmix.r_ordered() )
		{
			// success
			G.dp_flag = true;

			for(int i=0; i<(int)n_sample; ++i)
			{
				/*
				dp_cn[i] = G.gmix.assign_copynumber(var_depth[i]);
				if (dp_cn[i] >= 2)
				{
					dp_ns++;
				}
				*/

				// TODO: arbitrary
				if (var_depth[i] > 0.7 && var_depth[i] < 1.25)
				{
					dp_cn[i] = 2;
					dp_ns++;
				}
				else if (var_depth[i] > 1.4 &&  var_depth[i] < 1.8)
				{
					dp_cn[i] = 3;
					dp_ns ++;
				}
				else if (var_depth[i] >= 1.85)
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

			select_model(G.gmix2, means, D.dps[dp2_idx], D.dps[dp2_idx+1], G.MAX_P_OVERLAP);

			// 2-D genotyping
			if (G.gmix2.n_comp>1 && G.gmix2.r_ordered() )
			{
				G.dp2_flag = true;
				for(int i=0; i<(int)n_sample; ++i)
				{
	//                dp2_cn[i] = G.gmix2.assign_copynumber(D.dps[dp2_idx][i], D.dps[dp2_idx+1][i]);

					double dp1 = D.dps[dp2_idx][i];
					double dp2 = D.dps[dp2_idx+1][i];
					
					if (dp1 > 0.8 && dp1 < 1.25 && dp2>0.8 && dp2<1.25)
					{
						dp2_cn[i] = 2;
						dp2_ns++;
					}
					else if (dp1>1.4 && dp2 > 1.4 &&  dp2 < 1.8 && dp2<1.8)
					{
						dp2_cn[i] = 3;
						dp2_ns ++;
					}
					else if (dp1 > 1.85 && dp2 > 1.85)
					{
						dp2_cn[i] = round(var_depth[i] * 2.0);
						dp2_ns++;
					}
				}
			}
		}

		for(int i=0; i<n_sample; ++i)
		{
			if (dp_cn[i] <0 && dp2_cn[i]>=0)
				G.cn[i] = dp2_cn[i];
			else if (dp_cn[i] >= 0 && dp2_cn[i] < 0)
				G.cn[i] = dp_cn[i];
			else
				G.cn[i] = round((dp_cn[i] + dp2_cn[i])/2.0);
		}
	} // G.b_pre && G.b_post

    if (G.b_pre && G.b_post)
    {
        // If clustered, filter false positive variants
        for (int i=0;i<n_sample; ++i)
        {
          	if (abs(D.dps[0][i] - G.dp_pre_mean) > 2.0*G.dp_pre_std || abs(D.dps[1][i] - G.dp_post_mean) > 2.0* G.dp_post_std || abs(var_depth[i] - D.dps[1][i]) < 0.3 || abs(var_depth[i] - D.dps[0][i]) < 0.3)
            {
                G.cn[i] = -1;
                G.rp_cn[i] = -1;
                G.clip_cn[i] = -1;
            }
   			if (G.cn[i] < 0 && D.prepost_dp[i] > 0 && abs(var_depth[i] - D.prepost_dp[i]) < 0.15)
			{
				G.cn[i] = 2;
				G.gt[i] = 0;
			}
        }
    }

/*
    if (!G.pd_flag)
    {
        if (G.dp_flag && G.gmix.p_overlap > G.MAX_P_OVERLAP)
        {
            G.dp_flag = false;
        }
        if (G.dp2_flag && G.gmix2.p_overlap > G.MAX_P_OVERLAP)
        {
            G.dp2_flag = false;
        }
    }
    */

    int max_cn = 2;
    int min_cn = 2;
    for(int i=0; i<n_sample; ++i)
    {
        int sum_cn = 0;
        int cnt_cn = 0;
        // get consensus
        
        if (G.rp_geno_flag && G.rp_cn[i] >= 0)
        {
            sum_cn += G.rp_cn[i];
            cnt_cn ++;
        }
        if (G.clip_geno_flag && G.clip_cn[i] >= 0)
        {
            sum_cn += G.clip_cn[i];
            cnt_cn ++;
        }
        if (cnt_cn>0)
        {
            G.cn[i] = round(sum_cn / (double)cnt_cn);
        }
        
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

    G.ns = 0;
    G.ac = 0;
    for(int i=0; i<n_sample; ++i)
    {
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

    double callrate = (double)G.ns / n_sample;
    
	if ( (((G.dp_flag || G.dp2_flag) && G.clip_flag && G.read_flag) || G.clip_geno_flag || G.rp_geno_flag) && callrate>0.3 && G.ac > 0 && G.ac < G.ns*2)
        G.b_pass = true;

}
