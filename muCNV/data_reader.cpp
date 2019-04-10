//
//  multi_pileup.cpp
//  muCNV
//
//  Created by Goo Jun on 11/25/18.
//  Copyright © 2018 Goo Jun. All rights reserved.
//

#include "data_reader.h"
#include <math.h>
#include <algorithm>

bool ReadStat::del_support()
{
    int rp_cnt = n_pre_FR + n_post_FR;
    int sp_cnt = n_pre_split_in + n_post_split_in;
    int clip_cnt = n_pre_clip_in + n_post_clip_in;
    
    if (rp_cnt >=6   && rp_cnt + sp_cnt + clip_cnt > 10)
        return true;
		/*
    else if (sp_cnt > 6 && rp_cnt + sp_cnt + clip_cnt > 10)
        return true;
    else if (sp_cnt + rp_cnt + clip_cnt > 15)
        return true;
		*/
    
    return false;
}
bool ReadStat::dup_support()
{
    int rp_cnt = n_pre_RF + n_post_RF;
    int sp_cnt = n_pre_split_out + n_post_split_out;
    int clip_cnt = n_pre_clip_out + n_post_clip_out;
    
    if (rp_cnt >= 4  && rp_cnt + sp_cnt + clip_cnt > 10)
        return true;
		/*
    else if (sp_cnt > 6 && rp_cnt + sp_cnt + clip_cnt > 10)
        return true;
    else if (sp_cnt + rp_cnt + clip_cnt > 15)
        return true;
		*/
    
    return false;
}

bool ReadStat::inv_support()
{
    int rp_cnt = n_pre_FF + n_post_RR;
    
    // inversion can have both in and out clipping
    int sp_cnt = n_pre_split_out + n_post_split_out + n_pre_split_in + n_post_split_in;
    int clip_cnt = n_pre_clip_in + n_pre_clip_out + n_post_clip_in + n_post_clip_out;
    
    if (rp_cnt >= 5 && rp_cnt+sp_cnt+clip_cnt>=15)
        return true;

    return false;
}

bool ReadStat::ins_support()
{
    if (n_pre_INS > 3 && n_pre_INS + n_pre_clip_in + n_post_clip_in + n_pre_sp_missing > 10)
        return true; 
    if (n_pre_INS + n_pre_rp_missing + n_pre_clip_in + n_post_clip_in + n_pre_sp_missing > 10 && !del_support() && !dup_support() && !inv_support())
        return true;
    return false;
}

void median_filter(uint16_t* D, uint16_t* D_filt, int n_sample, int n_dp)
{
    for(int i=0; i<n_sample; ++i)
    {
        for(int j=2; j<n_dp-2; ++j) // first two and last two will be padded later
        {
            std::vector<uint16_t> buf (5,0);
            for(int k=0; k<5; ++k)
            {
                buf[k] = D[(j-2+k)*n_sample + i];
            }
            std::nth_element (buf.begin(), buf.begin()+2, buf.end());
            D_filt[j*n_sample+i] = buf[2];
        }
        D_filt[i] = D_filt[i+n_sample] = D_filt[i+n_sample*2]; // padding
        D_filt[(n_dp-1)*n_sample+i] = D_filt[(n_dp-2)*n_sample+i] = D_filt[(n_dp-3)*n_sample+i];
    }
	/*
	for(int i=0; i<n_dp;++i)
	{
		for(int j=0; j<n_sample;++j)
		{
			fprintf(stdout, "%.1f,%.1f\t", D[i*n_sample + j]/23.0, D_filt[i*n_sample +j]/32.0);
		}
		fprintf(stdout, "\n");
		fflush(stdout);
	}
	*/
}

int DataReader::load(std::vector<string>& base_names, std::vector<SampleStat> &stats, GcContent &gc, int chr)
{
    // base_names: list of base names for pileup/var/idx triples
    // return value: total number of samples
    
    n_sample_total = 0;
    n_pileup = (int) base_names.size();
	multi_idx.resize(n_pileup);

    int prev_n_var = 0;
    pileups.resize(n_pileup);
    var_files.resize(n_pileup);
    idx_files.resize(n_pileup);
    n_samples.resize(n_pileup);
    
    int idx_cnt = 1;
    chr_idx_rp.resize(gc.num_chr+1);
	if (chr>0)
	{
		for(int i=1;i<=gc.num_chr; ++i)
		{
			chr_idx_rp[i] = -1; // make everything invalid 
		}
		chr_idx_rp[chr] = idx_cnt;
		idx_cnt += ceil((double)gc.chr_size[chr] / 10000.0) ;
	}
	else
	{
		for(int i=1;i<=gc.num_chr; ++i)
		{
			chr_idx_rp[i] = idx_cnt; // chr_idx[1] contains the array index of pileup index of chr 1
			idx_cnt += ceil((double)gc.chr_size[i] / 10000.0) ;
		}
	}
    
    DMSG(idx_cnt << " indices should be in index file");
 
    for(int i=0; i<n_pileup ; ++i)
    {
        string pileup_name = base_names[i] + ".pileup";
        string var_name = base_names[i] + ".var";
        string idx_name = base_names[i] + ".idx";
        
        pileups[i].open(pileup_name, std::ios::in | std::ios::binary);
        var_files[i].open(var_name, std::ios::in | std::ios::binary);
        idx_files[i].open(idx_name, std::ios::in | std::ios::binary);
        
        // number of samples in each pileup
        pileups[i].read_int32(n_samples[i]);
        int32_t tmp;
        var_files[i].read_int32(tmp);
        if (n_samples[i] != tmp)
        {
            std::cerr << "Error! number of samples in pileup and vars do not match: " << base_names[i] << std::endl;
            exit(1);
        }
        n_sample_total += n_samples[i] ;

        for(int j=0; j<n_samples[i]; ++j)
        {
			char buf[256];
            pileups[i].read_sample_id(buf);
			sample_ids.push_back(std::string(buf));
        }
        
        for(int j=0; j<n_samples[i]; ++j)
        {
            SampleStat s;
            pileups[i].read_sample_stat(s);
            stats.push_back(s);
        }
        
        for(int j=0; j<n_samples[i]; ++j)
        {
            std::vector<double> gc_f;
            pileups[i].read_gc_factor(gc_f, gc.num_bin);

            /*
            std::vector<double> gcf10;

            //interpolate into 10-sub-bins
            gcf10.resize(gc_f.size()*10);

            for(int k=0; k< (int)gc_f.size()-1; ++k)
            {
                for(int l=0; l<10; ++l)
                {
                    // interpolate 1/gc_factor using linear interfolation and then inverse it
                    gcf10[k*10+l] = 1.0/((1.0/gc_f[k+1] - 1.0/gc_f[k]) *  (l/10.0) + 1.0/gc_f[k]);
                }
            }
            gcf10[(gc_f.size()-1) * 10] = gc_f[gc_f.size()-1];
            for(int k= (int) (gc_f.size()-1) * 10; k < gcf10.size(); ++k)
            {
                gcf10[k] = gc_f[gc_f.size()-1];
            }
                
            gc_factors.push_back(gcf10);
            */
            gc_factors.push_back(gc_f);
        }
        
        // number of variants
        var_files[i].read_int32(n_var);
        if (i>0 && n_var !=prev_n_var)
        {
            std::cerr << "Error! number of  variants do not match between " << base_names[i-1] << " and " << base_names[i] << std::endl;
            exit(1);
        }
        prev_n_var = n_var;
        
        multi_idx[i] = new uint64_t[idx_cnt];
        idx_files[i].read_uint64_multi(multi_idx[i], idx_cnt);
        idx_files[i].close(); // index file content is now on memory
        
        // idx offset?
    }

    chr_bytepos_dp100.resize(n_pileup);
    for(int i=0; i<n_pileup; ++i)
    {
        chr_bytepos_dp100[i].resize(gc.num_chr +1);
        chr_bytepos_dp100[i][0]= 0;
        chr_bytepos_dp100[i][1] = multi_idx[i][0];

		if (chr>0)
		{
			for(int c=2; c<=gc.num_chr; ++c)
			{
				chr_bytepos_dp100[i][chr] = -1; // make everything out of range
			}
			chr_bytepos_dp100[i][chr] = multi_idx[i][0]; // just this chromosome
		}
		else
		{
			for(int c=2; c<=gc.num_chr; ++c)
			{
				chr_bytepos_dp100[i][c] = chr_bytepos_dp100[i][c-1] + (gc.n_interval[c-1] * n_samples[i] * sizeof(uint16_t));
			}
		}
    }
    return n_sample_total;
}

int DataReader::read_depth100(sv& curr_sv, std::vector< std::vector<double> > &dvec_dp, GcContent& gc, bool b_dumpstat)
{
    // this information is not useful when sv length is short
    // process only for >200bp SVs (or at include least two full 100-bp intervals)
    
    bool b_medfilt = true; // whether to do median filtering or not
    
    int sample_idx = 0;
    
    int pre_start = (curr_sv.pos - 1500);
    int post_end = (curr_sv.end + 1500);
    
    if (pre_start<1) pre_start = 1 ;
    if (post_end > (int)gc.chr_size[curr_sv.chrnum]) post_end = gc.chr_size[curr_sv.chrnum];

    int idx_pre_start = pre_start / 100;
    int idx_post_end = (post_end/100);
    
    int idx_start = (curr_sv.pos / 100);
    int idx_end = (curr_sv.end / 100);
    
    int n_pre = idx_start - idx_pre_start;
    int n_post = idx_post_end - idx_end;
    
    // Read 'before start' depth and 'after end' depth
    // TODO: can be merged, or reordered to pre/center/post order for faster processing
    for(int i=0; i<n_pileup; ++i)
    {
        if (n_pre>2)
        {
            uint64_t start_byte = chr_bytepos_dp100[i][curr_sv.chrnum];
            start_byte += idx_pre_start * n_samples[i] * sizeof(uint16_t);

            int n_pre_by_sample = n_pre * n_samples[i];
            uint16_t *D = new uint16_t[n_pre_by_sample];
            
            pileups[i].seekg(start_byte);
            pileups[i].read_depth(D, n_pre_by_sample);
            
            std::vector<double> dp_sum (n_samples[i], 0);
            
            for(int j=0; j<n_pre; ++j)
            {
                for(int k=0; k<n_samples[i]; ++k)
                {
                   dp_sum[k] += correct_gc(gc, sample_idx+k, (double)D[j*n_samples[i] + k], curr_sv.chrnum, (pre_start + 50 + j*100) );
                }
            }
            for(int k=0; k<n_samples[i]; ++k)
            {
                dvec_dp[0][k + sample_idx] = (double) dp_sum[k] / n_pre / 32.0;
            }
            delete [] D;
        }
        if (n_post>2)
        {
            uint64_t start_byte = chr_bytepos_dp100[i][curr_sv.chrnum];
            start_byte += (idx_end + 1) * n_samples[i] * sizeof(uint16_t);

            int n_post_by_sample = n_post * n_samples[i];
            uint16_t *D = new uint16_t[n_post_by_sample];
            
            pileups[i].seekg(start_byte);
            pileups[i].read_depth(D, n_post_by_sample);
            
            std::vector<double> dp_sum (n_samples[i], 0);
            
            for(int j=0; j<n_post; ++j)
            {
                for(int k=0; k<n_samples[i]; ++k)
                {
                    dp_sum[k] += correct_gc(gc, sample_idx+k, (double)D[j*n_samples[i] + k], curr_sv.chrnum, (curr_sv.end + 150 + j*100) );
                }
            }
            for(int k=0; k<n_samples[i]; ++k)
            {
                dvec_dp[1][k + sample_idx] = (double) dp_sum[k] / n_post / 32.0;
            }
            delete [] D;
        }
        

        sample_idx += n_samples[i];
    }
    if (n_pre<=2)
    {
        for(int k=0; k<n_sample_total; ++k)
        {
            dvec_dp[0][k] = dvec_dp[1][k];
        }
    }
    else if (n_post <=2)
    {
        for(int k=0; k<n_sample_total; ++k)
        {
            dvec_dp[1][k] = dvec_dp[0][k];
        }
    }
    sample_idx = 0;
    
    int n_dp = idx_end - idx_start - 1;
    
    if (n_dp < 2)
    {
        for(int i=0; i<n_sample_total; ++i)
        {
            dvec_dp[2][i] = correct_gc(gc, i, dvec_dp[2][i], curr_sv.chrnum, (curr_sv.pos + curr_sv.end)/2);
        }
        return 0;
    }
    else if (n_dp <= 10)
        b_medfilt = false; // no median filtering
    
    for(int k=0; k<2; ++k)
    {
        std::vector<double> x(n_sample_total, 0);
        dvec_dp.push_back(x);
    }
    /*
    if (n_dp <= 40)
    {
        for(int k=0; k<2; ++k)
        {
            std::vector<double> x(n_sample_total, 0);
            dvec_dp.push_back(x);
        }
    }
    else
    {
        for(int k=0; k<4; ++k)
        {
            std::vector<double> x(n_sample_total, 0);
            dvec_dp.push_back(x);
        }
    }
     */

	if (n_dp <= 200)
	{
		for(int i=0; i<n_pileup; ++i)
		{
			uint64_t start_byte = chr_bytepos_dp100[i][curr_sv.chrnum];
			start_byte += (idx_start+1) * n_samples[i] * sizeof(uint16_t);
			
			int n_dp_by_sample = n_samples[i] * n_dp;
			
			uint16_t *D = new uint16_t[n_dp_by_sample];
			uint16_t *D_filt = new uint16_t[n_dp_by_sample];
			pileups[i].seekg(start_byte);
			pileups[i].read_depth(D, n_dp_by_sample);
            
			std::vector<unsigned> gc_sum (n_samples[i], 0);
						
			for(int j=0; j<n_dp; ++j)
			{
				for(int k=0; k<n_samples[i]; ++k)
				{
					double val = correct_gc(gc, sample_idx+k, (double)D[j*n_samples[i] + k], curr_sv.chrnum, (curr_sv.pos + 150 + j*100) );
                    gc_sum[k] += val;
                    D[j*n_samples[i] + k] = (uint16_t) val;
				}
			}

			if (b_medfilt)
				median_filter(D, D_filt, n_samples[i], n_dp);

			for(int k=0; k<n_samples[i]; ++k)
			{
				dvec_dp[2][k + sample_idx] = (double) gc_sum[k] / n_dp / 32.0;
			}

			std::vector<int> seg_starts;
			std::vector<int> seg_ends;

			if (n_dp <= 10)
			{
				seg_starts.push_back(0);
				seg_ends.push_back((int)n_dp/2);
				seg_starts.push_back(seg_ends[0]);
				seg_ends.push_back(n_dp);
			}
			else
			{
                seg_starts.push_back(2);
                seg_ends.push_back((int)n_dp/2);
                seg_starts.push_back(seg_ends[0]);
                seg_ends.push_back(n_dp-2);
                /*
                seg_starts.push_back(0);
                seg_ends.push_back((int)n_dp/4);
                seg_starts.push_back(seg_ends[0]);
                seg_ends.push_back((int)n_dp/2);
                seg_starts.push_back(seg_ends[1]);
                seg_ends.push_back((int)(n_dp*0.75));
                seg_starts.push_back(seg_ends[2]);
                seg_ends.push_back(n_dp);
                 */
			}
			for(int k=0; k<(int)seg_starts.size(); ++k)
			{
				std::vector<unsigned> dp_sum (n_samples[i], 0);
				std::vector<unsigned> dp_cnt (n_samples[i], 0);
				
				for(int j=seg_starts[k]; j<seg_ends[k]; ++j)
				{
					for(int m=0; m<n_samples[i]; ++m)
					{
						if (b_medfilt)
							dp_sum[m] += D_filt[j*n_samples[i] + m];
						else
							dp_sum[m] += D[j*n_samples[i] +m];
						dp_cnt[m] ++;
					}
				}
				for(int m=0; m<n_samples[i]; ++m)
					dvec_dp[k+3][m+sample_idx] = (double)dp_sum[m]/dp_cnt[m]/32.0 ;
			}
            
			sample_idx += n_samples[i];
			
			delete [] D;
			delete [] D_filt;
		}
	} // if n_dp <= 200
	else // n_dp >  200
	{
		for(int i=0; i<n_pileup; ++i)
		{
			// if n_inside < 200, make two
			// sample by sample sum
			std::vector<int> pos_starts;
            pos_starts.push_back(curr_sv.pos + (int)curr_sv.len / 8);
            pos_starts.push_back(curr_sv.pos + (int)curr_sv.len * 0.625);
            /*
			pos_starts.push_back(curr_sv.pos + 101);
			pos_starts.push_back(curr_sv.pos + (int)curr_sv.len / 4);
			pos_starts.push_back(curr_sv.pos + (int)curr_sv.len / 2);
			pos_starts.push_back(curr_sv.pos + (int)curr_sv.len*0.75);
*/
			std::vector<unsigned> gc_sum (n_samples[i], 0);
			n_dp = 50;
			int n_dp_by_sample = n_samples[i] * n_dp;
			uint16_t *D = new uint16_t[n_dp_by_sample];
			uint16_t *D_filt = new uint16_t[n_dp_by_sample];

			for (int seg=0; seg<2; ++seg)
			{
				uint64_t start_byte = multi_idx[i][0]; // This is the index position where dp100 record starts
				start_byte = chr_bytepos_dp100[i][curr_sv.chrnum];

				idx_start = (int)pos_starts[seg] / 100;
				start_byte += idx_start * n_samples[i] * sizeof(uint16_t);
				
				pileups[i].seekg(start_byte);
				pileups[i].read_depth(D, n_dp_by_sample);

				for(int j=0; j<n_dp; ++j)
				{
					for(int k=0; k<n_samples[i]; ++k)
					{
						D[j*n_samples[i] + k] = (uint16_t) correct_gc(gc, sample_idx+k, (double)D[j*n_samples[i] + k], curr_sv.chrnum, (pos_starts[seg] + 50 + j*100) );
						gc_sum[k] += D[j*n_samples[i] + k];
					}
				}

				median_filter(D, D_filt, n_samples[i], n_dp);

				std::vector<unsigned> dp_sum (n_samples[i], 0);
				
				for(int j=0; j<n_dp ; ++j)
				{
					for(int m=0; m<n_samples[i]; ++m)
					{
						dp_sum[m] += D_filt[j*n_samples[i] + m];
					}
				}

				for(int m=0; m<n_samples[i]; ++m)
					dvec_dp[seg+3][m+sample_idx] = (double)dp_sum[m]/n_dp/32.0;

			} // seg : 2

			for(int k=0; k<n_samples[i]; ++k)
			{
				dvec_dp[2][k + sample_idx] = (double) gc_sum[k] / (n_dp * 2)/ 32.0;
			}
			sample_idx += n_samples[i];

			delete [] D;
			delete [] D_filt;
			
		} // i : n_pileup
	}

    return 1;
}

void DataReader::read_pair_split(sv& curr_sv, std::vector<ReadStat>& rdstats, GcContent &gc)
{
    int sample_idx = 0;

	if (curr_sv.len < 50) return;
    // Starting breakpoint
    int start_pos1 = curr_sv.pos - 500;
    int end_pos1 = curr_sv.pos + 500;

    if (start_pos1<1) start_pos1 = 1;
	if (end_pos1 >= curr_sv.end) end_pos1 = curr_sv.end-1;

    int start_pos2 = curr_sv.end - 500;
    int end_pos2 = curr_sv.end + 500;

    if (start_pos2 <= curr_sv.pos) start_pos2 = curr_sv.pos+1;
	if (end_pos2 > (int)gc.chr_size[curr_sv.chrnum])  end_pos2 = gc.chr_size[curr_sv.chrnum]-1;
    
    uint64_t start_idx1 = chr_idx_rp[curr_sv.chrnum] + (int)start_pos1/10000;
    uint64_t end_idx1 = chr_idx_rp[curr_sv.chrnum] + (int)end_pos1/10000;

    uint64_t start_idx2 = chr_idx_rp[curr_sv.chrnum] + (int)start_pos2/10000;
    uint64_t end_idx2 = chr_idx_rp[curr_sv.chrnum] + (int)end_pos2/10000;

	bool b_overlap = (end_idx1 >= start_idx2); // in fact, > case should never happen

	uint64_t start_idx = start_idx1;
	uint64_t end_idx = (b_overlap) ? end_idx2 : end_idx1;

    for(int i=0; i<n_pileup; ++i)
    {
        pileups[i].seekg(multi_idx[i][start_idx]);

        for(uint64_t j=start_idx; j<=end_idx; ++j)
        {
            
            for(int k=0; k<n_samples[i]; ++k)
            {
//				printf("\nSample %d\n", sample_idx+k);
                uint32_t n_rp = 0;
                pileups[i].read_uint32(n_rp);
 
                for(uint32_t l=0; l<n_rp; ++l)
                {
                    readpair rp;
                    pileups[i].read_readpair(rp);
                    rp.chrnum = curr_sv.chrnum;
                    rp.selfpos += (j - chr_idx_rp[curr_sv.chrnum]) * 10000;
                    rp.matepos += (j - chr_idx_rp[curr_sv.chrnum]) * 10000;

					if (rp.selfpos >= start_pos1 && rp.selfpos <= end_pos1)
					{
						if (rp.matequal>0)
						{
							if (rp.matepos >= start_pos2 && rp.matepos <= end_pos2)
							{
                                //printf(">>>");
								if (rp.pairstr == 1 && rp.selfpos < rp.matepos)
									rdstats[sample_idx + k].n_pre_FR ++;
	 							else if (rp.pairstr == 2 && rp.selfpos <= rp.matepos)
									rdstats[sample_idx + k].n_pre_RF ++;
                                else if (rp.pairstr == 0 && rp.selfpos <= rp.matepos)
                                    rdstats[sample_idx + k].n_pre_FF ++;
                                
                            }
                            else if (rp.matepos >= start_pos1 && rp.matepos <= rp.selfpos + 300) // insert size 450 or below
                            {
                                if ((rp.pairstr == 1 && rp.selfpos < rp.matepos))
                                {
                                    rdstats[sample_idx + k].n_pre_INS ++;
                                }
                            }
						}
						else
                            rdstats[sample_idx + k].n_pre_rp_missing ++;
//                        printf("PRE_RP\t%d\t%d\t%d\t%u\t%d\n", rp.chrnum, rp.selfpos, rp.matepos, rp.matequal, rp.pairstr);

					}
 					if (b_overlap && rp.selfpos >= start_pos2 && rp.selfpos <= end_pos2)
					{
						if (rp.matequal > 0)
						{
                            if (rp.matepos >= start_pos1 && rp.matepos <= end_pos1)
                            {
                                //printf(">>>");

                                if (rp.pairstr == 2 && rp.selfpos > rp.matepos)
                                    rdstats[sample_idx + k].n_post_FR ++;
                                else if (rp.pairstr == 1 && rp.matepos <= rp.selfpos)
                                    rdstats[sample_idx + k].n_post_RF ++;
                                else if (rp.pairstr == 3 && rp.selfpos <= rp.matepos)
                                    rdstats[sample_idx + k].n_post_RR ++;
                            }
						}
						else
                            rdstats[sample_idx+k].n_post_rp_missing ++;
             //           printf("POST_RP\t%d\t%d\t%d\t%u\t%d\n", rp.chrnum, rp.selfpos, rp.matepos, rp.matequal, rp.pairstr);

					}
                }
				//printf("\n");
                uint32_t n_sp = 0;
                pileups[i].read_uint32(n_sp);
                for(uint32_t l=0; l<n_sp; ++l)
                {
                    splitread sp;
                    pileups[i].read_splitread(sp);
                    sp.chrnum = curr_sv.chrnum;
                    sp.pos += (j - chr_idx_rp[curr_sv.chrnum]) * 10000;
                    sp.sapos += (j - chr_idx_rp[curr_sv.chrnum]) * 10000;
                    
					if (sp.pos >=start_pos1 && sp.pos <= end_pos1)
                    {
                        if (sp.sapos >= start_pos2 && sp.sapos <= end_pos2)
                        {
                            //printf(">>>");

                            if (sp.pos < sp.sapos)
                            {
                                if (sp.firstclip <= 0 && sp.secondclip >=0)
                                {
                                    rdstats[sample_idx+k].n_pre_split_in ++;
                                }
                                else if (sp.firstclip >=0 && sp.secondclip <=0)
                                {
                                    rdstats[sample_idx+k].n_pre_split_out ++;
                                }
                            }

                        }
                        else if (sp.firstclip < 0 )
                        {
                            rdstats[sample_idx+k].n_pre_clip_in ++;
                        }
                        else if (sp.firstclip > 0 )
                        {
                            rdstats[sample_idx+k].n_pre_clip_out ++;
                        }
                        else if (sp.sapos == 0)
                        {
                            rdstats[sample_idx+k].n_pre_sp_missing ++;
                        }
                        //printf("PRE_SP\t%d\t%d\t%d\t%d\t%d\n", sp.chrnum, sp.pos, sp.sapos, sp.firstclip, sp.secondclip);

                    }
                    if (b_overlap && sp.pos >=start_pos2 && sp.pos <= end_pos2)
                    {
                        if (sp.sapos >= start_pos1 && sp.sapos <= end_pos1)
                        {
                            //printf(">>>");

                            if (sp.pos > sp.sapos)
                            {
                                if (sp.firstclip <= 0 && sp.secondclip >=0)
                                {
                                    rdstats[sample_idx+k].n_pre_split_out ++;
                                }
                                else if (sp.firstclip >=0 && sp.secondclip <=0)
                                {
                                    rdstats[sample_idx+k].n_pre_split_in ++;
                                }
                            }
                        }
                        else if (sp.firstclip > 0 )
                        {
                            rdstats[sample_idx+k].n_post_clip_in ++;
                        }
                        else if (sp.firstclip < 0 )
                        {
                            rdstats[sample_idx+k].n_post_clip_out ++;
                        }
                        else if (sp.sapos == 0)
                        {
                            rdstats[sample_idx+k].n_post_sp_missing ++;
                        }
                        //printf("POST_SP\t%d\t%d\t%d\t%d\t%d\n", sp.chrnum, sp.pos, sp.sapos, sp.firstclip, sp.secondclip);

                    }
                }
                
                uint32_t n_lclip = 0;
                pileups[i].read_uint32(n_lclip);
                for(uint32_t l=0; l<n_lclip; ++l)
                {
                    sclip myclip;
                    pileups[i].read_softclip(myclip);
                    myclip.chrnum = curr_sv.chrnum;
                    myclip.pos += (j - chr_idx_rp[curr_sv.chrnum]) * 10000;
                }
                
                uint32_t n_rclip = 0;
                pileups[i].read_uint32(n_rclip);
                for(uint32_t l=0; l<n_rclip; ++l)
                {
                    sclip myclip;
                    pileups[i].read_softclip(myclip);
                    myclip.chrnum = curr_sv.chrnum;
                    myclip.pos += (j - chr_idx_rp[curr_sv.chrnum]) * 10000;
                }
            }
        }

        sample_idx += n_samples[i];
    }
    
    if (b_overlap) return;
    
    
    // TODO: REFACTOR THIS
    // Non-overlapping case
    start_idx = start_idx2;
    end_idx = end_idx2;
    sample_idx = 0;
    
    for(int i=0; i<n_pileup; ++i)
    {
        pileups[i].seekg(multi_idx[i][start_idx]);
        
        for(uint64_t j=start_idx; j<=end_idx; ++j)
        {
            
            for(int k=0; k<n_samples[i]; ++k)
            {
                //printf("\nSample %d\n", sample_idx+k);
                uint32_t n_rp = 0;
                pileups[i].read_uint32(n_rp);
                
                for(uint32_t l=0; l<n_rp; ++l)
                {
                    readpair rp;
                    pileups[i].read_readpair(rp);
                    rp.chrnum = curr_sv.chrnum;
                    rp.selfpos += (j - chr_idx_rp[curr_sv.chrnum]) * 10000;
                    rp.matepos += (j - chr_idx_rp[curr_sv.chrnum]) * 10000;
                    
                    // THIS CANNOT HAPPEN BECAUSE ONLY + ISIZE READ PAIRS ARE WRITTEN
                    if (rp.selfpos >= start_pos2 && rp.selfpos <= end_pos2)
                    {
                        //printf(">>>");
                        if (rp.matequal > 0)
                        {
                            if (rp.matepos >= start_pos1 && rp.matepos <= end_pos1)
                            {
                                if (rp.pairstr == 2 && rp.selfpos > rp.matepos)
                                    rdstats[sample_idx + k].n_post_FR ++;
                                else if (rp.pairstr == 1 && rp.matepos <= rp.selfpos)
                                    rdstats[sample_idx + k].n_post_RF ++;
                                else if (rp.pairstr == 3 && rp.selfpos <= rp.matepos)
                                    rdstats[sample_idx + k].n_post_RR ++;
                            }
                        }
                        else
                            rdstats[sample_idx+k].n_post_rp_missing ++;
                        //printf("POST_RP\t%d\t%d\t%d\t%u\t%d\n", rp.chrnum, rp.selfpos, rp.matepos, rp.matequal, rp.pairstr);

                    }
                }
                //printf("\n");
                uint32_t n_sp = 0;
                pileups[i].read_uint32(n_sp);
                for(uint32_t l=0; l<n_sp; ++l)
                {
                    splitread sp;
                    pileups[i].read_splitread(sp);
                    sp.chrnum = curr_sv.chrnum;
                    sp.pos += (j - chr_idx_rp[curr_sv.chrnum]) * 10000;
                    sp.sapos += (j - chr_idx_rp[curr_sv.chrnum]) * 10000;
                    
      
                    if (sp.pos >=start_pos2 && sp.pos <= end_pos2)
                    {
                        if (sp.sapos >= start_pos1 && sp.sapos <= end_pos1)
                        {
                            //printf(">>>");
                            if (sp.pos > sp.sapos)
                            {
                                if (sp.firstclip <= 0 && sp.secondclip >=0)
                                {
                                    rdstats[sample_idx+k].n_post_split_out ++;
                                }
                                else if (sp.firstclip >=0 && sp.secondclip <=0)
                                {
                                    rdstats[sample_idx+k].n_pre_split_in ++;
                                }
                            }
                        }
                        else if (sp.sapos == 0)
                        {
                            rdstats[sample_idx+k].n_post_sp_missing ++;
                        }
                        //printf("POST_SP\t%d\t%d\t%d\t%d\t%d\n", sp.chrnum, sp.pos, sp.sapos, sp.firstclip, sp.secondclip);

                    }
                }
                
                uint32_t n_lclip = 0;
                pileups[i].read_uint32(n_lclip);
                for(uint32_t l=0; l<n_lclip; ++l)
                {
                    sclip myclip;
                    pileups[i].read_softclip(myclip);
                    myclip.chrnum = curr_sv.chrnum;
                    myclip.pos += (j - chr_idx_rp[curr_sv.chrnum]) * 10000;
                }
                
                uint32_t n_rclip = 0;
                pileups[i].read_uint32(n_rclip);
                for(uint32_t l=0; l<n_rclip; ++l)
                {
                    sclip myclip;
                    pileups[i].read_softclip(myclip);
                    myclip.chrnum = curr_sv.chrnum;
                    myclip.pos += (j - chr_idx_rp[curr_sv.chrnum]) * 10000;
                }
                
            }
        }
        
        sample_idx += n_samples[i];
    }
    
    
    return;
}

void DataReader::read_var_depth(int var_i, std::vector<double> &var_dp)
{
	// var_i : i-th var
	uint64_t bytepos = 0;
	int sample_idx = 0;
	uint16_t D[1000];  // max number of sample per each pileup, TEMPORARY

	for(int i=0; i<n_pileup; ++i)
	{
        // var_i : zero-index
        bytepos = 2*sizeof(int32_t) + 256*n_samples[i] + (var_i * sizeof(uint16_t) * n_samples[i]);
		var_files[i].seekg(bytepos);
		var_files[i].read_depth(D, n_samples[i]);
		for(int j=0; j<n_samples[i]; ++j)
		{
            var_dp[sample_idx + j] = (double)D[j]/32.0;
		}
        sample_idx += n_samples[i];
	}
    return;
}

// GC-correction of n-th sample at chr:pos
double DataReader::correct_gc(GcContent& gc, int n, double depth, int chr, int pos)
{
    int p = pos / gc.interval_dist;
    //    std::cerr << "pos " << pos << " p " << p << std::endl;
    int bin = gc.gc_array[chr][p];
    
    if (bin<(gc.num_bin-1) && gc_factors[n][bin]>0.1)
    {
        //        std::cerr << D << " at " << chr << ":" << pos << " is adjusted to " << D/gc_factor[bin] << " by gc Factor "<< gc_factor[bin] << std::endl;
        return depth * gc_factors[n][bin];
    }
    else
    {
        // Let's not make adjustment for bins with '255' value
        return depth;
    }
}

void DataReader::close()
{
    for(int i=0; i<n_pileup; ++i)
    {
        pileups[i].close();
        var_files[i].close();
        idx_files[i].close();
    }
}
