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
    int rp_cnt = n_rp[1]; // FR
    int clip_cnt = n_split_inward;
    
    if (rp_cnt > 4 || rp_cnt + clip_cnt > 7)
        return true;
    else
        return false;
}

bool ReadStat::dup_support()
{
    int rp_cnt = n_rp[2]; // RF
    int sp_cnt = n_split_outward;
    
    if (rp_cnt > 4  || rp_cnt + sp_cnt > 7)
        return true;
    else
        return false;
}

bool ReadStat::inv_support()
{
    int rp_cnt = n_rp[0] + n_rp[3];
    
    // inversion can have both in and out clipping
    int sp_cnt = n_split_inward + n_split_outward;
    
    if (rp_cnt > 4 && rp_cnt + sp_cnt> 7)
        return true;
    else
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
    // chr_idx_rp[] contains the array indices of the idxfile's indices
    // which points to the beginning of read pair blocks, for each chromosome
    
	if (chr>0) // if pileup files are split by chromosomes
	{
		for(int i=1;i<=gc.num_chr; ++i)
		{
			chr_idx_rp[i] = -1; // make everything invalid 
		}
		chr_idx_rp[chr] = idx_cnt;
        // if chr==20, chr_idx_rp[20] == 1;, chr_idx_rp[1..19], chr_idx_rp[21..24] == -1

		idx_cnt += ceil((double)gc.chr_size[chr] / 10000.0) ;
        // idx_cnt increases by number of 10kbp blocks in the single chromosome
	}
	else // if pileup files contains all chromosomes
	{
		for(int i=1;i<=gc.num_chr; ++i)
		{
			chr_idx_rp[i] = idx_cnt;
			idx_cnt += ceil((double)gc.chr_size[i] / 10000.0);
		}
        // idx_cnt increase by number of 10kbp blocks in all chromosomes, added up
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

void DataReader::adjust_gc_factor(GcContent& gc, std::vector<SampleStat>& stats, int chrnum)
{
    /*
    int sample_idx = 0;
    for(int i=0; i<n_pileup; ++i)
    {
        uint64_t start_byte = chr_bytepos_dp100[i][chrnum];
        
        int N = gc.n_interval[chrnum];
        
        int n_interval_by_sample = N * n_samples[i];
        
        uint16_t *D = new uint16_t[n_interval_by_sample];
        
        pileups[i].seekg(start_byte);
        pileups[i].read_depth(D, n_interval_by_sample);
        
        int64_t gc_sum[n_samples[i]][gc.num_bin];
        int64_t gc_cnt[n_samples[i]][gc.num_bin];
        double gc_avg[n_samples[i]][gc.num_bin];

        for(int j=0; j<n_samples[i]; ++j)
        {
            for(int k=0; k<gc.num_bin; ++k)
            {
                gc_sum[j][k] = 0;
                gc_cnt[j][k] = 0;
            }
        }
        for(int k=0; k<N; ++k)
        {
            double gcval = gc.get_gc_content(chrnum, k* gc.interval_dist, (k+1)*gc.interval_dist);
            for(int j=0; j<n_samples[i]; ++j)
            {
                int idx = (int) round(gcval*100);
                gc_sum[j][ idx ] += D[k*n_samples[i] + j]/32;
                gc_cnt[j][ idx ] ++;
            }
        }
        
        for(int j=0; j<n_samples[i]; ++j)
        {
            for(int k=2; k<gc.num_bin-2; ++k)
            {
                if (gc_cnt[j][k] + (gc_cnt[j][k-1] + gc_cnt[j][k+1])/2.0 + (gc_cnt[j][k-2]+gc_cnt[j][k+2])/4.0 > 100)
                {
                    gc_avg[j][k] = ((double)gc_sum[j][k] + 0.5*((double)gc_sum[j][k-1] + gc_sum[j][k+1]) + 0.25 * ((double)gc_sum[j][k-2] + gc_sum[j][k+2])) / (gc_cnt[j][k] +0.5*(gc_cnt[j][k-1]+gc_cnt[j][k+1]) + 0.25*(gc_cnt[j][k-2] + gc_cnt[j][k+2])) ;
                }
                else
                {
                    gc_avg[j][k] = 0;
                }
            }
            gc_avg[j][0] = gc_avg[j][1] = gc_avg[j][2];
            gc_avg[j][gc.num_bin-1] = gc_avg[j][gc.num_bin-2] = gc_avg[j][gc.num_bin-3];
        }
        
        for(int k=0;k<n_samples[i]; ++k)
        {
            for(int j=0; j<gc.num_bin; ++j)
            {
                fprintf(stderr, "sample %d, gc bin %d, avg dp %f\n", sample_idx + k, j, gc_avg[k][j]);
            }
        }
        for(int j=0; j<n_samples[i]; ++j)
        {
            fprintf(stderr, "sample %d, stat avg dp %f\n", sample_idx + j, stats[sample_idx+j].avg_dp);

            for(int k=0; k<gc.num_bin; ++k)
            {
                if (gc_avg[j][k]>0)
                {
                    gc_factors[sample_idx + j ][k] = stats[sample_idx+j].avg_dp / gc_avg[j][k];
                }
                else
                {
                    gc_factors[sample_idx + j][k] = 0;
                }
            }
        }
        sample_idx += n_samples[i];
    } */
    for(int i=0; i<n_sample_total; ++i)
    {
        for(int j=0; j<gc.num_bin; ++j)
        {
            fprintf(stderr, "sample %d, gc bin %d, gc factor %f\n", i, j, gc_factors[i][j]);
        }
    }
}

bool DataReader::read_depth100(sv& curr_sv, std::vector< std::vector<double> > &dvec_dp, GcContent& gc)
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
        return false;
    }
    else if (n_dp <= 10)
	{
        b_medfilt = false; // no median filtering
	}
    
	if (curr_sv.svtype == INV)
	{
		return false;
	}
    // deprecated 4/13/19
   // for(int k=0; k<2; ++k)
    //{
     //   std::vector<double> x(n_sample_total, 0);
      //  dvec_dp.push_back(x);
 //   }

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

                if (b_medfilt)
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

    return true;
}

bool DataReader::around_breakpoint(readpair &rp, sv &curr_sv)
{
    // It's already forced that matepos > selfpos
    
	switch(curr_sv.svtype)
	{
	case DEL:
		if (rp.selfpos >= curr_sv.pos - 500 && rp.selfpos <= curr_sv.pos + 100 && rp.matepos < curr_sv.end + 500 && rp.matepos >= curr_sv.end - 100 && rp.matequal > 0)
		{
			return true;
		}
		else
		{
			return false;
		}
		break;
	case DUP:
	case CNV:
		if (rp.selfpos >= curr_sv.pos - 100 && rp.selfpos < curr_sv.pos+500 && rp.matepos <= curr_sv.end + 100 && rp.matepos >= curr_sv.end -500 && rp.matequal > 0)
		{
			return true;
		}
		else
		{
			return false;
		}

		break;
	case INV:
		if ( (rp.selfpos >= curr_sv.pos - 500 && rp.selfpos < curr_sv.pos && rp.matepos < curr_sv.end  && rp.matepos >= curr_sv.end - 500 && rp.matepos > curr_sv.pos && rp.matequal > 30 && rp.pairstr == 0) || (rp.selfpos >= curr_sv.pos - 50 && rp.selfpos < curr_sv.pos + 500 && rp.selfpos < curr_sv.end && rp.matepos < curr_sv.end + 500 && rp.matepos >= curr_sv.end - 50 && rp.matequal > 30 && rp.pairstr == 3))
		{
			return true;
		}
		else
		{
			return false;
		}
		break;
    case BND:
    case INS:
        break;
	}
    return false;
}


void DataReader::read_pair_split(sv& curr_sv, std::vector<ReadStat>& rdstats, GcContent &gc, std::vector< std::vector<int> > &all_rps, std::vector<int> &all_lclips, std::vector<int> &all_rclips)
{
    int offset = 0;

//	if (curr_sv.len < 50) return; // clipped read works for short variants!
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

	bool b_overlap = (end_idx1 >= start_idx2);

	uint64_t start_idx = start_idx1;
	uint64_t end_idx = (b_overlap) ? end_idx2 : end_idx1;

    for(int repeat=0; repeat<2; ++repeat)
    {
        for(int batch_idx=0; batch_idx<n_pileup; ++batch_idx)
        {
            pileups[batch_idx].seekg(multi_idx[batch_idx][start_idx]);

            for(uint64_t j=start_idx; j<=end_idx; ++j)
            {
                for(int sample_idx=0; sample_idx<n_samples[batch_idx]; ++sample_idx)
                {
                    uint32_t n_rp = 0;
                    pileups[batch_idx].read_uint32(n_rp);
     
                    for(uint32_t l=0; l<n_rp; ++l)
                    {
                        readpair rp;
                        pileups[batch_idx].read_readpair(rp);
                        
                        rp.chrnum = curr_sv.chrnum;
                        
                        if (rp.matepos < rp.selfpos)
                        {
                            // swap them, always selfpos <= matepos
                            int pos = rp.selfpos;
                            rp.selfpos = rp.matepos;
                            rp.matepos = pos;
                            if (rp.pairstr == 1)
                            {
                                rp.pairstr = 2;
                            }
                            else if (rp.pairstr == 2)
                            {
                                rp.pairstr = 1;
                            }
                        }

                        if (around_breakpoint(rp, curr_sv))
                        {
                            int self_idx = (rp.selfpos - curr_sv.pos + 500 ) /10;
                            int mate_idx = (rp.matepos - curr_sv.end + 500 ) /10 + 100;
                            
                            if (self_idx <0 || self_idx>99 || mate_idx<100 || mate_idx>199)
                            {
                                // sanity check
                                std::cerr << "Something wrong with read pair indices selfidx " << self_idx << " mate idx " << mate_idx << std::endl;
                                exit(1);
                                
                            }
							rdstats[offset + sample_idx].n_rp[rp.pairstr] ++;
							
							int buf_start = self_idx-5;
							int buf_end = self_idx+5;
							if (buf_start<0) buf_start = 0;
							if (buf_end>99) buf_end = 99;
							
							for(int m=buf_start; m<=buf_end;++m)
							{
								// frequency in 100-bp window
								rdstats[offset + sample_idx].rp_seq[rp.pairstr][m]++;
								all_rps[rp.pairstr][m] ++;
							}
							buf_start = mate_idx - 5;
							buf_end = mate_idx + 5;
							if (buf_start<100) buf_start = 100;
							if (buf_end>199) buf_end = 199;
							for(int m=buf_start;m<=buf_end; ++m)
							{
								rdstats[offset + sample_idx].rp_seq[rp.pairstr][m]++;
								all_rps[rp.pairstr][m] ++;
							}
                            
                            PairSplit new_pair;
                            new_pair.positions.first = rp.selfpos;
                            new_pair.positions.second = rp.matepos;
                            // directions.first: indicate self forward (true) or self reverse (false)
                            // directions.second: indicate mate forward (true) or mate reverse (false)
                            // for deletions, should always be FR
                            // for duplications, RF
                            // inversions RR or FF
                            new_pair.directions.first = (rp.pairstr < 2); // pairstr 00b = FF, 01b = FR, 10b = RF, 11b = RR
                            new_pair.directions.second = (rp.pairstr == 0 || rp.pairstr == 2);
                            rdstats[offset+sample_idx].readpairs.push_back(new_pair);
                        }
                    }
                    //printf("\n");
                    uint32_t n_sp = 0;
                    pileups[batch_idx].read_uint32(n_sp);
                    for(uint32_t l=0; l<n_sp; ++l)
                    {
                        splitread sp;
                        pileups[batch_idx].read_splitread(sp);
                        
                        sp.chrnum = curr_sv.chrnum;
                        
                        // ||| : clipped bases
                        //
                        // 1. ----||||               ||||---- : DEL  , pos < sapos && firstclip < 0 && secondclip > 0, or all <> opposite
                        // 2. ||||----               ----|||| : DUP  , pos < sapos && firstclip > 0 && secondclip < 0, or all <> opposite
                        
                        // 3. ----||||               ----|||| : INV when strands are opposite, currently unavailable
                        // 4. ||||----               ||||---- : INV when strands are opposite, currently unavailable

                        /*
                        if (sp.pos > sp.sapos)
                        {
                            int pos = sp.pos;
                            sp.pos = sp.sapos;
                            sp.sapos = pos;
                            int clip = sp.firstclip;
                            sp.firstclip = sp.secondclip;
                            sp.secondclip = clip;
                            
                        }*/
                        
                        int break1 = sp.pos;
                        int break2 = sp.sapos;
                        
                        if (sp.firstclip < 0)
                        {
                            // TODO: 150 is hard coded
                            break1 += 151 + sp.firstclip;
                        }
                        if (sp.secondclip < 0)
                        {
                            break2 += 151 + sp.secondclip;
                            
                        }

                        // TODO: make split read arrays (+/- 100bp around sv start-end) a class
                        
                        if (break1 < break2)
                        {
                            if (sp.firstclip < 0 && sp.secondclip > 0) // DEL, inward split (pos)----||||      ||||(sapos)----
                            {
                                break2 = break2 - 1;

                                if (break1 >= curr_sv.pos - 100 && break1 < curr_sv.pos + 100 && break2 >= curr_sv.end-100 && break2 < curr_sv.end + 100)
                                {
                                    rdstats[offset+sample_idx].n_split_inward ++;
                                    rdstats[offset+sample_idx].sp_seq_in[break1 - curr_sv.pos + 100] ++;
                                    rdstats[offset+sample_idx].sp_seq_in[break2 - curr_sv.end + 300] ++;
                                    
                                    // std::cout << "break1 < break2, DEL, break1 " << break1 << " break2 " << break2 << std::endl;

                                    PairSplit new_split;
                                    new_split.positions.first = break1;
                                    new_split.positions.second = break2;
                                    new_split.directions.first = true;
                                    new_split.directions.second = false;
                                    rdstats[offset+sample_idx].splits.push_back(new_split);
                                }
                            }
                            else if (sp.firstclip > 0 && sp.secondclip < 0) // DUP, outward split ||||(pos)----    ----(sapos)||||
                            {
                                break2 = break2 - 1;
                                if (break1 >= curr_sv.pos - 100 && break1 < curr_sv.pos + 100 && break2 >= curr_sv.end-100 && break2 < curr_sv.end + 100)
                                {
                                    rdstats[offset+sample_idx].n_split_outward ++;
                                    rdstats[offset+sample_idx].sp_seq_out[break1 - curr_sv.pos + 100] ++;
                                    rdstats[offset+sample_idx].sp_seq_out[break2 - curr_sv.end + 300] ++;
                                    
                                    // std::cout << "break1 < break2, DUP, break1 " << break1 << " break2 " << break2 << std::endl;

                                    
                                    PairSplit new_split;
                                    new_split.positions.first = break1;
                                    new_split.positions.second = break2;
                                    new_split.directions.first = false;
                                    new_split.directions.second = true;
                                    rdstats[offset+sample_idx].splits.push_back(new_split);
                                }
                            }
                        }
                        else // break1 >= break2
                        {
                            if (sp.firstclip > 0 && sp.secondclip < 0) // DEL, inward split (sapos)----||||      ||||(pos)----
                            {
                                break2 = break2 -1 ;
                                if (break2 >= curr_sv.pos - 100 && break2 < curr_sv.pos + 100 && break1 >= curr_sv.end-100 && break1 < curr_sv.end + 100)
                                {
                                    rdstats[offset+sample_idx].n_split_inward ++;
                                    rdstats[offset+sample_idx].sp_seq_in[break2 - curr_sv.pos + 100] ++;
                                    rdstats[offset+sample_idx].sp_seq_in[break1 - curr_sv.end + 300] ++;
                                    
                                   // std::cout << "break2 < break1, DEL, break1 " << break1 << " break2 " << break2 << std::endl;
                                    
                                    PairSplit new_split;
                                    new_split.positions.first = break2;
                                    new_split.positions.second = break1;
                                    new_split.directions.first = true;
                                    new_split.directions.second = false;
                                    rdstats[offset+sample_idx].splits.push_back(new_split);
                                }
                            }
                            else if (sp.firstclip < 0 && sp.secondclip > 0) // DUP, outward split ||||(sapos)----    (pos)----||||
                            {
                                break2 = break2 - 1;
                                if (break2 >= curr_sv.pos - 100 && break2 < curr_sv.pos + 100 && break1 >= curr_sv.end-100 && break1 < curr_sv.end + 100)
                                {
                                    rdstats[offset+sample_idx].n_split_outward ++;
                                    rdstats[offset+sample_idx].sp_seq_out[break2 - curr_sv.pos + 100] ++;
                                    rdstats[offset+sample_idx].sp_seq_out[break1 - curr_sv.end + 300] ++;
                                    
                                    // std::cout << "break2 < break1, DUP, break1 " << break1 << " break2 " << break2 << std::endl;

                                    PairSplit new_split;
                                    new_split.positions.first = break2;
                                    new_split.positions.second = break1;
                                    new_split.directions.first = false;
                                    new_split.directions.second = true;
                                    rdstats[offset+sample_idx].splits.push_back(new_split);
                                }
                            }
                            
                        }
                    }
                    
                    uint32_t n_lclip = 0;
                    pileups[batch_idx].read_uint32(n_lclip);
                    for(uint32_t l=0; l<n_lclip; ++l)
                    {
                        sclip myclip;
                        pileups[batch_idx].read_softclip(myclip);
                        myclip.chrnum = curr_sv.chrnum;
                 
                        if (myclip.pos >= curr_sv.pos - 100 && myclip.pos < curr_sv.pos + 100)
                        {
                            int tmp_idx = myclip.pos - curr_sv.pos + 100;
                            
                            if ( tmp_idx <0 || tmp_idx > 199)
                            {
                                // sanity check
                                std::cerr << "Something wrong with clip idx  " << tmp_idx << std::endl;
                                exit(1);
                            }
                            
                            rdstats[offset+sample_idx].lclips[tmp_idx] ++;
                            rdstats[offset+sample_idx].n_lclip_start ++;
                            all_lclips[tmp_idx] ++;
                        }
                        if (myclip.pos >= curr_sv.end - 100 && myclip.pos < curr_sv.end + 100 )
                        {
                            int tmp_idx = myclip.pos - curr_sv.end + 300;
                            
                            if ( tmp_idx < 200 || tmp_idx > 399)
                            {
                                // sanity check
                                std::cerr << "Something wrong with clip idx end " << tmp_idx << std::endl;
                                exit(1);
                            }
                            
                            rdstats[offset+sample_idx].lclips[tmp_idx] ++;
                            rdstats[offset+sample_idx].n_lclip_end ++;
                            all_lclips[tmp_idx] ++;
                        }
                    }
                    
                    uint32_t n_rclip = 0;
                    pileups[batch_idx].read_uint32(n_rclip);
                    for(uint32_t l=0; l<n_rclip; ++l)
                    {
                        sclip myclip;
                        pileups[batch_idx].read_softclip(myclip);
                        myclip.chrnum = curr_sv.chrnum;

                        if (myclip.pos >= curr_sv.pos - 100 && myclip.pos < curr_sv.pos + 100)
                        {
                            int tmp_idx = myclip.pos - curr_sv.pos + 100;
                            
                            if ( tmp_idx <0 || tmp_idx > 199)
                            {
                                // sanity check
                                std::cerr << "Something wrong with clip idx  " << tmp_idx << std::endl;
                                exit(1);
                            }
                            rdstats[offset+sample_idx].rclips[tmp_idx] ++;
                            rdstats[offset+sample_idx].n_rclip_start ++;
                            all_rclips[tmp_idx] ++;
                        }
                        if (myclip.pos >= curr_sv.end - 100 && myclip.pos < curr_sv.end + 100 )
                        {
                            int tmp_idx = myclip.pos - curr_sv.end + 300;
                            
                            if ( tmp_idx < 200 || tmp_idx > 399)
                            {
                                // sanity check
                                std::cerr << "Something wrong with clip idx end " << tmp_idx << std::endl;
                                exit(1);
                            }
                            rdstats[offset+sample_idx].rclips[tmp_idx] ++;
                            rdstats[offset+sample_idx].n_rclip_end ++;
                            all_rclips[tmp_idx] ++;
                        }
                    }
                }
            }

            offset += n_samples[batch_idx];
        }
        
        if (b_overlap)
            return;
        
        // Prepare to read second readpair interval
        start_idx = start_idx2;
        end_idx = end_idx2;
        offset = 0;
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
    // TEMPORARY!!!
//    return depth;
    
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
