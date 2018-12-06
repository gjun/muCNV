//
//  multi_pileup.cpp
//  muCNV
//
//  Created by Goo Jun on 11/25/18.
//  Copyright © 2018 Goo Jun. All rights reserved.
//

#include "data_reader.h"
#include <math.h>

void medianfilter(const element* signal, element* result, int N)
{
    //   Move window through all elements of the signal
    for (int i = 2; i < N - 2; ++i)
    {
        //   Pick up window elements
        element window[5];
        for (int j = 0; j < 5; ++j)
            window[j] = signal[i - 2 + j];
        //   Order elements (only half of them)
        for (int j = 0; j < 3; ++j)
        {
            //   Find position of minimum element
            int min = j;
            for (int k = j + 1; k < 5; ++k)
                if (window[k] < window[min])
                    min = k;
            //   Put found minimum element in its place
            const element temp = window[j];
            window[j] = window[min];
            window[min] = temp;
        }
        //   Get result - the middle element
        result[i - 2] = window[2];
    }
}

int DataReader::load(std::vector<string>& base_names, std::vector<SampleStat> &stats, GcContent &gc)
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
    for(int i=1;i<=gc.num_chr; ++i)
    {
        chr_idx_rp[i] = idx_cnt; // chr_idx[1] contains the array index of pileup index of chr 1
        idx_cnt += ceil((double)gc.chr_size[i] / 10000.0) ;
    }
    
    std::cerr<< idx_cnt << " indices should be in index file" << std::endl;
 
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
        
        // idx offset?
    }

    chr_bytepos_dp100.resize(n_pileup);
    for(int i=0; i<n_pileup; ++i)
    {
        chr_bytepos_dp100[i].resize(gc.num_chr +1);
        chr_bytepos_dp100[i][0]= 0;
        chr_bytepos_dp100[i][1] = multi_idx[0][0];
        for(int c=2; c<=gc.num_chr; ++c)
        {
            chr_bytepos_dp100[i][c] = chr_bytepos_dp100[i][c-1] + (ceil((double)gc.chr_size[c] / 100.0) + 1) * n_samples[i] * sizeof(uint16_t);
        }
    }
    return n_sample_total;
}

int DataReader::read_depth100(sv& curr_sv, GcContent& gc)
{
    // this information is not useful when sv length is short
    // process only for >200bp SVs (or at include least two full 100-bp intervals)
    
    bool b_medfilt = true; // whether to do median filtering or not
    
    int n_inside =(int)curr_sv.end/100 - (int)curr_sv.pos/100;
    
    if (n_inside < 3)
        return -1;
    else if (n_inside < 5)
        b_medfilt = false; // no median filtering

    // now, let's forget edge detection and focus on the 'inside' portions
    // for small SVs ( <100kb ), it's okay to read the whole thing
    // for larger SVs, read two segments in the 'middle'
    
    if (curr_sv.len < 100000)
    {
        int startpos = 0;
        int endpos = 0;
        
        startpos = curr_sv.pos - 2000;
        
        if (startpos < 0)
            startpos = 1;
        
        endpos = curr_sv.end + 2000;
        if (endpos > (int) gc.chr_size[curr_sv.chrnum])
            endpos = (int) gc.chr_size[curr_sv.chrnum];
        
        int sample_idx = 0;
        
        int n_start = (startpos / 100);
        int n_end = (endpos / 100);
        int n_dp = n_end - n_start + 1;
        
        int x_start = (curr_sv.pos / 100) - n_start + 1;
        int x_end = (curr_sv.end + curr_sv.pos / 200) - n_start;
        int y_start = x_end+1;
        int y_end = (curr_sv.end/100) - n_start -1;
        
        for(int i=0; i<n_pileup; ++i)
        {
            uint64_t start_byte = multi_idx[i][0]; // This is the index position where dp100 record starts
            start_byte = chr_bytepos_dp100[i][curr_sv.chrnum];
            start_byte += n_start * n_samples[i] * sizeof(uint16_t);
            
            int n_dp_by_sample = n_samples[i] * n_dp;
            
            uint16_t *D = new uint16_t[n_dp_by_sample];
            uint16_t *D_filt = new uint16_t[n_dp_by_sample];
            
            pileups[i].seekg(start_byte);
            pileups[i].read_depth(D, n_dp_by_sample);
            
            // now do median filtering
            for(int j=0; j<n_samples[i]; ++j)
            {
                
            }
            // sample by sample sum
            for(int j=x_start; j<x_end; ++j)
            {
            
            }
            for(int j=y_start;j<y_end;++j)
            {
            
            }
            
            sample_idx += n_samples[i];
            
            delete [] D;
            delete [] D_filt;
        }
        
    }

    return 1;
}

void DataReader::read_pair_split(sv& curr_sv, std::vector<ReadStat>& rdstats, GcContent &gc)
{
    int sample_idx = 0;

    // Starting breakpoint
    int start_pos1 = curr_sv.pos - 500;
    int end_pos1 = curr_sv.pos + 500;

    if (start_pos1<1) start_pos1 = 1;
	if (end_pos1 >= curr_sv.end) end_pos1 = curr_sv.end-1;

    int start_pos2 = curr_sv.end - 500;
    int end_pos2 = curr_sv.end + 500;

    if (start_pos2 <= curr_sv.pos) start_pos2 = curr_sv.pos+1;
	if (end_pos2 > gc.chr_size[curr_sv.chrnum])  end_pos2 = gc.chr_size[curr_sv.chrnum]-1;
    
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
				printf("\nSample %d\n", sample_idx+k);
                uint32_t n_rp = 0;
                pileups[i].read_uint32(n_rp);
 
                for(uint32_t l=0; l<n_rp; ++l)
                {
                    readpair rp;
                    pileups[i].read_readpair(rp);
					if (rp.selfpos >= start_pos1 && rp.selfpos <= end_pos1)
					{
						if (rp.matequal>0)
						{
							if (rp.matepos >= start_pos2 && rp.matepos <= end_pos2)
							{
                                printf(">>>");
								if (rp.pairstr == 1 && rp.selfpos < rp.matepos)
									rdstats[sample_idx + k].n_pre_FR ++;
	 							else if (rp.pairstr == 2 && rp.selfpos <= rp.matepos)
									rdstats[sample_idx + k].n_pre_RF ++;
                                else if (rp.pairstr == 0 && rp.selfpos <= rp.matepos)
                                    rdstats[sample_idx + k].n_pre_FF ++;
                                
                            }
						}
						else
                            rdstats[sample_idx + k].n_pre_rp_missing ++;
                        printf("PRE_RP\t%d\t%d\t%d\t%u\t%d\n", rp.chrnum, rp.selfpos, rp.matepos, rp.matequal, rp.pairstr);

					}
 					if (b_overlap && rp.selfpos >= start_pos2 && rp.selfpos <= end_pos2)
					{
						if (rp.matequal > 0)
						{
                            if (rp.matepos >= start_pos1 && rp.matepos <= end_pos1)
                            {
                                printf(">>>");

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
                        printf("POST_RP\t%d\t%d\t%d\t%u\t%d\n", rp.chrnum, rp.selfpos, rp.matepos, rp.matequal, rp.pairstr);

					}
                }
				printf("\n");
                uint32_t n_sp = 0;
                pileups[i].read_uint32(n_sp);
                for(uint32_t l=0; l<n_sp; ++l)
                {
                    splitread sp;
                    pileups[i].read_splitread(sp);
					if (sp.pos >=start_pos1 && sp.pos <= end_pos1)
                    {
                        if (sp.sapos >= start_pos2 && sp.sapos <= end_pos2)
                        {
                            printf(">>>");

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
                        else if (sp.sapos == 0)
                        {
                            rdstats[sample_idx+k].n_pre_sp_missing ++;
                        }
                        printf("PRE_SP\t%d\t%d\t%d\t%d\t%d\n", sp.chrnum, sp.pos, sp.sapos, sp.firstclip, sp.secondclip);

                    }
                    if (b_overlap && sp.pos >=start_pos2 && sp.pos <= end_pos2)
                    {
                        if (sp.sapos >= start_pos1 && sp.sapos <= end_pos1)
                        {
                            printf(">>>");

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
                        else if (sp.sapos == 0)
                        {
                            rdstats[sample_idx+k].n_post_sp_missing ++;
                        }
                        printf("POST_SP\t%d\t%d\t%d\t%d\t%d\n", sp.chrnum, sp.pos, sp.sapos, sp.firstclip, sp.secondclip);

                    }
                }
                
            }
        }

        sample_idx += n_samples[i];
    }
    
    if (b_overlap) return;
    
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
                printf("\nSample %d\n", sample_idx+k);
                uint32_t n_rp = 0;
                pileups[i].read_uint32(n_rp);
                
                for(uint32_t l=0; l<n_rp; ++l)
                {
                    readpair rp;
                    pileups[i].read_readpair(rp);
     
                    if (rp.selfpos >= start_pos2 && rp.selfpos <= end_pos2)
                    {
                        printf(">>>");
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
                        printf("POST_RP\t%d\t%d\t%d\t%u\t%d\n", rp.chrnum, rp.selfpos, rp.matepos, rp.matequal, rp.pairstr);

                    }
                }
                printf("\n");
                uint32_t n_sp = 0;
                pileups[i].read_uint32(n_sp);
                for(uint32_t l=0; l<n_sp; ++l)
                {
                    splitread sp;
                    pileups[i].read_splitread(sp);
      
                    if (sp.pos >=start_pos2 && sp.pos <= end_pos2)
                    {
                        if (sp.sapos >= start_pos1 && sp.sapos <= end_pos1)
                        {
                            printf(">>>");
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
                        printf("POST_SP\t%d\t%d\t%d\t%d\t%d\n", sp.chrnum, sp.pos, sp.sapos, sp.firstclip, sp.secondclip);

                    }
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
    int p = pos*2 / gc.binsize;
    //    std::cerr << "pos " << pos << " p " << p << std::endl;
    int bin = gc.gc_array[chr][p];
    
    if (bin<20 && gc_factors[n][bin]>0.0001)
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
