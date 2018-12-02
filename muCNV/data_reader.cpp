//
//  multi_pileup.cpp
//  muCNV
//
//  Created by Goo Jun on 11/25/18.
//  Copyright © 2018 Goo Jun. All rights reserved.
//

#include "data_reader.h"
#include <math.h>

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

int DataReader::read_depth100(sv& curr_sv, std::vector< std::vector<double> > &dp100, GcContent& gc)
{
    // this information is not useful when sv length is short
    // process only for >300bp SVs
    
    if (curr_sv.len < 300 )
        return -1;
    
    // SV < 100kb : process all 100-bp intervals
    // SV >= 100kb : process 'around breakpoints' in original resolution and merge into 1kb blocks inside the SV

    if (curr_sv.len < 100000)
    {
        int startpos = 0;
        int endpos = 0;
        
        // put enough buffers before/after for median filtering
        startpos = curr_sv.pos - 2000;
        
        // return startpos to let the caller know where dp100 starts
        if (startpos < 0)
            startpos = 1;
        
        endpos = curr_sv.end + 2000;
        if (endpos > gc.chr_size[curr_sv.chrnum])
            endpos = gc.chr_size[curr_sv.chrnum];
    
        int sample_idx = 0;
        
        int n_start = (startpos / 100);
        int n_end = (endpos / 100);
        int n_dp = n_end - n_start + 1;
        
        for(int i=0; i<n_sample_total; ++i)
        {
            dp100[i].resize(n_dp);
        }
        
        for(int i=0; i<n_pileup; ++i)
        {
            uint64_t start_byte = multi_idx[i][0]; // This is the index position where dp100 record starts
			std::cerr << "startbyte " << start_byte << std::endl;
        
            start_byte = chr_bytepos_dp100[i][curr_sv.chrnum];
            
			std::cerr << "startbyte " << start_byte << std::endl;
            start_byte += n_start * n_samples[i] * sizeof(uint16_t);
			std::cerr << "startbyte " << start_byte << std::endl;
            
            int n_dp_by_sample = n_samples[i] * n_dp ;
            
            uint16_t *D = new uint16_t[n_dp_by_sample];
			std::cerr << "n_dp_by_sample " << n_dp_by_sample << std::endl;
            
            pileups[i].seekg(start_byte);
            pileups[i].read_depth(D, n_dp_by_sample);
			std::cerr << "read n_dp_by_sample " << n_dp_by_sample << std::endl;
         
            for(int j=0; j<n_dp; ++j)
            {
                for(int k=0; k<n_samples[i]; ++k)
                {
//					std::cerr << j << ", " << n_samples[i] << ", " << k << " : " << D[j*n_samples[i] + k ] << std::endl;
                    dp100[sample_idx + k][j] = (double)D[j*n_samples[i] + k] / 32.0;
                }
            }
            sample_idx += n_samples[i];

			std::cerr << "read n_dp_by_sample " << n_dp_by_sample << std::endl;
            delete [] D;
			std::cerr << "read n_dp_by_sample " << n_dp_by_sample << std::endl;
        }
        return startpos;
    }
    else if (curr_sv.len <= 1000000) 
    {
        int startpos = 0;
        int endpos = 0;
        
        // read start-2000 to start+2000 (100bp resolution)
        // read start+2000 to end+2000 (1kbp resolution)
        // read end-2000 to end+2000 (100bp resolution)
        
        return startpos;
    }
    else
    {
        int startpos = 0;
        int endpos = 0;
        // read start-2000 to start+2000 (10kbp resolution)
        // read start+2000 to end+2000 (10kbp resolution)
        // read end-2000 to end+2000 (10kbp resolution)
        return startpos;
    }


    return -1; // Should never reach here
    /*
    for(int i=0; i<n_pileup; ++i)
    {

            uint16_t dp100;
            for(int j=0;j<N;++j)
            {
                for(int i=0; i<n_sample; ++i)
                {
                    pup.read_depth(&dp100, 1);
                    dpsum[i] += dp100;
                    n_dp[i] +=1;
                }
            }
        }
    }
    
    }
    for(int i=0; i<n_sample; ++i)
    {
        SampleStat s;
        pup.read_sample_stat(s);
        printf("Sample %d, AVG DP: %f, STdev: %f, AVG ISIZE: %f, STdev: %f \n", i, s.avg_dp, s.std_dp, s.avg_isize, s.std_isize);
    }
    
    GcContent gc;
    gc.initialize(gc_file);
    
    for(int i=0; i<n_sample; ++i)
    {
        printf("GC-factors for sapmle %d:\n", i);
        std::vector<double> gc_factor (gc.num_bin);
        pup.read_gc_factor(gc_factor, gc.num_bin);
        for(int j=0; j<gc.num_bin; ++j)
        {
            printf("GC-bin %d: %f\n", j, gc_factor[j]);
        }
    }
    
    idx_file.read_uint64(curr_idx);
    printf("index position %d, tellg position %lu\n", (int)curr_idx, (unsigned long)pup.tellg());
    
    std::vector<uint64_t> dpsum (n_sample, 0);
    std::vector<uint64_t> n_dp (n_sample, 0);
    
    for(int c=1; c<=gc.num_chr; ++c)
    {
        int N = ceil((double)gc.chr_size[c] / 100.0) + 1;
        uint16_t dp100;
        for(int j=0;j<N;++j)
        {
            for(int i=0; i<n_sample; ++i)
            {
                pup.read_depth(&dp100, 1);
                dpsum[i] += dp100;
                n_dp[i] +=1;
            }
        }
    }
    for(int i=0; i<n_sample; ++i)
    {
        printf("Sample %d, average DP100: %d\n", i, (int)round((double)dpsum[i]/n_dp[i]/32.0));
    }
    
  
    
    pup.close();
    idx_file.close();
    
    int n_var = 0;
    
    var_file.read_int32(n_sample);
    var_file.read_int32(n_var);
    
    printf("Variant File, n_sample: %d, n_var : %d\n", n_sample, n_var);
    
    printf("Sample ID(s):");
    for(int i=0; i<n_sample; ++i)
    {
        char buf[256];
        
        var_file.read_sample_id(buf);
        printf("\t%s", buf);
    }
    printf("\n");
    
    for(int j=0;j<n_var;++j)
    {
        vec_sv[j].print();
        for(int i=0; i<n_sample; ++i)
        {
            uint16_t dp;
            var_file.read_depth(&dp, 1);
            printf("\t%f", (dp/32.0));
        }
        printf("\n");
    }
    var_file.close();
    */
}

void DataReader::read_pair_split(sv& curr_sv, std::vector< std::vector<readpair> > &dvec_rp, std::vector< std::vector<splitread> > &dvec_sp)
{
    int sample_idx = 0;
    
    // Starting breakpoint
    int startpos = curr_sv.pos - 500;
    int endpos = curr_sv.pos + 500;
    
    if (startpos<1) startpos = 1;
    if (endpos>=curr_sv.end) endpos = curr_sv.end;
    
    
    uint64_t start_idx = chr_idx_rp[curr_sv.chrnum] + (int)startpos/10000;
    uint64_t end_idx = chr_idx_rp[curr_sv.chrnum] + (int)endpos/10000;
    
    for(int i=0; i<n_pileup; ++i)
    {
        pileups[i].seekg(multi_idx[i][start_idx]);

        for(uint64_t j=start_idx; j<=end_idx; ++j)
        {
            
            for(int k=0; k<n_samples[i]; ++k)
            {
                uint32_t n_rp = 0;
                pileups[i].read_uint32(n_rp);
 
                for(int l=0; l<n_rp; ++l)
                {
                    readpair rp;
                    pileups[i].read_readpair(rp);
                    fprintf(stderr, "\t%d\t%d\t%d\t%u\t%d\n", rp.chrnum, rp.selfpos, rp.matepos, rp.matequal, rp.pairstr);
                }
                // check whether rp.selfpos is between startpos and endpos
                // if so, check strand matequal/matepos
                
                uint32_t n_sp = 0;
                pileups[i].read_uint32(n_sp);
                for(int l=0; l<n_sp; ++l)
                {
                    splitread sp;
                    pileups[i].read_splitread(sp);
                    fprintf(stderr, "\t%d\t%d\t%d\t%d\t%d\n", sp.chrnum, sp.pos, sp.sapos, sp.firstclip, sp.secondclip);
                    
                }
                
            }
        }

        sample_idx += n_samples[i];
    }
    
    // Ending breakpoint
    sample_idx = 0;
    startpos = curr_sv.end - 500;
    endpos = curr_sv.end + 500;
    
    if (startpos<curr_sv.pos) startpos = curr_sv.pos+1;
    
    start_idx = chr_idx_rp[curr_sv.chrnum] + (int)startpos/10000;
    end_idx = chr_idx_rp[curr_sv.chrnum] + (int)endpos/10000;
    
    for(int i=0; i<n_pileup; ++i)
    {
        pileups[i].seekg(multi_idx[i][start_idx]);

        for(uint64_t j=start_idx; j<=end_idx; ++j)
        {
            
            for(int k=0; k<n_samples[i]; ++k)
            {
                uint32_t n_rp = 0;
                pileups[i].read_uint32(n_rp);
                
                for(int l=0; l<n_rp; ++l)
                {
                    readpair rp;
                    pileups[i].read_readpair(rp);
                    printf("\t%d\t%d\t%d\t%u\t%d\n", rp.chrnum, rp.selfpos, rp.matepos, rp.matequal, rp.pairstr);
                }
                
                uint32_t n_sp = 0;
                pileups[i].read_uint32(n_sp);
                for(int l=0; l<n_sp; ++l)
                {
                    splitread sp;
                    pileups[i].read_splitread(sp);
                    printf("\t%d\t%d\t%d\t%d\t%d\n", sp.chrnum, sp.pos, sp.sapos, sp.firstclip, sp.secondclip);
                    
                }
                
            }
        }
        
        sample_idx += n_samples[i];
    }
    
    return;
}

void DataReader::read_var_depth(int k, std::vector<double> &var_dp)
{
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
