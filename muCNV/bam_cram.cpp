//
//  bamcram.cpp
//  muCNV
//
//  Created by Goo Jun on 11/29/17.
//  Copyright © 2017 Goo Jun. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <string>

#include "pileup.h"
#include "bam_cram.h"

#define IS_PROPERLYPAIRED(bam) (((bam)->core.flag&(BAM_FPAIRED|BAM_FPROPER_PAIR)) == (BAM_FPAIRED|BAM_FPROPER_PAIR) && !((bam)->core.flag&BAM_FUNMAP))

typedef enum { READ_UNKNOWN = 0, READ_1 = 1, READ_2 = 2 } readpart;

int get_cigar_clippos(std::string &cigar_str)
{
    int lclip = 0;
    int rclip = 0;
    int j=0;
    while(cigar_str[j] >= '0' && cigar_str[j] <='9' && j<(int)cigar_str.length())
        ++j;
    if (j==0 || j==(int)cigar_str.length())
    {
        lclip = 0; //something wrong
    }
    else if (cigar_str[j-1] == 'S')
    {
        lclip = atoi(cigar_str.substr(0,j).c_str());
    }
    
    j=(int)cigar_str.length()-1;
    
    if (cigar_str.back() == 'S')
    {
        --j;
        while(cigar_str[j] >= '0' && cigar_str[j] <='9' && j>=0)
            --j;
        if (j<0) // all numbers
        {
            rclip = atoi(cigar_str.substr(0,cigar_str.length()-1).c_str());
        }
        else
        {
            rclip = atoi(cigar_str.substr(j+1,cigar_str.length()-1).c_str());
        }
    }
    
    if (rclip >= lclip && rclip >= 10) // arbitrary cutoff, 15
        return -rclip;
    else if (lclip>rclip && lclip >= 10)
        return lclip;
    else
        return 0;
}

bool process_split(std::string &t, splitread &new_sp, int32_t tid, bool strand)
{
    int i=1;
    
    int32_t chr = 0;
    std::string cigar_str;
    std::string pos_str;

    if (t.substr(0,3) == "chr")
    {
        if (t[3]>='1' && t[3] <= '9')
        {
            if (t[4] == ',')
            {
                chr = t[3]-'0';
                i=5;
            }
            else if (t[4]>='0' && t[4]<='9' && t[5]==',')
            {
                chr = (t[3]-'0')*10 + (t[5]-'0');
                i=6;
            }
        }
        else if (t[3] == 'X' && t[4] == ',')
        {
            chr = 23;
            i=5;
        }
        else if (t[3] == 'Y' && t[4] == ',')
        {
            chr = 24;
            i=5;
        }
        else
        {
            return false;
        }
        
        if (chr != tid+1)
        {
            return false;
        }

		std::string::size_type c;
		std::string::size_type d;

		new_sp.sapos = stoi(t.substr(i), &c);
	//	DMSG("POS_STR: " << t.substr(i) << " SAPOS:" <<new_sp.sapos);
        
		if (new_sp.sapos < 1) return false;

		c=c+i+1;
		if (t[c] == '+')
		{
			if (!strand)
				return false;
		}
		else if (t[c] == '-')
		{
			if (strand)
				return false;
		}
		// DMSG("STRAND : " << t[c]);
		c+=2;

		d = t.find(',', c);
        if (d == std::string::npos) return false;
		cigar_str = t.substr(c,d-c);
	//	DMSG("CIGAR: " << cigar_str);
        
		if (cigar_str.length() == 0 ) return false;

        new_sp.secondclip = get_cigar_clippos(cigar_str);
	    return true;
    }
	else
	{
		return false;
	}
}

static int read_bam(void *data, bam1_t *b) // read level filters better go here to avoid pileup
{
    aux_t *aux = (aux_t*)data; // data in fact is a pointer to an auxiliary structure
    int ret;
    
    while (1)
    {
        ret = aux->iter ? sam_itr_next(aux->fp, aux->iter, b) : sam_read1(aux->fp, aux->hdr, b);
        
        if ( ret<0 ) break;
        if ( b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP) ) continue;
        if ( (int)b->core.qual < aux->min_mapQ ) continue;
		if ( aux->min_len && bam_cigar2qlen(b->core.n_cigar, bam_get_cigar(b)) < aux->min_len ) continue;
        if (in_centrome(b->core.tid+1, b->core.pos)) continue;
            
//		fprintf(stderr, "%s\ttid:%d\tpos:%d\tmtid:%d\tmtpos:%d",bam_get_qname(b), b->core.tid, b->core.pos, b->core.mtid, b->core.mpos);
        if ((b->core.flag & BAM_FPAIRED ) && !(b->core.flag & BAM_FSUPPLEMENTARY))
        {
            if (IS_PROPERLYPAIRED(b) && b->core.pos>0 && b->core.mtid == b->core.tid && b->core.isize>0)
            {
                // properly paired read, add to insert size distribution

                // get average isize statistics only from properly paired pairs
                // and also for isize>0 cases (do not double count)
                aux->sum_isz += b->core.isize;
                aux->sumsq_isz += (b->core.isize) * (b->core.isize);
                aux->n_isz += 1;
            }
            else 
            {
                // if not properly paired
                
                if (b->core.qual>=10)
                {
                    // add to read pair set
                    readpair new_rp;
                    new_rp.chrnum = b->core.tid+1;
					new_rp.matepos= 0;
					new_rp.matequal = 0;
                    new_rp.selfpos = b->core.pos;
                    new_rp.pairstr = (b->core.flag & BAM_FREVERSE) ? 2 : 0;

                    // Write only once for a read pair, record only (self-mate), not (mate-self)
                    if (b->core.tid == b->core.mtid && b->core.pos <= b->core.mpos)
                    {
                        uint8_t *aux_mq = bam_aux_get(b,"MQ");
						if (aux_mq != NULL)
                        {
                            new_rp.matequal = (int8_t) bam_aux2i(aux_mq);
                        }
                        else
                        {
                            new_rp.matequal = 0;
                        }
                        new_rp.matepos = b->core.mpos;
                        new_rp.pairstr += (b->core.flag & BAM_FMREVERSE) ? 1 : 0;
                        aux->n_rp++;
                        (*(aux->p_vec_rp)).push_back(new_rp);
                    }
                    else if (b->core.flag & BAM_FMUNMAP) // mate is unmapped, but self has good MQ  - insertion or inversion
                    {
                        // TODO: Check - Is this useful at all?
                        new_rp.matequal = 0;
                        new_rp.matepos = -1;
                        aux->n_rp++;
                        (*(aux->p_vec_rp)).push_back(new_rp);
					//	fprintf(stderr, "\tMateQualZero\n");
                    }
                }
            }
        }
        
        int16_t lclip = 0;
        int16_t rclip = 0;
        
        // Check whether there's soft clip for all reads with >1 cigar ops
        if ((b->core.n_cigar > 1) && b->core.qual >= 10) // min_MAPQ 10 for softclips
        {
            uint32_t *cigar  = bam_get_cigar(b);
            int ncigar = b->core.n_cigar;

            if (bam_cigar_op(cigar[0]) == BAM_CSOFT_CLIP)
            {
                lclip = bam_cigar_oplen(cigar[0]);
                sclip new_clip;
                new_clip.chrnum = b->core.tid + 1;
                new_clip.pos = b->core.pos + lclip;
                (*(aux->p_vec_lclip)).push_back(new_clip);
            }
            if (bam_cigar_op(cigar[ncigar-1]) == BAM_CSOFT_CLIP)
            {
                rclip = bam_cigar_oplen(cigar[ncigar-1]);
                sclip new_clip;
                new_clip.chrnum = b->core.tid + 1;
                new_clip.pos = bam_endpos(b) - rclip;
                (*(aux->p_vec_rclip)).push_back(new_clip);
            }
        }
        
		if ((b->core.flag & BAM_FSUPPLEMENTARY) && !in_centrome(b->core.tid+1, b->core.pos) )
		{
			uint8_t *aux_sa = bam_aux_get(b, "SA");
			char* p_sa;
			if (aux_sa && (p_sa = bam_aux2Z(aux_sa)))
			{
				std::string str_sa = std::string(p_sa);
                
				splitread new_sp;
				new_sp.chrnum = b->core.tid + 1;
				new_sp.pos = b->core.pos;
				new_sp.firstclip = 0;
				new_sp.secondclip = 0;

				int ncigar = b->core.n_cigar;
				uint32_t *cigar  = bam_get_cigar(b);

                if (bam_cigar_op(cigar[0]) == BAM_CSOFT_CLIP)
				{
					lclip = bam_cigar_oplen(cigar[0]);
				}
				if (bam_cigar_op(cigar[ncigar-1]) == BAM_CSOFT_CLIP)
				{
					rclip = bam_cigar_oplen(cigar[ncigar-1]);
				}
                
                if (rclip > lclip && rclip > 15) // arbitrary cutoff, 15
                    new_sp.firstclip = -rclip;
                else if (lclip>rclip && lclip > 15)
                    new_sp.firstclip = lclip;
                
                if (new_sp.firstclip != 0)
                // Process split read only if the read has soft-clipped ends
				{
					if (process_split(str_sa, new_sp, b->core.tid, !(b->core.flag & BAM_FREVERSE)) && !in_centrome(b->core.tid+1, new_sp.sapos))
					{
						aux->n_sp++;
						(*(aux->p_vec_sp)).push_back(new_sp);
					}
				}
			}
		}
        break;
    }
    return ret;
}


void BamCram::initialize_sequential(std::string &bname, GcContent &gc)
{
    data = (aux_t**)calloc(1, sizeof(aux_t*));
    data[0] = (aux_t*)calloc(1, sizeof(aux_t));
    data[0]->fp = hts_open(bname.c_str(), "r");
    
    data[0]->sum_isz = 0;
    data[0]->sumsq_isz = 0;
    data[0]->n_isz = 0;
    
    
    data[0]->min_mapQ = 1;
    data[0]->min_len = 0;
    data[0]->hdr = sam_hdr_read(data[0]->fp);;
    
    if (data[0]->hdr == NULL)
    {
        std::cerr << "Cannot open CRAM/BAM header" << std::endl;
        exit(1);
    }
    idx = NULL;

    gc_factor.resize(gc.num_bin);
    gc_count.resize(gc.num_bin);
    
    depth_interval.resize(gc.num_chr + 1);
    
    for(int i=1; i<=gc.num_chr; ++i)
    {
        depth_interval[i] = (uint16_t *) calloc(gc.n_interval[i], sizeof(uint16_t));
        
        DMSG("chr " << i << " bin size " << gc.n_interval[i]);
    }
}

void BamCram::read_depth_sequential(Pileup& pup, GcContent& gc, std::vector<breakpoint> &vec_bp, std::vector<sv> &vec_sv)
{
    // vec_bp should have been sorted beforehand
    int tid = -1, pos = -1;
    
    if (vec_bp[0].bptype > 0)
    {
        std::cerr << "Error: Merged interval's earliest position is SV-end, not SV-start." << std::endl;
        for(int i=0;i<20;++i)
        {
            std::cerr << vec_bp[i].chrnum << "\t" << vec_bp[i].pos << "\t" << vec_bp[i].bptype << std::endl;
        }
    }
    
    size_t nxt = 0;
   
    uint64_t sum_dp = 0;
    uint64_t sumsq_dp = 0;
    uint64_t n_dp = 0;
    
    std::vector<uint64_t> gc_sum;
    std::vector<uint64_t> gc_cnt;
    
    std::vector<int> dp_list;
    
    std::vector<double> dp_sum; // One interval (start-end) will add one entry on these
    std::vector<double> gc_dp_sum;
    std::vector<int> dp_cnt;

    data[0]->p_vec_rp = &(vec_rp);
    data[0]->p_vec_sp = &(vec_sp);
    data[0]->p_vec_rclip = &(vec_rclip);
    data[0]->p_vec_lclip = &(vec_lclip);
    
    data[0]->sum_isz = 0;
    data[0]->sumsq_isz = 0;
    data[0]->n_isz = 0;

	data[0]->n_rp = 0;
	data[0]->n_sp = 0;

    gc_sum.resize(gc.num_bin);
    gc_cnt.resize(gc.num_bin);
    
    for(int i=0;i<gc.num_bin;++i)
    {
        gc_sum[i] = 0;
        gc_cnt[i] = 0;
    }    
    
    int* n_plp = (int *) calloc(1, sizeof(int));;
    const bam_pileup1_t **plp = (const bam_pileup1_t **) calloc(1, sizeof(bam_pileup1_t*));

    bam_mplp_t mplp = bam_mplp_init(1, read_bam, (void**) data);
	bam_mplp_set_maxcnt(mplp, 64000);
    
    
    int sum100=0;
    int n100=0;
    
    int prev_chrnum = 1;
    int prev_pos = 1;
    
	std::cerr << "processing chr 1" << std::endl;
    
    while(bam_mplp_auto(mplp, &tid, &pos, n_plp, plp)>0 && tid < gc.num_chr)
    {
        int chrnum = tid+1;
        
        if (chrnum>prev_chrnum)
        {
            // clear stats
			std::cerr << "now processing chr " << chrnum << std::endl;
			std::cerr << "n_rp : " << data[0]->n_rp << " n_sp: " << data[0]->n_sp << std::endl;
            
            // Update 100-bp depth
            if (n100>0)
            {
               	int val = round((double) (sum100 * 32) / n100); // Now Depth100 stores (depth*32) as value
                if (val>65535) val=65535; // handle overflow, though unlikely
                
                //TODO: fix hard-coded 100
                depth_interval[prev_chrnum][prev_pos/100] = (uint16_t) val;
            }
			std::cerr << "n_isz : " << data[0]->n_isz << ", sum_isz : " << data[0]->sum_isz << ", sumsq_isz : " << data[0]->sumsq_isz << std::endl;

            // TODO: encapsulate this code
            // gc_array[chr][1] would contain 
			uint8_t bin = gc.gc_array[prev_chrnum][(int) prev_pos / gc.interval_dist];
			if (bin<gc.num_bin)
			{
				gc_sum[bin] += sum100;
				gc_cnt[bin] += n100;
			}
	
            sum100=0;
            n100=0;
            prev_chrnum = chrnum;
	
        }
        else if (pos/100 > prev_pos/100)
        {
            // Update 100-bp depth
            if (n100>0)
            {
                int val = round((double) (sum100 * 32.0) / n100);
                if (val>65535) val=65535; // handle overflow
                                          // TODO: fix hard-coded 100
                depth_interval[chrnum][prev_pos/100] = (uint16_t) val;
            }
            
            // TODO: encapsulate this code
			uint8_t bin = gc.gc_array[chrnum][(int)prev_pos / gc.interval_dist];
			if (bin<gc.num_bin)
			{
				gc_sum[bin] += sum100;
				gc_cnt[bin] += n100;
			}
	
            sum100=0;
            n100=0;
        }
        
        prev_pos = pos; // maybe use array idx instead of pos ?
        breakpoint curr_bp;
        
        curr_bp.chrnum = chrnum;
        curr_bp.pos = pos;
        
        while (nxt<vec_bp.size() && vec_bp[nxt] <= curr_bp)
        {
            if (vec_bp[nxt].bptype == 0)
            {
            	dp_list.push_back(vec_bp[nxt].idx);
            }
            else
            {                
                std::vector<int>::iterator dp_it = find(dp_list.begin(), dp_list.end(), vec_bp[nxt].idx) ;
				if (dp_it == dp_list.end())
				{
                    // Error: idx to be removed is not in dp_list
					int idx_nxt = vec_bp[nxt].idx;
					std::cerr << "Error, idx " << vec_bp[nxt].idx << " bptype " << vec_bp[nxt].bptype << " bp pos " << vec_bp[nxt].pos << ", sv " << vec_sv[idx_nxt].chrnum << ":" << vec_sv[idx_nxt].pos << "-" << vec_sv[idx_nxt].end << std::endl;
					std::cerr << " dp_list [";
					for(int ii=0;ii<(int)dp_list.size();++ii)
						std::cerr << dp_list[ii] << " " ;
					std::cerr << "]" << std::endl;
					exit(1);
				}
				else
                {
					*dp_it=dp_list.back();
					dp_list.pop_back();
				}
            }

            nxt++;
        }
        
        int m=0;
        for(int j=0;j<n_plp[0];++j)
        {
            const bam_pileup1_t *p = plp[0]+j;
            if (p->is_del || p->is_refskip )
                ++m;
        }
        int dpval = n_plp[0]-m;
        
        sum100 += dpval;
        n100 ++;
        
        // for whole-genome average
        sum_dp += dpval;
        sumsq_dp += dpval * dpval;
        n_dp += 1;
        
		for(size_t ii=0; ii<dp_list.size(); ++ii)
		{
			/*
			if (dp_list[ii]<0 || dp_list[ii]>=vec_sv.size())
			{
				std::cerr << "Error, " << dp_list[ii] << " is out of bound " << std::endl;
			}
			*/
			vec_sv[dp_list[ii]].dp_sum += dpval;
			vec_sv[dp_list[ii]].n_dp += 1;
		}
    }
    
    // TODO: refactor this..

    stat.avg_dp = (double) sum_dp / n_dp;
    stat.std_dp = sqrt(((double)sumsq_dp / n_dp - (stat.avg_dp * stat.avg_dp)));
    
    // Sort lclip & rclip
    // We can assume readpairs and splitreads are sorted, because they're added according to the read orders
    sort(vec_lclip.begin(), vec_lclip.end());
    sort(vec_rclip.begin(), vec_rclip.end());
    
    std::vector<double> gc_avg (gc.num_bin, 0);
    
    for(int i=2; i<gc.num_bin-2; ++i)
    {
        if (gc_cnt[i] + (gc_cnt[i-1] + gc_cnt[i+1])/2.0 + (gc_cnt[i-2]+gc_cnt[i+2])/4.0 > 100)
        {
            gc_avg[i] = ((double)gc_sum[i] + 0.5*((double)gc_sum[i-1] + gc_sum[i+1]) + 0.25 * ((double)gc_sum[i-2] + gc_sum[i+2])) / (gc_cnt[i] +0.5*(gc_cnt[i-1]+gc_cnt[i+1]) + 0.25*(gc_cnt[i-2] + gc_cnt[i+2])) ;
        }
        else
        {
            gc_avg[i] = -1;
        }
    }
    // padding at the end
    gc_avg[0] = gc_avg[1] = gc_avg[2];
    gc_avg[gc.num_bin-1] = gc_avg[gc.num_bin-2] = gc_avg[gc.num_bin-3];
    
    // Calculate GC-curve and save it with original DP to maximize information preservation instead of storing GC-corrected depths only
    for(int i=0;i<gc.num_bin;++i)
    {
        if (gc_avg[i]>0)
        {
            gc_factor[i] = stat.avg_dp / gc_avg[i];
 //           stat.avg_dp / ((double)gc_sum[i] / gc_cnt[i]) ; // multiplication factor
        }
        else
        {
            gc_factor[i] = 0;
        }
    }

	std::cerr << "n_rp : " << data[0]->n_rp << " n_sp: " << data[0]->n_sp << std::endl;
    stat.avg_isize = data[0]->sum_isz / (double)data[0]->n_isz;
    stat.std_isize = sqrt((double)data[0]->sumsq_isz /(double)data[0]->n_isz - ((double)stat.avg_isize * stat.avg_isize));
    
    std::cerr << "n_isz : " << data[0]->n_isz << ", sum_isz : " << data[0]->sum_isz << ", sumsq_isz : " << data[0]->sumsq_isz << std::endl;
	std::cerr << "avg_isz : " << stat.avg_isize << ", std_isize : " << stat.std_isize << std::endl;

    sam_itr_destroy(data[0]->iter);
    free(plp); free(n_plp);
    bam_mplp_destroy(mplp);
}


void BamCram::postprocess_depth(std::vector<sv> &vec_sv)
{
    // TODO: make this on-the-fly, not post-process, this is a remnant of 'normalized' and 'GC-corrected' version
    // Update average depth of each SV
    for(int i=0;i<(int)vec_sv.size(); ++i)
    {
        if (vec_sv[i].n_dp>0)
        {
            int val = round((double)(vec_sv[i].dp_sum * 32)/vec_sv[i].n_dp);
            if (val>65535)
                vec_sv[i].dp = 65535;
            else
                vec_sv[i].dp = (uint16_t) val;
        }
        else
        {
            vec_sv[i].dp = 0;
        }
    }
}
/*

void BamCram::initialize(std::string &bname)
{
    
    data = (aux_t**)calloc(1, sizeof(aux_t*));
    data[0] = (aux_t*)calloc(1, sizeof(aux_t));
    data[0]->fp = hts_open(bname.c_str(), "r");
    
    int rf = SAM_FLAG |  SAM_MAPQ | SAM_QUAL| SAM_POS | SAM_SEQ | SAM_CIGAR| SAM_TLEN | SAM_RNEXT | SAM_PNEXT;
    if (hts_set_opt(data[0]->fp, CRAM_OPT_REQUIRED_FIELDS, rf))
    {
        std::cerr << "Failed to set CRAM_OPT_REQUIRED_FIELDS value" << std::endl;
        exit(1);
    }
    if (hts_set_opt(data[0]->fp, CRAM_OPT_DECODE_MD, 0))
    {
        fprintf(stderr, "Failed to set CRAM_OPT_DECODE_MD value\n");
        exit(1);
    }
    
    data[0]->min_mapQ = 1;
    data[0]->min_len = 0;
    data[0]->hdr = sam_hdr_read(data[0]->fp);;
    
    if (data[0]->hdr == NULL)
    {
        std::cerr << "Cannot open CRAM/BAM header" << std::endl;
        exit(1);
    }
    
    hts_idx_t* tmp_idx = sam_index_load(data[0]->fp, bname.c_str());
    if (tmp_idx == NULL)
    {
        std::cerr << "Cannot open CRAM/BAM index" << std::endl;
        exit(1);
    }
    idx = tmp_idx;
}
*/
/*
 static readpart which_readpart(const bam1_t *b)
 {
 if ((b->core.flag & BAM_FREAD1) && !(b->core.flag & BAM_FREAD2)) {
 return READ_1;
 } else if ((b->core.flag & BAM_FREAD2) && !(b->core.flag & BAM_FREAD1)) {
 return READ_2;
 } else {
 return READ_UNKNOWN;
 }
 }
 */

// This function reads a BAM alignment from one BAM file.
/*
 static int read_bam_basic(void *data, bam1_t *b) // read level filters better go here to avoid pileup
 {
 aux_t *aux = (aux_t*)data; // data in fact is a pointer to an auxiliary structure
 int ret;
 
 while (1)
 {
 ret = aux->iter ? sam_itr_next(aux->fp, aux->iter, b) : sam_read1(aux->fp, aux->hdr, b);
 
 if ( ret<0 ) break;
 if ( b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP) ) continue;
 if ( (int)b->core.qual < aux->min_mapQ ) continue;
 // Nov 29, 2017, commented out
 //if ( aux->min_len && bam_cigar2qlen(b->core.n_cigar, bam_get_cigar(b)) < aux->min_len ) continue;
 //        std::cerr << bam_get_qname(b) << "/" << which_readpart(b) << " insert size : " << b->core.isize << std::endl;
 if (IS_PROPERLYPAIRED(b) && (b->core.qual > 20) && (b->core.tid == b->core.mtid) && b->core.mpos>0 )
 {
 if (b->core.isize < 10000 && b->core.isize>-10000)
 (*(aux->isz_list))[0].push_back(abs(b->core.isize));
 }
 break;
 }
 return ret;
 }
 */
/*
 // This function reads a BAM alignment from one BAM file.
 static int read_bam(void *data, bam1_t *b) // read level filters better go here to avoid pileup
 {
 aux_t *aux = (aux_t*)data; // data in fact is a pointer to an auxiliary structure
 int ret;
 
 while (1)
 {
 ret = aux->iter ? sam_itr_next(aux->fp, aux->iter, b) : sam_read1(aux->fp, aux->hdr, b);
 
 if ( ret<0 ) break;
 if ( b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP) ) continue;
 if ( (int)b->core.qual < aux->min_mapQ ) continue;
 // Nov 29, 2017, commented out
 //if ( aux->min_len && bam_cigar2qlen(b->core.n_cigar, bam_get_cigar(b)) < aux->min_len ) continue;
 // Filters: MQ>0, mate mapped to the same chr, isize>0, pos>0
 
 if ((int)b->core.qual>0 && b->core.tid == b->core.mtid && b->core.isize !=0 && b->core.mpos > 0)
 {
 if (((b->core.flag & BAM_FREVERSE) == BAM_FREVERSE) == ((b->core.flag & BAM_FMREVERSE) == BAM_FMREVERSE))
 {
 //            std::cerr << "REVERSED tid " << b->core.tid << " pos " << b->core.pos << " qual " <<(int)b->core.qual << " insert size " << b->core.isize << " mtid " << b->core.mtid <<  " mpos " << b->core.mpos << std::endl;
 
 for(set<int>::iterator it=(*(aux->isz_set)).begin(); it!=(*(aux->isz_set)).end(); ++it)
 {
 (*(aux->rev_isz_list))[*it].push_back(b->core.isize);
 //                    (*(aux->rev_pos_list))[*it].push_back(b->core.pos);
 }
 }
 else
 {
 if (! IS_PROPERLYPAIRED(b) )
 {
 //            std::cerr << "NOTPROPER tid " << b->core.tid << " pos " << b->core.pos << " qual " <<(int)b->core.qual << " insert size " << b->core.isize << " mtid " << b->core.mtid <<  " mpos " << b->core.mpos << std::endl;
 for(set<int>::iterator it=(*(aux->isz_set)).begin(); it!=(*(aux->isz_set)).end(); ++it)
 {
 (*(aux->isz_list))[*it].push_back(b->core.isize);
 //                    (*(aux->pos_list))[*it].push_back(b->core.pos);
 }
 }
 else
 {
 for(set<int>::iterator it=(*(aux->isz_set)).begin(); it!=(*(aux->isz_set)).end(); ++it)
 {
 (*(aux->isz_sum))[*it] += abs(b->core.isize);
 (*(aux->isz_cnt))[*it] += 1;
 }
 
 }
 }
 }
 break;
 }
 return ret;
 }
 */ // this version is commented out in Sep. 2018.

/*
// Get (overlapping) list of SV intervals, return average depth and GC-corrected average depth on intervals
void BamCram::read_depth(std::vector<sv> &m_interval, std::vector<std::string> &G )
{
	char reg[100];
	
	std::string chr = m_interval[0].chr;
	int chrnum = m_interval[0].chrnum;
	int startpos = 1;
	int endpos = m_interval[0].end;
	int n = (int) m_interval.size();
	int gap = med_isize-(avg_rlen/2);
	std::vector<breakpoint> bp;
	bp.resize(n*4);  // 4 checkpoints per interval

	for(int i=0; i<n; ++i)
	{
		// START - GAP
		bp[i*4].pos = m_interval[i].pos > gap ? m_interval[i].pos - gap : 1;
		bp[i*4].type = 0;
		bp[i*4].idx = i;
		
		// START
		bp[i*4+1].pos = m_interval[i].pos;
		bp[i*4+1].type = 1;
		bp[i*4+1].idx = i;
		
		// END
		bp[i*4+2].pos = m_interval[i].end;
		bp[i*4+2].type = 2;
		bp[i*4+2].idx = i;
		
		//END + GAP
		bp[i*4+3].pos = m_interval[i].end + gap;
		bp[i*4+3].type = 3;
		bp[i*4+3].idx = i;
		
		if (m_interval[i].end <= m_interval[i].pos )
		{
			std::cerr << "Error! interval end point " << m_interval[i].end << " is before start point " << m_interval[i].pos << std::endl;
		}
		
		if (m_interval[i].end + gap > endpos)
		{
			endpos = m_interval[i].end + gap;
		}
	}
	
	sort(bp.begin()+1, bp.end());

	startpos = bp[0].pos;
	endpos += 1;
	
	sprintf(reg, "chr%s:%d-%d", chr.c_str(), startpos, endpos);
//	sprintf(reg, "%s:%d-%d", chr.c_str(), startpos, endpos); // Check BAM/CRAM header for list of CHRs first?

	data[0]->iter = sam_itr_querys(idx, data[0]->hdr, reg);
	if (data[0]->iter == NULL)
	{
		std::cerr << reg << std::endl;
		std::cerr << "Can't parse region " << reg << std::endl;
		exit(1);
	}

	int tid, pos;
	
	if (bp[0].idx != 0)
	{
		std::cerr << "Error: Merged interval's first element is not the earliest." << std::endl;
		std::cerr << "first " << bp[0].idx  << " pos " << bp[0].pos << " second " << bp[1].idx << " pos " << bp[1].pos << std::endl;
	}
	if (bp[0].type > 1)
	{
		std::cerr << "Error: Merged interval's earliest position is SV-end, not SV-start." << std::endl;
	}
	int nxt = 1;
	
	// For average DP
	set<int> dp_set;
	std::vector<double> dp_sum (n,0);
	std::vector<double> gc_dp_sum (n,0);
	std::vector<int> dp_cnt (n,0);
	
	// For average Insert Size
	set<int> isz_set;

	std::vector< std::vector<int> > isz_list;
//	std::vector< std::vector<int> > pos_list;

	std::vector< std::vector<int> > rev_isz_list;
//	std::vector< std::vector<int> > rev_pos_list;

	std::vector<double> isz_sum(n,0);
	std::vector<int> isz_cnt(n,0);

	isz_list.resize(n);
//	pos_list.resize(n);
	rev_isz_list.resize(n);
//	rev_pos_list.resize(n);
	
	isz_set.insert(bp[0].idx);

	data[0]->isz_set = &isz_set;
	data[0]->isz_list = &isz_list;
//	data[0]->pos_list = &pos_list;
	data[0]->rev_isz_list = &rev_isz_list;
//	data[0]->rev_pos_list = &rev_pos_list;
	data[0]->isz_sum = &isz_sum;
	data[0]->isz_cnt = &isz_cnt;
	
	int* n_plp = (int *) calloc(1, sizeof(int));;
	const bam_pileup1_t **plp = (const bam_pileup1_t **) calloc(1, sizeof(bam_pileup1_t*));
	bam_mplp_t mplp = bam_mplp_init(1, read_bam, (void**) data);

//	const bam_pileup1_t *p;
	
	while(bam_mplp_auto(mplp, &tid, &pos, n_plp, plp)>0)
	{
		while (nxt < n*4 && pos>=bp[nxt].pos)
		{
			switch(bp[nxt].type)
			{
				case 0:
					//Pre-gap start
					isz_set.insert(bp[nxt].idx);
					nxt++;
					break;
				case 1:
					//Pre-gap end, interval start
					dp_set.insert(bp[nxt].idx);
//					std::cerr << "interval " << bp[nxt].idx << " : "  << m_interval[bp[nxt].idx].pos << "-" << m_interval[bp[nxt].idx].end << " has been inserted at " << pos <<std::endl;; 
					nxt++;
					break;
				case 2:
					//interval end, post-gap start
					dp_set.erase(bp[nxt].idx); 
//					std::cerr << "interval " << bp[nxt].idx << " : "  << m_interval[bp[nxt].idx].pos << "-" << m_interval[bp[nxt].idx].end << " has been removed at " << pos << std::endl; 
					nxt++;
					break;
				case 3:
					//Post-gap end
					isz_set.erase(bp[nxt].idx);
					nxt++;
					break;
				default:
					std::cerr << "Something Wrong" << std::endl;
					exit(1);
					break;
			}
		}
		int m=0;
		for(int j=0;j<n_plp[0];++j)
		{
			const bam_pileup1_t *p = plp[0]+j;
			if (p->is_del || p->is_refskip )
				++m;
		}
		double dpval = n_plp[0]-m;
		double gc_dpval = dpval;

		gc_dpval = gcCorrected(dpval, chrnum, pos);

		for(set<int>::iterator it=dp_set.begin(); it!=dp_set.end(); ++it)
		{
//			std::cerr<< "dpval for pos " << pos << " is " << dpval << std::endl;
			dp_sum[*it] += dpval;
			gc_dp_sum[*it] += gc_dpval;
			dp_cnt[*it]++;
		}
	}
	sam_itr_destroy(data[0]->iter);
	free(plp); free(n_plp);
	bam_mplp_destroy(mplp);

	for(int i=0;i<n;++i)
	{
		std::string &txt = G[i];

		char buf[100];

//		if (dp_cnt[i] == 0)
//		{
//			std::cerr << " pos " << m_interval[i].pos << "-" << m_interval[i].end << " has count 0." << std::endl;;
//		}
		double dp = (dp_cnt[i]>0) ? dp_sum[i]/(double)dp_cnt[i] : 0;
		double gc_dp = (dp_cnt[i]>0) ? gc_dp_sum[i]/(double)dp_cnt[i] : 0;

		sprintf(buf,  "%.1f:%.1f:", dp, gc_dp);

		txt = buf;

//		std::cerr << m_interval[i].chr <<  ":" << m_interval[i].pos << "-" << m_interval[i].end << "\t" << m_interval[i].len() << "\t"<<  m_interval[i].svtype << "\t" << X[i] << "\t" << GX[i] <<"\t";

		//process_readpair(m_interval[i], isz_list[i], pos_list[i], txt);
		process_readpair(m_interval[i], isz_list[i], txt);

		txt += ":";

//		process_readpair(m_interval[i], rev_isz_list[i], rev_pos_list[i], txt);
		process_readpair(m_interval[i], rev_isz_list[i], txt);

		txt += ":";

		if (isz_cnt[i] >0)
		{
			txt += to_string(isz_cnt[i]) + "," + to_string((int)((double)isz_sum[i] / (double)isz_cnt[i])) ;
		}
		else
		{
			txt += "0,.";
		}
//		std::cerr << std::endl;
//		Y[i] = (i_cnt[i]>0) ? i_sum[i]/(double)i_cnt[i] : 0;
	}
}*/

/*
 void BamCram::get_avg_depth()
 {
 std::vector<double> sums (GC.num_bin, 0);
 std::vector<double> cnts (GC.num_bin, 0);
 
 // For average Insert Size
 std::vector< std::vector<int> > i_list;
 
 i_list.resize(1);
 
 data[0]->isz_list = &i_list;
 
 std::vector< std::vector <double> > GCdata;
 GCdata.resize(GC.num_bin);
 
 for(int i=0; i<(int)GC.regions.size();++i)
 {
 //    while(rpos<chrlen[c])
 {
 char reg[100];
 sprintf(reg, "chr%d:%d-%d",GC.regions[i].chrnum, GC.regions[i].pos, GC.regions[i].end);
 
 data[0]->iter = sam_itr_querys(idx, data[0]->hdr, reg);
 if (data[0]->iter == NULL)
 {
 std::cerr << reg << std::endl;
 std::cerr << "Can't parse region" << std::endl;
 exit(1);
 }
 
 //            bam_plp_t plp = bam_plp_init(read_bam_basic, (void*) data);
 bam_mplp_t mplp = bam_mplp_init(1, read_bam_basic, (void**)data);
 
 //const bam_pileup1_t *p;
 const bam_pileup1_t **plp;
 
 int tid, pos;
 //            int n_plp;
 int *n_plp = (int *) calloc (1, sizeof(int));
 plp = (const bam_pileup1_t**) calloc(1, sizeof(bam_pileup1_t*));
 
 //            while((p = bam_plp_auto(plp, &tid, &pos, &n_plp))!=0)
 while(bam_mplp_auto(mplp, &tid, &pos, n_plp, plp) > 0)
 {
 //            std::cerr << "tid " << tid << " pos " << pos << std::endl;
 if (pos<GC.regions[i].pos || pos >GC.regions[i].end) continue;
 sums[GC.regions[i].gcbin] += n_plp[0];
 cnts[GC.regions[i].gcbin] += 1;
 GCdata[GC.regions[i].gcbin].push_back(n_plp[0]);
 }
 free(n_plp); free(plp);
 sam_itr_destroy(data[0]->iter);
 bam_mplp_destroy(mplp);
 }
 
 }
 
 double avg = 0;
 std::vector<double> meds (GC.num_bin,0);
 
 for(int i=0; i<GC.num_bin; ++i)
 {
 if (cnts[i] > 0)
 {
 avg += (sums[i]/cnts[i]) * GC.gc_dist[i];
 //        std::cerr << "Avg DP for GC bin " << i << " is " << sums[i]/cnts[i] << ", median is " ;
 sort(GCdata[i].begin(), GCdata[i].end());
 int m = (int)(GCdata[i].size()/2) -1;
 meds[i] = (GCdata[i][m] + GCdata[i][m+1]) /2.0;
 //std::cerr << meds[i] << std::endl;
 }
 }
 
 std::cerr << "Average Depth: " << avg << std::endl;
 
 med_isize = median(i_list[0]);
 
 double sum_is = 0, sumsq_is = 0;
 int cnt_is = 0;
 for(int i=0;i<(int)i_list[0].size();++i)
 {
 sum_is += i_list[0][i];
 sumsq_is += i_list[0][i] * i_list[0][i];
 cnt_is += 1;
 }
 //    std::cerr << "sum " << sum_is << " sumsq " << sumsq_is << " cnt " << cnt_is << std::endl;
 
 if (cnt_is>0)
 {
 avg_isize = sum_is / cnt_is;
 std_isize = sqrt( (sumsq_is / (double)cnt_is) - (avg_isize * avg_isize) );
 }
 else
 {
 avg_isize = 0;
 std_isize = 0;
 }
 
 gc_factor.resize(GC.num_bin);
 
 std::cerr << "Median Insert Size : " << med_isize << std::endl;
 std::cerr << "Average Insert Size : " << avg_isize << " (+/- " << std_isize << ")" << std::endl;
 
 //    std::cerr << "gc factor size : " << gc_factor.size() << std::endl;
 
 for(int i=0; i<GC.num_bin; ++i)
 {
 //        std::cerr << "meds[" << i << "]=" << meds[i]<< " avg=" << avg << " factor " << meds[i]/avg <<std::endl;
 gc_factor[i] = meds[i]/avg;
 
 //        std::cerr << "bin " << i << " factor " << gc_factor[i] << std::endl;
 }
 
 // TEMPORARY - HARD CODING for 20 BINS
 gc_factor[18] = gc_factor[17];
 gc_factor[19] = gc_factor[17];
 
 avg_dp = avg;
 avg_rlen = 150; //TEMPORARY - FIX!
 }
 */


//void BamCram::process_readpair(sv &currsv, std::vector<int> &isz_list, std::vector<int> &pos_list, std::string &txt)
/*
 void BamCram::process_readpair(sv &currsv, std::vector<int> &isz_list, std::string &txt)
 {
 
 std::vector<int> P_isz;
 std::vector<int> N_isz;
 //    std::vector<int> P_pos;
 //    std::vector<int> N_pos;
 for(int i=0;i<(int)isz_list.size();++i)
 {
 if (isz_list[i]>0 && isz_list[i] < 10*currsv.len )
 {
 P_isz.push_back(isz_list[i]);
 //            P_pos.push_back(pos_list[i]);
 }
 else if (isz_list[i] > -10 * currsv.len)
 {
 N_isz.push_back(isz_list[i]);
 //            N_pos.push_back(pos_list[i]);
 }
 }
 
 if (P_isz.size() > 0)
 {
 //txt += to_string(P_isz.size()) + "," + to_string(median(P_isz)) + "," + to_string(median(P_pos)+1 - currsv.pos );
 txt += to_string(P_isz.size()) + "," + to_string(median(P_isz)) ;
 }
 else
 {
 txt += ".";
 }
 txt += ":";
 
 if (N_isz.size() > 0)
 {
 txt += to_string(N_isz.size()) + "," + to_string(median(N_isz));
 }
 else
 {
 txt += ".";
 }
 }
 */
