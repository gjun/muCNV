//
//  bamcram.cpp
//  muCNV
//
//  Created by Goo Jun on 11/29/17.
//  Copyright © 2017 Goo Jun. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <set>
#include <unordered_set>
#include "muCNV.h"
#include <math.h>
#include <algorithm>
#include <string>
#include <queue>

#define IS_PROPERLYPAIRED(bam) (((bam)->core.flag&(BAM_FPAIRED|BAM_FPROPER_PAIR)) == (BAM_FPAIRED|BAM_FPROPER_PAIR) && !((bam)->core.flag&BAM_FUNMAP))

typedef enum { READ_UNKNOWN = 0, READ_1 = 1, READ_2 = 2 } readpart;

int get_cigar_clippos(string &cigar_str)
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
    
    if (rclip > lclip)
    {
        return -rclip;
    }
    else if (lclip>rclip)
    {
        return lclip;
    }
    return 0;
}

bool process_split(string &t, splitread &new_sp, int32_t tid, bool strand)
{
    int i=1;
    
    int32_t chr = 0;
    string cigar_str;
    string pos_str;

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
        
		/*
        string s;
        int c = i;
        
        //fprintf(stderr,"s = '%s', strlen(s) = %d, delims = '%s'\n",s,(int)strlen(s), delims);
        for(;t[i]!=',' && t[i]!=';' && t[i]!='\0';++i);

        if (t[i] ==';' || t[i] == '\0')
            return false;

        pos_str = std::string((char*)t+c, i-c);
        if (t[i+1] == '+')
        {
            if (!strand)
                return false;
        }
        else if (t[i+1] == '-')
        {
            if (strand)
                return false;
        }
        else
        {
            return false;
        }
        if (t[i+2] != ',')
            return false;
        i += 3;
        
        c=i;
        for(;t[i]!=',' && t[i]!=';' && t[i]!='\0';++i);
        if (t[i] ==';' || t[i] == '\0')
            return false;
		*/

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
		if (d == string::npos) return false;
		cigar_str = t.substr(c,d-c);
	//	DMSG("CIGAR: " << cigar_str);
		if (cigar_str.length() == 0 ) return false;
        // cigar_str = std::string((char*)t+c, i-c);
        // process pos and cigar here
//        new_sp.sapos = atoi(pos_str.c_str());
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

        if ((b->core.flag & BAM_FPAIRED ) && !(b->core.flag & BAM_FSUPPLEMENTARY) )
        {
            if (IS_PROPERLYPAIRED(b) && b->core.pos>0 && b->core.mtid == b->core.tid && b->core.isize!=0)
            {
				/*
				if (abs(b->core.isize) > 1000)
				{
					std::cerr << "isize " << b->core.isize << std::endl;
				}
				*/
                // get average isize statistics only from properly paired pairs
                aux->sum_isz += (b->core.isize > 0) ? b->core.isize : -b->core.isize;
                aux->sumsq_isz += (b->core.isize) * (b->core.isize);
                aux->n_isz += 1;

//				std::cerr << "n_isz : " << aux->n_isz << ", sum_isz : " << aux->sum_isz << ", sumsq_isz : " << aux->sumsq_isz << std::endl;
            }
            else 
            {
                if (b->core.qual>=10)
                {
                    // add to read pair set
                    readpair new_rp;
                    new_rp.chrnum = b->core.tid+1;
                    new_rp.selfpos = b->core.pos;
                    new_rp.pairstr = (b->core.flag & BAM_FREVERSE) ? 2 : 0;

                    if (b->core.tid == b->core.mtid)
                    {
                        uint8_t *aux_mq;
                        if ((aux_mq= bam_aux_get(b, "MQ")) != NULL)
                        {
                            new_rp.matequal = bam_aux2i(aux_mq);
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
                        new_rp.matequal =
                        new_rp.matepos = -1;
                        aux->n_rp++;
                        (*(aux->p_vec_rp)).push_back(new_rp);
                    }
                }
            }
        }
        
		if (b->core.flag & BAM_FSUPPLEMENTARY)
		{
			uint8_t *aux_sa = bam_aux_get(b, "SA");
			char* p_sa;
			if (aux_sa && (p_sa = bam_aux2Z(aux_sa)))
			{
				string str_sa = string(p_sa);

				splitread new_sp;
				new_sp.chrnum = b->core.tid + 1;
				new_sp.pos = b->core.pos;
				new_sp.firstclip = 0;
				new_sp.secondclip = 0;

				int ncigar = b->core.n_cigar;
				uint32_t *cigar  = bam_get_cigar(b);
				int16_t lclip = 0;
				int16_t rclip = 0;
				if (bam_cigar_op(cigar[0]) == BAM_CSOFT_CLIP)
				{
					lclip = bam_cigar_oplen(cigar[0]);
				}
				if (bam_cigar_op(cigar[ncigar-1]) == BAM_CSOFT_CLIP)
				{
					rclip = bam_cigar_oplen(cigar[0]);
				}
				if (rclip > lclip && rclip > 15) // arbitrary cutoff, 15
					new_sp.firstclip = -rclip;
				else if (lclip>rclip && lclip > 15)
					new_sp.firstclip = lclip;
				
				if (new_sp.firstclip != 0) // Process split read only if the read has soft-clipped ends
				{
					if (!process_split(str_sa, new_sp, b->core.tid, !(b->core.flag & BAM_FREVERSE)))
					{
						new_sp.sapos = 0;
						new_sp.secondclip = 0;
					}
					aux->n_sp++;
					(*(aux->p_vec_sp)).push_back(new_sp);
				}
			}
		}
        break;
    }
    return ret;
}


void bFile::initialize_sequential(string &bname)
{
    data = (aux_t**)calloc(1, sizeof(aux_t*));
    data[0] = (aux_t*)calloc(1, sizeof(aux_t));
    data[0]->fp = hts_open(bname.c_str(), "r");
    
    data[0]->sum_isz = 0;
    data[0]->sumsq_isz = 0;
    data[0]->n_isz = 0;
    
    //int rf = SAM_FLAG | SAM_RNAME | SAM_POS | SAM_MAPQ | SAM_CIGAR | SAM_SEQ | SAM_QUAL;
    int rf = SAM_FLAG |  SAM_MAPQ | SAM_QUAL| SAM_POS | SAM_SEQ | SAM_CIGAR| SAM_TLEN | SAM_RNEXT | SAM_PNEXT | SAM_AUX;
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

    idx = NULL;
}


void bFile::initialize(string &bname)
{

//	data = (aux_t*)calloc(1, sizeof(aux_t));
	data = (aux_t**)calloc(1, sizeof(aux_t*));
	data[0] = (aux_t*)calloc(1, sizeof(aux_t));
	data[0]->fp = hts_open(bname.c_str(), "r");
		
	//int rf = SAM_FLAG | SAM_RNAME | SAM_POS | SAM_MAPQ | SAM_CIGAR | SAM_SEQ | SAM_QUAL;
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
//	data->min_mapQ = 1; //filter out by MQ10
//	data->min_len = 0; // does not set minimum length
//	data->hdr = sam_hdr_read(data->fp);

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

double bFile::gcCorrected(double D, int chr, int pos)
{
    int p = pos*2 / GC.binsize;
//	std::cerr << "pos " << pos << " p " << p << std::endl;
	int bin = GC.gc_array[chr][p];

	if (bin<20 && gc_factor[bin]>0.0001)
	{
//		std::cerr << D << " at " << chr << ":" << pos << " is adjusted to " << D/gc_factor[bin] << " by gc Factor "<< gc_factor[bin] << std::endl;
		return D / gc_factor[bin];
	}
	else
	{
		// Let's not make adjustment for bins with '255' value
		return D;
	}
}


// Get (overlapping) list of SV intervals, return average depth and GC-corrected average depth on intervals
void bFile::read_depth_sequential(std::vector<breakpoint> &vec_bp, std::vector<sv> &vec_sv)
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
    
    gc_factor.resize(GC.num_bin);
    gc_sum.resize(GC.num_bin);
    gc_cnt.resize(GC.num_bin);
    for(int i=0;i<GC.num_bin;++i)
    {
        gc_sum[i] = 0;
        gc_cnt[i] = 0;
    }
    std::vector<int> dp_list;
    
    std::vector<double> dp_sum; // One interval (start-end) will add one entry on these
    std::vector<double> gc_dp_sum;
    std::vector<int> dp_cnt;

    data[0]->p_vec_rp = &vec_rp;
    data[0]->p_vec_sp = &vec_sp;
    
    data[0]->sum_isz = 0;
    data[0]->sumsq_isz = 0;
    data[0]->n_isz = 0;

	data[0]->n_rp = 0;
	data[0]->n_sp = 0;
    
    int* n_plp = (int *) calloc(1, sizeof(int));;
    const bam_pileup1_t **plp = (const bam_pileup1_t **) calloc(1, sizeof(bam_pileup1_t*));

    bam_mplp_t mplp = bam_mplp_init(1, read_bam, (void**) data);
	bam_mplp_set_maxcnt(mplp, 1000); // This is arbitrary, but we don't need > 1000 depth

    
    depth100.resize(GC.num_chr + 1);
	nbin_100.resize(GC.num_chr + 1);
	nbin_100[0] = 0;
    
    for(int i=1; i<=GC.num_chr; ++i)
    {
        nbin_100[i] = ceil((double)GC.chrSize[i] / 100.0) + 1 ;
        depth100[i] = (uint16_t *) calloc(nbin_100[i], sizeof(uint16_t));
    }
    
    int sum100=0;
    int n100=0;
    
    int prev_chrnum = 1;
    int prev_pos = 1;
    

	std::cerr << "processing chr 1" << std::endl;
    
    while(bam_mplp_auto(mplp, &tid, &pos, n_plp, plp)>0 && tid < GC.num_chr)
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
                depth100[prev_chrnum][prev_pos/100] = (uint16_t) val;
            }
			std::cerr << "n_isz : " << data[0]->n_isz << ", sum_isz : " << data[0]->sum_isz << ", sumsq_isz : " << data[0]->sumsq_isz << std::endl;

            sum100=0;
            n100=0;
            prev_chrnum = chrnum;
        }
        else if (pos/100 > prev_pos/100)
        {
            // Update 100-bp depth
            if (n100>0)
            {
                int val = round((double) (sum100 * 32) / n100);
                if (val>65535) val=65535; // handle overflow
                depth100[chrnum][pos/100] = (uint16_t) val;
            }

            sum100=0;
            n100=0;
        }
        
        prev_pos = pos; // maybe use array idx instead of pos ?
        breakpoint curr_bp;
        
        curr_bp.chrnum = chrnum;
        curr_bp.pos = pos;
        
        if (nxt>=vec_bp.size())
            break;
        
        while (vec_bp[nxt] <= curr_bp)
        {
            if (vec_bp[nxt].bptype == 0)
            {
//				std::cerr<<"adding bp idx " << vec_bp[nxt].idx << " pos " <<  vec_bp[nxt].pos;
            	dp_list.push_back(vec_bp[nxt].idx);
//				std::cerr << " dp_list [";
//				for(int ii=0;ii<dp_list.size();++ii)
//					std::cerr << dp_list[ii] << " " ;
//				std::cerr << "]" << std::endl;

            }
            else
            {                
                std::vector<int>::iterator it = find(dp_list.begin(), dp_list.end(), vec_bp[nxt].idx) ;
				if (it == dp_list.end())
				{
					int idx_nxt = vec_bp[nxt].idx;
					std::cerr << "Error, bptype " << vec_bp[nxt].bptype << " bp pos " << vec_bp[nxt].pos << ", sv " << vec_sv[idx_nxt].chrnum << ":" << vec_sv[idx_nxt].pos << "-" << vec_sv[idx_nxt].end << std::endl;

				}
				else
				{
//					dp_list.erase( it, it+1 );
			//		std::cerr<<"removing bp idx " << vec_bp[nxt].idx << " pos " <<  vec_bp[nxt].pos;
					*it=dp_list.back();
					dp_list.pop_back();
			//		std::cerr << " dp_list [";
			//		for(int ii=0;ii<dp_list.size();++ii)
			//			std::cerr << dp_list[ii] << " " ;
			//		std::cerr << "]" << std::endl;
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

        // double gc_dpval = dpval;
        // gc_dpval = gcCorrected(dpval, chrnum, pos);
        
        //sum100 += gc_dpval;
        sum100 += dpval;
        n100 ++;
        // for whole-genome average
        sum_dp += dpval;
        sumsq_dp += dpval * dpval;
        n_dp += 1;
        
        int bin = GC.gc_array[tid+1][(int)pos*2 / GC.binsize];
        gc_sum[bin] += dpval;
        gc_cnt[bin] += 1;
        
		for(std::vector<int>::iterator it=dp_list.begin(); it!=dp_list.end(); ++it)
		{
			if (*it <0 || *it>vec_sv.size())
			{
				std::cerr << "Error, " << *it << " is out of bound " << std::endl;
			}
			vec_sv[*it].dp_sum += dpval;
			vec_sv[*it].n_dp += 1;
		}
    }
    avg_dp = (double) sum_dp / n_dp;
    std_dp = sqrt(((double)sumsq_dp / n_dp - (avg_dp*avg_dp)));
    
    // Calculate GC-curve and save it with original DP to maximize information preservation instead of storing GC-corrected depths only
    for(int i=0;i<GC.num_bin;++i)
    {
        if (gc_cnt[i]>20)
        {
            gc_factor[i] = avg_dp / ((double)gc_sum[i] / gc_cnt[i]) ; // multiplication factor
        }
        else
        {
            gc_factor[i] = 0;
        }
    }

	std::cerr << "n_rp : " << data[0]->n_rp << " n_sp: " << data[0]->n_sp << std::endl;
    avg_isize = data[0]->sum_isz / (double)data[0]->n_isz;
    std_isize = sqrt((double)data[0]->sumsq_isz /(double)data[0]->n_isz - ((double)avg_isize*avg_isize));

	std::cerr << "n_isz : " << data[0]->n_isz << ", sum_isz : " << data[0]->sum_isz << ", sumsq_isz : " << data[0]->sumsq_isz << std::endl;
	std::cerr << "avg_isz : " << avg_isize << ", std_isize : " << std_isize << std::endl;

    sam_itr_destroy(data[0]->iter);
    free(plp); free(n_plp);
    bam_mplp_destroy(mplp);
}


void bFile::postprocess_depth(std::vector<sv> &vec_sv)
{

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

void bFile::write_pileup(string &sampID, std::vector<sv> &vec_sv)
{
    string pileup_name = sampID + ".pileup";
    string varfile_name = sampID + ".var";
    string idxfile_name = sampID + ".idx";

    size_t curr_pos = 0;
    
    std::ofstream pileupFile(pileup_name.c_str(), std::ios::out | std::ios::binary);
    std::ofstream idxFile(idxfile_name.c_str(), std::ios::out | std::ios::binary);

    int n_sample = 1;
    char pad[256] = {0};

    pileupFile.write(reinterpret_cast<char*>(&n_sample), sizeof(int));
    curr_pos += sizeof(int);

    // Sample ID (each with 256 bytes)
    if (sampID.length() > 255)
    {
        std::cerr << "Error, sample ID " << sampID << " is too long." << std::endl;
        exit(1);
    }
    pileupFile.write(sampID.c_str(), sampID.length());
    pileupFile.write(pad, 256-sampID.length());
	curr_pos += 256;

    // Write depth and isize stats
    pileupFile.write(reinterpret_cast<char*>(&avg_dp), sizeof(double));
    pileupFile.write(reinterpret_cast<char*>(&std_dp), sizeof(double));
    pileupFile.write(reinterpret_cast<char*>(&avg_isize), sizeof(double));
    pileupFile.write(reinterpret_cast<char*>(&std_isize), sizeof(double));

	curr_pos += sizeof(double) * 4;
    
    // Write GC curve
    for(int i=0;i<GC.num_bin;++i)
    {
        pileupFile.write(reinterpret_cast<char*>(&(gc_factor[i])), sizeof(double));
		curr_pos += sizeof(double);
    }

    // Write Index of var files (every chr offset, 1000-th variants)
//    std::cerr << "Sample " << sampID << ", header length " << curr_pos << std::endl;
    
    idxFile.write(reinterpret_cast<char*>(&curr_pos), sizeof(size_t)); // where SV DP starts

    // Write DP100
    for(int i=1; i<=GC.num_chr; ++i)
    {
        pileupFile.write(reinterpret_cast<char*>(depth100[i]), sizeof(uint16_t)*(nbin_100[i]));
		curr_pos += sizeof(uint16_t)*nbin_100[i];
    }

 //   std::cerr << "After DP100 written, curr_pos is at " << curr_pos << std::endl;

    int sp_idx = 0;
    int rp_idx = 0;
    int prev_sp = 0;
    int prev_rp = 0;
    
    int cnt_rp = 0;
    int cnt_sp = 0;
 
    for(int i=1;i<=GC.num_chr; ++i)
    {
        int N = ceil((double)GC.chrSize[i] / 10000.0) ;

        for(int j=1;j<=N;++j)
        {
			idxFile.write(reinterpret_cast<char*>(&curr_pos), sizeof(size_t)); // where each 10,000-bp interval starts;
            // RP
            while(rp_idx < (int)vec_rp.size() && vec_rp[rp_idx].chrnum == i && vec_rp[rp_idx].selfpos <= j*10000)
            {
                rp_idx ++;
            }
            uint16_t n_rp = (uint16_t) rp_idx - prev_rp;
            pileupFile.write(reinterpret_cast<char*>(&n_rp), sizeof(uint16_t));
            curr_pos += sizeof(uint16_t);
            
            for(int k=prev_rp; k<rp_idx; ++k)
            {
                pileupFile.write(reinterpret_cast<char*>(&(vec_rp[k].chrnum)), sizeof(int8_t));
                pileupFile.write(reinterpret_cast<char*>(&(vec_rp[k].selfpos)), sizeof(int32_t));
                pileupFile.write(reinterpret_cast<char*>(&(vec_rp[k].matepos)), sizeof(int32_t));
                pileupFile.write(reinterpret_cast<char*>(&(vec_rp[k].matequal)), sizeof(uint8_t));
                pileupFile.write(reinterpret_cast<char*>(&(vec_rp[k].pairstr)), sizeof(int8_t));
                curr_pos += sizeof(int8_t) + sizeof(uint32_t) + sizeof(uint32_t) + sizeof(uint8_t) + sizeof(int8_t);
                cnt_rp ++;
            }
            prev_rp = rp_idx;
            // SP
            while(sp_idx < (int)vec_sp.size() && vec_sp[sp_idx].chrnum == i && vec_sp[sp_idx].pos <= j*10000)
            {
                sp_idx ++;
            }
            uint16_t n_sp = (uint16_t) sp_idx - prev_sp;
            pileupFile.write(reinterpret_cast<char*>(&n_sp), sizeof(uint16_t));
            curr_pos += sizeof(uint16_t);
            
            for(int k=prev_sp; k<sp_idx; ++k)
            {
                pileupFile.write(reinterpret_cast<char*>(&(vec_sp[k].chrnum)), sizeof(int8_t));
                pileupFile.write(reinterpret_cast<char*>(&(vec_sp[k].pos)), sizeof(int32_t));
                pileupFile.write(reinterpret_cast<char*>(&(vec_sp[k].sapos)), sizeof(int32_t));
                pileupFile.write(reinterpret_cast<char*>(&(vec_sp[k].firstclip)), sizeof(int16_t));
                pileupFile.write(reinterpret_cast<char*>(&(vec_sp[k].secondclip)), sizeof(int16_t));
                curr_pos += sizeof(int8_t) + sizeof(uint32_t) + sizeof(uint32_t) + sizeof(int16_t) + sizeof(int16_t);
                cnt_sp ++;

            }
            prev_sp = sp_idx;
        }
       // fprintf(stderr, "\rChr %d, Index %d, cnt_rp %d, cnt_sp %d\n", i, (int)curr_pos, cnt_rp, cnt_sp);
    }
    pileupFile.close();
    idxFile.close();

    std::ofstream varFile(varfile_name.c_str(), std::ios::out | std::ios::binary);

    // Number of samples in this varFile (can include multiple samples)
    varFile.write(reinterpret_cast<char*>(&n_sample), sizeof(int));
    
    // Number of SVs in this varFile (can include multiple samples)
    int n_var = (int)vec_sv.size();
    varFile.write(reinterpret_cast<char*>(&n_var), sizeof(int));

    // Sample ID (each with 256 bytes)
    if (sampID.length() > 255)
    {
        std::cerr << "Error, sample ID " << sampID << " is too long." << std::endl;
        exit(1);
    }
    varFile.write(sampID.c_str(), sampID.length());
    varFile.write(pad, 256-sampID.length());
    
    for(int i=0;i<(int)vec_sv.size();++i)
    {
        varFile.write(reinterpret_cast<char*>(&(vec_sv[i].dp)), sizeof(uint16_t));
    }
    varFile.close();
}



void bFile::write_pileup_text(string &sampID, std::vector<sv> &vec_sv)
{
    // TODO: add write GC curve
    // TODO: write average 100-bp depth
    string fname = sampID + ".pileup.txt";
    FILE *fp = fopen(fname.c_str(), "w");
    
    fclose(fp);
    
    fname = sampID + ".var.txt";
    fp = fopen(fname.c_str(), "w");
    for(int i=0;i<(int)vec_sv.size();++i)
    {
        /*
        uint16_t n_rp = (uint16_t) vec_sv[i].vec_pair.size();
        uint16_t n_sp = (uint16_t) vec_sv[i].vec_split.size();
        */
        fprintf(fp, "%d", vec_sv[i].dp);
        /*
        fprintf(fp, "\t%d\t%d\t", n_rp, n_sp);

        for(int j=0;j<vec_sv[i].vec_pair.size();++j)
        {
            fprintf(fp, "%d,%d,%d,%d;", vec_sv[i].vec_pair[j].selfpos, vec_sv[i].vec_pair[j].matepos, vec_sv[i].vec_pair[j].selfstr, vec_sv[i].vec_pair[j].matestr);
        }
        fprintf(fp,"\t");
        for(int j=0;j<vec_sv[i].vec_split.size();++j)
        {
            fprintf(fp, "%d,%d,%d,%d;", vec_sv[i].vec_split[j].pos, vec_sv[i].vec_split[j].sapos, vec_sv[i].vec_split[j].firstclip, vec_sv[i].vec_split[j].secondclip);
        }
         */
        fprintf(fp, "\n");
         
    }
    fclose(fp);
}

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
void bFile::read_depth(std::vector<sv> &m_interval, std::vector<string> &G )
{
	char reg[100];
	
	string chr = m_interval[0].chr;
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
		string &txt = G[i];

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
 void bFile::get_avg_depth()
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


//void bFile::process_readpair(sv &currsv, std::vector<int> &isz_list, std::vector<int> &pos_list, string &txt)
/*
 void bFile::process_readpair(sv &currsv, std::vector<int> &isz_list, string &txt)
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
