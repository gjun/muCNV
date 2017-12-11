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
#include "muCNV.h"
#include <algorithm>

// GRCh38
// TODO: Take FASTA and fai file input and match MD5 with BAMs? Or just use whatever given?

breakpoint::breakpoint()
{
	pos = 0;
	type = 0;
}

static int readcount;

bool breakpoint::operator < (const breakpoint& b) const
{
	return (pos<b.pos);
}

bool breakpoint::operator == (const breakpoint& b) const
{
	return (pos == b.pos);
}

typedef enum { READ_UNKNOWN = 0, READ_1 = 1, READ_2 = 2 } readpart;

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

// This function reads a BAM alignment from one BAM file.
static int read_bam(void *data, bam1_t *b) // read level filters better go here to avoid pileup
{
	aux_t *aux = (aux_t*)data; // data in fact is a pointer to an auxiliary structure
	int ret;

	while (1)
	{
		ret = aux->iter ? sam_itr_next(aux->fp, aux->iter, b) : sam_read1(aux->fp, aux->hdr, b);
		
		// readcount++;
		// cerr << "insert size : " << b->core.isize << endl;
	

		if ( ret<0 ) break;
		if ( b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP) ) continue;
		if ( (int)b->core.qual < aux->min_mapQ ) continue;

		// Nov 29, 2017, commented out
		//if ( aux->min_len && bam_cigar2qlen(b->core.n_cigar, bam_get_cigar(b)) < aux->min_len ) continue;
//		readcount++;
//		cerr << bam_get_qname(b) << "/" << which_readpart(b) << " insert size : " << b->core.isize << endl;
		if (b->core.isize > 0)
		{
			for(set<int>::iterator it=(*(aux->fwd_set)).begin(); it!=(*(aux->fwd_set)).end(); ++it)
			{
				(*(aux->isz_sum))[*it] += b->core.isize;
				(*(aux->isz_cnt))[*it]++;
			}
		}
		if(b->core.isize<0)
		{
			for(set<int>::iterator it=(*(aux->rev_set)).begin(); it!=(*(aux->rev_set)).end(); ++it)
			{
				(*(aux->isz_sum))[*it] -= b->core.isize;
				(*(aux->isz_cnt))[*it]++;
			}
		}
		break;
	}

	
	return ret;
}

void bFile::initialize(string &bname)
{
	data = (aux_t*)calloc(1, sizeof(aux_t));
	data->fp = hts_open(bname.c_str(), "r");
		
	//int rf = SAM_FLAG | SAM_RNAME | SAM_POS | SAM_MAPQ | SAM_CIGAR | SAM_SEQ | SAM_QUAL;
	// Nov30, 2017, removed CIGAR flag
	int rf = SAM_FLAG | SAM_QNAME | SAM_AUX | SAM_POS | SAM_MAPQ | SAM_SEQ | SAM_QUAL | SAM_TLEN ;
	if (hts_set_opt(data->fp, CRAM_OPT_REQUIRED_FIELDS, rf))
	{
		cerr << "Failed to set CRAM_OPT_REQUIRED_FIELDS value" << endl;
		exit(1);
	}
	if (hts_set_opt(data->fp, CRAM_OPT_DECODE_MD, 0))
	{
		fprintf(stderr, "Failed to set CRAM_OPT_DECODE_MD value\n");
		exit(1);
	}
	data->min_mapQ = 1; //filter out only MQ0
	data->min_len = 0; // does not set minimum length
	data->hdr = sam_hdr_read(data->fp);
	if (data->hdr == NULL)
	{
		cerr << "Cannot open CRAM/BAM header" << endl;
		exit(1);
	}
	
	hts_idx_t* tmp_idx = sam_index_load(data->fp, bname.c_str());
	if (tmp_idx == NULL)
	{
		cerr << "Cannot open CRAM/BAM index" << endl;
		exit(1);
	}
	idx = tmp_idx;
}

void bFile::get_avg_depth(double &avg, double &GCavg, double &avg_insert)
{
	// TODO : SAMPLE regions of whole genome
	// TODO : Get GC Content estimates and calculate GC-corrected average depth
	
	double sum = 0;
	double cnt = 0;
	for(int c=20;c<21;++c)
	{
		int rpos=100000;
	//	while(rpos<chrlen[c])
		{
			char reg[100];
			sprintf(reg, "chr%d:%d-%d",c, rpos, rpos+100);
			//			cerr << reg << endl;
			

			data->iter = sam_itr_querys(idx, data->hdr, reg);
			if (data->iter == NULL)
			{
				cerr << reg << endl;
				cerr << "Can't parse region" << endl;
				exit(1);
			}
			bam_plp_t plp = bam_plp_init(read_bam, (void*) data);
			const bam_pileup1_t *p;
			
			int tid, pos;
			int n_plp;
			
			while((p = bam_plp_auto(plp, &tid, &pos, &n_plp))!=0)
			{
				//			cerr << "tid " << tid << " pos " << pos << endl;
				if (pos<rpos || pos >=rpos+10) continue;
				sum += n_plp;
				cnt ++;
			}
			delete p;
			bam_plp_destroy(plp);
			rpos+=1000000;
		}
	}

	avg = sum/cnt;
	GCavg = avg;
	avg_insert = 500; //TEMPORARY - FIX!
}


// Get (overlapping) list of SV intervals, return average depth and GC-corrected average depth on intervals
void bFile::read_depth(vector<sv> &m_interval, vector<double> &X, vector<double> &GX, vector<double> &Y)
{
	char reg[100];
	
	string chr = m_interval[0].chr;
	int startpos = 1;
	int endpos = m_interval[0].end;
	int n = (int) m_interval.size();
	vector<breakpoint> bp;
	bp.resize(n*4);  // 4 checkpoints per interval
	// Make list of all breakpoints, to identify the transition points for calculating 'average depth'
	// Breakpoints + read pair info (GAP = (AVG INSERT SIZE - AVG READ LEN) GAP bp before START, GAP bp after END point)
	
	// TODO : CHECK CIGAR to find SOFT CLIP with Breakpoint-Overlapping Reads? -- maybe not very complicated, but not sure how to incorporate into clustering model
	for(int i=0; i<n; ++i)
	{
		// START - GAP
		bp[i*4].pos = m_interval[i].pos > 300 ? m_interval[i].pos - 300 : 1;
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
		bp[i*4+3].pos = m_interval[i].end + 300;
		bp[i*4+3].type = 3;
		bp[i*4+3].idx = i;
		
		if (m_interval[i].end <= m_interval[i].pos )
		{
			cerr << "Error! interval end point " << m_interval[i].end << " is before start point " << m_interval[i].pos << endl;
		}
		
		if (m_interval[i].end + 300 > endpos)
		{
			endpos = m_interval[i].end + 300;
		}
	}
	
	startpos = bp[0].pos;
	
	sort(bp.begin()+1, bp.end());
	
	//sprintf(reg, "chr%s:%d-%d", chr.c_str(), startpos, endpos);
	sprintf(reg, "%s:%d-%d", chr.c_str(), startpos, endpos); // Check BAM/CRAM header for list of CHRs first?
	data->iter = sam_itr_querys(idx, data->hdr, reg);
	if (data->iter == NULL)
	{
		cerr << reg << endl;
		cerr << "Can't parse region " << reg << endl;
		exit(1);
	}

	int tid, pos;
	int n_plp;
	
	if (bp[0].idx != 0)
	{
		cerr << "Error: Merged interval's first element is not the earliest." << endl;
		cerr << "first " << bp[0].idx  << " pos " << bp[0].pos << " second " << bp[1].idx << " pos " << bp[1].pos << endl;
	}
	if (bp[0].type > 1)
	{
		cerr << "Error: Merged interval's earliest position is SV-end, not SV-start." << endl;
	}
	int nxt = 1;
	
	// For average DP
	set<int> dp_set;
	vector<double> dp_sum (n,0);
	vector<int> dp_cnt (n,0);
	
	// For average Insert Size
	set<int> f_set;
	set<int> r_set;
	vector<double> i_sum(n,0);
	vector<int> i_cnt(n,0);
	f_set.insert(bp[0].idx);
	
	data->fwd_set = &f_set;
	data->rev_set = &r_set;
	data->isz_sum = &i_sum;
	data->isz_cnt = &i_cnt;
	
	bam_plp_t plp = bam_plp_init(read_bam, (void*) data);
	const bam_pileup1_t *p;
	
	while((p = bam_plp_auto(plp, &tid, &pos, &n_plp))!=0)
	{
		if (pos<startpos || pos >=endpos) continue;
		
		if (pos>=bp[nxt].pos)
		{
			switch(bp[nxt].type)
			{
				case 0:
					//Pre-gap start
					f_set.insert(bp[nxt].idx);
					nxt++;
					break;
				case 1:
					//Pre-gap end, interval start
					f_set.erase(bp[nxt].idx);  // TODO, error checking : should check whether it exists in the set
					dp_set.insert(bp[nxt].idx);
					nxt++;
					break;
				case 2:
					//interval end, post-gap start
					r_set.insert(bp[nxt].idx);
					dp_set.erase(bp[nxt].idx); 
					nxt++;
					break;
				case 3:
					//Post-gap end
					r_set.erase(bp[nxt].idx);
					break;
			}
		}
		for(set<int>::iterator it=dp_set.begin(); it!=dp_set.end(); ++it)
		{
			int m=0;
		
			for(int j=0;j<n_plp;++j)
			{
				if ((p+j)->is_del || (p+j)->is_refskip )
					++m;
			}
			dp_sum[*it] += n_plp-m;
			dp_cnt[*it]++;
		}
	}
	sam_itr_destroy(data->iter);
	bam_plp_destroy(plp);

	for(int i=0;i<n;++i)
	{
		X[i] = (dp_cnt[i]>0) ? dp_sum[i]/(double)dp_cnt[i] : 0;
		GX[i] = X[i]; // TODO: do GC correction
	}
	for(int i=0;i<n;++i)
	{
		Y[i] = (i_cnt[i]>0) ? i_sum[i]/(double)i_cnt[i] : 0;
	}
}

double bFile::read_pair(sv &interval)
{
	char reg[100];
	
	int isize = 400; // TEMPORARY - get value from average insert size & read lengh?
	int readlen = 150;
	int margin = isize - readlen/2;
	int sum = 0;
	int cnt = 0;
	
	string chr = interval.chr;
	int start1 = interval.pos > margin ? interval.pos - margin : 1;
	int end1 = interval.pos > readlen ? interval.pos - readlen : 1;

	int start2 = interval.end;
	int end2 = interval.end + margin;
	
	sprintf(reg, "chr%s:%d-%d", chr.c_str(), start1, end1);
	
	data->iter = sam_itr_querys(idx, data->hdr, reg);
	bam1_t *b;
	b=bam_init1();
	
	if (data->iter == NULL)
	{
		cerr << reg << endl;
		cerr << "Can't parse region " << reg << endl;
		exit(1);
	}
	int result;
	while((result = sam_itr_next(data->fp , data->iter, b)) >= 0)
	{
		if (b->core.flag & (BAM_FDUP)) continue;
		if (b->core.isize > 0)
		{
			sum += b->core.isize;
			cnt ++ ;
		}

	}
	hts_itr_destroy(data->iter);
	
	//sprintf(reg, "chr%s:%d-%d", chr.c_str(), start2, end2);
	sprintf(reg, "%s:%d-%d", chr.c_str(), start2, end2);

	data->iter = sam_itr_querys(idx, data->hdr, reg);
	if (data->iter == NULL)
	{
		cerr << reg << endl;
		cerr << "Can't parse region " << reg << endl;
		exit(1);
	}
	while((result = sam_itr_next(data->fp , data->iter, b)) >= 0)
	{
		if (b->core.flag & (BAM_FDUP)) continue;
		if (b->core.isize < 0)
		{
			sum -= b->core.isize;
			cnt ++;
		}
	}
	bam_destroy1(b);
	double ret = (cnt>0) ? (double)sum/(double)cnt : 0;
	return ret;
}

int gcContent::gcbin(int chr, int pos)
{
	return 0;
}

double gcContent::gcFactor(int, int)
{
	// chr, pos, return factor for GC correction
	return 1;
}

void gcContent::initialize(string &gcFile)
{// filename for GC content, populate all vectors
	ifstream inFile(gcFile.c_str(), std::ios::in | std::ios::binary);
	
	// read number of chrs
	inFile.read(reinterpret_cast <char *> (&num_chr), sizeof(uint16_t));
	
	// read size of chrs
	chrSize.resize(num_chr);
	for(uint16_t i=0;i<num_chr;++i)
	{
		inFile.read(reinterpret_cast <char *> (&chrSize[i]), sizeof(uint16_t));
	}
	// read size of GC bins
	inFile.read(reinterpret_cast <char *> (&binsize), sizeof(uint16_t));

	// read number of GC bins
	inFile.read(reinterpret_cast <char *> (&num_bin), sizeof(uint16_t));

	// read number of intervals per GC bin
	inFile.read(reinterpret_cast <char *> (&int_per_bin), sizeof(uint16_t));

	// read intervals (chr, start, end)
	
	// read GC content for each (bin size)-bp interval for each chromosome
	

//	vector< vector<sv> > regions; // Double array to store list of regions for each GC bin -- non-overlapping, so let's just be it out-of-order
//	vector< vector<int> > gc_array; // Array to store GC content for every 400-bp (?) interval
//	vector<double> gc_dist; // Array to store proportion of Ref genome for each GC content bin

}

