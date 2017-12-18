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
#include <math.h>
#include <algorithm>
#include <string>

// GRCh38
// TODO: Take FASTA and fai file input and match MD5 with BAMs? Or just use whatever given?

breakpoint::breakpoint()
{
	pos = 0;
	type = 0;
}

bool breakpoint::operator < (const breakpoint& b) const
{
	return (pos<b.pos);
}

bool breakpoint::operator == (const breakpoint& b) const
{
	return (pos == b.pos);
}

typedef enum { READ_UNKNOWN = 0, READ_1 = 1, READ_2 = 2 } readpart;

void readmagic(ifstream &F)
{
	char buf[100];
	F.read(reinterpret_cast<char *>(buf), 7);
	buf[7] = '\0';
	if (strcmp(buf, "mCNVMGC"))
	{
		cerr << "Error: GC content file is corrupted." << endl;
		exit(1);
	}
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
static int read_bam(void *data, bam1_t *b) // read level filters better go here to avoid pileup
{
	aux_t *aux = (aux_t*)data; // data in fact is a pointer to an auxiliary structure
	int ret;

	while (1)
	{
		ret = aux->iter ? sam_itr_next(aux->fp, aux->iter, b) : sam_read1(aux->fp, aux->hdr, b);
		
		// cerr << "insert size : " << b->core.isize << endl;

		if ( ret<0 ) break;
		if ( b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP) ) continue;
		if ( (int)b->core.qual < aux->min_mapQ ) continue;

		// Nov 29, 2017, commented out
		//if ( aux->min_len && bam_cigar2qlen(b->core.n_cigar, bam_get_cigar(b)) < aux->min_len ) continue;
//		cerr << bam_get_qname(b) << "/" << which_readpart(b) << " insert size : " << b->core.isize << endl;
		if (b->core.isize > 0)
		{
			for(set<int>::iterator it=(*(aux->fwd_set)).begin(); it!=(*(aux->fwd_set)).end(); ++it)
			{
				(*(aux->isz_list))[*it].push_back(b->core.isize);
//				(*(aux->isz_cnt))[*it]++;
			}
		}
		if(b->core.isize<0)
		{
			for(set<int>::iterator it=(*(aux->rev_set)).begin(); it!=(*(aux->rev_set)).end(); ++it)
			{

				(*(aux->isz_list))[*it].push_back(0-b->core.isize);
//				(*(aux->isz_sum))[*it] -= b->core.isize;
//				(*(aux->isz_cnt))[*it]++;
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

void bFile::get_avg_depth()
{
	// TODO : SAMPLE regions of whole genome
	// TODO : Get GC Content estimates and calculate GC-corrected average depth
	
	vector<double> sums (GC.num_bin, 0);
	vector<double> cnts (GC.num_bin, 0);

	// For average Insert Size
	set<int> f_set;
	set<int> r_set;
	vector< vector<int> > i_list;

	i_list.resize(1);

	f_set.insert(0);
	r_set.insert(0);

	data->fwd_set = &f_set;
	data->rev_set = &r_set;
	data->isz_list = &i_list;

	vector< vector <double> > GCdata;
	GCdata.resize(GC.num_bin);

	for(int i=0; i<GC.regions.size();++i)
	{
	//	while(rpos<chrlen[c])
		{
			char reg[100];
			sprintf(reg, "chr%d:%d-%d",GC.regions[i].chrnum, GC.regions[i].pos, GC.regions[i].end);

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
				if (pos<GC.regions[i].pos || pos >GC.regions[i].end) continue;
				sums[GC.regions[i].gcbin] += n_plp;
				cnts[GC.regions[i].gcbin] += 1;
				GCdata[GC.regions[i].gcbin].push_back(n_plp);
			}
			sam_itr_destroy(data->iter);
			bam_plp_destroy(plp);
		}
	}

	double avg = 0;
	vector<double> meds (GC.num_bin,0);

	for(int i=0; i<GC.num_bin; ++i)
	{
		if (cnts[i] > 0)
		{
			avg += (sums[i]/cnts[i]) * GC.gc_dist[i];
//			cerr << "Avg DP for GC bin " << i << " is " << sums[i]/cnts[i] << ", median is " ;
			sort(GCdata[i].begin(), GCdata[i].end());
			int m = (int)(GCdata[i].size()/2) -1;
			meds[i] = (GCdata[i][m] + GCdata[i][m+1]) /2.0;
//			cerr << meds[i] << endl;
		}
	}

	cerr << "Average Depth: " << avg << endl;

	sort(i_list[0].begin(), i_list[0].end());
	int m = (int)(i_list[0].size()/2) -1;
	avg_isize = (i_list[0][m] + i_list[0][m+1])/2.0; //TEMPORARY - FIX!

	cerr << "Median Insert Size : " << avg_isize << endl;
	gc_factor.resize(GC.num_bin);

	for(int i=0; i<GC.num_bin; ++i)
	{
		gc_factor[i] = meds[i] / avg;
		cerr << "bin " << i << " factor " << gc_factor[i] << endl;
	}

	avg_dp = avg;
	avg_rlen = 150; //TEMPORARY - FIX!
}


double bFile::gcCorrected(double D, int chr, int pos)
{
	int p = (int)round(pos / 200.0);
	int bin = GC.gc_array[chr][p];
//	if (bin>17) bin=17;

	if (bin<20 && gc_factor[bin]>0)
	{
//		cerr << D << " at " << chr << ":" << pos << " is adjusted to " << D/gc_factor[bin] << " by gc Factor "<< gc_factor[bin] << endl;
		return D / gc_factor[bin];
	}
	else
	{
		return D;
	}
}

// Get (overlapping) list of SV intervals, return average depth and GC-corrected average depth on intervals
void bFile::read_depth(vector<sv> &m_interval, vector<double> &X, vector<double> &GX, vector< vector<int> > &Y)
{
	char reg[100];
	
	string chr = m_interval[0].chr;
	int chrnum = m_interval[0].chrnum;
	int startpos = 1;
	int endpos = m_interval[0].end;
	int n = (int) m_interval.size();
	int gap = avg_isize-avg_rlen;
	vector<breakpoint> bp;
	bp.resize(n*4);  // 4 checkpoints per interval
	// Make list of all breakpoints, to identify the transition points for calculating 'average depth'
	// Breakpoints + read pair info (GAP = (AVG INSERT SIZE - AVG READ LEN) GAP bp before START, GAP bp after END point)
	
	// TODO : CHECK CIGAR to find SOFT CLIP with Breakpoint-Overlapping Reads? -- maybe not very complicated, but not sure how to incorporate into clustering model
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
			cerr << "Error! interval end point " << m_interval[i].end << " is before start point " << m_interval[i].pos << endl;
		}
		
		if (m_interval[i].end + gap > endpos)
		{
			endpos = m_interval[i].end + gap;
		}
	}
	
	sort(bp.begin()+1, bp.end());

	startpos = bp[0].pos;
	
	sprintf(reg, "chr%s:%d-%d", chr.c_str(), startpos, endpos);
//	sprintf(reg, "%s:%d-%d", chr.c_str(), startpos, endpos); // Check BAM/CRAM header for list of CHRs first?
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
	vector<double> gc_dp_sum (n,0);
	vector<int> dp_cnt (n,0);
	
	// For average Insert Size
	set<int> f_set;
	set<int> r_set;
	vector< vector<int> > i_list;
	i_list.resize(n);

//	vector<double> i_sum(n,0);
//	vector<int> i_cnt(n,0);
	f_set.insert(bp[0].idx);
	
	data->fwd_set = &f_set;
	data->rev_set = &r_set;
	data->isz_list = &i_list;
	/*
	data->isz_sum = &i_sum;
	data->isz_cnt = &i_cnt;
	*/
	
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
		int m=0;
		for(int j=0;j<n_plp;++j)
		{
			if ((p+j)->is_del || (p+j)->is_refskip )
				++m;
		}
		double dpval = n_plp-m;
		double gc_dpval = gcCorrected(n_plp-m, chrnum, pos);

		for(set<int>::iterator it=dp_set.begin(); it!=dp_set.end(); ++it)
		{
			dp_sum[*it] += dpval;
			gc_dp_sum[*it] += gc_dpval;
			dp_cnt[*it]++;
		}
	}
	sam_itr_destroy(data->iter);
	bam_plp_destroy(plp);

	for(int i=0;i<n;++i)
	{
		X[i] = (dp_cnt[i]>0) ? dp_sum[i]/(double)dp_cnt[i] : 0;
		GX[i] = (dp_cnt[i]>0) ? gc_dp_sum[i]/(double)dp_cnt[i] : 0;
	}
	for(int i=0;i<n;++i)
	{
		for(int j=0;j<i_list[i].size();++j)
		{

			if ((i_list[i][j] - avg_isize) < m_interval[i].len()/2)
			{
				Y[i][0]++;
			}
			else if ((i_list[i][j] - avg_isize) < m_interval[i].len()*1.5)
			{
				Y[i][1]++;
			}
			else
			{
				Y[i][2]++;
			}
		}
//		Y[i] = (i_cnt[i]>0) ? i_sum[i]/(double)i_cnt[i] : 0;
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

void gcContent::initialize(string &gcFile)
{// filename for GC content, populate all vectors
	ifstream inFile(gcFile.c_str(), std::ios::in | std::ios::binary);
	if (!inFile.good())
	{
		cerr << "Error: cannot open GC file."<< endl;
		exit(1);
	}
	
	readmagic(inFile);
	// read number of chrs
	inFile.read(reinterpret_cast <char *> (&num_chr), sizeof(uint8_t));
	
	cerr << (int)num_chr << " chromosomes in GC content file." << endl;

	// read size of chrs
	chrSize.resize(num_chr);
	for(int i=0;i<num_chr;++i)
	{
		inFile.read(reinterpret_cast <char *> (&chrSize[i]), sizeof(uint32_t));
//		cerr << "Chr " << i << " size: " << chrSize[i] <<  endl;
	}

	// read size of GC-interval bin
	inFile.read(reinterpret_cast <char *> (&binsize), sizeof(uint16_t));
	cerr << "Bin size: " << (int) binsize << endl;

	// read number of GC bins
	inFile.read(reinterpret_cast <char *> (&num_bin), sizeof(uint16_t));
	cerr << "Num_bin : " << num_bin << endl;

	// read number of total intervals
	inFile.read(reinterpret_cast <char *> (&total_bin), sizeof(uint16_t));
	cerr << "Total bin : " << total_bin << endl;

	regions.resize(total_bin);

	readmagic(inFile);

	gc_array.resize(num_chr);

	// read intervals (chr, start, end)
	for(int i=0; i<total_bin; ++i)
	{
		uint8_t c;
		uint32_t pos1, pos2;
		uint8_t gc;
		inFile.read(reinterpret_cast<char *>(&c), sizeof(uint8_t));
		inFile.read(reinterpret_cast<char *>(&pos1), sizeof(uint32_t));
		inFile.read(reinterpret_cast<char *>(&pos2), sizeof(uint32_t));
		inFile.read(reinterpret_cast<char *>(&gc), sizeof(uint8_t));

		regions[i].chrnum = c;
		regions[i].chr = "chr" + to_string(c);
		regions[i].pos = pos1;
		regions[i].end = pos2;
		regions[i].gcbin = gc;
	}

	readmagic(inFile);
	// read GC content for each (bin size)-bp interval for each chromosome
	// cerr << "Currnet position : " << inFile.tellg() << endl;
	for(int i=0; i<num_chr; ++i)
	{
		int N = ceil((chrSize[i] / (double)binsize)*2.0) ;
		gc_array[i] = (uint8_t *) calloc(N, sizeof(uint8_t));
		inFile.read(reinterpret_cast<char *>(gc_array[i]), sizeof(uint8_t)*N);
		readmagic(inFile);

		if (!inFile.good())
		{
			cerr << "Cannot finish reading GC content file." <<endl;
			exit(1);
		}
//		cerr << "Chr " << i << " GC content array loaded for "<<  N << " segments. " <<  endl;
	}
//	cerr << "Currnet position : " << inFile.tellg() << endl;

	for(int i=0;i<num_bin;++i)
	{
		if (!inFile.good())
		{
			cerr << "Cannot finish reading GC content file." <<endl;
			exit(1);
		}
		double v;
		inFile.read(reinterpret_cast<char *>(&v), sizeof(double));
//		cerr << "Bin " << i << " GC content proportion: "<< v << endl;
		gc_dist.push_back(v);
	}
	readmagic(inFile);
	inFile.close();
}

