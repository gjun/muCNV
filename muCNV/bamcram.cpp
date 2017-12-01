//
//  bamcram.cpp
//  muCNV
//
//  Created by Goo Jun on 11/29/17.
//  Copyright © 2017 Goo Jun. All rights reserved.
//

#include <stdio.h>
#include <set>
#include "muCNV.h"

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

int chrlen[26] = {0, 248956422,242193529,198295559,190214555,181538259,170805979,159345973,145138636,138394717,133797422,135086622,133275309,114364328,107043718,101991189,90338345,83257441,80373285,58617616,64444167,46709983,50818468,156040895,57227415, 16569};

// This function reads a BAM alignment from one BAM file.
static int read_bam(void *data, bam1_t *b) // read level filters better go here to avoid pileup
{
	aux_t *aux = (aux_t*)data; // data in fact is a pointer to an auxiliary structure
	int ret;
	while (1)
	{
		ret = aux->iter? sam_itr_next(aux->fp, aux->iter, b) : sam_read1(aux->fp, aux->hdr, b);
		if ( ret<0 ) break;
		if ( b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP) ) continue;
		if ( (int)b->core.qual < aux->min_mapQ ) continue;
		
		// Nov 29, 2017, commented out
		//if ( aux->min_len && bam_cigar2qlen(b->core.n_cigar, bam_get_cigar(b)) < aux->min_len ) continue;
		break;
	}
	return ret;
}

void bFile::initialize(string &bname)
{
//	n = (int)bnames.size();
	
//	data = (aux_t**)calloc(n, sizeof(aux_t*));
	data = (aux_t*)calloc(1, sizeof(aux_t));
	

	data = (aux_t*)calloc(1, sizeof(aux_t));
	data->fp = hts_open(bname.c_str(), "r");
		
	//int rf = SAM_FLAG | SAM_RNAME | SAM_POS | SAM_MAPQ | SAM_CIGAR | SAM_SEQ | SAM_QUAL;
	// Nov30, 2017, removed CIGAR flag
	int rf = SAM_FLAG | SAM_RNAME | SAM_POS | SAM_MAPQ | SAM_SEQ | SAM_QUAL;
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
		while(rpos<chrlen[c])
		{
			char reg[100];
			// TEMPORARY!! Use correct CHR name
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
void bFile::read_depth(vector<sv> &m_interval, vector<double> &X, vector<double> &GX)
{
	char reg[100];
	//	sprintf(reg, "%d:%d-%d",interval.chr, interval.pos, interval.end);
	
	// TEMPORARY!! Use correct CHR name
	// TODO: Change SV.chr to string
	
	int chr = m_interval[0].chr;
	int startpos = m_interval[0].pos;
	int endpos = m_interval[0].end;
	int n = (int) m_interval.size();
	vector<breakpoint> bp;
	bp.resize(n*2);
	
	for(int i=0; i<n; ++i)
	{
		bp[i*2].pos = m_interval[i].pos;
		bp[i*2].type = false;
		bp[i*2].idx = i;
		bp[i*2+1].pos = m_interval[i].end;
		bp[i*2+1].type = true;
		bp[i*2+1].idx = i;
		
		if (m_interval[i].end > endpos)
		{
			endpos = m_interval[i].end;
		}
	}
	
	sort(bp.begin(), bp.end());
	// make sure there's no SVs witn pos == end
	
	sprintf(reg, "chr%d:%d-%d", chr, startpos, endpos);

	data->iter = sam_itr_querys(idx, data->hdr, reg);
	if (data->iter == NULL)
	{
		cerr << reg << endl;
		cerr << "Can't parse region" << reg << endl;
		exit(1);
		
	}

	bam_plp_t plp = bam_plp_init(read_bam, (void*) data);
	const bam_pileup1_t *p;
	
	int tid, pos;
	int n_plp;
	vector<double> sum (n,0);
	vector<int> cnt (n,0);
	set<int> currset;
	currset.insert(bp[0].idx); // Should be always 0 because m_interval is sorted by position.
	// What if there're SVs with the same starting positions?
	
	if (bp[0].idx != 0)
	{
		cerr << "Error: Merged interval's first element is not the earliest." << endl;
	}
	if (bp[0].type)
	{
		cerr << "Error: Merged interval's earliest position is SV-end, not SV-start." << endl;
	}
	int next = 1;
	
	while((p = bam_plp_auto(plp, &tid, &pos, &n_plp))!=0)
	{
		//			cerr << "tid " << tid << " pos " << pos << endl;
		if (pos<startpos || pos >=endpos) continue;
		if (pos>=bp[next].pos)
		{
			if (bp[next].type)
			{
				//SV-end
				currset.erase(bp[next].idx);
				next++;
			}
			else
			{
				//SV-start
				currset.insert(bp[next].idx);
				next++;
			}
		}
		for(set<int>::iterator it=currset.begin(); it!=currset.end(); ++it)
		{
			sum[*it] += n_plp;
			cnt[*it]++;
		}
	}

	for(int i=0;i<n;++i)
	{
		X[i] = (cnt[i]>0) ? sum[i]/(double)cnt[i] : 0;
		GX[i] = X[i]; // TODO: do GC correction
	}
}

void bFile::read_pair(vector<sv> &m_interval, vector<double>& X)
{
}
