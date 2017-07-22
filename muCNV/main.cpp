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

#include <algorithm>

// #include "pFile.h"
#include "muCNV.h"

// TCLAP headers
#include "tclap/CmdLine.h"
#include "tclap/Arg.h"

#include "sam.h"

using namespace std;

bool bUseGL = false;
double P_THRESHOLD = 0.9;
double BE_THRESHOLD = 0.01;
double RO_THRESHOLD = 0.8;

typedef struct {     // auxiliary data structure
	samFile *fp;     // the file handle
	bam_hdr_t *hdr;  // the file header
	hts_itr_t *iter; // NULL if a region not specified
	int min_mapQ, min_len; // mapQ filter; length filter
} aux_t;


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
		if ( aux->min_len && bam_cigar2qlen(b->core.n_cigar, bam_get_cigar(b)) < aux->min_len ) continue;
		break;
	}
	return ret;
}

int main(int argc, char** argv)
{
	cout << "muCNV 0.5 -- Multi-sample genotyping of CNVs from depth data" << endl;
	cout << "(c) 2017 Goo Jun" << endl << endl;
	cerr.setf(ios::showpoint);

	bool bVerbose;
	string sInFile;
	string sEventFile;
	string sFamFile;
	string sIntervalFile;
	string sOutPrefix;
	string sChr;

	vector<string> sampleIDs;
	vector<string> sampleDirs;
	map<string, unsigned> hIdSex;

	vector<double> AvgDepth;
	//	vector<interval_t> del_intervals;
	//	vector<interval_t> dup_intervals;
	vector<Interval> del_intervals;
	vector<Interval> dup_intervals;
	vector<Interval> all_intervals;

	srand((unsigned int)time(NULL));

	// Parsing command-line arguments
	try 
	{
		TCLAP::CmdLine cmd("Command description message", ' ', "0.01");
		TCLAP::ValueArg<string> argIn("i","index","Input index file (sample ID, candidate VCF, BAM/CRAM)",false,"","string");
		TCLAP::ValueArg<string> argOut("o","out","Prefix for output filename",false,"muCNV","string");
		
		//		TCLAP::ValueArg<string> argInterval("n","interval","File containing list of candidate intervals",false,"","string");
		TCLAP::ValueArg<double> argPos("p","posterior","(Optional) Posterior probability threshold",false,0.9,"double");
		TCLAP::ValueArg<double> argBE("b","bayes","(Optional) Bayes error threshold",false,0.01,"double");

		TCLAP::ValueArg<double> argRO("r","reciprocal","(Optional) Reciprocal overlap threshold to merge candidate intervals (default: 0.8)",false,0.8,"double");
		//		TCLAP::SwitchArg switchGL("g","gl","Use likelihood instead of posterior to call",cmd,false);
		TCLAP::SwitchArg switchVerbose("v","verbose","Turn on verbose mode",cmd,false);

		cmd.add(argIn);
		cmd.add(argOut);
		cmd.add(argPos);

		//		cmd.add(argInterval);
		cmd.add(argBE);
		cmd.add(argRO);
		cmd.parse(argc, argv);

        // Index file:
		// three columns (tab separated)
		// sample ID, VCF with canddiate intervals, sequence file (BAM/CRAM)
        // SampleID, DepthFile(bgzipped, tabixed)
		
        
		sInFile = argIn.getValue();
		//sIntervalFile = argInterval.getValue();
		sOutPrefix = argOut.getValue();

		P_THRESHOLD = argPos.getValue();
		BE_THRESHOLD = argBE.getValue();
		RO_THRESHOLD = argRO.getValue();

        bVerbose = switchVerbose.getValue();

	}
	catch (TCLAP::ArgException &e) 
	{
		cerr << "Error: " << e.error() << " for arg " << e.argId() << endl;
		abort();
	}
	string fname = "/Users/gjun/data/cram/NA12878.fragment.cram";

	
	aux_t **data;
	data = (aux_t**)calloc(1, sizeof(aux_t*));
	data[0] = (aux_t*)calloc(1, sizeof(aux_t));
	data[0]->fp = hts_open(fname.c_str(), "r");

	
	int count=0;
	int rf = SAM_FLAG | SAM_RNAME | SAM_POS | SAM_MAPQ | SAM_CIGAR | SAM_SEQ | SAM_QUAL;
	
	if (hts_set_opt(data[0]->fp, CRAM_OPT_REQUIRED_FIELDS, rf))
	{
		cerr << "Failed to set CRAM_OPT_REQUIRED_FIELDS value" << endl;
		exit(1);
	}
	if (hts_set_opt(data[0]->fp, CRAM_OPT_DECODE_MD, 0))
	{
		fprintf(stderr, "Failed to set CRAM_OPT_DECODE_MD value\n");
		return
		1;
	}
	data[0]->min_mapQ = 20;
	data[0]->min_len = 0;
	data[0]->hdr = sam_hdr_read(data[0]->fp);
	if (data[0]->hdr == NULL)
	{
		cerr << "Cannot open CRAM/BAM header" << endl;
		exit(1);
	}
	
	
	hts_idx_t *idx = sam_index_load(data[0]->fp, fname.c_str());
	if (idx == NULL)
	{
		cerr << "Cannot open CRAM/BAM index" << endl;
		exit(1);
	}
	
	hts_itr_t *iter = sam_itr_querys(idx, data[0]->hdr, "20:6000000-6001000");
	hts_idx_destroy(idx);
	if (iter == NULL)
	{
		cerr << "Can't parse region" << endl;
		exit(1);
	}

	bam_mplp_t mplp = bam_mplp_init(1, read_bam, (void**)data);
	bam1_t *aln = bam_init1();

	int *n_plp = (int*)calloc(1, sizeof(int));
	const bam_pileup1_t **plp = (const bam_pileup1_t**)calloc(1, sizeof(bam_pileup1_t*));
	int tid, pos;
	int last_pos = -1;
	bam_hdr_t *h = data[0]->hdr;
	
	while(bam_mplp_auto(mplp, &tid, &pos, n_plp, plp)>0)
	{
		if (pos<6000000 || pos >=6000100) continue;
		cout << h->target_name[tid] << "\t" << pos +1;
		int j,m = 0;
		for(j=0;j<n_plp[0]; ++j)
		{
			const bam_pileup1_t *p = plp[0] + j;
			if (p->is_del || p->is_refskip) ++m;
			else if (bam_get_qual(p->b)[p->qpos] < 20) ++m;
		}
		cout << "\t" << n_plp[0] - m << endl;
		
	}
	
	free(n_plp);
	free(plp);
	bam_mplp_destroy(mplp);
	
	//while(sam_itr_next(myFile, iter, aln)>=0)
	{
		//	cout << myHeader->target_name[aln->core.tid];
		//	cout << aln->core.pos << endl;
		//count++;
	}
	
	//	cout << "Number of reads : " << count << endl;

/*

	unsigned n_sample = 0;
	unsigned n_del = 0;
	unsigned n_dup = 0;

	vector<string> depthFiles;

	vector<Interval> del_events;
	vector<Interval> dup_events;

	vector<string> eventFiles;

    // Read sample index
    SampleList samples;
    samples.readIndex(sInFile);
    
//	readIndex(sInFile, sampleIDs, sampleDirs, eventFiles, depthFiles);
    
//	n_sample = (unsigned) sampleIDs.size();
    
    if (bVerbose) cerr << samples.n_sample << " samples identified\n"<<endl;
  
	AvgDepth.resize(n_sample, 0);

	for(unsigned i=0; i<n_sample;++i)
	{
		// Read average depth per sample
		AvgDepth[i] = getAvgDepth(sampleIDs[i], sampleDirs[i]);
		if (bVerbose) cerr << "\rSample " << i+1 << " / "<< n_sample << ": " << sampleIDs[i] << ", average Depth " << AvgDepth[i];

		// Read CNV files
		string sSegmentFile = getCNVsegmentFileName(sampleIDs[i], sampleDirs[i]);
		readInterval(sSegmentFile, 1, del_events, dup_events);

		// Read event files
		readInterval(eventFiles[i], 2, del_events, dup_events);
	}

	if (bVerbose) 
	{
		cerr << endl << del_events.size() << " deletion and ";
		cerr << dup_events.size() << " duplication candidate intervals identified" <<  endl;
	}

	if (sIntervalFile != "")
	{
		readInterval(sIntervalFile, 3, del_events, dup_events);
	}

	// Sort intervals, remove duplicates, and cluster overlapping intervals
	sort(del_events.begin(), del_events.end(), compareIntervals);

	// To do: remove duplicates before clustering -- will make things faster, but results should be the same
	//		remove_duplicates(del_events);
	clusterIntervals(del_events, del_intervals);

	sort(dup_events.begin(), dup_events.end(), compareIntervals);
	//		remove_duplicates(dup_events);
	clusterIntervals(dup_events, dup_intervals);

	// merge dup_intervals and del_intervals 
	all_intervals = del_intervals;
	uint32_t n_sv = (uint32_t) all_intervals.size();

//	n_del = del_intervals.size();
//	n_dup = dup_intervals.size();

	// Store average depth per interval per sample
	vector< vector<double> > X(n_del, vector<double>(n_sample,0));   
	vector< vector<double> > Y(n_dup, vector<double>(n_sample,0));   

	// Editing: read full length from one interval at a time, and do breakpoint refinement // Breakpoint refinement will be implemented in later version.
	//
//
//   Read  2 x interval (0.5 before, 0.5 after) from all samples
//   Do clustering on the original or smallest (? confident ) interval
//   If passed BIC and initial Bayes error criteria, Refine breakpoint
//   Normalize each read depth by average depth across samples ( how to adjust for allele count ? )
//   Re-do clustering, evaluate Bayes error

//	readDepth(depthFiles, del_intervals, X, AvgDepth);
//	readDepth(depthFiles, dup_intervals, Y, AvgDepth);

	uint64_t max_svlen = 0;
	for(unsigned i=0;i<n_sv;++i)
	{
		uint64_t L = all_intervals[i].len();
		// Temporary, for quick run
		if (L>1e6)
		{
			all_intervals.erase(all_intervals.begin()+i);
		}
		else if (max_svlen < L )
		{
			max_svlen = L;
		}
	}

	// Reserve memory space to read all intervals -- up to 2 * max_svlen for each sample
	vector< vector<double> > D(n_sample, vector<double> (1,0));
	vector<double>  AvgD(n_sample, 0);

	if (bVerbose)
	{
		cerr << n_sv << "intervals with max SVLEN " << max_svlen << endl;
	}

	for(unsigned i=0;i<n_sample;++i)
	{
		D[i].resize(max_svlen + 2000);
	}

	string sVcfFile = sOutPrefix + ".vcf" ;

	FILE *vcfFile;
	vcfFile = fopen(sVcfFile.c_str(), "wt");

	write_vcfheader(vcfFile, sampleIDs);

	// To Do : make a class for VCF writer, and add methods to write header and SV entries

	for(unsigned i=0;i<n_sv;++i)
	{
		// Read depth for i-th interval
	//TMPTMP	readDepth(depthFiles, all_intervals[i],  D, AvgD, AvgDepth);

		// Do clustering 
		if (all_intervals[i].sv_type == "DEL")
		{
//TMPTMP			call_deletion(all_intervals[i], AvgD, AvgDepth, vcfFile);
		}
		else
		{
//			call_duplication(all_intervals[i], vcfFile);
//		Test 
		}

		// Refine breakpoints
		
		// Re-do clustering

		// Write to VCF
	}

	call_deletions(X, AvgDepth, sampleIDs,  del_intervals, vcfFile);
	call_duplications(Y, AvgDepth, sampleIDs, dup_intervals, vcfFile);

	fclose(vcfFile);
*/
	return 0;
}
