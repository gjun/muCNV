//
//  sv.cpp
//  muCNV
//
//  Created by Goo Jun on 7/23/17.
//  Copyright © 2017 Goo Jun. All rights reserved.
//

#include "muCNV.h"
#include <algorithm>

bool sv::operator < (const sv& s) const
{
	if (chrnum==s.chrnum)
	{
		if (pos==s.pos)
		{
			return(end<s.end);
		}
		else
		{
			return(pos<s.pos);
		}
	}
	else
	{
		return(chrnum<s.chrnum);
	}
}

int find_start(vector<sv> &L, int pos)
{
	// Assume list is sorted according to pos
	int left = 0;
	int right = L.size();
	int idx = (left+right)/2;

	while(right > left)
	{
		if (L[idx].pos > pos)
		{
			right = idx;
		}
		else if (L[idx].pos < pos)
		{
			left = idx+1;
		}
		else if (L[idx].pos == pos)
		{
			while(idx>0 && L[--idx].pos == pos);
			return idx;
		}
		idx = (left+right)/2;
	}
	return idx;
}

sv::sv()
{
	svtype = "";
	chr = "";
	chrnum = -1;
	pos = -1;
//	ci_pos.first = 0;
//	ci_pos.second = 0;
//	ci_end.first = 0;
//	ci_end.second = 0;
	
}

bool sv::operator == (const sv& s) const
{
	return (chrnum == s.chrnum && pos == s.pos && end==s.end && svtype == s.svtype);
}

void pick_sv_from_merged(sv &picked, vector<sv> &merged)
{

	if (merged.size() == 1) // trivial case
	{
		picked = merged[0];
		return;
	}
	
	int max_supp = 0;
	// find ones with best supp
	for(int i=0;i<(int)merged.size(); ++i)
	{
		if (merged[i].supp > max_supp)
			max_supp = merged[i].supp;
	}

	vector<sv> max_svs;
	vector<int> starts;
	vector<int> ends;
	for(int i=0;i<(int)merged.size(); ++i)
	{
		if (merged[i].supp == max_supp)
		{
			max_svs.push_back(merged[i]);
			starts.push_back(merged[i].pos);
			ends.push_back(merged[i].end);
		}
	}
	int med_start = median(starts);
	int med_end = median(ends);

	int best_idx = 0;
	double best_dist = DBL_MAX;

	for(int i=0; i<(int)max_svs.size(); ++i)
	{
		double D = (med_start-max_svs[i].pos)*(med_start-max_svs[i].pos) + (med_end-max_svs[i].end)*(med_end-max_svs[i].end);
		if (D<best_dist)
		{
			best_idx = i;
			best_dist = D;
		}
	}
	picked = max_svs[best_idx];

//	new_sv.ci_pos.first = pos[0] - new_sv.pos;
//	new_sv.ci_pos.second = pos[pos.size()-1] - new_sv.pos;
//	new_sv.ci_end.first = end[0] - new_sv.end;
//	new_sv.ci_end.second = end[pos.size()-1] - new_sv.end;
}
