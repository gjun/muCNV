//
//  sv.cpp
//  muCNV
//
//  Created by Goo Jun on 7/23/17.
//  Copyright © 2017 Goo Jun. All rights reserved.
//

#include "sv.h"
#include "common.h"
#include <algorithm>

svType svTypeNum(int t)
{
    switch(t)
    {
        case 0:
            return DEL;
            break;
        case 1:
            return DUP;
            break;
        case 2:
            return INV;
            break;
        case 3:
            return CNV;
            break;
        case 4:
            return INS;
            break;
        case 5:
            return BND;
            break;
    }
    return DEL;
}

std::string svTypeName(svType t)
{
    std::string S[6] = {"DEL", "DUP", "INV", "CNV", "INS", "BND"};
    return S[t];
}

breakpoint::breakpoint()
{
    chrnum = 0;
    pos = 0;
    bptype = 0;
    idx = 0;
}

bool breakpoint::operator < (const breakpoint& b) const
{
    return (chrnum<b.chrnum || (chrnum == b.chrnum && pos<b.pos));
}

bool breakpoint::operator == (const breakpoint& b) const
{
    return (chrnum == b.chrnum && pos == b.pos);
}

bool breakpoint::operator <= (const breakpoint& b) const
{
	return (chrnum<b.chrnum || (chrnum==b.chrnum && pos<=b.pos));
//    return (*this <b || *this ==b);
}

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

int find_start(std::vector<sv> &L, int pos)
{
	// Assume list is sorted according to pos
	int left = 0;
	int right = (int) L.size();
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
	svtype = DEL;
	//chr = "";
	chrnum = -1;
	pos = -1;
    supp=0;
    n_dp = 0;
    dp_sum = 0;
    dp = 0;
//	ci_pos.first = 0;
//	ci_pos.second = 0;
//	ci_end.first = 0;
//	ci_end.second = 0;
}

bool sv::operator == (const sv& s) const
{
	return (chrnum == s.chrnum && pos == s.pos && end==s.end && svtype == s.svtype);
}

void sv::print(void)
{
    printf("%d:%d-%d_%s", chrnum, pos, end, svTypeName(svtype).c_str());
}


void pick_sv_from_merged(sv &picked, std::vector<sv> &merged)
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

	std::vector<sv> max_svs;
	std::vector<int> starts;
	std::vector<int> ends;
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

double RO(sv &x, sv &y)
{
    double l = 0;
    double L = 0;
    
    if (y.pos > x.end)
    {
        return 0;
    }
    else if (x.pos > y.end)
    {
        return 0;
    }
    
    bool bStart = (x.pos < y.pos);
    bool bEnd = (x.end < y.end);
    
    l = (bEnd ? x.end : y.end) -  (bStart ? y.pos : x.pos);
    L = (bEnd ? y.end : x.end) - (bStart ? x.pos : y.pos);
    return (l/L);
}

int find_overlap_sv(sv &S , std::vector<sv>& dels)
{
    int idx = 0;
    double max_RO = 0;
    double max_idx = 0;
    
    if (S.pos > 5000000)
    {
        idx = find_start(dels, S.pos - 5000000);
    }
    while(dels[idx].pos < S.end && idx<(int)dels.size())
    {
        if (dels[idx].end > S.pos)
        {
            double r = RO(dels[idx], S);
            if (r>max_RO)
            {
                max_RO = r;
                max_idx = idx;
            }
        }
        idx++;
    }
    if (max_RO>0.6)
    {
        return max_idx;
    }
    else
    {
        return -1;
    }
}

bool in_centrome(int chrnum, int pos)
{
    // GRCh38 Centromere coordinates , hardcoded
    int centro[24][2] = {
        {121700000, 125100000},
        {91800000, 96000000},
        {87800000, 94000000},
        {48200000, 51800000},
        {46100000, 51400000},
        {58500000, 62600000},
        {58100000, 62100000},
        {43200000, 47200000},
        {42200000, 45500000},
        {38000000, 41600000},
        {51000000, 55800000},
        {33200000, 37800000},
        {16500000, 18900000},
        {16100000, 18200000},
        {17500000, 20500000},
        {35300000, 38400000},
        {22700000, 27400000},
        {15400000, 21500000},
        {24200000, 28100000},
        {25700000, 30400000},
        {10900000, 13000000},
        {13700000, 17400000},
        {58100000, 63800000},
        {10300000, 10600000}
    };
    if (pos >= centro[chrnum-1][0] && pos <= centro[chrnum-1][1])
        return true;
    else
        return false;
}

bool in_centrome(sv &S)
{
    // GRCh38 Centromere coordinates , hardcoded
    int centro[24][2] = {
        {121700000, 125100000},
        {91800000, 96000000},
        {87800000, 94000000},
        {48200000, 51800000},
        {46100000, 51400000},
        {58500000, 62600000},
        {58100000, 62100000},
        {43200000, 47200000},
        {42200000, 45500000},
        {38000000, 41600000},
        {51000000, 55800000},
        {33200000, 37800000},
        {16500000, 18900000},
        {16100000, 18200000},
        {17500000, 20500000},
        {35300000, 38400000},
        {22700000, 27400000},
        {15400000, 21500000},
        {24200000, 28100000},
        {25700000, 30400000},
        {10900000, 13000000},
        {13700000, 17400000},
        {58100000, 63800000},
        {10300000, 10600000}
    };
    //  int hetero[2] = {51078349, 54425074};
    
    //  if (S.pos >= centro[S.chrnum-1][0]-300000 && S.pos <= centro[S.chrnum-1][1]+300000)
    if (S.pos >= centro[S.chrnum-1][0] && S.pos <= centro[S.chrnum-1][1])
        return true;
    if (S.end >= centro[S.chrnum-1][0] && S.end <= centro[S.chrnum-1][1])
        return true;
    if (S.pos <= centro[S.chrnum-1][0] && S.end >= centro[S.chrnum-1][0])
        return true;
    if (S.pos <= centro[S.chrnum-1][1] && S.end >= centro[S.chrnum-1][1])
        return true;
    //  if (S.chrnum == 7  && S.pos >= hetero[0] && S.pos <= hetero[1]+1000000)
    //      return true;
    //  if (S.chrnum == 7  && S.end >= hetero[0] && S.end <= hetero[1])
    //      return true;
    return false;
}
