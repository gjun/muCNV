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

sv::sv()
{
	svtype = "";
	chr = "";
	chrnum = -1;
	pos = -1;
	ci_pos.first = 0;
	ci_pos.second = 0;
	ci_end.first = 0;
	ci_end.second = 0;
	
}

bool sv::operator == (const sv& s) const
{
	return (chrnum == s.chrnum && pos == s.pos && end==s.end && svtype == s.svtype);
}

void pick_sv_from_merged(vector<sv> &sv_list, vector<sv> &merged)
{
	if (merged.size() == 1) // trivial case
	{
		sv_list.push_back(merged[0]);
		return;
	}
	
	// let's do the position only
	
	vector<double> pos, end;

	string svtype = "DEL";
	sv new_sv;
	
	for(int i=0;i<(int)merged.size();++i)
	{
		pos.push_back(merged[i].pos);
		end.push_back(merged[i].end);
		if (merged[i].svtype != svtype)
		{
			svtype = "CNV";
		}
	}
	sort(pos.begin(), pos.end());
	sort(end.begin(), end.end());
	

	new_sv.svtype = svtype;
	new_sv.chr = merged[0].chr;
	new_sv.chrnum = merged[0].chrnum;
	if (pos.size()%2)
	{
		int idx =  floor(pos.size()/2.0);
		new_sv.pos = pos[idx];
		new_sv.end = end[idx];
	}
	else
	{
		int idx =  floor(pos.size()/2.0);
		new_sv.pos = floor((pos[idx]+pos[idx-1])/2.0);
		new_sv.end =  floor((end[idx]+end[idx-1])/2.0);

	}
	if (new_sv.pos >= new_sv.end)
	{
		new_sv.pos = pos[0];
		new_sv.end = end[end.size()-1];
	}
	new_sv.ci_pos.first = pos[0] - new_sv.pos;
	new_sv.ci_pos.second = pos[pos.size()-1] - new_sv.pos;
	new_sv.ci_end.first = end[0] - new_sv.end;
	new_sv.ci_end.second = end[pos.size()-1] - new_sv.end;
	
	sv_list.push_back(new_sv);
}
