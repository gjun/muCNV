//
//  sv.cpp
//  muCNV
//
//  Created by Goo Jun on 7/23/17.
//  Copyright © 2017 Goo Jun. All rights reserved.
//

#include "muCNV.h"

sv::sv()
{
	svtype = "NA";
	source = "NA";
	chr = -1;
	pos = -1;
	end = -1;
	ci_pos.first = 0;
	ci_pos.second = 0;
	ci_end.first = 0;
	ci_end.second = 0;
}

bool sv::operator < (const sv& s) const
{
	if (chr==s.chr)
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
		return(chr<s.chr);
	}
}
