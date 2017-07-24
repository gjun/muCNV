//
//  sv.cpp
//  muCNV
//
//  Created by Goo Jun on 7/23/17.
//  Copyright © 2017 Goo Jun. All rights reserved.
//

#include "muCNV.h"

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

bool sv::operator == (const sv& s) const
{
	return (chr == s.chr && pos == s.pos && end==s.end && svtype == s.svtype);
}
