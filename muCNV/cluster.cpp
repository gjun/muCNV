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

#include "muCNV.h"
#include <algorithm>

extern double RO_THRESHOLD;

double RO(sv x, sv y)
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

void cluster_svs(vector<sv> &candidates , vector< vector<sv> > &merged_candidates)
{

	int curr = 0;

	// 1. sort intervals
	std::sort(candidates.begin(), candidates.end());
	candidates.erase( std::unique( candidates.begin(), candidates.end() ), candidates.end() );
	
	// find a block of sv intervals with overlap
	while(curr<candidates.size())
	{
		int block_end = candidates[curr].end;
		
		int last_idx = curr;
		
		while(++last_idx < candidates.size() && candidates[last_idx].pos < block_end)
		{
			if (block_end<candidates[last_idx].end)
			{
				block_end = candidates[last_idx].end;
			}
		}

		int n = last_idx - curr;
		double D[n][n];
		double max_RO = 0;
		int clusters[n];
		int max_i = 0;
		int max_j = 0;
		
		for(int i=0;i<n;++i)
		{
			clusters[i] = i;
		}
		
		for(int i=0;i<n;++i)
		{
			for(int j=i+1;j<n;++j)
			{
				D[i][j] = D[j][i] = RO(candidates[curr+i], candidates[curr+j]);
				if (D[i][j] > max_RO)
				{
					max_RO = D[i][j];
					max_i = i;
					max_j = j;
				}
			}
			D[i][i] = 0;
		}

		while(max_RO>RO_THRESHOLD)
		{
			// Merge clusters
			clusters[max_j] = clusters[max_i];
			// Update distances
			for(int i=0;i<n;++i)
			{
				if (i != max_i)
				{
					D[max_i][i] = D[i][max_i] = (D[max_i][i] > D[max_j][i]) ? D[max_i][i] : D[max_j][i];
				}
			}
			
			for(int i=0;i<n;++i)
			{
				D[max_j][i] = D[i][max_j] = 0;
			}
			
			max_RO = 0;
			for(int i=0;i<n;++i)
			{
				for(int j=i+1;j<n;++j)
				{
					if (D[i][j] > max_RO)
					{
						max_RO = D[i][j];
						max_i = i;
						max_j = j;
					}
				}
			}
		}
		
		for(int i=0;i<n;++i)
		{
			vector<sv> t;
			for(int j=0; j<n; ++j)
			{
				if (clusters[j] == i)
				{
					t.push_back(candidates[j+curr]);
				}
			}
			if (t.size()>0)
			{
				merged_candidates.push_back(t); // Maybe not the best way with unnecessary copies
			}
		}
		curr=last_idx;
	}
}
