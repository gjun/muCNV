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

//void merge_svs(vector<sv> &candidates , vector< vector<sv> > &merged_candidates)
void merge_svs(vector<sv> &candidates , vector<int> &idxs)
{
	
	int curr = 0;
	
	// 1. sort intervals
	std::sort(candidates.begin(), candidates.end());
	candidates.erase( std::unique( candidates.begin(), candidates.end() ), candidates.end() );
	
	// find a block of sv intervals with overlap
	while(curr<(int)candidates.size())
	{
		int block_end = candidates[curr].end;
		int last_idx = curr;
		idxs.push_back(curr);
		
		int cnt = 0;
		while(++last_idx < (int)candidates.size() && candidates[last_idx].chr == candidates[curr].chr &&  candidates[last_idx].pos < block_end && cnt<500) // MAX candidates in a single interval: 200
		{
			if (block_end<candidates[last_idx].end)
			{
				block_end = candidates[last_idx].end;
			}
			cnt++;
		}
		
		/*
		vector<sv> t;
		for(int i=curr; i<last_idx; ++i)
		{
			t.push_back(candidates[i]);
		}
		merged_candidates.push_back(t);
		*/

		curr=last_idx;
	}
}

void cluster_svs(vector<sv> &candidates , vector< vector<sv> > &merged_candidates)
{

	int curr = 0;

	// 1. sort intervals
	std::sort(candidates.begin(), candidates.end());
	candidates.erase( std::unique( candidates.begin(), candidates.end() ), candidates.end() );

	// find a block of sv intervals with overlap
	while(curr<(int)candidates.size())
	{
		//		cerr << "curr : " << curr<< endl;
		int block_end = candidates[curr].end;

		//		cerr << "block_end : " << block_end << endl;
		int last_idx = curr;

		while(++last_idx < (int)candidates.size() && candidates[last_idx].chr == candidates[curr].chr &&  candidates[last_idx].pos < block_end)
		{
			//			cerr<< "pos : " << candidates[last_idx].pos  << ", end : " << candidates[last_idx].end << endl;
			if (block_end<candidates[last_idx].end)
			{
				block_end = candidates[last_idx].end;

				//				cerr << "block_end : " << block_end << endl;
			}
		}

		///		cerr << "last_idx : " << last_idx << endl; // TEMP : code reached up to here

		int n = last_idx - curr;
		if (n>2)
		{
			vector< vector <double> >  D; // This might be problematic with dense intervals -- change it to malloc or std::vector
			D.resize(n);
			for(int i=0;i<n;++i)
			{
				D[i].resize(n);
			}

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
					// TODO:: Do not store all candidate intervals, just hold the median one
					//     :: Maybe should be processed by a preprocessor -- not to deal with it here
				}
			}

		}
		else
		{
			vector<sv> t;
			t.push_back(candidates[curr]);

			merged_candidates.push_back(t);
		}
		curr=last_idx;
	}
}
