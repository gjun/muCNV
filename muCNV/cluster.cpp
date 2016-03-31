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
#include <set>

extern double RO_THRESHOLD;

double RO(interval_t x, interval_t y)
{
	double l = 0;
	double L = 0;

	if (y.first > x.second)
	{
		return 0;
	}
	else if (x.first > y.second)
	{
		return 0;
	}

	bool bStart = (x.first < y.first);
	bool bEnd = (x.second < y.second);

	l = (bEnd ? x.second : y.second) -  (bStart ? y.first : x.first);
	L = (bEnd ? y.second : x.second) - (bStart ? x.first : y.first);
	return (l/L);
}

string merge_sources(const char* s1, const char *s2) 
{
	vector<string> tokens1;
	vector<string> tokens2;
	set<string> sources;
	string ret;

	pFile::tokenizeLine(s1, ",", tokens1);
	pFile::tokenizeLine(s2, ",", tokens2);

	for(unsigned i=0;i<tokens1.size();++i)
	{
		sources.insert(tokens1[i]);
	}
	for(unsigned i=0;i<tokens2.size();++i)
	{
		sources.insert(tokens2[i]);
	}

	set<string>::iterator it=sources.begin();
	ret = *it;

	++it;
	for(;it!=sources.end();++it)
	{
		ret = ret + "," + (*it);
	}
	return ret;
}

void clusterIntervals(vector<Interval> &events, vector<Interval> &intervals) 
{
	// Dummy last interval
	Interval last_event;

	last_event.T.first = events.back().T.second + 10000;
	last_event.T.second =  last_event.T.first + 20000;

	// This should be popped back 
	events.push_back(last_event);

	Interval curr(events.front());

	unsigned k0= 0;

	for(unsigned k=1; k<events.size(); ++k)
	{
		if (events[k].T.first > curr.T.second)
		{
//			cerr << "Cluster detected : " << k-k0 << " events in " << curr.first << " to " << curr.second << endl;

			vector<Interval> cluster(events.begin()+k0, events.begin()+k);

			double max_RO = 1;
			unsigned k1=0, k2=0;

			if (k-k0>1)
			{
				vector<double> col_max(cluster.size(),0);
				vector<double> row_max(cluster.size(),0);
				vector<unsigned> col_idx(cluster.size(),0);
				vector<unsigned> row_idx(cluster.size(),0);

				max_RO = -1;
				// Initiate similarity matrix
				for(unsigned m=0; m<cluster.size()-1; ++m)
				{
					for(unsigned n=m+1; n<cluster.size(); ++n)
					{
						double r = RO(cluster[m].T, cluster[n].T);
						if (r > max_RO)
						{
							max_RO = r;
							k1 = m;
							k2 = n;
						}
						if (r > row_max[m])
						{
							row_max[m] = r;
							row_idx[m] = n;
						}
						if (r > col_max[n])
						{
							col_max[n] = r;
							col_idx[n] = m;
						}
					}
				}

				while(max_RO>RO_THRESHOLD && cluster.size()>1)
				{
					// Merge intervals, always the latter one is deleted, since intervals are sorted
					cluster[k1].T.second = max(cluster[k1].T.second, cluster[k2].T.second);

					// merge source, cipos, ciend
					int32_t posdiff = (int32_t)(cluster[k1].T.first - cluster[k2].T.first);
					int32_t enddiff = (int32_t)(cluster[k1].T.second - cluster[k2].T.second);
					/*
					1) enddiff <0 
					|------------------|
					      |---------------------|
					1: ci_end1.first + enddiff
					2: ci_end2.first 
					new_ci_end.first  = min(1,2)
					3: ci_end1.second + enddiff
					4: ci_end2.second 
					new_ci_end.second = max(3,4)

					2) enddiff > 0     
				 	|--------------------|
						|-------------|
					1: ci_end1.first + enddiff
					2: ci_end2.first 
					new_ci_end.first  = min(1,2)

					3: ci_end1.second + enddiff
					4: ci_end2.second 
					new_ci_end.second = max(3,4)
					*/
					cluster[k1].ci_pos.first = min(cluster[k1].ci_pos.first, posdiff + cluster[k2].ci_pos.first);
					cluster[k1].ci_pos.second = max(cluster[k1].ci_pos.second, posdiff + cluster[k2].ci_pos.second);

					cluster[k1].ci_end.first = min(cluster[k1].ci_end.first + enddiff, cluster[k2].ci_end.first);
					cluster[k1].ci_end.second = max(cluster[k1].ci_end.second, enddiff + cluster[k2].ci_end.second);

					cluster[k1].source = merge_sources(cluster[k1].source.c_str(), cluster[k2].source.c_str());

					cluster.erase(cluster.begin()+k2);
					col_max.erase(col_max.begin()+k2);
					col_idx.erase(col_idx.begin()+k2);
					row_max.erase(row_max.begin()+k2);
					row_idx.erase(row_idx.begin()+k2);

					unsigned prev_k1 = k1;
					unsigned prev_k2 = k2;
					double max_COL = 0;
					double max_ROW = 0;
					vector<unsigned> col_ties, row_ties;

					for(unsigned m=0;m<cluster.size();++m)
					{
						if (m==prev_k1)
						{
							row_max[m] = 0;
							row_idx[m] = 0;
						}
						else if (row_idx[m] == prev_k2 || row_idx[m] == prev_k1)
						{
							row_max[m] = 0;
							for(unsigned n=m+1;n<cluster.size();++n)
							{
								double r2 = RO(cluster[m].T, cluster[n].T);
								if (r2>row_max[m])
								{
									row_idx[m] = n; 
									row_max[m] = r2;
								}
							}
						}
						else if (row_idx[m] > prev_k2) 
						{
							--row_idx[m];
						}
					}

					for(unsigned m=1;m<cluster.size(); ++m)
					{
						if (m==prev_k1)
						{
							col_max[m] = 0;
							col_idx[m] = 0;
						}
						else if (col_idx[m] == prev_k2 || col_idx[m] == prev_k1)
						{
							col_max[m] = 0;
							for(unsigned n=0;n<m;++n)
							{
								double r2 = RO(cluster[m].T, cluster[n].T);

								if (r2>col_max[m])
								{
									col_idx[m] = n; 
									col_max[m] = r2;
								}
							}
						}
						else if (col_idx[m] > prev_k2)
						{
							--col_idx[m];
						}
					}
					for(unsigned m=0; m<cluster.size(); ++m)
					{
						if (m!=prev_k1)
						{
							double r = RO(cluster[m].T, cluster[prev_k1].T);

							if (m<prev_k1)
							{
								if (r>row_max[m])
								{
									row_max[m] = r;
									row_idx[m] = prev_k1;
								}
								if (r > col_max[prev_k1])
								{
									col_max[prev_k1] = r;
									col_idx[prev_k1] = m;
								}
							}
							else
							{
								if (r>row_max[prev_k1])
								{
									row_max[prev_k1] = r;
									row_idx[prev_k1] = m;
								}
								if (r > col_max[m])
								{
									col_max[m] = r;
									col_idx[m] = prev_k1;
								}
		
							}

							if (row_max[m] > max_ROW)
							{
								max_ROW = row_max[m];
								k1 = m;
								row_ties.clear();
							}
							else if (row_max[m] == max_ROW)
							{
								if (row_ties.empty())
								{
									row_ties.push_back(k1);
								}
								row_ties.push_back(m);
							}

							if (col_max[m] > max_COL)
							{
								max_COL = col_max[m];
								k2 = m;
								col_ties.clear();
							}
							else if (col_max[m] == max_COL)
							{
								if (col_ties.empty())
								{
									col_ties.push_back(k2);
								}
								col_ties.push_back(m);
							}
						}
					}

					if (row_max[prev_k1] > max_ROW)
					{
						max_ROW = row_max[prev_k1];
						k1 = prev_k1;
						row_ties.clear();
					}
					else if (row_max[prev_k1] == max_ROW)
					{
						if (row_ties.empty())
						{
							row_ties.push_back(k1);
						}
						row_ties.push_back(prev_k1);
					}

					if (col_max[prev_k1] > max_COL)
					{
						max_COL = col_max[prev_k1];
						k2 = prev_k1;
						col_ties.clear();
					}
					else if (col_max[prev_k1] == max_COL)
					{
						if (col_ties.empty())
						{
							col_ties.push_back(k2);
						}
						col_ties.push_back(prev_k1);
					}
	
					if (!row_ties.empty() || !col_ties.empty())
					{
						if (row_ties.empty())
						{
							row_ties.push_back(k1);
						}
						else if (col_ties.empty())
						{
							col_ties.push_back(k2);
						}

						double max_RO = 0;
						for(unsigned m=0;m<row_ties.size();++m)
						{
							for(unsigned n=0;n<col_ties.size();++n)
							{
								if (row_ties[m]<col_ties[n])
								{
									double r = RO(cluster[row_ties[m]].T, cluster[col_ties[n]].T);
									if (r>max_RO)
									{
										max_RO = r;
										k1 = row_ties[m];
										k2 = col_ties[n];
									}
								}
							}
						}
					}
					max_RO = max_ROW;
				}
			}

			intervals.insert(intervals.end(), cluster.begin(), cluster.end());
			k0 = k;
			curr.T.first = events[k].T.first;
			curr.T.second = events[k].T.second;
		}
		else if (curr.T.second < events[k].T.second)
		{
			curr.T.second = events[k].T.second;
		}
	}

}

/*
void clusterIntervals(vector<interval_t> &events, vector<interval_t> &intervals) 
{

	// Dummy last interval
	interval_t last_event;

	last_event.first = events.back().second + 10000;
	last_event.second =  last_event.first + 20000;

	events.push_back(last_event);

	interval_t curr(events.front());

	unsigned k0= 0;

	for(unsigned k=1; k<events.size(); ++k)
	{
		if (events[k].first > curr.second)
		{
//			cerr << "Cluster detected : " << k-k0 << " events in " << curr.first << " to " << curr.second << endl;

			vector<interval_t> cluster(events.begin()+k0, events.begin()+k);
			vector<interval_t> minmax(events.begin()+k0, events.begin()+k);

			double max_RO = 1;
			unsigned k1=0, k2=0;

			if (k-k0>1)
			{
				vector<double> col_max(cluster.size(),0);
				vector<double> row_max(cluster.size(),0);
				vector<unsigned> col_idx(cluster.size(),0);
				vector<unsigned> row_idx(cluster.size(),0);

				max_RO = -1;
				// Initiate similarity matrix
				for(unsigned m=0; m<cluster.size()-1; ++m)
				{
					for(unsigned n=m+1; n<cluster.size(); ++n)
					{
						double r = RO(cluster[m], cluster[n]);
						if (r > max_RO)
						{
							max_RO = r;
							k1 = m;
							k2 = n;
						}
						if (r > row_max[m])
						{
							row_max[m] = r;
							row_idx[m] = n;
						}
						if (r > col_max[n])
						{
							col_max[n] = r;
							col_idx[n] = m;
						}
					}
				}

				while(max_RO>RO_THRESHOLD && cluster.size()>1)
				{
					// Merge intervals
					cluster[k1].second = (cluster[k1].second > cluster[k2].second) ? cluster[k1].second : cluster[k2].second;
					minmax[k1].second = (cluster[k1].second > cluster[k2].second) ? cluster[k2].second : cluster[k1].second;
					minmax[k1].first = (cluster[k1].first > cluster[k2].first) ? cluster[k1].first : cluster[k2].first;

					cluster.erase(cluster.begin()+k2);
					col_max.erase(col_max.begin()+k2);
					col_idx.erase(col_idx.begin()+k2);
					row_max.erase(row_max.begin()+k2);
					row_idx.erase(row_idx.begin()+k2);

					unsigned prev_k1 = k1;
					unsigned prev_k2 = k2;
					double max_COL = 0;
					double max_ROW = 0;
					vector<unsigned> col_ties, row_ties;

					for(unsigned m=0;m<cluster.size();++m)
					{
						if (m==prev_k1)
						{
							row_max[m] = 0;
							row_idx[m] = 0;
						}
						else if (row_idx[m] == prev_k2 || row_idx[m] == prev_k1)
						{
							row_max[m] = 0;
							for(unsigned n=m+1;n<cluster.size();++n)
							{
								double r2 = RO(cluster[m], cluster[n]);
								if (r2>row_max[m])
								{
									row_idx[m] = n; 
									row_max[m] = r2;
								}
							}
						}
						else if (row_idx[m] > prev_k2) 
						{
							--row_idx[m];
						}
					}

					for(unsigned m=1;m<cluster.size(); ++m)
					{
						if (m==prev_k1)
						{
							col_max[m] = 0;
							col_idx[m] = 0;
						}
						else if (col_idx[m] == prev_k2 || col_idx[m] == prev_k1)
						{
							col_max[m] = 0;
							for(unsigned n=0;n<m;++n)
							{
								double r2 = RO(cluster[m], cluster[n]);

								if (r2>col_max[m])
								{
									col_idx[m] = n; 
									col_max[m] = r2;
								}
							}
						}
						else if (col_idx[m] > prev_k2)
						{
							--col_idx[m];
						}
					}
					for(unsigned m=0; m<cluster.size(); ++m)
					{
						if (m!=prev_k1)
						{
							double r = RO(cluster[m], cluster[prev_k1]);

							if (m<prev_k1)
							{
								if (r>row_max[m])
								{
									row_max[m] = r;
									row_idx[m] = prev_k1;
								}
								if (r > col_max[prev_k1])
								{
									col_max[prev_k1] = r;
									col_idx[prev_k1] = m;
								}
							}
							else
							{
								if (r>row_max[prev_k1])
								{
									row_max[prev_k1] = r;
									row_idx[prev_k1] = m;
								}
								if (r > col_max[m])
								{
									col_max[m] = r;
									col_idx[m] = prev_k1;
								}
		
							}

							if (row_max[m] > max_ROW)
							{
								max_ROW = row_max[m];
								k1 = m;
								row_ties.clear();
							}
							else if (row_max[m] == max_ROW)
							{
								if (row_ties.empty())
								{
									row_ties.push_back(k1);
								}
								row_ties.push_back(m);
							}

							if (col_max[m] > max_COL)
							{
								max_COL = col_max[m];
								k2 = m;
								col_ties.clear();
							}
							else if (col_max[m] == max_COL)
							{
								if (col_ties.empty())
								{
									col_ties.push_back(k2);
								}
								col_ties.push_back(m);
							}
						}
					}

					if (row_max[prev_k1] > max_ROW)
					{
						max_ROW = row_max[prev_k1];
						k1 = prev_k1;
						row_ties.clear();
					}
					else if (row_max[prev_k1] == max_ROW)
					{
						if (row_ties.empty())
						{
							row_ties.push_back(k1);
						}
						row_ties.push_back(prev_k1);
					}

					if (col_max[prev_k1] > max_COL)
					{
						max_COL = col_max[prev_k1];
						k2 = prev_k1;
						col_ties.clear();
					}
					else if (col_max[prev_k1] == max_COL)
					{
						if (col_ties.empty())
						{
							col_ties.push_back(k2);
						}
						col_ties.push_back(prev_k1);
					}
	
					if (!row_ties.empty() || !col_ties.empty())
					{
						if (row_ties.empty())
						{
							row_ties.push_back(k1);
						}
						else if (col_ties.empty())
						{
							col_ties.push_back(k2);
						}

						double max_RO = 0;
						for(unsigned m=0;m<row_ties.size();++m)
						{
							for(unsigned n=0;n<col_ties.size();++n)
							{
								if (row_ties[m]<col_ties[n])
								{
									double r = RO(cluster[row_ties[m]], cluster[col_ties[n]]);
									if (r>max_RO)
									{
										max_RO = r;
										k1 = row_ties[m];
										k2 = col_ties[n];
									}
								}
							}
						}
					}
					max_RO = max_ROW;
				}
			}

			intervals.insert(intervals.end(), cluster.begin(), cluster.end());
			k0 = k;
			
			curr.first = events[k].first;
			curr.second = events[k].second;
		}
		else if (curr.second < events[k].second)
		{
			curr.second = events[k].second;
		}
	}
}
*/

