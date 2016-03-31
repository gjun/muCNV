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

extern double P_THRESHOLD;
extern double BE_THRESHOLD;

void call_duplications(vector<vector<double> > &X, vector<double> &AvgDepth, vector<string> &sampleIDs, vector<Interval> &intervals, FILE* vcfFile)
{
	unsigned n_interval = intervals.size();
	unsigned n_sample = sampleIDs.size();

	vector<unsigned short> GT(n_sample,0);
	vector<unsigned> GQ(n_sample, 0);

	// For each candidate region, run EM
	// Run EM with 2 and 3 components, compare likelihoods with the 1-gaussian model, apply BIC

	unsigned write_cnt = 0;
	for(unsigned k=0; k<n_interval;++k)
	{
		vector<Gaussian> Comps(1);

		double min_BIC = DBL_MAX;

		double M = mean(X[k]);
		double S = stdev(X[k], M);

		Comps[0].Mean = M;
		Comps[0].Stdev = S; 
		Comps[0].Alpha = 1;

		min_BIC = BIC(X[k], Comps);
//		cerr << "1 comp, bic " << min_BIC << endl;
		unsigned n_comp = 1;

		// Check up to 30 -- should be enough? -- let's go for 10 now
		// It seems like too many components make probs too small hence failing

		for(unsigned i=2; i<4; ++i)
		{
			Gaussian* pG = new Gaussian;
			Comps.push_back(*pG);

			for(unsigned j=0; j<i; ++j)
			{
				Comps[j].Mean = 1.0 + (double)j*0.5;
				Comps[j].Stdev = 0.1;
				Comps[j].Alpha = 1;
			}
			// Run EM with i components
			EM(X[k], Comps);

			double bic = BIC(X[k], Comps);

			if (bic<min_BIC)
			{
				min_BIC = bic;
				n_comp = i;
			}
		}

		Comps.resize(n_comp);

		// Determine whether this is a polymorphic interval or not, make calls

//		cerr << "n_comp: " << n_comp  << endl;

		for(unsigned j=0;j<n_comp;++j)
		{
			Comps[j].Mean = 1.0 + (double)j*0.5;
			Comps[j].Stdev = 0.1;
			Comps[j].Alpha = 1;
		}
		EM(X[k], Comps);
		double BE = BayesError(Comps);
		if (BE<BE_THRESHOLD)
		{
			unsigned NS=0;
			unsigned AC = classify(X[k], GT, GQ, NS, Comps);
			if (AC>0)
			{
				++write_cnt;
				write_vcf_dup(vcfFile, write_cnt, GT, GQ, intervals[k], AC, NS, X[k], AvgDepth, Comps, true);
			}
		}
		else
		{
			unsigned NS=0;
			unsigned AC = classify(X[k], GT, GQ, NS, Comps);
			if (AC>0)
			{
				++write_cnt;
				write_vcf_dup(vcfFile, write_cnt, GT, GQ, intervals[k], AC, NS, X[k], AvgDepth, Comps, false);
			}
	
		}
	}
}

/*
void call_duplications(vector<vector<double> > &X, vector<double> &AvgDepth, vector<string> &sampleIDs, vector<Interval> &intervals, FILE* vcfFile)
{
	unsigned n_interval = intervals.size();
	unsigned n_sample = sampleIDs.size();

	vector<unsigned short> GT(n_sample,0);
	vector<unsigned> GQ(n_sample, 0);

	// For each candidate region, run EM
	// Run EM with 2 and 3 components, compare likelihoods with the 1-gaussian model, apply BIC

	unsigned write_cnt = 0;
	for(unsigned k=0; k<n_interval;++k)
	{
		vector<Gaussian> Comps(1);
		double min_BIC = DBL_MAX;
		double min_ERR = DBL_MAX;
		double M = mean(X[k]);
		double S = stdev(X[k], M);

		Comps[0].Mean = M;
		Comps[0].Stdev = S; 
		Comps[0].Alpha = 1;

		min_BIC = BIC(X[k], Comps);
//		cerr << "1 comp, bic " << min_BIC << endl;
		unsigned n_comp = 1;

		// Check up to 30 -- should be enough? -- let's go for 10 now
		// It seems like too many components make probs too small hence failing

		for(unsigned i=2; i<4; ++i)
		{
			Gaussian* pG = new Gaussian;
			Comps.push_back(*pG);
			double bic = DBL_MAX;
			double err = DBL_MAX;
			for(unsigned m=0;m<20;++m)
			{
				for(unsigned j=0; j<i; ++j)
				{
					// Randomly pick initial means ~10 times -- will take some time but worth it
					int pick = rand()%n_sample;
					Comps[j].Mean = X[k][pick];
					Comps[j].Stdev = 0.1;
					Comps[j].Alpha = 1;
				}
				// Run EM with i components
				EM(X[k], Comps);

				double val = BIC(X[k], Comps);

				if (val<bic)
				{
					bic = val;
					err = BayesError(Comps);
				}
			}

//			cerr << i << " comps, bic " << bic << ", BE = " << BayesError(Comps) <<  endl;

//			if (bic<min_BIC && err< 0.1)
			if (err < min_ERR)
			{
				min_BIC = bic;
				min_ERR = err;
				n_comp = i;
			}
		}
		if (min_ERR>0.1)
		{
			n_comp = 1;
		}
		Comps.resize(n_comp);

		// Determine whether this is a polymorphic interval or not, make calls

//		cerr << "n_comp: " << n_comp  << endl;
		// Two-component is selected
		if (n_comp>1)
		{
			for(unsigned m=0; m<20; ++m)
			{
				for(unsigned j=0;j<n_comp;++j)
				{
					int pick = rand()%n_sample;
					Comps[j].Mean = X[k][pick];
					Comps[j].Stdev = 0.1;
					Comps[j].Alpha = 1;
				}
				EM(X[k], Comps);
			}
			double BE = BayesError(Comps);
			if (BE<0.1)
			{
				unsigned NS=0;
				unsigned AC = classify(X[k], GT, GQ, NS, Comps);


			//	if (AC>0)
				{
					++write_cnt;
					write_vcf_dup(vcfFile, write_cnt, GT, GQ, intervals[k], AC, NS, X[k], AvgDepth, Comps);
				}
			}
		}
		
	}
}
*/



void EM(vector<double>& x, vector<Gaussian>& Comps)
{
	unsigned n_sample = x.size();
	unsigned n_comp = Comps.size();
	unsigned n_iter = 30;

	for(unsigned i=0; i<n_iter; ++i)
	{
		vector<double> sum (n_comp,0);
		vector<double> sum_pr (n_comp,0);
		vector<double> sum_err (n_comp,0);

		// E step
		for(unsigned j=0; j<n_sample; ++j)
		{
			double sum_p = 0;
			vector<double> pr(n_comp, 0);
			for(unsigned m=0;m<n_comp;++m)
			{
				pr[m] = Comps[m].Alpha * normpdf(x[j], Comps[m]);
				sum_p += pr[m];
			}

			if (sum_p > 1e-100) // if the value is an outlier, exclude it from calculations
			{
				for(unsigned m=0;m<n_comp;++m)
				{
					pr[m] /= sum_p;
					sum[m] += pr[m] * x[j];
					sum_err[m] += pr[m] * (x[j] - Comps[m].Mean)*(x[j]-Comps[m].Mean);
					sum_pr[m] += pr[m];
				}
			}
		}

		// M step
		for(unsigned m=0;m<n_comp;++m)
		{
			Comps[m].Mean = sum[m]/sum_pr[m];
			Comps[m].Stdev = sqrt(sum_err[m] / (sum_pr[m] ));
			Comps[m].Alpha = sum_pr[m] / n_sample;
		}
	}
}

unsigned classify(vector<double>& x, vector<unsigned short>& GT, vector<unsigned>& GQ, unsigned& ns, vector<Gaussian>& Comps)
{
	unsigned n_sample = x.size();
	unsigned n_comp = Comps.size();
	unsigned ac = 0;
	ns = 0;

	for(unsigned j=0; j<n_sample; ++j)
	{
		double max_pos = -1.0 * DBL_MAX;
		double max_lk = -1.0 * DBL_MAX;
		double lk[n_comp];
		double pos[n_comp];
		double sum_pos = 0;
		double sum_lk = 0;
		unsigned short gt=0;

		for(unsigned m=0; m<n_comp; ++m)
		{
			//posterior
			lk[m]  = normpdf(x[j], Comps[m] );
			pos[m] = Comps[m].Alpha * lk[m];

			sum_pos += pos[m];
			sum_lk += lk[m];

			if (lk[m]>max_lk)
			{
				max_lk = lk[m];
			}
			if (pos[m]>max_pos)
			{
				gt = round(Comps[m].Mean*2.0);
				max_pos = pos[m];
			}
		}

		if (gt<2 || gt>4)
		{
			// non-biallelic, let's handle this later
			return 0;
		}

		double max_p = max_pos;
		double sum_p = sum_pos;

		if (max_p / sum_p > P_THRESHOLD && sum_p > 1e-100)
		{
			GT[j] = gt;
			if (gt>2)
			{
				ac += (gt - 2);
				++ns;
			}
			double p_err = (sum_p-max_p) / sum_p;

			if (p_err > 8.4e-12)
			{
				GQ[j] = (unsigned)floor(-10.0 * log10((sum_p - max_p) / (sum_p) ));

				if (GQ[j]>255)
				{
					GQ[j] = 255;
					//cerr << "error: max_p is " << max_p << " and sum_p is " << sum_p << ", GQ is " << GQ[j] << endl;
				}
			}
			else
			{
				GQ[j] = 255;
			}
		}
		else
		{
			GT[j] = 0;
			if (sum_p>1e-100)
			{
				double p_err = (sum_p-max_p) / sum_p;

				if (p_err > 8.4e-12)
				{
					GQ[j] = (unsigned)floor(-10.0 * log10 ((sum_p - max_p) / (sum_p) ));

					if (GQ[j]>255)
					{
						GQ[j] = 255;
				//		cerr << "error: max_p is " << max_p << " and sum_p is " << sum_p << ", GQ is " << GQ[j] << endl;
					}
				}
				else
				{
					GQ[j] = 255;
				}
			}
			else
			{
				GQ[j] = 0;
			}
		}
	}
	return ac;
}

/*
unsigned classify(vector<double>& x, vector<unsigned short>& GT, vector<unsigned>& GQ, unsigned& ns, vector<Gaussian>& Comps)
{
	unsigned n_sample = x.size();
	unsigned n_comp = Comps.size();
	unsigned ac = 0;
	ns = 0;

	for(unsigned j=0; j<n_sample; ++j)
	{
		double max_pos = -1.0 * DBL_MAX;
		double max_lk = -1.0 * DBL_MAX;
		double lk[n_comp];
		double pos[n_comp];
		double sum_pos = 0;
		double sum_lk = 0;
		unsigned short gt=0;

		for(unsigned m=0; m<n_comp; ++m)
		{
			//posterior
			lk[m]  = normpdf(x[j], Comps[m] );
			pos[m] = Comps[m].Alpha * lk[m];

			sum_pos += pos[m];
			sum_lk += lk[m];

			if (lk[m]>max_lk)
			{
				max_lk = lk[m];
			}
			if (pos[m]>max_pos)
			{
				gt = round(Comps[m].Mean*2.0)+1;
				max_pos = pos[m];
			}
		}

		double max_p = max_pos;
		double sum_p = sum_pos;

		if (max_p / sum_p > P_THRESHOLD && sum_p > 1e-100)
		{
			GT[j] = gt;
			if (gt!=0)
			{
				++ac;
				++ns;
			}
			double p_err = (sum_p-max_p) / sum_p;

			if (p_err > 8.4e-12)
			{
				GQ[j] = (unsigned)floor(-10.0 * log10((sum_p - max_p) / (sum_p) ));

				if (GQ[j]>255)
				{
					GQ[j] = 255;
					//cerr << "error: max_p is " << max_p << " and sum_p is " << sum_p << ", GQ is " << GQ[j] << endl;
				}
			}
			else
			{
				GQ[j] = 255;
			}
		}
		else
		{
			GT[j] = 0;
			if (sum_p>1e-100)
			{
				double p_err = (sum_p-max_p) / sum_p;

				if (p_err > 8.4e-12)
				{
					GQ[j] = (unsigned)floor(-10.0 * log10 ((sum_p - max_p) / (sum_p) ));

					if (GQ[j]>255)
					{
						GQ[j] = 255;
				//		cerr << "error: max_p is " << max_p << " and sum_p is " << sum_p << ", GQ is " << GQ[j] << endl;
					}
				}
				else
				{
					GQ[j] = 255;
				}
			}
			else
			{
				GQ[j] = 0;
			}
		}
	}
	return ac;
}
*/
