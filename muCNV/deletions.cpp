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

extern bool bUseGL;
extern double P_THRESHOLD;
extern double BE_THRESHOLD;


bool ordered(vector<Gaussian> &C)
{
	for(unsigned i=0; i<C.size()-1;++i)
	{
		if (C[i].Mean<=C[i+1].Mean)
			return false;
	}
	return true;
}

#ifdef NEVER //TMPTMP
void call_deletion(Interval &interval, vector<double> X, vector<double> &AvgDepth, FILE* vcfFile)
{
	unsigned n_sample = X.size();

	unsigned short GT;
	vector<unsigned> GL(3,255);
	unsigned GQ;
	unsigned write_cnt = 0; // TMPTMP

	vector<Gaussian> C1(1); // 1-component model
	vector<Gaussian> C2(2); // 2-component model
	vector<Gaussian> C3(3); // 3-component model

	// For each candidate region, run EM
	// Run EM with 2 and 3 components, compare likelihoods with the 1-gaussian model, apply BIC
	bool bFlip = false;

	double BIC1, BIC2, BIC3;

	// Single-component model
	C1[0].Mean = mean(X);
	C1[0].Stdev = stdev(X, C1[0].Mean);
	C1[0].Alpha = 1;

	BIC1 = BIC(X[k], C1);

	// Two-component model
	C2[0].Mean = 1;
	C2[1].Mean = 0.5;
	C2[0].Stdev = C2[1].Stdev = 0.1;
	C2[0].Alpha = C2[1].Alpha = 0.5;

	// Three-component model
	C3[0].Mean = 1;
	C3[1].Mean = 0.5;
	C3[2].Mean = 0;
	C3[0].Stdev = C3[1].Stdev = C3[2].Stdev = 0.1;
	C3[0].Alpha = C3[1].Alpha = C3[2].Alpha = 1.0/3.0;

	EM(X, C3, bFlip);
	BIC3 = BIC(X, C3);

	// check which one is major
	if (C3[2].Mean < C3[0].Mean && C3[2].Alpha > C3[0].Alpha)
	{

		bFlip = true;
		C2[0].Mean = 0.5;
		C2[1].Mean = 0;

		EM(X, C2, bFlip);
	}
	else
	{
		// Two-component EM and BIC
		EM(X, C2, bFlip);
	}
	BIC2 = BIC(X, C2);

	// Determine whether this is a polymorphic interval or not, make calls
	if (BIC1>BIC2 && BIC3>BIC2)
	{
		// 2 cluster
		double BE = BayesError(C2);
		if (ordered(C2) && BE< BE_THRESHOLD)
		{
			unsigned NS=0;
			unsigned AC = classify(X, GT, GL, GQ, NS, C2, bFlip);

	//		cerr << "AC: " << AC << endl;
			if (AC>0)
			{
				++write_cnt;
				write_vcf(vcfFile, write_cnt, GT, GL, GQ, interval, AC, NS, X, AvgDepth, C2, BE, true);
			}
		}
		else if (ordered(C2))
		{
			unsigned NS=0;
			unsigned AC = classify(X[k], GT, GL, GQ, NS, C2, bFlip);

	//		cerr << "AC: " << AC << endl;
			if (AC>0)
			{
				++write_cnt;
				write_vcf(vcfFile, write_cnt, GT, GL, GQ, interval, AC, NS, X, AvgDepth, C2, BE, false);
			}

		}
		
		// Two-component is selected
	//	printCluster(C2);
	}
	else if (BIC3<BIC1 && BIC3<BIC2)
	{
		// 3 cluster
		double BE = BayesError(C3);
		if (ordered(C3) && BE< BE_THRESHOLD)
		{
			unsigned NS=0;
			unsigned AC = classify(X, GT, GL, GQ, NS, C3, bFlip);

			if (AC>0)
			{
				++write_cnt;
				write_vcf(vcfFile, write_cnt, GT, GL, GQ, interval, AC, NS, X, AvgDepth, C3, BE, true);
			}
		}
		else if (ordered(C3))
		{
			unsigned NS=0;
			unsigned AC = classify(X[k], GT, GL, GQ, NS, C3, bFlip);

			if (AC>0)
			{
				++write_cnt;
				write_vcf(vcfFile, write_cnt, GT, GL, GQ, interval, AC, NS, X, AvgDepth, C3, BE, false);
			}
		}

	//	printCluster(C3);
	}
}
#endif //TMPTMP

void call_deletions(vector<vector<double> > &X, vector<double> &AvgDepth, vector<string> &sampleIDs, vector<Interval> &intervals, FILE* vcfFile)
{
	unsigned n_interval = intervals.size();
	unsigned n_sample = sampleIDs.size();

	vector<unsigned short> GT(n_sample,0);
	vector< vector<unsigned> > GL(n_sample, vector<unsigned>(3,255));
	vector<unsigned> GQ(n_sample, 0);

	vector<Gaussian> C1(1); // 1-component model
	vector<Gaussian> C2(2); // 2-component model
	vector<Gaussian> C3(3); // 3-component model

	// For each candidate region, run EM
	// Run EM with 2 and 3 components, compare likelihoods with the 1-gaussian model, apply BIC

	unsigned write_cnt = 0;
	for(unsigned k=0; k<n_interval;++k)
	{
		bool bFlip = false;

		double BIC1, BIC2, BIC3;

		// Single-component model
		C1[0].Mean = mean(X[k]);
		C1[0].Stdev = stdev(X[k], C1[0].Mean);
		C1[0].Alpha = 1;

		BIC1 = BIC(X[k], C1);

		// Two-component model
		C2[0].Mean = 1;
		C2[1].Mean = 0.5;
		C2[0].Stdev = C2[1].Stdev = 0.1;
		C2[0].Alpha = C2[1].Alpha = 0.5;

		// Three-component model
		C3[0].Mean = 1;
		C3[1].Mean = 0.5;
		C3[2].Mean = 0;
		C3[0].Stdev = C3[1].Stdev = C3[2].Stdev = 0.1;
		C3[0].Alpha = C3[1].Alpha = C3[2].Alpha = 1.0/3.0;

		EM(X[k], C3, bFlip);
		BIC3 = BIC(X[k], C3);

		// check which one is major
		if (C3[2].Mean < C3[0].Mean && C3[2].Alpha > C3[0].Alpha)
		{

			bFlip = true;
			C2[0].Mean = 0.5;
			C2[1].Mean = 0;

			EM(X[k], C2, bFlip);
		}
		else
		{
			// Two-component EM and BIC
			EM(X[k], C2, bFlip);
		}
		BIC2 = BIC(X[k], C2);

		// Determine whether this is a polymorphic interval or not, make calls
		if (BIC1>BIC2 && BIC3>BIC2)
		{
			// 2 cluster
			double BE = BayesError(C2);
			if (ordered(C2) && BE< BE_THRESHOLD)
			{
				unsigned NS=0;
				unsigned AC = classify(X[k], GT, GL, GQ, NS, C2, bFlip);

		//		cerr << "AC: " << AC << endl;
				if (AC>0)
				{
					++write_cnt;
					write_vcf(vcfFile, write_cnt, GT, GL, GQ, intervals[k], AC, NS, X[k], AvgDepth, C2, BE, true);
				}
			}
			else if (ordered(C2))
			{
				unsigned NS=0;
				unsigned AC = classify(X[k], GT, GL, GQ, NS, C2, bFlip);

		//		cerr << "AC: " << AC << endl;
				if (AC>0)
				{
					++write_cnt;
					write_vcf(vcfFile, write_cnt, GT, GL, GQ, intervals[k], AC, NS, X[k], AvgDepth, C2, BE, false);
				}

			}
			
			// Two-component is selected
		//	printCluster(C2);
		}
		else if (BIC3<BIC1 && BIC3<BIC2)
		{
			// 3 cluster
			double BE = BayesError(C3);
			if (ordered(C3) && BE< BE_THRESHOLD)
			{
				unsigned NS=0;
				unsigned AC = classify(X[k], GT, GL, GQ, NS, C3, bFlip);

				if (AC>0)
				{
					++write_cnt;
					write_vcf(vcfFile, write_cnt, GT, GL, GQ, intervals[k], AC, NS, X[k], AvgDepth, C3, BE, true);
				}
			}
			else if (ordered(C3))
			{
				unsigned NS=0;
				unsigned AC = classify(X[k], GT, GL, GQ, NS, C3, bFlip);

				if (AC>0)
				{
					++write_cnt;
					write_vcf(vcfFile, write_cnt, GT, GL, GQ, intervals[k], AC, NS, X[k], AvgDepth, C3, BE, false);
				}
			}
	
		//	printCluster(C3);
		}
	
	}
}
void EM(vector<double>& x, vector<Gaussian>& C, bool bFlip)
{
	unsigned n_sample = x.size();
	unsigned n_comp = C.size();

	unsigned n_iter = 30;
	// pseudo-counts
	unsigned p_count= 2;
	// pseudo-means
	double p_val[3] = {1.0, 0.5, 0};

	if (bFlip && n_comp==2) 
	{
		// if fitting for (0, 1) copies, set pseudo-means accordingly
		p_val[0] = 0.5;
		p_val[1] = 0;
	}

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
				pr[m] = C[m].Alpha * normpdf(x[j], C[m]);
				if (m==2)
				{
					// half-normal scaler
					pr[m] *= 2.0;
				}
				else if (bFlip && n_comp==2 &&  m==1)
				{
					pr[m] *=2.0;
				}
				sum_p += pr[m];
			}

			if (sum_p > 1e-100) // if the value is an outlier, exclude it from calculations
			{
				for(unsigned m=0;m<n_comp;++m)
				{
					pr[m] /= sum_p;
					sum[m] += pr[m] * x[j];
					sum_err[m] += pr[m] * (x[j] - C[m].Mean)*(x[j]-C[m].Mean);
					sum_pr[m] += pr[m];
				}
			}
		}

		// Add pseudo-count values
		for(unsigned m=0;m<n_comp; ++m)
		{
			sum[m] += p_val[m] * p_count;
			sum_err[m] += (p_val[m] - C[m].Mean)*(p_val[m] - C[m].Mean) * p_count;
			sum_pr[m] += p_count;
		}

		// M step
		for(unsigned m=0;m<n_comp;++m)
		{
			C[m].Mean = sum[m]/sum_pr[m];
			C[m].Stdev = sqrt(sum_err[m] / (sum_pr[m] ));
			C[m].Alpha = sum_pr[m] / (n_sample + p_count*n_comp);
			if (m==2)
			{
				// half-normal with mean 0 and sigma
				C[m].Mean = 0;
			}
			else if (bFlip && n_comp==2 && m==1)
			{
				C[m].Mean = 0;
			}

		}
	}
}

void conEM(vector<double>& x, vector<double>& M, vector<double>& S, vector<double>& A)
{
	unsigned n_sample = x.size();
	unsigned n_comp = M.size();

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
				pr[m] = A[m] * normpdf(x[j], M[m], S[m]);
				if (m==2)
				{
					// half-normal scaler
					pr[m] *= 2.0;
				}
				sum_p += pr[m];
			}

			if (sum_p > 1e-100) // if the value is an outlier, exclude it from calculations
			{
				for(unsigned m=0;m<n_comp;++m)
				{
					pr[m] /= sum_p;
					sum[m] += pr[m] * x[j];
					sum_err[m] += pr[m] * (x[j] - M[m])*(x[j]-M[m]);
					sum_pr[m] += pr[m];
				}
			}
		}

		// M step

		M[0] = (sum[0] + sum[1]*2.0) / (sum_pr[0] + sum_pr[1]);
		S[0] = sqrt(sum_err[0] / sum_pr[0]);
		A[0] = sum_pr[0] / (n_sample);


		M[1] = M[0]/2.0;
		S[1] = sqrt(sum_err[1] / sum_pr[1]);
		A[1] = sum_pr[1] / (n_sample);

		if (n_comp>2)
		{
			// half-normal with mean 0 and sigma
			M[2] = 0;
			S[2] = sqrt(sum_err[2] / sum_pr[2]);
			A[2] = sum_pr[2] / (n_sample);
		}


	}
}

unsigned classify(vector<double>& x, vector<unsigned short>& GT, vector< vector<unsigned> >& GL, vector<unsigned>& GQ, unsigned& ns, vector<Gaussian>& C, bool bFlip)
{
	unsigned n_sample = x.size();
	unsigned n_comp = C.size();
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
			lk[m]  = normpdf(x[j], C[m]);
			pos[m] = C[m].Alpha * lk[m];

			sum_pos += pos[m];
			sum_lk += lk[m];

			if (lk[m]>max_lk)
			{
				if (bUseGL) gt = m+1;
				max_lk = lk[m];
			}
			if (pos[m]>max_pos)
			{
				if (!bUseGL) gt = m+1;
				max_pos = pos[m];
			}
		}

		if (bFlip && n_comp==2)
		{
			++gt;
		}
		for(unsigned m=0; m<n_comp; ++m)
		{
			if (sum_lk>1e-100)
			{
				if (lk[m] < 8.4e-12)
				{
					GL[j][m] = 255;
				}
				else
				{
					GL[j][m] = -10.0*log10(lk[m]/sum_lk);
					if (GL[j][m]>255)
					{
						cerr << "error: lk is " << lk[m] << " and sum_lk is " << sum_lk << ", GL is " << GL[j][m] << endl;
					}
				}
			}
			else
			{
				GL[j][m] = 255;
			}
		}

		if (n_comp==2)
		{
			if (bFlip)
			{
				GL[j][0] = 255;
			}
			else
			{
				GL[j][2] = 255;
			}
		}

		double max_p = (bUseGL) ? max_lk : max_pos;
		double sum_p = (bUseGL) ? sum_lk : sum_pos;

		if (max_p / sum_p > P_THRESHOLD && sum_p > 1e-100)
		{
			GT[j] = gt;
			if (gt>0)
			{
				ac += (gt-1);
				++ns;
			}
			double p_err = (sum_p-max_p) / sum_p;

			if (p_err > 8.4e-12)
			{
				GQ[j] = (unsigned)floor(-10.0 * log10((sum_p - max_p) / (sum_p) ));

				if (GQ[j]>255)
				{
					cerr << "error: max_p is " << max_p << " and sum_p is " << sum_p << ", GQ is " << GQ[j] << endl;
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
						cerr << "error: max_p is " << max_p << " and sum_p is " << sum_p << ", GQ is " << GQ[j] << endl;
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
