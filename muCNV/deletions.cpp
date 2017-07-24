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
