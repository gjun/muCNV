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

extern double BE_THRESHOLD;
extern double P_THRESHOLD;

int gtype::assign(double x, vector<Gaussian> &C)
{
	int n_comp = (int) C.size();
	double p[n_comp];
	double max_P = -1;
	double max_R = -1;
	int ret = -1;

	for(int i=0;i<n_comp; ++i)
	{
// Major allele dominates when sample size is large 
//		p[i] = C[i].Alpha * C[i].pdf(x);
		p[i] = C[i].pdf(x);
		
		if (p[i] > max_P)
		{
			max_P = p[i];
			ret = i;
		}
	}
	for(int i=0;i<n_comp; ++i)
	{
		if (ret != i)
		{
			double R = p[i] / max_P;
			if (R>max_R)
			{
				max_R = R;
			}
		}
	}
	
	ret = round(C[ret].Mean * 2);
	
	if (max_R > P_THRESHOLD)
	{
		return -1;
	}
	
	return ret;
}


void svgeno::print(sv &S, svdata &D, string &ln)
{
	ln = S.chr;
	ln += "\t" + to_string(S.pos) + "\t" + S.svtype + "_" + S.chr + ":" + to_string(S.pos) + "-" + to_string(S.end) + "\t.\t<" + S.svtype + ">\t.\t";

	string info = "SVTYPE=" + S.svtype + ";END=" + to_string(S.end) + ";SVLEN=" + to_string(S.len) + ";NCLUS=" + to_string(C.size()) + ";P_OVERLAP=" + to_string(BE_dp);


	if (S.svtype == "INV")
	{
		if (cnv_pos_flag)
		{
			info += ";RF_POS";
		}
		if (cnv_neg_flag)
		{
			info += ";RF_NEG";
		}
		if (inv_pos_flag)
		{
			info += ";FF";
		}
		if (inv_neg_flag)
		{
			info += ";RR";
		}
	}
	else
	{
		if (dp_flag)
		{
			info += ";DP";
		}
		if (pos_flag)
		{
			info += ";POS";
		}
		if (neg_flag)
		{
			info += ";NEG";
		}
		info += ";MEAN=" + to_string(Comps[0].Mean);
		for(unsigned i=1;i<Comps.size(); ++i)
		{
			info += "," + to_string(Comps[i].Mean);
		}

		info += ";STDEV=" + to_string(Comps[0].Stdev);
		for(unsigned i=1;i<Comps.size(); ++i)
		{
			info += "," + to_string(Comps[i].Stdev);
		}

		info += ";FRAC=" + to_string(Comps[0].Alpha);
		for(unsigned i=1;i<Comps.size(); ++i)
		{
			info += "," + to_string(Comps[i].Alpha);
		}
	}

	info +=";AC=" + to_string(ac) + ";NS=" + to_string(ns) + ";AF=" + to_string((double)ac/(double)(2.0*ns));   
	
	if (b_biallelic)
	{
		if (S.svtype == "INV")
		{
			info += "\tGT:PR:DP:FP:FN:IP:IN";
		}
		else
		{
			info += "\tGT:CN:ND:DP:FP:FN";
		}
	}

	string filt = "FAIL";
	if (dp_flag || pos_flag || neg_flag)
	{
		if ((ac>0 &&  ac < (2*n_sample) && ns>(n_sample *0.5)) && ((dp_flag && p_overlap<0.2) ||  (p_overlap < 0.3 && (neg_flag || pos_flag)) || (p_overlap < 0.4 && neg_flag && pos_flag)))
		{
			b_pass = true;
			filt = "PASS";
		}
	}
	else
	{
		filt = "FAIL;DUMP;
	}
	
	ln += filt + "\t" + info;
	
	for(int i=0;i<n_sample;++i)
	{
		switch(gt[i])
		{
			case 0:
				ln += "\t0/0:2";
				break;
			case 1:
				ln += "\t0/1:1";
				break;
			case 2:
				ln += "\t1/1:0";
				break;
			default:
				ln += "\t.:.";
				break;
		}
		ln += ":" + to_string(D.norm_dp[i]).substr(0,4) + ":" + to_string((int)D.dp[i]) + ":" + to_string((int)D.n_cnv_pos[i]) + "," + to_string((int)D.cnv_pos[i]) + ":" + to_string((int)D.n_cnv_neg[i]) + "," + to_string((int)D.cnv_neg[i]);
	}
}

void gtype::copyComps(vector<Gaussian> &C, vector<Gaussian> &C0)
{
	C.clear();
	C.resize(C0.size());
	for(unsigned j=0;j<C.size(); ++j)
	{
		C[j].set(C0[j].Mean, C0[j].Stdev);
		C[j].Alpha = C0[j].Alpha;
	}
}


void gtype::call_del(sv &S, svdata& dt, svgeno &G)
{
	int n_sample=dt.n;
	double bic1, bic2, bic2_rev, bic3;

	vector<Gaussian> C1(1); // 1-component model
	vector<Gaussian> C2(2); // 2-component model
	vector<Gaussian> C2_rev(2); // 2-component model, with reversed
	vector<Gaussian> C3(3); // 3-component model
	
	// For each candidate region, run EM
	// Fit Gaussian mixture models with 1, 2, and 3 components, compare BIC
	C1[0].estimate(dt.norm_dp);
	C1[0].Alpha = 1;
	
	double min_bic = DBL_MAX;

	// P_Overlap
	G.P_overlap = 1;
	G.b_biallelic = true; //deletions 
	
	// One-component model
	bic1 = BIC(dt.norm_dp, C1);
	min_bic = bic1;
	copyComps(G.Comps, C1);
	
	// Two-component model
	C2[0].set(1, 0.1);
	C2[1].set(0.5, 0.1);
	C2[0].Alpha = C2[1].Alpha = 0.5;
	EM(dt.norm_dp, C2);
	bic2 = BIC(dt.norm_dp, C2);

//	if (bic2 < min_bic && C2[0].Mean > C2[1].Mean && C2[0].Mean > 0.75 && C2[1].Mean > 0.3 && C2[1].Mean < 0.75)
	if (bic2 < min_bic) // now we do not compare mean
	{
		min_bic = bic2;
		dp_flag = true;
		BE_dp = BayesError(C2);
		copyComps(G.Comps, C2);
	}

	// Three-component model
	C3[0].set(1, 0.1);
	C3[1].set(0.5, 0.1);
	C3[2].set(0, 0.1);
	C3[0].Alpha = C3[1].Alpha = C3[2].Alpha = 1.0/3.0;
	EM(dt.norm_dp, C3);
	bic3 = BIC(dt.norm_dp, C3);

//	if (bic3 < min_bic &&  C3[0].Mean > 0.75 && C3[0].Mean > (C3[1].Mean+0.35) && C3[1].Mean > (C3[2].Mean+0.35))
	if (bic3 < min_bic)
	{
		min_bic = bic3;
		dp_flag = true;
		BE_dp = BayesError(C3);
		copyComps(G.Comps, C3);
	}

	vector<int> pos_gt(n_sample, 0);
	vector<int> neg_gt(n_sample, 0);

	if (S.len > 50) // 100?
	{
		for(int i=0;i<n_sample;++i)
		{
			if (dt.norm_cnv_pos[i] > 0.5 && dt.norm_cnv_pos[i] < 2 && dt.n_cnv_pos[i] > 2)
			{
				pos_gt[i] = 1;
				G.pos_flag = true;
			}
			if (dt.norm_cnv_neg[i] > 0.5 && dt.norm_cnv_neg[i] < 2 && dt.n_cnv_neg[i] > 2)
			{
				neg_gt[i] = 1;
				G.neg_flag = true;
			}
		}
	}
		
//	if (G.dp_flag || G.pos_flag || G.neg_flag)
	{
		if (!G.dp_flag)
		{
			// Depth-based clustering failed, force 2 or 3 component model
			if (bic2 < bic3)
			{
				copyComps(G.Comps, C2);
			}
			else
			{
				copyComps(G.Comps, C3);
			}
			G.p_overlap = BayesError(C);
		}

		for(int i=0;i<n_sample;++i)
		{
			G.cn[i] = assign(X[i], C);

			// To prevent high depth ones are classified as double deletion due to long tail
			if ((dt.norm_dp[i]*2 - cn) > 1 || (dt.norm_dp[i]*2 - cn) < -1)
			{
				G.cn[i] = -1;
			}
			switch(G.cn[i])
			{
				case -1:
					if (pos_gt[i] && neg_gt[i])
					{
						if ( 0.2 <= dt.norm_dp[i] && dt.norm_dp[i]< 1)
						{
							G.gt[i] = 1;
							G.cn[i] = 1;
						}
						else if (dt.norm_dp[i] < 0.2)
						{
							G.gt[i] = 2;
							G.cn[i] = 0;
						}
					}
					else if (pos_gt[i] || neg_gt[i])
					{
						if ( 0.3 <= dt.norm_dp[i] && dt.norm_dp[i] <= 0.7)
						{
							G.gt[i] = 1;
							G.cn[i] = 1;
						}
						else if (dt.norm_dp[i] < 0.2)
						{
							G.gt[i] = 2;
							G.cn[i] = 0;
						}
					}
					break;
				case 0:
					G.gt[i] = 2;
					G.cn[i] = 0;
					break;
				case 1:
					G.gt[i] = 1;
					G.cn[i] = 1;
					break;
				case 2:
					G.gt[i] = 0;
					G.cn[i] = 2;
					if (pos_gt[i] && neg_gt[i] && dt.norm_dp[i] < 1)
					{
						G.gt[i] = 1;
						G.cn[i] = 1;
					}
					else if ((pos_gt[i] || neg_gt[i]) && dt.norm_dp[i] <=0.7 )
					{
						G.gt[i] = 1;
						G.cn[i] = 1;
					}
					break;
			}
			if (gt[i] >=0 )
			{
				G.ac += gt[i];
				G.ns += 1;
			}
		}
	}
}

void gtype::call_cnv(sv &S, svdata& D, svgeno& G)
{
	int n_sample = D.n;

	vector<Gaussian> C(1);
	double min_BIC = DBL_MAX;

	unsigned n_comp = 1;
	
	vector<int> cnt(10,0);
	// Count CN for each bin
	for(int i=0;i<n_sample;++i)
	{
		int idx = round(D.norm_dp[i]*2.0);
		if (idx>9) idx = 9;

		cnt[idx]++;
	}

	sort(cnt.begin(), cnt.end(), std::greater<int>());

	C[0].Mean = cnt[0]/2.0;
	C[0].Stdev = stdev(D.norm_dp, C[0].Mean);
	min_BIC = BIC(D.norm_dp, C);
	copyComps(G.Comps, C);

	for(int i=1; i<10 && cnt[i]>0; ++i)
	{
		Gaussian *pG = new Gaussian;
		C.push_back(*pG);

		for(unsigned j=0; j<i; ++j)
		{
			C[j].Mean = (double)cnt[j]*0.5;
			C[j].Stdev = 0.1;
			C[j].Alpha = 1.0/(i+1);
		}
		EM(D.norm_dp, C);

		double bic = BIC(D.norm_dp, C); 
		if (bic < min_BIC)
		{
			n_comp = i+1;
			min_BIC = bic;
			copyComps(G.Comps, C);
			G.dp_flag = true;
		}
	}
	G.p_overlap = BayesError(G.Comps);

	if (S.len > 50) // 100?
	{
		for(int i=0;i<n_sample;++i)
		{
			if (round(D.norm_dp[i]*2.0) != 2 && D.n_cnv_pos[i] > 2)
			{
				G.pos_gt[i] = 1;
				G.pos_flag = true;
			}
			if (round(D.norm_dp[i]*2.0) != 2 && D.n_cnv_neg[i] > 2)
			{
				G.neg_gt[i] = 1;
				G.neg_flag = true;
			}
		}
	}
		
	// if any of three clustering meets criteria
//	if (dp_flag || pos_flag || neg_flag)
	{
		for(int i=0;i<n_sample;++i)
		{
			
			G.cn[i] = assign(D.norm_dp[i], G.Comps);

			if (G.cn[i] == -1)
			{
				if ((G.pos_gt[i] && G.neg_gt[i]) || (G.pos_gt[i] && D.n_cnv_pos[i]>4) || (G.neg_gt[i] && D.n_cnv_neg[i]>4))
				{
					cn[i] = round(D.norm_dp[i]*2.0);
				}
			}
		}
		
		G.b_biallelic = true;
		for(int i=0;i<n_sample;++i)
		{
			if (G.cn[i]>4 || G.cn[i] < 2)
			{
				G.b_biallelic = false; // Let's handle other biallelic case later
			}
			if (G.cn[i]>=0)
			{
				G.ns++;
				if (G.cn[i] != 2)
				{
					G.ac++;
				}
			}
		}
	}
}


void gtype::call_inv(sv &S, svdata& D, svgeno& G)
{
	int n_sample = D.n;

	// do clustering with positive insert size - if successful, starting breakpoint is accurate
	vector<Gaussian> C1(1);
	vector<Gaussian> C2(2);
	
	C1[0].estimate(D.norm_cnv_pos);
	C1[0].Alpha=1;
	
	C2[0].set(0,0.1);
	C2[1].set(1,0.1);
	C2[0].Alpha = C2[1].Alpha = 0.5;

	EM(D.norm_cnv_pos, C2);

	if (BIC(D.norm_cnv_pos, C2) < BIC(D.norm_cnv_pos, C1))
	{
		G.cnv_pos_flag = true;
	}
	
	C1[0].estimate(D.norm_cnv_neg);
	C1[0].Alpha=1;
	
	C2[0].set(0,0.1);
	C2[1].set(1,0.1);
	C2[0].Alpha = C2[1].Alpha = 0.5;
	
	EM(D.norm_cnv_neg, C2);
	
	if (BIC(D.norm_cnv_neg, C2) <BIC(D.norm_cnv_neg, C1))
	{
		G.cnv_neg_flag = true;
	}

	for(int i=0;i<n_sample; ++i)
	{
		if (D.n_inv_pos[i] > 2)
		{
			G.inv_pos_flag = true;
		}
		if (D.n_inv_neg[i] > 2)
		{
			G.inv_neg_flag = true;
		}
	}
	
	// if any of three clustering meets criteria
	if (cnv_pos_flag || cnv_neg_flag || inv_pos_flag || inv_neg_flag)
	{
		vector<int> support(n_sample, 0);
		
		format_output(s, ln);
		
		string filt = "FAIL";
		string info = "SVTYPE=" + s.svtype + ";END=" + to_string(s.end) + ";SVLEN=" + to_string(s.len) ;
	
		for(int i=0;i<n_sample;++i)
		{
			int cnt = 0;
			if (cnv_pos_flag && dt.norm_cnv_pos[i] > 0.5)
			{
				cnt++;
			}
			if (cnv_neg_flag && dt.norm_cnv_neg[i] > 0.5)
			{
				cnt++;
			}
			if (inv_pos_flag)
			{
				cnt++;
			}
			if (inv_neg_flag)
			{
				cnt++;
			}
			j
			if (cnt>0 && (inv_pos_flag||inv_neg_flag) && dt.norm_readcount[i] < 0.5)
			{
				// 1/1
				gt += "\t1/1";
				ns ++;
				ac +=2;
			}
			else if (cnt >0 && (inv_pos_flag||inv_neg_flag))
			{
				gt += "\t0/1";
				ns ++;
				ac +=1;
			}
			else if ((inv_pos_flag || inv_neg_flag) || cnt>1)
			{
				gt += "\t./.";
			}
			else
			{
				gt += "\t0/0";
				ns ++;
				ac +=0;
			}
			gt += ":" + to_string(dt.norm_readcount[i]).substr(0,4) + ":" + to_string((int)dt.dp[i]) + ":" + to_string((int)dt.cnv_pos[i]) + ":" + to_string((int)dt.cnv_neg[i]) + ":" + to_string((int)dt.inv_pos[i]) + ":" + to_string((int)dt.inv_neg[i]);

		}
		
		info +=";AC=" + to_string(ac) + ";NS=" + to_string(ns) + ";AF=" + to_string((double)ac/(double)(2.0*ns))  + "\tGT:PR:DP:FP:FN:IP:IN";

		if (ac>0)
		{
			filt = "PASS";
		}
		
		ln += filt + "\t" + info + "\t" + gt;
	}
	else
	{
		ln = "";
	}
}

// EM for general CNVs
void gtype::EM(vector<double>& x, vector<Gaussian>& Comps)
{
	unsigned n_sample = (unsigned) x.size();
	unsigned n_comp = (unsigned) Comps.size();
	unsigned n_iter = 15;
	
	unsigned p_count= 2;
	double p_val[n_comp];
	int zeroidx = -1;
	
	for(unsigned i=0; i<n_comp; ++i)
	{
		p_val[i] = Comps[i].Mean;
		if (p_val[i] <1e-10) // zero
		{
			zeroidx = i;
		}
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
				pr[m] = Comps[m].Alpha * normpdf(x[j], Comps[m]);
				if (zeroidx == (int)m )
				{
					pr[m] *= 2.0;
				}
				sum_p += pr[m];
			}
			
			if (sum_p > 1e-30) // if the value is an outlier, exclude it from calculations
			{
				for(unsigned m=0;m<n_comp;++m)
				{
					pr[m] /= sum_p;
//					sum[m] += pr[m] * x[j];
					sum_err[m] += pr[m] * (x[j] - Comps[m].Mean)*(x[j]-Comps[m].Mean);
					sum_pr[m] += pr[m];
				}
			}
		}
		
		// Add pseudo-count values
		for(unsigned m=0; m<n_comp; ++m)
		{
//			sum[m] += p_val[m] * p_count;
			
			sum_err[m] += (p_val[m] - Comps[m].Mean)*(p_val[m]-Comps[m].Mean) * p_count;
			sum_pr[m] += p_count;
		}
		
		// M step
		for(unsigned m=0;m<n_comp;++m)
		{
// TEMPORARY: NOT UPDATING MEANS
//			Comps[m].Mean = sum[m]/sum_pr[m];
			Comps[m].Stdev = sqrt(sum_err[m] / (sum_pr[m] ));

			Comps[m].Alpha = sum_pr[m] /( n_sample + n_comp*p_count);
			
			if (zeroidx == (int)m)
			{
				Comps[m].Mean = 0;
			}
		}
	}
}


// EM for general CNVs
void gtype::fit(vector<double>& x, vector<Gaussian>& Comps)
{
	unsigned n_sample = (unsigned) x.size();
	unsigned n_comp = (unsigned) Comps.size();
	
	// Just update the variance, do not change mean

	vector<double> sum (n_comp,0);
	vector<double> sum_pr (n_comp,0);
	vector<double> sum_err (n_comp,0);
		
	for(unsigned j=0; j<n_sample; ++j)
	{
		double sum_p = 0;
		vector<double> pr(n_comp, 0);
		for(unsigned m=0;m<n_comp;++m)
		{
			pr[m] = Comps[m].Alpha * normpdf(x[j], Comps[m]);
			if (Comps[m].Mean < 1e-10)
			{
				pr[m] *= 2.0;
			}
			sum_p += pr[m];
		}
			
		if (sum_p > 1e-30) // if the value is an outlier, exclude it from calculations
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
	
	for(unsigned m=0;m<n_comp;++m)
	{
		if (sum_pr[m] > 1e-10)
			Comps[m].Stdev = sqrt(sum_err[m] / (sum_pr[m] ));

		Comps[m].Alpha = sum_pr[m] /( n_sample);
	}
}


