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

gtype::gtype()
{
	min_bic = 0;
	bUseGL = false;
	dFlag = false;
	p_overlap = 0;
}

int gtype::assign(double x, vector<Gaussian> &C)
{
	int n_comp = (int) C.size();
	double p[n_comp];
	double max_P = -1;
	double max_R = -1;
	int ret = -1;

	for(int i=0;i<n_comp; ++i)
	{
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


void gtype::format_output(sv &s, string &ln)
{
	ln = s.chr;
	
	ln += "\t" + to_string(s.pos) + "\t" + s.svtype + "_" + s.chr + ":" + to_string(s.pos) + "-" + to_string(s.end) + "\t.\t<" + s.svtype + ">\t.\t";
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


void gtype::call_del(sv &s, svdata& dt, string &ln)
{
	int n_sample=dt.n;
	double bic1, bic2, bic2_rev, bic3;

	vector<Gaussian> C;
	vector<Gaussian> C1(1); // 1-component model
	vector<Gaussian> C2(2); // 2-component model
	vector<Gaussian> C2_rev(2); // 2-component model, with reversed
	vector<Gaussian> C3(3); // 3-component model
	
	vector<double> &X = dt.norm_dp;
	
	// For each candidate region, run EM
	// Run EM with 2 and 3 components, compare likelihoods with the 1-gaussian model, apply BIC

	
	C1[0].estimate(X);
	C1[0].Alpha = 1;
	
	double min_bic = DBL_MAX;
	bool dp_flag = false;
	double BE_dp = 1;
	
	// Single-component model
	bic1 = BIC(X, C1);
	if (bic1 < min_bic)
	{
		min_bic = bic1;
		copyComps(C, C1);
	}
	
	// Two-component model, af < 0.5
	C2[0].set(1, 0.1);
	C2[1].set(0.5, 0.1);
	C2[0].Alpha = C2[1].Alpha = 0.5;
	//EM(X, C2);
	fit(X, C2);
	bic2 = BIC(X, C2);
//	if (bic2 < min_bic && C2[0].Mean > C2[1].Mean && C2[0].Mean > 0.75 && C2[1].Mean > 0.35 && C2[1].Mean < 0.65)
	if (bic2 < min_bic) // now we do not compare mean
	{
		min_bic = bic2;
		dp_flag = true;
		BE_dp = BayesError(C2);
		copyComps(C, C2);
	}
	// Two-component model, af > 0.5
	/*
	C2_rev[0].set(0.5, 0.1);
	C2_rev[1].set(0, 0.1);
	C2_rev[0].Alpha = C2_rev[1].Alpha = 0.5;
	EM(X, C2_rev);
	bic2_rev = BIC(X, C2_rev);
	if (bic2_rev < min_bic && (C2_rev[0].Mean < 0.7 && C2_rev[0].Mean >0.3) && (C2_rev[1].Mean < 0.2))
	{
		min_bic = bic2_rev;
	//	dp_flag = true;
		BE_dp = BayesError(C2_rev);
		copyComps(C, C2_rev);
	}
	*/
	// Three-component model
	C3[0].set(1, 0.1);
	C3[1].set(0.5, 0.1);
	C3[2].set(0, 0.1);
	C3[0].Alpha = C3[1].Alpha = C3[2].Alpha = 1.0/3.0;
//	EM(X, C3);
	fit(X, C3);
	bic3 = BIC(X, C3);
//	if (bic3 < min_bic && C3[0].Mean > (C3[1].Mean+0.3) && C3[1].Mean > (C3[2].Mean+0.3))
	if (bic3 < min_bic)
	{
		min_bic = bic3;
		dp_flag = true;
		BE_dp = BayesError(C3);
		copyComps(C, C3);
	}
	// TODO: Check whether clusters are in order, after each clustering
	bool pos_flag = false;
	bool neg_flag = false;

//	if (s.len > 50)
	{
		// do clustering with positive insert size - if successful, starting breakpoint is accurate
		vector<Gaussian> C_pos1(1);
		C_pos1[0].estimate(dt.norm_cnv_pos);
		C_pos1[0].Alpha=1;
		
		vector<Gaussian> C_pos2(2);
		C_pos2[0].set(0,0.1);
		C_pos2[1].set(1,0.1);
		C_pos2[0].Alpha = C_pos2[1].Alpha = 0.5;
		fit(dt.norm_cnv_pos, C_pos2);
//		EM(dt.norm_cnv_pos, C_pos2);

//		if (BIC(dt.norm_cnv_pos, C_pos2) <BIC(dt.norm_cnv_pos, C_pos1) && C_pos2[1].Mean > 0.5 &&  C_pos2[1].Mean < 2 && C_pos2[1].Alpha > 4.5/(n_sample + 4.0))
		if (BIC(dt.norm_cnv_pos, C_pos2) <BIC(dt.norm_cnv_pos, C_pos1))
		{
			pos_flag = true;
		}
		
		// do clustering with negative insert size - if successful, ending breakpoint is accurate
		vector<Gaussian> C_neg1(1);
		C_neg1[0].estimate(dt.norm_cnv_neg);
		C_neg1[0].Alpha = 1;
		
		vector<Gaussian> C_neg2(2);
		C_neg2[0].set(0,0.1);
		C_neg2[1].set(1,0.1);
		C_neg2[0].Alpha = C_neg2[1].Alpha = 0.5;
		
	//	EM(dt.norm_cnv_neg, C_neg2);
		fit(dt.norm_cnv_neg, C_neg2);
		
//		if (BIC(dt.norm_cnv_neg, C_neg2) <BIC(dt.norm_cnv_neg, C_neg1) && C_neg2[1].Mean > 0.5 && C_neg2[1].Mean < 2 && C_neg2[1].Alpha > 4.5/(n_sample + 4.0))
		if (BIC(dt.norm_cnv_neg, C_neg2) <BIC(dt.norm_cnv_neg, C_neg1))
		{
			neg_flag = true;
		}
	}

	
	// if any of three clustering meets criteria
	if (dp_flag || pos_flag || neg_flag)
	{
		if (!dp_flag)
		{
			// Depth-based clustering failed, force 3-component model
			if (bic2 < bic3)
			{
				C.clear();
				copyComps(C, C2);
			}
			else
			{
				C.clear();
				copyComps(C, C3);
			}
			BE_dp = BayesError(C);
		}

		vector<int> gt(n_sample, 0);

		format_output(s, ln);
		
		string filt = "FAIL";
		string info = "SVTYPE=" + s.svtype + ";END=" + to_string(s.end) + ";SVLEN=" + to_string(s.len) + ";NCLUS=" + to_string(C.size()) + ";P_OVERLAP=" + to_string(BE_dp);
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
		
		info += ";MEAN=" + to_string(C[0].Mean);
		for(unsigned i=1;i<C.size(); ++i)
		{
			info += "," + to_string(C[i].Mean);
		}
		info += ";STDEV=" + to_string(C[0].Stdev);
		for(unsigned i=1;i<C.size(); ++i)
		{
			info += "," + to_string(C[i].Stdev);
		}
		info += ";FRAC=" + to_string(C[0].Alpha);
		for(unsigned i=1;i<C.size(); ++i)
		{
			info += "," + to_string(C[i].Alpha);
		}
		
		int ac = 0;
		int ns = 0;

/*
		for(int m=0; m<(int)C.size(); ++m)
		{
			C[m].Stdev = 0.1;
		}
		
		*/
		for(int i=0;i<n_sample;++i)
		{
			int cn;
			bool pos_gt = false, neg_gt = false;
			
			cn = assign(X[i], C);
			if (pos_flag && dt.norm_cnv_pos[i] > 0.5 && dt.norm_cnv_pos[i] < 2 && s.len>50 && dt.n_cnv_pos[i] > 1)
			{
				pos_gt = true;
			}
			if (neg_flag && dt.norm_cnv_neg[i] > 0.5 && dt.norm_cnv_neg[i] < 2 && s.len>50 && dt.n_cnv_neg[i] > 1)
			{
				neg_gt = true;
			}
			
			gt[i] = -1;
			switch(cn)
			{
				case -1:
					if (pos_gt && neg_gt)
					{
						if ( 0.2 <= dt.norm_dp[i] && dt.norm_dp[i]< 0.8)
						{
							gt[i] = 1;
						}
						else if (dt.norm_dp[i] < 0.2)
						{
							gt[i] = 2;
						}
					}
					else if (pos_gt || neg_gt)
					{
						if ( 0.2 <= dt.norm_dp[i] && dt.norm_dp[i] < 0.7)
						{
							gt[i] = 1;
						}
						else if (dt.norm_dp[i] < 0.2)
						{
							gt[i] = 2;
						}
					}
					break;
				case 0:
					gt[i] = 2;
					break;
				case 1:
					gt[i] = 1;
					break;
				case 2:
					if (pos_gt && neg_gt)
					{
						if (dt.norm_dp[i] < 0.8)
						{
							gt[i] = 1;
						}
						else
						{
							gt[i] = 0;
						}
					}
					else if (pos_gt || neg_gt)
					{
						if (dt.norm_dp[i] < 0.75)
						{
							gt[i] = 1;
						}
						else
						{
							gt[i] = 0;
						}
					}
					else
					{
						gt[i] = 0;
					}
					break;
			}
			if (gt[i] >=0 )
			{
				ac += gt[i];
				ns += 1;
			}
		}
		
		info +=";AC=" + to_string(ac) + ";NS=" + to_string(ns) + ";AF=" + to_string((double)ac/(double)(2.0*ns))  + "\tGT:CN:ND:DP:FP:FN";
		
		if ((ac>0 &&  ac < (2*n_sample) && ns>(n_sample *0.5)) && ((dp_flag && BE_dp<BE_THRESHOLD) ||  (BE_dp < BE_THRESHOLD *2 && (neg_flag || pos_flag)) || (BE_dp < 0.5&& neg_flag && pos_flag)))
		{
			filt = "PASS";
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
			ln += ":" + to_string(dt.norm_dp[i]).substr(0,4) + ":" + to_string((int)dt.dp[i]) + ":" + to_string((int)dt.cnv_pos[i]) + ":" + to_string((int)dt.cnv_neg[i]);
		}
	}
	else
	{
		ln = "";
	}

}

void gtype::call_cnv(sv &s, svdata& dt, string &ln)
{
	int n_sample = dt.n;
	vector<double> &X = dt.norm_dp;
	
	vector<Gaussian> Comps(1);
	
	double min_BIC = DBL_MAX;
	vector<double> bic(10,0);
	
	Comps[0].estimate(X);
	Comps[0].Alpha = 1;
	
	min_BIC = bic[0] = BIC(X, Comps);
	//		cerr << "1 comp, bic " << min_BIC << endl;
	unsigned n_comp = 1;
	
	// Check up to 30 -- should be enough? -- let's go for 10 now
	// It seems like too many components make probs too small hence failing
	
	for(unsigned i=2; i<10; ++i)
	{
		Gaussian* pG = new Gaussian;
		Comps.push_back(*pG);
		
		for(unsigned j=0; j<i; ++j)
		{
			Comps[j].Mean = (double)j*0.5;
			Comps[j].Stdev = 0.1;
			Comps[j].Alpha = 1;
		}
		// Run EM with i components
//		EM(X, Comps);
		fit(X, Comps);
		
		bic[i-1] = BIC(X, Comps);
		
		if (bic[i-1]<min_BIC)
		{
			min_BIC = bic[i-1];
			n_comp = i;
		}
	}
	
	Comps.clear();
	Comps.resize(n_comp);
	
	for(unsigned j=0;j<n_comp;++j)
	{
		Comps[j].Mean = (double)j*0.5;
		Comps[j].Stdev = 0.1;
		Comps[j].Alpha = 1;
	}

//	EM(X, Comps);
	fit(X, Comps);
	
	bool dp_flag = false;
	double BE = 1;

	if (n_comp>1)
	{
		BE = BayesError(Comps);
		dp_flag = (BE<BE_THRESHOLD);
	}

	bool pos_flag = false;
	bool neg_flag = false;
	
//	if (s.len > 50)
	{
		// do clustering with positive insert size - if successful, starting breakpoint is accurate
		vector<Gaussian> C_pos1(1);
		C_pos1[0].estimate(dt.norm_cnv_pos);
		C_pos1[0].Alpha=1;
		
		vector<Gaussian> C_pos2(2);
		C_pos2[0].set(0,0.1);
		C_pos2[1].set(1,0.1);
		C_pos2[0].Alpha = C_pos2[1].Alpha = 0.5;
		//EM(dt.norm_cnv_pos, C_pos2);
		fit(dt.norm_cnv_pos, C_pos2);
		
//		if (BIC(dt.norm_cnv_pos, C_pos2) <BIC(dt.norm_cnv_pos, C_pos1) && C_pos2[1].Mean > 0.5 &&  C_pos2[1].Mean < 2 && C_pos2[1].Alpha > 4.5/(n_sample + 4.0))
		if (BIC(dt.norm_cnv_pos, C_pos2) <BIC(dt.norm_cnv_pos, C_pos1))
		{
			pos_flag = true;
		}
		
		// do clustering with negative insert size - if successful, ending breakpoint is accurate
		
		vector<Gaussian> C_neg1(1);
		C_neg1[0].estimate(dt.norm_cnv_neg);
		C_neg1[0].Alpha = 1;
		
		vector<Gaussian> C_neg2(2);
		C_neg2[0].set(0,0.1);
		C_neg2[1].set(1,0.1);
		C_neg2[0].Alpha = C_neg2[1].Alpha = 0.5;
		
		//EM(dt.norm_cnv_neg, C_neg2);
		fit(dt.norm_cnv_neg, C_neg2);
		
//		if (BIC(dt.norm_cnv_neg, C_neg2) <BIC(dt.norm_cnv_neg, C_neg1) && C_neg2[1].Mean > 0.5 && C_neg2[1].Mean < 2 && C_neg2[1].Alpha > 4.5/(n_sample + 4.0))
		if (BIC(dt.norm_cnv_neg, C_neg2) <BIC(dt.norm_cnv_neg, C_neg1))
		{
			neg_flag = true;
		}
	}
	
	// if any of three clustering meets criteria
	if (dp_flag || pos_flag || neg_flag)
	{
		vector<int> cn(n_sample, 0);
		
		format_output(s, ln);
		
		string filt = "FAIL";
		string info = "SVTYPE=" + s.svtype + ";END=" + to_string(s.end) + ";SVLEN=" + to_string(s.len) + ";NCLUS=" + to_string(Comps.size()) + ";P_OVERLAP=" + to_string(BE);
		
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
		
		int ac = 0;
		int ns = 0;
		
		for(int i=0;i<n_sample;++i)
		{
			bool pos_gt = false, neg_gt = false;
			
			cn[i] = assign(X[i], Comps);
			if (pos_flag && dt.n_cnv_pos[i]>1)
			{
				pos_gt = true;
			}
			if (neg_flag && dt.n_cnv_neg[i] >1)
			{
				neg_gt = true;
			}
			
			if (cn[i] == -1)
			{
				if (pos_gt || neg_gt)
				{
					cn[i] = round(dt.norm_dp[i]*2.0);
				}
			}
		}
		
		bool b_biallelic = true;
		for(int i=0;i<n_sample;++i)
		{
			if (cn[i]>4 || cn[i] ==1 || cn[i] == 0)
			{
				b_biallelic = false; // Let's handle other biallelic case (e.g. CN=2,3,5)
			}
		}
		string gt = "";
		
		if (b_biallelic)
		{

			for(int i=0; i<n_sample; ++i )
			{
				switch(cn[i])
				{
					case 2:
						gt += "\t0/0:2";
						ns++;
						break;
					case 3:
						gt += "\t0/1:3";
						ac++;
						ns++;
						break;
					case 4:
						gt += "\t1/1:4";
						ac+=2;
						ns++;
						break;
					case -1:
						gt += "\t.:.";
						break;
				}
				gt += ":" + to_string(dt.norm_dp[i]).substr(0,4) + ":" + to_string((int)dt.dp[i]) + ":" + to_string((int)dt.cnv_pos[i]) + ":" + to_string((int)dt.cnv_neg[i]);

			}
			info +=";AC=" + to_string(ac) + ";NS=" + to_string(ns) + ";AF=" + to_string((double)ac/(double)(2.0*ns))  + "\tGT:CN:ND:DP:FP:FN";
		}
		else
		{
			for(int i=0; i<n_sample; ++i )
			{
				if (cn[i] >= 0)
				{
					gt += "\t" + to_string(cn[i]);
					ns ++;
					if (cn[i]>2)
					{
						ac++;
					}
				}
				else
				{
					gt += "\t.";
				}
				gt += ":" + to_string(dt.norm_dp[i]).substr(0,4) + ":" + to_string((int)dt.dp[i]) + ":" + to_string((int)dt.cnv_pos[i]) + ":" + to_string((int)dt.cnv_neg[i]);
			}
			info +=";AC=" + to_string(ac) + ";NS=" + to_string(ns) + ";AF=" + to_string((double)ac/(double)(ns))  + "\tCN:ND:DP:FP:FN";
		}
				
		
		if ((ac>0 &&  ac < (2*n_sample) && ns>(n_sample *0.5)) && ((dp_flag && BE<BE_THRESHOLD) ||  (BE < BE_THRESHOLD *2 && (neg_flag || pos_flag)) || (BE < 0.5 && neg_flag && pos_flag)))
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


void gtype::call_inv(sv &s, svdata& dt, string &ln)
{
	int n_sample = dt.n;
	
	bool cnv_pos_flag = false;
	bool cnv_neg_flag = false;

	bool inv_pos_flag = false;
	bool inv_neg_flag = false;
	
	
	// do clustering with positive insert size - if successful, starting breakpoint is accurate
	vector<Gaussian> C1(1);
	vector<Gaussian> C2(2);
	
	C1[0].estimate(dt.norm_cnv_pos);
	C1[0].Alpha=1;
	
	C2[0].set(0,0.1);
	C2[1].set(1,0.1);
	C2[0].Alpha = C2[1].Alpha = 0.5;

	//EM(dt.norm_cnv_pos, C2);
	fit(dt.norm_cnv_pos, C2);
	
	//if (BIC(dt.norm_cnv_pos, C2) <BIC(dt.norm_cnv_pos, C1) && C2[1].Mean > 0.5 &&  C2[1].Mean < 3 && C2[1].Alpha > 4.5/(n_sample + 4.0))
	if (BIC(dt.norm_cnv_pos, C2) <BIC(dt.norm_cnv_pos, C1))
	{
		cnv_pos_flag = true;
	}
	
	C1[0].estimate(dt.norm_cnv_neg);
	C1[0].Alpha=1;
	
	C2[0].set(0,0.1);
	C2[1].set(1,0.1);
	C2[0].Alpha = C2[1].Alpha = 0.5;
	
//	EM(dt.norm_cnv_neg, C2);
	fit(dt.norm_cnv_neg, C2);
	
	//if (BIC(dt.norm_cnv_neg, C2) <BIC(dt.norm_cnv_neg, C1) && C2[1].Mean > 0.5 &&  C2[1].Mean < 3 && C2[1].Alpha > 4.5/(n_sample + 4.0))
	if (BIC(dt.norm_cnv_neg, C2) <BIC(dt.norm_cnv_neg, C1))
	{
		cnv_neg_flag = true;
	}
	
	C1[0].estimate(dt.norm_inv_pos);
	C1[0].Alpha=1;
	
	C2[0].set(0,0.1);
	C2[1].set(1,0.1);
	C2[0].Alpha = C2[1].Alpha = 0.5;
	
	//EM(dt.norm_inv_pos, C2);
	fit(dt.norm_inv_pos, C2);
	
	//if (BIC(dt.norm_inv_pos, C2) <BIC(dt.norm_inv_pos, C1) && C2[1].Mean > 0.5 &&  C2[1].Mean < 3 && C2[1].Alpha > 4.5/(n_sample + 4.0))
	if (BIC(dt.norm_inv_pos, C2) <BIC(dt.norm_inv_pos, C1))
	{
		cnv_pos_flag = true;
	}
	
	/*
	C1[0].estimate(dt.norm_inv_neg);
	C1[0].Alpha=1;
	
	C2[0].set(0,0.1);
	C2[1].set(1,0.1);
	C2[0].Alpha = C2[1].Alpha = 0.5;
	
	//EM(dt.norm_inv_neg, C2);
	fit(dt.norm_inv_neg, C2);
	
//	if (BIC(dt.norm_inv_neg, C2) <BIC(dt.norm_inv_neg, C1) && C2[1].Mean > 0.5 &&  C2[1].Mean < 3 && C2[1].Alpha > 4.5/(n_sample + 4.0))
	if (BIC(dt.norm_inv_neg, C2) <BIC(dt.norm_inv_neg, C1))
		*/
	for(int i=0;i<n_sample; ++i)
	{
		if (dt.n_inv_neg[i] >0 && dt.n_inv_pos[i] >0)
			cnv_neg_flag = true;
	}
	
	// if any of three clustering meets criteria
	if (cnv_pos_flag || cnv_neg_flag || inv_pos_flag || inv_neg_flag)
	{
		vector<int> support(n_sample, 0);
		
		format_output(s, ln);
		
		string filt = "FAIL";
		string info = "SVTYPE=" + s.svtype + ";END=" + to_string(s.end) + ";SVLEN=" + to_string(s.len) ;
		
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
		
		int ac = 0;
		int ns = 0;
		
		string gt="";

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
	unsigned n_iter = 20;
	
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
					sum[m] += pr[m] * x[j];
					sum_err[m] += pr[m] * (x[j] - Comps[m].Mean)*(x[j]-Comps[m].Mean);
					sum_pr[m] += pr[m];
				}
			}
		}
		
		// Add pseudo-count values
		for(unsigned m=0; m<n_comp; ++m)
		{
			sum[m] += p_val[m] * p_count;
			
			sum_err[m] += (p_val[m] - Comps[m].Mean)*(p_val[m]-Comps[m].Mean) * p_count;
			sum_pr[m] += p_count;
		}
		
		// M step
		for(unsigned m=0;m<n_comp;++m)
		{

			Comps[m].Mean = sum[m]/sum_pr[m];
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


