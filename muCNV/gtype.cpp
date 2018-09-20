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

bool ordered(vector<Gaussian> &C)
{
	for(unsigned i=0; i<C.size()-1;++i)
	{
		if (C[i].Mean<=C[i+1].Mean)
			return false;
	}
	return true;
}

bool r_ordered(vector<Gaussian> &C)
{
	for(unsigned i=0; i<C.size()-1;++i)
	{
		if (C[i].Mean>=C[i+1].Mean)
			return false;
	}
	return true;
}

void copy(Gaussian &x, Gaussian &y)
{
	x.set(y.Mean, y.Stdev);
	x.Alpha = y.Alpha;
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
// Major allele might domiate Alpha when sample size is large 
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

	int up = ceil(x*2.0);
	int down = floor(x*2.0);

	if (ret != up && ret != down)
	{
		return -1;
	}
	
	if (max_R > 0.1) //   P_THRESHOLD)
	{
		return -1;
	}
	
	return ret;
}

void svgeno::print(sv &S, svdata &D, string &ln, vector<double>& wt)
{
	ln = to_string(S.chrnum);
	ln += "\t" + to_string(S.pos) + "\t" + svTypeName(S.svtype) + "_" + to_string(S.chrnum) + ":" + to_string(S.pos) + "-" + to_string(S.end) + "\t.\t<" + svTypeName(S.svtype) + ">\t.\t";

	if (b_pass)
	{
		ln += "PASS\t";
	}
	else if (b_dump)
	{
		ln += "FAIL;DUMP\t";
	}
	else
	{
		ln += "FAIL\t";
	}
	
	ln += "SVTYPE=" + string(svTypeName(S.svtype)) + ";END=" + to_string(S.end) + ";SVLEN=" + to_string(S.len) + ";AC=" + to_string(ac) + ";NS=" + to_string(ns) + ";AF=";
	if (ns>0)
		ln+=to_string((double)ac/(double)(2.0*ns));
	else
		ln+="0";
	
	ln+= info;

	bool bic_flag = false ;
	if (Comps.size()>0)
	{
		ln+= ";NCLUS=" + to_string(Comps.size()) + ";P_OVERLAP=" + to_string(p_overlap);

		ln+= ";BIC_DP="+ to_string(bic[0]) + "," + to_string(bic[1]) + "," + to_string(bic[2]);
		for(unsigned i=3;i<Comps.size(); ++i)
		{
			ln += "," + to_string(bic[i]);
			if (bic[i]<bic[0])
			{
				bic_flag = true;
			}
		}
	}

	if (bic_flag)
		ln += ";BIC";

	if (dp_flag)
		ln += ";DP";
	if (pos_flag)
		ln += ";POS";
	if (neg_flag)
		ln += ";NEG";

	if (Comps.size()>0)
	{
		ln += ";MEAN=" + to_string(Comps[0].Mean);

		for(unsigned i=1;i<Comps.size(); ++i)
			ln += "," + to_string(Comps[i].Mean);

		ln += ";STDEV=" + to_string(Comps[0].Stdev);
		for(unsigned i=1;i<Comps.size(); ++i)
			ln += "," + to_string(Comps[i].Stdev);

		ln += ";FRAC=" + to_string(Comps[0].Alpha);
		for(unsigned i=1;i<Comps.size(); ++i)
			ln += "," + to_string(Comps[i].Alpha);
	}

	if (b_biallelic)
		ln += "\tGT:CN:ND:DP:FP:FN:WT";
	//	ln += "\tGT";
	else
		ln += "\tCN:ND:DP:FP:FN";
//		ln += "\tCN";
	
	for (int i=0; i<n_sample; ++i)
	{
		if (b_biallelic)
		{
			switch(gt[i])
			{
				case 0:
					ln += "\t0/0";
					break;
				case 1:
					ln += "\t0/1";
					break;
				case 2:
					ln += "\t1/1";
					break;
				default:
					ln += "\t.";
					break;
			}
			ln+= ":";
		}
		else
		{
			ln += "\t";
		}

		if (S.svtype == INV)
		{
			if (cn[i]<0)
			{
				ln += to_string(cn[i]) + ":" + to_string(D.norm_readcount[i]).substr(0,4) + ":" + to_string((int)D.dp[i]) + ":" +  to_string(D.n_inv_pos[i]) + "," + to_string((int)D.inv_pos[i]) + ":" + to_string(D.n_inv_neg[i]) + "," + to_string((int)D.inv_neg[i]);
			}
			else
			{
				ln += ".:" + to_string(D.norm_readcount[i]).substr(0,4) + ":" + to_string((int)D.dp[i]) + ":" +  to_string(D.n_inv_pos[i]) + "," + to_string((int)D.inv_pos[i]) + ":" + to_string(D.n_inv_neg[i]) + "," + to_string((int)D.inv_neg[i]) ;
			}
		}
		else
		{
			if (cn[i] <0)
			{
				ln +=  ".:";
			}
			else 
			{
				ln += to_string(cn[i])+":";
			}
			ln += to_string(D.norm_dp[i]).substr(0,4) + ":" + to_string((int)D.dp[i]) + ":" + to_string((int)D.n_cnv_pos[i]) + "," + to_string((int)D.cnv_pos[i]) + ":" + to_string((int)D.n_cnv_neg[i]) + "," + to_string((int)D.cnv_neg[i]) + ":" + to_string((int)wt[i]); 
		}
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


void gtype::call_tmp(sv &S, svdata& D, svgeno &G, vector<double> &avg_isz, vector<double>&std_isz, vector<double> &wt)
{
	int n_sample=D.n;

	double sum_dp = 0;
	int n_dp = 0;
	int n_pos = 0;
	int n_neg = 0;
	int n_other = 0;

	double sum_var_pos = 0;
	double sum_var_neg = 0;

	double sum_other_pos = 0;
	double sum_other_neg = 0;

	G.b_biallelic = true;
	G.p_overlap = 0; 

	vector<Gaussian> C(2);
	C[0].set(1,0.1);
	C[1].set(0.5,0.1);
	copyComps(G.Comps, C);

	for(int i=0; i<n_sample; ++i)
	{
		G.cn[i] = 2;
		if (wt[i]>0)
		{
			G.gt[i] = 1;
			n_dp ++;

			if (S.svtype == INV)
			{
				sum_dp += D.norm_readcount[i];

				if (D.n_inv_pos[i]>2 && D.norm_inv_pos[i] < 1.3 && D.norm_inv_pos[i]> 0.7)
				{
					sum_var_pos += D.inv_pos[i];
					n_pos ++;
				}
				if (D.n_inv_neg[i]>2 && D.norm_inv_neg[i] < 1.3 && D.norm_inv_neg[i]> 0.7)
				{
					sum_var_neg += D.inv_neg[i];
					n_neg ++;
				}
			}
			else
			{
				sum_dp += D.norm_dp[i];

				if (D.n_cnv_pos[i]>2 && D.norm_cnv_pos[i] <1.3 && D.norm_cnv_pos[i]>0.7)
				{
					sum_var_pos += D.cnv_pos[i];
					n_pos ++;
				}
				if (D.n_cnv_neg[i]>2 && D.norm_cnv_neg[i]<1.3 && D.norm_cnv_neg[i]>0.7)
				{
					sum_var_neg += D.cnv_neg[i];
					n_neg ++;
				}
			}
		}
		else
		{
			G.gt[i] = 0;
			if (S.svtype == INV)
			{
				sum_other_pos += D.inv_pos[i];
				sum_other_neg += D.inv_neg[i];
			}
			else
			{
				sum_other_pos += D.cnv_pos[i];
				sum_other_neg += D.cnv_neg[i];
			}
			n_other ++;
		}
	}

	double avg_var_dp = (n_dp>0) ? sum_dp/n_dp : 1;

	double avg_var_pos = (n_pos>0) ? sum_var_pos/n_pos : 0;
	double avg_var_neg = (n_neg>0) ? sum_var_neg/n_neg : 0;

	double avg_other_pos = (n_other>0) ? sum_other_pos/n_other : 0;
	double avg_other_neg = (n_other>0) ? sum_other_neg/n_other : 0;

	G.info = ";POS=(" + to_string((int)avg_other_pos) + "," + to_string((int)avg_var_pos) + ")";
	G.info += ";NEG=(" + to_string((int)avg_other_neg) + "," + to_string((int)avg_var_neg) + ")";
	C[0].set(1,0.1);
	C[1].set(avg_var_dp, 0.1);
	copyComps(G.Comps, C);

	G.ac = 0;
	G.ns = 0;
	G.b_pass = true;
	G.b_dump = false;

	for(int i=0; i<n_sample; ++i)
	{
		if (G.gt[i] == 0)
		{
			if (n_pos>0 && n_neg > 0 && S.len>300) 
			{
				double dpos1, dneg1, dpos2, dneg2;
				if (S.svtype == INV)
				{
					dpos1 = (D.inv_pos[i] - avg_var_pos);
					dneg1 = (D.inv_neg[i] - avg_var_neg);
					dpos2 = (D.inv_pos[i] - avg_other_pos);
					dneg2 = (D.inv_neg[i] - avg_other_neg);
				}
				else
				{
					dpos1 = (D.cnv_pos[i] - avg_var_pos);
					dneg1 = (D.cnv_neg[i] - avg_var_neg);
					dpos2 = (D.cnv_pos[i] - avg_other_pos);
					dneg2 = (D.cnv_neg[i] - avg_other_neg);
				}


				if ((dpos1*dpos1 + dneg1*dneg1) < (dpos2*dpos2 + dneg2*dneg2))
				{
					G.gt[i] = -1;
				}
				else
				{
					G.ns++;
				}
			}
			/*
			else if (S.svtype == "INV" &&  abs(D.norm_readcount[i]-avg_var_dp) < abs(D.norm_readcount[i] - 1))
			{
				G.gt[i] = -1;
			}
			*/
			else if (S.svtype != INV &&  abs(D.norm_dp[i]-avg_var_dp) < abs(D.norm_dp[i] - 1))
			{
				G.gt[i] = -1;
			}
			else
			{
				G.ns ++;
			}
		}
		else
		{
			if (S.svtype == DEL && D.norm_dp[i]<0.15)
			{
				G.gt[i] = 2;
				G.ac += 2;
			}
			else if (S.svtype == DUP && D.norm_dp[i]>1.8)
			{
				G.gt[i] = 2;
				G.ac += 2;
			}
			else if (S.svtype == INV && D.norm_readcount[i] < 0.15)
			{
				G.gt[i] = 2;
				G.ac += 2;
			}
			else
			{
				G.ac ++;
			}
			G.ns ++;
		}
	}
	if (G.ns < 0.2*n_sample)
	{
		G.b_pass = false;
	}
	if (S.svtype == DEL && avg_var_dp > 1)
	{
		G.b_pass = false;
	}
	if (S.svtype == DUP && avg_var_dp < 1)
	{
		G.b_pass = false;
	}
}


void gtype::call_del_tmp(sv &S, svdata& D, svgeno &G, vector<double> &avg_isz, vector<double>&std_isz, vector<double> &wt)
{
	int n_sample=D.n;

	G.b_biallelic = true;
	G.p_overlap = 0; 

	vector<double> X;
	vector<double> Y;
	vector<double> X2;
	vector<double> P_X;
	vector<double> N_X;
	vector<double> P_Y;
	vector<double> N_Y;

	X.reserve(n_sample);
	X2.reserve(n_sample);
	Y.reserve(n_sample);

	P_X.reserve(n_sample);
	N_X.reserve(n_sample);

	P_Y.reserve(n_sample);
	N_Y.reserve(n_sample);

	for(int i=0; i<n_sample; ++i)
	{
		G.cn[i] = 2;
		if (wt[i]>0)
		{
			if (D.norm_dp[i] < 0.2)
			{
				X2.push_back(D.norm_dp[i]);
			}
			else
			{
				X.push_back(D.norm_dp[i]);
			}

			if (S.len>150)
			{
				if (D.n_cnv_pos[i]>2 && D.norm_cnv_pos[i]<1.3 && D.norm_cnv_pos[i]>0.7)
				{
					P_X.push_back(D.norm_cnv_pos[i]);
				}
				if (D.n_cnv_neg[i]>2 && D.norm_cnv_neg[i]<1.3 && D.norm_cnv_neg[i]>0.7)
				{
					N_X.push_back(D.norm_cnv_neg[i]);
				}
			}
		}
	}

	vector<Gaussian> C1(1);
	vector<Gaussian> C2(2);
	vector<Gaussian> C3(3);

	C1[0].estimate(D.norm_dp);
	C1[0].Alpha = 1;
	copyComps(G.Comps, C1);

	C2[1].estimate(X);

	if (C2[1].Mean > 0.8 && P_X.size() == 0 && N_X.size() == 0)
	{
		// no good signes in both depth & read pair
		G.b_pass = false;
		return;
	}

	for(int i=0;i<n_sample; ++i)
	{
		if (wt[i]==0 && abs(D.norm_dp[i] - 1) < abs(D.norm_dp[i] - C2[1].Mean))
		{
			if (S.len<=150 || D.n_cnv_pos[i]==0 || D.n_cnv_neg[i]==0 || D.norm_cnv_pos[i] <0.7 || D.norm_cnv_pos[i]>1.3)
			{
				Y.push_back(D.norm_dp[i]);
				P_Y.push_back(D.norm_cnv_pos[i]);
				N_Y.push_back(D.norm_cnv_neg[i]);
			}
		}
	}

	C2[0].estimate(Y);
	C2[0].Alpha = ((double)Y.size()) / (X.size() + Y.size());
	C2[1].Alpha = ((double)X.size()) / (X.size() + Y.size());

	// Fit 1 or 2 comp
	G.bic[0] = BIC(D.norm_dp, C1);
	G.bic[1] = BIC(D.norm_dp, C2);

	if (G.bic[1] < G.bic[0])
	{
		G.dp_flag = true;

		// ToDo : 3-component, construct 3-components from the beginning based on depth
		G.p_overlap = BayesError(C2);
		copyComps(G.Comps, C2);
	}
	if (X2.size()>0)
	{
		C3[0].estimate(Y);
		C3[1].estimate(X);
		C3[2].estimate(X2);

		C3[0].Alpha = ((double)Y.size()) / (X.size() + Y.size() + X2.size());
		C3[1].Alpha = ((double)X.size()) / (X.size() + Y.size() + X2.size());
		C3[2].Alpha = ((double)X2.size()) / (X.size() + Y.size() + X2.size());
		G.bic[2] = BIC(D.norm_dp, C3);

		if (G.bic[2] < G.bic[1] && G.bic[2] < G.bic[0])
		{
			G.dp_flag = true;

			G.p_overlap = BayesError(C3);
			copyComps(G.Comps, C3);
		}
	}

	vector<Gaussian> P1(1);
	vector<Gaussian> P2(2);
	if (P_X.size() > 0)
	{

		P1[0].estimate(D.cnv_pos);
		P1[0].Alpha = 1;

		P2[0].estimate(P_Y);
		P2[1].estimate(P_X);
		P2[0].Alpha = (double) P_Y.size() / (P_X.size() + P_Y.size());
		P2[1].Alpha = (double) P_X.size() / (P_X.size() + P_Y.size());

		if (BIC(D.norm_cnv_pos, P2) < BIC(D.norm_cnv_pos, P1))
		{
			G.pos_flag = true;
			G.info = ";POS=(" + to_string((int)(P2[0].Mean*S.len)) + "," + to_string((int)(P2[1].Mean*S.len)) + ")";
		}
	}

	vector<Gaussian> N1(1);
	vector<Gaussian> N2(2);
	if (N_X.size() > 0)
	{
		N1[0].estimate(D.cnv_neg);
		N1[0].Alpha = 1;

		N2[0].estimate(N_Y);
		N2[1].estimate(N_X);
		N2[0].Alpha = (double) N_Y.size() / (N_X.size() + N_Y.size());
		N2[1].Alpha = (double) N_X.size() / (N_X.size() + N_Y.size());

		if (BIC(D.norm_cnv_neg, N2) < BIC(D.norm_cnv_neg, N1))
		{
			G.neg_flag = true;
			G.info += ";NEG=(" + to_string((int)(N2[0].Mean*S.len)) + "," + to_string((int)(N2[1].Mean*S.len)) + ")";
		}
	}

	G.ac = 0;
	G.ns = 0;
	G.b_pass = false;
	G.b_dump = true;

	if (G.dp_flag || G.pos_flag || G.neg_flag)
	{
		bool b_hom = false;
		for(int i=0; i<n_sample; ++i)
		{
			int pos_gt = -1;
			int neg_gt = -1;
			G.cn[i] = -1;
			if (G.dp_flag)
			{
				G.cn[i] = assign(D.norm_dp[i], G.Comps);
			}

			if (G.pos_flag && D.n_cnv_pos[i] > 1 && D.norm_dp[i] < 0.9)
			{
				double p0 = P2[0].pdf(D.norm_cnv_pos[i]);
				double p1 = P2[1].pdf(D.norm_cnv_pos[i]);

				if (p1 > 5.0 * p0)
				{
					pos_gt = 1;
				}
				else if (p0 > p1)
				{
					pos_gt = 0;
				}
			}

			if (G.neg_flag && D.n_cnv_neg[i] > 1 && D.norm_dp[i] < 0.9)
			{
				double p0 = N2[0].pdf(D.norm_cnv_neg[i]);
				double p1 = N2[1].pdf(D.norm_cnv_neg[i]);
				if (p1 > 5.0 * p0)
				{
					neg_gt = 1;
				}
				else if (p0 > p1)
				{
					neg_gt = 0;
				}
			}

			if (G.cn[i] == 0 || (pos_gt==1 && neg_gt==1 && D.norm_dp[i] < 0.2))
			{
				G.gt[i] = 2;
				G.ac+=2;
				G.ns++;
				b_hom = true;
			}
			else if (G.cn[i] == 1 || (pos_gt==1 && neg_gt==1))
			{
				G.gt[i] = 1;
				G.ac++;
				G.ns++;
			}
			else if (G.cn[1] == 2 || (pos_gt==0 && neg_gt==0))
			{
				G.gt[i] = 0;
				G.ns++;
			}
		}

		G.b_dump = false;

		if (G.ac > 0.95*G.ns && b_hom == false)
		{
			G.b_pass=  false;
		}
		else if (G.ns > 0.5*n_sample && G.ac > 0 && G.ac < (G.ns*2.0))
		{
			G.b_pass = true;
		}
	}
	else
	{
		G.b_dump = true;
	}
	
}


void gtype::call_dup_tmp(sv &S, svdata& D, svgeno &G, vector<double> &avg_isz, vector<double>&std_isz, vector<double> &wt)
{
	int n_sample=D.n;

	G.b_biallelic = true;
	G.p_overlap = 0; 

	vector<double> X;
	vector<double> Y;
	vector<double> X2;
	vector<double> X3;
	vector<double> X4;
	vector<double> X5;
	vector<double> X6;
	vector<double> P_X;
	vector<double> N_X;
	vector<double> P_Y;
	vector<double> N_Y;

	X.reserve(n_sample);
	X2.reserve(n_sample);
	X3.reserve(n_sample);
	X4.reserve(n_sample);
	X5.reserve(n_sample);
	X6.reserve(n_sample);

	Y.reserve(n_sample);

	P_X.reserve(n_sample);
	N_X.reserve(n_sample);

	P_Y.reserve(n_sample);
	N_Y.reserve(n_sample);

	for(int i=0; i<n_sample; ++i)
	{
		G.cn[i] = 2;
		if (wt[i]>0)
		{
			if (D.norm_dp[i] > 3.8)
			{
				X6.push_back(D.norm_dp[i]);
			}
			else if (D.norm_dp[i] > 3.3)
			{
				X5.push_back(D.norm_dp[i]);
			}
			else if (D.norm_dp[i] > 2.8)
			{
				X4.push_back(D.norm_dp[i]);
			}
			else if (D.norm_dp[i] > 2.3)
			{
				X3.push_back(D.norm_dp[i]);
			}
			else if (D.norm_dp[i] > 1.8)
			{
				X2.push_back(D.norm_dp[i]);
			}
			else
			{
				X.push_back(D.norm_dp[i]);
			}

			if (S.len>500)
			{
				if (D.n_cnv_pos[i]>2 && D.norm_cnv_pos[i]<1.3 && D.norm_cnv_pos[i]>0.7)
				{
					P_X.push_back(D.norm_cnv_pos[i]);
				}
				if (D.n_cnv_neg[i]>2 && D.norm_cnv_neg[i]<1.3 && D.norm_cnv_neg[i]>0.7)
				{
					N_X.push_back(D.norm_cnv_neg[i]);
				}
			}
		}
	}

	vector<Gaussian> C1(1);
	vector<Gaussian> C2(2);
	vector<Gaussian> C3(3);
	vector<Gaussian> C4(4);
	vector<Gaussian> C5(5);
	vector<Gaussian> C6(6);
	vector<Gaussian> C7(7);

	C1[0].estimate(D.norm_dp);
	C1[0].Alpha = 1;
	copyComps(G.Comps, C1);

	C2[1].estimate(X);

	if (C2[1].Mean > 0.8 && P_X.size() == 0 && N_X.size() == 0)
	{
		// no good signes in both depth & read pair
		G.b_pass = false;
		return;
	}

	for(int i=0;i<n_sample; ++i)
	{
		if (wt[i]==0 && abs(D.norm_dp[i] - 1) < abs(D.norm_dp[i] - C2[1].Mean))
		{
			if (S.len<=500 || D.n_cnv_pos[i]==0 || D.n_cnv_neg[i]==0 || D.norm_cnv_pos[i] <0.7 || D.norm_cnv_pos[i]>1.3)
			{
				Y.push_back(D.norm_dp[i]);
				P_Y.push_back(D.norm_cnv_pos[i]);
				N_Y.push_back(D.norm_cnv_neg[i]);
			}
		}
	}

	C2[0].estimate(Y);
	C2[0].Alpha = ((double)Y.size()) / (X.size() + Y.size());
	C2[1].Alpha = ((double)X.size()) / (X.size() + Y.size());

	// Fit 1 or 2 comp
	G.bic[0] = BIC(D.norm_dp, C1);
	G.bic[1] = BIC(D.norm_dp, C2);
//	double be2 = BayesError(C2);

	if (G.bic[1] < G.bic[0])
	{
		G.dp_flag = true;

		// ToDo : 3-component, construct 3-components from the beginning based on depth
		G.p_overlap = BayesError(C2);
		copyComps(G.Comps, C2);
	}
	if (X2.size()>0)
	{
		C3[0].estimate(Y);
		C3[1].estimate(X);
		C3[2].estimate(X2);

		C3[0].Alpha = ((double)Y.size()) / (X.size() + Y.size() + X2.size());
		C3[1].Alpha = ((double)X.size()) / (X.size() + Y.size() + X2.size());
		C3[2].Alpha = ((double)X2.size()) / (X.size() + Y.size() + X2.size());
		G.bic[2] = BIC(D.norm_dp, C3);

		if (G.bic[2] < G.bic[1] && G.bic[2] < G.bic[0])
		{
			G.dp_flag = true;

			G.p_overlap = BayesError(C3);
			copyComps(G.Comps, C3);
		}
	}

	vector<Gaussian> P1(1);
	vector<Gaussian> P2(2);
	if (P_X.size() > 0)
	{

		P1[0].estimate(D.cnv_pos);
		P1[0].Alpha = 1;

		P2[0].estimate(P_Y);
		P2[1].estimate(P_X);
		P2[0].Alpha = (double) P_Y.size() / (P_X.size() + P_Y.size());
		P2[1].Alpha = (double) P_X.size() / (P_X.size() + P_Y.size());

		if (BIC(D.norm_cnv_pos, P2) < BIC(D.norm_cnv_pos, P1))
		{
			G.pos_flag = true;
			G.info = ";POS=(" + to_string((int)(P2[0].Mean*S.len)) + "," + to_string((int)(P2[1].Mean*S.len)) + ")";
		}
	}

	vector<Gaussian> N1(1);
	vector<Gaussian> N2(2);
	if (N_X.size() > 0)
	{
		N1[0].estimate(D.cnv_neg);
		N1[0].Alpha = 1;

		N2[0].estimate(N_Y);
		N2[1].estimate(N_X);
		N2[0].Alpha = (double) N_Y.size() / (N_X.size() + N_Y.size());
		N2[1].Alpha = (double) N_X.size() / (N_X.size() + N_Y.size());

		if (BIC(D.norm_cnv_neg, N2) < BIC(D.norm_cnv_neg, N1))
		{
			G.neg_flag = true;
			G.info += ";NEG=(" + to_string((int)(N2[0].Mean*S.len)) + "," + to_string((int)(N2[1].Mean*S.len)) + ")";
		}
	}

	G.ac = 0;
	G.ns = 0;
	G.b_pass = false;
	G.b_dump = true;

	if (G.dp_flag || G.pos_flag || G.neg_flag)
	{
		bool b_hom = false;
		for(int i=0; i<n_sample; ++i)
		{
			int pos_gt = -1;
			int neg_gt = -1;
			G.cn[i] = -1;
			if (G.dp_flag)
			{
				G.cn[i] = assign(D.norm_dp[i], G.Comps);
			}

			if (G.pos_flag && D.n_cnv_pos[i] > 1 && D.norm_dp[i] < 0.9)
			{
				double p0 = P2[0].pdf(D.norm_cnv_pos[i]);
				double p1 = P2[1].pdf(D.norm_cnv_pos[i]);

				if (p1 > 5.0 * p0)
				{
					pos_gt = 1;
				}
				else if (p0 > p1)
				{
					pos_gt = 0;
				}
			}

			if (G.neg_flag && D.n_cnv_neg[i] > 1 && D.norm_dp[i] < 0.9)
			{
				double p0 = N2[0].pdf(D.norm_cnv_neg[i]);
				double p1 = N2[1].pdf(D.norm_cnv_neg[i]);
				if (p1 > 5.0 * p0)
				{
					neg_gt = 1;
				}
				else if (p0 > p1)
				{
					neg_gt = 0;
				}
			}

			if (G.cn[i] == 0 || (pos_gt==1 && neg_gt==1 && D.norm_dp[i] < 0.2))
			{
				G.gt[i] = 2;
				G.ac+=2;
				G.ns++;
				b_hom = true;
			}
			else if (G.cn[i] == 1 || (pos_gt==1 && neg_gt==1))
			{
				G.gt[i] = 1;
				G.ac++;
				G.ns++;
			}
			else if (G.cn[1] == 2 || (pos_gt==0 && neg_gt==0))
			{
				G.gt[i] = 0;
				G.ns++;
			}
		}

		G.b_dump = false;

		if (G.ac > 0.95*G.ns && b_hom == false)
		{
			G.b_pass=  false;
		}
		else if (G.ns > 0.5*n_sample && G.ac > 0 && G.ac < (G.ns*2.0))
		{
			G.b_pass = true;
		}
	}
	else
	{
		G.b_dump = true;
	}
	
}


void gtype::call_del(sv &S, svdata& D, svgeno &G, vector<double> &avg_isz, vector<double>&std_isz, vector<double> &wt)
{
	int n_sample=D.n;

	vector<Gaussian> P1(1);
	vector<Gaussian> P2(2);

	vector<Gaussian> N1(1);
	vector<Gaussian> N2(2);

//	cerr << S.svtype << "_" << S.chr << ":" << S.pos << "-" << S.end << endl;

	vector<int> posneg_gt(n_sample, -1);

	for(int i=0;i<n_sample; ++i)
	{
		if (D.norm_dp[i]>1.4)
			wt[i] = 0;
	}

	if (S.len > 150)  
	{
		vector<int> pos_i(n_sample, 0);
		vector<int> neg_i(n_sample, 0);;

		double sum_pos0 = 0, sum_pos1 = 0;
		double sum_neg0 = 0, sum_neg1 = 0;
		double sumsq_pos0 = 0, sumsq_pos1 = 0;
		double sumsq_neg0 = 0, sumsq_neg1 = 0;

		int n_pos0 = 0, n_pos1 = 0;
		int n_neg0 = 0, n_neg1 = 0;;

		for(int i=0;i<n_sample; ++i)
		{
			double max_err = S.len / 10.0 ; 
			if (max_err > 1000) 
			{
				max_err = 1000;
			}
			else if (max_err < 75)
			{
				max_err = 75;
			}

			double isz_min = S.len + avg_isz[i] - max_err;
			double isz_max = S.len + avg_isz[i] + max_err;
			
			if (D.n_cnv_pos[i] > 2 && D.cnv_pos[i] > isz_min && D.cnv_pos[i] < isz_max && D.norm_dp[i]<0.75)
			{
				pos_i[i] = 1;
				sum_pos1 += D.n_cnv_pos[i] * D.cnv_pos[i];
				n_pos1 += D.n_cnv_pos[i];
			}
			else
			{
				pos_i[i] = 0;
				sum_pos0 += D.cnv_pos[i];
				n_pos0 ++;
			}
			if (D.n_cnv_neg[i] > 2 && D.cnv_neg[i] > isz_min && D.cnv_neg[i] < isz_max && D.norm_dp[i]<0.75)
			{
				sum_neg1 += D.n_cnv_neg[i] * D.cnv_neg[i];
				neg_i[i] = 1;
				n_neg1 += D.n_cnv_neg[i];
			}
			else
			{
				sum_neg0 += D.cnv_neg[i];
				neg_i[i] = 0;
				n_neg0 ++;
			}
			if (pos_i[i] && neg_i[i])
			{
				wt[i] += (D.n_cnv_pos[i] + D.n_cnv_neg[i]) / 2.0;
			}
		}
		if (n_pos1 > 0)
		{
			P1[0].estimate(D.cnv_pos);
			P1[0].Alpha = 1;

			double m0 = sum_pos0 / n_pos0;

			// pseudocount;
			sum_pos1 += S.len + S.len + 800;
			n_pos1 +=2; 
			double m1 = sum_pos1 / n_pos1;

			for(int i=0; i<n_sample; ++i)
			{
				if (pos_i[i])
				{
					sumsq_pos1 += D.n_cnv_pos[i] * (D.cnv_pos[i]-m1) * (D.cnv_pos[i]-m1) ;
				}
				else
				{
					sumsq_pos0 += (D.cnv_pos[i] -m0 ) * (D.cnv_pos[i] - m0);
				}
			}
			sumsq_pos1 += (S.len+300 - m1)*(S.len+300 -m1) + (S.len+500-m1)*(S.len+500-m1);

			double stdev0 = sqrt(sumsq_pos0 / n_pos0);
			double stdev1 = sqrt(sumsq_pos1 / n_pos1);

			P2[0].set(m0 , stdev0) ;
			P2[0].Alpha = (double)n_pos0 / (double)(n_pos0 + n_pos1) ;

			P2[1].set(m1 , stdev1) ;
			P2[1].Alpha = (double)n_pos1 / (double)(n_pos0 + n_pos1) ;

			G.info += ";POS=" + to_string(P2[0].Mean) + "(" + to_string(P2[0].Stdev) + ")"  + "," + to_string(P2[1].Mean) + "(" + to_string(P2[1].Stdev) + ")";

			double bic_p1 = BIC(D.cnv_pos, P1, wt);
			double bic_p2 = BIC(D.cnv_pos, P2, wt);

//			double be_p1 = BayesError(P2);
			if (bic_p1 > bic_p2)
			{
				G.pos_flag = true;
			}
			G.info += ";BIC_P=" + to_string(bic_p1) +  "," + to_string(bic_p2);
//			G.info += ";BE_P=" + to_string(be_p1);
		}
		if (n_neg1 > 0)
		{
			N1[0].estimate(D.cnv_neg);
			N1[0].Alpha = 1;

			double m0 = sum_neg0 / n_neg0;

			// pseudocount;
			sum_neg1 += S.len + S.len + 800;
			n_neg1 +=2; 
			double m1 = sum_neg1 / n_neg1;

			for(int i=0; i<n_sample; ++i)
			{
				if (neg_i[i])
				{
					sumsq_neg1 += D.n_cnv_neg[i] * (D.cnv_neg[i]-m1) * (D.cnv_neg[i]-m1) ;
				}
				else
				{
					sumsq_neg0 += (D.cnv_neg[i] -m0 ) * (D.cnv_neg[i] - m0);
				}
			}
			sumsq_neg1 += (S.len+300 - m1)*(S.len+300 -m1) + (S.len+500-m1)*(S.len+500-m1);
			double stdev0 = sqrt(sumsq_neg0 / n_neg0);
			double stdev1 = sqrt(sumsq_neg1 / n_neg1);

			N2[0].set(m0 , stdev0) ;
			N2[0].Alpha = (double)n_neg0 / (double)(n_neg0 + n_neg1) ;

			N2[1].set(m1 , stdev1) ;
			N2[1].Alpha = (double)n_neg1 / (double)(n_neg0 + n_neg1) ;

			G.info += ";NEG=" + to_string(N2[0].Mean) + "(" + to_string(N2[0].Stdev) + ")"  + "," + to_string(N2[1].Mean) + "(" + to_string(N2[1].Stdev) + ")";

			double bic_n1 = BIC(D.cnv_neg, N1, wt);
			double bic_n2 = BIC(D.cnv_neg, N2, wt);
			if (bic_n1 > bic_n2)
			{
				G.neg_flag = true;
			}
			G.info += ";BIC_N=" + to_string(bic_n1) +  "," + to_string(bic_n2);
		}

		if (G.pos_flag && G.neg_flag)
		{
			for(int i=0; i<n_sample; ++i)
			{
				double p0 = P2[0].pdf(D.cnv_pos[i]);
				double p1 = P2[1].pdf(D.cnv_pos[i]);

				double n0 = N2[0].pdf(D.cnv_neg[i]);
				double n1 = N2[1].pdf(D.cnv_neg[i]);

				if ((p1 > 2*p0 && D.n_cnv_pos[i]>1) && (n1 > 2*n0 && D.n_cnv_neg[i] >1) && D.norm_dp[i] < 0.85) // Add depth constraint
				{
					posneg_gt[i] = 1;
				}
				else if (p0 >p1 && n0 >n1)
				{
					posneg_gt[i] = 0;
					// default : missing
				}
			}
		}
	} //S.len > 150

	double be2, be3;

	vector<Gaussian> C1(1); // 1-component model
	vector<Gaussian> C2(2); // 2-component model
	vector<Gaussian> C3(3); // 3-component model
	
	// Fit Gaussian mixture models with 1, 2, and 3 components, compare BIC
	C1[0].estimate(D.norm_dp);
	C1[0].Alpha = 1;
	
	double min_bic = DBL_MAX;

	// P_Overlap
	G.p_overlap = 1;
	G.b_biallelic = true; //deletions 
	
	// One-component model
	G.bic[0] = BIC(D.norm_dp, C1, wt);
	min_bic = G.bic[0]; 
	
	// Two-component model
	C2[0].set(1, 0.01);
	C2[1].set(0.5, 0.01);
	C2[0].Alpha = C2[1].Alpha = 0.5;
	EM(D.norm_dp, wt, C2);
	//fit(D.norm_dp, C2);

	copyComps(G.Comps, C1);

	G.bic[1] = BIC(D.norm_dp, C2, wt);
	G.dp_flag = false;
	be2 = BayesError(C2);

	if (G.bic[1] < min_bic || be2<0.2)
	{
		if (ordered(C2) && C2[0].Mean > 0.8 && C2[0].Mean < 1.3 && C2[1].Mean<0.7 && C2[1].Mean > 0.3)
		{
			G.dp_flag = true;
			min_bic = G.bic[1];
			copyComps(G.Comps, C2);
			G.p_overlap = be2;
		}
	}

	// Three-component model
	C3[0].set(1, 0.01);
	C3[1].set(0.5, 0.01);
	C3[2].set(0.05, 0.01);
	C3[0].Alpha = C3[1].Alpha = C3[2].Alpha = 1.0/3.0;
	EM(D.norm_dp, wt, C3);

	//fit(D.norm_dp, C3);

	G.bic[2] = BIC(D.norm_dp, C3, wt);
	be3 = BayesError(C3);

	if ((G.bic[2] < min_bic || be3<0.2) && ordered(C3) && be3<G.p_overlap+0.2)
	{
		min_bic = G.bic[2];
		G.dp_flag = true;
		copyComps(G.Comps, C3);
		G.p_overlap = be3;
	}

	vector<int> dp_gt(n_sample, -1);

	if ((G.dp_flag && G.p_overlap<0.35) || (G.dp_flag && (G.pos_flag||G.neg_flag) && G.p_overlap < 0.4 ) || (G.dp_flag && G.pos_flag && G.neg_flag && G.p_overlap<0.45))
	{
		for(int i=0;i<n_sample;++i)
		{
			switch(assign(D.norm_dp[i], G.Comps))
			{
			case 0:
				dp_gt[i] = 2;
				break;
			case 1:
				dp_gt[i] = 1;
				break;
			case 2:
				dp_gt[i] = 0;
				break;
			}
		}
	}
	
	for(int i=0;i<n_sample;++i)
	{
		if (dp_gt[i] == 2 || (posneg_gt[i] == 1 && D.norm_dp[i]<0.2))
		{
			G.gt[i] = 2;
			G.cn[i] = 0;
			G.ac +=2;
			G.ns +=1;
		}
		else if (dp_gt[i] == 1 || posneg_gt[i] == 1)
		{
			G.gt[i] = 1;
			G.cn[i] = 1;
			G.ac += 1;
			G.ns +=1;
		}
		else if (dp_gt[i] == 0 || posneg_gt[i] == 0)
		{
			G.gt[i] = 0;
			G.cn[i] = 2;

			G.ns += 1;
		}
	}

	if (G.dp_flag || G.pos_flag || G.neg_flag)
	{
		G.b_dump = false;
		if (G.Comps.size() == 2 && G.Comps[1].Alpha > 0.8)
		{
			// This is unlikely. If HET is highly common, there should be HOMALT.
			G.b_pass = false;
		}
		else if ( (G.ac>0) &&  (G.ac < (2*G.ns)) && (G.ns>(n_sample*0.5)))
		{
			G.b_pass = true;
		}
	}
	else
	{
		G.b_dump = true;
	}
}

void gtype::call_cnv(sv &S, svdata& D, svgeno& G, vector<double> &avg_isz, vector<double> &std_isz, vector<double> &wt)
{
	int n_sample = D.n;

	vector<Gaussian> P1(1);
	vector<Gaussian> P2(2);
	
	vector<Gaussian> N1(1);
	vector<Gaussian> N2(2);
	
	vector<int> posneg_gt(n_sample, -1);
	for(int i=0;i<n_sample; ++i)
	{
		if (D.norm_dp[i]<0.6)
			wt[i] = 0;
	}



	if (S.len > 500)
	{
		vector<int> pos_i(n_sample, 0);
		vector<int> neg_i(n_sample, 0);;
		
		double sum_pos0 = 0, sum_pos1 = 0;
		double sum_neg0 = 0, sum_neg1 = 0;
		double sumsq_pos0 = 0, sumsq_pos1 = 0;
		double sumsq_neg0 = 0, sumsq_neg1 = 0;
		
		int n_pos0 = 0, n_pos1 = 0;
		int n_neg0 = 0, n_neg1 = 0;;
		
		for(int i=0;i<n_sample; ++i)
		{
			double max_err = S.len / 10.0 ;
			if (max_err > 1000)
			{
				max_err = 1000;
			}
			else if (max_err < 75)
			{
				max_err = 75;
			}
			
			double isz_min = S.len - avg_isz[i] - max_err;
			double isz_max = S.len - avg_isz[i] + max_err;
			
			if (D.n_cnv_pos[i] > 2 && D.cnv_pos[i] > isz_min && D.cnv_pos[i] < isz_max && D.norm_dp[i]>1.3)
			{
				pos_i[i] = 1;
				sum_pos1 += D.n_cnv_pos[i] * D.cnv_pos[i];
				n_pos1 += D.n_cnv_pos[i];
			}
			else
			{
				pos_i[i] = 0;
				sum_pos0 += D.cnv_pos[i];
				n_pos0 ++;
			}
			if (D.n_cnv_neg[i] > 2 && D.cnv_neg[i] > isz_min && D.cnv_neg[i] < isz_max && D.norm_dp[i]>1.3)
			{
				sum_neg1 += D.n_cnv_neg[i] * D.cnv_neg[i];
				neg_i[i] = 1;
				n_neg1 += D.n_cnv_neg[i];
			}
			else
			{
				sum_neg0 += D.cnv_neg[i];
				neg_i[i] = 0;
				n_neg0 ++;
			}
			if (pos_i[i] && neg_i[i] )
			{
				wt[i] += (D.n_cnv_pos[i] + D.n_cnv_neg[i]) / 2.0;
			}
		}
		if (n_pos1 > 0)
		{
			P1[0].estimate(D.cnv_pos);
			P1[0].Alpha = 1;
			
			double m0 = sum_pos0 / n_pos0;
			
			// pseudocount;
			sum_pos1 += S.len + S.len - 800;
			n_pos1 +=2;
			double m1 = sum_pos1 / n_pos1;
			
			for(int i=0; i<n_sample; ++i)
			{
				if (pos_i[i])
				{
					sumsq_pos1 += D.n_cnv_pos[i] * (D.cnv_pos[i]-m1) * (D.cnv_pos[i]-m1) ;
				}
				else
				{
					sumsq_pos0 += (D.cnv_pos[i] -m0 ) * (D.cnv_pos[i] - m0);
				}
			}
			sumsq_pos1 += (S.len-300 - m1)*(S.len-300 -m1) + (S.len-500-m1)*(S.len-500-m1);
			
			double stdev0 = sqrt(sumsq_pos0 / n_pos0);
			double stdev1 = sqrt(sumsq_pos1 / n_pos1);
			
			P2[0].set(m0 , stdev0) ;
			P2[0].Alpha = (double)n_pos0 / (double)(n_pos0 + n_pos1) ;
			
			P2[1].set(m1 , stdev1) ;
			P2[1].Alpha = (double)n_pos1 / (double)(n_pos0 + n_pos1) ;
			
			G.info += ";POS=" + to_string(P2[0].Mean) + "(" + to_string(P2[0].Stdev) + ")"  + "," + to_string(P2[1].Mean) + "(" + to_string(P2[1].Stdev) + ")";
			
			double bic_p1 = BIC(D.cnv_pos, P1, wt);
			double bic_p2 = BIC(D.cnv_pos, P2, wt);
			double be_p1 = BayesError(P2);
			if (bic_p1 > bic_p2)
			{
				G.pos_flag = true;
			}
			G.info += ";BIC_P=" + to_string(bic_p1) +  "," + to_string(bic_p2);
			G.info += ";BE_P=" + to_string(be_p1);
		}
		if (n_neg1 > 0)
		{
			N1[0].estimate(D.cnv_neg);
			N1[0].Alpha = 1;
			
			double m0 = sum_neg0 / n_neg0;
			
			// pseudocount;
			sum_neg1 += S.len + S.len - 800;
			n_neg1 +=2;
			double m1 = sum_neg1 / n_neg1;
			
			for(int i=0; i<n_sample; ++i)
			{
				if (neg_i[i])
				{
					sumsq_neg1 += D.n_cnv_neg[i] * (D.cnv_neg[i]-m1) * (D.cnv_neg[i]-m1) ;
				}
				else
				{
					sumsq_neg0 += (D.cnv_neg[i] -m0 ) * (D.cnv_neg[i] - m0);
				}
			}
			sumsq_neg1 += (S.len-300 - m1)*(S.len-300 -m1) + (S.len-500-m1)*(S.len-500-m1);
			double stdev0 = sqrt(sumsq_neg0 / n_neg0);
			double stdev1 = sqrt(sumsq_neg1 / n_neg1);
			
			N2[0].set(m0 , stdev0) ;
			N2[0].Alpha = (double)n_neg0 / (double)(n_neg0 + n_neg1) ;
			
			N2[1].set(m1 , stdev1) ;
			N2[1].Alpha = (double)n_neg1 / (double)(n_neg0 + n_neg1) ;
			
			G.info += ";NEG=" + to_string(N2[0].Mean) + "(" + to_string(N2[0].Stdev) + ")"  + "," + to_string(N2[1].Mean) + "(" + to_string(N2[1].Stdev) + ")";
			

			double bic_n1 = BIC(D.cnv_neg, N1, wt);
			double bic_n2 = BIC(D.cnv_neg, N2, wt);
			double be_n1 = BayesError(N2);
		
			if (bic_n1 > bic_n2)
			{
				G.neg_flag = true;
			}

			G.info += ";BIC_N=" + to_string(bic_n1) +  "," + to_string(bic_n2);
			G.info += ";BE_N=" + to_string(be_n1);
		}
		
		if (G.pos_flag && G.neg_flag)
		{
			for(int i=0; i<n_sample; ++i)
			{
				double p0 = P2[0].pdf(D.cnv_pos[i]);
				double p1 = P2[1].pdf(D.cnv_pos[i]);
				
				double n0 = N2[0].pdf(D.cnv_neg[i]);
				double n1 = N2[1].pdf(D.cnv_neg[i]);
				
				if (((p1 > 2*p0 && D.n_cnv_pos[i]>1) && (n1 > 2*n0 && D.n_cnv_neg[i] >1)) && D.norm_dp[i]>1.2)
				{
					posneg_gt[i] = 1;
				}
				else if (p0 >p1 && n0 >n1)
				{
					posneg_gt[i] = 0;
					// default : missing
				}
			}
		}
	} //S.len > 500

	vector<Gaussian> C(1);
	double min_BIC = DBL_MAX;

	unsigned n_comp = 1;

	C[0].Mean = 1;
	C[0].Stdev = stdev(D.norm_dp, C[0].Mean);
	G.bic[0] = BIC(D.norm_dp, C, wt);
	min_BIC = G.bic[0];
	copyComps(G.Comps, C);

	bool keep_going = true;
	G.p_overlap = 1;
	for(int i=1;i<10 && keep_going; ++i)
	{
		Gaussian *pG = new Gaussian;
		C.push_back(*pG);

		for(int j=0; j<=i; ++j)
		{
			C[j].Mean = 0.5 * j + 1.0;
			C[j].Stdev = 0.02;
			C[j].Alpha = 1.0/(i+1);
		}
//		C[0].Stdev = 0.1; // More variance in CN2 only
		
		EM(D.norm_dp, wt, C);
		//fit(D.norm_dp, C);

		G.bic[i] = BIC(D.norm_dp, C, wt);
		double be = BayesError(C);
		
		keep_going = false;

		if (G.bic[i] < min_BIC && C[i].Alpha > 0.5/(double)n_sample && r_ordered(C) && be < G.p_overlap+0.2) // not significantly increase Bayes error
		{
			n_comp = i+1;

			min_BIC = G.bic[i];

			copyComps(G.Comps, C);

			G.dp_flag = true;
			G.p_overlap = be;

			keep_going = true;
		}
	}

	vector<int> dp_cn(n_sample, -1);
	
	if ((G.dp_flag && G.p_overlap<0.35) || (G.dp_flag && (G.pos_flag||G.neg_flag) && G.p_overlap < 0.4) || (G.dp_flag && G.pos_flag && G.neg_flag && G.p_overlap<0.45))
	{
		for(int i=0;i<n_sample;++i)
		{
			dp_cn[i] = assign(D.norm_dp[i], G.Comps);
		}
	}

	G.b_biallelic = (n_comp < 4);

	for(int i=0;i<n_sample;++i)
	{
		if (G.b_biallelic)
		{
			if (dp_cn[i] == 4 || (posneg_gt[i] == 1 && D.norm_dp[i]>1.8))
			{
				G.gt[i] = 2;
				G.cn[i] = 4;
				G.ac +=2;
				G.ns +=1;
			}
			else if (dp_cn[i] == 3 || posneg_gt[i] == 1)
			{
				G.gt[i] = 1;
				G.cn[i] = 3;
				G.ac += 1;
				G.ns +=1;
			}
			else if (dp_cn[i] == 2 || posneg_gt[i] == 0)
			{
				G.gt[i] = 0;
				G.cn[i] = 2;
				G.ns += 1;
			}
			else
			{
				G.gt[i] = -1;
				G.cn[i] = -1;
			}
		}
		else
		{
			if (G.cn[i]>=0)
			{
				G.ns++;
				if (G.cn[i] >  2)
				{
					G.ac++;
				}
			}
		}
	}

	if (G.dp_flag || G.pos_flag || G.neg_flag)
	{
		G.b_dump = false;

		if (n_comp ==2 && G.Comps[1].Alpha > 0.8)
		{
			// This is unlikely
			G.b_pass = false;
		}
		else if ( (G.ac>0) &&  (G.ac < (2*G.ns)) && (G.ns>(n_sample*0.5)))
		{
			G.b_pass = true;
		}
	}
	else
	{
		G.b_dump = true;
	}
}


void gtype::call_inv(sv &S, svdata& D, svgeno& G, vector<double>& avg_isz, vector<double> &std_isz)
{
	int n_sample = D.n;

	double bic2, bic3;
	double be2, be3;
	double max_overlap = 0.3;

//	vector<Gaussian> C1(1); // 1-component model
	vector<Gaussian> C2(2); // 2-component model
	vector<Gaussian> C3(3); // 3-component model
	
	// For each candidate region, run EM
	// Fit Gaussian mixture models with 1, 2, and 3 components, compare BIC
//	C1[0].estimate(D.norm_readcount);
//	C1[0].Alpha = 1;
	
//	double min_bic = DBL_MAX;

	// P_Overlap
	G.p_overlap = 1;
	G.b_biallelic = true; //deletions 
	
	// One-component model
//	bic1 = BIC(D.norm_readcount, C1);
//	min_bic = bic1;
//	copyComps(G.Comps, C1);
	
	// Two-component model
	C2[0].set(1.1, 0.1);
	C2[1].set(0.6, 0.1);
	C2[0].Alpha = C2[1].Alpha = 0.5;
	EM(D.norm_readcount, C2);
	bic2 = BIC(D.norm_readcount, C2);

//	if (bic2 < min_bic && C2[0].Mean > C2[1].Mean && C2[0].Mean > 0.75 && C2[1].Mean > 0.3 && C2[1].Mean < 0.75)
//	if (bic2 < min_bic) // now we do not compare mean
		be2 = BayesError(C2);
		if (be2 < max_overlap && C2[0].Alpha > 2.5/(n_sample+2.0) && C2[1].Alpha > 2.5/(n_sample+2.0))
		{
//			G.p_overlap = be2;
//			min_bic = bic2;
			G.dp_flag = true;
//			copyComps(G.Comps, C2);
		}

	// Three-component model
	C3[0].set(1, 0.1);
	C3[1].set(0.5, 0.1);
	C3[2].set(0, 0.1);
	C3[0].Alpha = C3[1].Alpha = C3[2].Alpha = 1.0/3.0;
	EM(D.norm_dp, C3);
	bic3 = BIC(D.norm_dp, C3);

//		min_bic = bic3;
		be3 = BayesError(C3);
		if (be3 < max_overlap && C3[0].Alpha > 2.5/(n_sample+2.0) && C3[1].Alpha > 2.5/(n_sample+2.0) && C3[2].Alpha > 2.5/(n_sample+2.0))
		{
//			G.p_overlap = be3;
			G.dp_flag = true;
//			copyComps(G.Comps, C3);
		}

	vector<int> pos_gt(n_sample, 0);
	vector<int> neg_gt(n_sample, 0);

	for(int i=0;i<n_sample;++i)
	{
//			if (D.norm_cnv_pos[i] > 0.5 && D.norm_cnv_pos[i] < 2 && D.n_cnv_pos[i] > 2)
		double min_isz = S.len - (std_isz[i]*3);
		double max_isz = S.len + (std_isz[i]*3);


		if (D.n_cnv_pos[i]>0 && D.n_cnv_neg[i]>0 && D.cnv_pos[i] > min_isz && D.cnv_pos[i] < max_isz && D.cnv_neg[i] < max_isz && D.cnv_neg[i] > min_isz)
		{
			pos_gt[i] = 1;
			neg_gt[i] = 1;
			G.pos_flag = true;
			G.neg_flag = true;
		}

		if (D.n_inv_pos[i]>1 && D.inv_pos[i] > min_isz && D.inv_pos[i] < max_isz)
		{
			pos_gt[i] = 1;
			G.pos_flag = true;
		}
		if (D.n_inv_neg[i]>1 && D.inv_neg[i] > min_isz && D.inv_neg[i] < max_isz)
		{
			neg_gt[i] = 1;
			G.neg_flag = true;
		}
	}

//	if (!G.dp_flag)
	{
		// Depth-based clustering failed, force 2 or 3 component model
		if (bic2 < bic3)
			copyComps(G.Comps, C2);
		else
			copyComps(G.Comps, C3);
		G.p_overlap = BayesError(G.Comps);
	}

	for(int i=0;i<n_sample;++i)
	{
		G.cn[i] = assign(D.norm_readcount[i], G.Comps);
		G.gt[i] = -1;
		
		// To prevent high depth ones are classified as double deletion due to long tail
		if ((D.norm_readcount[i]*2 - G.cn[i]) > 1 || (D.norm_readcount[i]*2 - G.cn[i]) < -1)
		{
			G.cn[i] = -1;
		}
		switch(G.cn[i])
		{
			case -1:
				if (pos_gt[i] && neg_gt[i])
				{
					if ( 0.2 <= D.norm_readcount[i] && D.norm_readcount[i]< 0.9)
					{
						G.gt[i] = 1;
						G.cn[i] = 1;
					}
					else if (D.norm_readcount[i] < 0.2)
					{
						G.gt[i] = 2;
						G.cn[i] = 0;
					}
				}
				else if (pos_gt[i] || neg_gt[i])
				{
					if ( 0.3 <= D.norm_readcount[i] && D.norm_readcount[i] <= 0.7)
					{
						G.gt[i] = 1;
						G.cn[i] = 1;
					}
					else if (D.norm_readcount[i] < 0.2)
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
				if (pos_gt[i] && neg_gt[i] && D.norm_readcount[i] < 0.9)
				{
					G.gt[i] = 1;
					G.cn[i] = 1;
				}
				else if ((pos_gt[i] || neg_gt[i]) && D.norm_readcount[i] <=0.7 )
				{
					G.gt[i] = 1;
					G.cn[i] = 1;
				}
				break;
		}
		if (G.gt[i] >=0 )
		{
			G.ac += G.gt[i];
			G.ns += 1;
		}
	}
	if (G.dp_flag || G.pos_flag || G.neg_flag)
	{
		G.b_dump = false;
		G.b_pass = false;
		if ((G.ac>0 &&  G.ac < (2*n_sample) && G.ns>(n_sample *0.5)) && ((G.dp_flag && G.p_overlap<0.2) ||  (G.p_overlap < 0.3 && (G.neg_flag || G.pos_flag)) || (G.p_overlap < 0.4 && G.neg_flag && G.pos_flag)))
		{
			G.b_pass = true;
		}
	}
	else
	{
		G.b_pass = false;
		G.b_dump = true;
	}

}

// EM with weights
void gtype::EM(vector<double>& x, vector<double> &w, vector<Gaussian>& Comps)
{
	unsigned n_sample = (unsigned) x.size();
	unsigned n_comp = (unsigned) Comps.size();
	unsigned n_iter = 10;
	
	unsigned p_count= 1;
	double p_val[n_comp];
	int zeroidx = -1;
	
	for(unsigned i=0; i<n_comp; ++i)
	{
		p_val[i] = Comps[i].Mean;
		/*
		if (Comps[i].Mean < 0.1)
		{
			zeroidx = i;
		}
		*/
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
				/*
				if (zeroidx == (int)m )
				{
					pr[m] *= 2.0;
				}
				*/
				sum_p += pr[m];
			}
			
			if (sum_p > 1e-30) // if the value is an outlier, exclude it from calculations
			{
				for(unsigned m=0;m<n_comp;++m)
				{
					pr[m] /= sum_p;
					sum[m] += pr[m] * x[j] * w[j];
					sum_err[m] += pr[m] * w[j] * (x[j] - Comps[m].Mean)*(x[j]-Comps[m].Mean);
					sum_pr[m] += pr[m] * w[j];
				}
			}
		}
		
		double sumsum = 0;
		// Add pseudo-count values
		for(unsigned m=0; m<n_comp; ++m)
		{
			sum[m] += p_val[m] * p_count;
			
			sum_err[m] += (p_val[m] - Comps[m].Mean)*(p_val[m]-Comps[m].Mean) * p_count;
			sum_pr[m] += p_count;
			sumsum += sum_pr[m];
		}
		
		// M step
		for(unsigned m=0;m<n_comp;++m)
		{
			Comps[m].Mean = sum[m]/sum_pr[m];
			Comps[m].Stdev = sqrt(sum_err[m] / sum_pr[m]) ;
			Comps[m].Alpha = sum_pr[m] / sumsum;
			if (Comps[m].Stdev < 1e-10)
				Comps[m].Stdev = 0.005;

//			cerr << "\t(" << Comps[m].Mean << "," << Comps[m].Stdev  << " : " << Comps[m].Alpha << ")" ;
			/*
			if (zeroidx == (int)m)
			{
				Comps[m].Mean = 0;
			}
			*/
		}
//		cerr << endl;
	}
}


// EM for general CNVs without weights
void gtype::EM(vector<double>& x, vector<Gaussian>& Comps)
{
	unsigned n_sample = (unsigned) x.size();
	unsigned n_comp = (unsigned) Comps.size();
	unsigned n_iter = 10;
	
	unsigned p_count= 2;
	double p_val[n_comp];
	int zeroidx = -1;
	
	for(unsigned i=0; i<n_comp; ++i)
	{
		p_val[i] = Comps[i].Mean;
		if (Comps[i].Mean < 0.1)
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
			
			if (sum_p > 1e-12) // if the value is an outlier, exclude it from calculations
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
			Comps[m].Stdev = sqrt(sum_err[m] / sum_pr[m]) ;
			Comps[m].Alpha = sum_pr[m] /( n_sample + n_comp*p_count);
			//Comps[m].Alpha = (sum_pr[m] - p_count)/( n_sample);

			// cerr << "\t(" << Comps[m].Mean << "," << Comps[m].Stdev  << " : " << Comps[m].Alpha << ")" ;
			if (zeroidx == (int)m)
			{
				Comps[m].Mean = 0;
			}
		}
		//cerr << endl;
	}
}


// EM with means fixed multiples of base
void gtype::fit(vector<double>& x, vector<Gaussian>& Comps)
{
	int n_sample = (int) x.size();
	int n_comp = (int) Comps.size();
	int n_iter = 10;
	
	int p_count= 2;
	double p_val[n_comp];
	int zeroidx = -1;
	
	int factor[n_comp];
	double base = 0.5;

	for(int i=0; i<n_comp; ++i)
	{
		p_val[i] = Comps[i].Mean;
		factor[i] = round(Comps[i].Mean / base);
		if (p_val[i] <0.1) // zero
		{
			zeroidx = i;
		}
	}
	
	for(int i=0; i<n_iter; ++i)
	{
		vector<double> sum (n_comp,0);
		vector<double> sum_pr (n_comp,0);
		vector<double> sum_err (n_comp,0);
		
		// E step
		double sum_b = 0;
		for(int j=0; j<n_sample; ++j)
		{
			double sum_p = 0;
			vector<double> pr(n_comp, 0);
			for(int m=0;m<n_comp;++m)
			{
				pr[m] = Comps[m].Alpha * Comps[m].pdf(x[j]);
				sum_p += pr[m];
			}
			
			if (sum_p > 1e-12) // if the value is an outlier, exclude it from calculations
			{
				for(int m=0;m<n_comp;++m)
				{
					pr[m] /= sum_p;
					if (m!= zeroidx)
					{
						sum_b +=  pr[m] *x[j] / (double)factor[m];
					}
//					sum[m] += pr[m] * x[j];
					sum_err[m] += pr[m] * (x[j] - Comps[m].Mean)*(x[j]-Comps[m].Mean);
					sum_pr[m] += pr[m];
				}
			}
		}
		
		// Add pseudo-count values
		for(int m=0; m<n_comp; ++m)
		{
//			sum[m] += p_val[m] * p_count;
			sum_err[m] += (p_val[m] - Comps[m].Mean)*(p_val[m]-Comps[m].Mean) * p_count;
			sum_pr[m] += p_count;
		}
		
		// M step
		double  sum_nonzero = 0;
		for(int m=0;m<n_comp;++m)
		{
//			Comps[m].Mean = sum[m] / sum_pr[m];
			if (m!=zeroidx)
				sum_nonzero += sum_pr[m];

			Comps[m].Stdev = sqrt(sum_err[m] / (sum_pr[m] ));
			Comps[m].Alpha = sum_pr[m] /( n_sample + n_comp*p_count);
		}
		base = sum_b / sum_nonzero;

		for(int m=0;m<n_comp;++m)
		{
			Comps[m].Mean = base * factor[m];
			//cerr << "\t(" << Comps[m].Mean << "," << Comps[m].Stdev  << " : " << Comps[m].Alpha << ")" ;
		}
//		cerr << endl;
	}

	double sum_p = 0;
	for(int m=0;m<n_comp;++m)
	{
		Comps[m].Alpha  -= p_count / ((double)n_sample + n_comp * p_count);
		if (Comps[m].Alpha < 0) Comps[m].Alpha = 0;
		sum_p += Comps[m].Alpha;
	}
//	cerr << "Final: ";
	for(int m=0;m<n_comp;++m)
	{
		Comps[m].Alpha  = Comps[m].Alpha / sum_p;
//		cerr << "\t(" << Comps[m].Mean << "," << Comps[m].Stdev  << " : " << Comps[m].Alpha << ")" ;
	}
//	cerr << endl;

}


