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

void svgeno::print(sv &S, svdata &D, string &ln)
{
	ln = S.chr;
	ln += "\t" + to_string(S.pos) + "\t" + S.svtype + "_" + S.chr + ":" + to_string(S.pos) + "-" + to_string(S.end) + "\t.\t<" + S.svtype + ">\t.\t";

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
	
	ln += "SVTYPE=" + S.svtype + ";END=" + to_string(S.end) + ";SVLEN=" + to_string(S.len) + ";AC=" + to_string(ac) + ";NS=" + to_string(ns) + ";AF=";
	if (ns>0)
		ln+=to_string((double)ac/(double)(2.0*ns));
	else
		ln+="0";
	
	ln+= info;
	ln+= ";NCLUS=" + to_string(Comps.size()) + ";P_OVERLAP=" + to_string(p_overlap);


	bool bic_flag = false ;
	ln+= ";bic="+ to_string(bic[0]) + "," + to_string(bic[1]) + "," + to_string(bic[2]);
	for(unsigned i=3;i<Comps.size(); ++i)
	{
		ln += "," + to_string(bic[i]);
		if (bic[i]<bic[0])
		{
			bic_flag = true;
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

	ln += ";MEAN=" + to_string(Comps[0].Mean);
	for(unsigned i=1;i<Comps.size(); ++i)
		ln += "," + to_string(Comps[i].Mean);

	ln += ";STDEV=" + to_string(Comps[0].Stdev);
	for(unsigned i=1;i<Comps.size(); ++i)
		ln += "," + to_string(Comps[i].Stdev);

	ln += ";FRAC=" + to_string(Comps[0].Alpha);
	for(unsigned i=1;i<Comps.size(); ++i)
		ln += "," + to_string(Comps[i].Alpha);

	if (b_biallelic)
		ln += "\tGT:CN:ND:DP:FP:FN";
	else
		ln += "\tCN:ND:DP:FP:FN";
	
	for (int i=0; i<n_sample; ++i)
	{
		if (b_biallelic)
		{
			switch(gt[i])
			{
				case 0:
					ln += "\t0/0:";
					break;
				case 1:
					ln += "\t0/1:";
					break;
				case 2:
					ln += "\t1/1:";
					break;
				default:
					ln += "\t.:";
					break;
			}
		}
		else
		{
			ln += "\t";
		}

		if (S.svtype == "INV")
		{
			if (cn[i]<0)
			{
				ln += to_string(cn[i]) + ":" + to_string(D.norm_readcount[i]).substr(0,4) + ":" + to_string((int)D.dp[i]) + ":" + to_string(D.n_cnv_pos[i]) + "," + to_string((int)D.cnv_pos[i]) + ":" + to_string(D.n_cnv_neg[i]) + "," + to_string((int)D.cnv_neg[i]) + ":" +  to_string(D.n_inv_pos[i]) + "," + to_string((int)D.inv_pos[i]) + ":" + to_string(D.n_inv_neg[i]) + "," + to_string((int)D.inv_neg[i]);
			}
			else
			{
				ln += ".:" + to_string(D.norm_readcount[i]).substr(0,4) + ":" + to_string((int)D.dp[i]) + ":" + to_string(D.n_cnv_pos[i]) + "," + to_string((int)D.cnv_pos[i]) + ":" + to_string(D.n_cnv_neg[i]) + "," + to_string((int)D.cnv_neg[i]) + ":" +  to_string(D.n_inv_pos[i]) + "," + to_string((int)D.inv_pos[i]) + ":" + to_string(D.n_inv_neg[i]) + "," + to_string((int)D.inv_neg[i]) ;
			}
		}
		else
		{
			if (cn[i] <0)
			{
				ln +=  ".:" + to_string(D.norm_dp[i]).substr(0,4) + ":" + to_string((int)D.dp[i]) + ":" + to_string((int)D.n_cnv_pos[i]) + "," + to_string((int)D.cnv_pos[i]) + ":" + to_string((int)D.n_cnv_neg[i]) + "," + to_string((int)D.cnv_neg[i]) +  ":" +  to_string(D.n_inv_pos[i]) + "," + to_string((int)D.inv_pos[i]) + ":" + to_string(D.n_inv_neg[i]) + "," + to_string((int)D.inv_neg[i]) + ":" + to_string((int)D.n_isz[i]); 
			}
			else
			{
				ln += to_string(cn[i]) + ":" + to_string(D.norm_dp[i]).substr(0,4) + ":" + to_string((int)D.dp[i]) + ":" + to_string((int)D.n_cnv_pos[i]) + "," + to_string((int)D.cnv_pos[i]) + ":" + to_string((int)D.n_cnv_neg[i]) + "," + to_string((int)D.cnv_neg[i])  + ":" +  to_string(D.n_inv_pos[i]) + "," + to_string((int)D.inv_pos[i]) + ":" + to_string(D.n_inv_neg[i]) + "," + to_string((int)D.inv_neg[i]) + ":" + to_string((int)D.n_isz[i]); 
			}
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


void gtype::call_del(sv &S, svdata& D, svgeno &G, vector<double> &avg_isz, vector<double>&std_isz, vector<double> &wt)
{
	int n_sample=D.n;

	vector<Gaussian> P1(1);
	vector<Gaussian> P2(2);

	vector<Gaussian> N1(1);
	vector<Gaussian> N2(2);

	vector<int> posneg_gt(n_sample, -1);

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
			double max_err = S.len / 7.0 ; 
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

/*
			double bic_p1 = BIC(D.cnv_pos, P1);
			double bic_p2 = BIC(D.cnv_pos, P2);
			double be_p1 = BayesError(P2);
			if (bic_p1 > bic_p2)
				*/
			{
				G.pos_flag = true;
			}
//			G.info += ";BIC_P=" + to_string(bic_p1) +  "," + to_string(bic_p2);
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

/*
			if (BIC(D.cnv_neg, N1) > BIC(D.cnv_neg, N2))*/
			{
				G.neg_flag = true;
			}
		}

		if (G.pos_flag && G.neg_flag)
		{
			for(int i=0; i<n_sample; ++i)
			{
				double p0 = P2[0].pdf(D.cnv_pos[i]);
				double p1 = P2[1].pdf(D.cnv_pos[i]);

				double n0 = N2[0].pdf(D.cnv_neg[i]);
				double n1 = N2[1].pdf(D.cnv_neg[i]);

				if ((p1 > 2*p0 && D.n_cnv_pos[i]>1) || (n1 > 2*n0 && D.n_cnv_neg[i] >1))
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
	G.bic[0] = BIC(D.norm_dp, C1);
	min_bic = G.bic[0]; 
	
	// Two-component model
	C2[0].set(1, 0.1);
	C2[1].set(0.5, 0.05);
	C2[0].Alpha = C2[1].Alpha = 0.5;
	EM(D.norm_dp, wt, C2);
	//fit(D.norm_dp, C2);

	copyComps(G.Comps, C1);

	G.bic[1] = BIC(D.norm_dp, C2);
	G.dp_flag = false;
	be2 = BayesError(C2);

	if (G.bic[1] < min_bic)
	{
		if (ordered(C2) && C2[0].Mean > 0.8 && C2[0].Mean < 1.3 && C2[1].Mean<0.7 && C2[1].Mean > 0.3)
			G.dp_flag = true;
		min_bic = G.bic[1];
		copyComps(G.Comps, C2);
		G.p_overlap = be2;
	}

	// Three-component model
	C3[0].set(1, 0.1);
	C3[1].set(0.5, 0.05);
	C3[2].set(0, 0.05);
	C3[0].Alpha = C3[1].Alpha = C3[2].Alpha = 1.0/3.0;
	EM(D.norm_dp, wt, C3);

	//fit(D.norm_dp, C3);

	G.bic[2] = BIC(D.norm_dp, C3);
	be3 = BayesError(C3);

	if (G.bic[2] < min_bic && ordered(C3) && be3<0.4 && be3<be2+0.2) 
	{
		min_bic = G.bic[2];
		G.dp_flag = true;
		copyComps(G.Comps, C3);
		G.p_overlap = be3;
	}

	vector<int> dp_gt(n_sample, -1);

	if ((G.dp_flag && G.p_overlap<0.2) || (G.dp_flag && (G.pos_flag||G.neg_flag) && G.p_overlap < 0.3 ) || (G.dp_flag && G.pos_flag && G.neg_flag && G.p_overlap<0.4))
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
		if (dp_gt[i] == 2)
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
		if ( (G.ac>0) &&  (G.ac < (2*G.ns)) && (G.ns>(n_sample*0.5)))
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

	vector<Gaussian> C(1);
	double min_BIC = DBL_MAX;

	unsigned n_comp = 1;

	C[0].Mean = 1;
	C[0].Stdev = stdev(D.norm_dp, C[0].Mean);
	G.bic[0] = BIC(D.norm_dp, C);
	min_BIC = G.bic[0];
	copyComps(G.Comps, C);

/*
	vector<int> cnt(20,0);
	// Count CN for each bin
	for(int i=0;i<n_sample;++i)
	{
		int k = round(D.norm_dp[i]*2.0);
		if (k>4) 
		{
			if (k>10) k=10;
			cnt[k]++;
		}
	}

	vector<size_t> idx = sort_indexes(cnt);
	for (int i=1; i<=2 ; ++i)
	{
		Gaussian *pG = new Gaussian;
		C.push_back(*pG);
		for(int j=0;j<=i;++j)
		{
			C[j].Mean = 1.0 + j*0.5;
			C[j].Stdev = 0.1;
			C[j].Alpha = 1.0/(i+1);
		}
		EM(D.norm_dp, C);
		double bic = BIC(D.norm_dp, C);
		G.bic[i] = bic;

		bool mean_flag = false;
		for(int j=0;j<=i; ++j)
		{
			if (abs(C[j].Mean  - (1.0 + j*0.5)) > 0.2)
			{
				mean_flag = true;
			}
		}

		if (bic < min_BIC && BayesError(C)<max_overlap && !mean_flag)
		{
			n_comp = i+1;
			min_BIC = bic;
			copyComps(G.Comps, C);
			G.dp_flag = true;
		}
	}

	if (n_comp == 3)
	{
		for(int i=0; i<10 && cnt[idx[i]]>0; ++i)
		{
			Gaussian *pG = new Gaussian;
			C.push_back(*pG);

			for(unsigned j=0;j<3; ++j)
			{
				C[j].Mean = 1.0 + j*0.5;
				C[j].Stdev = 0.1;
				C[j].Alpha = 1.0/(i+4);
			}
			for(unsigned j=0; j<=i; ++j)
			{
				C[j+3].Mean = (double)idx[j]*0.5;
				C[j+3].Stdev = 0.1;
				C[j+3].Alpha = 1.0/(i+4);
			}
			EM(D.norm_dp, C);

			double bic = BIC(D.norm_dp, C); 
			G.bic[i+3] = bic;

			bool mean_flag = false;
			for(int j=0;j<3; ++j)
			{
				if (abs(C[j].Mean  - (1.0 + j*0.5)) > 0.2)
				{
					mean_flag = true;
				}
			}

			for(int j=0;j<=i; ++j)
			{
				if (abs(C[j+3].Mean  - idx[j]*0.5 ) > 0.2)
				{
					mean_flag = true;
				}
			}


			if (bic < min_BIC && BayesError(C)<max_overlap && !mean_flag) // More stringent for higher copy numbers
			{
				n_comp = i+4;
				min_BIC = bic;
				copyComps(G.Comps, C);
				G.dp_flag = true;
			}
		}
	}
	*/

	bool keep_going = true;
	for(int i=1;i<10 && keep_going; ++i)
	{
		Gaussian *pG = new Gaussian;
		C.push_back(*pG);

		for(int j=0; j<=i; ++j)
		{
			C[j].Mean = 0.5 * j + 1.0;
			C[j].Stdev = 0.05;
			C[j].Alpha = 1.0/(i+1);
		}

		EM(D.norm_dp, C);
		//fit(D.norm_dp, C);

		G.bic[i] = BIC(D.norm_dp, C);
		keep_going = false;
		if (G.bic[i] < min_BIC && C[i].Alpha > 0.5/(double)n_sample) // not significantly increase Bayes error
		{
			n_comp = i+1;

			min_BIC = G.bic[i];
			copyComps(G.Comps, C);
			G.dp_flag = true;
			keep_going = true;
		}
	}

	G.p_overlap = BayesError(G.Comps);

	vector<int> pos_gt(n_sample, 0);
	vector<int> neg_gt(n_sample, 0);

	if (S.len >500)// 100?
	{
		for(int i=0;i<n_sample;++i)
		{
			double min_isz = S.len - avg_isz[i] - (std_isz[i]*2);
			double max_isz = S.len + (std_isz[i]*2);

			if (D.n_cnv_pos[i]>1 && D.n_cnv_neg[i]>1 && D.cnv_pos[i] > min_isz && D.cnv_pos[i] < max_isz && D.cnv_neg[i] < max_isz && D.cnv_neg[i] > min_isz)
			{
				pos_gt[i] = 1;
				neg_gt[i] = 1;
				G.pos_flag = true;
				G.neg_flag = true;
			}

			if (D.n_cnv_pos[i]>2 && D.cnv_pos[i] > min_isz && D.cnv_pos[i] < max_isz)
			{
				pos_gt[i] = 1;
				G.pos_flag = true;
			}
			if (D.n_cnv_neg[i]>2 && D.cnv_neg[i] > min_isz && D.cnv_neg[i] < max_isz)
			{
				neg_gt[i] = 1;
				G.neg_flag = true;
			}
		}
	}

	// if any of three clustering meets criteria
//	if (dp_flag || pos_flag || neg_flag)
	{
		for(int i=0;i<n_sample;++i)
		{
			G.gt[i] = -1;
			G.cn[i] = assign(D.norm_dp[i], G.Comps);

/*
			if (G.cn[i] == -1  && D.norm_dp[i]>1.25)
			{
				//if ((pos_gt[i] && neg_gt[i]) || (pos_gt[i] && D.n_cnv_pos[i]>3) || (neg_gt[i] && D.n_cnv_neg[i]>3))
				if ((pos_gt[i] && neg_gt[i]))
				{
					G.cn[i] = round(D.norm_dp[i]*2.0);
				}
			}
			if (G.cn[i] == 2 && D.norm_dp[i] > 1.25) 
			{
				//if ((pos_gt[i] && neg_gt[i]) || (pos_gt[i] && D.n_cnv_pos[i]>3) || (neg_gt[i] && D.n_cnv_neg[i]>3))
				if ((pos_gt[i] && neg_gt[i]))
				{
					G.cn[i] = round(D.norm_dp[i]*2.0);
				//	G.cn[i] = 3;  // guess there's a duplication
				}
			}
			*/
		}
		
		G.b_biallelic = (n_comp < 4);
		for(int i=0;i<n_sample;++i)
		{
			if( G.b_biallelic)
			{
				switch(G.cn[i])
				{
				case 2:
					G.gt[i] = 0;
					G.ns ++;
					break;
				case 3:
					G.gt[i] = 1;
					G.ns ++;
					G.ac += 1;
					break;
				case 4:
					G.gt[i] = 2;
					G.ns ++;
					G.ac += 2;
					break;
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
	}
	if (G.dp_flag || G.pos_flag || G.neg_flag)
	{
		G.b_dump = false;
		//if ((G.ac>0 &&  G.ac < (2*n_sample) && G.ns>(n_sample *0.5)) && ((G.dp_flag && G.p_overlap<0.25) ||  (G.p_overlap < 0.35 && (G.neg_flag || G.pos_flag)) || (G.p_overlap < 0.45 && G.neg_flag && G.pos_flag)))
		if (n_comp ==2 && G.Comps[1].Alpha > 0.9)
		{
			// This is unlikely
			G.b_pass = false;
		}
		else if ((G.ac>0 &&  G.ac < (2*n_sample) && G.ns>(n_sample *0.5)) && (n_comp>1) && ((G.dp_flag && G.p_overlap<0.2) || (G.p_overlap<0.3 && G.dp_flag && G.neg_flag && G.pos_flag)))
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
	
	unsigned p_count= 10;
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
				sum_p += pr[m] * w[j];
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



// EM for general CNVs
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
			
			if (sum_p > 1e-30) // if the value is an outlier, exclude it from calculations
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


