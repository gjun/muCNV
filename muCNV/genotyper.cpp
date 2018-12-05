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

#include "genotyper.h"

bool ordered(std::vector<Gaussian> &C)
{
	for(unsigned i=0; i<C.size()-1;++i)
	{
		if (C[i].Mean<=C[i+1].Mean)
			return false;
	}
	return true;
}

bool r_ordered(std::vector<Gaussian> &C)
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

int Genotyper::assign(double x, std::vector<Gaussian> &C)
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
	
	if (max_R > 0.1)
	{
		return -1;
	}
	
	return ret;
}

void svgeno::print(sv &S, svdata &D, std::string &ln, std::vector<double>& wt)
{
	ln = std::to_string(S.chrnum);
	ln += "\t" + std::to_string(S.pos) + "\t" + svTypeName(S.svtype) + "_" + std::to_string(S.chrnum) + ":" + std::to_string(S.pos) + "-" + std::to_string(S.end) + "\t.\t<" + svTypeName(S.svtype) + ">\t.\t";

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
	
    ln += "SVTYPE=" + std::string(svTypeName(S.svtype)) + ";END=" + std::to_string(S.end) + ";SVLEN=" + std::to_string(S.len) + ";AC=" + std::to_string(ac) + ";NS=" + std::to_string(ns) + ";AF=";
	if (ns>0)
		ln+=std::to_string((double)ac/(double)(2.0*ns));
	else
		ln+="0";
	
	ln+= info;

	bool bic_flag = false ;
	if (Comps.size()>0)
	{
		ln+= ";NCLUS=" + std::to_string(Comps.size()) + ";P_OVERLAP=" + std::to_string(p_overlap);

		ln+= ";BIC_DP="+ std::to_string(bic[0]) + "," + std::to_string(bic[1]) + "," + std::to_string(bic[2]);
		for(unsigned i=3;i<Comps.size(); ++i)
		{
			ln += "," + std::to_string(bic[i]);
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
		ln += ";MEAN=" + std::to_string(Comps[0].Mean);

		for(unsigned i=1;i<Comps.size(); ++i)
			ln += "," + std::to_string(Comps[i].Mean);

		ln += ";STDEV=" + std::to_string(Comps[0].Stdev);
		for(unsigned i=1;i<Comps.size(); ++i)
			ln += "," + std::to_string(Comps[i].Stdev);

		ln += ";FRAC=" + std::to_string(Comps[0].Alpha);
		for(unsigned i=1;i<Comps.size(); ++i)
			ln += "," + std::to_string(Comps[i].Alpha);
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
				ln += std::to_string(cn[i]) + ":" + std::to_string(D.norm_readcount[i]).substr(0,4) + ":" + std::to_string((int)D.dp[i]) + ":" +  std::to_string(D.n_inv_pos[i]) + "," + std::to_string((int)D.inv_pos[i]) + ":" + std::to_string(D.n_inv_neg[i]) + "," + std::to_string((int)D.inv_neg[i]);
			}
			else
			{
				ln += ".:" + std::to_string(D.norm_readcount[i]).substr(0,4) + ":" + std::to_string((int)D.dp[i]) + ":" +  std::to_string(D.n_inv_pos[i]) + "," + std::to_string((int)D.inv_pos[i]) + ":" + std::to_string(D.n_inv_neg[i]) + "," + std::to_string((int)D.inv_neg[i]) ;
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
				ln += std::to_string(cn[i])+":";
			}
			ln += std::to_string(D.norm_dp[i]).substr(0,4) + ":" + std::to_string((int)D.dp[i]) + ":" + std::to_string((int)D.n_cnv_pos[i]) + "," + std::to_string((int)D.cnv_pos[i]) + ":" + std::to_string((int)D.n_cnv_neg[i]) + "," + std::to_string((int)D.cnv_neg[i]) + ":" + std::to_string((int)wt[i]); 
		}
	}
}

void Genotyper::copyComps(std::vector<Gaussian> &C, std::vector<Gaussian> &C0)
{
	C.clear();
	C.resize(C0.size());
	for(unsigned j=0;j<C.size(); ++j)
	{
		C[j].set(C0[j].Mean, C0[j].Stdev);
		C[j].Alpha = C0[j].Alpha;
	}
}

void Genotyper::call_del(sv &S, svdata& D, svgeno &G, std::vector<double> &avg_isz, std::vector<double>&std_isz, std::vector<double> &wt)
{
	int n_sample=D.n;

	std::vector<Gaussian> P1(1);
	std::vector<Gaussian> P2(2);

	std::vector<Gaussian> N1(1);
	std::vector<Gaussian> N2(2);

//	cerr << S.svtype << "_" << S.chr << ":" << S.pos << "-" << S.end << endl;

	std::vector<int> posneg_gt(n_sample, -1);

	for(int i=0;i<n_sample; ++i)
	{
		if (D.norm_dp[i]>1.4)
			wt[i] = 0;
	}

	if (S.len > 150)  
	{
		std::vector<int> pos_i(n_sample, 0);
		std::vector<int> neg_i(n_sample, 0);;

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

			G.info += ";POS=" + std::to_string(P2[0].Mean) + "(" + std::to_string(P2[0].Stdev) + ")"  + "," + std::to_string(P2[1].Mean) + "(" + std::to_string(P2[1].Stdev) + ")";

			double bic_p1 = BIC(D.cnv_pos, P1, wt);
			double bic_p2 = BIC(D.cnv_pos, P2, wt);

//			double be_p1 = BayesError(P2);
			if (bic_p1 > bic_p2)
			{
				G.pos_flag = true;
			}
			G.info += ";BIC_P=" + std::to_string(bic_p1) +  "," + std::to_string(bic_p2);
//			G.info += ";BE_P=" + std::to_string(be_p1);
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

			G.info += ";NEG=" + std::to_string(N2[0].Mean) + "(" + std::to_string(N2[0].Stdev) + ")"  + "," + std::to_string(N2[1].Mean) + "(" + std::to_string(N2[1].Stdev) + ")";

			double bic_n1 = BIC(D.cnv_neg, N1, wt);
			double bic_n2 = BIC(D.cnv_neg, N2, wt);
			if (bic_n1 > bic_n2)
			{
				G.neg_flag = true;
			}
			G.info += ";BIC_N=" + std::to_string(bic_n1) +  "," + std::to_string(bic_n2);
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

	std::vector<Gaussian> C1(1); // 1-component model
	std::vector<Gaussian> C2(2); // 2-component model
	std::vector<Gaussian> C3(3); // 3-component model
	
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

	std::vector<int> dp_gt(n_sample, -1);

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

void Genotyper::call_cnv(sv &S, svdata& D, svgeno& G, std::vector<double> &avg_isz, std::vector<double> &std_isz, std::vector<double> &wt)
{
	int n_sample = D.n;

	std::vector<Gaussian> P1(1);
	std::vector<Gaussian> P2(2);
	
	std::vector<Gaussian> N1(1);
	std::vector<Gaussian> N2(2);
	
	std::vector<int> posneg_gt(n_sample, -1);
	for(int i=0;i<n_sample; ++i)
	{
		if (D.norm_dp[i]<0.6)
			wt[i] = 0;
	}



	if (S.len > 500)
	{
		std::vector<int> pos_i(n_sample, 0);
		std::vector<int> neg_i(n_sample, 0);;
		
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
			
			G.info += ";POS=" + std::to_string(P2[0].Mean) + "(" + std::to_string(P2[0].Stdev) + ")"  + "," + std::to_string(P2[1].Mean) + "(" + std::to_string(P2[1].Stdev) + ")";
			
			double bic_p1 = BIC(D.cnv_pos, P1, wt);
			double bic_p2 = BIC(D.cnv_pos, P2, wt);
			double be_p1 = BayesError(P2);
			if (bic_p1 > bic_p2)
			{
				G.pos_flag = true;
			}
			G.info += ";BIC_P=" + std::to_string(bic_p1) +  "," + std::to_string(bic_p2);
			G.info += ";BE_P=" + std::to_string(be_p1);
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
			
			G.info += ";NEG=" + std::to_string(N2[0].Mean) + "(" + std::to_string(N2[0].Stdev) + ")"  + "," + std::to_string(N2[1].Mean) + "(" + std::to_string(N2[1].Stdev) + ")";
			

			double bic_n1 = BIC(D.cnv_neg, N1, wt);
			double bic_n2 = BIC(D.cnv_neg, N2, wt);
			double be_n1 = BayesError(N2);
		
			if (bic_n1 > bic_n2)
			{
				G.neg_flag = true;
			}

			G.info += ";BIC_N=" + std::to_string(bic_n1) +  "," + std::to_string(bic_n2);
			G.info += ";BE_N=" + std::to_string(be_n1);
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

	std::vector<Gaussian> C(1);
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

	std::vector<int> dp_cn(n_sample, -1);
	
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



