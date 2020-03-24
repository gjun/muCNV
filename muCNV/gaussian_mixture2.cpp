//
//  gaussian_mixture2.cpp
//  Project: muCNV
//  Created by: Goo Jun (goo.jun@gmail.com)

#include "gaussian_mixture2.h"
#include <math.h>
#include <iomanip>
#include <sstream>

GaussianMixture2::GaussianMixture2(std::vector<double> &m, std::vector<double> &s)
{
	if (m.size() != s.size())
	{
		std::cerr << "Error, mean and std. dev vector sizes are different" << std::endl;
		exit(1);
	}
	n_comp = (int) m.size();

	Comps.resize(n_comp);
	for(int i=0; i< n_comp; ++i)
	{
		Comps[i].Mean[0] = m[i];
		Comps[i].Mean[1] = m[i];

		Comps[i].Cov[0] = s[i];
		Comps[i].Cov[3] = s[i];

		Comps[i].Alpha = 1.0/n_comp;
		if (Comps[i].Mean[0] + Comps[i].Mean[1] < 0.1)
			zeroidx = i;
	}
}

GaussianMixture2::GaussianMixture2(std::vector<double> &m1, std::vector<double> &m2, std::vector<double> &s)
{
	if (m1.size() != s.size() || m1.size() != m2.size())
	{
		std::cerr << "Error, mean and std. dev vector sizes are different" << std::endl;
		exit(1);
	}
	n_comp = (int) m1.size();

	Comps.resize(n_comp);
	for(int i=0; i< n_comp; ++i)
	{
		Comps[i].Mean[0] = m1[i];
		Comps[i].Mean[1] = m2[i];

		Comps[i].Cov[0] = s[i];
		Comps[i].Cov[3] = s[i];

		Comps[i].Alpha = 1.0/n_comp;
		if (Comps[i].Mean[0] + Comps[i].Mean[1] < 0.1)
			zeroidx = i;
	}
}

GaussianMixture2& GaussianMixture2::operator = (const GaussianMixture2& gmix)
{
	n_comp = (int) gmix.n_comp;

	Comps.resize(gmix.n_comp);

	for(int i=0; i<n_comp; ++i)
	{
		Comps[i].Mean[0] = gmix.Comps[i].Mean[0];
		Comps[i].Mean[1] = gmix.Comps[i].Mean[1];

		for(int j=0; j<4; ++j)
		{
			Comps[i].Cov[j] = gmix.Comps[i].Cov[j];
		}
        Comps[i].update();
		Comps[i].Alpha = gmix.Comps[i].Alpha;
	}
	bic = gmix.bic;
    
	p_overlap = gmix.p_overlap;
	zeroidx = gmix.zeroidx;

	return *this;
}


void GaussianMixture2::estimate(std::vector<double> &x, std::vector<double> &y, std::vector<int> &z, int n)
{
    // TODO: clean current components?

    n_comp = n;
    int n_sample = (int)x.size();
    if (n_sample != (int)y.size())
    {
        std::cerr<< "Error, GaussianMixture2 has different size of X and Y" << std::endl;
        exit(1);
    }        
    if (n_sample != (int)z.size())
    {
        std::cerr<< "Error, GaussianMixture2 has different size of data and label" << std::endl;
        exit(1);
    }

	std::vector<double> sum_x (n_comp, 0);
   	std::vector<double> sum_y (n_comp, 0);
	std::vector<double> sum_xx (n_comp, 0);
	std::vector<double> sum_xy (n_comp, 0);
	std::vector<double> sum_yy (n_comp, 0);
    std::vector<int> cnt(n_comp, 0);
	
	for(int i=0; i<n_sample; ++i)
	{
		if (z[i] >=0 && z[i] < n_comp)
		{
			sum_x[z[i]] +=x[i];
			sum_y[z[i]] +=y[i];
			sum_xx[z[i]] +=x[i]*x[i];
			sum_xy[z[i]] +=x[i]*y[i];
			sum_yy[z[i]] +=y[i]*y[i];
			cnt[z[i]]++;
		}
	}
	
    for(int i=0; i<n_comp; ++i)
	{
        if (cnt[i]>0)
	    {
            Comps[i].Mean[0] = sum_x[i]/(double)cnt[i];
            Comps[i].Mean[1] = sum_y[i]/(double)cnt[i];
            Comps[i].Cov[0] = sum_xx[i]/(double)cnt[i] - (Comps[i].Mean[0]*Comps[i].Mean[0]);
            Comps[i].Cov[1] = sum_xy[i]/(double)cnt[i] - (Comps[i].Mean[0]*Comps[i].Mean[1]);
            Comps[i].Cov[2] = Comps[i].Cov[1];
            Comps[i].Cov[3] = sum_yy[i]/(double)cnt[i] - (Comps[i].Mean[1]*Comps[i].Mean[1]);
            Comps[i].N = cnt[i];
            Comps[i].update();

    	}
    }

    bic = BIC(x, y);
}


// 2-D EM for general without weights
void GaussianMixture2::KM2(std::vector<double>& x, std::vector<double> &y, bool b_mahalanobis)
{
	int n_sample = (int) x.size();
	int n_iter = 15;

	// pseudo-counts
	int p_count= 5;
	double p_val[n_comp][2];

	if (n_comp == 1)
	{
		Comps[0].estimate(x, y);
		Comps[0].Alpha = 1;
		bic = BIC(x, y);
		p_overlap = 0;
		return;
	}

	for(int i=0; i<n_comp; ++i)
	{
		p_val[i][0] = Comps[i].Mean[0];
		p_val[i][1] = Comps[i].Mean[1];
	}
	std::vector<double> member (n_sample, 0);

	for(int i=0; i<n_iter; ++i)
	{

		std::vector<double> sum_x (n_comp,0);
		std::vector<double> sum_y (n_comp,0);
		std::vector<int> cnt (n_comp,0);

		// Hard assignment of samples to centroids (E step)
		for(int j=0; j<n_sample; ++j)
		{
			double min_dist = DBL_MAX;
			for(int m=0; m<n_comp; ++m)
			{
				double dx = (x[j] - Comps[m].Mean[0]);
				double dy = (y[j] - Comps[m].Mean[1]);

				double dist = 0;

				// Mahalanobis dist
				if (b_mahalanobis)
				{
					dist = (dx * Comps[m].Prc[0] + dy * Comps[m].Prc[2]) * dx + (dx * Comps[m].Prc[1] + dy * Comps[m].Prc[3]) * dy;
				}
				else
				{
					dist = dx*dx + dy*dy;
				}
				if (dist < min_dist)
				{
					min_dist = dist;
					member[j] = m;
				}
			}
			sum_x[member[j]] += x[j];
			sum_y[member[j]] += y[j];
			cnt[member[j]] ++;
		}

		// M step
		// Update means
		for(int m=0;m<n_comp;++m)
		{
			Comps[m].Mean[0] = (sum_x[m] + p_val[m][0]*p_count )/ (cnt[m] + p_count);
			Comps[m].Mean[1] = (sum_y[m] + p_val[m][1]*p_count )/ (cnt[m] + p_count);
		}

		if (b_mahalanobis)
		{
			std::vector<double> sum_xx(n_comp, 0);
			std::vector<double> sum_yy(n_comp, 0);
			std::vector<double> sum_xy(n_comp, 0);
			for(int j=0; j<n_sample; ++j)
			{
				double dx = Comps[member[j]].Mean[0] - x[j];
				double dy = Comps[member[j]].Mean[1] - y[j];
				sum_xx[member[j]] += dx * dx;
				sum_yy[member[j]] += dy * dy;
				sum_xy[member[j]] += dx * dy;
			}

			for(int m=0;m<n_comp;++m)
			{
				// Wishart prior : S = [0.01 0; 0 0.01]
				Comps[m].Cov[0] = (sum_xx[m] + 0.01*p_count) / (cnt[m] + p_count - 1);
				Comps[m].Cov[1] = sum_xy[m] / (cnt[m] + p_count - 1);
				Comps[m].Cov[2] = Comps[m].Cov[1];
				Comps[m].Cov[3] = (sum_yy[m] + 0.01*p_count) / (cnt[m] + p_count - 1);

				Comps[m].update();

				Comps[m].Alpha = (double)cnt[m] / n_sample;
			}
		}
	}

	if (!b_mahalanobis)
	{
		std::vector<double> sum_xx(n_comp, 0);
		std::vector<double> sum_yy(n_comp, 0);
		std::vector<double> sum_xy(n_comp, 0);
		std::vector<int> cnt(n_comp, 0);
		for(int j=0; j<n_sample; ++j)
		{
			double dx = Comps[member[j]].Mean[0] - x[j];
			double dy = Comps[member[j]].Mean[1] - y[j];
			sum_xx[member[j]] += dx * dx;
			sum_yy[member[j]] += dy * dy;
			sum_xy[member[j]] += dx * dy;

			cnt[member[j]] ++;
		}

		for(int m=0;m<n_comp;++m)
		{
			// Wishart prior : S = [0.01 0; 0 0.01]
			Comps[m].Cov[0] = (sum_xx[m] + 0.01*p_count) / (cnt[m] + p_count - 1);
			Comps[m].Cov[1] = sum_xy[m] / (cnt[m] + p_count - 1);
			Comps[m].Cov[2] = Comps[m].Cov[1];
			Comps[m].Cov[3] = (sum_yy[m] + 0.01*p_count) / (cnt[m] + p_count - 1);

			Comps[m].update();

			Comps[m].Alpha = (double)cnt[m] / n_sample;
		}
	}

	double llk = 0;
	for(int j=0; j<n_sample; ++j)
	{
		double l = 0;
		for(int m=0; m<n_comp; ++m)
		{
			l += Comps[m].Alpha * Comps[m].pdf(x[j], y[j]);
		}
		if (l>0 && !isnan(l))
		{
			llk += log(l);
		}
	}

	bic = -2.0 * llk + 5*n_comp*log(n_sample);
	p_overlap = BayesError();
	//    std::cerr << "BIC: " << bic << ", P_OVERLAP: " << p_overlap << std::endl;
}

// 2-D EM for general without weights
void GaussianMixture2::EM2(std::vector<double>& x, std::vector<double> &y)
{
	// Let's not consider half-normal distribution -- for now

	int n_sample = (int) x.size();
	int n_iter = 15;

	// pseudo-counts
	double p_count= 0.1;

	double p_val[n_comp][2];

	if (n_comp == 1)
	{
		Comps[0].estimate(x, y);
		Comps[0].Alpha = 1;
		bic = BIC(x, y);
		p_overlap = 0;
		return;
	}

	for(int i=0; i<n_comp; ++i)
	{
		p_val[i][0] = Comps[i].Mean[0];
		p_val[i][1] = Comps[i].Mean[1];

		if (Comps[i].Mean[0] < 0.01 && Comps[i].Mean[1] <0.01)
		{
			zeroidx = i;
		}
	}

	// pseudo-means
	//    double p_val[3] = {1.0, 0.5, 0};

	for(int i=0; i<n_iter; ++i)
	{
		std::vector<double> sum_x (n_comp,0);
		std::vector<double> sum_y (n_comp,0);

		std::vector<double> sum_e_xx (n_comp,0);
		std::vector<double> sum_e_xy (n_comp,0);
		std::vector<double> sum_e_yy (n_comp,0);

		std::vector<double> sum_pr (n_comp,0);

		double pr[n_comp][n_sample];
		bool b_include[n_sample];

		// E step
		for(int j=0; j<n_sample; ++j)
		{
			b_include[j] = true;
			double sum_p = 0;
			for(int m=0;m<n_comp;++m)
			{
				pr[m][j] = Comps[m].Alpha * Comps[m].pdf(x[j], y[j]);
                pr[m][j] *= 4.0;

				sum_p += pr[m][j];
			}
			if (sum_p < 1e-10) b_include[j] = false;

			if (b_include[j]) // if the value is an outlier, exclude it from calculations
			{
				for(int m=0;m<n_comp;++m)
				{
					pr[m][j] /= sum_p;
					sum_x[m] += pr[m][j] * x[j];
					sum_y[m] += pr[m][j] * y[j];

					sum_pr[m] += pr[m][j];
				}
			}
		}
		double sumsum = 0;
		// Add pseudo-count values
		for(int m=0; m<n_comp; ++m)
		{
            if ( m!= zeroidx)
            {
                sum_x[m] += p_val[m][0] * p_count;
                sum_y[m] += p_val[m][1] * p_count;
                sum_pr[m] += p_count;
            }
			sumsum += sum_pr[m];
		}

		// M step

		// Update means
		for(int m=0;m<n_comp;++m)
		{
			if (m != zeroidx)
			{
				Comps[m].Mean[0] = sum_x[m]/sum_pr[m];
				Comps[m].Mean[1] = sum_y[m]/sum_pr[m];
			}
		}

		// Update covariance
		for(int j=0; j<n_sample; ++j)
		{
			if (b_include[j]) // if the value is an outlier, exclude it from calculations
			{
				for(int m=0;m<n_comp;++m)
				{
					double ex = x[j] - Comps[m].Mean[0];
					double ey = y[j] - Comps[m].Mean[1];

					sum_e_xx[m] += pr[m][j] * ex * ex;
					sum_e_xy[m] += pr[m][j] * ex * ey;
					sum_e_yy[m] += pr[m][j] * ey * ey;
				}
			}
		}

		for(int m=0;m<n_comp;++m)
		{
            double ex = p_val[m][0] - Comps[m].Mean[0];
            double ey = p_val[m][1] - Comps[m].Mean[1];

            sum_e_xx[m] += ex * ex;
            sum_e_xy[m] += ex * ey;
            sum_e_yy[m] += ey * ey;

            // Wishart prior : S = [0.01 0; 0 0.01]
            if (sum_pr[m]>1e-10)
            {
                Comps[m].Cov[0] = (sum_e_xx[m]  + 0.01 * p_count) / sum_pr[m];
                Comps[m].Cov[1] = Comps[m].Cov[2] = sum_e_xy[m] / sum_pr[m];
                Comps[m].Cov[3] = (sum_e_yy[m] + 0.01 * p_count) / sum_pr[m];
            }
            Comps[m].update();

            Comps[m].Alpha = sum_pr[m] / sumsum;
		}
		//print();
	}
	double llk = 0;
	for(int j=0; j<n_sample; ++j)
	{
		double l = 0;
		for(int m=0; m<n_comp; ++m)
		{
            l += Comps[m].Alpha * Comps[m].pdf(x[j], y[j]);
		}
		if (l>0 && !isnan(l))
		{
			llk += log(l);
		}
	}

	bic = -2.0 * llk + (5*n_comp-1.0)*log(n_sample);
	p_overlap = BayesError();
	//	std::cerr << "BIC: " << bic << ", P_OVERLAP: " << p_overlap << std::endl;
}


// 2-D EM for general without weights
void GaussianMixture2::EM2_select(std::vector<double>& x, std::vector<double> &y, std::vector<bool> &mask)
{
	// Let's not consider half-normal distribution -- for now

	int n_sample = (int) x.size();
	int n_iter = 15;

	// pseudo-counts
	double p_count = 0.2;

	double p_val[n_comp][2];

    DDPRINT("Starting clustering with %d components\n",  n_comp );

#ifdef DDEBUG
	print(stderr);
#endif

	if (n_comp == 1)
	{
		Comps[0].estimate_select(x, y, mask);
		Comps[0].Alpha = 1;
		
		bic = BIC_select(x, y, mask);
		p_overlap = 0;
#ifdef DDEBUG
        print(stderr);
#endif
	    DDMSG("BIC: " << bic << ", P_OVERLAP: N/A");
		return;
	}

	for(int i=0; i<n_comp; ++i)
	{
		p_val[i][0] = Comps[i].Mean[0];
		p_val[i][1] = Comps[i].Mean[1];

		if (Comps[i].Mean[0] < 0.01 && Comps[i].Mean[1] <0.01)
		{
			zeroidx = i;
		}
	}

	// pseudo-means
	//    double p_val[3] = {1.0, 0.5, 0};

	for(int i=0; i<n_iter; ++i)
	{
		std::vector<double> sum_x (n_comp,0);
		std::vector<double> sum_y (n_comp,0);

		std::vector<double> sum_e_xx (n_comp,0);
		std::vector<double> sum_e_xy (n_comp,0);
		std::vector<double> sum_e_yy (n_comp,0);

		std::vector<double> sum_pr (n_comp,0);

		double pr[n_comp][n_sample];
		bool b_include[n_sample];

		// E step
		for(int j=0; j<n_sample; ++j)
		{
			b_include[j] = false;
			if (!mask[j])
				continue;
			
			b_include[j] = true;

			double sum_p = 0;
			for(int m=0;m<n_comp;++m)
			{
				pr[m][j] = Comps[m].Alpha * Comps[m].pdf(x[j], y[j]);
 //               pr[m][j] *= 4.0;

				sum_p += pr[m][j];
			}
//			if (sum_p < 1e-10) b_include[j] = false;

			if (b_include[j]) // if the value is an outlier, exclude it from calculations
			{
				for(int m=0;m<n_comp;++m)
				{
					pr[m][j] /= sum_p;
					sum_x[m] += pr[m][j] * x[j];
					sum_y[m] += pr[m][j] * y[j];

					sum_pr[m] += pr[m][j];
				}
			}
		}
		double sumsum = 0;
		// Add pseudo-count values
		for(int m=0; m<n_comp; ++m)
		{
            if ( m!= zeroidx)
            {
                sum_x[m] += p_val[m][0] * p_count;
                sum_y[m] += p_val[m][1] * p_count;
                sum_pr[m] += p_count;
            }
			sumsum += sum_pr[m];
		}

		// M step

		// Update means
		for(int m=0;m<n_comp;++m)
		{
			if (m != zeroidx)
			{
				Comps[m].Mean[0] = sum_x[m]/sum_pr[m];
				Comps[m].Mean[1] = sum_y[m]/sum_pr[m];
			}
		}

		// Update covariance
		for(int j=0; j<n_sample; ++j)
		{
			if (mask[j] && b_include[j]) // if the value is an outlier, exclude it from calculations
			{
				for(int m=0;m<n_comp;++m)
				{
					double ex = x[j] - Comps[m].Mean[0];
					double ey = y[j] - Comps[m].Mean[1];

					sum_e_xx[m] += pr[m][j] * ex * ex;
					sum_e_xy[m] += pr[m][j] * ex * ey;
					sum_e_yy[m] += pr[m][j] * ey * ey;
				}
			}
		}

		for(int m=0;m<n_comp;++m)
		{
            double ex = p_val[m][0] - Comps[m].Mean[0];
            double ey = p_val[m][1] - Comps[m].Mean[1];

            sum_e_xx[m] += ex * ex;
            sum_e_xy[m] += ex * ey;
            sum_e_yy[m] += ey * ey;

            // Wishart prior : S = [0.01 0; 0 0.01]
            if (sum_pr[m]>1e-10)
            {
                Comps[m].Cov[0] = (sum_e_xx[m]  + 0.01 * p_count) / sum_pr[m];
                Comps[m].Cov[1] = Comps[m].Cov[2] = sum_e_xy[m] / sum_pr[m];
                Comps[m].Cov[3] = (sum_e_yy[m] + 0.01 * p_count) / sum_pr[m];
            }
            Comps[m].update();

            Comps[m].Alpha = sum_pr[m] / (sumsum);
		}
#ifdef DDEBUG
		print(stderr);
#endif
	}
	double llk = 0;
	int n_mask = 0;
	for(int j=0; j<n_sample; ++j)
	{
		if (!mask[j])
			continue;
		double l = 0;
		for(int m=0; m<n_comp; ++m)
		{
            l += Comps[m].Alpha * Comps[m].pdf(x[j], y[j]);
		}
		if (l>0 && !isnan(l))
		{
			llk += log(l);
		}
		n_mask ++;
	}

	if (n_mask>0)
	{
		bic = -2.0 * llk + (5*n_comp-1.0)*log((double)n_mask);
	}
	p_overlap = BayesError();
	std::cerr << "BIC: " << bic << ", P_OVERLAP: " << p_overlap << std::endl;
}

bool GaussianMixture2::ordered()
{
    if (Comps[0].Mean[0] < 0.8 || Comps[0].Mean[1] > 1.2 || Comps[0].Mean[1] < 0.8 || Comps[0].Mean[1] > 1.2 )
		return false;

    if ((Comps[1].Mean[0] < 0.4 && Comps[1].Mean[1] < 0.4) || (Comps[1].Mean[0] > 0.6 && Comps[1].Mean[1] > 0.6))
		return false;

	for (int i=0; i<n_comp-1; ++i)
	{
		if (Comps[i].Mean[0] + Comps[i].Mean[1] - Comps[i+1].Mean[0] - Comps[i+1].Mean[1] < 0.6)
        {
			return false;
        }
	}


    if (Comps[1].Alpha > 0.5 && (n_comp < 3 || Comps[2].Alpha < 0.01))
    {
        return false;
    }
	return true;
}

bool GaussianMixture2::r_ordered()
{
    
    if (Comps[0].Mean[0] < 0.8 || Comps[0].Mean[1] > 1.2 || Comps[0].Mean[1] < 0.8 || Comps[0].Mean[1] > 1.2 )
        return false;

    if ((Comps[1].Mean[0] < 1.4 && Comps[1].Mean[1] < 1.4))
		return false;

    /*
	if (Comps[0].Mean[0] + Comps[0].Mean[1] < 1.4 || Comps[0].Mean[0] + Comps[0].Mean[1] > 2.6 )
		return false;

	if (Comps[1].Mean[0] + Comps[1].Mean[1] < 2.6 )
		return false;
     */
    for (int i=0; i<n_comp-1; ++i)
    {
        if ( Comps[i+1].Mean[0] + Comps[i+1].Mean[1] - Comps[i].Mean[0] - Comps[i].Mean[1] < 0.6)
            return false;
    }
	return true;
}


int GaussianMixture2::assign_copynumber(double x, double y)
{
	double p[n_comp];
	double max_P = -1;
	double max_R = -1;
	int ret = -1;

	for(int i=0;i<n_comp; ++i)
	{
		// Major allele might domiate Alpha when sample size is large
		//        p[i] = C[i].Alpha * C[i].pdf(x);
		p[i] = Comps[i].pdf(x, y);

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
    
	// ret = round(Comps[ret].Mean[0] + Comps[ret].Mean[1]); // TODO: what if only one of the dimensions cluster correctly? (0, 0.5, 1) + (1, 1, 1) = (1 1.5 2) 

	// int up = ceil(x+y);
	// int down = floor(x+y);

	// if (ret != up && ret != down)
	// {
	// 	return -1;
	// }

	if (max_R > 0.1)
	{
		return -1;
	}

	return ret;
}

int GaussianMixture2::assign_dpcnt_copynumber(double x, double y)
{
	double p[n_comp];
	double max_P = -1;
	double max_R = -1;
	int ret = -1;

	for(int i=0;i<n_comp; ++i)
	{
		p[i] = Comps[i].pdf(x, y);

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
    
	
	if (ret >= 0 && ret < (int)Comps.size())
	{
		ret = 2 - ret;
		if (round(x*2.0) != ret)
		{
			if (x>1.0 && ret == 0 && y<0.05)
			{
				ret = 2;
			}
			else
			{
				ret = -1;
			}
		}
	}

	// ret = round(Comps[ret].Mean[0] + Comps[ret].Mean[1]); // TODO: what if only one of the dimensions cluster correctly? (0, 0.5, 1) + (1, 1, 1) = (1 1.5 2) 

	// int up = ceil(x+y);
	// int down = floor(x+y);

	// if (ret != up && ret != down)
	// {
	// 	return -1;
	// }

	if (max_R > 0.1)
	{
		return -1;
	}

	return ret;
}


double GaussianMixture2::BIC(std::vector<double>& x, std::vector<double>& y)
{
	int n_sample = (int)x.size();

	double ret = 0;
	double llk = 0;

	for(int j=0; j<n_sample; ++j)
	{
		double l = 0;
		for(int m=0; m<n_comp; ++m)
		{
            l += Comps[m].Alpha * Comps[m].pdf(x[j], y[j]);
		}
		if (l>0)
		{
			llk += log(l);
		}
	}

	ret = -2.0 * llk +  (5 * n_comp - 1.0) *log(n_sample);

	return ret;
}


double GaussianMixture2::BIC_select(std::vector<double>& x, std::vector<double>& y, std::vector<bool> &mask)
{
	int n_sample = 0;

	double ret = 0;
	double llk = 0;

	for(int j=0; j<(int)x.size(); ++j)
	{
		if (mask[j])
		{ 
			double l = 0;
			for(int m=0; m<n_comp; ++m)
			{
				l += Comps[m].Alpha * Comps[m].pdf(x[j], y[j]);
			}
			if (l>0)
			{
				llk += log(l);
			}
			n_sample ++;
		}
	}

	if (n_sample > 0)
	{
		ret = -2.0 * llk +  (5 * n_comp - 1.0) *log(n_sample);
	}

	return ret;
}

double GaussianMixture2::BayesError()
{
	double min_d = DBL_MAX;

	if (n_comp <2 )
	{
		return DBL_MAX;
	}

	// Get minimum distance between all pairs
	for(int i=0; i<n_comp-1; ++i)
	{
		for(int j=i+1; j<n_comp; ++j)
		{
			double C[4], P[4];
			for(int k=0;k<4;++k)
			{
				C[k] = (Comps[i].Cov[k] + Comps[j].Cov[k])/2.0;
			}
			double D = det(C);
			/*
			   while (D<1e-20)
			   {
			   C[0] += 1e-10;
			   C[1] += 1e-10;
			   D = det(C);
			   }
			 */
			P[0] = C[3]/D;
			P[3] = C[0]/D;
			P[1] = P[2] = -1.0*C[1]/D;
			double m1 = Comps[i].Mean[0] - Comps[j].Mean[0];
			double m2 = Comps[i].Mean[1] - Comps[j].Mean[1];

			double d = 0.125*(m1 * P[0] + m2 * P[2])*m1 + (m1 * P[1] + m2 * P[3])*m2 + 0.5*log( D / sqrt(Comps[i].Det*Comps[j].Det));

    		if (d<min_d)
			{
				min_d = d;
			}
		}
	}
	return exp(-1.0*min_d);
}

void GaussianMixture2::print(FILE *fp)
{
	for(int i=0; i<n_comp; ++i)
	{
        fprintf(fp, "Comp %d (%f,%f; %f,%f,%f), %f\n", i, Comps[i].Mean[0], Comps[i].Mean[1], Comps[i].Cov[0], Comps[i].Cov[1], Comps[i].Cov[3], Comps[i].Alpha);
    }
}

std::string GaussianMixture2::print_str()
{
	std::stringstream str;
	str << "N=" << n_comp << "|(";
	str << std::fixed << std::setprecision(4) << Comps[0].Mean[0] << "," << Comps[0].Mean[1] << ":" << Comps[0].Cov[0] << ":" << Comps[0].Cov[1]  << ":" << Comps[0].Cov[2] << ":" << Comps[0].Alpha;
    for(int i=1; i<n_comp; ++i)
	{			
		str << std::fixed << std::setprecision(4) << Comps[i].Mean[0] << "," << Comps[i].Mean[1] << ":" << Comps[i].Cov[0] << ":" << Comps[i].Cov[1]  << ":" << Comps[i].Cov[2] << ":" << Comps[i].Alpha;
	}
	str << std::fixed << std::setprecision(2) << "|BIC=" << bic;
	str << std::fixed << std::setprecision(4) << "|OVERLAP=" << p_overlap; 
	return str.str();
}
