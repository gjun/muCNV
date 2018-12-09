//
//  mix_gaussian.cpp
//  muCNV
//
//  Created by Goo Jun on 12/3/18.
//  Copyright © 2018 Goo Jun. All rights reserved.
//

#include "gaussian_mixture.h"
#include <math.h>

// 2-D EM without weights
// EM with weights
GaussianMixture::GaussianMixture(std::vector<double> &m, std::vector<double> &s)
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
        Comps[i].Mean = m[i];
        Comps[i].Stdev = s[i];
        Comps[i].Alpha = 1.0/n_comp;
		if (m[i] < 0.01)
			zeroidx = i;
    }
}

GaussianMixture& GaussianMixture::operator = (const GaussianMixture& gmix)
{
    n_comp = (int) gmix.n_comp;
    Comps.resize(gmix.n_comp);
    
    for(int i=0; i<n_comp; ++i)
    {
        Comps[i].Mean = gmix.Comps[i].Mean;
        Comps[i].Stdev = gmix.Comps[i].Stdev;
        Comps[i].Alpha = gmix.Comps[i].Alpha;
    }
    bic = gmix.bic;
    p_overlap = gmix.p_overlap;
	zeroidx = gmix.zeroidx;
    
    return *this;
}

void GaussianMixture::print()
{
	for(int i=0; i<n_comp; ++i)
	{
		std::cerr << "Comp " << i <<  " (" << Comps[i].Mean << ", " << Comps[i].Stdev << "), " << Comps[i].Alpha << std::endl;
	}
}

void GaussianMixture::EM(std::vector<double>& x)
{
    int n_sample = (int) x.size();
    int n_iter = 10;
    
    int p_count= 1;
    double p_val[n_comp];
    
    if (n_comp == 1)
    {
        Comps[0].estimate(x);
        Comps[0].Alpha = 1;
        bic = BIC(x);
        p_overlap = 0;
        return;
    }
    for(int i=0; i<n_comp; ++i)
    {
        p_val[i] = Comps[i].Mean;
    }
    
    for(int i=0; i<n_iter; ++i)
    {
        std::vector<double> sum (n_comp,0);
        std::vector<double> sum_pr (n_comp,0);
        std::vector<double> sum_err (n_comp,0);
        
        double pr[n_comp][n_sample];
        bool b_include[n_sample];
        
        // E step
        for(int j=0; j<n_sample; ++j)
        {
            double sum_p = 0;
            for(int m=0;m<n_comp;++m)
            {
                pr[m][j] = Comps[m].Alpha * normpdf(x[j], Comps[m]);
                if (zeroidx == (int)m)
                    pr[m][j] *= 2.0; // Truncated normal.
                sum_p += pr[m][j];
            }
            
            if (sum_p < 1e-10)
                b_include[j] = false; // if the value is an outlier, exclude it from calculations
            else
                b_include[j] = true;
            
            if (b_include[j])
            {
                for(int m=0;m<n_comp;++m)
                {
                    pr[m][j] /= sum_p;
                    sum[m] += pr[m][j] * x[j];
                    sum_pr[m] += pr[m][j];
                }
            }

        }
        
        double sumsum = 0;
        // Add pseudo-count values
        for(int m=0; m<n_comp; ++m)
        {
            sum[m] += p_val[m] * p_count;
            sum_pr[m] += p_count;
            sumsum += sum_pr[m];
        }
        // M step
        for(int m=0;m<n_comp;++m)
        {
            if (m != zeroidx)   // Mean of '0' in truncated normal should always be 0
                Comps[m].Mean = sum[m]/sum_pr[m];
        }
        
        for(int j=0; j<n_sample; ++j)
        {
            for(int m=0; m<n_comp; ++m)
            {
                sum_err[m] += pr[m][j] * (x[j] - Comps[m].Mean)*(x[j]-Comps[m].Mean);
            }
        }
        
        // Pseudo-count for variance -- commented out
        /*
        for(int m=0; m<n_comp; ++m)
        {
            sum_err[m] += (p_val[m] - Comps[m].Mean)*(p_val[m]-Comps[m].Mean) * p_count;
        }
        */
        for(int m=0; m<n_comp; ++m)
        {
			Comps[m].Stdev = sqrt(sum_err[m] / sum_pr[m]) ;
			Comps[m].Alpha = sum_pr[m] / sumsum;
        }
        //            cerr << "\t(" << Comps[m].Mean << "," << Comps[m].Stdev  << " : " << Comps[m].Alpha << ")";
            //        cerr << endl;
//		print();
    }
    
	double llk = 0;
    for(int j=0; j<n_sample; ++j)
    {
        double l = 0;
        for(int m=0; m<n_comp; ++m)
        {
			if (m != zeroidx)
				l += Comps[m].Alpha * normpdf(x[j], Comps[m]);
			else
				l += 2.0 * Comps[m].Alpha * normpdf(x[j], Comps[m]);
        }
        if (l>0)
        {
            llk += log(l);
        }
    }
    
    bic = -2.0 * llk +  2*n_comp*log(n_sample);
    p_overlap = BayesError();
}


int GaussianMixture::assign_copynumber(double x)
{
    double p[n_comp];
    double max_P = -1;
    double max_R = -1;
    int ret = -1;
    
    for(int i=0;i<n_comp; ++i)
    {
        // Major allele might domiate Alpha when sample size is large
        //        p[i] = C[i].Alpha * C[i].pdf(x);
        p[i] = Comps[i].pdf(x);
        
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
    
    ret = round(Comps[ret].Mean * 2);
    
    int up = ceil(x*2.0);
    int down = floor(x*2.0);
    
    if (ret != up && ret != down)
    {
        return -1;
    }
    
    if (max_R >= 0.2)
    {
        return -1;
    }
    
    return ret;
}


bool GaussianMixture::ordered()
{
    for(int i=0; i<n_comp-1;++i)
    {
        if (Comps[i].Mean<=Comps[i+1].Mean)
            return false;
    }
    return true;
}

bool GaussianMixture::r_ordered()
{
    for(int i=0; i<(int)Comps.size()-1;++i)
    {
        if (Comps[i].Mean>=Comps[i+1].Mean)
            return false;
    }
    return true;
}


double GaussianMixture::BIC(std::vector<double>& x)
{
    int n_sample = (int)x.size();

    double ret = 0;
	double llk = 0;
    
    for(int j=0; j<n_sample; ++j)
    {
        double l = 0;
        for(int m=0; m<n_comp; ++m)
        {
			if (m != zeroidx)
				l += Comps[m].Alpha * normpdf(x[j], Comps[m]);
			else
				l += 2.0 * Comps[m].Alpha * normpdf(x[j], Comps[m]);
        }
        if (l>0)
        {
            llk += log(l);
        }
    }
    
    ret = -2.0 * llk +  2*n_comp*log(n_sample);
    
    return ret;
}


double GaussianMixture::BayesError()
{
    // Returns maximum Bhattacharyya coefficient (== exp(-D) ) between components
    double min_d = DBL_MAX;
    
    if (n_comp <2 )
    {
        return 1;
    }
    
    // Get minimum distance between all pairs
    for(int i=0; i<n_comp-1; ++i)
    {
        for(int j=i+1; j<n_comp; ++j)
        {
            double s = (Comps[i].Stdev*Comps[i].Stdev + Comps[j].Stdev*Comps[j].Stdev); // sigma_p^2 + sigma_q^2
            double d = (Comps[i].Mean-Comps[j].Mean)*(Comps[i].Mean-Comps[j].Mean)/(4.0*s) + 0.5*log( 0.5 * s / (Comps[i].Stdev*Comps[j].Stdev));
            if (d<min_d)
            {
                min_d = d;
            }
        }
    }
    return exp(-1.0*min_d);
}



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
		if (Comps[i].Mean[0] + Comps[i].Mean[1] < 0.01)
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
            Comps[i].Cov[j] = gmix.Comps[i].Cov[j];
        Comps[i].Alpha = gmix.Comps[i].Alpha;
    }
    bic = gmix.bic;
    p_overlap = gmix.p_overlap;
	zeroidx = gmix.zeroidx;
    
    return *this;
}


// 2-D EM for general without weights
void GaussianMixture2::KM2(std::vector<double>& x, std::vector<double> &y)
{
    // Let's not consider half-normal distribution -- for now
    
    int n_sample = (int) x.size();
    int n_iter = 20;
    
    // pseudo-counts
    int p_count= 5;
    
    if (n_comp == 1)
    {
        Comps[0].estimate(x, y);
        Comps[0].Alpha = 1;
        bic = BIC(x, y);
        p_overlap = 0;
        return;
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
                double dist = dx*dx + dy*dy;
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
            Comps[m].Mean[0] = sum_x[m] / cnt[m];
            Comps[m].Mean[1] = sum_y[m] / cnt[m];
        }
        print();

    }
    
    std::vector<double> sum_xx(n_comp, 0);
    std::vector<double> sum_yy(n_comp, 0);
    std::vector<double> sum_xy(n_comp, 0);
    std::vector<int> cnt(n_comp, 0);
    for(int j=0; j<n_sample; ++j)
    {
        double dx = Comps[member[j]].Mean[0] - x[j];
        double dy = Comps[member[j]].Mean[1] - y[j];
        cnt[member[j]] ++;
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
    print();

    double llk = 0;
    for(int j=0; j<n_sample; ++j)
    {
        double l = 0;
        for(int m=0; m<n_comp; ++m)
        {
            if (m != zeroidx)
                l += Comps[m].Alpha * Comps[m].pdf(x[j], y[j]);
            else
                l += Comps[m].Alpha * Comps[m].pdf(x[j], y[j]) * 4.0;
        }
        if (l>0 && !isnan(l))
        {
            llk += log(l);
            if (isnan(llk))
            {
                std::cerr << "isnan!" << std::endl;
            }
        }
    }
    
    bic = -2.0 * llk + 5*n_comp*log(n_sample);
    p_overlap = BayesError();
    std::cerr << "BIC: " << bic << ", P_OVERLAP: " << p_overlap << std::endl;
}

// 2-D EM for general without weights
void GaussianMixture2::EM2(std::vector<double>& x, std::vector<double> &y)
{
    // Let's not consider half-normal distribution -- for now
    
    int n_sample = (int) x.size();
    int n_iter = 20;
    
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
				if (m == zeroidx)
				{
					pr[m][j] *= 4.0;
				}
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
        // Add pseudo-count values
        for(int m=0; m<n_comp; ++m)
        {
            sum_x[m] += p_val[m][0] * p_count;
            sum_y[m] += p_val[m][1] * p_count;
            sum_pr[m] += p_count;
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
			if (m != zeroidx)
			{
				double ex = p_val[m][0] - Comps[m].Mean[0];
				double ey = p_val[m][1] - Comps[m].Mean[1];
				
				sum_e_xx[m] += ex * ex;
				sum_e_xy[m] += ex * ey;
				sum_e_yy[m] += ey * ey;
				
				// Wishart prior : S = [0.01 0; 0 0.01]
				Comps[m].Cov[0] = (sum_e_xx[m]  + 0.001 * p_count) / sum_pr[m];
				Comps[m].Cov[1] = Comps[m].Cov[2] = sum_e_xy[m] / sum_pr[m];
				Comps[m].Cov[3] = (sum_e_yy[m] + 0.001 * p_count) / sum_pr[m];
				Comps[m].update();
			}
            
            Comps[m].Alpha = sum_pr[m] / (n_sample + n_comp*p_count);
        }
		print();
    }
	double llk = 0;
    for(int j=0; j<n_sample; ++j)
    {
        double l = 0;
        for(int m=0; m<n_comp; ++m)
        {
			if (m != zeroidx)
				l += Comps[m].Alpha * Comps[m].pdf(x[j], y[j]);
			else
				l += Comps[m].Alpha * Comps[m].pdf(x[j], y[j]) * 4.0;
        }
        if (l>0 && !isnan(l))
        {
            llk += log(l);
			if (isnan(llk))
			{
				std::cerr << "isnan!" << std::endl;
			}
        }
    }
    
    bic = -2.0 * llk + 5*n_comp*log(n_sample);
    p_overlap = BayesError();
	std::cerr << "BIC: " << bic << ", P_OVERLAP: " << p_overlap << std::endl;
}


bool GaussianMixture2::ordered()
{
    for(int i=0; i<n_comp-1;++i)
    {
        if (Comps[i].Mean[0] + Comps[i].Mean[1] <= Comps[i+1].Mean[0] + Comps[i+1].Mean[1])
            return false;
    }
    return true;
}

bool GaussianMixture2::r_ordered()
{
    for(int i=0; i<(int)Comps.size()-1;++i)
    {
        if (Comps[i].Mean[0]+Comps[i].Mean[1] >= Comps[i+1].Mean[0]+Comps[i+1].Mean[1])
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
    
    ret = round(Comps[ret].Mean[0] + Comps[ret].Mean[1]); // TODO: what if only one of the dimensions cluster correctly? (0, 0.5, 1) + (1, 1, 1) = (1 1.5 2) 
    
    int up = ceil(x+y);
    int down = floor(x+y);
    
    if (ret != up && ret != down)
    {
		std::cerr << std::endl << "ret " << ret << " up " << up << " down " << down << std::endl;
		std::cerr << "x " << x << " y " << y << std::endl;
        return -1;
    }
    
    if (max_R > 0.2)
    {
		std::cerr << "max_R" << max_R <<  std::endl;
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
			if (m != zeroidx)
				l += Comps[m].Alpha * Comps[m].pdf(x[j], y[j]);
			else
				l += 4.0*Comps[m].Alpha * Comps[m].pdf(x[j], y[j]);
        }
        if (l>0)
        {
            llk += log(l);
        }
    }
    
    ret = -2.0 * llk +  5*n_comp*log(n_sample);
    
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

void GaussianMixture2::print()
{
	for(int i=0; i<n_comp; ++i)
	{
		std::cerr << "Comp " << i <<  " (" << Comps[i].Mean[0] << "," << Comps[i].Mean[1] << ";" << Comps[i].Cov[0] << "," << Comps[i].Cov[1] <<"." << Comps[i].Cov[3]  << "), " << Comps[i].Alpha << std::endl;
	}
}

