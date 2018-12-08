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
    llk = gmix.llk;
    p_overlap = gmix.p_overlap;
    
    return *this;
}

void GaussianMixture::EM(std::vector<double>& x)
{
    unsigned n_sample = (unsigned) x.size();
    unsigned n_iter = 12;
    
    unsigned p_count= 1;
    double p_val[n_comp];
    int zeroidx = -1;
    
    if (n_comp == 1)
    {
        Comps[0].estimate(x);
        Comps[0].Alpha = 1;
        bic = BIC(x);
        p_overlap = 0;
        return;
    }
    for(unsigned i=0; i<n_comp; ++i)
    {
        p_val[i] = Comps[i].Mean;
        if (Comps[i].Mean < 0.01)
        {
            zeroidx = i;
        }
    }
    
    for(unsigned i=0; i<n_iter; ++i)
    {
        std::vector<double> sum (n_comp,0);
        std::vector<double> sum_pr (n_comp,0);
        std::vector<double> sum_err (n_comp,0);
        
        double pr[n_comp][n_sample];
        bool b_include[n_sample];
        
        // E step
        for(unsigned j=0; j<n_sample; ++j)
        {
            double sum_p = 0;
            for(unsigned m=0;m<n_comp;++m)
            {
                pr[m][j] = Comps[m].Alpha * normpdf(x[j], Comps[m]);
                if (zeroidx == (int)m)
                    pr[m][j] *= 2.0; // Truncated normal.
                sum_p += pr[m][j];
            }
            
            if (sum_p < 1e-30)
                b_include[j] = false; // if the value is an outlier, exclude it from calculations
            else
                b_include[j] = true;
            
            if (b_include[j])
            {
                for(unsigned m=0;m<n_comp;++m)
                {
                    pr[m][j] /= sum_p;
                    sum[m] += pr[m][j] * x[j];
                    sum_pr[m] += pr[m][j];
                }
            }

        }
        
        double sumsum = 0;
        // Add pseudo-count values
        for(unsigned m=0; m<n_comp; ++m)
        {
            sum[m] += p_val[m] * p_count;
            sum_pr[m] += p_count;
            sumsum += sum_pr[m];
        }
        // M step
        for(unsigned m=0;m<n_comp;++m)
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
        for(unsigned m=0; m<n_comp; ++m)
        {
            sum_err[m] += (p_val[m] - Comps[m].Mean)*(p_val[m]-Comps[m].Mean) * p_count;
        }
        */
        for(unsigned m=0; m<n_comp; ++m)
        {
            Comps[m].Stdev = sqrt(sum_err[m] / sum_pr[m]) ;
            Comps[m].Alpha = sum_pr[m] / sumsum;
            if (Comps[m].Stdev < 1e-10)
                Comps[m].Stdev = 0.005;
        }
        //            cerr << "\t(" << Comps[m].Mean << "," << Comps[m].Stdev  << " : " << Comps[m].Alpha << ")";
            //        cerr << endl;
    }
    
    for(int j=0; j<n_sample; ++j)
    {
        double l = 0;
        for(int m=0; m<n_comp; ++m)
        {
            l += Comps[m].Alpha * normpdf(x[j], Comps[m]);
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
    
    if (max_R > 0.2)
    {
        return -1;
    }
    
    return ret;
}


bool GaussianMixture::ordered()
{
    for(unsigned i=0; i<n_comp-1;++i)
    {
        if (Comps[i].Mean<=Comps[i+1].Mean)
            return false;
    }
    return true;
}

bool GaussianMixture::r_ordered()
{
    for(unsigned i=0; i<Comps.size()-1;++i)
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
    
    for(int j=0; j<n_sample; ++j)
    {
        double l = 0;
        for(int m=0; m<n_comp; ++m)
        {
            l += Comps[m].Alpha * normpdf(x[j], Comps[m]);
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
    for(unsigned i=0; i<n_comp-1; ++i)
    {
        for(unsigned j=i+1; j<n_comp; ++j)
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


// 2-D EM with weights
// EM with weights
void GaussianMixture::wEM(std::vector<double>& x, std::vector<double> &w)
{
    unsigned n_sample = (unsigned) x.size();
    unsigned n_comp = (unsigned) Comps.size();
    unsigned n_iter = 12;
    
    unsigned p_count= 1;
    double p_val[n_comp];
    //    int zeroidx = -1;
    
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
        std::vector<double> sum (n_comp,0);
        std::vector<double> sum_pr (n_comp,0);
        std::vector<double> sum_err (n_comp,0);

        
        // E step
        for(unsigned j=0; j<n_sample; ++j)
        {
            double sum_p = 0;
            std::vector<double> pr(n_comp, 0);
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
            
            //            cerr << "\t(" << Comps[m].Mean << "," << Comps[m].Stdev  << " : " << Comps[m].Alpha << ")" ;
            /*
             if (zeroidx == (int)m)
             {
             Comps[m].Mean = 0;
             }
             */
        }
        //        cerr << endl;
    }
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
    llk = gmix.llk;
    p_overlap = gmix.p_overlap;
    
    return *this;
}


// 2-D EM for general without weights
void GaussianMixture2::EM2(std::vector<double>& x, std::vector<double> &y)
{
    // Let's not consider half-normal distribution -- for now
    
    unsigned n_sample = (unsigned) x.size();
    unsigned n_iter = 12;
    
    // pseudo-counts
    unsigned p_count= 5;
    double p_val[n_comp][2];
    
    if (n_comp == 1)
    {
        Comps[0].estimate(x, y);
        Comps[0].Alpha = 1;
        bic = BIC(x, y);
        p_overlap = 0;
        return;
    }
    
    for(unsigned i=0; i<n_comp; ++i)
    {
        p_val[i][0] = Comps[i].Mean[0];
        p_val[i][1] = Comps[i].Mean[1];
    }
    
    // pseudo-means
    //    double p_val[3] = {1.0, 0.5, 0};
    
    for(unsigned i=0; i<n_iter; ++i)
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
        for(unsigned j=0; j<n_sample; ++j)
        {
            b_include[j] = true;
            double sum_p = 0;
            for(unsigned m=0;m<n_comp;++m)
            {
                pr[m][j] = Comps[m].Alpha * Comps[m].pdf(x[j], y[j]);
                sum_p += pr[m][j];
            }
            if (sum_p < 1e-30) b_include[j] = false;
            
            if (b_include[j]) // if the value is an outlier, exclude it from calculations
            {
                for(unsigned m=0;m<n_comp;++m)
                {
                    pr[m][j] /= sum_p;
                    sum_x[m] += pr[m][j] * x[j];
                    sum_y[m] += pr[m][j] * y[j];
                    
                    sum_pr[m] += pr[m][j];
                }
            }
        }
        // Add pseudo-count values
        for(unsigned m=0; m<n_comp; ++m)
        {
            sum_x[m] += p_val[m][0] * p_count;
            sum_y[m] += p_val[m][1] * p_count;
            sum_pr[m] += p_count;
        }
        
        // M step
        
        // Update means
        for(unsigned m=0;m<n_comp;++m)
        {
            Comps[m].Mean[0] = sum_x[m]/sum_pr[m];
            Comps[m].Mean[1] = sum_y[m]/sum_pr[m];
        }
        
        // Update covariance
        for(unsigned j=0; j<n_sample; ++j)
        {
            if (b_include[j]) // if the value is an outlier, exclude it from calculations
            {
                for(unsigned m=0;m<n_comp;++m)
                {
                    double ex = x[j] - Comps[m].Mean[0];
                    double ey = y[j] - Comps[m].Mean[1];
                    
                    sum_e_xx[m] += pr[m][j] * ex * ex;
                    sum_e_xy[m] += pr[m][j] * ex * ey;
                    sum_e_yy[m] += pr[m][j] * ey * ey;
                }
            }
        }
        
        for(unsigned m=0;m<n_comp;++m)
        {
            double ex = p_val[m][0] - Comps[m].Mean[0];
            double ey = p_val[m][1] - Comps[m].Mean[1];
            
            sum_e_xx[m] += ex * ex;
            sum_e_xy[m] += ex * ey;
            sum_e_yy[m] += ey * ey;
            
            Comps[m].Cov[0] = sum_e_xx[m] / sum_pr[m];
            Comps[m].Cov[1] = Comps[m].Cov[2] = sum_e_xy[m] / sum_pr[m];
            Comps[m].Cov[3] = sum_e_yy[m] / sum_pr[m];
            
            Comps[m].Alpha = sum_pr[m] / (n_sample + n_comp*p_count);
            Comps[m].update();
        }
    }
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
    
    bic = -2.0 * llk +  5*n_comp*log(n_sample);
    p_overlap = BayesError();
}


bool GaussianMixture2::ordered()
{
    for(unsigned i=0; i<n_comp-1;++i)
    {
        if (Comps[i].Mean[0] + Comps[i].Mean[1] <= Comps[i+1].Mean[0] + Comps[i+1].Mean[1])
            return false;
    }
    return true;
}

bool GaussianMixture2::r_ordered()
{
    for(unsigned i=0; i<Comps.size()-1;++i)
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
    
    int up = ceil(x*2.0);
    int down = floor(x*2.0);
    
    if (ret != up && ret != down)
    {
        return -1;
    }
    
    if (max_R > 0.2)
    {
        return -1;
    }
    
    return ret;
}


double GaussianMixture2::BIC(std::vector<double>& x, std::vector<double>& y)
{
    int n_sample = (int)x.size();
    
    double ret = 0;
    
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
            if (D>1e-10)
            {
                P[0] = C[3]/D;
                P[3] = C[0]/D;
                P[1] = P[2] = -1.0*C[1]/D;
            }
            else
            {
                // for singular cases, just use Cov = [0.01 0; 0 0.01]
                P[0]=100;
                P[1]=0;
                P[2]=0;
                P[3]=100;
            }
            
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
