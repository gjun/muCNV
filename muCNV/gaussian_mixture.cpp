//
//  mix_gaussian.cpp
//  muCNV
//
//  Created by Goo Jun on 12/3/18.
//  Copyright © 2018 Goo Jun. All rights reserved.
//

#include "gaussian_mixture.h"
#include <math.h>

// 2-D EM with weights
// EM with weights
void GaussianMixture::EM(std::vector<double>& x, std::vector<double> &w)
{
    unsigned n_sample = (unsigned) x.size();
    unsigned n_comp = (unsigned) Comps.size();
    unsigned n_iter = 10;
    
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



// 2-D EM for general without weights
void GaussianMixture2::EM2(std::vector<double>& x, std::vector<double> &y)
{
    // Let's not consider half-normal distribution -- for now
    
    unsigned n_sample = (unsigned) x.size();
    unsigned n_iter = 20;
    
    // pseudo-counts
    unsigned p_count= 5;
    double p_val[n_comp][2];
    
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
}



