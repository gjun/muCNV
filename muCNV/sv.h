//
//  sv.h
//  muCNV
//
//  Created by Goo Jun on 11/13/18.
//  Copyright © 2018 Goo Jun. All rights reserved.
//

#ifndef sv_h
#define sv_h

#include <string>
#include "gaussian.h"

enum svType {DEL=0, DUP=1, INV=2, CNV=3, INS=4, BND=5};
std::string svTypeName(svType);
svType svTypeNum(int);

class sv
{
public:
    svType svtype;
    //    string source;
    //    string chr;
    int chrnum;
    int pos;
    int end;
    int len;
    int supp;
    //    pair<int,int> ci_pos;
    //    pair<int,int> ci_end;
    void get_len()
    {
        len = end - pos + 1;
    };
    void print(void);
    bool operator < (const sv&) const;
    bool operator == (const sv&) const;
    
    uint64_t dp_sum;
    int n_dp;
    uint16_t dp;
    
    sv();
};

class breakpoint
{
public:
    uint8_t chrnum;
    int bptype; // 0 : pos-gap, 1: pos, 2: pos+gap, 3:end-gap, 4: end, 5: end+gap
    
    int pos;
    int idx; // SV_id
    bool operator < (const breakpoint&) const;
    bool operator == (const breakpoint&) const;
    bool operator <= (const breakpoint&) const;
    
    breakpoint();
};


class svdata
{
public:
    int n;
    std::vector<double> dp;
    std::vector<double> isz;
    
    std::vector<double> cnv_pos;
    std::vector<double> cnv_neg;
    std::vector<double> inv_pos;
    std::vector<double> inv_neg;
    
    std::vector<int> n_isz;
    std::vector<int> n_cnv_pos;
    std::vector<int> n_cnv_neg;
    std::vector<int> n_inv_pos;
    std::vector<int> n_inv_neg;
    
    std::vector<double> norm_dp;
    std::vector<double> norm_readcount;
    std::vector<double> norm_cnv_pos;
    std::vector<double> norm_cnv_neg;
    std::vector<double> norm_inv_pos;
    std::vector<double> norm_inv_neg;
    
    void set_size(int num)
    {
        n = num;
        dp.resize(n);
        isz.resize(n);
        
        cnv_pos.resize(n);
        cnv_neg.resize(n);
        inv_pos.resize(n);
        inv_neg.resize(n);
        
        n_isz.resize(n);
        n_cnv_pos.resize(n);
        n_cnv_neg.resize(n);
        n_inv_pos.resize(n);
        n_inv_neg.resize(n);
        
        norm_dp.resize(n);
        norm_cnv_pos.resize(n);
        norm_cnv_neg.resize(n);
        norm_inv_pos.resize(n);
        norm_inv_neg.resize(n);
        norm_readcount.resize(n);
    };
    
    void normalize(sv &interval, std::vector<double> &avg_depth, std::vector<double> &avg_isize)
    {
        for(int i=0;i<n;++i)
        {
            norm_dp[i] = dp[i] / avg_depth[i];
            if (interval.svtype == DEL)
            {
                norm_cnv_pos[i] = (cnv_pos[i] - avg_isize[i] ) / interval.len;
                norm_cnv_neg[i] = (cnv_neg[i] - avg_isize[i] ) / interval.len;
                norm_inv_pos[i] = (inv_pos[i] - avg_isize[i] ) / interval.len;
                norm_inv_neg[i] = (inv_neg[i] - avg_isize[i] ) / interval.len;
            }
            else if (interval.svtype == DUP || interval.svtype == CNV)
            {
                norm_cnv_pos[i] = (cnv_pos[i] + avg_isize[i]) / interval.len;
                norm_cnv_neg[i] = (cnv_neg[i] + avg_isize[i]) / interval.len;
                norm_inv_pos[i] = (inv_pos[i] + avg_isize[i]) / interval.len;
                norm_inv_neg[i] = (inv_neg[i] + avg_isize[i]) / interval.len;
            }
            else if (interval.svtype == INV )
            {
                norm_cnv_pos[i] = (cnv_pos[i] ) / avg_isize[i];
                norm_cnv_neg[i] = (cnv_neg[i] ) / avg_isize[i];
                norm_inv_pos[i] = (inv_pos[i] ) / interval.len;
                norm_inv_neg[i] = (inv_neg[i] ) / interval.len;
            }
            // READLEN fixed to 150 : later!!
            norm_readcount[i] = (double)n_isz[i] * 150.0 / (interval.len + 2.0*(avg_isize[i]-75)) / avg_depth[i];
        }
    };
    
    void print(sv &interval)
    {
    }
    
};


double RO(sv &, sv &);

void pick_sv_from_merged(sv &, std::vector<sv> &);
int find_overlap_sv(sv& , std::vector<sv>&);
int find_start(std::vector<sv> &, int );

bool in_centrome(sv &);
bool in_centrome(int, int);

#endif /* sv_h */
