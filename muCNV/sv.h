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
#endif /* sv_h */
