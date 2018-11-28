//
//  base_file.hpp
//  muCNV
//
//  Created by Goo Jun on 11/25/18.
//  Copyright © 2018 Goo Jun. All rights reserved.
//

#ifndef base_file_hpp
#define base_file_hpp

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cstdint>

class BaseFile
{
public:
    void open(std::string &, std::ios_base::openmode);
    void close();
    
    int write_int32(int32_t);
    int write_uint32(uint32_t);
    int write_uint64(uint64_t);
    int write_sample_id(std::string &);

    int read_int32(int32_t &);
    int read_uint32(uint32_t &);
    int read_uint64(uint64_t &);
    int read_uint64_multi(uint64_t *, int);
    int read_sample_id(char *);
    size_t tellg() { return fs.tellg(); };
    bool good() { return fs.good(); };
    std::istream& seekg(size_t p) { return fs.seekg(p); };
    
protected:
    std::fstream fs;
    std::string filename;
};

#endif /* base_file_hpp */
