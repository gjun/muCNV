//
//  base_file.cpp
//  muCNV
//
//  Created by Goo Jun on 11/25/18.
//  Copyright © 2018 Goo Jun. All rights reserved.
//

#include "base_file.h"

void BaseFile::open(std::string &fname, std::ios_base::openmode mode)
{
    filename = fname;
    fs.open(filename.c_str(), mode);
}

void BaseFile::close()
{
    if (fs.is_open())
    {
        fs.close();
    }
    else
    {
        std::cerr << "Error, " << filename << " is not open" << std::endl;
    }
}

int BaseFile::write_int32(int32_t number)
{
    fs.write(reinterpret_cast<char*>(&number), sizeof(int32_t));
    return(sizeof(int32_t));
}

int BaseFile::write_uint32(uint32_t number)
{
    fs.write(reinterpret_cast<char*>(&number), sizeof(uint32_t));
    return(sizeof(uint32_t));
}

int BaseFile::write_uint64(uint64_t number)
{
    fs.write(reinterpret_cast<char*>(&number), sizeof(uint64_t));
    return(sizeof(uint64_t));
}

int BaseFile::write_sample_id(std::string &sample_id)
{
    char pad[256] = {0};
    
    // Write sample ID (each with 256 bytes)
    if (sample_id.length() > 255)
    {
        std::cerr << "Error, sample ID " << sample_id << " is too long." << std::endl;
        exit(1);
    }
    
    fs.write(sample_id.c_str(), sample_id.length());
    fs.write(pad, 256-sample_id.length());
    
    return 256;
}

int BaseFile::read_int32(int32_t & number)
{
    fs.read(reinterpret_cast<char*>(&number), sizeof(int32_t));
    return(sizeof(int32_t));
}

int BaseFile::read_uint32(uint32_t & number)
{
    fs.read(reinterpret_cast<char*>(&number), sizeof(uint32_t));
    return(sizeof(uint32_t));
}

int BaseFile::read_uint64(uint64_t & number)
{
    fs.read(reinterpret_cast<char*>(&number), sizeof(uint64_t));
    return(sizeof(uint64_t));
}

int BaseFile::read_uint64_multi(uint64_t* array, int N)
{
    fs.read(reinterpret_cast<char*>(array), sizeof(uint64_t)*N);
    return(sizeof(uint64_t)*N);
}


int BaseFile::read_sample_id(char* buf)
{
    // read sample ID (each with 256 bytes)
    fs.read(buf, 256);
    return 256;
}
