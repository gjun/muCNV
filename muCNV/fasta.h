//
//  fasta.h
//  muCNV
//
//  Created by Goo Jun on 12/23/18.
//  Copyright © 2018 Goo Jun. All rights reserved.
//

#ifndef fasta_h
#define fasta_h

#include <string>
#include <iostream>
#include <fstream>
#include <map>

struct fastaChr
{
    size_t size;
    size_t length;
    size_t startpos;
    unsigned basesPerLine;
    unsigned bytesPerLine;
};

class fasta
{
private:
    std::ifstream fastaFile;
    std::map<std::string, fastaChr> faidx;
    std::vector<std::string> chrs;
    unsigned curr_column;
    unsigned columnwidth;
    unsigned columnbytes;
public:
    int load(const char*, const char*);;
    int seek(std::string, size_t);
    size_t chrlen(std::string) ;
    int read(int64_t, std::string &);
    void printIndex();
    ~fasta();
};

fasta::~fasta()
{
    if (fastaFile.is_open())
    {
        fastaFile.close();
    }
}

size_t fasta::chrlen(std::string chr)
{
    return faidx[chr].length;
}

int fasta::load(const char *fastaFileName, const char *faidxFileName)
{
    std::ifstream faidxFile;
    
    // Load FAIDX
    
    faidxFile.open(faidxFileName, std::ifstream::in);
    if (!faidxFile.is_open())
    {
        std::cerr << "Error opening faidx file: " << faidxFileName << std::endl;
        exit(1);
    }
    
    while(faidxFile.good())
    {
        std::string chr;
        fastaChr T;
        
        faidxFile >> chr;
        faidxFile >> T.length;
        faidxFile >> T.startpos;
        faidxFile >> T.basesPerLine;
        faidxFile >> T.bytesPerLine;
        faidx[chr] = T;
        chrs.push_back(chr);
    }
    faidxFile.close();
    
    curr_column = 0;
    columnwidth = 0;
    columnbytes = 0;
    
    // Load FASTA
    fastaFile.open(fastaFileName, std::ifstream::in);
    if (!fastaFile.is_open())
    {
        std::cerr << "Error opening fasta file: " << fastaFileName <<std::endl;
        exit(1);
    }
    
    return 0;
}


int fasta::seek(std::string chr, size_t pos)
{
    if (pos>faidx[chr].length)
    {
        pos = faidx[chr].length;
    }
    //    std::cout << "chr " << chr << " pos " << pos << std::endl;
    
    pos--;
    
    columnwidth = faidx[chr].basesPerLine;
    columnbytes = faidx[chr].bytesPerLine;
    size_t n_line = (size_t)(pos/columnwidth);
    //    std::cout << "line number " << n_line << std::endl;
    curr_column = pos%columnwidth;
    std::streampos fastaPtr = (std::streampos)(faidx[chr].startpos + n_line * faidx[chr].bytesPerLine + curr_column);
    
    //    std::cout << "file position  " << fastaPtr << std::endl;
    
    fastaFile.seekg(fastaPtr, std::ios::beg);
    
    return 0;
}

int fasta::read(int64_t len, std::string &contig)
{
    contig = "";
    char *buf = new char [256];
    
    //    printf("Curr column %d, column width %d, len %d\n", curr_column, columnwidth, len);
    // To do : chromosome length
    while(len>0)
    {
        unsigned readlen = columnwidth-curr_column;
        //        printf("readlen %d \n", readlen);
        if (readlen > len)
        {
            readlen = (unsigned int) len;
            //            printf("readlen shortened to %d \n", readlen);
            fastaFile.read(reinterpret_cast <char *> (buf), readlen);
            buf[readlen] = 0;
            contig += buf;
            //            printf("buf: %s\n", buf);
            
            curr_column += readlen;
        }
        else
        {
            fastaFile.read(reinterpret_cast <char *> (buf), readlen);
            buf[readlen] = 0;
            //            printf("buf: %s\n", buf);
            contig += buf;
            
            fastaFile.read(reinterpret_cast <char *> (buf), columnbytes-columnwidth);
            buf[columnbytes-columnwidth] = 0;
            //            printf("null buf: %s\n", buf);
            
            curr_column = 0;
        }
        
        len -= readlen;
        //        printf("len shortened to %d \n", len);
    }
    
    return 0;
}

void fasta::printIndex()
{
    // Just for sanity check
    for(int i=1;i<23;++i)
    {
        std::string chr = "chr" + std::to_string(i);;
        std::cout << chr << "\t" << faidx[chr].length << "\t" << faidx[chr].startpos << "\t" << faidx[chr].basesPerLine << "\t" << faidx[chr].bytesPerLine << "\n";
    }
}

#endif /* fasta_h */
