#ifndef __COMMON_H__
#define __COMMON_H__

#include <stdio.h>

#include "sv.h"

void split(const char*, const char*, std::vector<std::string>&);
template <class T> void vprint(std::vector<T> &);
int median(std::vector<int> &);
double average(std::vector<double> &);

#endif // __COMMON_H__
