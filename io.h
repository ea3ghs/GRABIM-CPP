#ifndef IO_H
#define IO_H
#include "matchingnetwork.h"
#include <iostream>
#include <fstream>

using namespace std;

class IO
{
public:
    IO();
    int exportGNUplot(GRABIM_Result, string, int);
};

#endif // IO_H
