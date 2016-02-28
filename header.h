#ifndef HEADER_H
#define HEADER_H

#include <iostream>
#include <math.h>
#include <nlopt.h>
#include <armadillo>


const double c0 = 299792458;//Speed of light (m/s)
const double Z0 = 50;//Ohm. Reference impedance

using namespace arma;

typedef struct DeviceData {
    cx_mat S;
    cx_mat Z;
    vec freq;

} DeviceData;


DeviceData ReadS2PFile(std::string);


#endif // HEADER_H
