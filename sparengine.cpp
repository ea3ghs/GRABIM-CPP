#include "sparengine.h"

SparEngine::SparEngine()
{
}


// Returns the S matrix at a given frequency
cx_mat SparEngine::getSparams(rowvec x, cx_double zs, cx_double zl, double f, std::string topology)
{
    cx_mat ABCD = getABCDmatrix(x, f, topology);
    cx_mat S;
    S << -1 << -1 << endr << -1 << -1 << endr;
    //Convert ABCD to S parameters
    S(0,0) = (ABCD(0,0)*zl+ABCD(0,1)-ABCD(1,0)*conj(zs)*zl-ABCD(1,1)*conj(zs))/(ABCD(0,0)*zl+ABCD(0,1)+ABCD(1,0)*zs*zl+ABCD(1,1)*zs);
    S(0,1) = (2.*(ABCD(0,0)*ABCD(1,1)-ABCD(0,1)*ABCD(1,0))*sqrt(real(zs)*real(zl)))/(ABCD(0,0)*zl+ABCD(0,1)+ABCD(1,0)*zs*zl+ABCD(1,1)*zs);
    S(1,0) = (2.*sqrt(real(zs)*real(zl)))/(ABCD(0,0)*zl+ABCD(0,1)+ABCD(1,0)*zs*zl+ABCD(1,1)*zs);
    S(1,1) = (-ABCD(0,0)*conj(zl)+ABCD(0,1)-ABCD(1,0)*conj(zl)*zs+ABCD(1,1)*zs)/(ABCD(0,0)*zl+ABCD(0,1)+ABCD(1,0)*zs*zl+ABCD(1,1)*zs);
    return S;
}


// Returns the ABCD matrix at a given frequency
cx_mat SparEngine::getABCDmatrix(rowvec x, double f, std::string topology)
{
    int element;
    double w = 2*datum::pi*f;
    double beta = w/c0;
    cx_double gamma = cx_double(0, beta);
    cx_mat ABCD, ABCD_t;
    ABCD << 1 << 0 << endr << 0 << 1 << endr;

    unsigned int i, k;

    for (i = 0, k=0; i < topology.length(); i++, k++)
    {
        element = atoi(topology.substr(i,1).c_str());
        switch(element)
        {
        case 0: ABCD_t << 1. << cx_double(0,w*x.at(k)) << endr << 0 << 1. << endr;
            break;
        case 1: ABCD_t << 1. << cx_double(0,-1/(w*x.at(k))) << endr << 0 << 1. << endr;
            break;
        case 2: ABCD_t << 1. << 0 << endr << cx_double(0,-1./(w*x.at(k))) << 1. << endr;
            break;
        case 3: ABCD_t << 1. << 0 << endr << cx_double(0, w*x.at(k)) << 1. << endr;
            break;
        case 4: ABCD_t << cosh(gamma*x.at(k+1)) << x.at(k)*sinh(gamma*x.at(k+1)) << endr << sinh(gamma*x.at(k+1))/x(k) << cosh(gamma*x.at(k+1)) << endr;
            k++;//It involves two parameters, so we need to skip the next index
            break;
        case 5: ABCD_t << 1. << 0 << endr << (tanh(gamma*x.at(k+1)))/x.at(k) << 1. << endr;
            k++;
            break;
        case 6: ABCD_t << 1. << 0 << endr << 1./(x.at(k)*tanh(gamma*x.at(k+1))) << 1. << endr;
            k++;
            break;
        default: return -ABCD.eye();
        }

        ABCD = ABCD*ABCD_t;
    }
    return ABCD;
}
