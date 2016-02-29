#ifndef MATCHINGNETWORK_H
#define MATCHINGNETWORK_H
#include <nlopt.hpp>
#include <armadillo>
#include <string.h>
#include <iostream>
#include <queue>

using namespace arma;

const double c0 = 299792458;//Speed of light (m/s)

typedef struct DeviceData {
    cx_mat S;
    cx_mat Z;
    vec freq;

} DeviceData;

typedef struct GRABIM_Result {
    rowvec x_grid_search;
    double grid_val;
    rowvec x_nlopt;
    double nlopt_val;

} GRABIM_Result;


class MatchingNetwork
{
public:
    MatchingNetwork();
    int SetSourceImpedance(std::string);
    int SetSourceImpedance(cx_vec, vec);
    int SetLoadImpedance(std::string);
    int SetLoadImpedance(cx_vec, vec);
    int SetInitialPivot(rowvec);
    rowvec GetInitialPivot();
    int SetMatchingBand(double, double, int);
    int SetTopology(std::string);
    GRABIM_Result RunGRABIM();
    int SetMaxIterGridSearch(int);
    int SetThreshold(double);
    double GetThreshold();
    int SetVerbose(bool);
    double CandidateEval(rowvec);
    int SetNLoptAlg(nlopt::algorithm);
    nlopt::algorithm GetNLoptAlg();

private:
    DeviceData LoadS2PData(std::string);
    rowvec GridSearch();
    rowvec LocalOptimiser(rowvec);
    int ResampleImpedances();
    cx_mat getSparams(rowvec, cx_double, cx_double, double);
    mat GeneratingMatrix(int);
    rowvec x_ini;
    vec f_matching;
    cx_vec ZS, ZL;
    cx_vec ZS_matching, ZL_matching;
    vec fS, fL;
    vec f_analysis;
    bool verbose;
    std::string topology;
    int Grid_MaxIter;
    double MatchingThreshold;
    nlopt::algorithm NLoptAlgo;
};

#endif // MATCHINGNETWORK_H
