#include <iostream>
#include "matchingnetwork.h"
#include "io.h"

using namespace std;

int main()
{
    MatchingNetwork MatchingObject;

    ///// MATCHING SETTINGS ////////
    double f1 = 1e9;
    double f2 = 2e9;
    int N = 20;//Number of samples to match in the matching band

    std::string sourcefile = "./s2p_files/f531p84g0.s2p";
    std::string loadfile = "./s2p_files/TQP7M9101.s2p";

    rowvec x_ini, x_opt;


    MatchingObject.SetMatchingBand(f1, f2, N);

    // SET SOURCE AND LOAD IMPEDANCES
    MatchingObject.SetSourceImpedance(sourcefile);
    MatchingObject.SetLoadImpedance(loadfile);
    vec freq = linspace(5e8, 10e9, 200);

    /*vec ZSr = 100*ones(freq.n_rows,1);
    vec ZSi = 0*ones(freq.n_rows,1);
    vec ZLr = 50*ones(freq.n_rows,1);
    vec ZLi = 0*ones(freq.n_rows,1);

    cx_vec ZS(ZSr, ZSi);
    cx_vec ZL(ZLr, ZLi);

    MatchingObject.SetSourceImpedance(ZS, freq);
    MatchingObject.SetLoadImpedance(ZL, freq);*/

    MatchingObject.SetTopology("131202");

    x_ini  << 1e-12 << 1e-12 << 1e-12 <<1e-9 << 1e-9 << 1e-9;
    MatchingObject.SetInitialPivot(x_ini);
    MatchingObject.SetMaxIterGridSearch(1000);
    MatchingObject.SetThreshold(-30);
    MatchingObject.SetNLoptAlg(nlopt::algorithm::LN_NELDERMEAD);
    MatchingObject.SetObjectiveFunction(ObjectiveFunction::NINF_POWERTRANS);
    GRABIM_Result R = MatchingObject.RunGRABIM();
    IO io;
    io.exportGNUplot(R, "GRABIM.dat");
    std::cout << "GRID SEARCH: S11_max = "<< R.grid_val << "dB <= " << R.x_grid_search << std::endl;
    std::cout << "NLOPT: S11_max = "<< R.nlopt_val << "dB <= " <<  R.x_nlopt << std::endl;

}



