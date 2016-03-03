#include <iostream>
#include "matchingnetwork.h"
#include "io.h"

using namespace std;

int main(int argc, char *argv[])
{
    MatchingNetwork MatchingObject;

    if (argc == 3)
    {

    }

    ///// MATCHING SETTINGS ////////
    double f1 = 1e9;
    double f2 = 3e9;
    int N = 20;//Number of samples to match in the matching band

    std::string sourcefile = "./s2p_files/TQP7M9101.s2p";
    std::string loadfile = "./s2p_files/TQP7M9102.s2p";

    rowvec x_ini, x_opt;


    MatchingObject.SetMatchingBand(f1, f2, N);

    // SET SOURCE AND LOAD IMPEDANCES
    //MatchingObject.SetSourceImpedance(sourcefile);
    //MatchingObject.SetLoadImpedance(loadfile);
    vec freq = linspace(5e8, 6e9, 200);

    vec ZSr = 100*ones(freq.n_rows,1);
    vec ZSi = 0*ones(freq.n_rows,1);
    vec ZLr = 50*ones(freq.n_rows,1);
    vec ZLi = 0*ones(freq.n_rows,1);
    cx_vec ZS(ZSr, ZSi);
    cx_vec ZL(ZLr, ZLi);
    MatchingObject.SetSourceImpedance(ZS, freq);
    MatchingObject.SetLoadImpedance(ZL, freq);

    MatchingObject.SetTopology("-1");

   // x_ini  << 100 << 0.05 << 100 <<0.05 << 100 <<0.05;
    MatchingObject.SetInitialPivot(x_ini);
    MatchingObject.SetMaxIterGridSearch(1000);
    MatchingObject.SetThreshold(-30);
    MatchingObject.SetNLoptAlg(nlopt::algorithm::LN_PRAXIS);
    MatchingObject.SetObjectiveFunction(ObjectiveFunction::NINF_S11dB);
    GRABIM_Result R = MatchingObject.RunGRABIM();
    IO io;
    io.exportGNUplot(R, "GRABIM.dat");
    std::cout << "GRID SEARCH: S11_max = "<< R.grid_val << "dB <= " << R.x_grid_search << std::endl;
    std::cout << "NLOPT: S11_max = "<< R.nlopt_val << "dB <= " <<  R.x_nlopt << std::endl;

}



