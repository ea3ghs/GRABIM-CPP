#include <iostream>
#include "matchingnetwork.h"
#include "io.h"

using namespace std;

cx_double check_string(char *arg)
{
    string arg_str = arg;
    int index = arg_str.find_first_of("j");
    if (index != -1)//Then it is a single impedance
    {
        double zreal = atof(arg_str.substr(0, index-1).c_str());
        double zimag = atof(arg_str.substr(index+1).c_str());
        return cx_double(zreal, zimag);
    }
    else return cx_double(-1, -1);
}

int main(int argc, char *argv[])
{
    MatchingNetwork MatchingObject;

    std::string sourcefile;
    std::string loadfile;
    MatchingObject.SetTopology("-1");
    cx_double ZS, ZL;

    vec freq = linspace(5e8, 10e9, 200);//Set the frequency vector por plotting the final results

    int N = 20;//Number of samples to match in the matching band

    // SOURCE, LOAD IMPEDANCES AND MATCHING SETTINGS
    if (argc == 5)
    {
       MatchingObject.SetMatchingBand(atof(argv[3]), atof(argv[4]), N);
       ZS = check_string(argv[1]);
       ZL = check_string(argv[2]);
       if(ZS.real() == -1) MatchingObject.SetSourceImpedance(argv[1]);
       if(ZL.real() == -1) MatchingObject.SetLoadImpedance(argv[2]);
       MatchingObject.SetNLoptAlg(nlopt::LN_NELDERMEAD);
    }
    if (argc == 6)
    {
       MatchingObject.SetMatchingBand(atof(argv[3]), atof(argv[4]), N);
       ZS = check_string(argv[1]);
       ZL = check_string(argv[2]);
       if(ZS.real() == -1) MatchingObject.SetSourceImpedance(argv[1]);
       if(ZL.real() == -1) MatchingObject.SetLoadImpedance(argv[2]);
       string fifth_arg = argv[5];
       if(fifth_arg.substr(0,2).compare("NL"))
       {
           MatchingObject.SetTopology(argv[5]);
       }
       else
       {
           //The 5th argument is the NLopt algorithm
           if (fifth_arg.compare("NLOPT_LN_PRAXIS")) MatchingObject.SetNLoptAlg(nlopt::LN_PRAXIS);
           else if (fifth_arg.compare("NLOPT_LN_NELDERMEAD")) MatchingObject.SetNLoptAlg(nlopt::LN_NELDERMEAD);
           else if (fifth_arg.compare("NLOPT_LN_SBPLX")) MatchingObject.SetNLoptAlg(nlopt::LN_SBPLX);
           else if (fifth_arg.compare("NLOPT_LN_NELDERMEAD")) MatchingObject.SetNLoptAlg(nlopt::LN_COBYLA);
           else if (fifth_arg.compare("NLOPT_LN_NELDERMEAD")) MatchingObject.SetNLoptAlg(nlopt::LN_BOBYQA);
           else if (fifth_arg.compare("NLOPT_LN_NELDERMEAD")) MatchingObject.SetNLoptAlg(nlopt::LN_AUGLAG);
           //GLOBAL OPT ALGORITHMS
           else if (fifth_arg.compare("NLOPT_LN_NELDERMEAD")) MatchingObject.SetNLoptAlg(nlopt::GN_ESCH);
           else if (fifth_arg.compare("NLOPT_LN_NELDERMEAD")) MatchingObject.SetNLoptAlg(nlopt::GN_ISRES);
           else if (fifth_arg.compare("NLOPT_LN_NELDERMEAD")) MatchingObject.SetNLoptAlg(nlopt::GD_STOGO);
       }
    }


    ///// MATCHING SETTINGS ////////
    rowvec x_ini, x_opt;

    if (ZS.real() != -1)
    {
        vec ZSr = ZS.real()*ones(freq.n_rows,1);
        vec ZSi = ZS.imag()*ones(freq.n_rows,1);
        cx_vec zs(ZSr, ZSi);
        MatchingObject.SetSourceImpedance(zs,freq);
    }

    if (ZL.real() != -1)
    {
        vec ZLr = ZL.real()*ones(freq.n_rows,1);
        vec ZLi = ZL.imag()*ones(freq.n_rows,1);
        cx_vec zl(ZLr, ZLi);
        MatchingObject.SetLoadImpedance(zl,freq);
    }


   // x_ini  << 100 << 0.05 << 100 <<0.05 << 100 <<0.05;
    MatchingObject.SetInitialPivot(x_ini);
    MatchingObject.SetMaxIterGridSearch(1000);
    MatchingObject.SetThreshold(-30);
    MatchingObject.SetObjectiveFunction(ObjectiveFunction::NINF_S11dB);
    GRABIM_Result R = MatchingObject.RunGRABIM();
    IO io;
    io.exportGNUplot(R, "GRABIM.dat");
    std::cout << "GRID SEARCH: S11_max = "<< R.grid_val << "dB <= " << R.x_grid_search << std::endl;
    std::cout << "NLOPT: S11_max = "<< R.nlopt_val << "dB <= " <<  R.x_nlopt << std::endl;

}



