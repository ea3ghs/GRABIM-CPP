#include <iostream>
#include "matchingnetwork.h"
#include "io.h"

using namespace std;

//This function parses the source and load impedance entered from command line.
cx_double check_string(char *arg)
{
    string arg_str = arg;
    int index = arg_str.find_first_of("j");
    int sign = 1;
    double zreal, zimag;
    if (index != -1)//Then it is a single impedance
    {
        zreal = atof(arg_str.substr(0, index-1).c_str());
        if (!arg_str.substr(index-1, 1).compare("-")) sign = -1;
        zimag = sign*atof(arg_str.substr(index+1).c_str());
        return cx_double(zreal, zimag);
    }
    else
    {
        zreal = atof(arg);
        if (zreal > 0)return cx_double(zreal, 0);
        else
            return cx_double(-1, -1);
    }
}

int main(int argc, char *argv[])
{
    MatchingNetwork MatchingObject;
    cx_double ZS, ZL;

    //DEFAULT VALUES
    MatchingObject.SetTopology("-1");//By default, no topology is asigned
    int plot = 1;
    string output_data = "GRABIM.dat";
    vec freq = linspace(1e6, 100e6, 200);//Set the frequency vector por plotting the final results
    int N = 20;//Number of samples to match in the matching band

    // SETTINGS
    if (argc == 9)
    {

printf("\n\n\n\n!!\n!! GRABIM-CPP rev20160512\n!!\n!! https://github.com/andresmmera/GRABIM-CPP\n");
printf("!! [SRC:%s]--[MatchingNetwork:%s]--[LOAD:%s]\n",argv[1],argv[5],argv[2]);
printf("!! BW: %s..%s\n",argv[3],argv[4]);
printf("!! algorithm: %s\n",argv[6]);

        MatchingObject.SetMatchingBand(atof(argv[3]), atof(argv[4]), N);
        ZS = check_string(argv[1]);
        ZL = check_string(argv[2]);

        if(ZS.real() == -1) MatchingObject.SetSourceImpedance(argv[1]);
        if(ZL.real() == -1) MatchingObject.SetLoadImpedance(argv[2]);
        MatchingObject.SetTopology(argv[5]);
        string NL_ALG = argv[6];
        //NLopt algorithm
        if (NL_ALG.compare("NLOPT_LN_PRAXIS")) MatchingObject.SetNLoptAlg(nlopt::LN_PRAXIS);
        else if (NL_ALG.compare("NLOPT_LN_NELDERMEAD")) MatchingObject.SetNLoptAlg(nlopt::LN_NELDERMEAD);
        else if (NL_ALG.compare("NLOPT_LN_SBPLX")) MatchingObject.SetNLoptAlg(nlopt::LN_SBPLX);
        else if (NL_ALG.compare("NLOPT_LN_COBYLA")) MatchingObject.SetNLoptAlg(nlopt::LN_COBYLA);
        else if (NL_ALG.compare("NLOPT_LN_BOBYQA")) MatchingObject.SetNLoptAlg(nlopt::LN_BOBYQA);
        else if (NL_ALG.compare("NLOPT_LN_AUGLAG")) MatchingObject.SetNLoptAlg(nlopt::LN_AUGLAG);
        //GLOBAL (genetic...) OPT ALGORITHMS
        else if (NL_ALG.compare("NLOPT_GN_ESCH")) MatchingObject.SetNLoptAlg(nlopt::GN_ESCH);
        else if (NL_ALG.compare("NLOPT_GN_ISRES")) MatchingObject.SetNLoptAlg(nlopt::GN_ISRES);
        else if (NL_ALG.compare("NLOPT_GD_STOGO")) MatchingObject.SetNLoptAlg(nlopt::GD_STOGO);

        plot = atoi(argv[7]);//Display results?
        output_data = argv[8];//Output data for GNUplot
    }
    else//Something went wrong...
    {
        cout << "Invalid arguments:" << endl;
        cout << "./GRABIM-CPP <source-impedance> <load-impedance> <topology> <NL-algorithm> <enable-display> <output-data-file>";
        return -1;
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

    //Grid search settings
    MatchingObject.SetMaxIterGridSearch(1000);//Maximum number of iterations for the grid search algorithm
    MatchingObject.SetThreshold(-30);//It specifies the mininum S11 [dB] required, typically,
                                     //S11 < 10 dB is considered as valid for common applications
    MatchingObject.SetObjectiveFunction(ObjectiveFunction::NINF_S11dB);//Sets the kind of objective functions. Everything seems to
                                                                       //suggest that NINF_S11dB gives the best results

    GRABIM_Result R = MatchingObject.RunGRABIM();
    IO io;io.exportGNUplot(R, output_data, plot);
    //std::cout << "GRID SEARCH: S11_max = "<< R.grid_val << "dB <= " << R.x_grid_search;
    //std::cout << "NLOPT: S11_max = "<< R.nlopt_val << "dB <= " <<  R.x_nlopt;
    
}



