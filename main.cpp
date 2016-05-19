#include "GRABIM.h"
#include "io.h"
#include<cctype>

using namespace std;


cx_double getComplexImpedanceFromText(char *arg)
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

    if (argc == 7)//Check seed setup
    {
        if (!strcmp(argv[5], "--set-seed"))
        {
            arma_rng::set_seed(atof(argv[6]));
        }
        else
        {
            arma_rng::set_seed_random();  // set the seed to a random value
        }
    }


    string SourceFile = argv[1];
    string LoadFile = argv[2];

    bool ZSisConstant = isdigit(SourceFile.at(0));
    bool ZLisConstant = isdigit(LoadFile.at(0));

    //Before starting the matching engine, we must ensure that the impedance data is already loaded
    if (SourceFile.empty() && (!ZSisConstant))
    {
        cerr << "Please specify a Touchstone file for the source" << endl;
        return 0;
    }

    if (LoadFile.empty()&& (!ZLisConstant))
    {
        cerr << "Please specify a Touchstone file for the load" << endl;
        return 0;
    }

    //Check filename extension
    int formatSource=-1, formatLoad=-1;;
    if (SourceFile.find(".s1p")) formatSource = 0;
    if (LoadFile.find(".s1p")) formatLoad = 0;


    if ((formatSource != 0) && (!ZSisConstant))
    {
        cerr << "Unknown source impedace" << endl;
        return 0;
    }

    if ((formatLoad != 0)  && (!ZLisConstant))
    {
        cerr << "Unknown load impedance" << endl;
        return 0;
    }

    //Impedance data paths were already specified, let's proceed to bring the S-par data into memory
    IO inout_operations;

    if (!ZSisConstant)//Read source impedance from file
    {
        inout_operations.loadS1Pdata(SourceFile, SOURCE);//s1p
    }
    else//Set constant source impedance
    {
       cx_double zs_temp;
       char * text = argv[1];
       zs_temp = getComplexImpedanceFromText(text);
       if (zs_temp.real() == -1)//Check if the input value is correct
       {
           cerr << "The input given for the source impedance is not valid" << endl;
           return 0;
       }
       inout_operations.set_constant_ZS_vs_freq(zs_temp);
    }

    if (!ZLisConstant)
    {
       inout_operations.loadS1Pdata(LoadFile, LOAD);//s1p
    }
    else
    {
        cx_double zl_temp;
        char * text = argv[2];
        zl_temp = getComplexImpedanceFromText(text);

        if (zl_temp.real() == -1)//Check if the input value is correct
        {
            cerr << "The input given for the load impedance is not valid" << endl;
            return 0;
        }
        inout_operations.set_constant_ZL_vs_freq(zl_temp);

    }

    //Check frequency specifications
    double fmatching_min = atof(argv[3]);
    double fmatching_max = atof(argv[4]);

    if ((fmatching_min == -1) || (fmatching_max == -1))
    {
        cerr << "Incorrect frequency settings" << endl;
        return 0;
    }
    else//Everything correct... lets set frequency
    {
        inout_operations.set_matching_band(fmatching_min, fmatching_max);

        //Check if the specified frequencies lie with the s1p/s2p data
        inout_operations.ResampleImpedances();//Force data update
        if (fmatching_min < inout_operations.getFrequency().min())//The lower freq is not present at s1p/s2p
        {
            cerr <<"One of the impedance data files does not contain the specified lower frequency" << endl;
            return 0;
        }
        if (fmatching_max > inout_operations.getFrequency().max())//The maximum freq is not present at s1p/s2p
        {
            cerr <<"One of the impedance data files does not contain the specified upper frequency" << endl;
            return 0;
        }
    }


    string TopoScript_path = "predefined_topologies";
    string QucsSchPath = "GRABIM_result.sch";

    GRABIM MatchingObject;
    // Impedance and frequency settings
    MatchingObject.SetSourceImpedance(inout_operations.getSourceImpedance());
    MatchingObject.SetLoadImpedance(inout_operations.getLoadImpedance());
    MatchingObject.SetFrequency(inout_operations.getFrequency());
    MatchingObject.setTopoScript(TopoScript_path);

    MatchingObject.SetTopology("-1");

    inout_operations.setLocalOptimiser(nlopt::LN_NELDERMEAD);
    inout_operations.set_qucs_sch_path(QucsSchPath);


    GRABIM_Result R = MatchingObject.RunGRABIM();//Runs GRABIM. Well, this is not exactly the algorithm
    // detailed at [1] but from my point of view is functional and easier to code...
    //Notes:
    // 1) The candidate vector is not in log units. I do not see a good reason for doing so. Maybe I am missing something important
    // 2) Frequency is not in rad/s.
    // 3) The objective function is the magnitude of S11 expressed in dB. log(x) functions usually have strong
    // gradients so it seem to suggest that this is good for derivative free opt algorithms
    // 4) This code takes advantage from NLopt derivative-free local optimisers. NLopt is easy to link and it works
    // fine. Despite the fact that the Nelder-Mead algorithm does not guarantee convergence (among other problems), it leads to achieve a good local
    // (probably, global) optimum. This is caused by the fact that the matching network should be as simple as possible => few elements => xk \in R^N, where
    // N is typically < 6. Even N=6 is a big number, please consider that matching networks are tight to physical constraints in practice, so, the larger the
    // network, the harder the 'tuning'.

    inout_operations.PrintNetwork_StandardOutput(R);

    (ZSisConstant) ? R.source_path = "" : R.source_path = SourceFile;
    (ZLisConstant) ? R.load_path = "": R.load_path = LoadFile;

    R.QucsVersion = "0.0.19";


    string GNUplot_path = "GRABIM.dat";

    cout << "Finished: GRABIM has successfully finished." << endl;
    cout << "Please execute: 'gnuplot plotscript' to take a look to the results." << endl;
    cout << "A new Qucs schematic has been generated at "<< inout_operations.get_qucs_sch_path() << endl;

    inout_operations.exportGNUplot(R, GNUplot_path);
    inout_operations.ExportQucsSchematic(R);
}



