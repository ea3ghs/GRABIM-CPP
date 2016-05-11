GRABIM_Result MatchingNetwork::RunGRABIM()
{
    GRABIM_Result Res;
    double lambda4 = c0/(4.*mean(f_matching));
    if (!topology.compare("-1"))//The user did not entered any specific network, so it seems
    {//reasonable to try some typical wideband matching network
        rowvec Vopt, Vaux;
        string candidate;
        double gridtest, opttopo;

        cout << "Grid search result for some common networks..." << endl;
        topology = "313202";//PiCPiL
        x_ini << 1e-12 << 1e-12 << 1e-12 << 1e-9 << 1e-9 << 1e-9;
        Vaux = GridSearch();
        gridtest = CandidateEval(Vaux);
        candidate = topology;
        opttopo = gridtest;
        Vopt = Vaux;
        cout << "PiCPiL: S11_min = " <<  gridtest << " dB" << endl;


        topology = "202313";//PiLPiC
        x_ini << 1e-9 << 1e-9 << 1e-9 << 1e-12 << 1e-12 << 1e-12;
        Vaux = GridSearch();
        gridtest = CandidateEval(Vaux);
        cout << "PiLPiC: S11_min = " <<  gridtest << " dB" << endl;
        if (gridtest < opttopo - 0.5)//It worths to change the topology if the new one improves the result significantly
        {
            candidate = topology;
            opttopo = gridtest;
            Vopt = Vaux;
        }

/*
        topology = "444";//3 cascaded lambda/4 sections
        double meanZS = mean(real(ZS_matching));
        double meanZL = mean(real(ZL_matching));

        if (meanZS < meanZL)
        {
            x_ini << 1.1*meanZS << lambda4<< .5*(meanZS+meanZL) << lambda4 << .9*meanZL << lambda4;
        }
        else
        {
            x_ini << 0.9*meanZS << lambda4<< .5*(meanZS+meanZL) << lambda4 << 1.1*meanZL << lambda4;
        }
        Vaux = GridSearch();
        gridtest = CandidateEval(Vaux);
        cout << "3-sections lambda/4: S11_min = " <<  gridtest << " dB" << endl;
        if (gridtest < opttopo - 0.5)
        {
            candidate = topology;
            opttopo = gridtest;
            Vopt = Vaux;
        }

        topology = "546";//Open circuit stub + Transmission line + Short circuited stub
        double meanZ = .5*(meanZS+meanZL);
        x_ini << meanZ << lambda4 << meanZ << lambda4<< meanZ << lambda4;
        Vaux = GridSearch();
        gridtest = CandidateEval(Vaux);
        cout << "OC stub + TL + SC stub: S11_min = " <<  gridtest << " dB" << endl;
        if (gridtest < opttopo - 0.5)
        {
            candidate = topology;
            opttopo = gridtest;
            Vopt = Vaux;
        }

        topology = "645";//Short circuited stub + Transmission line + Open circuit stub
        x_ini << meanZ << lambda4 << meanZ << lambda4<< meanZ << lambda4;
        Vaux = GridSearch();
        gridtest = CandidateEval(Vaux);
        cout << "SC stub + TL + OC stub: S11_min = " <<  gridtest << " dB" << endl;
        if (gridtest < opttopo - 0.5)
        {
            candidate = topology;
            opttopo = gridtest;
            Vopt = Vaux;
        }

*/

        topology = "020313";
        x_ini << 1e-9 << 1e-9 << 1e-9 << 1e-12 << 1e-12 << 1e-12;
        Vaux = GridSearch();
        gridtest = CandidateEval(Vaux);
        cout << "TeeLPiC: S11_min = " <<  gridtest << " dB" << endl;
        if (gridtest < opttopo - 0.5)
        {
            candidate = topology;
            opttopo = gridtest;
            Vopt = Vaux;
        }


        topology = "313020";
        x_ini  << 1e-12 << 1e-12 << 1e-12 << 1e-9 << 1e-9 << 1e-9;
        Vaux = GridSearch();
        gridtest = CandidateEval(Vaux);
        cout << "PiCTeeL: S11_min = " <<  gridtest << " dB" << endl;
        if (gridtest < opttopo - 0.5)
        {
            candidate = topology;
            opttopo = gridtest;
            Vopt = Vaux;
        }



        topology = "03030";
        x_ini << 1e-9 << 1e-12 << 1e-9 << 1e-12 << 1e-9;
        Vaux = GridSearch();
        gridtest = CandidateEval(Vaux);
        cout << "Lowpass LC ladder: S11_min = " <<  gridtest << " dB" << endl;
        if (gridtest < opttopo - 0.5)
        {
            candidate = topology;
            opttopo = gridtest;
            Vopt = Vaux;
        }

/*
        topology = "43434";
        x_ini << 75 << lambda4 << 1e-12 << 75 << lambda4  << 1e-12 << 75 << lambda4 ;
        Vaux = GridSearch();
        gridtest = CandidateEval(Vaux);
        cout << "Lowpass TL-C ladder: S11_min = " <<  gridtest << " dB" << endl;
        if (gridtest < opttopo - 0.5)
        {
            candidate = topology;
            opttopo = gridtest;
            Vopt = Vaux;
        }

*/

        topology = "01230123";//BPP3
        x_ini << 1e-9 << 1e-12 << 1e-9 << 1e-12<< 1e-9 << 1e-12<< 1e-9 << 1e-12;
        Vaux = GridSearch();
        gridtest = CandidateEval(Vaux);
        cout << "BPP3: S11_min = " <<  gridtest << " dB" << endl;
        if (gridtest < opttopo - 0.5)
        {
            candidate = topology;
            opttopo = gridtest;
            Vopt = Vaux;
        }

        topology = "23012301";//BPS3
        x_ini << 1e-9 << 1e-12 << 1e-9 << 1e-12<< 1e-9 << 1e-12<< 1e-9 << 1e-12;
        Vaux = GridSearch();
        gridtest = CandidateEval(Vaux);
        cout << "BPS3: S11_min = " <<  gridtest << " dB" << endl;
        if (gridtest < opttopo - .5)
        {
            candidate = topology;
            opttopo = gridtest;
            Vopt = Vaux;
        }


        topology = "131020131";
        x_ini << 1e-12 << 1e-12 << 1e-12 << 1e-9 << 1e-9 << 1e-9 << 1e-12 << 1e-12 << 1e-12;
        Vaux = GridSearch();
        gridtest = CandidateEval(Vaux);
        cout << "TeeCTeeLTeeC: S11_min = " <<  gridtest << " dB" << endl;
        if (gridtest < opttopo - .5)
        {
            candidate = topology;
            opttopo = gridtest;
            Vopt = Vaux;
        }


        topology = candidate;
        Res.x_grid_search = Vopt;

    }
    else//It builds an initial vector using default values 1nH for inductances, 1 pF for capacitors, 100 Ohm & lambda/4 for transmission lines
    {
        queue <double> XINI;
        for (unsigned int i = 0; i< topology.size();i++)
        {
            if ((!topology.substr(i,1).compare("0"))||(!topology.substr(i,1).compare("2"))) XINI.push(1e-9);
            if ((!topology.substr(i,1).compare("1"))||(!topology.substr(i,1).compare("3"))) XINI.push(1e-12);
            if ((!topology.substr(i,1).compare("4"))||(!topology.substr(i,1).compare("5"))||(!topology.substr(i,1).compare("6")))XINI.push(100),XINI.push(lambda4);
        }
        x_ini = ones(1, XINI.size());
        for (unsigned int i = 0; i < x_ini.size();i++)
        {
            x_ini.at(i) = XINI.front();
            XINI.pop();
        }
        Res.x_grid_search = GridSearch();
    }

    Res.grid_val = CandidateEval(Res.x_grid_search);
    Res.x_nlopt = LocalOptimiser(Res.x_grid_search);
    Res.nlopt_val = CandidateEval(Res.x_nlopt);
    Res.ZS = ZS;
    Res.ZL = ZL;
    Res.f_analysis = f_analysis;

    //Initialize S param vectors
    Res.S11_gridsearch=cx_vec(f_analysis,f_analysis);
    Res.S21_gridsearch=cx_vec(f_analysis,f_analysis);
    Res.S12_gridsearch=cx_vec(f_analysis,f_analysis);
    Res.S22_gridsearch=cx_vec(f_analysis,f_analysis);

    Res.S11_nlopt=cx_vec(f_analysis,f_analysis);
    Res.S21_nlopt=cx_vec(f_analysis,f_analysis);
    Res.S12_nlopt=cx_vec(f_analysis,f_analysis);
    Res.S22_nlopt=cx_vec(f_analysis,f_analysis);


    //Generate S parameter results
    cx_mat S_gridsearch, S_nlopt;
    S_gridsearch << 1<<1<<endr << 1<<1<< endr;
    S_nlopt << 1<<1<<endr << 1<<1<< endr;

    for (unsigned int i = 0; i < f_analysis.n_rows; i++)
    {
        // Grid search
        S_gridsearch = getSparams(Res.x_grid_search, ZS.at(i), ZL.at(i), f_analysis.at(i));
        Res.S11_gridsearch.at(i) = S_gridsearch(0,0);
        Res.S21_gridsearch.at(i) = S_gridsearch(1,0);
        Res.S12_gridsearch.at(i) = S_gridsearch(0,1);
        Res.S22_gridsearch.at(i) = S_gridsearch(1,1);

        // NLopt
        S_nlopt = getSparams(Res.x_nlopt, ZS.at(i), ZL.at(i), f_analysis.at(i));
        Res.S11_nlopt.at(i) = S_nlopt(0,0);
        Res.S21_nlopt.at(i) = S_nlopt(1,0);
        Res.S12_nlopt.at(i) = S_nlopt(0,1);
        Res.S22_nlopt.at(i) = S_nlopt(1,1);
    }



    return Res;
}
