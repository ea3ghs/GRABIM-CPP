#include "matchingnetwork.h"

MatchingNetwork::MatchingNetwork()
{
    verbose  = false;
}

// It sets if the program should show detailed information about the iterations or not
int MatchingNetwork::SetVerbose(bool vb)
{
    verbose = vb;
    return 0;
}

// It loads a s2p file and calculates the output impedance
int MatchingNetwork::SetSourceImpedance(std::string sourcepath)
{
    DeviceData DATA = LoadS2PData(sourcepath);
    if (DATA.freq.at(0) == -1) return -1;//Prevents from corrupted file format
    fS = DATA.freq;
    ZS = DATA.Z.col(3);
    if (!ZL.empty())ResampleImpedances();
    return 0;
}

// It sets the source impedance vs frequency
int MatchingNetwork::SetSourceImpedance(cx_vec zs, vec freq)
{
    if (zs.n_rows != freq.n_rows) return -1;
    ZS = zs;
    ZS_matching = zs;
    fS = freq;
    if (!ZL.empty())ResampleImpedances();
    return 0;
}

// It loads a s2p file and calculates the input impedance
int MatchingNetwork::SetLoadImpedance(cx_vec zl, vec freq)
{
    if (zl.n_rows != freq.n_rows) return -1;
    ZL = zl;
    ZL_matching = ZL;
    fL = freq;
    if (!ZS.empty())ResampleImpedances();
    return 0;
}

// It sets the load impedance vs frequency
int MatchingNetwork::SetLoadImpedance(std::string loadpath)
{
    DeviceData DATA = LoadS2PData(loadpath);
    if (DATA.freq.at(0) == -1) return -1;//Prevents from corrupted file format
    fL = DATA.freq;
    ZL = DATA.Z.col(0);
    if (!ZS.empty())ResampleImpedances();
    return 0;
}

// Sets the first point of the grid search
int MatchingNetwork::SetInitialPivot(rowvec x)
{
    x_ini = x;
    return 0;
}

// Returns the first point of the grid search
rowvec MatchingNetwork::GetInitialPivot()
{
    return x_ini;
}

//Sets the portion of the spectrum where matching is required
int MatchingNetwork::SetMatchingBand(double f1, double f2, int N)
{
    f_matching = linspace(f1, f2, N);
    return 0;
}

// This function sets the ladder arrangement of the matching network
int MatchingNetwork::SetTopology(std::string s)
{
    topology = s;
    return 0;
}

string tolower(string str)
{
    char c;
    for (unsigned int i =0; i < str.length();i++)
    {
        c=str.at(i);
        str.at(i) = tolower(c);
    }
    return str;
}


string RemoveBlankSpaces(string line)
{
    //Remove consecutive repeated blank spaces and space at the beginning
    //Sometimes, the fields may be separated by \t...
    for (unsigned int i = 0; i< line.length(); i++)
    {
        if (i == 0)//Remove first space
        {
            if((!line.substr(0,1).compare(" "))||(!line.substr(0,1).compare("\t")))
            {
                line.erase(0, 1);
                i--;
            }
            continue;
        }
        if (((!line.substr(i-1, 1).compare(" "))||(!line.substr(i-1,1).compare("\t")))&&((!line.substr(i, 1).compare(" "))||(!line.substr(i,1).compare("\t"))))
        {
            line.erase(i, 1);
            i--;
        }
    }
    return line;
}

//Loads a s2p data file and parses its contents. These data will be used to set source/load impedances
DeviceData MatchingNetwork::LoadS2PData(std::string filepath)
{
    DeviceData DATA;
    std::ifstream s2pfile(filepath.c_str());//Tries to open the data file.
    if(!s2pfile.is_open())//The data file cannot be opened => error
    {
        DATA.freq.at(0) = -1;
        return DATA;
    }

    std::string line;
    double freq_scale = 1;
    double Zref = 50;

    std::getline(s2pfile, line);
    while(line.compare(0, 1, "#"))//Looking for # field
    {
        std::getline(s2pfile, line);
    }

    line = tolower(line);
    //Freq scale
    if (line.find("ghz") != -1)
    {
        freq_scale = 1e9;
    }
    else
    {
        if (line.find("mhz") != -1)
        {
            freq_scale = 1e6;
        }
        else
        {
            if ((line.find("khz") != -1))
            {
                freq_scale = 1e3;
            }
        }
    }


    //Get the impedance at which the S params were measured

    int Rindex = line.find_last_of("r");
    Rindex = line.find_first_not_of(" ", Rindex);
    Zref = atof(line.substr(Rindex+1).c_str());
    int is_indB = line.find("db");


    while( getline(s2pfile, line) )
    {//Looking for the start of the raw data

        line = RemoveBlankSpaces(line);

        if ((!line.compare(0,1, "!"))|| (line.length() == 1)) continue;
        else break;


    }

    //DATA beginning.
    //At this point, the number of frequency samples is not known, so it's better to
    //push data into queues and then arrange it into armadillo structures
    std::queue <double> frequency, S11M, S11A, S21M, S21A, S12M, S12A, S22M, S22A;
    unsigned int qsize=0;

    do
    {
        line = RemoveBlankSpaces(line);

        if (line.empty()|| (line.length()==1))break;
        if (line.at(0) == '!') break;//Comment

        //Frequency
        int index = line.find_first_of(" ");
        if (index == -1)index = line.find_first_of("\t");
        frequency.push(atof(line.substr(0,index).c_str()));
        line.erase(0, index+1);


        index = line.find_first_of(" ");
        if (index == -1)index = line.find_first_of("\t");
        S11M.push(atof(line.substr(0,index).c_str()));
        line.erase(0, index+1);

        index = line.find_first_of(" ");
        if (index == -1)index = line.find_first_of("\t");
        S11A.push(atof(line.substr(0,index).c_str()));
        line.erase(0, index+1);

        index = line.find_first_of(" ");
        if (index == -1)index = line.find_first_of("\t");
        S21M.push(atof(line.substr(0,index).c_str()));
        line.erase(0, index+1);

        index = line.find_first_of(" ");
        if (index == -1)index = line.find_first_of("\t");
        S21A.push(atof(line.substr(0,index).c_str()));
        line.erase(0, index+1);

        index = line.find_first_of(" ");
        if (index == -1)index = line.find_first_of("\t");
        S12M.push(atof(line.substr(0,index).c_str()));
        line.erase(0, index+1);

        index = line.find_first_of(" ");
        if (index == -1)index = line.find_first_of("\t");
        S12A.push(atof(line.substr(0,index).c_str()));
        line.erase(0, index+1);

        index = line.find_first_of(" ");
        if (index == -1)index = line.find_first_of("\t");
        S22M.push(atof(line.substr(0,index).c_str()));
        line.erase(0, index+1);

        index = line.find_first_of(" ");
        if (index == -1)index = line.find_first_of("\t");
        if (index == -1)index = line.find_first_of("\r");//End of the line
        S22A.push(atof(line.substr(0,index).c_str()));
        qsize++;
    }while (std::getline(s2pfile, line));

    vec freq(qsize);
    cx_mat S(qsize, 4);
    cx_mat Z(qsize, 4);
    double phi, S11m, S11a, S21m, S21a, S12m, S12a, S22m, S22a;
    for (unsigned int i = 0; i < qsize; i++)
    {
        freq(i) = freq_scale*frequency.front();
        frequency.pop();

        S11m = S11M.front();
        S11a = S11A.front();
        S11M.pop();
        S11A.pop();
        if (is_indB != -1) S11m = pow(10, .05*S11m);
        phi = (datum::pi/180)*S11a;
        S(i, 0) = cx_double(S11m,0)*cx_double(cos(phi), sin(phi));

        S21m = S21M.front();
        S21a = S21A.front();
        S21M.pop();
        S21A.pop();
        if (is_indB != -1) S21m = pow(10, .05*S21m);
        phi = (datum::pi/180)*S21a;
        S(i, 1) = cx_double(S21m,0)*cx_double(cos(phi), sin(phi));


        S12m = S12M.front();
        S12a = S12A.front();
        S12M.pop();
        S12A.pop();
        if (is_indB != -1) S12m = pow(10, .05*S12m);
        phi = (datum::pi/180)*S12a;
        S(i, 2) = cx_double(S12m,0)*cx_double(cos(phi), sin(phi));



        S22m = S22M.front();
        S22a = S22A.front();
        S22M.pop();
        S22A.pop();
        if (is_indB != -1) S22m = pow(10, .05*S22m);
        phi = (datum::pi/180)*S22a;
        S(i, 3) = cx_double(S22m,0)*cx_double(cos(phi), sin(phi));

        cx_double K = cx_double(Zref,0)/((1.-S(i,0))*(1.-S(i,3))-S(i,1)*S(i,2));
        Z(i, 0) = K*((1.+S(i, 0))*(1.-S(i,3))+S(i,1)*S(i,2));//Z11
        Z(i, 1) = K*2.*S(i,1);//Z21
        Z(i, 2) = K*2.*S(i,2);//Z12
        Z(i, 3) = K*((1.-S(i, 0))*(1.+S(i,3))+S(i,1)*S(i,2));//Z22

        //cout << freq.at(i) << " " << Z.row(i) << endl;

    }


    DATA.freq = freq;
    DATA.S = S;
    DATA.Z = Z;

    return DATA;
}

// Load and source impedances may be sampled at different frequencies. It is essential to resample them
// using the same frequency basis. This requires interpolation of complex data. It would be desirable to use
// spline or cubic interpolation, but it seems that they are not implemented in Armadillo
int MatchingNetwork::ResampleImpedances()
{
    //Check if the inputs lie in the s2p measurements
    double fmin, fmax;
    (fS.min() > fL.min()) ? fmin = fS.min() : fmin = fL.min();
    (fS.max() > fL.max()) ? fmax = fL.max() : fmax = fS.max();

    int N = f_matching.n_rows;
    if (fmin > f_matching.min())
    {
        std::cout << "Warning: Data not found for "<< f_matching.min() << ". The lowest frequency was set to "<< fmin << std::endl;
        f_matching  = linspace(fmin,f_matching.max(), N);
    }
    if (fmax < f_matching.max())
    {
        std::cout << "Warning: Data not found for "<< f_matching.max() << ". The highest frequency was set to "<< fmax << std::endl;
        f_matching  = linspace(f_matching.min(), fmax, N);
    }

    double N_density = N/(f_matching.max() - f_matching.min());
    f_analysis = linspace(fmin, fmax, floor(N_density*(fmax-fmin)));


    //Impedance interpolation. This way both vector have the same frequency basis.
    vec ZS_r(real(ZS)), ZS_i(imag(ZS)), ZL_r(real(ZL)), ZL_i(imag(ZL));
    vec ZS_inter_R, ZS_inter_I, ZL_inter_R, ZL_inter_I;
    interp1(fS, ZS_r, f_analysis, ZS_inter_R);//Armadillo lacks of spline interpolation
    interp1(fS, ZS_i, f_analysis, ZS_inter_I);
    interp1(fL, ZL_r, f_analysis, ZL_inter_R);
    interp1(fL, ZL_i, f_analysis, ZL_inter_I);

    ZS = cx_vec(ZS_inter_R, ZS_inter_I);
    ZL = cx_vec(ZL_inter_R, ZL_inter_I);

    //Select impedance values at the matching band
    vec V = abs(f_analysis - f_matching.min());
    uvec v_find = find(V < V.min()+1e-3);
    int index1 = v_find.at(0);
    V = abs(f_analysis - f_matching.max());
    v_find = find(V < V.min()+1e-3);
    int index2 = v_find.at(0);

    ZS_matching = ZS.rows(index1, index2);
    ZL_matching = ZL.rows(index1, index2);

    //Now, it's needed to u-pdate the matching band according to f_analysis values
    f_matching = f_analysis.rows(index1, index2);
    return 0;
}


// Returns the S matrix at a given frequency
cx_mat MatchingNetwork::getSparams(rowvec x, cx_double zs, cx_double zl, double f)
{
    cx_mat ABCD = getABCDmatrix(x, f);
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
cx_mat MatchingNetwork::getABCDmatrix(rowvec x, double f)
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

// Sets the maximum number of iterations for the grid search engine
int MatchingNetwork::SetMaxIterGridSearch(int max_iter)
{
    Grid_MaxIter = max_iter;
    return 0;
}

// Sets the threshold of what the program consideres a good matching
int MatchingNetwork::SetThreshold(double th)
{
    MatchingThreshold = th;
    return 0;
}

// Returns the matching threshold
double MatchingNetwork::GetThreshold()
{
    return MatchingThreshold;
}


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
        x_ini << 5e-12 << 5e-12 << 5e-12 << 5e-9 << 5e-9 << 5e-9;
        Vaux = GridSearch();
        gridtest = CandidateEval(Vaux);
        candidate = topology;
        opttopo = gridtest;
        Vopt = Vaux;
        cout << "PiCPiL: S11_min = " <<  gridtest << " dB" << endl;


        topology = "202313";//PiLPiC
        x_ini << 5e-9 << 5e-9 << 5e-9 << 5e-12 << 5e-12 << 5e-12;
        Vaux = GridSearch();
        gridtest = CandidateEval(Vaux);
        cout << "PiLPiC: S11_min = " <<  gridtest << " dB" << endl;
        if (gridtest < opttopo - 0.5)//It worths to change the topology if the new one improves the result significantly
        {
            candidate = topology;
            opttopo = gridtest;
            Vopt = Vaux;
        }

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


        topology = "020313";
        x_ini << 5e-9 << 5e-9 << 5e-9 << 5e-12 << 5e-12 << 5e-12;
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
        x_ini  << 5e-12 << 5e-12 << 5e-12 << 5e-9 << 5e-9 << 5e-9;
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
        x_ini << 5e-9 << 5e-12 << 5e-9 << 5e-12 << 5e-9;
        Vaux = GridSearch();
        gridtest = CandidateEval(Vaux);
        cout << "Lowpass LC ladder: S11_min = " <<  gridtest << " dB" << endl;
        if (gridtest < opttopo - 0.5)
        {
            candidate = topology;
            opttopo = gridtest;
            Vopt = Vaux;
        }

        topology = "43434";
        x_ini << 75 << lambda4 << 5e-12 << 75 << lambda4  << 5e-12 << 75 << lambda4 ;
        Vaux = GridSearch();
        gridtest = CandidateEval(Vaux);
        cout << "Lowpass TL-C ladder: S11_min = " <<  gridtest << " dB" << endl;
        if (gridtest < opttopo - 0.5)
        {
            candidate = topology;
            opttopo = gridtest;
            Vopt = Vaux;
        }


        topology = "01230123";//BPP3
        x_ini << 5e-9 << 5e-12 << 5e-9 << 5e-12<< 5e-9 << 5e-12<< 5e-9 << 5e-12;
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
        x_ini << 5e-9 << 5e-12 << 5e-9 << 5e-12<< 5e-9 << 5e-12<< 5e-9 << 5e-12;
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
        x_ini << 5e-12 << 5e-12 << 5e-12 << 5e-9 << 5e-9 << 5e-9 << 5e-12 << 5e-12 << 5e-12;
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



rowvec MatchingNetwork::GridSearch()
{
    vec delta_k, best, step;
    int dim = x_ini.n_cols;
    rowvec f_xkq(pow(2,dim)+1);
    rowvec xk(x_ini);
    mat xkq = ones(pow(2, dim)+1, dim);
    mat best_pivot = ones(4, dim);
    delta_k << 1. << 1./4 << 1./16 << 1./64 << endr;
    best << 100<< 100<<100<<100<< endr;
    double fxk = CandidateEval(xk);
    double best_candidate;
    best.at(0) = fxk;
    mat C = GeneratingMatrix(dim);
    unsigned int niter = 0;


    /*
    //Initial guess. Sometimes, the initial point may be far away from the optimum, so it
    // seems to be reasonable to pick a handful of random vectors nearby the inital point and
    // choose the best one.

    int N_initial_guess  = 100;
    mat XKQ = ones(N_initial_guess, xk.n_cols);
    vec FXK = ones(N_initial_guess);
    for (int k = 0; k <N_initial_guess; k++)
    {
        XKQ.row(k) = xk %(randu(1, xk.n_cols)*3);
        FXK.at(k) =CandidateEval(XKQ.row(k));
    }
    uvec v_find_guess = find(FXK < FXK.min()+1e-6);
    int imin = v_find_guess.at(0);
    if (FXK.min() < fxk)
    {
        xk = XKQ.row(imin);
        fxk = CandidateEval(xk);
    }*/

    // Factorial search
    for (int i = 0; i<4; i++, niter++)
    {
        if ((best.min() < MatchingThreshold)&&(niter > Grid_MaxIter))
        {
            return xk;
        }

        for (int q = 0; q<=pow(2, dim);q++)
        {
            step = delta_k.at(i)*C.col(q);
            xkq.row(q) = xk%(1. + step.t()*.99);
            f_xkq.at(q) = CandidateEval(xkq.row(q));
        }
        best_candidate = f_xkq.min();
        uvec v_find = find(f_xkq < best_candidate+1e-6);
        int imin = v_find.at(0);

        if (verbose)std::cout << "BEST CANDIDATE "<< i << ": " << xkq.row(imin) << " => " << best_candidate << std::endl;


        if (best_candidate < best.min())
        {
            best << best_candidate << 0 << 0 << 0 << endr;//Prepare best for new hypercube cycle
            xk = xkq.row(imin);//Set new pivot
            xk = InspectCandidate(xk);
            i = -1;//New pivot
            if (verbose) std::cout << "New pivot: " << xk << std::endl;
            continue;
        }
        else
        {
            best.at(i) = best_candidate;
            best_pivot.row(i) = xk;
        }

        if (i > 0)//Checks if the i-th cycle improved the result with regard to the (i-1)-th cycle
        {
            if ((abs(best.at(i)-best.at(i-1))<1e-3)||(best(i)>best.at(i-1)))//The current cycle did not improve the result
            {
                uvec v_find = find(best < best.min()+1e-6);

                //Suggestion for improvement. When GRABIM gets stuck, sometimes it helps to try some
                //random vectors near the current pivot. The distribution was chosen to be U(-0.5, 0.5)
                int Ncatchup  = 100;
                mat XKQ = ones(Ncatchup, xk.n_cols);
                vec FXK = ones(Ncatchup);
                for (int k = 0; k <Ncatchup; k++)
                {
                    XKQ.row(k) = best_pivot.row(0) %(randu(1, xk.n_cols)+.5);
                    FXK.at(k) =CandidateEval(XKQ.row(k));
                }
                uvec v_find_catchup = find(FXK < FXK.min()+1e-6);
                int imin = v_find_catchup.at(0);
                if (FXK.min() < best.min())
                {
                    xk = XKQ.row(imin);
                    xk = InspectCandidate(xk);
                    if (verbose)std::cout << "=> " << xk << ": " << FXK.min()<< std::endl;
                    i=-1;
                }
            }

        }
    }
    return xk;
}


// This function checks whether the candidate vector contains some
// irrelevant value or not. Unsignificant values may prevent the algorithm
// to find a better or more realistic solution. In this sense, shunt elements
// whose impedance exceeds 4kOhm are interpreted as open circuits, and conversely,
// series components whose impedance is below 1 Ohm are treated as short circuits
rowvec MatchingNetwork::InspectCandidate(rowvec xk)
{
    double impedance, fmax = f_matching.max(), fmin = f_matching.min();
    unsigned int element;
    double wmax =2*datum::pi*fmax, wmin =2*datum::pi*fmin;
    for (unsigned int i = 0; i < topology.length(); i++)
    {
        element = atoi(topology.substr(i,1).c_str());
        switch(element)
        {
        case 0://Series inductor
            impedance = wmax*xk.at(i);
            if (impedance < 1)
            {
                xk.at(i) = 1e-30;
            }
            if (xk.at(i) < 0)//Inductance must be > 0, so it seems that a capacitor would do a better job
            {
                xk.at(i) = 1./(2e3*wmin);//High impedance capacitor
                topology[i] = '1';//Series capacitor
                cout << "Warning: The selected topology leads to L < 0" <<endl;
                cout << "Warging: Topology changed. The new topology is: " << topology << endl;
            }
            continue;
        case 2://Shunt inductor
            impedance = wmin*xk.at(i);
            if (impedance > 4e3)
            {
                xk.at(i) = 100;
            }
            if (xk.at(i) < 0)//Inductance must be > 0, so it seems that a capacitor would do a better job
            {
                xk.at(i) = 1./(2e3*wmin);//High impedance capacitor
                topology[i]='3';//Shunt capacitor
                cout << "Warning: The selected topology leads to L < 0" <<endl;
                cout << "Warging: Topology changed. The new topology is: " << topology << endl;
            }

            continue;
        case 1://Series capacitor
            impedance = 1./(wmin*xk.at(i));
            if (impedance < 1)
            {
                xk.at(i) = 100;
            }
            if (xk.at(i) < 0)//C must be > 0, so it seems that an inductor would do a better job
            {
                xk.at(i) = 5./(wmax);//Low impedance inductance
                topology[i] = '0';//Series inductance
                cout << "Warning: The selected topology leads to C < 0" <<endl;
                cout << "Warging: Topology changed. The new topology is: " << topology << endl;
            }
            continue;
        case 3://Shunt capacitor
            impedance = 1./(wmax*xk.at(i));
            if (impedance > 4e3)
            {xk.at(i) = 1e-30;}
            if (xk.at(i) < 0)//C must be > 0, so it seems that an inductor would do a better job
            {
                xk.at(i) = 5./(wmax);//Low impedance inductance
                topology[i] = '2';//Shunt inductance
                cout << "Warning: The selected topology leads to C < 0" <<endl;
                cout << "Warging: Topology changed. The new topology is: " << topology << endl;
            }
            continue;

        }
    }
  return xk;
}


mat MatchingNetwork::GeneratingMatrix(int dim)
{
    mat C = zeros(dim, pow(2, dim)+1);
    for (int i = 0; i <dim; i++)
    {
        mat V =  ones(1, pow(2, (i+1)));
        V.cols(0, pow(2, i)-1) = ones(1, pow(2, i));
        V.cols(pow(2,i), pow(2, i+1)-1) = -ones(1, pow(2, i));
        C.row(i).cols(0, pow(2, dim)-1) = repmat(V, 1., pow(2, dim-i-1));
    }
    return C;
}

double MatchingNetwork::CandidateEval(rowvec x)
{
    double fobj = -1e3;
    cx_mat S, ABCD;
    for (unsigned int i = 0; i < f_matching.n_rows; i++)
    {

        if (ObjFun == ObjectiveFunction::NINF_S11dB)
        {
            S = getSparams(x, ZS_matching.at(i), ZL_matching.at(i), f_matching.at(i));
            if (abs(S(0,0)) > fobj) fobj = abs(S(0,0));
        }
        if (ObjFun == ObjectiveFunction::NINF_POWERTRANS)
        {
            ABCD = getABCDmatrix(x, f_matching.at(i));
            fobj = CalcInvPowerTransfer(ABCD, ZS_matching.at(i), ZL_matching.at(i));
        }
    }
    if (ObjFun == ObjectiveFunction::NINF_S11dB)fobj = 20*log10(fobj);//|grad{log(x)}| > |grad{x}| when x < 1;
    return fobj;
}


typedef struct NLoptData {
    vec f_matching;
    cx_vec ZS_matching, ZL_matching;
    std::string topology;
    ObjectiveFunction objf;

} NLoptData;

double myfunc(const std::vector<double> &x, std::vector<double> &grad, void *n)
{
    MatchingNetwork M;
    rowvec x_ = ones(1,x.size());
    NLoptData * N = (NLoptData *) n;
    M.SetTopology(N->topology);
    M.SetMatchingBand(N->f_matching.min(), N->f_matching.max(), N->f_matching.n_rows);
    M.SetSourceImpedance(N->ZS_matching, N->f_matching);
    M.SetLoadImpedance(N->ZL_matching, N->f_matching);
    M.SetObjectiveFunction(N->objf);
    for (unsigned int i = 0; i < x.size(); i++) x_.at(i) = x[i];

    double eval = M.CandidateEval(x_);
    std::cout<< eval << " <= " << x_<< std::endl;
    return eval;
}

rowvec MatchingNetwork::LocalOptimiser(rowvec x_grid)
{
    int dim = x_grid.n_cols;
    nlopt::opt opt(NLoptAlgo, dim);
    NLoptData n;
    n.f_matching = f_matching;
    n.topology = topology;
    n.ZL_matching = ZL_matching;
    n.ZS_matching = ZS_matching;
    n.objf = ObjFun;
    opt.set_min_objective(myfunc, &n);
    opt.set_stopval(this->GetThreshold()-10);
    opt.set_maxeval(1e6);
    opt.set_xtol_rel(1e-4);

    //Bounds
    std::vector<double> lb(dim), ub(dim);
    unsigned int i, k;
    for (i = 0, k=0; i < topology.length(); i++, k++)
    {
        int element = atoi(topology.substr(i,1).c_str());
        if (element < 4)
        {
            lb[k] = 1e-20;
            ub[k] = 1e-6;
        }
        else
        {
            lb[k] = 5;   lb[k+1] = 1e-6;
            ub[k] = 1e4; ub[k+1] = 10;
            k++;
        }
    }


    // opt.set_lower_bounds(lb);
    // opt.set_upper_bounds(ub);

    std::vector<double> x(dim);
    for (int i = 0; i < dim; i++)x[i] = x_grid.at(i);
    double minf;

    cout << endl << "Local optimiser output:" << endl;
    nlopt::result result = opt.optimize(x, minf);
    std::cout << "NLopt status: " << std::endl;

    switch(result)
    {
    //Success
    case 1: std::cout << "Success!" <<std::endl;
        break;
    case 2: std::cout << "Goal achieved" <<std::endl;
        break;
    case 3: std::cout << "The algorithm reached convergence (relative)" <<std::endl;
        break;
    case 4: std::cout << "The algorithm reached convergence (absolute)" <<std::endl;
        break;
    case 5: std::cout << "Maximum number of evaluations reached" <<std::endl;
        break;
    case 6: std::cout << "Maximum time reached" <<std::endl;
        break;
        //Failure
    case -1: std::cout << "NLopt failed!" <<std::endl;
        break;
    case -2: std::cout << "Invalid arguments" <<std::endl;
        break;
    case -3: std::cout << "There is not enough available memory" <<std::endl;
        break;
    case -4: std::cout << "Forced termination" <<std::endl;
        break;
    }

    rowvec x_nlopt = ones(1, dim);
    for (int i = 0; i < dim; i++)x_nlopt.at(i) = x[i];

    return x_nlopt;

}


int MatchingNetwork::SetNLoptAlg(nlopt::algorithm NLA)
{
    NLoptAlgo = NLA;
    return 0;
}

nlopt::algorithm MatchingNetwork::GetNLoptAlg()
{
    return NLoptAlgo;
}


int MatchingNetwork::SetObjectiveFunction(ObjectiveFunction of)
{
    ObjFun = of;
    return 0;
}

ObjectiveFunction MatchingNetwork::GetObjectiveFunction()
{
    return ObjFun;
}

// Inverse power transfer from ABCD matrix
// [1] Eq (6.2.1), (6.2.2)
double MatchingNetwork::CalcInvPowerTransfer(cx_mat ABCD, cx_double ZS, cx_double ZL)
{
    cx_double p = real(ZS)*real(ZL) - imag(ZS)*imag(ZL);
    cx_double q = imag(ZS)*real(ZL) + imag(ZL)*real(ZS);
    cx_double a = ABCD(0,0)*real(ZL) - ABCD(1,0)*q + ABCD(1,1)*real(ZS);
    cx_double b = ABCD(0,1) + ABCD(1,0)*p + ABCD(1,1)*imag(ZS) + ABCD(0,0)*imag(ZL);
    return abs((a*a + b*b)/(4*real(ZS)*real(ZL)));
}
