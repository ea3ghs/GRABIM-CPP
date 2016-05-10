#include "io.h"

IO::IO()
{
    Nsamples = 30;
}

double string_to_double( const std::string& s )
{
    std::istringstream i(s);
    double x;
    if (!(i >> x))
        return 0;
    return x;
}


//Export data to display results with GNUplot
int IO::exportGNUplot(GRABIM_Result Res, string filepath)
{
    ofstream GNUplotExport;
    GNUplotExport.open (filepath, std::ofstream::out);
    if(!GNUplotExport.is_open()) return -1;
    GNUplotExport << "#GRABIM data" << endl;
    GNUplotExport << "#freq S11(Grid search) S11(Local optimiser) S21(Grid search) S21(Local optimiser) " << endl;
    for (unsigned int i = 0; i < Res.freq.n_rows; i++)
    {//Writes frequency S11_grid[db] S11_optimiser[db] and S21[dB]
        GNUplotExport << Res.freq.at(i)*1e-9 << " " << 20*log10(abs(Res.S11_gridsearch.at(i))) << " "
                      << 20*log10(abs(Res.S21_gridsearch.at(i))) << " " << 20*log10(abs(Res.S11_nlopt.at(i)))
                      << " " << 20*log10(abs(Res.S21_nlopt.at(i))) << " " << real(Res.S11_nlopt.at(i))
                      << " " << imag(Res.S11_nlopt.at(i)) <<endl;
    }
    GNUplotExport.close();

    //Update plotscript
    ifstream plotscriptfile("plotscript");
    if (!plotscriptfile.is_open())
    {
        cout << "Something went wrong with file export..." << endl;
        return -1;
    }
    string s = "source = \"" + filepath + "\"\n", line;//Path to the data dile
    getline(plotscriptfile, line);
    while (getline(plotscriptfile, line))s+=line+"\n";
    plotscriptfile.close();

    ofstream plotscript_write("plotscript");//Updates the name of the source data
    plotscript_write << s;
    plotscript_write.flush();
    plotscript_write.close();

    return 0;
}


// Extends std::tolower(int c) capability to st::string arguments
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

//Removes consecutive blank spaces in a std::string
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


int IO::loadS1Pdata(std::string filepath, terminal Port, bool isRawImpedance)
{
    std::ifstream s2pfile(filepath.c_str());//Tries to open the data file.
    if(!s2pfile.is_open())//The data file cannot be opened => error
    {
        return -1;
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
    freq_scale = getS2PfreqScale(line);


    //Get the impedance at which the S params were measured

    int Rindex = line.find_last_of("r");
    Rindex = line.find_first_not_of(" ", Rindex);
    Zref = atof(line.substr(Rindex+1).c_str());
    bool is_indB = (line.find("db") != -1);
    bool RI = (line.find("ma") == -1);


    while( getline(s2pfile, line) )
    {//Looking for the start of the raw data

        line = RemoveBlankSpaces(line);

        if ((!line.compare(0,1, "!"))|| (line.length() == 1)) continue;
        else break;


    }

    //DATA beginning.
    //At this point, the number of frequency samples is not known, so it's better to
    //push data into queues and then arrange it into armadillo structures
    std::queue <double> frequency, S11M, S11A;
    unsigned int qsize=0;

    do
    {
        line = RemoveBlankSpaces(line);
        if (line.empty()|| (line.length()==1))break;
        if (line.at(0) == '!') break;//Comment

        //Frequency
        int index = line.find_first_of(" ");
        if (index == -1)index = line.find_first_of("\t");
        frequency.push(string_to_double(line.substr(0,index)));
        line.erase(0, index+1);


        index = line.find_first_of(" ");
        if (index == -1)index = line.find_first_of("\t");
        S11M.push(string_to_double(line.substr(0,index)));
        line.erase(0, index+1);

        index = line.find_first_of(" ");
        if (index == -1)index = line.find_first_of("\t");
        S11A.push(string_to_double(line.substr(0,index)));

        qsize++;
    }while (std::getline(s2pfile, line));

    vec freq(qsize);
    cx_mat S(qsize, 1);
    cx_mat Z(qsize, 1);
    double phi, S11m, S11a;
    for (unsigned int i = 0; i < qsize; i++)
    {
        freq(i) = freq_scale*frequency.front();
        frequency.pop();

        S11m = S11M.front();
        S11a = S11A.front();
        S11M.pop();
        S11A.pop();
        if (is_indB) S11m = pow(10, .05*S11m);
        phi = (datum::pi/180)*S11a;

        if (RI)
        {
            S(i, 0) = cx_double(S11m, S11a);
        }
        else
        {
            S(i, 0) = cx_double(S11m,0)*cx_double(cos(phi), sin(phi));
        }

        if (isRawImpedance)
        {//The data file contains impedance data
            Z(i, 0) =  S(i, 0);
        }
        else
        {//The data file contains s1p data
        Z(i, 0) = Zref*((1.+S(i, 0))/(1.-S(i, 0)));//Z
        }

    }

    if (Port == SOURCE)//Set source port properties
    {
        fS = freq;
        ZS = Z.col(0);
        if (!ZL.empty())ResampleImpedances();
    }

    if (Port == LOAD)//Set load port properties
    {
        fL = freq;
        ZL = Z.col(0);;
        if (!ZS.empty())ResampleImpedances();
    }
    return 0;
}


// Load and source impedances may be sampled at different frequencies. It is essential to resample them
// using the same frequency basis. This requires interpolation of complex data. It would be desirable to use
// spline or cubic interpolation, but it seems that they are not implemented in Armadillo
int IO::ResampleImpedances()
{
    //Check whether the vectors are already resampled or not
    if ((ZS.n_elem == ZL.n_elem) && (ZS.n_elem == freq.n_elem))
    {
        return 0;
    }

    double fmin, fmax;
    /*We have to discern the following cases
    1: ZS constant, ZL constant versus frequency
    2: ZS s2p/s1p, ZL constant vs freq
    3: ZS constant vs freq, ZS s2p/s1p
    4: ZS s2p/s1p vs ZL s2p/s1p
    */

    if ((ZS.n_elem == 1) && (ZL.n_elem == 1)) //ZS constant, ZL constant vs frequency
    {
        freq = linspace(fmatching_min, fmatching_max, Nsamples);//Available freqs
        cx_double zs_temp = ZS(0,0);
        cx_double zl_temp = ZL(0,0);
        //Create ZS and ZL vectors with the same length than freq
        ZS.resize(freq.size());
        ZS.ones();
        ZS = ZS*zs_temp;

        ZL.resize(freq.size());
        ZL.ones();
        ZL = ZL*zl_temp;
        return 0;
    }

    if ((ZS.n_elem == 1) && (ZL.n_elem != 1)) //ZS constant vs ZL s2p/s1p
    {
        fS = fL;
        freq = fL;
        cx_double zs_temp = ZS(0,0);
        ZS.resize(ZL.size());
        ZS.ones();
        ZS = ZS*zs_temp;
        fmin = freq.min();
        fmax = freq.max();
    }

    if ((ZS.n_elem != 1) && (ZL.n_elem == 1)) //ZS s2p/s1p vs ZL constant
    {
        fL = fS;
        freq = fL;
        cx_double zl_temp = ZL(0,0);
        ZL.resize(ZS.size());
        ZL.ones();
        ZL = ZL*zl_temp;
        fmin = freq.min();
        fmax = freq.max();
    }
    else
    {
      if ((ZS.n_elem != 1) && (ZL.n_elem != 1)) //ZS s2p/s1p vs ZL s2p/s1p
      {
          //Define vector of available freqs
          (fS.min() > fL.min()) ? fmin = fS.min() : fmin = fL.min();
          (fS.max() > fL.max()) ? fmax = fL.max() : fmax = fS.max();

          //int N;//Number of points for linear interpolation
          //(ZS.n_elem > ZL.n_elem) ? N = 2*ZS.n_elem : N = 2*ZL.n_elem;
      }
    }

        freq = linspace(fmin, fmax, Nsamples);//Available freqs

        //Impedance interpolation. This way both vector have the same frequency basis.
        vec ZS_r(real(ZS)), ZS_i(imag(ZS)), ZL_r(real(ZL)), ZL_i(imag(ZL));
        vec ZS_inter_R, ZS_inter_I, ZL_inter_R, ZL_inter_I;
        interp1(fS, ZS_r, freq, ZS_inter_R);//Armadillo lacks of spline interpolation
        interp1(fS, ZS_i, freq, ZS_inter_I);
        interp1(fL, ZL_r, freq, ZL_inter_R);
        interp1(fL, ZL_i, freq, ZL_inter_I);

        ZS = cx_vec(ZS_inter_R, ZS_inter_I);
        ZL = cx_vec(ZL_inter_R, ZL_inter_I);
        return 0;

}


int IO::loadS2Pdata(std::string filepath, terminal Port)
{
    std::ifstream s2pfile(filepath.c_str());//Tries to open the data file.
    if(!s2pfile.is_open())//The data file cannot be opened => error
    {
        return -1;
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
    freq_scale = getS2PfreqScale(line);

    //Get the impedance at which the S params were measured

    int Rindex = line.find_last_of("r");
    Rindex = line.find_first_not_of(" ", Rindex);
    Zref = atof(line.substr(Rindex+1).c_str());
    int is_indB = line.find("db");
    bool RI = (line.find("ma") == -1);

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
        S11M.push(string_to_double(line.substr(0,index)));
        line.erase(0, index+1);

        index = line.find_first_of(" ");
        if (index == -1)index = line.find_first_of("\t");
        S11A.push(string_to_double(line.substr(0,index)));
        line.erase(0, index+1);

        index = line.find_first_of(" ");
        if (index == -1)index = line.find_first_of("\t");
        S21M.push(string_to_double(line.substr(0,index)));
        line.erase(0, index+1);

        index = line.find_first_of(" ");
        if (index == -1)index = line.find_first_of("\t");
        S21A.push(string_to_double(line.substr(0,index)));
        line.erase(0, index+1);

        index = line.find_first_of(" ");
        if (index == -1)index = line.find_first_of("\t");
        S12M.push(string_to_double(line.substr(0,index)));
        line.erase(0, index+1);

        index = line.find_first_of(" ");
        if (index == -1)index = line.find_first_of("\t");
        S12A.push(string_to_double(line.substr(0,index)));
        line.erase(0, index+1);

        index = line.find_first_of(" ");
        if (index == -1)index = line.find_first_of("\t");
        S22M.push(string_to_double(line.substr(0,index)));
        line.erase(0, index+1);

        index = line.find_first_of(" ");
        if (index == -1)index = line.find_first_of("\t");
        if (index == -1)index = line.find_first_of("\r");//End of the line
        S22A.push(string_to_double(line.substr(0,index)));
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

        if (RI)
        {
            S(i, 0) = cx_double(S11m, S11a);
        }
        else
        {
            phi = (datum::pi/180)*S11a;
            S(i, 0) = cx_double(S11m,0)*cx_double(cos(phi), sin(phi));
        }


        S21m = S21M.front();
        S21a = S21A.front();
        S21M.pop();
        S21A.pop();
        if (is_indB != -1) S21m = pow(10, .05*S21m);

        if (RI)
        {
            S(i, 1) = cx_double(S21m, S21a);
        }
        else
        {
            phi = (datum::pi/180)*S21a;
            S(i, 1) = cx_double(S21m,0)*cx_double(cos(phi), sin(phi));
        }


        S12m = S12M.front();
        S12a = S12A.front();
        S12M.pop();
        S12A.pop();
        if (is_indB != -1) S12m = pow(10, .05*S12m);

        if (RI)
        {
            S(i, 2) = cx_double(S12m, S12a);
        }
        else
        {
            phi = (datum::pi/180)*S12a;
            S(i, 2) = cx_double(S12m,0)*cx_double(cos(phi), sin(phi));
        }

        S22m = S22M.front();
        S22a = S22A.front();
        S22M.pop();
        S22A.pop();
        if (is_indB != -1) S22m = pow(10, .05*S22m);

        if (RI)
        {
            S(i, 3) = cx_double(S22m, S22a);
        }
        else
        {
            phi = (datum::pi/180)*S22a;
            S(i, 3) = cx_double(S22m,0)*cx_double(cos(phi), sin(phi));
        }


        cx_double K = cx_double(Zref,0)/((1.-S(i,0))*(1.-S(i,3))-S(i,1)*S(i,2));
        Z(i, 0) = K*((1.+S(i, 0))*(1.-S(i,3))+S(i,1)*S(i,2));//Z11
        Z(i, 1) = K*2.*S(i,1);//Z21
        Z(i, 2) = K*2.*S(i,2);//Z12
        Z(i, 3) = K*((1.-S(i, 0))*(1.+S(i,3))+S(i,1)*S(i,2));//Z22

        //cout << freq.at(i) << " " << Z.row(i) << endl;

    }

    if (Port == SOURCE)//Set source port properties
    {
        fS = freq;
        ZS = Z.col(3);//The source port is the output of a prior stage
        ResampleImpedances();
    }

    if (Port == LOAD)//Set load port properties
    {
        fL = freq;
        ZL = Z.col(0);//The load port is the input of the following stage
        ResampleImpedances();
    }

    return 0;
}



cx_mat IO::getSourceImpedance()
{
    ResampleImpedances();
    return ZS;
}
cx_mat IO::getLoadImpedance()
{
    ResampleImpedances();
    return ZL;
}

vec IO::getFrequency()
{
    ResampleImpedances();
    return freq;
}


void IO::set_constant_ZS_vs_freq(cx_double zs)
{
    ZS = zs;
}


void IO::set_constant_ZL_vs_freq(cx_double zl)
{
    ZL = zl;
}


void IO::set_matching_band(double fmin, double fmax)
{
    fmatching_min = fmin;
    fmatching_max = fmax;
}

// Set load/source matching impedace
int IO::setMatchingImpedances()
{
    //Select impedance values at the matching band
    int index1 = getFreqIndex(fmatching_min);
    int index2 = getFreqIndex(fmatching_max);

    ZS_matching = ZS.rows(index1, index2);
    ZL_matching = ZL.rows(index1, index2);

    //Now, it's needed to update the matching band according to f_analysis values
    f_matching = freq.rows(index1, index2);
    return 0;
}


// Gets the index of a given frequency
int IO::getFreqIndex(double f1)
{
    vec V = abs(freq - f1);
    uvec v_find = find(V < V.min()+1e-3);
    return v_find.at(0);
}

//Sets local optimisation algorithm
void IO::setLocalOptimiser(nlopt::algorithm NL)
{
    LocalOptAlgo = NL;
}


nlopt::algorithm IO::getLocalOptimiser()
{
    return LocalOptAlgo;
}


//Get freq scale from a string line
double IO::getS2PfreqScale(string line)
{
    if (line.find("ghz") != std::string::npos)
    {
        return 1e9;
    }
    else
    {
        if (line.find("mhz") != std::string::npos)
        {
            return 1e6;
        }
        else
        {
            if ((line.find("khz") != std::string::npos))
            {
                return 1e3;
            }
            else
            {
                return 1;//Hz
            }
        }
    }
}


// This function exports the best topology found to a Qucs schematic
int IO::ExportQucsSchematic(GRABIM_Result R)
{
    QString wirestr = "";
    QString componentstr = "";
    QString paintingstr = "";
    int x_pos = 0;
    SchematicParser(R, x_pos, componentstr, wirestr, paintingstr);
    CreateSchematic(componentstr, wirestr, paintingstr, R.QucsVersion);
    return 0;
}



// Given a string code of inductors, capacitors and transmission lines, it generates the Qucs network. Notice that the schematic is split into
// three part: components, wires and paintings, all of them are passed by reference.
int IO::SchematicParser(GRABIM_Result R, int & x_pos, QString & componentstr, QString & wirestr, QString & paintingstr)
{
    QString component;
    int x_series = 120, x_shunt = 20;//x-axis spacing depending on whether the component is placed in a series or shunt configuration
    //Clear input strings (just in case)
    componentstr = "";
    wirestr = "";
    paintingstr = "";

    // Schematic code
    // 0: Port 1 + S parameter simulation
    // 1: Port 1
    // 2: Port 2
    // 3: Port 1, Port 2 and S parameter simulation


    if ((R.source_path.isEmpty()) && (abs(ZS.at(0).imag()) < 1e-3) && (ZS.at(0).real() > 1e-3))
    {//Conventional term
        componentstr += QString("<Pac P1 1 %2 -30 18 -26 0 1 \"1\" 1 \"%1 Ohm\" 1 \"0 dBm\" 0 \"1 GHz\" 0>\n").arg(ZS.at(0).real()).arg(x_pos);
        componentstr += QString("<GND * 1 %1 0 0 0 0 0>\n").arg(x_pos);

        wirestr += QString("<%1 -60 %1 -120>\n").arg(x_pos);
        wirestr += QString("<%1 -120 %2 -120>\n").arg(x_pos).arg(x_pos+60);
        x_pos +=60;
    }
    else
    {//Place a S-param file
        if (R.source_path.contains(".s1p"))
        {
           componentstr += QString("<SPfile X1 1 %1 -120 -26 -67 0 0 \"%2\" 1 \"rectangular\" 0 \"linear\" 0 \"open\" 0 \"1\" 0>").arg(x_pos).arg(R.source_path);
           componentstr += QString("<GND * 1 %1 -150 0 0 0 0>").arg(x_pos);
           x_pos = 30;
           wirestr += QString("<%1 -120 %2 -120>\n").arg(x_pos).arg(x_pos+60);
           x_pos +=60;
        }
        if (R.source_path.contains(".s2p"))
        {
           componentstr += QString("<SPfile X1 1 %1 -120 -26 -67 0 0 \"%2\" 1 \"rectangular\" 0 \"linear\" 0 \"open\" 0 \"2\" 0>").arg(x_pos).arg(R.source_path);
           componentstr += QString("<GND * 1 %1 -150 0 0 0 0>").arg(x_pos);
           x_pos = 30;
           wirestr += QString("<%1 -120 %2 -120>\n").arg(x_pos).arg(x_pos+60);
           x_pos +=60;
        }
    }



    //Add the frequency range for the S-param simulation
    double freq_start = std::max(0., freq.min()-5e8);
    double freq_stop = freq.max()+5e8;
    componentstr += QString("<.SP SP1 1 0 100 0 67 0 0 \"lin\" 1 \"%2Hz\" 1 \"%3Hz\" 1 \"300\" 1 \"no\" 0 \"1\" 0 \"2\" 0>\n").arg((freq_start)).arg((freq_stop));
    componentstr += QString("<Eqn Eqn1 1 200 100 -28 15 0 0 \"S11_dB=dB(S[1,1])\" 1 \"yes\" 0>\n");
    componentstr += QString("<Eqn Eqn1 1 200 200 -28 15 0 0 \"S21_dB=dB(S[2,1])\" 1 \"yes\" 0>\n");


    // The string format is as follows: "XX<value>;XX<value2>;...XX<valueN>;"
    // where XX, YY, ZZ define the type of component and its configuration.
    //    0: Series inductance
    //    1: Series capacitor
    //    2: Shunt inductance
    //    3: Shunt capacitor
    //    4: Series transmission line
    //    5: Open stub
    //    6: Short circuited stub

    int value_index = 0;
    for (int i = 0; i < R.topology.count(); i++)
    {
        component = R.topology.at(i);

        if (!component.compare("0"))//Series inductor
        {
            componentstr += QString("<L L1 1 %1 -120 -26 10 0 0 \"%2\" 1 "" 0>\n").arg(x_pos+60).arg(R.x_nlopt.at(value_index));
            wirestr += QString("<%1 -120 %2 -120 "" 0 0 0 "">\n").arg(x_pos).arg(x_pos+30);
            wirestr += QString("<%1 -120 %2 -120 "" 0 0 0 "">\n").arg(x_pos+90).arg(x_pos+x_series);
            x_pos += x_series;
            value_index++;
        }
        else if (!component.compare("1"))//Series capacitor
        {
            componentstr += QString("<C C1 1 %1 -120 -26 17 0 0 \"%2\" 1 "" 0>\n").arg(x_pos+60).arg(R.x_nlopt.at(value_index));
            wirestr += QString("<%1 -120 %2 -120 "" 0 0 0 "">\n").arg(x_pos).arg(x_pos+30);
            wirestr += QString("<%1 -120 %2 -120 "" 0 0 0 "">\n").arg(x_pos+90).arg(x_pos+x_series);
            x_pos += x_series;
            value_index++;
        }
        else if (!component.compare("2"))//Shunt inductor
        {
            componentstr += QString("<GND * 1 %1 0 0 0 0 0>\n").arg(x_pos);
            componentstr += QString("<L L1 1 %1 -30 -26 2 0 1 \"%2\" 1 "" 0>\n").arg(x_pos).arg(R.x_nlopt.at(value_index));
            wirestr += QString("<%1 -60 %1 -120 "" 0 0 0 "">\n").arg(x_pos);
            wirestr += QString("<%1 -120 %2 -120 "" 0 0 0 "">\n").arg(x_pos).arg(x_pos+x_shunt);
            x_pos += x_shunt;
            value_index++;
        }
        else if (!component.compare("3"))//Shunt capacitor
        {
            componentstr += QString("<GND * 1 %1 0 0 0 0 0>\n").arg(x_pos);
            componentstr += QString("<C C1 1 %1 -30 -26 17 0 1 \"%2\" 1 "" 0>\n").arg(x_pos).arg(R.x_nlopt.at(value_index));
            wirestr += QString("<%1 -60 %1 -120 "" 0 0 0 "">\n").arg(x_pos);
            wirestr += QString("<%1 -120 %2 -120 "" 0 0 0 "">\n").arg(x_pos).arg(x_pos+x_shunt);
            x_pos += x_shunt;
            value_index++;
        }
        else if (!component.compare("4"))//Transmission line
        {
            componentstr += QString("<TLIN Line1 1 %3 -120 -26 20 0 0 \"%1\" 1 \"%2\" 1 \"0 dB\" 0 \"26.85\" 0>\n").arg(R.x_nlopt.at(value_index)).arg(R.x_nlopt.at(value_index+1)).arg(x_pos+60);
            wirestr += QString("<%1 -120 %2 -120 "" 0 0 0 "">\n").arg(x_pos).arg(x_pos+30);
            wirestr += QString("<%1 -120 %2 -120 "" 0 0 0 "">\n").arg(x_pos+90).arg(x_pos+x_series);
            x_pos += x_series;
            value_index+=2;
        }
        else if (!component.compare("5"))//Open stub
        {
            componentstr += QString("<TLIN Line1 1 %3 -60 -26 20 0 1 \"%1\" 1 \"%2\" 1 \"0 dB\" 0 \"26.85\" 0>\n").arg(R.x_nlopt.at(value_index)).arg(R.x_nlopt.at(value_index+1)).arg(x_pos);
            wirestr += QString("<%1 -90 %1 -120 "" 0 0 0 "">\n").arg(x_pos);
            wirestr += QString("<%1 -120 %2 -120 "" 0 0 0 "">\n").arg(x_pos).arg(x_pos+x_shunt);
            x_pos += x_shunt;
            value_index+=2;
        }
        else if (!component.compare("6"))//Short circuited stub
        {
            componentstr += QString("<TLIN Line1 1 %3 -60 -26 20 0 1 \"%1\" 1 \"%2\" 1 \"0 dB\" 0 \"26.85\" 0>\n").arg(R.x_nlopt.at(value_index)).arg(R.x_nlopt.at(value_index+1)).arg(x_pos);
            componentstr += QString("<GND * 1 %1 -30 0 0 0 0>\n").arg(x_pos);
            wirestr += QString("<%1 -90 %1 -120 "" 0 0 0 "">\n").arg(x_pos);
            wirestr += QString("<%1 -120 %2 -120 "" 0 0 0 "">\n").arg(x_pos).arg(x_pos+x_shunt);
            x_pos += x_shunt;
            value_index+=2;
        }

    }


    if ((R.load_path.isEmpty()) && (abs(ZL.at(0).imag()) < 1e-3) && (ZL.at(0).real() > 1e-3))
    {//Conventional term
        componentstr += QString("<Pac P1 1 %2 -30 18 -26 0 1 \"1\" 1 \"%1 Ohm\" 1 \"0 dBm\" 0 \"1 GHz\" 0>\n").arg(ZL.at(0).real()).arg(x_pos);
        componentstr += QString("<GND * 1 %1 0 0 0 0 0>\n").arg(x_pos);

        wirestr += QString("<%1 -60 %1 -120>\n").arg(x_pos);
    }
    else
    {//Place a S-param file
        if (R.source_path.contains(".s1p"))
        {
           componentstr += QString("<SPfile X1 1 %1 -120 -26 -67 0 0 \"%2\" 1 \"rectangular\" 0 \"linear\" 0 \"open\" 0 \"1\" 0>").arg(x_pos).arg(R.load_path);
           componentstr += QString("<GND * 1 %1 -150 0 0 0 0>").arg(x_pos);
        }
        if (R.source_path.contains(".s2p"))
        {
           componentstr += QString("<SPfile X1 1 %1 -120 -26 -67 0 0 \"%2\" 1 \"rectangular\" 0 \"linear\" 0 \"open\" 0 \"2\" 0>").arg(x_pos).arg(R.load_path);
           componentstr += QString("<GND * 1 %1 -150 0 0 0 0>").arg(x_pos);
        }
    }

    return 0;
}


//-----------------------------------------------------------------------------
// Given the components, wires and paintings, it creates the schematic and copies on the clipboard
bool IO::CreateSchematic(QString components, QString wires, QString paintings, QString QucsVersion)
{
    //Header
    QString Schematic = QString("<Qucs Schematic ") + QucsVersion + QString(">\n");

    //Add components
    Schematic += "<Components>\n";
    Schematic+=components;
    Schematic += "</Components>\n";

    //Add wires
    Schematic+= "<Wires>\n";
    Schematic += wires;
    Schematic+= "</Wires>\n";

    //Add paintings
    Schematic += "<Paintings>\n";
    Schematic += paintings;
    Schematic += "</Paintings>\n";

    QApplication::clipboard()->setText(Schematic, QClipboard::Clipboard);
    return true;
}
