#define MAXNETWORKS 39
struct _networks { char topology[16]; char name[24]; } NETWORKS[MAXNETWORKS]={
"303","PCLC",
"03","SLP1",
"0303","SLP2",
"030303","SLP3",
"30","PLP1",
"3030","PLP2",
"303030","PLP3",
"030","SLCL",
"12","SHP1",
"121","SCLC",
"1212","SHP2",
"121212","SHP3",
"21","PHP1",
"2121","PHP2",
"212121","PHP3",
"231023","PBP3",
"2310","PBP2",
"0123","SBP2",
"013201","SBP3",
"32023","BP+L+BP",
"32132","BP+C+BP",
"01301","AAA",
"213","TRAFOC",
"3023","SCALA",
"1231","RES",
//"2130","PHPLP",
"3013013","BBB",
"3102103","CCC",
"202313","PIL+PIC",
"020131","TEL+TEC",
"313202","PIC+PIL",
"131020","TEC+TEL",
"23102310","PBP4",
"01320123","SBP4",
"020303020","TEL+PCLC+TEL",
"020131020","TEL+TEC+TEL",
"2130","PHPLP",   //revisar spice
"213021","PHPLPHP",
"213021303","PHPLPHPLP",
"21303","HPLPPI"
};

GRABIM_Result MatchingNetwork::RunGRABIM()
{
GRABIM_Result Res;
double lambda4 = c0/(4.*mean(f_matching));

if(!topology.compare("-1"))
    {
    //The user did not entered any specific network, so it seems
    //reasonable to try some typical wideband matching network
    printf("!! running GridSearch() with some common networks...\n");
    printf("!! legend  0:L-- 1:C-- 2:L// 3:C// 4:TL-- 5:stub+oc 6:stub+sc\n");
    rowvec Vopt, Vaux;
    string candidate;
    double gridtest, opttopo=99;

    for(int n=0;n!=MAXNETWORKS;n++)
        {
        printf("!! %2d | %-16s %-24s ",n,NETWORKS[n].topology,NETWORKS[n].name);

        topology = NETWORKS[n].topology;
        ////inicializo componentes
        queue <double> XINI;
        for (unsigned int i = 0; i< topology.size();i++)
            {
            if ((!topology.substr(i,1).compare("0"))||(!topology.substr(i,1).compare("2"))) XINI.push(1e-9);
            if ((!topology.substr(i,1).compare("1"))||(!topology.substr(i,1).compare("3"))) XINI.push(1e-12);
            if ((!topology.substr(i,1).compare("4"))||(!topology.substr(i,1).compare("5"))||(!topology.substr(i,1).compare("6")))XINI.push(100),XINI.push(lambda4);
            }
        x_ini = ones(1, XINI.size());
        for (unsigned int i = 0; i < x_ini.size();i++)
            {x_ini.at(i) = XINI.front();XINI.pop();}
        ////
        Vaux = GridSearch();
        gridtest = CandidateEval(Vaux);
        printf(" %6.2f dB",gridtest);
        if (gridtest < opttopo - .5)  //quedarme con la de menos inductores?
            {
            candidate = topology;
            opttopo = gridtest;
            Res.x_grid_search = Vaux;
            printf("     MINS11 found!");
            }
        printf("\n");
        }
    topology=candidate;
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
    cout<<"!! running GridSearch() with network "<< topology << endl;
    Res.x_grid_search = GridSearch();
    }

cout<<"!! running LocalOptimiser() with network "<< topology << endl;

//FASE2
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

std::cout << "!! GRIDSEARCH  S11<"<< Res.grid_val  << "dB\tval:" << Res.x_grid_search;
std::cout << "!! NLOPTSEARCH S11<"<< Res.nlopt_val << "dB\tval:" << Res.x_nlopt;

printf("\n+----SRC-----+");
for(unsigned int i=0;i<topology.size();i++)
    {
                                      cout<<"\n|            |  ";
if(!topology.substr(i,1).compare("0"))cout<<"\n|            L  "<<Res.x_nlopt[i]*1E9 << "nH";
if(!topology.substr(i,1).compare("1"))cout<<"\n|            C  "<<Res.x_nlopt[i]*1E12<< "pF";
if(!topology.substr(i,1).compare("2"))cout<<"\n+-----L------+  "<<Res.x_nlopt[i]*1E9 << "nH";
if(!topology.substr(i,1).compare("3"))cout<<"\n+-----C------+  "<<Res.x_nlopt[i]*1E12<< "pF";
if(!topology.substr(i,1).compare("4"))cout<<"\n|            T  "<<
                                            "\n|            l  "<<Res.x_nlopt[i];
if(!topology.substr(i,1).compare("5"))cout<<"\n|      oc+stub  "<<Res.x_nlopt[i];
if(!topology.substr(i,1).compare("6"))cout<<"\n|      sc+stub  "<<Res.x_nlopt[i];
                                      cout<<"\n|            |  ";
    }
printf("\n+----LOAD----+\n");
return Res;
}



/*

    //Runs GRABIM. Well, this is not exactly the algorithm
    // detailed at [1] but from my point of view is functional and easier to code...
    //Diferences:
    // 1) The candidate vector is not in log units. I do not see a good reason for doing so. Maybe I am missing something important
    // 2) Frequency is not in rad/s.
    // 3) The objective function is the magnitude of S11 expressed in dB. log(x) functions usually have strong
    // gradients so it seem to suggest that this is good for derivative free opt algorithms
    // 4) This code takes advantage from NLopt derivative-free local optimisers. NLopt is easy to link and it works
    // fine. Despite the fact that the Nelder-Mead algorithm does not guarantee convergence (among other problems), it leads to achieve a good local
    // (probably, global) optimum. This is caused by the fact that the matching network should be as simple as possible => few elements => xk \in R^N, where
    // N is typically < 6. Even N=6 is a big number, please consider that matching networks are tight to physical constraints in practice, so, the larger the
    // network, the harder the 'tuning'.

*/