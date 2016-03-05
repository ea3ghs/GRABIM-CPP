#include "io.h"

IO::IO()
{
}

//Export data to display results with GNUplot
int IO::exportGNUplot(GRABIM_Result Res, string filepath, int plot)
{
    ofstream GNUplotExport;
    GNUplotExport.open (filepath, std::ofstream::out);
    if(!GNUplotExport.is_open()) return -1;
    GNUplotExport << "#GRABIM data" << endl;
    GNUplotExport << "#freq S11(Grid search) S11(Local optimiser) S21(Grid search) S21(Local optimiser) " << endl;
    for (int i = 0; i < Res.f_analysis.n_rows; i++)
    {
        GNUplotExport << Res.f_analysis.at(i)*1e-9 << " " << 20*log10(abs(Res.S11_gridsearch.at(i))) << " "
                      << 20*log10(abs(Res.S21_gridsearch.at(i))) << " " << 20*log10(abs(Res.S11_nlopt.at(i)))
                      << " " << 20*log10(abs(Res.S21_nlopt.at(i))) << endl;
    }
    GNUplotExport.close();

    if (plot)//Display the result using gnuplot
    {
      string command;
      command = "gnuplot plotscript";
      system(command.c_str());
    }
    return 0;
}
