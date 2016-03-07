# GRABIM-CPP
Grid search algorithm for wideband matching network synthesis + Local optimiser

Reference:
[1] Broadband direct-coupled and RF matching networks. Thomas R. Cuthbert, 1999

Requirements:

* Armadillo: http://arma.sourceforge.net/download.html
* NLopt: http://ab-initio.mit.edu/wiki/index.php/NLopt
* GNUplot: http://gnuplot.info/download.html

How to build:

* From QtCreator: just open GRABIM-CPP.pro
* From terminal: g++ main.cpp matchingnetwork.cpp io.cpp -I(NLoptLibPath)include -L(NLoptLibPath)lib -lnlopt -lm -std=c++11

where (NLoptLibPath) is the directory where NLopt is installed (by default, /usr/local/)


How to run:

./GRABIM-CPP (path-to-s2p-file-source) (path-to-s2p-file-load) f1 f2 topology opt_algorithm
./GRABIM-CPP (source-impedance) (load-impedances) f1 f2 topology opt_algorithm


where:

* f1 and f2 are the lowest and the highest frequency where matching is required.
* (source-impedance) and (load-impedances) represent constant and complex impedances (i.e. the must be something like 50+j10, 75+j0,... )
* topology defines the kind of matching network. The user can build any ladder circuit using the following coding:
  - 0: Series inductance
  - 1: Series capacitor
  - 2: Parallel inductance
  - 3: Parallel capacitor
  - 4: Transmission line
  - 5: Open stub
  - 6: Short-circuited stub 

By default, if topology = "-1", the program tries a few well known wideband (bandpass and lowpass) structures.

Examples:

* ./GRABIM-CPP 75+j10 13+j30 2e9 4e9 -1 NLOPT_LN_NELDERMEAD 1 test1_Nelder
* ./GRABIM-CPP 10+j10 20-j30 1e9 4e9 -1 NLOPT_LN_NELDERMEAD 1 test2_Nelder
* ./GRABIM-CPP 10-j10 50+j0 5e9 7e9 -1 NLOPT_LN_SBPLX 1 test1_Subplex





