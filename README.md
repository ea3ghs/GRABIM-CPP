# GRABIM-CPP
Grid search algorithm for wideband matching network synthesis + Local optimiser

Requirements:

Armadillo: http://arma.sourceforge.net/download.html
NLopt: http://ab-initio.mit.edu/wiki/index.php/NLopt
GNUplot: http://gnuplot.info/download.html

How to build:

* From QtCreator: just open GRABIM-CPP.pro
* From terminal: g++ main.cpp matchingnetwork.cpp io.cpp -I(NLoptLibPath)include -L(NLoptLibPath)lib -lnlopt -lm -std=c++11

where (NLoptLibPath) is the directory where NLopt is installed (by default, /usr/local/)
