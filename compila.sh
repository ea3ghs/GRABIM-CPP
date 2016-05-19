clear

rm *~;ls -l *.c* *.h*;rm a.out;
g++ -g -o a.out *.cpp -lnlopt -lm -larmadillo -std=c++11 
#clang++-3.8  -o a.out *.cpp -lnlopt -lm -larmadillo -std=c++11 

ls -l a.out
rm *.png temporal*

./a.out 75+j10 13+j30 2e9 4e9 -1 NLOPT_LN_NELDERMEAD  1 1pngLN_NELDERMEAD 
#gdb --args ./a.out 75+j10 13+j30 2e9 4e9 -1 NLOPT_LN_NELDERMEAD 1 xx
