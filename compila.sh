clear

rm *~;ls -l *.c* *.h*;rm a.out;
g++ -g -o a.out GRABIM.cpp main.cpp sparengine.cpp io.cpp -lnlopt -lm -larmadillo -std=c++11
#clang++-3.8

./a.out 6L6.s1p 50 24900e3 28500e3 1212 ;gnuplot plotscript
./a.out 6L6.s1p 50 14000e3 28500e3 1212 ;gnuplot plotscript

#./a.out 75+j10 13+j30 2e9 4e9 123456
#./a.out 75+j10 13+j30 2e9 4e9  
#gdb --args ./a.out 75+j10 13+j30 2e9 4e9 
