# Program for running Project2_f_Lanczos.cpp

import os

os.system('g++ Project2_f_Lanczos.cpp -o Project2_f_Lanczos.o -O2 -I /Users/Henrik/FYS4150\ -\ Computational\ physics/armadillo-5.500.2/include -lblas -llapack')
os.system('./Project2_f_Lanczos.o')