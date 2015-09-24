# Program for running Jacobi's rotation algorithm

import os

os.system('g++ Project2_c_eigenvalues.cpp -o Project2_c_eigenvalues.o -O2 -I /Users/Henrik/FYS4150\ -\ Computational\ physics/armadillo-5.500.2/include -lblas -llapack')
os.system('./Project2_c_eigenvalues.o')