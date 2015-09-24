# Program to compile and run Project2_d.cpp

import os

os.system('g++ Project2_d.cpp -o Project2_d.o -O2 -I /Users/Henrik/FYS4150\ -\ Computational\ physics/armadillo-5.500.2/include -lblas -llapack')
os.system('./Project2_d.o')

from math import *
import numpy as np
import matplotlib.pyplot as plt
omega_r = 5.0

def read_rho_wavefunctions(filename):
    infile = open(filename, 'r')
    # Elements to be read in file:
    rho = []; u1 = []; u2 = []; u3 = [];
    # Read lines except for the first one:
    lines = infile.readlines()[1:]
    for line in lines:
        words = line.split()
        rho.append(float(words[0]))
        u1.append(float(words[1]))
        u2.append(float(words[2]))
        u3.append(float(words[3]))
    infile.close()
    return rho, u1, u2, u3

# Fetching data by a call on read_x_u_v for three different n:
rho, u1, u2, u3 = read_rho_wavefunctions('Wavefunction_relative.txt')

# Plotting commands to look at the wave functions:
plt.rcParams.update({'font.size': 14})
fig, ax = plt.subplots(1)
ax.plot(rho,u1,'r-',label='$\mid u_0(\\rho) \mid^2$')
ax.plot(rho,u2,'b-',label='$\mid u_1(\\rho) \mid^2$')
ax.plot(rho,u3,'g-',label='$\mid u_2(\\rho) \mid^2$')
ax.set_xlabel('Radial coordinate $\\rho$')
ax.set_ylabel('Relative radial wavefunction: $\mid u(\\rho) \mid^2$')
ax.legend(loc='upper right',fancybox='True')
ax.set_title('Relative radial wave function with $\omega_r =$ %.2f.' % omega_r)
ax.grid()
plt.show()