# Program for plotting the wave functions found 
# with Armadillo library and comparing them to
# exact results.

import os

os.system('g++ Project2_a_plot_eigenfunctions.cpp -o Project2_a_plot_eigenfunctions.o -O2 -I /Users/Henrik/FYS4150\ -\ Computational\ physics/armadillo-5.500.2/include -lblas -llapack')
os.system('./Project2_a_plot_eigenfunctions.o')

from math import *
import numpy as np
import matplotlib.pyplot as plt

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

# Exact results are defined in some functions:
def analytical_ground_state1(rho):
	return 4/sqrt(pi)*rho**2*exp(-rho**2)

def analytical_ground_state2(rho):
	return 8/(3*sqrt(pi))*rho**2*(1.5-rho**2)**2*exp(-rho**2)

def analytical_ground_state3(rho):
	return 8.0/(15*sqrt(pi))*rho**2*(15.0/4-5*rho**2+rho**4)**2*exp(-rho**2)

# Fetching data by a call on read_x_u_v for three different n:
rho, u1, u2, u3 = read_rho_wavefunctions('Wavefunction.txt')
analytical1 = []; analytical2 = []; analytical3 = [];
for i in range(len(rho)):
	rho_val = rho[i]
	analytical1.append(analytical_ground_state1(rho_val))
	analytical2.append(analytical_ground_state2(rho_val))
	analytical3.append(analytical_ground_state3(rho_val))
	#print "Analytical1: ", analytical_ground_state1(rho_val), " Numerical1: ", u1[i]

# Plotting commands:
plt.rcParams.update({'font.size': 14})
fig, ax = plt.subplots(1)
ax.plot(rho,u1,'r-',label='$\mid u_0(\\rho) \mid^2$, numerical')
ax.plot(rho,u2,'b-',label='$\mid u_1(\\rho) \mid^2$, numerical')
ax.plot(rho,u3,'g-',label='$\mid u_2(\\rho) \mid^2$, numerical')
ax.plot(rho,analytical1,'r+',label='$\mid u_0(\\rho) \mid^2$, analytical')
ax.plot(rho,analytical2,'b+',label='$\mid u_1(\\rho) \mid^2$, analytical')
ax.plot(rho,analytical3,'g+',label='$\mid u_2(\\rho) \mid^2$, analytical')
ax.set_xlabel('Radial coordinate $\\rho$')
ax.set_ylabel('Radial wavefunction: $\mid u(\\rho) \mid^2$')
ax.legend(loc='upper right',fancybox='True')
ax.set_title('Radial wavefunction.')
ax.grid()
plt.show()