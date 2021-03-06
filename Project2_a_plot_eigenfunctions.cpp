/*
**     Project2: a)
**     A program to plot the normalized eigenfunctions.
**     They are compared to exact wave functions.
*/

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <armadillo>
#include <cstring>
#include <time.h>

using namespace std;
using namespace arma;
ofstream ofile;

int main() {
    
    // Initial variables: 
    int n = 350;
    double rho_min, rho_max, h, e;
    double calculation_time_arma;
    string outfilename;

    rho_min = 0.0;
    rho_max = 5.0;
    outfilename = "Wavefunction.txt";

    h = (rho_max - rho_min)/n;
    e = -1.0/(h*h);
    
    vec d(n-1);
    vec rho(n+1);
    mat A(n-1,n-1);

    for (int i=0; i <= n; i++) { rho(i) = rho_min + i*h; }
    // Now rho(0) = rho_min and rho(n) = rho_max.
    for (int i=0; i < n-1; i++) { d(i) = 2/(h*h) + rho(i+1)*rho(i+1); }
    // Now d(0) is evaluated at rho(1) etc.

    // Initializing matrix A:
    for (int i=0; i < n-1; i++) {
   	  A(i,i) = d(i);
   	  if (i != n-2) { 
   	 	A(i,i+1) = e;
   	  	A(i+1,i) = e;
   	  }
    }

    // Calculate computation time and print out eigenvalues to check:
    clock_t start_arma, finish_arma;
    start_arma = clock();
    mat eigvec;
    vec eigval;
    eig_sym(eigval, eigvec, A);
    finish_arma = clock();
    calculation_time_arma = (finish_arma - start_arma)/(double)CLOCKS_PER_SEC;

    cout << endl;
    cout << "Standard armadillo solver. Time: " << calculation_time_arma << " s." << endl;
    for (int i= 0; i < 3; i++) {
    	cout << eigval(i) << endl;
    }

    // Eigenfunctions are normed by the use of a simple trapezoidal integration:
    double norm1, norm2, norm3;
    vec V1 = eigvec.col(0);
    vec V2 = eigvec.col(1);
    vec V3 = eigvec.col(2);

    vec U_square_normed1(n-1);
    vec U_square_normed2(n-1);
    vec U_square_normed3(n-1);

    for (int i = 0; i < n-1; i++) {
        double V1_i2 = V1(i)*V1(i);
        double V2_i2 = V2(i)*V2(i);
        double V3_i2 = V3(i)*V3(i);

        norm1 += V1_i2; norm2 += V2_i2; norm3 += V3_i2;
        U_square_normed1(i) = V1_i2;
        U_square_normed2(i) = V2_i2;
        U_square_normed3(i) = V3_i2;
    }
    norm1 *= h; norm2 *= h; norm3 *= h;
    U_square_normed1 /= norm1;
    U_square_normed2 /= norm2;
    U_square_normed3 /= norm3;

    // Normed results are written to file and read by a python script:
    ofile.open(outfilename);
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << "       rho:         u1             u2             u3:" << endl;
    // Loop over all n producing table of time used:
    for (int i=0; i<n-1; i++) {
      double U_val1 = U_square_normed1(i);
      double U_val2 = U_square_normed2(i);
      double U_val3 = U_square_normed3(i);
      ofile << setw(15) << setprecision(8) << rho(i+1);
      ofile << setw(15) << setprecision(8) << U_val1;
      ofile << setw(15) << setprecision(8) << U_val2;
      ofile << setw(15) << setprecision(8) << U_val3 << endl;
    }
    ofile.close();

    return 0;
}
