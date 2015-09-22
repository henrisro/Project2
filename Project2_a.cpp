/*
**     Project2: a) and b)
**     A first brute force implementation of Jacobi's rotation algorithm.
**     The method is applied to solving Shrödinger's equation
**     with a harmonic oscillator potenitial for a single electron.
**     Computation time is compared to Armadillo's standard solver and written
**     to txt-file "Table_of_results.txt".
*/

#include <iostream>
#include <cmath>
#include <armadillo>
#include <fstream>
#include <iomanip>
#include <cstring>
#include <time.h>

using namespace std;
using namespace arma;
ofstream ofile;

// Function to find indices of off-diagonal element in (symmetric) matrix A with
// the largest absolute value. Indices k and l are returned.
double off_diag_abs(mat A, int * k, int * l, int n) {
	double max_val = 0.0;
	double A_ij;
    for (int i = 0; i < n; i++) {
    	for (int j = i+1; j < n; j++) {
    		A_ij = fabs(A(i,j));
            // Found new element with larger absolute value:
            if (A_ij > max_val) { max_val = A_ij; *k = i; *l = j; }
    	}
    }
    return max_val;
}


// Function to make one iteration of Jacobi method.
// Rotation is defined by the indices k and l (returned by the function off_diag_abs).
// Matrix A of dimension n is rotated. Matrix R contains eigenvector basis and is also rotated.
// Implementation of this is very close to that found in the lecture notes of FYS4150 (2015).
void Jacobi_rotate(mat & A, mat & R, int k, int l, int n) {
    double t, tau, s, c;

    if ( A(k,l) != 0.0 ) {
      tau = ( A(l,l) - A(k,k) )/( 2*A(k,l) );

      // Ensure that we get the smallest absolute value of rotation angle:
      if (tau >= 0.0)  { t = 1.0/(fabs(tau) + sqrt(1.0 + tau*tau)); }
      else { t = -1.0/(fabs(tau) + sqrt(1.0 + tau*tau)); }
    
      c = 1.0/sqrt(1.0 + t*t);
      s = t*c;
    }
    else {
        // If A(k,l) = 0, no need for rotation:
    	c = 1.0; s = 0.0;
    }

    // Calculate all rotations by hand to avoid multiplication of large matrices:
    double R_ik, R_il, A_ik, A_il;
    double A_kk = A(k,k);
    double A_ll = A(l,l);

    // Rotation of elements defined by k and l:
    A(k,k) = A_kk*c*c - 2*A(k,l)*c*s + A_ll*s*s;
    A(l,l) = A_ll*c*c + 2*A(k,l)*c*s + A_kk*s*s;
    A(k,l) = 0.0;
    A(l,k) = 0.0;

    // Other elements affected by rotation:
    for (int i = 0; i < n; i++) {
    	if (i != k && i != l) {
            A_ik = A(i,k);
            A_il = A(i,l);

    		A(i,k) = A_ik*c - A_il*s;
    		A(k,i) = A(i,k);
    		A(i,l) = A_il*c + A_ik*s;
    		A(l,i) = A(i,l);
    	}
    	R_ik = R(i,k);
        R_il = R(i,l);
        R(i,k) = R_ik*c - R_il*s;
        R(i,l) = R_il*c + R_ik*s;
    }
    return;
}

int main() {
    
    // Initial variables:
    int n, counter;
    int max_counter = 1E6;
    double eps = 1.0E-8;
    double rho_min, rho_max, h, e;
    double calculation_time_Jacobi, calculation_time_arma;
    string outfilename;
    rowvec N;

    // Consider the following choice of rho_max:
    rho_min = 0.0;
    rho_max = 5.0;

    // Comparison will be made for the following system sizes:
    N  << 50 << 100 << 150 << 200 << 250 << 300 << 350;
    
    outfilename = "Table_of_results.txt";
    ofile.open(outfilename);
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << "  rho_max set to: " << rho_max << " and tolerance in Jacobi method set to: " << eps << endl;
    ofile << "   n:         Eigenvalues:           T(Armadillo):     T(Jacobi):   #Transformations:" << endl;
    // Loop over different system sizes:
    for (int j = 0; j < 7; j++) {
    
    n = N(j);
    counter = 0;

    h = (rho_max - rho_min)/n;
    e = -1.0/(h*h);
    
    // Initializing matrices and vectors that will be used:
    vec d(n-1);
    vec rho(n+1);
    vec eigenvalues(n-1);
    mat A(n-1,n-1);
    A.zeros();
    mat R(n-1,n-1);
    // Making R equal to identity initially.
    // Corrsponds to initial choice of basis:
    R.eye();

    for (int i=0; i <= n; i++) { rho(i) = rho_min + i*h; }
    // Now rho(0) = rho_min and rho(n) = rho_max.
    for (int i=0; i < n-1; i++) { d(i) = 2/(h*h) + rho(i+1)*rho(i+1); }
    // Now d(0) is evaluated at rho(1) etc.

    // Initiating matrix A:
    for (int i=0; i < n-1; i++) {
   	  A(i,i) = d(i);
   	  if (i != n-2) { 
   	 	A(i,i+1) = e;
   	  	A(i+1,i) = e;
   	  }
    }
    
    // Make a copy of A that will be used for Armadillo library method:
    mat B = A;
    
    // Calculate computation time and eigenvectors with Jacobi method:
    clock_t start_Jacobi, finish_Jacobi;
    int k, l;
    double diag_maximum = off_diag_abs(A, &k, &l, n-1);
    start_Jacobi = clock();
    while (diag_maximum > eps && counter < max_counter) {
        diag_maximum = off_diag_abs(A, &k, &l, n-1);
        Jacobi_rotate(A, R, k, l, n-1);        
    	counter++;
    }
    finish_Jacobi = clock();
    calculation_time_Jacobi = (finish_Jacobi - start_Jacobi)/(double)CLOCKS_PER_SEC;
    // Extract diagonal elements (eigenvalues) and sort them:
    for (int i = 0; i < n-1; i++) { eigenvalues(i) = A(i,i); }
    eigenvalues = sort(eigenvalues);


    // Calculate computation time and eigenvectors with standard Armadillo library:
    clock_t start_arma, finish_arma;
    start_arma = clock();
    mat eigvec;
    vec eigval;
    eig_sym(eigval, eigvec, B);
    finish_arma = clock();
    calculation_time_arma = (finish_arma - start_arma)/(double)CLOCKS_PER_SEC;

    // Write results to outfile:
    ofile << setw(5) << n << "    ";
    ofile << "{" << setprecision(6) << eigenvalues(0) << ", ";
    ofile << setprecision(6) << eigenvalues(1) << ", ";
    ofile << setprecision(6) << eigenvalues(2) << "}";
    ofile << setw(12) << setprecision(6) << calculation_time_arma << " s.";
    ofile << setw(12) << setprecision(6) << calculation_time_Jacobi << " s.";
    ofile << setw(10) << counter << endl;

    }
    ofile.close();
    return 0;
}
