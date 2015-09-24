/*
**     Project2: c)
**     A program to print out the
**     eigenvalues for comparison.
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

// Approximate formula derived in reference article:
double article_eigenvalue_approximate(int m, double omega_r) {
    double epsilon = 3.0/2.0*pow(omega_r/2.0, 2.0/3.0) + sqrt(3.0)*omega_r*(m+0.5);
    return 2.0*epsilon;
}

int main() {
    
    // Initial variables: ///////
    // rho_max and n are meant to be adjusted.
    int n = 350;
    double rho_min, rho_max, h, e;
    double calculation_time_arma;
    double omega_r = 5.0;
    string outfilename;

    rho_min = 0.0;
    rho_max = 2.0;
    /////////////////////////////

    h = (rho_max - rho_min)/n;
    e = -1.0/(h*h);
    
    vec d(n+1);
    vec rho(n+1);
    mat A(n-1,n-1);

    for (int i=0; i <= n; i++) { rho(i) = rho_min + i*h; }
    // Now rho(0) = rho_min and rho(n) = rho_max.
    for (int i=0; i < n-1; i++) { d(i) = 2.0/(h*h) + omega_r*omega_r*rho(i+1)*rho(i+1) + 1.0/rho(i+1); }
    // Now d(0) is evaluated at rho(1) etc.

    // Initiating matrix A:
    for (int i=0; i < n-1; i++) {
      A(i,i) = d(i);
      if (i != n-2) { 
        A(i,i+1) = e;
        A(i+1,i) = e;
      }
    }

    mat eigvec;
    vec eigval;
    eig_sym(eigval, eigvec, A);
    double temp_eigenval;

    // Print out information to screen:
    cout << endl;
    cout << "Ground state information with oscillator frequency: " << omega_r << endl;
    temp_eigenval = article_eigenvalue_approximate(0,omega_r);
    cout << "Eigenvalue (groundtate) found numerically: " << eigval(0) << endl;
    cout << "Article formula: " << temp_eigenval << endl;
    cout << "Relative error: " << (eigval(0)-temp_eigenval)/temp_eigenval << endl;
    cout << endl;

    return 0;
}
