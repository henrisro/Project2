/*
**     Project2: f)
**     A first implementation of Lanczos' algorithm.
**     The function Lanczos() takes in matrix A and 
**     transforms it to a tridiagonal matrix which is
**     returned. If one puts q1(0) = 1 (and the rest equal
**     to zero), this produces basically the same matrix as
**     with Householder's algorithm.
*/

#include <iostream>
#include <cmath>
#include <armadillo>
#include <fstream>
#include <iomanip>
#include <cstring>

using namespace std;
using namespace arma;
ofstream ofile;

// Computes and returns the 2-norm of a vector v:
double norm_by_hand(vec v, int n) {
    double sum_norm = 0.0;
    for (int i=0; i < n; i++) {
        sum_norm += v(i)*v(i);
    }
    sum_norm = sqrt(sum_norm);
    return sum_norm;
}

// Transform matrix to tridiagonal form with Lanczos'
// iterative algorithm and later extract its eigenvalues by
// some other method.
void Lanczos(mat & A, int n) {
    // tolerance in beta values:
    double btol = 1E-5;
    // initiate vectors and matrices used in the loop:
    int kmax = n;
    int k = 0;
    vec alpha = zeros(kmax);
    vec beta = zeros(kmax);
    vec q_k = zeros(n);
    vec q1 = zeros(n);
    //q1(0) = 1;
    // Generate random normed vector:
    q1 = randu<vec>(n);
    q1 /= norm_by_hand(q1,n);
    vec r = q1;
    double b = 1;

    // Loop over k:
    while (b > btol && k < kmax) {
        k += 1;
        vec qkm1 = q_k;
        q_k = r/b;
        vec Aqk = A*q_k;
        // Update alpha element:
        alpha(k-1) = dot(trans(q_k),Aqk);
        r = Aqk - q_k*alpha(k-1) - qkm1*b;
        b = norm_by_hand(r,n);
        // Update beta element:
        beta(k-1) = b;
    }

    // Overwrite matrix A with the newly computed tridiagonal one:
    A = zeros(n,n);
    for (int i=0; i < n; i++) {
        A(i,i) = alpha(i);
        if (i != n-1) { 
            A(i+1,i) = beta(i);
            A(i,i+1) = beta(i);
        }
    }
    return;
}

int main() {
    // Initial variables:
    int n = 50;
    string outfilename;
    mat A(n,n, fill::randu);
    // Produce an initially symmetric matrix:
    mat B = trans(A);
    A = 0.5*(B+A);

    // Calculate eigenvalues with Armadillo library directly on A:
    mat eigvec0;
    vec eigval0;
    eig_sym(eigval0, eigvec0, A);

    // Make A tridiagonal to check if the matrix has
    // eigenvalues close to that of the original matrix:
    Lanczos(A,n);
    mat eigvec;
    vec eigval;
    eig_sym(eigval, eigvec, A);

    outfilename = "Table_of_results_Lanczos.txt";
    ofile.open(outfilename);
    ofile << setiosflags(ios::showpoint | ios::uppercase);

    // Write results to txt file:
    ofile << "Eigenvalues found with Armadillo: " << endl;
    ofile << "On matrix A:    On matrix T:   Relative difference:" << endl;
    for (int i=0; i < n; i++) {
        ofile << setprecision(6) << setw(10) << eigval0(i) << "       ";
        ofile << setprecision(6) << setw(10) << eigval(i) << "         ";
        ofile << setw(10) << fabs((eigval0(i) - eigval(i))/eigval0(i)) << endl;
    }
    ofile.close();
    return 0;
}
