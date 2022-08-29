// c++ -o 1e.exe task1e.cpp -larmadillo -llapack -lblas -std=c++11

#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>
#include <iomanip>
#include <armadillo>

using namespace std;
using namespace arma;

double f(double x, double h) {
        return h*h*100*exp(-10*x);
}

double exact(double x) {
    return 1.0-(1-exp(-10))*x-exp(-10*x);
}


void LU_decomp(int n, double *a, double *b, double *c, double *g) {
    
    mat A = zeros<mat>(n, n);
    vec b_vector = zeros<vec>(n);
    
    for(int i = 0; i < n; i++) {
        b_vector(i) = g[i];       
        for(int j = 0; j < n; j++) {
            if(i == j) {
                A(i, j) = b[i];
            }else if(j == i+1) {
                A(i, j) = a[i];
            }else if(j == i-1) {
                A(i,j) = c[i];
            }    
        }                        
    }

    auto start = chrono::high_resolution_clock::now();        
    vec x = solve(A, b_vector);        
    auto finish = chrono::high_resolution_clock::now();
    cout << "Time used to compute the LU-decomposition: " << chrono::duration_cast<chrono::nanoseconds>(finish-start).count()/pow(10,9) << " seconds" << endl;
}
    
int main(int argc, char* argv[]) {    
    int n;

    if(argc == 0) {
          cout << "Error: n missing. " << endl;
          exit(1);
           
    } else {

        n = atoi(argv[1]);
   }
               
    double h = 1.0/(n);
    double *g = new double[n+1];
    double *x = new double[n+1];

    double *a = new double[n+1];
    double *b = new double[n+1];
    double *c = new double[n+1];
    
    double *b_tilde = new double[n+1];
    double *g_tilde = new double[n+1];
      
    for(int i = 0; i < n+1; i++) {       
        a[i] = -1.0;
        b[i] = 2.0;
        c[i] = -1.0;
        x[i] = i*h;
        g[i] = f(x[i], h);
    }
    
    b_tilde[0] = b[0];
    g_tilde[0] = g[0];
    
    LU_decomp(n, a, b, c, g);

    delete[] a; delete[] b; delete[] c; delete[] b_tilde; 
    delete[] g_tilde, delete[] g, delete[] x;
    
    return 0;
}
