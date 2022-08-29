//c++ -o 1b.exe task1b.cpp -std=c++11

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <iomanip>
using namespace std;

double f(double x, double h) {
        return h*h*100*exp(-10*x);
}

double exact(double x) {
    return 1.0-(1-exp(-10))*x-exp(-10*x);
}

void gen_alg(int n, double h, double *b, double *a, double *c, double *b_tilde, double *g_tilde, double *g, string name) {
    double *v = new double[n+1];

    //Forward substitution loop     
    for(int i = 1; i < n+1; i++) {
        b_tilde[i] = b[i] - (a[i] *c[i-1])/b_tilde[i-1];
        g_tilde[i] = g[i] - (a[i]/b_tilde[i-1] * g_tilde[i-1]);
    } 

    //Backwards substitution loop
    v[n] = g_tilde[n]/b[n];
    for(int i = n-1; i >= 1; i--) {
        v[i] = (g_tilde[i] - c[i] * v[i+1]/b_tilde[i]);
    }   
    
    ofstream fileout;
    string filename = name + ".txt";
    fileout.open(filename);
    fileout << n << endl;
    
    for(int i = 0; i < n; i++) {
        fileout << i*h << ", " << v[i] << endl;    
    }  
    fileout.close();  
}    
    
int main(int argc, char* argv[]) {    
    string outputfile;
    int n;

    if(argc == 1) {
          cout << "Error: Filename not found. " << endl;
          exit(1);
    
    } else if(argc < 3) {
        cout << "Error: n missing." << endl;
        exit(1);
           
    } else {
        outputfile = argv[1];
        n = atoi(argv[2]);
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

    gen_alg(n, h, b, a, c, b_tilde, g_tilde, g, outputfile);

    delete[] a; delete[] b; delete[] c; delete[] b_tilde; 
    delete[] g_tilde, delete[] g, delete[] x;
    
    return 0;
}





