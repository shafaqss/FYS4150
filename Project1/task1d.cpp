//c++ -o 1c.exe task1c.cpp -std=c++11
/*This program finds the relative error, to run the program
please include the value of n on the command line, example
of execution

        ./task1d.out 100
*/
#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>
#include <iomanip>
using namespace std;

double f(double x, double h) {
        return h*h*100*exp(-10*x);
}

double exact(double x) {
    return 1.0-(1-exp(-10))*x-exp(-10*x);
}

void spec_algo(int n, double h, double *g_tilde, double *g, double *analytic) {
    double *v = new double[n+1];
    double *b_tilde_s = new double[n+1];

    for(int i = 0; i <= n; i++) {
        b_tilde_s[i] = (i+2.0)/(i+1.0);
    }

    //Forward substitution loop
    for(int i = 1; i < n+1; i++) {
        g_tilde[i] = g[i] + g_tilde[i-1]/b_tilde_s[i-1];
    }
    v[n] = (g_tilde[n])/2.0;

    //Backwards substitution loop
    for(int i = n-1; i >= 1; i--) {
        v[i] = (g_tilde[i] + v[i+1])/(b_tilde_s[i]);
    }

    //relative Error
    double error[n+1];
    for (int i=1;i<n+1; i++) {
      error[i] = log10(abs((v[i]-analytic[i])/analytic[i]));
    }
    //finding the max value of the relative error
    double max_value = error[0];
    for (int i=0; i<n; i++){
      if (abs(error[i])>abs(max_value)){
        max_value = error[i];
      }
    }

    cout << "n = "<< n<<endl;
    cout<< "log10(h): "<<log10(h)<<endl;
    cout<<"The maximum relative errror is: "<< max_value<<endl;
}

int main(int argc, char* argv[]) {
    int n;

    if(argc == 1) {
          cout << "Error: n missing. " << endl;
          exit(1);

    } else {
        n = atoi(argv[1]);
   }

    double h = 1.0/(n);
    double *g = new double[n+1];
    double *x = new double[n+1];
    double *u_analytical = new double[n+1];

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
        u_analytical[i] = exact(x[i]);
    }

    b_tilde[0] = b[0];
    g_tilde[0] = g[0];

    spec_algo(n, h, g_tilde, g, u_analytical);

    delete[] a; delete[] b; delete[] c; delete[] b_tilde;
    delete[] g_tilde, delete[] g, delete[] x; delete[] u_analytical;

    return 0;
}
