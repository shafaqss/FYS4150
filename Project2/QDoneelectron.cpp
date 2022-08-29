//Quantum dots in 3-dimensions, one electron
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <armadillo>
#include <vector>
#include <iomanip>

using namespace std;
using namespace arma;

mat create_matrix(double &rho_max, int n){
  /*Creating the matrix A*/
  double h = rho_max/(n+1); // step size
  mat A = mat(n,n); A.fill(0.0);
  double d = 2/(h*h);  //diagonal element
  double e = -1/(h*h); // nondiagonal element

  for (int i=0; i<n ; i++){ //filling the matrix
    double rho = (i+1)*h;
    A(i,i) = d + (rho*rho);
    if (i<n-1){
      A(i,i+1)=e;
      A(i+1,i)=e;
    }
  }
  return A;
}

double max_offdiagonal(mat &A, int &k, int &l, int n){
  /*Returns the non-diagonal element with the maximum absolute value
  assuming that the input matrix A is symmetric*/
  double current_max = 0;
  for (int i=0; i<n; i++){
    for (int j=i+1; j<n; j++){
      double aij = fabs(A(i,j));
      if (aij > current_max){
        current_max = aij;
        k = i;  l = j;
      }
    }
  }
  return current_max;
}

void jacobi_rotation(mat &A, mat &R, int &k, int &l, int n){
  /*Performs the jacobi rotation in  in the kth column and lth row*/
  double tau, t, c, s;
  tau = (A(l,l) - A(k,k)) / (2*A(k,l));
  if (tau >= 0) {
      t = 1.0/(tau + sqrt(1 + tau*tau));
  } else {
      t = -1.0/(-tau + sqrt(1 + tau*tau));
  }
  c = 1.0/(sqrt(1 + (t*t))); //cosine element
  s = c*t;                   //sine element
  if (A(k,l) == 0) {
      c = 1.0; //cosine
      s = 0.0; //sine
  }
  // defning the matrix elements
  double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
  // changing the matrix elements with indices k and l
  a_kk = A(k, k);
  a_ll = A(l, l);
  A(k ,k) = a_kk*c*c - 2*A(k, l)*c*s + a_ll*s*s;
  A(l, l) = a_ll*c*c + 2*A(k, l)*c*s + a_kk*s*s;
  A(l, k) = 0.0;
  A(k, l) = 0.0;
  // changing the remaining elemenst
  for (int i=0; i<n; i++){
    if ((i != k) && (i != l)) {
        a_ik = A(i, k);
        a_il = A(i, l);
        A(i, k) = a_ik*c - a_il*s;
        A(k, i) = A(i, k);
        A(i, l) = a_il*c + a_ik*s;
        A(l, i) = A(i, l);
    }

    // now we compute the new eigenvectors
    r_ik = R(i, k);
    r_il = R(i, l);
    R(i, k) = c*r_ik - s*r_il;
    R(i, l) = c*r_il + s*r_ik;
  }
}

void jacobi_method(mat &A, mat &R, int n){
  /*performing the jacobi algorithm to find the eigenvalues*/
  double eps = 1.0E-8;
  int maximum_iterations = n*n*n; //setting an upper bound for the iterations
  int k = 0;
  int l = 0;
  int iterations = 0; // counter variable
  double max_nondiag = max_offdiagonal(A, k, l, n);

  while ((max_nondiag > eps) && (iterations <= maximum_iterations)){
    jacobi_rotation(A, R, k, l, n);
    max_nondiag = max_offdiagonal(A, k, l, n);
    iterations ++;

  }
  cout << "number of iterations are = " << iterations<<endl;
}

double eigenvalue_analytic(double rho_max, int n,int j){
  /*return the analytical eigenvalues for the matrix*/
  //j is the number of the eigenvalue
  double pi = acos(-1.0);
  double h = rho_max/(n+1);
  double dd = 2/(h*h); //diagonal element
  double aa = -1/(h*h); //nondiagonal element
  double exact = dd + (2*aa*cos(((j+1)*pi)/(n+1))); //eigenvalue calculation
  return exact;
}

//TEST FUNCTIONS
void test_eigenvalues(){
  /*test to see if the jacobi method computes the correct
  eigenvalues, testing using a simple 3x3 matrix*/
  double tol = 1.0E-8;
  int n = 3;
  mat T;
  mat R; R.eye(n,n);
  T   << 6  << -2 << -1 << endr
      << -2 << 6 << -1 << endr
      << -1 <<-1 << 5 << endr;
  vec eigval(n);
  eig_sym(eigval, T);
  jacobi_method(T,R,n);
  vec computed_eigval = T.diag();
  computed_eigval = sort(computed_eigval, "ascend");

  for (int i=0; i<n; i++){
    if (abs(computed_eigval(i)-eigval(i)) > tol){
      cout<<"Error: The computed eigenvalues are incorrect!"<<endl;
    }}
}

void test_max_offdiagonal(){
  /*testing the max_offdiagonal() function to make sure the maximum
  offdiagonal elemnet is found. Note however that this function
  assumes the input matrix is symmetric*/
  double tol = 1.0E-8;
  int n = 4;
  int k = 0;
  int l = 0;
  mat T;
  T   << 6  << -2 << -1<< 21<<endr
      << -2 << 6 <<  8 << 3 <<endr
      << -1 << 8 << -1 << 5 <<endr
      << 21 << 3 << 5  << 8 <<endr;
  double comp_max = max_offdiagonal(T,k,l,n);
  if (abs(comp_max-21)>tol){
    cout << "Error:The maximum diagonal function has failed!"<<endl;
  }
}

int main(int argc, const char **argv) {
  //test_eigenvalues();
  //test_max_offdiagonal();
  int n = 5; double rho_max = 8;
  cout << "Matrix of size, n = " << n << endl;
  mat A,R; R.eye(n,n); // create the matrices
  A = create_matrix(rho_max,n);

  vec eigvals(n);// finding eigenvlaues usig buit in armadillo function
  eig_sym(eigvals,A);
  for (int i=0; i<n; i++){
    cout <<"The eigenvalue is "<< eigvals(i)<<endl;
  }
  jacobi_method(A,R,n); // eigenvalues using jacobis method


return 0;
}
