// [[Rcpp::interfaces(r,cpp)]]

#include <Rcpp.h>
using namespace Rcpp;

double* tridiagonalSolving(double* a, double* b, double* c, double* d, int n) {
  double* e = new double[n];
  double* f = new double[n];
  double* u = new double[n];
  
  //Forward steps
  double e_prev = 0.0;
  double f_prev = 0.0;
  for(int i=0;i<n;i++) {
    e[i] = c[i]/(b[i] - a[i]*e_prev);
    f[i] = (d[i] - a[i]*f_prev)/(b[i] - a[i]*e_prev);
    // Rcout<<i<< " "<< e[i]<< " "<< f[i]<<"\n";
    e_prev = e[i];
    f_prev = f[i];
  }
  //Backward steps
  u[n-1] = f[n-1];
  for(int i = (n - 2);i>=0;i--) {
    u[i] = f[i] - e[i]*u[i + 1];
  }  
  delete[] e;
  delete[] f;
  return(u);
}

