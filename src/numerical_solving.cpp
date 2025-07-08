// [[Rcpp::interfaces(r,cpp)]]

#include <Rcpp.h>
using namespace Rcpp;

// a,b,c,d are input vectors
// e,f are vectors used internally
// sol is the output
void tridiagonalSolving(NumericVector a, NumericVector b, NumericVector c, NumericVector d,
                        NumericVector e, NumericVector f, NumericVector sol) {
  int n = b.size();
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
  sol[n-1] = f[n-1];
  for(int i = (n - 2);i>=0;i--) {
    sol[i] = f[i] - e[i]*sol[i + 1];
  }
}

