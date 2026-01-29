// [[Rcpp::interfaces(r,cpp)]]
#include <RcppArmadillo.h>
#include "modelInput_c.h"
#include "communication_structures_c.h"
using namespace Rcpp;

Rcpp::NumericMatrix copyNumericMatrix_c(arma::mat comm, int rows, int cols) {
  NumericMatrix out(rows, cols);
  for(int r=0;r<rows;r++) {
    for(int c=0;c<cols; c++) {
      out(r, c) = comm(r, c);
    }
  }
  return(out);
}
