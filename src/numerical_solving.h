#include <RcppArmadillo.h>

#ifndef NUMERICAL_SOLVING_H
#define NUMERICAL_SOLVING_H
#endif
using namespace Rcpp;

void tridiagonalSolving(NumericVector a, NumericVector b, NumericVector c, NumericVector d,
                                 NumericVector e, NumericVector f, NumericVector sol);
