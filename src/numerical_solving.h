#include <Rcpp.h>

#ifndef NUMERICAL_SOLVING_H
#define NUMERICAL_SOLVING_H
#endif
using namespace Rcpp;

NumericVector tridiagonalSolving(NumericVector a, NumericVector b, NumericVector c, NumericVector d);