#include <Rcpp.h>

#ifndef NUMERICAL_SOLVING_H
#define NUMERICAL_SOLVING_H
#endif
using namespace Rcpp;

NumericVector tridiagonalSolving(double* a, double* b, double* c, double* d);