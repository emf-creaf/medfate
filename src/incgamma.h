#include <Rcpp.h>

#ifndef INCGAMMA_H
#define INCGAMMA_H
#endif
using namespace Rcpp;

NumericVector incgam(double a, double x);
double invincgam(double a, double p, double q);