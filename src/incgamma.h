#include <Rcpp.h>

#ifndef INCGAMMA_H
#define INCGAMMA_H
#endif
using namespace Rcpp;

double* incgam(double a, double x);
double invincgam(double a, double p, double q);
double errorfunction(double x, bool erfcc, bool expo);