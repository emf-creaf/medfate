#include <Rcpp.h>

#ifndef DECOMPOSITION_H
#define DECOMPOSITION_H
#endif
using namespace Rcpp;

double litterMetabolicFraction(double ligninPercent, double Nmass);
double annualLitterDecompositionRate(double AET, double lignin);