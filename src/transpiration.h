#include <Rcpp.h>

#ifndef TRANSPIRATION_H
#define TRANSPIRATION_H
#endif
using namespace Rcpp;

List profitMaximization(List supplyFunction, DataFrame photosynthesisFunction, int type, double Gwmin, double Gmax, double kleafmax = NA_REAL);
