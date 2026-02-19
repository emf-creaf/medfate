#include <RcppArmadillo.h>
using namespace Rcpp;

#ifndef GROWTH_DAY_H
#define GROWTH_DAY_H

void growthDay_private(List internalCommunication, List x, NumericVector meteovec, 
                       double latitude, double elevation, double slope, double aspect,
                       double solarConstant, double delta, 
                       double runon = 0.0, Nullable<NumericVector> lateralFlows = R_NilValue, double waterTableDepth = NA_REAL, 
                       bool verbose = false);

#endif
