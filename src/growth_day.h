#include <RcppArmadillo.h>

#ifndef GROWTH_DAY_H
#define GROWTH_DAY_H
#endif
using namespace Rcpp;

List growthDay(List x, CharacterVector date, NumericVector meteovec, 
               double latitude, double elevation, double slope = NA_REAL, double aspect = NA_REAL,  
               double runon = 0.0, Nullable<NumericVector> lateralFlows = R_NilValue, double waterTableDepth = NA_REAL, 
               bool modifyInput = true);

void growthDay_private(List internalCommunication, List x, NumericVector meteovec, 
                       double latitude, double elevation, double slope, double aspect,
                       double solarConstant, double delta, 
                       double runon = 0.0, Nullable<NumericVector> lateralFlows = R_NilValue, double waterTableDepth = NA_REAL, 
                       bool verbose = false);
