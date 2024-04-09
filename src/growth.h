#include <Rcpp.h>

#ifndef GROWTH_H
#define GROWTH_H
#endif
using namespace Rcpp;

List growthDay(List x, CharacterVector date, NumericVector meteovec, 
               double latitude, double elevation, double slope = NA_REAL, double aspect = NA_REAL,  
               double runon = 0.0, Nullable<NumericVector> lateralFlows = R_NilValue, double waterTableDepth = NA_REAL, 
               bool modifyInput = true);