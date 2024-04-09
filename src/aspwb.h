#include <Rcpp.h>

#ifndef ASPWB_H
#define ASPWB_H
#endif
using namespace Rcpp;

List aspwbInput(double crop_factor, List control, List soil);

List aspwb_day_internal(List x, NumericVector meteovec, 
                        double elevation, double slope, double aspect, 
                        double runon =  0.0, Nullable<NumericVector> lateralFlows = R_NilValue, double waterTableDepth = NA_REAL,
                        bool verbose = false);

List aspwb_day(List x, CharacterVector date, NumericVector meteovec, 
               double latitude, double elevation, double slope = NA_REAL, double aspect = NA_REAL,
               double runon =  0.0, Nullable<NumericVector> lateralFlows = R_NilValue, double waterTableDepth = NA_REAL,
               bool modifyInput = true);