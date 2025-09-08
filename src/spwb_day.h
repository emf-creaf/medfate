#include <Rcpp.h>

#ifndef SPWB_DAY_H
#define SPWB_DAY_H
#endif
using namespace Rcpp;

void fccsHazard(NumericVector fireHazard, List x, NumericVector meteovec, List transpOutput, double slope);

List spwbDay(List x, CharacterVector date, NumericVector meteovec,
             double latitude, double elevation, double slope = NA_REAL, double aspect = NA_REAL,  
             double runon = 0.0, Nullable<NumericVector> lateralFlows = R_NilValue, double waterTableDepth = NA_REAL,
             bool modifyInput = true);
void spwbDay_basic(List internalCommunication, List x, NumericVector meteovec, 
                   double elevation, double slope, double aspect,
                   double runon = 0.0, Nullable<NumericVector> lateralFlows = R_NilValue, double waterTableDepth = NA_REAL, 
                   bool verbose = false);
void spwbDay_advanced(List internalCommunication, List x, NumericVector meteovec, 
                      double latitude, double elevation, double slope, double aspect,
                      double solarConstant, double delta, 
                      double runon = 0.0, Nullable<NumericVector> lateralFlows = R_NilValue, double waterTableDepth = NA_REAL, 
                      bool verbose = false);
