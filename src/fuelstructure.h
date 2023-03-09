#include <Rcpp.h>

#ifndef FUELSTRUCTURE_H
#define FUELSTRUCTURE_H
#endif
using namespace Rcpp;

List fuelLiveStratification(List x, DataFrame SpParams, double gdd = NA_REAL, 
                            double heightProfileStep = 10.0, double maxHeightProfile = 5000.0, double bulkDensityThreshold = 0.05);
DataFrame FCCSproperties(List object, DataFrame SpParams, NumericVector cohortFMC = NumericVector::create(), 
                         double gdd = NA_REAL, 
                         double heightProfileStep = 10.0, double maxHeightProfile = 5000, double bulkDensityThreshold = 0.05,
                         String depthMode = "crownaverage");

double layerFuelAverageParameter(double minHeight, double maxHeight, NumericVector cohortParameter, NumericVector cohortLoading, NumericVector H, NumericVector CR);
NumericVector woodyFuelProfile(NumericVector z, List x, DataFrame SpParams,  
                               double gdd = NA_REAL);
