#include <Rcpp.h>

#ifndef FUELSTRUCTURE_H
#define FUELSTRUCTURE_H
#endif
using namespace Rcpp;

List fuelLiveStratification(List x, DataFrame SpParams, double gdd = NA_REAL, String mode = "MED",
                            double heightProfileStep = 10.0, double maxHeightProfile = 5000.0, double bulkDensityThreshold = 0.05);
DataFrame FCCSproperties(List object, double ShrubCover, double CanopyCover, DataFrame SpParams, NumericVector cohortFMC = NumericVector::create(), 
                         double gdd = NA_REAL, String mode = "MED", 
                         double heightProfileStep = 10.0, double maxHeightProfile = 5000, double bulkDensityThreshold = 0.05);
// List fuelStructure(List object, DataFrame SpParams, DataFrame FuelModelParams, 
//                    double gdd = NA_REAL, String mode = "MED", 
//                    double heightProfileStep = 10.0, double maxHeightProfile = 5000, double bulkDensityThreshold = 0.05, bool useModelForLive = false);
NumericVector woodyFuelProfile(NumericVector z, List x, DataFrame SpParams,  
                               double gdd = NA_REAL, String mode = "MED");
