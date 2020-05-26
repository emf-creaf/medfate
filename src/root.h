#include <Rcpp.h>

#ifndef ROOT_H
#define ROOT_H
#endif
using namespace Rcpp;


NumericVector ldrRS_one(double Z50, double Z95, NumericVector d);
NumericVector conicRS_one(double Z, NumericVector d);
NumericMatrix conicDistribution(NumericVector Z, NumericVector d);
NumericMatrix ldrDistribution(NumericVector Z50, NumericVector Z95, NumericVector d);

double specificRootSurfaceArea(double specificRootLength, double rootTissueDensity);
double fineRootArea(double vgrhizo_kmax, double leafArea);
double fineRootBiomass(double vgrhizo_kmax, double leafArea, 
                       double specificRootLength, double rootTissueDensity);
double fineRootSoilVolume(double fineRootBiomass, double specificRootLength, double rootLengthDensity = 10.0);

NumericVector rootLengths(NumericVector v, NumericVector d, double depthWidthRatio = 1.0);
NumericVector xylemConductanceProportions(NumericVector v, NumericVector d, double depthWidthRatio = 1.0);
List horizontalProportions(NumericMatrix V, NumericVector poolProportions, double LAIcell,
                           double poolOverlapFactor);
List horizontalProportionsAdvanced(NumericVector poolProportions, NumericVector VolInd, 
                                   NumericVector N, NumericMatrix V, 
                                   NumericVector d, NumericVector bulkDensity);