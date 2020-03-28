#include <Rcpp.h>

#ifndef ROOT_H
#define ROOT_H
#endif
using namespace Rcpp;


NumericVector ldrRS_one(double Z50, double Z95, NumericVector d);
NumericVector conicRS_one(double Z, NumericVector d);
NumericMatrix conicDistribution(NumericVector Z, NumericVector d);
NumericMatrix ldrDistribution(NumericVector Z50, NumericVector Z95, NumericVector d);
NumericVector rootLengths(NumericVector v, NumericVector d, double depthWidthRatio = 1.0);
NumericVector xylemConductanceProportions(NumericVector v, NumericVector d, double depthWidthRatio = 1.0);
List horizontalProportions(NumericMatrix V, NumericVector LAIlive, double poolOverlapFactor);