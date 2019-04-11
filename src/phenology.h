#include <Rcpp.h>

#ifndef PHENOLOGY_H
#define PHENOLOGY_H
#endif
using namespace Rcpp;

NumericVector gdd(IntegerVector DOY, NumericVector Temp, double Tbase = 5.0, double cum = 0.0);

double leafDevelopmentStatus(double Sgdd, double gdd);

NumericVector leafDevelopmentStatus(NumericVector Sgdd, double gdd);

void updateLeaves(List x, double doy, double tmean, double wind, double Tbase = 5.0);