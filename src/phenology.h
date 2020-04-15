#include <Rcpp.h>

#ifndef PHENOLOGY_H
#define PHENOLOGY_H
#endif
using namespace Rcpp;

NumericVector gdd(IntegerVector DOY, NumericVector Temp, double Tbase = 5.0, double cum = 0.0);

double leafDevelopmentStatus(double Sgdd, double gdd);
double leafSenescenceStatus(double Ssen, double sen);

NumericVector leafDevelopmentStatus(NumericVector Sgdd, NumericVector gdd);
NumericVector leafSenescenceStatus(NumericVector Ssen, NumericVector sen);

void updateLeaves(List x, int doy, double photoperiod, double tmean, double wind);