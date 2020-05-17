#include <Rcpp.h>

#ifndef PHENOLOGY_H
#define PHENOLOGY_H
#endif
using namespace Rcpp;

NumericVector gdd(IntegerVector DOY, NumericVector Temp, double Tbase = 5.0, double cum = 0.0);

double leafDevelopmentStatus(double Sgdd, double gdd, double unfoldingDD = 300.0);
bool leafSenescenceStatus(double Ssen, double sen);

NumericVector leafDevelopmentStatus(NumericVector Sgdd, NumericVector gdd, double unfoldingDD = 300.0);
LogicalVector leafSenescenceStatus(NumericVector Ssen, NumericVector sen);

void updatePhenology(List x, int doy, double photoperiod, double tmean);
void updateLeaves(List x, double wind, bool fromGrowthModel);