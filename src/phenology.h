#include <RcppArmadillo.h>

#ifndef PHENOLOGY_H
#define PHENOLOGY_H
using namespace Rcpp;

void updatePhenology(List x, int doy, double photoperiod, double tmean);
void updateLeaves(List x, double wind, bool fromGrowthModel);
#endif
