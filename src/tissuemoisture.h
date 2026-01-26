#include <RcppArmadillo.h>

#ifndef TISSUEMOISTURE_H
#define TISSUEMOISTURE_H
using namespace Rcpp;

double sapwoodWaterCapacity(double Al2As, double height, NumericVector V, NumericVector L, double wd);
NumericVector plantWaterContent(List x);
#endif
