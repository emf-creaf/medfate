#include <RcppArmadillo.h>

#ifndef LIGHTEXTINCTION_ADVANCED_H
#define LIGHTEXTINCTION_ADVANCED_H
using namespace Rcpp;

NumericVector leafAngleBetaParameters(double leafAngle, double leafAngleSD);

#endif
