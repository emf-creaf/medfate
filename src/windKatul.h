#include <RcppArmadillo.h>

#ifndef WINDKATUL_H
#define WINDKATUL_H
using namespace Rcpp;

void windCanopyTurbulence_inner(DataFrame output, NumericVector zmid, NumericVector LAD, double canopyHeight,
                                double u, double windMeasurementHeight, String model);

#endif
