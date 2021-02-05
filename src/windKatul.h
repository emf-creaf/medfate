#include <Rcpp.h>

#ifndef WINDKATUL_H
#define WINDKATUL_H
#endif
using namespace Rcpp;

DataFrame windCanopyTurbulenceModel(NumericVector zm, NumericVector Cx, double hm, double d0, double z0,
                                     String model = "k-epsilon");

DataFrame windCanopyTurbulence(NumericVector zmid, NumericVector LAD, double canopyHeight,
                                double u2m, String model = "k-epsilon");