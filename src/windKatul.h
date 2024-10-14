#include <Rcpp.h>

#ifndef WINDKATUL_H
#define WINDKATUL_H
#endif
using namespace Rcpp;

NumericMatrix windCanopyTurbulenceModel_inner(NumericVector zm, NumericVector Cx, double hm, double d0, double z0,
                                              String model = "k-epsilon");

void windCanopyTurbulence_inner(DataFrame output, NumericVector zmid, NumericVector LAD, double canopyHeight,
                                double u, double windMeasurementHeight = 200, String model = "k-epsilon");