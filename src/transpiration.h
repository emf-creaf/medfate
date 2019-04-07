#include <Rcpp.h>

#ifndef TRANSPIRATION_H
#define TRANSPIRATION_H
#endif
using namespace Rcpp;

List profitMaximization(List supplyFunction, DataFrame photosynthesisFunction, int type, double Gwmin, double Gmax, double kleafmax = NA_REAL);
List transpGranier(List x, List soil, double tday, double pet);
List transpSperry(List x, List soil, double tmin, double tmax, double rhmin, double rhmax, double rad, double wind, 
                  double latitude, double elevation, double solarConstant, double delta, double prec,
                  bool verbose = false, int stepFunctions = NA_INTEGER);