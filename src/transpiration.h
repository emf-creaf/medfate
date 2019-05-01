#include <Rcpp.h>

#ifndef TRANSPIRATION_H
#define TRANSPIRATION_H
#endif
using namespace Rcpp;

List profitMaximization(List supplyFunction, DataFrame photosynthesisFunction, int type, double Gwmin, double Gmax, double kleafmax = NA_REAL);
List transpirationGranier(List x, List soil, double tday, double pet, bool modifyInput = true);
List transpirationSperry(List x, List soil, double tmin, double tmax, double rhmin, double rhmax, double rad, double wind, 
                  double latitude, double elevation, double slope, double aspect,
                  double solarConstant, double delta, double prec,
                  double canopyEvaporation = 0.0, double snowMelt = 0.0, double soilEvaporation = 0.0,
                  bool verbose = false, int stepFunctions = NA_INTEGER, bool modifyInput = true);