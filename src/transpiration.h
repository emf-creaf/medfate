#include <Rcpp.h>

#ifndef TRANSPIRATION_H
#define TRANSPIRATION_H
#endif
using namespace Rcpp;

List profitMaximization(List supplyFunction, DataFrame photosynthesisFunction, double Gswmin, double Gswmax);
List transpirationBasic(List x, NumericVector meteovec,
                        double elevation, bool modifyInput = true,
                        String transpirationMode = "Granier");
List transpirationAdvanced(List x, NumericVector meteovec,
                  double latitude, double elevation, double slope, double aspect,
                  double solarConstant, double delta,
                  double canopyEvaporation = 0.0, double snowMelt = 0.0, double soilEvaporation = 0.0,
                  bool verbose = false, int stepFunctions = NA_INTEGER, 
                  bool modifyInput = true, String transpirationMode = "Sperry");