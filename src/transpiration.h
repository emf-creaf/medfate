#include <Rcpp.h>

#ifndef TRANSPIRATION_H
#define TRANSPIRATION_H
#endif
using namespace Rcpp;

List profitMaximization(List supplyFunction, DataFrame photosynthesisFunction, double Gswmin, double Gswmax, 
                        double gainModifier = 1.0, double costModifier = 1.0, String costWater = "dEdP");
List transpirationGranier(List x, NumericVector meteovec, 
                          bool modifyInput = true);
List transpirationSperry(List x, NumericVector meteovec,
                  double latitude, double elevation, double slope, double aspect,
                  double solarConstant, double delta,
                  double canopyEvaporation = 0.0, double snowMelt = 0.0, double soilEvaporation = 0.0,
                  bool verbose = false, int stepFunctions = NA_INTEGER, 
                  bool modifyInput = true);