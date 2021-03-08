#include <Rcpp.h>

#ifndef TRANSPIRATION_H
#define TRANSPIRATION_H
#endif
using namespace Rcpp;

List profitMaximization(List supplyFunction, DataFrame photosynthesisFunction, double Gswmin, double Gswmax, 
                        double gainModifier = 1.0, double costModifier = 1.0, String costWater = "dEdP");
List transpirationGranier(List x, double tday, double pet, 
                          bool modifyInputX = true, bool modifyInputSoil = true);
List transpirationSperry(List x, double tmin, double tmax, 
                  double tminPrev, double tmaxPrev, double tminNext, double rhmin, double rhmax, double rad, double wind, 
                  double latitude, double elevation, double slope, double aspect,
                  double solarConstant, double delta, double prec,
                  double canopyEvaporation = 0.0, double snowMelt = 0.0, double soilEvaporation = 0.0,
                  bool verbose = false, int stepFunctions = NA_INTEGER, 
                  bool modifyInputX = true, bool modifyInputSoil = true);