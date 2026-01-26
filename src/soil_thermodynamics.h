#include <RcppArmadillo.h>

#ifndef SOIL_THERMODYNAMICS_H
#define SOIL_THERMODYNAMICS_H
#endif
using namespace Rcpp;

NumericVector layerThermalConductivity(NumericVector sand, NumericVector clay, 
                                       NumericVector W, NumericVector Theta_SAT, NumericVector Theta_FC,
                                       NumericVector Temp);
NumericVector temperatureChange_inner(List SEBcommunication, NumericVector widths, NumericVector Temp,
                                      NumericVector sand, NumericVector clay, 
                                      NumericVector W, NumericVector Theta_SAT, NumericVector Theta_FC,
                                      double Gdown, double tstep);
