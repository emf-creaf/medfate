#include <Rcpp.h>

#ifndef SOIL_THERMODYNAMICS_H
#define SOIL_THERMODYNAMICS_H
#endif
using namespace Rcpp;

NumericVector layerThermalConductivity(NumericVector sand, NumericVector clay, NumericVector W, NumericVector Theta_FC);
NumericVector temperatureChange(NumericVector dVec, NumericVector Temp,
                                NumericVector sand, NumericVector clay, 
                                NumericVector W, NumericVector Theta_FC,
                                double Gdown, double tstep);