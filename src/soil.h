#include <Rcpp.h>

#ifndef SOIL_H
#define SOIL_H
#endif
using namespace Rcpp;

double theta2psi(double clay, double sand, double theta, double om = NA_REAL);
double psi2theta(double clay, double sand, double psi, double om = NA_REAL);
NumericVector waterFC(List soil);

String soilUSDAType(double clay, double sand);
NumericVector vanGenuchtenParams(String soilType);
List soil(List SoilParams, NumericVector W = NumericVector::create(1.0,1.0,1.0));

NumericVector layerthermalconductivity(NumericVector sand, NumericVector clay, NumericVector W);
NumericVector soilthermalconductivity(List soil);
NumericVector soilTemperatureChange(NumericVector dVec, NumericVector Temp,
                                    NumericVector sand, NumericVector clay, 
                                    NumericVector W, NumericVector Theta_FC,
                                    double Gdown);