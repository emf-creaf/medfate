#include <Rcpp.h>

#ifndef SOIL_H
#define SOIL_H
#endif
using namespace Rcpp;

double theta2psi(double clay, double sand, double theta);
double psi2theta(double clay, double sand, double psi);
String soilUSDAType(double clay, double sand);
NumericVector vanGenuchtenParams(String soilType);
List soil(List SoilParams, NumericVector W = NumericVector::create(1.0,1.0,1.0));

