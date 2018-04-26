#include <Rcpp.h>

#ifndef SOIL_H
#define SOIL_H
#endif
using namespace Rcpp;

double theta2psiSaxton(double clay, double sand, double theta, double om = NA_REAL);
double psi2thetaSaxton(double clay, double sand, double psi, double om = NA_REAL);
double theta2psiVanGenuchten(double n, double alpha, double theta_res, double theta_sat, double theta);
double psi2thetaVanGenuchten(double n, double alpha, double theta_res, double theta_sat, double psi);

NumericVector waterFC(List soil, String model = "SX");
NumericVector thetaFC(List soil, String model = "SX");
NumericVector theta(List soil, String model="SX");
NumericVector psi(List soil, String model="SX");

NumericVector soilthermalconductivity(List soil, String model = "SX");
NumericVector soilthermalcapacity(List soil, String model = "SX");

String soilUSDAType(double clay, double sand);
NumericVector vanGenuchtenParams(String soilType);
List soil(List SoilParams, NumericVector W = NumericVector::create(1.0,1.0,1.0));

NumericVector layerthermalconductivity(NumericVector sand, NumericVector clay, NumericVector W, NumericVector Theta_FC);
NumericVector soilTemperatureChange(NumericVector dVec, NumericVector Temp,
                                    NumericVector sand, NumericVector clay, 
                                    NumericVector W, NumericVector Theta_FC,
                                    double Gdown);