#include <Rcpp.h>

#ifndef SOIL_H
#define SOIL_H
#endif
using namespace Rcpp;

double theta2psiSaxton(double clay, double sand, double theta, double om = NA_REAL);
double psi2thetaSaxton(double clay, double sand, double psi, double om = NA_REAL);
double theta2psiVanGenuchten(double n, double alpha, double theta_res, double theta_sat, double theta);
double psi2thetaVanGenuchten(double n, double alpha, double theta_res, double theta_sat, double psi);

NumericVector waterSAT(List soil, String model = "SX");
NumericVector thetaSAT(List soil, String model = "SX");
NumericVector waterFC(List soil, String model = "SX");
NumericVector thetaFC(List soil, String model = "SX");
NumericVector waterWP(List soil, String model = "SX");
NumericVector thetaWP(List soil, String model = "SX");
NumericVector theta(List soil, String model="SX");
NumericVector psi(List soil, String model="SX");
double waterTableDepth(List soil, String model = "SX");

NumericVector thermalConductivity(List soil, String model = "SX");
NumericVector thermalCapacity(List soil, String model = "SX");

String USDAType(double clay, double sand);

NumericVector vanGenuchtenParamsCarsel(String soilType);
NumericVector vanGenuchtenParamsToth(double clay, double sand, double om, double bd, bool topsoil);
  
List soil(List SoilParams, String VG_PTF = "Carsel", NumericVector W = NumericVector::create(1.0,1.0,1.0));

NumericVector layerThermalConductivity(NumericVector sand, NumericVector clay, NumericVector W, NumericVector Theta_FC);
NumericVector temperatureChange(NumericVector dVec, NumericVector Temp,
                                NumericVector sand, NumericVector clay, 
                                NumericVector W, NumericVector Theta_FC,
                                double Gdown);