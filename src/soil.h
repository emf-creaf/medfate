#include <Rcpp.h>

#ifndef SOIL_H
#define SOIL_H
#endif
using namespace Rcpp;

const double cmdTOmmolm2sMPa = 655.2934;
const double cmTOMPa = 0.00009804139432;
const double mTOMPa = 0.009804139432;

CharacterVector layerNames(int nlayers);

double saturatedConductivitySaxton(double clay, double sand, double bd, double om = NA_REAL, bool mmol = true);
  
double theta2psiSaxton(double clay, double sand, double theta, double om = NA_REAL);
double psi2thetaSaxton(double clay, double sand, double psi, double om = NA_REAL);
double theta2psiVanGenuchten(double n, double alpha, double theta_res, double theta_sat, double theta);
double psi2thetaVanGenuchten(double n, double alpha, double theta_res, double theta_sat, double psi);
double psi2kVanGenuchten(double ksat, double n, double alpha, double theta_res, double theta_sat, double psi);
double psi2kVanGenuchtenMicropores(double k_b, double n, double alpha, double theta_res, double theta_sat, 
                                   double psi, double psi_b);
double psi2DVanGenuchten(double k_sat, double n, double alpha, double theta_res, double theta_sat, 
                         double psi);
double psi2cVanGenuchten(double n, double alpha, double theta_res, double theta_sat, double psi);

NumericVector waterExtractable(DataFrame soil, String model="SX", double minPsi = -5.0);
NumericVector waterSAT(DataFrame soil, String model = "SX");
NumericVector thetaSAT(DataFrame soil, String model = "SX");
NumericVector waterFC(DataFrame soil, String model = "SX");
NumericVector thetaFC(DataFrame soil, String model = "SX");
NumericVector waterWP(DataFrame soil, String model = "SX");
NumericVector thetaWP(DataFrame soil, String model = "SX");
NumericVector water(DataFrame soil, String model="SX");
NumericVector waterPsi(DataFrame soil, double psi, String model="SX");
NumericVector theta(DataFrame soil, String model="SX");
NumericVector psi(DataFrame soil, String model="SX");
NumericVector conductivity(DataFrame soil, String model="SX");
NumericVector capacitance(DataFrame soil, String model="SX");

NumericVector psi2thetasoil(DataFrame soil, NumericVector psi, String model="SX");
  
double saturatedWaterDepth(DataFrame soil, String model = "SX");


String USDAType(double clay, double sand);

NumericVector vanGenuchtenParamsCarsel(String soilType);
NumericVector vanGenuchtenParamsToth(double clay, double sand, double om, double bd, bool topsoil);
NumericVector campbellParamsClappHornberger(String soilType);

DataFrame soilInit(DataFrame x, String VG_PTF = "Carsel");
