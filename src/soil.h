#include <Rcpp.h>

#ifndef SOIL_H
#define SOIL_H

struct SoilInit {
  std::vector<double> widths;
  std::vector<double> clay;
  std::vector<double> sand;
  std::vector<double> om;
  std::vector<double> nitrogen;
  std::vector<double> ph;
  std::vector<double> bd;
  std::vector<double> rfc;
  std::vector<double> macro;
  std::vector<double> Ksat;
  std::vector<double> VG_alpha;
  std::vector<double> VG_n;
  std::vector<double> VG_theta_res;
  std::vector<double> VG_theta_sat;
  std::vector<double> W;
  std::vector<double> psi;
  std::vector<double> Temp;
  int nlayers;
};

using namespace Rcpp;

/**
 * Conversion factor from conductivity in cm·day-1 to molH20·m-2·MPa-1·s-1
 *  1 day = 86400 sec
 *  1 mol H20 = 18.01528 g
 */
const double cmdTOmmolm2sMPa = 655.2934;//100.0/(18.01528*86400.0*0.00009804139432); 
/**
 * Conversion factor from cm to MPa
 */
const double cmTOMPa = 0.00009804139432;
/**
 * Conversion factor from m to MPa
 */
const double mTOMPa = 0.009804139432;//1/9.804139*0.000001; 

CharacterVector layerNames(int nlayers);

double theta2psiSaxton(double clay, double sand, double theta, double om);
double psi2thetaSaxton(double clay, double sand, double psi, double om);
double theta2psiVanGenuchten(double n, double alpha, double theta_res, double theta_sat, double theta);
double psi2thetaVanGenuchten(double n, double alpha, double theta_res, double theta_sat, double psi);
double psi2kVanGenuchten(double ksat, double n, double alpha, double theta_res, double theta_sat, double psi);
double psi2kVanGenuchtenMicropores(double k_b, double n, double alpha, double theta_res, double theta_sat, 
                                   double psi, double psi_b);
double psi2DVanGenuchten(double k_sat, double n, double alpha, double theta_res, double theta_sat, 
                         double psi);
double psi2cVanGenuchten(double n, double alpha, double theta_res, double theta_sat, double psi);

NumericVector waterExtractable(DataFrame soil, String model, double minPsi);
NumericVector waterSAT(DataFrame soil, String model);
NumericVector thetaSAT(DataFrame soil, String model);
NumericVector waterFC(DataFrame soil, String model);
NumericVector thetaFC(DataFrame soil, String model);
NumericVector waterWP(DataFrame soil, String model);
NumericVector thetaWP(DataFrame soil, String model);
NumericVector water(DataFrame soil, String model);
NumericVector waterPsi(DataFrame soil, double psi, String model);
NumericVector theta(DataFrame soil, String model);
NumericVector psi(DataFrame soil, String model);
NumericVector conductivity(DataFrame soil, String model, bool mmol);
NumericVector capacitance(DataFrame soil, String model);

NumericVector psi2thetasoil(DataFrame soil, NumericVector psi, String model);
  
double saturatedWaterDepth(DataFrame soil, String model);


String USDAType(double clay, double sand);

NumericVector vanGenuchtenParamsCarsel(String soilType);
NumericVector vanGenuchtenParamsToth(double clay, double sand, double om, double bd, bool topsoil);
NumericVector campbellParamsClappHornberger(String soilType);

DataFrame soilInit(DataFrame x, String VG_PTF);

#endif