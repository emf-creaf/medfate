#include <Rcpp.h>

#ifndef HYDROLOGY_H
#define HYDROLOGY_H
#endif
using namespace Rcpp;

double soilEvaporationAmount(double DEF,double PETs, double Gsoil);
NumericVector herbaceousTranspiration(double pet, double LherbSWR, double herbLAI, 
                                      List soil, String soilFunctions, bool modifySoil = true);
double soilEvaporation(List soil, String soilFunctions, double pet, double LgroundSWR,
                       bool modifySoil = true);


double infiltrationBoughton(double input, double Ssoil);
double infitrationGreenAmpt(double t, double Psi_w, double Ksat, double theta_sat, double theta_dry);
double infiltrationAmount(double rainfallInput, double rainfallIntensity, List soil, String soilFunctions, 
                          String model = "Green-Ampt");
NumericVector infiltrationRepartition(double I, NumericVector dVec, NumericVector macro, 
                                      double a = -0.005, double b = 3.0);

double rainfallIntensity(int doy, double prec, double Rconv = 5.6, double Rsyn = 1.5);

double interceptionGashDay(double Rainfall, double Cm, double p, double ER=0.05);
double interceptionLiuDay(double Rainfall, double Cm, double p, double ER=0.05);

double snowMelt(double tday, double rad, double LgroundSWR, double elevation);

NumericVector soilWaterInputs(List soil, String soilFunctions, String interceptionMode,
                              double prec, double rainfallIntensity,
                              double pet, double tday, double rad, double elevation,
                              double Cm, double LgroundPAR, double LgroundSWR, 
                              double runon = 0.0,
                              bool snowpack = true, bool modifySoil = true);
double soilFlows(List soil, NumericVector sourceSink, int nsteps = 24,
                 bool modifySoil = true);