#include <Rcpp.h>

#ifndef HYDROLOGY_H
#define HYDROLOGY_H
#endif
using namespace Rcpp;

double soilEvaporationAmount(double DEF,double PETs, double Gsoil);
NumericVector herbaceousTranspiration(double pet, double LherbSWR, double herbLAI, 
                                      DataFrame soil, String soilFunctions, bool modifySoil = true);
double soilEvaporation(DataFrame soil, double snowpack, 
                       String soilFunctions, double pet, double LgroundSWR,
                       bool modifySoil = true);


double infiltrationBoughton_c(double input, double Ssoil);
double infitrationGreenAmpt_c(double t, double Psi_w, double Ksat, double theta_sat, double theta_dry);

double infiltrationAmount(double rainfallInput, double rainfallIntensity, DataFrame soil, 
                          String soilFunctions, String model = "GreenAmpt1911", double K_correction = 1.0);
NumericVector infiltrationRepartition(double I, NumericVector widths, NumericVector macro, 
                                      double a = -0.005, double b = 3.0);

double rainfallIntensity(int month, double prec, NumericVector rainfallIntensityPerMonth);
double rainfallIntensity_c(int month, double prec, std::vector<double> rainfallIntensityPerMonth);
double interceptionGashDay_c(double Rainfall, double Cm, double p, double ER=0.05);
double interceptionLiuDay_c(double Rainfall, double Cm, double p, double ER=0.05);
double snowMelt_c(double tday, double rad, double LgroundSWR, double elevation);

NumericVector waterInputs(List x, 
                          double prec, double rainfallIntensity,
                          double pet, double tday, double rad, double elevation,
                          double Cm, double LgroundPAR, double LgroundSWR, 
                          bool modifyInput = true);
NumericVector soilWaterBalance_inner(List SWBcommunication, DataFrame soil, String soilFunctions, 
                               double rainfallInput, double rainfallIntensity, double snowmelt, NumericVector sourceSink, 
                               double runon = 0.0, Nullable<NumericVector> lateralFlows = R_NilValue, double waterTableDepth = NA_REAL,
                               String infiltrationMode = "GreenAmpt1911", double infiltrationCorrection = 5.0, String soilDomains = "buckets", 
                               int nsteps = 24, int max_nsubsteps = 3600, bool modifySoil = true);