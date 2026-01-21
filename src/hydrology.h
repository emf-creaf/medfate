#include <Rcpp.h>

#ifndef HYDROLOGY_H
#define HYDROLOGY_H
#endif
using namespace Rcpp;

NumericVector herbaceousTranspiration(double pet, double LherbSWR, double herbLAI, 
                                      DataFrame soil, String soilFunctions, bool modifySoil = true);
double soilEvaporation(DataFrame soil, double snowpack, 
                       String soilFunctions, double pet, double LgroundSWR,
                       bool modifySoil = true);

double infiltrationAmount(double rainfallInput, double rainfallIntensity, DataFrame soil, 
                          String soilFunctions, String model = "GreenAmpt1911", double K_correction = 1.0);
NumericVector infiltrationRepartition(double I, NumericVector widths, NumericVector macro, 
                                      double a = -0.005, double b = 3.0);

double rainfallIntensity(int month, double prec, NumericVector rainfallIntensityPerMonth);

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