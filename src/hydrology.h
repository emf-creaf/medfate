#include <Rcpp.h>

#ifndef HYDROLOGY_H
#define HYDROLOGY_H
#endif
using namespace Rcpp;

double soilEvaporationAmount(double DEF,double PETs, double Gsoil);
double herbaceousTranspiration(double pet, double LherbSWR, double herbLAI, 
                               List soil, String soilFunctions, bool modifySoil = true);
NumericVector soilEvaporation(List soil, String soilFunctions, double pet, double LgroundSWR,
                              bool modifySoil = true);


double infiltrationAmount(double input, double Ssoil);
NumericVector infiltrationRepartition(double I, NumericVector dVec, NumericVector macro, 
                                      double a = -0.005, double b = 3.0);

double interceptionGashDay(double Precipitation, double Cm, double p, double ER=0.05);

double erFactor(int doy, double pet, double prec, double Rconv = 5.6, double Rsyn = 1.5);

double snowMelt(double tday, double rad, double LgroundSWR, double elevation);

NumericVector soilWaterInputs(List soil, String soilFunctions, double prec, double er, double tday, double rad, double elevation,
                              double Cm, double LgroundPAR, double LgroundSWR, 
                              double runon = 0.0,
                              bool snowpack = true, bool modifySoil = true);
NumericVector soilInfiltrationPercolation(List soil, String soilFunctions, 
                                          double waterInput,
                                          bool modifySoil = true);
double soilFlows(List soil, NumericVector sourceSink, int nsteps = 24, 
                 String lowerBoundary = "free",
                 bool modifySoil = true);