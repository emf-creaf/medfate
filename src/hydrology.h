#include <Rcpp.h>

#ifndef HYDROLOGY_H
#define HYDROLOGY_H
#endif
using namespace Rcpp;

double soilEvaporationDay(double DEF,double PETs, double Gsoil);
double infiltrationDay(double input, double Ssoil);
double interceptionGashDay(double Precipitation, double Cm, double p, double ER=0.05);

double erFactor(int doy, double pet, double prec, double Rconv = 5.6, double Rsyn = 1.5);
NumericVector soilEvaporation(List soil, String soilFunctions, double pet, double LgroundSWR,
                              bool modifySoil = true);
NumericVector infiltrationRepartition(double I, NumericVector dVec, NumericVector macro);
NumericVector verticalInputs(List soil, String soilFunctions, double prec, double er, double tday, double rad, double elevation,
                             double Cm, double LgroundPAR, double LgroundSWR, 
                             double runon = 0.0,
                             bool snowpack = true, bool drainage = true, bool modifySoil = true);