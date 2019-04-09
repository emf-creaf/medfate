#include <Rcpp.h>

#ifndef HYDROLOGY_H
#define HYDROLOGY_H
#endif
using namespace Rcpp;

double soilEvaporationAmount(double DEF,double PETs, double Gsoil);
NumericVector soilEvaporation(List soil, String soilFunctions, double pet, double LgroundSWR,
                              bool modifySoil = true);

double infiltrationAmount(double input, double Ssoil);
NumericVector infiltrationRepartition(double I, NumericVector dVec, NumericVector macro);

double interceptionGashDay(double Precipitation, double Cm, double p, double ER=0.05);

double erFactor(int doy, double pet, double prec, double Rconv = 5.6, double Rsyn = 1.5);
NumericVector verticalInputs(List soil, String soilFunctions, double prec, double er, double tday, double rad, double elevation,
                             double Cm, double LgroundPAR, double LgroundSWR, 
                             double runon = 0.0,
                             bool snowpack = true, bool drainage = true, bool modifySoil = true);