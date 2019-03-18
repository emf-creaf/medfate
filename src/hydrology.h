#include <Rcpp.h>

#ifndef HYDROLOGY_H
#define HYDROLOGY_H
#endif
using namespace Rcpp;

NumericVector er(IntegerVector DOY, double ERconv=0.05, double ERsyn = 0.2);
double soilevaporation(double DEF,double PETs, double Gsoil);
double infiltrationDay(double input, double Ssoil);
double interceptionGashDay(double Precipitation, double Cm, double p, double ER=0.05);
NumericVector infiltrationRepartition(double I, NumericVector dVec, NumericVector macro);
