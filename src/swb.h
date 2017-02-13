#include <Rcpp.h>

#ifndef SWB_H
#define SWB_H
#endif
using namespace Rcpp;

NumericVector er(IntegerVector DOY, double ERconv=0.05, double ERsyn = 0.2);
NumericVector gdd(IntegerVector DOY, NumericVector Temp, double Tbase = 5.0);
double soilevaporation(double DEF,double PETs, double Gsoil);
double infiltrationDay(double NetPrec, double Ssoil);

void checkswbInput(List x, List soil, String transpirationMode);

List swbDay1(List x, List soil, double tday, double pet, double rain, double er, double runon=0.0, bool verbose=false);
List swbDay2(List x, List soil, double tmin, double tmax, double rhmin, double rhmax, double rad, double wind, 
             double latitude, double elevation, double slope, double aspect, double delta, 
             double rain, double er, double runon=0.0, bool verbose = false);