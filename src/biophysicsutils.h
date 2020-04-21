#include <Rcpp.h>

#ifndef BIOPHYSICS_UTILS_H
#define BIOPHYSICS_UTILS_H
#endif
using namespace Rcpp;

IntegerVector date2doy(CharacterVector dateStrings);
IntegerVector date2photoperiod(CharacterVector dateStrings, double latitude);


double leafTemperature(double absRad, double airTemperature, double u, double E,  double leafWidth = 0.01);
double temperatureDiurnalPattern(double t, double tmin, double tmax, 
                                 double tminPrev, double tmaxPrev, double tminNext, double daylength);
double radiationDiurnalPattern(double t, double daylength);
double irradianceToPhotonFlux(double I, double lambda = 546.6507);

double waterDynamicViscosity(double temp);