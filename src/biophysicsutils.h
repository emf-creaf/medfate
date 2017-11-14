#include <Rcpp.h>

#ifndef BIOPHYSICS_UTILS_H
#define BIOPHYSICS_UTILS_H
#endif
using namespace Rcpp;


double leafTemperature(double absRad, double airTemperature, double u, double E,  double leafWidth = 0.01);
double temperatureDiurnalPattern(double t, double tmin, double tmax, double daylength);
double radiationDiurnalPattern(double t, double daylength);
double irradianceToPhotonFlux(double I, double lambda = 546.6507);