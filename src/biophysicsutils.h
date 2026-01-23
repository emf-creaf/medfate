#include <Rcpp.h>

#ifndef BIOPHYSICS_UTILS_H
#define BIOPHYSICS_UTILS_H
using namespace Rcpp;

IntegerVector date2doy(CharacterVector dateStrings);
NumericVector date2photoperiod(CharacterVector dateStrings, double latitude);

#endif
