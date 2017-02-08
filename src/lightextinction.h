#include <Rcpp.h>

#ifndef LIGHTEXTINCTION_H
#define LIGHTEXTINCTION_H
#endif
using namespace Rcpp;

double availableLight(double h, NumericVector H, NumericVector LAI, NumericVector k, NumericVector CR);
NumericVector parcohortC(NumericVector H, NumericVector LAI, NumericVector k, NumericVector CR);
NumericVector parcohort(IntegerVector SP, NumericVector H, NumericVector CR, NumericVector LAI, DataFrame SpParams);
NumericVector parheight(NumericVector heights, IntegerVector SP, NumericVector H, NumericVector CR, NumericVector LAI, DataFrame SpParams);
NumericVector swrheight(NumericVector heights, IntegerVector SP, NumericVector H, NumericVector CR, NumericVector LAI, DataFrame SpParams);
NumericVector cohortAbsorbedSWRFraction(NumericVector z, NumericVector LAI, NumericVector H, NumericVector CR, NumericVector kPAR);
NumericVector cohortAbsorbedSWRFraction(NumericVector LAI, NumericVector H, NumericVector CR, NumericVector kPAR);
