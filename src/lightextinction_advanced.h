#include <RcppArmadillo.h>

#ifndef LIGHTEXTINCTION_ADVANCED_H
#define LIGHTEXTINCTION_ADVANCED_H
using namespace Rcpp;

NumericVector leafAngleBetaParameters(double leafAngle, double leafAngleSD);

List instantaneousLightExtinctionAbsortion(NumericMatrix LAIme, NumericMatrix LAImd, NumericMatrix LAImx, 
                                           NumericVector p, NumericVector q, NumericVector ClumpingIndex, 
                                           NumericVector alphaSWR, NumericVector gammaSWR,
                                           DataFrame ddd, int ntimesteps = 24, double trunkExtinctionFraction = 0.1);

void longwaveRadiationSHAW_inner(List internalLWR, NumericMatrix LAIme, NumericMatrix LAImd, NumericMatrix LAImx, 
                                 double LWRatm, double Tsoil, NumericVector Tair, double trunkExtinctionFraction = 0.1);


#endif
