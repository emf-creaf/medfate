#include <Rcpp.h>

#ifndef LIGHTEXTINCTION_ADVANCED_H
#define LIGHTEXTINCTION_ADVANCED_H
#endif
using namespace Rcpp;

NumericVector leafAngleBetaParameters(double leafAngle, double leafAngleSD);

NumericVector layerDirectIrradianceFraction(NumericMatrix LAIme, NumericMatrix LAImd,NumericMatrix LAImx, 
                                            NumericVector kb, NumericVector ClumpingIndex, 
                                            NumericVector alpha, NumericVector gamma, double trunkExtinctionFraction = 0.1);

NumericMatrix layerDiffuseIrradianceFraction(NumericMatrix LAIme, NumericMatrix LAImd,NumericMatrix LAImx, 
                                             NumericMatrix K, NumericVector ClumpingIndex, NumericVector ZF,
                                             NumericVector alpha, NumericVector gamma, double trunkExtinctionFraction = 0.1);

List cohortSunlitShadeAbsorbedRadiation(double Ib0, double Id0,
                                        NumericMatrix LAIme, NumericMatrix LAImd, NumericMatrix LAImx,
                                        NumericVector kb, NumericMatrix K, NumericVector ClumpingIndex, NumericVector ZF, 
                                        NumericVector alpha, NumericVector gamma, double trunkExtinctionFraction = 0.1);

NumericVector layerSunlitFraction(NumericMatrix LAIme, NumericMatrix LAImd, 
                                  NumericVector kb, NumericVector ClumpingIndex);

List instantaneousLightExtinctionAbsortion(NumericMatrix LAIme, NumericMatrix LAImd, NumericMatrix LAImx, 
                                           NumericVector p, NumericVector q, NumericVector ClumpingIndex, 
                                           NumericVector alphaSWR, NumericVector gammaSWR,
                                           DataFrame ddd, int ntimesteps = 24, double trunkExtinctionFraction = 0.1);

void longwaveRadiationSHAW_inner(List internalLWR, NumericMatrix LAIme, NumericMatrix LAImd, NumericMatrix LAImx, 
                                 double LWRatm, double Tsoil, NumericVector Tair, double trunkExtinctionFraction = 0.1);