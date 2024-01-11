#include <Rcpp.h>

#ifndef LIGHTEXTINCTION_H
#define LIGHTEXTINCTION_H
#endif
using namespace Rcpp;

double availableLight(double h, NumericVector H, NumericVector LAI, NumericVector k, NumericVector CR);
NumericVector parcohortC(NumericVector H, NumericVector LAI_expanded, NumericVector LAI_dead, NumericVector k, NumericVector CR);
NumericVector parcohort(IntegerVector SP, NumericVector H, NumericVector CR, NumericVector LAI, DataFrame SpParams);
NumericVector parheight(NumericVector heights, IntegerVector SP, NumericVector H, NumericVector CR, NumericVector LAI, DataFrame SpParams);
NumericVector swrheight(NumericVector heights, IntegerVector SP, NumericVector H, NumericVector CR, NumericVector LAI, DataFrame SpParams);
NumericVector cohortAbsorbedSWRFraction(NumericVector z, NumericVector LAI_expanded, NumericVector LAI_dead, NumericVector H, NumericVector CR, NumericVector kPAR);

NumericVector layerIrradianceFraction(NumericMatrix LAIme, NumericMatrix LAImd, NumericMatrix LAImx, NumericVector k, NumericVector alpha, double trunkExtinctionFraction = 0.1);

double groundIrradianceFraction(NumericMatrix LAIme, NumericMatrix LAImd, NumericMatrix LAImx, NumericVector k, NumericVector alpha, double trunkExtinctionFraction = 0.1);

List cohortSunlitShadeAbsorbedRadiation(double Ib0, double Id0, NumericVector Ibf, NumericVector Idf,
                                   NumericMatrix LAIme, NumericMatrix LAImd, 
                                   NumericVector kb,  NumericVector kd, NumericVector alpha, NumericVector gamma);

NumericVector layerSunlitFraction(NumericMatrix LAIme, NumericMatrix LAImd, NumericVector kb);

List instantaneousLightExtinctionAbsortion(NumericMatrix LAIme, NumericMatrix LAImd, NumericMatrix LAImx, 
                                           NumericVector LeafAngle, NumericVector kPAR, NumericVector alphaSWR, NumericVector gammaSWR,
                                           DataFrame ddd, int ntimesteps = 24, double trunkExtinctionFraction = 0.1);

List longwaveRadiationSHAW(NumericMatrix LAIme, NumericMatrix LAImd, NumericMatrix LAImx, 
                           double LWRatm, double Tsoil, NumericVector Tair, double trunkExtinctionFraction = 0.1);