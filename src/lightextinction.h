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
NumericVector cohortAbsorbedSWRFraction(NumericVector LAI_expanded, NumericVector LAI_dead, NumericVector H, NumericVector CR, NumericVector kPAR);

NumericVector layerIrradianceFraction(NumericMatrix LAIme, NumericMatrix LAImd, NumericVector k, NumericVector alpha);
NumericVector layerIrradianceFractionBottomUp(NumericMatrix LAIme, NumericMatrix LAImd, NumericVector k, NumericVector alpha);

double groundIrradianceFraction(NumericMatrix LAIme, NumericMatrix LAImd, NumericVector k, NumericVector alpha);

List cohortSunlitShadeAbsorbedRadiation(double Ib0, double Id0, NumericVector Ibf, NumericVector Idf, double beta,
                                   NumericMatrix LAIme, NumericMatrix LAImd, 
                                   NumericVector kb,  NumericVector kd, NumericVector alpha, double gamma);

NumericVector layerSunlitFraction(NumericMatrix LAIme, NumericMatrix LAImd, NumericVector kb);

List instantaneousLightExtinctionAbsortion(NumericMatrix LAIme, NumericMatrix LAImd, NumericVector kPAR,
                                           DataFrame ddd, NumericVector LWR_diffuse,
                                           int ntimesteps = 24, String canopyMode= "sunshade");