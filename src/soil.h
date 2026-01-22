#include <Rcpp.h>
#include "soil_c.h"

#ifndef SOIL_H
#define SOIL_H

using namespace Rcpp;

CharacterVector layerNames(int nlayers);

NumericVector waterExtractable(DataFrame soil, String model, double minPsi);
NumericVector waterSAT(DataFrame soil, String model);
NumericVector thetaSAT(DataFrame soil, String model);
NumericVector waterFC(DataFrame soil, String model);
NumericVector thetaFC(DataFrame soil, String model);
NumericVector waterWP(DataFrame soil, String model);
NumericVector thetaWP(DataFrame soil, String model);
NumericVector water(DataFrame soil, String model);
NumericVector waterPsi(DataFrame soil, double psi, String model);
NumericVector theta(DataFrame soil, String model);
NumericVector psi(DataFrame soil, String model);
NumericVector conductivity(DataFrame soil, String model, bool mmol);
NumericVector capacitance(DataFrame soil, String model);

NumericVector psi2thetasoil(DataFrame soil, NumericVector psi, String model);
  
double saturatedWaterDepth(DataFrame soil, String model);


String USDAType(double clay, double sand);

NumericVector vanGenuchtenParamsCarsel(String soilType);
NumericVector vanGenuchtenParamsToth(double clay, double sand, double om, double bd, bool topsoil);
NumericVector campbellParamsClappHornberger(String soilType);

DataFrame soilInit(DataFrame x, String VG_PTF);
Soil soilDataFrameToStructure(DataFrame x, String model);

#endif
