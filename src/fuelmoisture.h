#include <RcppArmadillo.h>

#ifndef FUELMOISTURE_H
#define FUELMOISTURE_H
#endif
using namespace Rcpp;

DataFrame deadFuelMoisture(NumericVector m0, NumericVector airTemp, NumericVector airHumidity, NumericVector fuelRadiation, NumericVector fuelWindSpeed, NumericVector effRain, NumericVector rainDuration);
double canopyLiveFuelMoisture(double canopyBaseHeight, double canopyTopHeight, NumericVector cohortFMC, NumericVector cohortLoading, NumericVector H, NumericVector CR);
double fuelbedLiveFuelMoisture(double fuelbedHeight, NumericVector cohortFMC, NumericVector cohortLoading, NumericVector H, NumericVector CR);
