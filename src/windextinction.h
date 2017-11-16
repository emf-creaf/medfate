#include <Rcpp.h>

#ifndef WINDEXTINCTION_H
#define WINDEXTINCTION_H
#endif
using namespace Rcpp;

double windSpeedAtCanopyHeight(double wind20H, double canopyHeight);
double windSpeedAtHeightOverCanopy(double z, double wind20H, double canopyHeight);
double windSpeedMassmanExtinction(double z, double wind20H, double LAIc, double canopyHeight);
double unshelteredMidflameWindSpeed(double wind20H, double fuelBedHeight);
double shelteredMidflameWindSpeed(double wind20H, double crownFillProportion, double topCanopyHeight);
NumericVector windExtinctionProfile(NumericVector z, double wind20H, double LAIc, double canopyHeight);
NumericVector windExtinctionCohort(NumericVector H, NumericVector CR, double wind20H, double LAIc, double canopyHeight);
double aerodynamicResistance(double canopyHeight, double wind);