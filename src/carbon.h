#include <Rcpp.h>

#ifndef CARBON_H
#define CARBON_H
#endif
using namespace Rcpp;

const double sucroseMolarWeight = 342.3; //g*mol-1
const double starchMolarWeight = 162.1406; //g*mol-1
const double leafCperDry = 0.3; //g C · g dry-1
const double rootCperDry = 0.4959;

NumericVector carbonCompartments(double SA, double LAI, double H, double Z, double N, double SLA, double WoodDensity, double WoodC);