#include <Rcpp.h>

#ifndef GROWTH_H
#define GROWTH_H
#endif
using namespace Rcpp;

NumericVector carbonCompartments(double SA, double LAI, double H, double Z, 
                                 double N, double SLA, double WoodDens, double WoodC);