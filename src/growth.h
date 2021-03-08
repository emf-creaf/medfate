#include <Rcpp.h>

#ifndef GROWTH_H
#define GROWTH_H
#endif
using namespace Rcpp;

NumericVector carbonCompartments(double SA, double LAI, double H, double Z, 
                                 double N, double SLA, double WoodDens, double WoodC);

List growthDay(List x, CharacterVector date, double tmin, double tmax, 
               double rhmin, double rhmax, double rad, double wind, 
               double latitude, double elevation, double slope, double aspect,  
               double prec, double runon=0.0);
List growthDay1(List x, double tday, double pet, double prec, double er, double runon=0.0, 
                double rad = NA_REAL, double elevation = NA_REAL, bool verbose=false);
List growthDay2(List x, double tmin, double tmax, double tminPrev, double tmaxPrev, double tminNext, 
                double rhmin, double rhmax, double rad, double wind, 
                double latitude, double elevation, double slope, double aspect,
                double solarConstant, double delta, 
                double prec, double pet, double er, double runon=0.0, bool verbose = false);