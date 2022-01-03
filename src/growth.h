#include <Rcpp.h>

#ifndef GROWTH_H
#define GROWTH_H
#endif
using namespace Rcpp;

List growthDay(List x, CharacterVector date, double tmin, double tmax, 
               double rhmin, double rhmax, double rad, double wind, 
               double latitude, double elevation, double slope, double aspect,  
               double prec, double runon=0.0);
List growthDay1(List x, NumericVector meteovec, 
                double elevation = NA_REAL, 
                double runon=0.0, bool verbose = false);
List growthDay2(List x, NumericVector meteovec, 
                double latitude, double elevation, double slope, double aspect,
                double solarConstant, double delta, 
                double runon=0.0, bool verbose = false);