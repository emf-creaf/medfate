#include <Rcpp.h>

#ifndef GROWTH_H
#define GROWTH_H
#endif
using namespace Rcpp;

List growthDay(List x, CharacterVector date, double tmin, double tmax, 
               double rhmin, double rhmax, double rad, double wind, 
               double latitude, double elevation, double slope, double aspect,  
               double prec, double runon=0.0, bool modifyInput = true);