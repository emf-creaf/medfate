#include <Rcpp.h>

#ifndef GROWTH_H
#define GROWTH_H
#endif
using namespace Rcpp;

List growthDay(List x, CharacterVector date, NumericVector meteovec, 
               double latitude, double elevation, double slope, double aspect,  
               double runon=0.0, bool modifyInput = true);