#include <Rcpp.h>

#ifndef FIREBEHAVIOUR_H
#define FIREBEHAVIOUR_H
#endif
using namespace Rcpp;

double criticalFirelineIntensity(double CBH, double M);

List FCCSbehaviour(DataFrame FCCSpropsSI,
                   NumericVector MliveSI = NumericVector::create(90, 90, 60), 
                   NumericVector MdeadSI = NumericVector::create(6, 6, 6, 6, 6), 
                   double slope = 0.0, double windSpeedSI = 11.0);
