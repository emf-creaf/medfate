#include <RcppArmadillo.h>

#ifndef FIREBEHAVIOUR_H
#define FIREBEHAVIOUR_H
using namespace Rcpp;

List FCCSbehaviour(DataFrame FCCSpropsSI,
                   NumericVector MliveSI, 
                   NumericVector MdeadSI, 
                   double slope, 
                   double windSpeedSI);

#endif
