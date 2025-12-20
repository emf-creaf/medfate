#include <Rcpp.h>

#ifndef DECOMPOSITION_H
#define DECOMPOSITION_H
#endif
using namespace Rcpp;

double litterMetabolicFraction(double ligninPercent, double Nmass);
double annualLitterDecompositionRate(double AET, double lignin);

double DAYCENTInner(List commDecomp,
                    DataFrame structuralLitter, NumericVector CENTURYPools,
                    DataFrame paramsDecomposition,
                    NumericVector baseAnnualRates, double annualTurnoverRate,
                    double sand, double clay, double soilTemperature, double soilMoisture, double soilPH, 
                    double soilO2 = 1.0, double cultfac = 1.0,
                    double tstep = 1.0);
