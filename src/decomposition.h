#include <Rcpp.h>

#ifndef DECOMPOSITION_H
#define DECOMPOSITION_H
#endif
using namespace Rcpp;

double litterMetabolicFraction(double ligninPercent, double Nmass);
double annualLitterDecompositionRate(double AET, double lignin);

void addLeafLitter(String species_litter, double leaf_litter, 
                   DataFrame structuralLitter, 
                   DataFrame paramsLitterDecomposition,
                   NumericVector internalSOM);

void addFineRootLitter(String species_litter, double fineroot_litter, 
                       DataFrame structuralLitter, 
                       DataFrame paramsLitterDecomposition,
                       NumericVector internalSOM);
  
void addSmallBranchLitter(String species_litter, double smallbranch_litter, 
                          DataFrame structuralLitter);

double DAYCENTInner(List commDecomp,
                    DataFrame structuralLitter, NumericVector CENTURYPools,
                    DataFrame paramsLitterDecomposition,
                    NumericVector baseAnnualRates, double annualTurnoverRate,
                    double sand, double clay, double soilTemperature, double soilMoisture, double soilPH, 
                    double soilO2 = 1.0, double cultfac = 1.0,
                    double tstep = 1.0);
