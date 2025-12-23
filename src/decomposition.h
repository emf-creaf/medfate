#include <Rcpp.h>

#ifndef DECOMPOSITION_H
#define DECOMPOSITION_H
#endif
using namespace Rcpp;

double litterMetabolicFraction(double ligninPercent, double Nmass);
double annualLitterDecompositionRate(double AET, double lignin);

void addLeafTwigLitter(String species_litter, double leaf_litter, double twig_litter,
                       DataFrame litter, 
                       DataFrame paramsLitterDecomposition,
                       NumericVector internalSOC);

void addSmallBranchLitter(String species_litter, double smallbranch_litter, 
                          DataFrame litter);

void addLargeWoodLitter(String species_litter, double largewood_litter, 
                        DataFrame litter);

void addCoarseRootLitter(String species_litter, double coarsewood_litter, 
                         DataFrame litter);

void addFineRootLitter(String species_litter, double fineroot_litter, 
                       DataFrame litter, 
                       DataFrame paramsLitterDecomposition,
                       NumericVector SOC);
  
double DAYCENTInner(List commDecomp,
                    DataFrame snags, DataFrame litter, NumericVector SOC,
                    DataFrame paramsLitterDecomposition,
                    NumericVector baseAnnualRates, double annualTurnoverRate,
                    double airTemperature, double airRelativeHumidity, 
                    double sand, double clay, double soilTemperature, double soilMoisture, double soilPH, 
                    double soilO2 = 1.0, double cultfac = 1.0,
                    double tstep = 1.0);
  