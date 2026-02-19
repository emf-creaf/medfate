#include <RcppArmadillo.h>

#ifndef DECOMPOSITION_H
#define DECOMPOSITION_H

using namespace Rcpp;

const int LITDECOMPCOM_TRANSFER_SURFACE_ACTIVE = 0;
const int LITDECOMPCOM_TRANSFER_SURFACE_SLOW = 1;
const int LITDECOMPCOM_TRANSFER_SOIL_ACTIVE = 2;
const int LITDECOMPCOM_TRANSFER_SOIL_SLOW = 3;
const int LITDECOMPCOM_FLUX_RESPIRATION = 4;

const int SNAGDECOMPCOM_TRANSFER_SURFACE_ACTIVE = 0;
const int SNAGDECOMPCOM_TRANSFER_SURFACE_SLOW = 1;
const int SNAGDECOMPCOM_FLUX_RESPIRATION = 2;

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
                    double soilO2, double cultfac,
                    double tstep);

#endif