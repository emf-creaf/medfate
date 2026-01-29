#include "RcppArmadillo.h"

#ifndef FORESTUTILS_C_H
#define FORESTUTILS_C_H
using namespace Rcpp;

double leafAreaProportion_c(double z1, double z2, double zmin, double zmax);

void updateLAIdistributionVectors_c(arma::mat& LAIdist, 
                                    const std::vector<double>& z, 
                                    const std::vector<double>& LAI, 
                                    const std::vector<double>& H, 
                                    const std::vector<double>& CR);

#endif
