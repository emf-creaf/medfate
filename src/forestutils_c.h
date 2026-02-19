#include "RcppArmadillo.h"

#ifndef FORESTUTILS_C_H
#define FORESTUTILS_C_H
using namespace Rcpp;

std::vector<std::string> cohortType_c(const std::vector<std::string> IDs);

double leafAreaProportion_c(double z1, double z2, double zmin, double zmax);

void updateLAIdistributionVectors_c(arma::mat& LAIdist, 
                                    const std::vector<double>& z, 
                                    const std::vector<double>& LAI, 
                                    const std::vector<double>& H, 
                                    const std::vector<double>& CR);

std::vector<double> largerTreeBasalArea_c(const std::vector<double>& N, const std::vector<double>& dbh, 
                                          double self_include_prop = 0.5);

std::vector<double> treeCrownRatioAllometric_c(std::vector<double>& N, std::vector<double>& dbh, std::vector<double>& H, 
                                               const std::vector<double>& Acw, const std::vector<double>& Bcw,
                                               const std::vector<double>& Acr, const std::vector<double>& B1cr, 
                                               const std::vector<double>& B2cr, 
                                               const std::vector<double>& B3cr,
                                               const std::vector<double>& C1cr, const std::vector<double>& C2cr);

#endif
