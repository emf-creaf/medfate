#include <RcppArmadillo.h>

#ifndef LIGHTEXTINCTION_BASIC_C_H
#define LIGHTEXTINCTION_BASIC_C_H

double availableLight_c(double h, 
                        const std::vector<double>& H, 
                        const std::vector<double>& LAI_expanded, 
                        const std::vector<double>& LAI_dead, 
                        const std::vector<double>& k, 
                        const std::vector<double>& CR);
void parcohortC_c(std::vector<double>& PARcohort,
                  const std::vector<double>& H, 
                  const std::vector<double>& LAI_expanded, 
                  const std::vector<double>& LAI_dead, 
                  const std::vector<double>& k, 
                  const std::vector<double>& CR);
void cohortAbsorbedSWRFraction_c(std::vector<double>& SWRfraction, 
                                 const arma::mat& LAIme, 
                                 const arma::mat& LAImd, 
                                 const std::vector<double>& kSWR);
#endif
