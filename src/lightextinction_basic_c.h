#include <RcppArmadillo.h>

#ifndef LIGHTEXTINCTION_BASIC_C_H
#define LIGHTEXTINCTION_BASIC_C_H

struct AbsorbedSWR_COMM {
  std::vector<double> fi; 
  std::vector<double> rem;
  arma::mat fij; 
  AbsorbedSWR_COMM(size_t numCohorts = 0, size_t ncanlayers = 0) : fi(ncanlayers), rem(ncanlayers), fij(ncanlayers, numCohorts) {}
};
  
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
                                 AbsorbedSWR_COMM& AbSWRcomm,
                                 const arma::mat& LAIme, 
                                 const arma::mat& LAImd, 
                                 const std::vector<double>& kSWR);
#endif
