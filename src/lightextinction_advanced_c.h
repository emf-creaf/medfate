#include <RcppArmadillo.h>

#ifndef LIGHTEXTINCTION_ADVANCED_C_H
#define LIGHTEXTINCTION_ADVANCED_C_H

struct CohortSunlitShadeAbsorbedRadiation_RESULT {
  arma::mat I_sunlit;
  arma::mat I_shade;
  
  CohortSunlitShadeAbsorbedRadiation_RESULT(size_t ncanlayers, size_t numCohorts) : 
    I_sunlit(ncanlayers, numCohorts, arma::fill::zeros),
    I_shade(ncanlayers, numCohorts, arma::fill::zeros){}
};
Rcpp::List copyCohortSunlitShadeAbsorvedRadiationResult_c(CohortSunlitShadeAbsorbedRadiation_RESULT& res);

double leafAngleCDF_c(double leafAngle, double p, double q);
double directionalExtinctionCoefficient_c(double p, double q, double solarElevation);

void layerDirectIrradianceFraction_c(std::vector<double>& Ifraction,
                                     const arma::mat& LAIme, const arma::mat& LAImd, const arma::mat& LAImx, 
                                     const std::vector<double>& kb, const std::vector<double>& ClumpingIndex, 
                                     const std::vector<double>& alpha, const std::vector<double>& gamma, 
                                     double trunkExtinctionFraction = 0.1);
void layerDiffuseIrradianceFraction_c(arma::mat& Ifraction,
                                      const arma::mat& LAIme, const arma::mat& LAImd, const arma::mat& LAImx, 
                                      const arma::mat& K, const std::vector<double>& ClumpingIndex, const std::vector<double>& ZF,
                                      const std::vector<double>& alpha, const std::vector<double>& gamma, 
                                      double trunkExtinctionFraction = 0.1);

void cohortSunlitShadeAbsorbedRadiation_c(CohortSunlitShadeAbsorbedRadiation_RESULT& res,
                                          double Ib0, double Id0,
                                          const arma::mat& LAIme, const arma::mat& LAImd, const arma::mat& LAImx,
                                          const std::vector<double>& kb, const arma::mat& K, const std::vector<double>& ClumpingIndex, const std::vector<double>& ZF, 
                                          const std::vector<double>& alpha, const std::vector<double>& gamma, double trunkExtinctionFraction = 0.1);

void layerSunlitFraction_c(std::vector<double>& fSL,
                           const arma::mat& LAIme, const arma::mat& LAImd, 
                           const std::vector<double>& kb, const std::vector<double>& ClumpingIndex);
#endif
