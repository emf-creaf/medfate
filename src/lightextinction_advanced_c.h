#include <RcppArmadillo.h>
#include "radiation_c.h"

#ifndef LIGHTEXTINCTION_ADVANCED_C_H
#define LIGHTEXTINCTION_ADVANCED_C_H


struct InstantaneousMultilayerAbsortion_RESULT {
  std::vector<arma::mat> PAR_SL;
  std::vector<arma::mat> PAR_SH;
  std::vector<arma::mat> SWR_SL;
  std::vector<arma::mat> SWR_SH;
  
  InstantaneousMultilayerAbsortion_RESULT(size_t numCohorts = 0, size_t ncanlayers = 0, size_t ntimesteps = 0) : 
    PAR_SL(ntimesteps),
    PAR_SH(ntimesteps),
    SWR_SL(ntimesteps),
    SWR_SH(ntimesteps) {
      for(size_t n=0;n<ntimesteps;n++) {
        PAR_SL[n] = arma::mat(ncanlayers, numCohorts);
        PAR_SH[n] = arma::mat(ncanlayers, numCohorts);
        SWR_SL[n] = arma::mat(ncanlayers, numCohorts);
        SWR_SH[n] = arma::mat(ncanlayers, numCohorts);
      }
   }
};
struct InstantaneousSunlitShadeAbsortion_RESULT {
  std::vector<std::vector<double>> PAR_SL;
  std::vector<std::vector<double>> PAR_SH;
  std::vector<std::vector<double>> SWR_SL;
  std::vector<std::vector<double>> SWR_SH;
  
  InstantaneousSunlitShadeAbsortion_RESULT(size_t numCohorts = 0, size_t ntimesteps = 0) : 
    PAR_SL(ntimesteps),
    PAR_SH(ntimesteps),
    SWR_SL(ntimesteps),
    SWR_SH(ntimesteps) {
    for(size_t n=0;n<ntimesteps;n++) {
      PAR_SL[n] = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
      PAR_SH[n] = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
      SWR_SL[n] = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
      SWR_SH[n] = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    }
  }
};
struct InstantaneousLightExtinctionAbsortion_RESULT {
  std::vector<std::vector<double>> fsunlit;
  InstantaneousMultilayerAbsortion_RESULT multilayer;
  InstantaneousSunlitShadeAbsortion_RESULT sunshade;
  std::vector<double> SWR_can;
  std::vector<double> SWR_soil;
  std::vector<double> gbf;
  std::vector<double> gdf;
  InstantaneousLightExtinctionAbsortion_RESULT(size_t numCohorts = 0, size_t ncanlayers = 0, size_t ntimesteps = 0) : 
    fsunlit(ntimesteps),
    multilayer(numCohorts, ncanlayers, ntimesteps),
    sunshade(numCohorts, ntimesteps){
    SWR_can = std::vector<double>(ntimesteps, medfate::NA_DOUBLE);
    SWR_soil = std::vector<double>(ntimesteps, medfate::NA_DOUBLE);
    gbf = std::vector<double>(ntimesteps, medfate::NA_DOUBLE);
    gdf = std::vector<double>(ntimesteps, medfate::NA_DOUBLE);
    for(size_t n=0;n<ntimesteps;n++) {
      fsunlit[n] = std::vector<double>(ncanlayers, medfate::NA_DOUBLE);
    }
  }
};
Rcpp::List copyInstantaneousLightExtinctionAbsortionResult_c(const InstantaneousLightExtinctionAbsortion_RESULT& res);

struct LongWaveRadiationLayer_RESULT {
  std::vector<double> Lup;
  std::vector<double> Ldown;
  std::vector<double> Lnet;
  std::vector<double> tau;
  std::vector<double> sumTauComp;
  
  LongWaveRadiationLayer_RESULT(size_t ncanlayers) {
    Lup = std::vector<double>(ncanlayers, medfate::NA_DOUBLE);
    Ldown = std::vector<double>(ncanlayers, medfate::NA_DOUBLE);
    Lnet = std::vector<double>(ncanlayers, medfate::NA_DOUBLE);
    tau = std::vector<double>(ncanlayers, medfate::NA_DOUBLE);
    sumTauComp = std::vector<double>(ncanlayers, medfate::NA_DOUBLE);
  }
};
struct LongWaveRadiation_RESULT {
  LongWaveRadiationLayer_RESULT LWR_layer;
  double Ldown_ground;
  double Lup_ground;
  double Lnet_ground;
  double Ldown_canopy;
  double Lup_canopy;
  double Lnet_canopy;
  arma::mat Lnet_cohort_layer;
  
  LongWaveRadiation_RESULT(size_t ncanlayers, size_t numCohorts) : LWR_layer(ncanlayers), Lnet_cohort_layer(ncanlayers, numCohorts){}
};
Rcpp::List copyLongWaveRadiationResult_c(const LongWaveRadiation_RESULT& res);

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

void cohortSunlitShadeAbsorbedRadiation_c(arma::mat& I_sunlit, arma::mat& I_shade,
                                          double Ib0, double Id0,
                                          const arma::mat& LAIme, const arma::mat& LAImd, const arma::mat& LAImx,
                                          const std::vector<double>& kb, const arma::mat& K, const std::vector<double>& ClumpingIndex, const std::vector<double>& ZF, 
                                          const std::vector<double>& alpha, const std::vector<double>& gamma, double trunkExtinctionFraction = 0.1);

void layerSunlitFraction_c(std::vector<double>& fSL,
                           const arma::mat& LAIme, const arma::mat& LAImd, 
                           const std::vector<double>& kb, const std::vector<double>& ClumpingIndex);

void instantaneousLightExtinctionAbsortion_c(InstantaneousLightExtinctionAbsortion_RESULT& res,
                                             const arma::mat& LAIme, const arma::mat& LAImd, const arma::mat& LAImx, 
                                             const std::vector<double>& p, const std::vector<double>& q, const std::vector<double>& ClumpingIndex, 
                                             const std::vector<double>& alphaSWR, const std::vector<double>& gammaSWR,
                                             const DirectDiffuseDay_RESULT& ddd, int ntimesteps, double trunkExtinctionFraction);

void longwaveRadiationSHAW_inner_c(LongWaveRadiation_RESULT& res, 
                                   const arma::mat& LAIme, const arma::mat& LAImd, const arma::mat& LAImx, 
                                   double LWRatm, double Tsoil, const std::vector<double>& Tair, double trunkExtinctionFraction);
#endif
