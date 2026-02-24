#include <RcppArmadillo.h>
#include "medfate.h"
#include "modelInput_c.h"

#ifndef DECOMPOSITION_C_H
#define DECOMPOSITION_C_H

const int DECOMPCOM_SURFACE_METABOLIC = 0;
const int DECOMPCOM_SOIL_METABOLIC = 1;
const int DECOMPCOM_SURFACE_ACTIVE = 2;
const int DECOMPCOM_SOIL_ACTIVE = 3;
const int DECOMPCOM_SURFACE_SLOW = 4;
const int DECOMPCOM_SOIL_SLOW = 5;
const int DECOMPCOM_SOIL_PASSIVE = 6;


//communication structures
struct SnagDecomposition_COMM {
  double transfer_surface_active, transfer_surface_slow, flux_respiration;
  SnagDecomposition_COMM() {
    transfer_surface_active = transfer_surface_slow = flux_respiration = 0.0;
  }
};
struct LitterDecomposition_COMM {
  double transfer_surface_active, transfer_surface_slow, transfer_soil_active, transfer_soil_slow, flux_respiration;
  LitterDecomposition_COMM() {
    transfer_surface_active = transfer_surface_slow = transfer_soil_active = transfer_soil_slow = flux_respiration = 0.0;
  }
};
// Structure with the following matrices/vectors:
//   \itemize{
//     \item{\code{xi}: Environmental scalar matrix.}
//     \item{\code{A}: Carbon transfer matrix.} 
//     \item{\code{pathf}: Fractional carbon flow from pool j to pool i.} 
//     \item{\code{respf}: Fractional respiration loss for carbon flow from pool j to pool i.} 
// }
struct Decomposition_COMM {
  SnagDecomposition_COMM sdo;
  LitterDecomposition_COMM ldo;
  std::vector<double> xi, K;
  arma::mat A, pathf, respf;
  double Kmix, K_s21;
  Decomposition_COMM(size_t npool) {
    Kmix = K_s21 = 0.0;
    K = std::vector<double>(npool, 0.0);
    xi = std::vector<double>(npool, 0.0);
    A = arma::mat(npool, npool);
    pathf = arma::mat(npool, npool);
    respf = arma::mat(npool, npool);
    for(int i = 0; i< (int) npool; i++) {
      for(int j = 0; j< (int) npool; j++) {
        A(i,j) = 0.0;
        pathf(i,j) = 0.0;
        respf(i,j) = 0.0;
      }
    }
  }
  Decomposition_COMM() {
    Decomposition_COMM(7);
  }
};
//low-level functions
double litterMetabolicFraction_c(double ligninPercent, double Nmass);
double annualLitterDecompositionRate_c(double AET, double lignin);
double pHEffect_c(double x, const std::string& pool);
double moistureEffect_c(double sand, double clay, double soilMoisture);
double temperatureEffect_c(double soilTemperature);

//increase litter pools
void addLeafTwigLitter_c(std::string& species_litter, double leaf_litter, double twig_litter,
                         InternalLitter& litter, 
                         LitterDecompositionParams& paramsLitterDecomposition,
                         InternalSOC& SOC);
void addSmallBranchLitter_c(std::string& species_litter, double smallbranch_litter, 
                            InternalLitter& litter);
void addLargeWoodLitter_c(std::string& species_litter, double largewood_litter, 
                          InternalLitter& litter);
void addCoarseRootLitter_c(std::string& species_litter, double coarsewood_litter, 
                           InternalLitter& litter);
void addFineRootLitter_c(std::string& species_litter, double fineroot_litter, 
                         InternalLitter& litter, 
                         LitterDecompositionParams& paramsLitterDecomposition,
                         InternalSOC& SOC);

//Update DAYCENT parameters
void updateBaseRates_c(Decomposition_COMM& DECcomm,
                       DecompositionAnnualBaseRates& baseAnnualRates, double annualTurnoverRate);
void updateDecompositionRateScalars_c(Decomposition_COMM& DECcomm, 
                                      double sand, double clay,
                                      double soilTemperature, double soilMoisture, double soilPH, 
                                      double soilO2, double cultfac);
void updateCarbonTransferMatrices_c(Decomposition_COMM& DECcomm, 
                                    double sand, double clay, double soilO2);


void DAYCENTsnagsInner_c(SnagDecomposition_COMM& sdo,
                         InternalSnags& snags, LitterDecompositionParams& paramsLitterDecomposition,
                         DecompositionAnnualBaseRates& baseAnnualRates,
                         double airTemperature, double airRelativeHumidity,
                         double tstep);

void DAYCENTlitterInner_c(LitterDecomposition_COMM& ldo, 
                          InternalLitter& litter, LitterDecompositionParams& paramsLitterDecomposition,
                          DecompositionAnnualBaseRates& baseAnnualRates,
                          double sand, double clay, double soilTemperature, double soilMoisture, double soilPH, 
                          double soilO2 = 1.0, double cultfac = 1.0,
                          double tstep = 1.0);

double DAYCENTInner_c(Decomposition_COMM& DECcomm,
                      InternalSnags& snags, InternalLitter& litter, InternalSOC& SOC,
                      LitterDecompositionParams& paramsLitterDecomposition,
                      DecompositionAnnualBaseRates& baseAnnualRates, double annualTurnoverRate,
                      double airTemperature, double airRelativeHumidity, 
                      double sand, double clay, double soilTemperature, double soilMoisture, double soilPH, 
                      double soilO2 = 1.0, double cultfac = 1.0,
                      double tstep = 1.0);
#endif
  