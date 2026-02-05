#include <RcppArmadillo.h>
#include "medfate.h"
#include "control_c.h"
#include "modelInput_c.h"

#ifndef INNER_SPERRY_C_H
#define INNER_SPERRY_C_H

struct SperryNetwork {
  SperryWBParams sperryParams;

  std::vector<double> psisoil;
  std::vector<double> krhizomax;
  std::vector<double> nsoil;
  std::vector<double> alphasoil;
  std::vector<double> krootmax;
  double rootc;
  double rootd;
  double kstemmax;
  double stemc;
  double stemd;
  double kleafmax;
  double kleafapomax;
  double kleafsymp;
  double leafc;
  double leafd;
  double PLCstem;
  double PLCleaf;

};

void initSperryNetwork_inner_c(SperryNetwork& network,
                               int c,
                               const InternalWater& internalWater, 
                               const TranspirationParams& paramsTranspiration, 
                               const WaterStorageParams& paramsWaterStorage,
                               const std::vector<double>& VCroot_kmax, 
                               const std::vector<double>& VGrhizo_kmax,
                               const std::vector<double>& psiSoil, 
                               const std::vector<double>& VG_n, 
                               const std::vector<double>& VG_alpha,
                               const ControlParameters& control,
                               double sapFluidityDay = 1.0);
#endif