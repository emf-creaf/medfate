#include <RcppArmadillo.h>
#include "medfate.h"
#include "control_c.h"
#include "photosynthesis_c.h"
#include "modelInput_c.h"
#include "transpiration_advanced_c.h"

#ifndef INNER_SPERRY_C_H
#define INNER_SPERRY_C_H


struct NetworkSteadyState {
  double E;
  std::vector<double> ERoot;
  std::vector<double> ERhizo;
  double psiLeaf;
  double psiStem;
  double psiRootCrown;
  std::vector<double> psiRhizo;
  std::vector<double> x;
  NetworkSteadyState(int nlayers, double E) {
    E = E;
    ERhizo = std::vector<double>(nlayers, medfate::NA_DOUBLE);
    ERoot = std::vector<double>(nlayers, medfate::NA_DOUBLE);
    psiRhizo = std::vector<double>(nlayers, medfate::NA_DOUBLE);
    x = std::vector<double>(nlayers+1, medfate::NA_DOUBLE);
  }
};
struct SupplyFunction {
  std::vector<double> E;
  std::vector<double> dEdP;
  arma::mat ERhizo;
  arma::mat psiRhizo;
  std::vector<double> psiRoot;
  std::vector<double> psiStem;
  std::vector<double> psiLeaf;
  SupplyFunction(int nlayers, int maxNsteps) {
    E = std::vector<double>(maxNsteps, medfate::NA_DOUBLE);
    dEdP = std::vector<double>(maxNsteps, medfate::NA_DOUBLE);
    ERhizo = arma::mat(maxNsteps,nlayers);
    psiRhizo = arma::mat(maxNsteps,nlayers);
    psiRoot = std::vector<double>(maxNsteps, medfate::NA_DOUBLE);
    psiStem = std::vector<double>(maxNsteps, medfate::NA_DOUBLE);
    psiLeaf = std::vector<double>(maxNsteps, medfate::NA_DOUBLE);
  }
  SupplyFunction() {
    SupplyFunction(0,0);
  }
};

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
  
  SupplyFunction supply;
};

struct ProfitMaximization {
  PhotoFunction photosynthesisFunction;
  double Profit;
  int iMaxProfit;
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

void fillSupplyFunctionNetwork_c(SperryNetwork&  hydraulicNetwork,
                                 double minFlow = 0.0, double pCrit = 0.001);

void innerSperry_c(ModelInput& x,
                   SperryNetwork* networks, 
                   InnerTranspirationInput_COMM& input, 
                   AdvancedTranspiration_RESULT& output, int n, double tstep, 
                   int stepFunctions = -1);
#endif