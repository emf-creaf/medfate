#include "RcppArmadillo.h"
#include "lowlevel_structures_c.h"
#include "transpiration_basic_c.h"
#include "lightextinction_advanced_c.h"
#include "soil_thermodynamics_c.h"
#include "modelInput_c.h"
#include "radiation_c.h"
#include "windKatul_c.h"

#ifndef TRANSPIRATION_ADVANCED_C_H
#define TRANSPIRATION_ADVANCED_C_H


// ----------------------------------------------------------------------------
// Advanced Transpiration Output Structures
// ----------------------------------------------------------------------------
struct PlantsAdvancedTranspiration_RESULT {
  std::vector<double> LAI;
  std::vector<double> LAIlive;
  std::vector<double> FPAR;
  std::vector<double> Extraction;
  std::vector<double> Transpiration;
  std::vector<double> GrossPhotosynthesis;
  std::vector<double> NetPhotosynthesis;
  std::vector<double> RootPsi;
  std::vector<double> StemPsi;
  std::vector<double> StemPLC;
  std::vector<double> LeafPLC;
  std::vector<double> LeafPsiMin;
  std::vector<double> LeafPsiMax;
  std::vector<double> dEdP;
  std::vector<double> DDS;
  std::vector<double> StemRWC;
  std::vector<double> LeafRWC;
  std::vector<double> LFMC;
  std::vector<double> WaterBalance;
  
  PlantsAdvancedTranspiration_RESULT(size_t numCohorts = 0) {
    LAI = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    LAIlive = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    FPAR = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    Extraction = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    Transpiration = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    GrossPhotosynthesis = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    NetPhotosynthesis = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    RootPsi = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    StemPsi = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    StemPLC = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    LeafPLC = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    LeafPsiMin = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    LeafPsiMax = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    dEdP = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    DDS = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    StemRWC = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    LeafRWC = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    LFMC = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    WaterBalance = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
  }
};
Rcpp::DataFrame copyPlantAdvancedTranspirationResult_c(const PlantsAdvancedTranspiration_RESULT& plants, ModelInput& x);


struct PlantsAdvancedTranspirationInst_RESULT {
  arma::mat E, Ag, An, dEdP, RootPsi, StemPsi, LeafPsi;
  arma::mat StemSympPsi, LeafSympPsi, StemPLC, LeafPLC, StemRWC, LeafRWC, StemSympRWC, LeafSympRWC;
  arma::mat PWB;
  
  PlantsAdvancedTranspirationInst_RESULT(size_t numCohorts, size_t ntimesteps) : 
    E(numCohorts, ntimesteps),
    Ag(numCohorts, ntimesteps),
    An(numCohorts, ntimesteps),
    dEdP(numCohorts, ntimesteps),
    RootPsi(numCohorts, ntimesteps),
    StemPsi(numCohorts, ntimesteps),
    LeafPsi(numCohorts, ntimesteps),
    StemSympPsi(numCohorts, ntimesteps),
    LeafSympPsi(numCohorts, ntimesteps),
    StemPLC(numCohorts, ntimesteps),
    LeafPLC(numCohorts, ntimesteps),
    StemRWC(numCohorts, ntimesteps),
    LeafRWC(numCohorts, ntimesteps),
    StemSympRWC(numCohorts, ntimesteps),
    LeafSympRWC(numCohorts, ntimesteps),
    PWB(numCohorts, ntimesteps) {
    
  }
};
Rcpp::List copyPlantAdvancedTranspirationInstResult_c(const PlantsAdvancedTranspirationInst_RESULT& plants, ModelInput& x);

struct LeafAdvancedTranspiration_RESULT {
  std::vector<double> LeafPsiMin;
  std::vector<double> LeafPsiMax;
  std::vector<double> GSWMin;
  std::vector<double> GSWMax;
  std::vector<double> TempMin;
  std::vector<double> TempMax;
  
  LeafAdvancedTranspiration_RESULT(size_t numCohorts = 0) {
    LeafPsiMin = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    LeafPsiMax = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    GSWMin = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    GSWMax = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    TempMin = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    TempMax = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
  }
};
Rcpp::DataFrame copyLeafAdvancedTranspirationResult_c(const LeafAdvancedTranspiration_RESULT& leaf, ModelInput& x);


struct LeafAdvancedTranspirationInst_RESULT {
  arma::mat LAI;
  arma::mat Vmax298;
  arma::mat Jmax298;
  arma::mat Abs_SWR;
  arma::mat Abs_PAR;
  arma::mat Net_LWR;
  arma::mat Ag;
  arma::mat An;
  arma::mat Ci;
  arma::mat E;
  arma::mat Gsw;
  arma::mat VPD;
  arma::mat Temp;
  arma::mat Psi;
  
  LeafAdvancedTranspirationInst_RESULT(size_t numCohorts, size_t ntimesteps) : 
    LAI(numCohorts, ntimesteps),
    Vmax298(numCohorts, ntimesteps),
    Jmax298(numCohorts, ntimesteps),
    Abs_SWR(numCohorts, ntimesteps),
    Abs_PAR(numCohorts, ntimesteps),
    Net_LWR(numCohorts, ntimesteps),
    Ag(numCohorts, ntimesteps),
    An(numCohorts, ntimesteps),
    Ci(numCohorts, ntimesteps),
    E(numCohorts, ntimesteps),
    Gsw(numCohorts, ntimesteps),
    VPD(numCohorts, ntimesteps),
    Temp(numCohorts, ntimesteps),
    Psi(numCohorts, ntimesteps)
    { }
};
Rcpp::List copyLeafAdvancedTranspirationInstResult_c(const LeafAdvancedTranspirationInst_RESULT& leaf_inst, ModelInput& x);

struct EnergyBalance_RESULT {
  std::vector<double> SolarHour;
  std::vector<double> Tatm;
  std::vector<double> Tcan;
  arma::mat SoilTemperature;
  std::vector<double> SWRcan;
  std::vector<double> LWRcan;
  std::vector<double> LEVcan;
  std::vector<double> LEFsnow;
  std::vector<double> Hcan;
  std::vector<double> Ebalcan;
  std::vector<double> Hcansoil;
  std::vector<double> LEVsoil;
  std::vector<double> SWRsoil;
  std::vector<double> LWRsoil;
  std::vector<double> Ebalsoil;
  arma::mat TemperatureLayers;
  arma::mat VaporPressureLayers;
  
  EnergyBalance_RESULT(size_t nlayers = 0, size_t ncanlayers = 0, size_t ntimesteps = 0) : 
    SoilTemperature(ntimesteps, nlayers, arma::fill::zeros),
    TemperatureLayers(ntimesteps, ncanlayers, arma::fill::zeros),
    VaporPressureLayers(ntimesteps, ncanlayers, arma::fill::zeros)
    {
    SolarHour = std::vector<double>(ntimesteps, medfate::NA_DOUBLE);
    Tatm = std::vector<double>(ntimesteps, medfate::NA_DOUBLE);
    Tcan = std::vector<double>(ntimesteps, medfate::NA_DOUBLE);
    SWRcan = std::vector<double>(ntimesteps, medfate::NA_DOUBLE);
    LWRcan = std::vector<double>(ntimesteps, medfate::NA_DOUBLE);
    LEVcan = std::vector<double>(ntimesteps, medfate::NA_DOUBLE);
    LEFsnow = std::vector<double>(ntimesteps, medfate::NA_DOUBLE);
    Hcan = std::vector<double>(ntimesteps, medfate::NA_DOUBLE);
    Ebalcan = std::vector<double>(ntimesteps, medfate::NA_DOUBLE);
    Hcansoil = std::vector<double>(ntimesteps, medfate::NA_DOUBLE);
    LEVsoil = std::vector<double>(ntimesteps, medfate::NA_DOUBLE);
    SWRsoil = std::vector<double>(ntimesteps, medfate::NA_DOUBLE);
    LWRsoil = std::vector<double>(ntimesteps, medfate::NA_DOUBLE);
    Ebalsoil = std::vector<double>(ntimesteps, medfate::NA_DOUBLE);
};
};
Rcpp::List copyEnergyBalanceResult_c(const EnergyBalance_RESULT& EBres, ModelInput& x);


struct AdvancedTranspiration_RESULT {
  size_t numCohorts;
  size_t nlayers;
  size_t ncanlayers;
  size_t ntimesteps;
  
  // Stand-level (4 fields)
  StandBasicTranspiration_RESULT stand;
  
  // Plants data frame
  PlantsAdvancedTranspiration_RESULT plants;
  
  PlantsAdvancedTranspirationInst_RESULT plants_inst;
  
  LeafAdvancedTranspiration_RESULT sunlit, shade;
  
  LeafAdvancedTranspirationInst_RESULT sunlit_inst, shade_inst;
  
  EnergyBalance_RESULT energy;
  
  // Extraction matrices
  arma::mat extraction;
  arma::mat extractionInst;
  std::vector<arma::mat> extractionPools;

  arma::mat rhizoPsi;
  
  CanopyTurbulence_RESULT canopyTurbulence;
  
  DirectDiffuseDay_RESULT directDiffuseDay;

  InstantaneousLightExtinctionAbsortion_RESULT lightExtinctionAbsortion;
  
  std::vector<LongWaveRadiation_RESULT> lwrExtinction;
  
  AdvancedTranspiration_RESULT(size_t numCohorts = 0, size_t nlayers = 0, size_t ncanlayers = 0, size_t ntimesteps = 0) : 
    numCohorts(numCohorts),
    nlayers(nlayers),
    ncanlayers(ncanlayers),
    ntimesteps(ntimesteps),
    plants(numCohorts), 
    plants_inst(numCohorts, ntimesteps),
    sunlit(numCohorts),
    shade(numCohorts),
    sunlit_inst(numCohorts, ntimesteps),
    shade_inst(numCohorts, ntimesteps),
    energy(nlayers, ncanlayers, ntimesteps),
    extraction(numCohorts, nlayers, arma::fill::zeros),
    extractionInst(nlayers, ntimesteps, arma::fill::zeros),
    extractionPools(numCohorts),
    rhizoPsi(numCohorts, nlayers, arma::fill::zeros),
    canopyTurbulence(ncanlayers),
    directDiffuseDay(ntimesteps),
    lightExtinctionAbsortion(numCohorts, ncanlayers, ntimesteps),
    lwrExtinction(ntimesteps, LongWaveRadiation_RESULT(ncanlayers, numCohorts)) {
    for(size_t c = 0; c < numCohorts; c++) {
      extractionPools[c] = arma::mat(numCohorts, nlayers, arma::fill::zeros);
    }
  }
};
Rcpp::List copyAdvancedTranspirationResult_c(const AdvancedTranspiration_RESULT& BTres, ModelInput& x);



// ----------------------------------------------------------------------------
// Advanced Transpiration Communication Structures
// ----------------------------------------------------------------------------
struct AdvancedTranspiration_COMM {
  
  arma::mat RHOPCohDyn;
  CanopyTurbulenceModel_RESULT canopyTurbulenceModel;
  SoilEnergyBalance_COMM SEBcomm;
  
  AdvancedTranspiration_COMM(size_t numCohorts = 0, size_t nlayers = 0, size_t ncanlayers = 0, size_t ntimesteps = 0) : 
    RHOPCohDyn(numCohorts, nlayers),
    canopyTurbulenceModel(ncanlayers),
    SEBcomm(nlayers) {}
};


struct InnerTranspirationInput_COMM {
  double Patm;
  std::vector<double> zWind;
  double f_dry;
  std::vector<int> iLayerCohort, iLayerSunlit, iLayerShade;
  std::vector<int> iPMSunlit, iPMShade;
  std::vector<int> nlayerscon;
  arma::Mat<uint8_t> layerConnected;
  std::vector<arma::Mat<uint8_t>> layerConnectedPools;
  std::vector<double> psiSoil;
  arma::mat psiSoilM;
  arma::mat KunsatM;
  
  InnerTranspirationInput_COMM(size_t numCohorts = 0, size_t nlayers = 0, size_t ncanlayers = 0) {
    zWind = std::vector<double>(ncanlayers, medfate::NA_DOUBLE);
    iLayerCohort = std::vector<int>(numCohorts, medfate::NA_INTEGER);
    iLayerSunlit = std::vector<int>(numCohorts, medfate::NA_INTEGER);
    iLayerShade = std::vector<int>(numCohorts, medfate::NA_INTEGER);
    iPMSunlit = std::vector<int>(numCohorts, 0); //Initial values set to closed stomata
    iPMShade = std::vector<int>(numCohorts, 0);
    nlayerscon = std::vector<int>(numCohorts, 0);
    layerConnected = arma::Mat<uint8_t>(numCohorts, nlayers);
    layerConnectedPools = std::vector<arma::Mat<uint8_t>>(numCohorts);
    for(size_t c = 0; c<numCohorts; c++) {
      layerConnectedPools[c] = arma::Mat<uint8_t>(numCohorts, nlayers);
    }
    psiSoil = std::vector<double>(nlayers, medfate::NA_DOUBLE);
    psiSoilM = arma::mat(numCohorts, nlayers);
    KunsatM = arma::mat(numCohorts, nlayers);
  }
};


void transpirationAdvanced_c(AdvancedTranspiration_RESULT& ATres, AdvancedTranspiration_COMM& ATcomm, ModelInput& x, 
                             const WeatherInputVector& meteovec,
                             const double latitude, double elevation, double slope, double aspect, 
                             const double solarConstant, const double delta,
                             const double canopyEvaporation, const double snowMelt, const double soilEvaporation, const double herbTranspiration, 
                             int stepFunctions);


#endif
