#include "RcppArmadillo.h"
#include "communication_structures_c.h"
#include "transpiration_basic_c.h"
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
  // Stand-level (4 fields)
  StandBasicTranspiration_RESULT stand;
  
  // Plants data frame
  PlantsAdvancedTranspiration_RESULT plants;
  
  LeafAdvancedTranspiration_RESULT sunlit, shade;
  
  EnergyBalance_RESULT energy;
  
  // Extraction matrices
  arma::mat extraction;
  arma::mat extractionInst;
  std::vector<arma::mat> extractionPools;

  arma::mat rhizoPsi;
  
  CanopyTurbulence_RESULT canopyTurbulence;
  
  DirectDiffuseDay_RESULT directDiffuseDay;
  
  AdvancedTranspiration_RESULT(size_t numCohorts = 0, size_t nlayers = 0, size_t ncanlayers = 0, size_t ntimesteps = 0) : 
    plants(numCohorts), 
    sunlit(numCohorts),
    shade(numCohorts),
    energy(nlayers, ncanlayers, ntimesteps),
    extraction(numCohorts, nlayers, arma::fill::zeros),
    extractionInst(nlayers, ntimesteps, arma::fill::zeros),
    extractionPools(numCohorts),
    rhizoPsi(numCohorts, nlayers, arma::fill::zeros),
    canopyTurbulence(ncanlayers),
    directDiffuseDay(ntimesteps){
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
  
  AbsorbedSWR_COMM AbSWRcomm;
  std::vector<double> CohASWRF;
  std::vector<double> Tmax;
  std::vector<double> TmaxCoh;
  arma::mat RHOPCohDyn;
  CanopyTurbulenceModel_RESULT canopyTurbulenceModel;
  
  AdvancedTranspiration_COMM(size_t numCohorts = 0, size_t ncanlayers = 0, size_t nlayers= 0) : 
    AbSWRcomm(numCohorts, ncanlayers), 
    CohASWRF(numCohorts),
    Tmax(numCohorts),
    TmaxCoh(numCohorts), 
    RHOPCohDyn(numCohorts, nlayers),
    canopyTurbulenceModel(ncanlayers){}
};

void transpirationAdvanced_c(AdvancedTranspiration_RESULT& ATres, AdvancedTranspiration_COMM& ATcomm, ModelInput& x, 
                             const WeatherInputVector& meteovec,
                             const double latitude, double elevation, double slope, double aspect, 
                             const double solarConstant, const double delta,
                             const double canopyEvaporation, const double snowMelt, const double soilEvaporation, const double herbTranspiration, 
                             int stepFunctions);


#endif
