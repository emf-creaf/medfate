#include "lowlevel_structures_c.h"
#include "spwb_day_c.h"
#include "carbon_c.h"
#include "decomposition_c.h"

#ifndef GROWTH_DAY_C_H
#define GROWTH_DAY_C_H

struct LabileCarbonBalance_RESULT {
  std::vector<double> GrossPhotosynthesis;
  std::vector<double> MaintenanceRespiration;
  std::vector<double> GrowthCosts;
  std::vector<double> RootExudation;
  std::vector<double> LabileCarbonBalance;
  std::vector<double> SugarLeaf;
  std::vector<double> StarchLeaf;
  std::vector<double> SugarSapwood; 
  std::vector<double> StarchSapwood;
  std::vector<double> SugarTransport;
  LabileCarbonBalance_RESULT(size_t numCohorts = 0) {
    GrossPhotosynthesis = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    MaintenanceRespiration = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    GrowthCosts = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    RootExudation = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    LabileCarbonBalance = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    SugarLeaf = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    StarchLeaf = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    SugarSapwood = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    StarchSapwood = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    SugarTransport = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
  }
};
Rcpp::DataFrame copyLabileCarbonBalanceResult_c(const LabileCarbonBalance_RESULT& LCBres, ModelInput& x);

struct LabileCarbonBalanceInst_RESULT {
  arma::mat GrossPhotosynthesis;
  arma::mat MaintenanceRespiration;
  arma::mat GrowthCosts;
  arma::mat RootExudation;
  arma::mat LabileCarbonBalance;
  arma::mat SugarLeaf;
  arma::mat StarchLeaf;
  arma::mat SugarSapwood; 
  arma::mat StarchSapwood;
  arma::mat SugarTransport;
  LabileCarbonBalanceInst_RESULT(size_t numCohorts = 0, size_t ntimesteps = 0) {
    GrossPhotosynthesis = arma::mat(numCohorts, ntimesteps);
    MaintenanceRespiration = arma::mat(numCohorts, ntimesteps);
    GrowthCosts = arma::mat(numCohorts, ntimesteps);
    RootExudation = arma::mat(numCohorts, ntimesteps);
    LabileCarbonBalance = arma::mat(numCohorts, ntimesteps);
    SugarLeaf = arma::mat(numCohorts, ntimesteps);
    StarchLeaf = arma::mat(numCohorts, ntimesteps);
    SugarSapwood = arma::mat(numCohorts, ntimesteps);
    StarchSapwood = arma::mat(numCohorts, ntimesteps);
    SugarTransport = arma::mat(numCohorts, ntimesteps);
  }
};
Rcpp::DataFrame copyLabileCarbonBalanceInstResult_c(const LabileCarbonBalanceInst_RESULT& LCBres, ModelInput& x);


struct PlantBiomassBalance_RESULT {
  std::vector<double> InitialDensity;
  std::vector<double> InitialSapwoodBiomass;
  std::vector<double> InitialStructuralBiomass;
  std::vector<double> StructuralBiomassBalance; 
  std::vector<double> StructuralBiomassChange;
  std::vector<double> InitialLabileBiomass;
  std::vector<double> LabileBiomassBalance;
  std::vector<double> LabileBiomassChange;
  std::vector<double> InitialLivingPlantBiomass;
  std::vector<double> InitialPlantBiomass;
  std::vector<double> PlantBiomassBalance;
  std::vector<double> PlantBiomassChange;
  std::vector<double> MortalityBiomassLoss;
  std::vector<double> InitialCohortBiomass;
  std::vector<double> CohortBiomassBalance;
  std::vector<double> CohortBiomassChange;
  
  PlantBiomassBalance_RESULT(size_t numCohorts = 0) {
    InitialDensity = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    InitialSapwoodBiomass = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    InitialStructuralBiomass = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    StructuralBiomassBalance = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    StructuralBiomassChange = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    InitialLabileBiomass = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    LabileBiomassBalance = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    LabileBiomassChange = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    InitialLivingPlantBiomass = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    InitialPlantBiomass = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    PlantBiomassBalance = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    PlantBiomassChange = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    MortalityBiomassLoss = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    InitialCohortBiomass = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    CohortBiomassBalance = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    CohortBiomassChange = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
  }
};
Rcpp::DataFrame copyPlantBiomassBalanceResult_c(const PlantBiomassBalance_RESULT& PBBres, ModelInput& x);

struct PlantStructure_RESULT {
  std::vector<double> LeafBiomass;
  std::vector<double> SapwoodBiomass;
  std::vector<double> FineRootBiomass;
  std::vector<double> LeafArea;
  std::vector<double> SapwoodArea;
  std::vector<double> FineRootArea;
  std::vector<double> HuberValue;
  std::vector<double> RootAreaLeafArea;
  std::vector<double> DBH;
  std::vector<double> Height;
  PlantStructure_RESULT(size_t numCohorts = 0) {
    LeafBiomass = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    SapwoodBiomass = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    FineRootBiomass = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    LeafArea = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    SapwoodArea = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    FineRootArea = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    HuberValue = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    RootAreaLeafArea = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    DBH = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    Height = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
  }
};
Rcpp::DataFrame copyPlantStructureResult_c(const PlantStructure_RESULT& PSres, ModelInput& x);

struct GrowthMortality_RESULT {
  std::vector<double> SAgrowth;
  std::vector<double> LAgrowth;
  std::vector<double> FRAgrowth;
  std::vector<double> StarvationRate;
  std::vector<double> DessicationRate;
  std::vector<double> MortalityRate;  
  GrowthMortality_RESULT(size_t numCohorts = 0) {
    SAgrowth = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    LAgrowth = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    FRAgrowth = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    StarvationRate = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    DessicationRate = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    MortalityRate = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
  }
};
Rcpp::DataFrame copyGrowthMortalityResult_c(const GrowthMortality_RESULT& GMres, ModelInput& x);

struct GROWTH_RESULT {
  
  LabileCarbonBalance_RESULT LCBres;
  PlantBiomassBalance_RESULT PBBres;
  PlantStructure_RESULT PSres;
  GrowthMortality_RESULT GMres;

  StandCB_RESULT standCB;  
  
  GROWTH_RESULT(size_t numCohorts) : 
    LCBres(numCohorts),
    PBBres(numCohorts),
    PSres(numCohorts),
    GMres(numCohorts){}
  virtual ~GROWTH_RESULT() = default;
};
Rcpp::List copyGROWTHResult_c(GROWTH_RESULT& GROWTHres, ModelInput& x);

struct BasicGROWTH_RESULT : GROWTH_RESULT {
  BasicSPWB_RESULT BSPWBres;
  
  BasicGROWTH_RESULT(BasicSPWB_RESULT SPWBres, size_t numCohorts) : 
    GROWTH_RESULT(numCohorts),
    BSPWBres(SPWBres){}
};
Rcpp::List copyBasicGROWTHResult_c(BasicGROWTH_RESULT& GROWTHres, ModelInput& x);

struct AdvancedGROWTH_RESULT : GROWTH_RESULT {
  AdvancedSPWB_RESULT ASPWBres;
  LabileCarbonBalanceInst_RESULT LCBInstres;
  
  AdvancedGROWTH_RESULT(AdvancedSPWB_RESULT SPWBres, size_t numCohorts, size_t ntimesteps) : 
  GROWTH_RESULT(numCohorts),
  ASPWBres(SPWBres),
  LCBInstres(numCohorts, ntimesteps) {}
};
Rcpp::List copyAdvancedGROWTHResult_c(AdvancedGROWTH_RESULT& GROWTHres, ModelInput& x);


struct InitialFinalCarbonCompartments{
  CarbonCompartments ccIni_g_ind;
  CarbonCompartments ccFin_g_ind;
  CarbonCompartments ccIni_gC_m2;
  CarbonCompartments ccFin_gC_m2;
  
  InitialFinalCarbonCompartments(size_t numCohorts) :
    ccIni_g_ind(numCohorts),
    ccFin_g_ind(numCohorts),
    ccIni_gC_m2(numCohorts),
    ccFin_gC_m2(numCohorts) {}
};

struct GROWTHCommunicationStructures {
  WBCommunicationStructures WBcomm;
  InitialFinalCarbonCompartments initialFinalCC;
  Decomposition_COMM DECcomm;
  
  GROWTHCommunicationStructures(size_t numCohorts, size_t nlayers, size_t ncanlayers, size_t ntimesteps) : 
    WBcomm(numCohorts, nlayers, ncanlayers, ntimesteps),
    initialFinalCC(numCohorts),
    DECcomm(7) {}
};

double dailyMortalityProbability_c(double stressValue, double stressThreshold);
double phloemFlow_c(double psiUpstream, double psiDownstream,
                    double concUpstream, double concDownstream,
                    double temp, double k_f, double nonSugarConc);
double qResp_c(double Tmean);

void growthDay_private_c(GROWTH_RESULT& GROWTHres, GROWTHCommunicationStructures& GROWTHcomm, ModelInput& x, 
                         WeatherInputVector meteovec, 
                         double latitude, double elevation, double slope, double aspect,
                         double solarConstant, double delta, 
                         const double runon, 
                         const std::vector<double>& lateralFlows, const double waterTableDepth);

void growthDay_inner_c(GROWTH_RESULT& GROWTHres, GROWTHCommunicationStructures& GROWTHcomm, ModelInput& x, 
                       std::string date,
                       WeatherInputVector meteovec, 
                       double latitude, double elevation, double slope, double aspect,
                       const double runon, 
                       const std::vector<double>& lateralFlows, const double waterTableDepth);
#endif
