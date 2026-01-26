#include "medfate.h"
#include "Rcpp.h"

#ifndef CONTROL_H
#define CONTROL_H

struct FillMissing {
  bool fillMissingRootParams;
  bool fillMissingSpParams;
  bool fillMissingWithGenusParams;
};

struct Results {
  bool subdailyResults;
  bool standResults;
  bool soilResults;
  bool soilPoolResults;
  bool snowResults;
  bool plantResults;
  bool labileCarbonBalanceResults;
  bool plantStructureResults;
  bool growthMortalityResults;
  bool decompositionPoolResults;
  bool leafResults;
  bool temperatureResults;
  bool fireHazardResults;
};

struct WeatherParams {
  double windMeasurementHeight;
  double defaultWindSpeed;
  double defaultCO2;
  std::vector<double> defaultRainfallIntensityPerMonth;
};

struct PhenoParams {
  bool leafPhenology;
  double unfoldingDD;
};

struct CommonWBParams {
  double verticalLayerSize;
  bool bareSoilEvaporation;
  bool unlimitedSoilWater;
  std::string interceptionMode;
  std::string infiltrationMode;
  std::string stemCavitationRecovery;
  std::string leafCavitationRecovery;
  bool segmentedXylemVulnerability;
  double cavitationRecoveryMaximumRate;
  bool truncateRootDistribution;
  double fullRhizosphereOverlapConductivity;
  std::string soilFunctions;
  std::string VG_PTF;
  int max_nsubsteps_soil;
  double infiltrationCorrection;
};

struct BasicWBParams {
  double hydraulicRedistributionFraction;
};

struct AdvancedWBParams {
  int ndailysteps;
  int nsubsteps_canopy;
  bool taper;
  bool multiLayerBalance;
  bool sapFluidityVariation;
  double TPhase_gmin;
  double Q10_1_gmin;
  double Q10_2_gmin;
  double rootRadialConductance;
  double averageFracRhizosphereResistance;
  double thermalCapacityLAI;
  double boundaryLayerSize;
  bool sunlitShade;
};

struct SperryWBParams {
  bool leafCavitationEffects;
  bool stemCavitationEffects;
  int maxNsteps;
  int ntrial;
  double psiTol;
  double ETol;
};

struct SureauWBParams {
  std::string stomatalSubmodel;
  bool plantCapacitance;
  bool cavitationFlux;
  bool leafCuticularTranspiration;
  bool stemCuticularTranspiration;
  double C_SApoInit;
  double C_LApoInit;
  double k_SSym;
  double fractionLeafSymplasm;
  double gs_NightFrac;
  double JarvisPAR;
  double fTRBToLeaf;
};

struct FireHazardParams {
  double standardWind;
  double standardDFMC;
  std::string lfmcComponent;
};

struct EquilibriumOsmoticConcentration {
  double leaf;
  double sapwood;
};

struct ConstructionCosts {
  double leaf;
  double sapwood;
  double fineroot;
};

struct SenescenceRates {
  double sapwood;
  double fineroot;
};

struct MaximumRelativeGrowthRates {
  double leaf;
  double cambium;
  double sapwood;
  double fineroot;
};

struct GrowthParams {
  bool subdailyCarbonBalance;
  bool allowDessication;
  bool allowStarvation;
  bool sinkLimitation;
  bool shrubDynamics;
  bool herbDynamics;
  std::string allocationStrategy;
  double phloemConductanceFactor;
  double nonSugarConcentration;
  EquilibriumOsmoticConcentration equilibriumOsmoticConcentration;
  double minimumRelativeStarchForGrowth;
  ConstructionCosts constructionCosts;
  SenescenceRates senescenceRates;
  MaximumRelativeGrowthRates maximumRelativeGrowthRates;
};

struct MortalityParams {
  std::string mortalityMode;
  double mortalityBaselineRate;
  double mortalityRelativeSugarThreshold;
  double mortalityRWCThreshold;
};
struct RecruitmentParams {
  double recrTreeDBH;
  double recrTreeDensity;
  double ingrowthTreeDBH;
  double ingrowthTreeDensity;  
};
struct DecompositionAnnualBaseRates {
  double SurfaceMetabolic;
  double SoilMetabolic;
  double Leaves;
  double FineRoots;
  double Twigs;
  double SmallBranches;
  double LargeWood;
  double CoarseRoots;
  double SurfaceActive;
  double SoilActive;
  double SurfaceSlow;
  double SoilSlow;
  double SoilPassive;
};

struct DecompositionParams {
  
  DecompositionAnnualBaseRates annualBaseRates;
  double decompositionAnnualTurnoverRate;
};
// Create struct for control parameters paralleling the output of R function defaultControl()
struct ControlParameters {

  bool verbose;
  
  FillMissing fillMissing;
  
  Results results;
  
  std::string transpirationMode;
  std::string soilDomains;
  std::string rhizosphereOverlap;
  
  WeatherParams weather;
  PhenoParams phenology;
  CommonWBParams commonWB;
  BasicWBParams basicWB;
  AdvancedWBParams advancedWB;
  SperryWBParams sperry;
  SureauWBParams sureau;
  FireHazardParams fireHazard;
  GrowthParams growth;
  MortalityParams mortality;
  RecruitmentParams recruitment;
  DecompositionParams decomposition;
  
  ControlParameters();
  ControlParameters(Rcpp::List x);
};

#endif


