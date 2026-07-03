#include "medfate.h"
#include "RcppArmadillo.h"

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

struct PhenologyControlParams {
  bool leafPhenology;
  double unfoldingDD;
};

struct MistletoeParams {
  double kPAR = 0.5;
  double g = 0.8;
  double LeafWidth = 1.0;
  double Tmax_LAI = 0.134;
  double Tmax_LAIsq = -0.006;
  double Gsw_P50_Baldocchi = -2.5;
  double Gsw_slope_Baldocchi = 30.0;
  double Gsw_AC_slope_Baldocchi = 8.0;
  double Beta_p =  1.907817;
  double Beta_q = 1.289641;
  double ClumpingIndex = 0.75;
  double alphaSWR = 0.7;
  double gammaSWR = 0.14;
  double Vmax298 = 80.0;
  double Jmax298 = 120.0;
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
  bool cavitationInducedDefoliation;
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

struct GrowthControlParams {
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
  double BelowgroundMetabolic;
  double Leaves;
  double FineRoots;
  double Twigs;
  double SmallBranches;
  double LargeWood;
  double CoarseRoots;
  double SurfaceActive;
  double BelowgroundActive;
  double SurfaceSlow;
  double BelowgroundSlow;
  double BelowgroundPassive;
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
  PhenologyControlParams phenology;
  CommonWBParams commonWB;
  BasicWBParams basicWB;
  AdvancedWBParams advancedWB;
  SperryWBParams sperry;
  SureauWBParams sureau;
  MistletoeParams mistletoe;
  FireHazardParams fireHazard;
  GrowthControlParams growth;
  MortalityParams mortality;
  RecruitmentParams recruitment;
  DecompositionParams decomposition;
  
  ControlParameters();
  ControlParameters(Rcpp::List x);
};

#endif


