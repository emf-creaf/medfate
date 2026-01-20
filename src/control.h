#include "medfate.h"

#ifndef CONTROL_H
#define CONTROL_H


// struct SperryNumericParams {
//   int maxNsteps;
//   int ntrial;
//   double psiTol;
//   double ETol;
// };

// Create struct for control parameters paralleling the output of R function defaultControl()
struct ControlParameters {
  bool fillMissingRootParams;
  bool fillMissingSpParams;
  bool fillMissingWithGenusParams;
  bool verbose;
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
  double fireHazardStandardWind;
  double fireHazardStandardDFMC;
  
  std::string transpirationMode;
  std::string soilDomains;
  std::string rhizosphereOverlap;
  bool truncateRootDistribution;
  double fullRhizosphereOverlapConductivity;
  
  std::string soilFunctions;
  std::string VG_PTF;
  
  int ndailysteps;
  int max_nsubsteps_soil;
  
  double defaultWindSpeed;
  double defaultCO2 = 386; // #ppm
  std::vector<double> defaultRainfallIntensityPerMonth;
  bool leafPhenology;
  bool bareSoilEvaporation;
  bool unlimitedSoilWater;
  std::string interceptionMode;
  std::string infiltrationMode;
  double infiltrationCorrection;
  double unfoldingDD;
  double verticalLayerSize;
  double windMeasurementHeight;
  bool segmentedXylemVulnerability;
  std::string stemCavitationRecovery;
  std::string leafCavitationRecovery;
  std::string lfmcComponent;
  
  //spwb with granier
  double hydraulicRedistributionFraction;
    
  // spwb with sperry/sureau
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
  double cavitationRecoveryMaximumRate;
  bool sunlitShade;
  
  // SperryNumericParams snp;
  bool leafCavitationEffects;
  bool stemCavitationEffects;
    
  // spwb with sureau
  std::string stomatalSubmodel;
  bool plantCapacitance;
  bool cavitationFlux;
  bool leafCuticularTranspiration;
  bool stemCuticularTranspiration;
  double C_SApoInit = 2.0e-05;
  double C_LApoInit = 1.0e-05;
  double k_SSym = 0.26;
  double fractionLeafSymplasm = 0.5;
  double gs_NightFrac = 0.05;
  double JarvisPAR = 0.003;
  double fTRBToLeaf = 0.8;
    
  // growth/mortality
  bool subdailyCarbonBalance;
  bool allowDessication;
  bool allowStarvation;
  bool sinkLimitation;
  bool shrubDynamics;
  bool herbDynamics;
  std::string allocationStrategy;
  double phloemConductanceFactor;
  double nonSugarConcentration;
    // equilibriumOsmoticConcentration = list(leaf = 0.8, sapwood = 0.6),  # (Paljakka et al. 2017)
  double minimumRelativeStarchForGrowth;
        // constructionCosts = list(leaf = 1.5, 
        //                          sapwood = 1.47, 
        //                          fineroot = 1.30), #  g gluc · g dw -1
    // senescenceRates = list(sapwood = 0.000135, # day-1 Equivalent to annual 4.8% 1-(1-0.048)^(1.0/365)
    //                          fineroot = 0.001897231), #day-1 Equivalent to annual 50% 1-(1-0.5)^(1.0/365)
    //   maximumRelativeGrowthRates = list(leaf = 0.09, # m2 leaf ·cm-2 sapwood· day-1
    //                                       cambium = 0.0025, # cm2 sapwood ·cm-1 cambium· day-1
    //                                       sapwood = 0.002, # cm2 sapwood ·cm-2 sapwood· day-1
    //                                       fineroot = 0.1), # g dw · g dw -1 · day -1
    std::string mortalityMode;
    double mortalityBaselineRate;
    double mortalityRelativeSugarThreshold;
    double mortalityRWCThreshold;
    double recrTreeDBH;
    double recrTreeDensity;
    double ingrowthTreeDBH;
    double ingrowthTreeDensity;
      // decompositionAnnualBaseRates = c("SurfaceMetabolic" = 8.0, 
      //                                  "SoilMetabolic" = 18.5, 
      //                                  "Leaves" = 2.0, 
      //                                  "FineRoots" = 4.9,
      //                                  "Twigs" =  1.8, 
      //                                  "SmallBranches" =  1.5, 
      //                                  "LargeWood" = 0.02, 
      //                                  "CoarseRoots" = 0.1, 
      //                                  "SurfaceActive" = 6.0, 
      //                                  "SoilActive" = 11.0, 
      //                                  "SurfaceSlow" = 0.08, 
      //                                  "SoilSlow" = 0.4, 
      //                                  "SoilPassive" = 0.0033),
      double decompositionAnnualTurnoverRate;
      
      bool allowSeedBankDynamics;
      bool allowRecruitment;
      bool allowResprouting;
      std::string recruitmentMode;
      bool removeEmptyCohorts;
      double minimumTreeCohortDensity;
      double minimumShrubCohortCover;
      bool dynamicallyMergeCohorts;
      bool keepCohortsWithObsID;
      // seedRain = NULL,
      double seedProductionTreeHeight;
      double seedProductionShrubHeight;
      double probRecr;
      double minTempRecr;
      double minMoistureRecr;
      double minFPARRecr;
      double recrAge;
      double recrTreeHeight;
      double recrShrubCover;
      double recrShrubHeight;
};


#endif
