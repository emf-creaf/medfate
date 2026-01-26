#include <RcppArmadillo.h>
#include "control.h"
using namespace Rcpp;

/**
 * Implementation of ControlParameters class
 */
ControlParameters::ControlParameters(){
  
}
ControlParameters::ControlParameters(List x) {
  verbose = as<bool>(x["verbose"]);
  transpirationMode = as<std::string>(x["transpirationMode"]);
  soilDomains = as<std::string>(x["soilDomains"]);
  rhizosphereOverlap = as<std::string>(x["rhizosphereOverlap"]);
  
  fillMissing.fillMissingRootParams = as<bool>(x["fillMissingRootParams"]);
  fillMissing.fillMissingSpParams = as<bool>(x["fillMissingSpParams"]);
  fillMissing.fillMissingWithGenusParams = as<bool>(x["fillMissingWithGenusParams"]);
  
  results.subdailyResults = as<bool>(x["subdailyResults"]);
  results.standResults = as<bool>(x["standResults"]);
  results.soilResults = as<bool>(x["soilResults"]);
  results.soilPoolResults = as<bool>(x["soilPoolResults"]);
  results.snowResults = as<bool>(x["snowResults"]);
  results.plantResults = as<bool>(x["plantResults"]);
  results.labileCarbonBalanceResults = as<bool>(x["labileCarbonBalanceResults"]);
  results.plantStructureResults = as<bool>(x["plantStructureResults"]);
  results.growthMortalityResults = as<bool>(x["growthMortalityResults"]);
  results.decompositionPoolResults = as<bool>(x["decompositionPoolResults"]);
  results.leafResults = as<bool>(x["leafResults"]);
  results.temperatureResults = as<bool>(x["temperatureResults"]);
  results.fireHazardResults = as<bool>(x["fireHazardResults"]);
  
  fireHazard.standardWind = as<double>(x["fireHazardStandardWind"]);
  fireHazard.standardDFMC = as<double>(x["fireHazardStandardDFMC"]);
  
  
  weather.windMeasurementHeight = as<double>(x["windMeasurementHeight"]);
  weather.defaultWindSpeed = as<double>(x["defaultWindSpeed"]);
  weather.defaultCO2 = as<double>(x["defaultCO2"]);
  weather.defaultRainfallIntensityPerMonth = as< std::vector<double> >(x["defaultRainfallIntensityPerMonth"]);
  
  phenology.leafPhenology = as<bool>(x["leafPhenology"]);
  phenology.unfoldingDD = as<double>(x["unfoldingDD"]);
  
  fireHazard.lfmcComponent = as<std::string>(x["lfmcComponent"]);

  commonWB.truncateRootDistribution = as<bool>(x["truncateRootDistribution"]);
  commonWB.fullRhizosphereOverlapConductivity = as<double>(x["fullRhizosphereOverlapConductivity"]);
  commonWB.soilFunctions = as<std::string>(x["soilFunctions"]);
  commonWB.VG_PTF = as<std::string>(x["VG_PTF"]);
  commonWB.max_nsubsteps_soil = as<int>(x["max_nsubsteps_soil"]);
  commonWB.infiltrationCorrection = as<double>(x["infiltrationCorrection"]);
  commonWB.verticalLayerSize = as<double>(x["verticalLayerSize"]);
  commonWB.bareSoilEvaporation = as<bool>(x["bareSoilEvaporation"]);
  commonWB.unlimitedSoilWater = as<bool>(x["unlimitedSoilWater"]);
  commonWB.interceptionMode = as<std::string>(x["interceptionMode"]);
  commonWB.infiltrationMode = as<std::string>(x["infiltrationMode"]);
  commonWB.cavitationRecoveryMaximumRate = as<double>(x["cavitationRecoveryMaximumRate"]);
  commonWB.stemCavitationRecovery = as<std::string>(x["stemCavitationRecovery"]);
  commonWB.leafCavitationRecovery = as<std::string>(x["leafCavitationRecovery"]);
  commonWB.segmentedXylemVulnerability = as<bool>(x["segmentedXylemVulnerability"]);

  basicWB.hydraulicRedistributionFraction = as<double>(x["hydraulicRedistributionFraction"]);
  
  advancedWB.nsubsteps_canopy = as<int>(x["nsubsteps_canopy"]);
  advancedWB.ndailysteps = as<int>(x["ndailysteps"]);
  advancedWB.taper = as<bool>(x["taper"]);
  advancedWB.multiLayerBalance = as<bool>(x["multiLayerBalance"]);
  advancedWB.sapFluidityVariation = as<bool>(x["sapFluidityVariation"]);
  advancedWB.TPhase_gmin = as<double>(x["TPhase_gmin"]);
  advancedWB.Q10_1_gmin = as<double>(x["Q10_1_gmin"]);
  advancedWB.Q10_2_gmin = as<double>(x["Q10_2_gmin"]);
  advancedWB.rootRadialConductance = as<double>(x["rootRadialConductance"]);
  advancedWB.averageFracRhizosphereResistance = as<double>(x["averageFracRhizosphereResistance"]);
  advancedWB.thermalCapacityLAI = as<double>(x["thermalCapacityLAI"]);
  advancedWB.boundaryLayerSize = as<double>(x["boundaryLayerSize"]);
  advancedWB.sunlitShade = as<bool>(x["sunlitShade"]);

  
  sperry.leafCavitationEffects = as<bool>(x["leafCavitationEffects"]);
  sperry.stemCavitationEffects = as<bool>(x["stemCavitationEffects"]);
  List numericParams = x["numericParams"];
  sperry.maxNsteps = as<int>(numericParams["maxNsteps"]);
  sperry.ntrial = as<int>(numericParams["ntrial"]);
  sperry.psiTol = as<double>(numericParams["psiTol"]);
  sperry.ETol = as<double>(numericParams["ETol"]);
  
  sureau.stomatalSubmodel = as<std::string>(x["stomatalSubmodel"]);
  sureau.plantCapacitance = as<bool>(x["plantCapacitance"]);
  sureau.cavitationFlux = as<bool>(x["cavitationFlux"]);
  sureau.leafCuticularTranspiration = as<bool>(x["leafCuticularTranspiration"]);
  sureau.stemCuticularTranspiration = as<bool>(x["stemCuticularTranspiration"]);
  sureau.C_SApoInit = as<double>(x["C_SApoInit"]);
  sureau.C_LApoInit = as<double>(x["C_LApoInit"]);
  sureau.k_SSym = as<double>(x["k_SSym"]);
  sureau.fractionLeafSymplasm = as<double>(x["fractionLeafSymplasm"]);
  sureau.gs_NightFrac = as<double>(x["gs_NightFrac"]);
  sureau.JarvisPAR = as<double>(x["JarvisPAR"]);
  sureau.fTRBToLeaf = as<double>(x["fTRBToLeaf"]);
  
  growth.subdailyCarbonBalance = as<bool>(x["subdailyCarbonBalance"]);
  growth.allowDessication = as<bool>(x["allowDessication"]);
  growth.allowStarvation = as<bool>(x["allowStarvation"]);
  growth.sinkLimitation = as<bool>(x["sinkLimitation"]);
  growth.shrubDynamics = as<bool>(x["shrubDynamics"]);
  growth.herbDynamics = as<bool>(x["herbDynamics"]);
  growth.allocationStrategy = as<std::string>(x["allocationStrategy"]);
  growth.phloemConductanceFactor = as<double>(x["phloemConductanceFactor"]);
  growth.nonSugarConcentration = as<double>(x["nonSugarConcentration"]);
  
  List eqOsmConcs = x["equilibriumOsmoticConcentration"];
  growth.equilibriumOsmoticConcentration.leaf = as<double>(eqOsmConcs["leaf"]);
  growth.equilibriumOsmoticConcentration.sapwood = as<double>(eqOsmConcs["sapwood"]);
  growth.minimumRelativeStarchForGrowth = as<double>(x["minimumRelativeStarchForGrowth"]);
  List constructionCosts = x["constructionCosts"];
  growth.constructionCosts.leaf = as<double>(constructionCosts["leaf"]);
  growth.constructionCosts.sapwood = as<double>(constructionCosts["sapwood"]);
  growth.constructionCosts.fineroot = as<double>(constructionCosts["fineroot"]);
  List senescenceRates = x["senescenceRates"];
  growth.senescenceRates.sapwood = as<double>(senescenceRates["sapwood"]);
  growth.senescenceRates.fineroot = as<double>(senescenceRates["fineroot"]);
  List maximumRelativeGrowthRates = x["maximumRelativeGrowthRates"];
  growth.maximumRelativeGrowthRates.leaf = as<double>(maximumRelativeGrowthRates["leaf"]);
  growth.maximumRelativeGrowthRates.cambium = as<double>(maximumRelativeGrowthRates["cambium"]);
  growth.maximumRelativeGrowthRates.sapwood = as<double>(maximumRelativeGrowthRates["sapwood"]);
  growth.maximumRelativeGrowthRates.fineroot = as<double>(maximumRelativeGrowthRates["fineroot"]);
  
  mortality.mortalityMode = as<std::string>(x["mortalityMode"]);
  mortality.mortalityBaselineRate = as<double>(x["mortalityBaselineRate"]);
  mortality.mortalityRelativeSugarThreshold = as<double>(x["mortalityRelativeSugarThreshold"]);
  mortality.mortalityRWCThreshold = as<double>(x["mortalityRWCThreshold"]);
  
  recruitment.recrTreeDBH = as<double>(x["recrTreeDBH"]);
  recruitment.recrTreeDensity = as<double>(x["recrTreeDensity"]);
  recruitment.ingrowthTreeDBH = as<double>(x["ingrowthTreeDBH"]);
  recruitment.ingrowthTreeDensity = as<double>(x["ingrowthTreeDensity"]);

  NumericVector decompRates = x["decompositionAnnualBaseRates"];  
  decomposition.decompositionAnnualTurnoverRate = as<double>(x["decompositionAnnualTurnoverRate"]);
  decomposition.annualBaseRates.SurfaceMetabolic = decompRates["SurfaceMetabolic"];
  decomposition.annualBaseRates.SoilMetabolic = decompRates["SoilMetabolic"];
  decomposition.annualBaseRates.Leaves = decompRates["Leaves"];
  decomposition.annualBaseRates.FineRoots = decompRates["FineRoots"];
  decomposition.annualBaseRates.Twigs = decompRates["Twigs"];
  decomposition.annualBaseRates.SmallBranches = decompRates["SmallBranches"];
  decomposition.annualBaseRates.LargeWood = decompRates["LargeWood"];
  decomposition.annualBaseRates.CoarseRoots = decompRates["CoarseRoots"];
  decomposition.annualBaseRates.SurfaceActive = decompRates["SurfaceActive"];
  decomposition.annualBaseRates.SurfaceSlow = decompRates["SurfaceSlow"];
  decomposition.annualBaseRates.SoilActive = decompRates["SoilActive"];
  decomposition.annualBaseRates.SoilSlow = decompRates["SoilSlow"];
  decomposition.annualBaseRates.SoilPassive = decompRates["SoilPassive"];  
}

// [[Rcpp::export(.testControlListToStructure)]]
int testControlListToStructure(List x) {
  ControlParameters ctl = ControlParameters(x);
  // Rcout << sureau.fTRBToLeaf << "\n";
  return(sizeof(ctl));
}
