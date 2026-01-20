#include <Rcpp.h>
#include "control.h"
using namespace Rcpp;

// Copies Rcpp control list to ControlParameters
ControlParameters controlListToStructure(List x) {
  ControlParameters ctl;
  ctl.verbose = as<bool>(x["verbose"]);
  ctl.transpirationMode = as<std::string>(x["transpirationMode"]);
  ctl.soilDomains = as<std::string>(x["soilDomains"]);
  ctl.rhizosphereOverlap = as<std::string>(x["rhizosphereOverlap"]);
  
  ctl.fillMissing.fillMissingRootParams = as<bool>(x["fillMissingRootParams"]);
  ctl.fillMissing.fillMissingSpParams = as<bool>(x["fillMissingSpParams"]);
  ctl.fillMissing.fillMissingWithGenusParams = as<bool>(x["fillMissingWithGenusParams"]);
  
  ctl.results.subdailyResults = as<bool>(x["subdailyResults"]);
  ctl.results.standResults = as<bool>(x["standResults"]);
  ctl.results.soilResults = as<bool>(x["soilResults"]);
  ctl.results.soilPoolResults = as<bool>(x["soilPoolResults"]);
  ctl.results.snowResults = as<bool>(x["snowResults"]);
  ctl.results.plantResults = as<bool>(x["plantResults"]);
  ctl.results.labileCarbonBalanceResults = as<bool>(x["labileCarbonBalanceResults"]);
  ctl.results.plantStructureResults = as<bool>(x["plantStructureResults"]);
  ctl.results.growthMortalityResults = as<bool>(x["growthMortalityResults"]);
  ctl.results.decompositionPoolResults = as<bool>(x["decompositionPoolResults"]);
  ctl.results.leafResults = as<bool>(x["leafResults"]);
  ctl.results.temperatureResults = as<bool>(x["temperatureResults"]);
  ctl.results.fireHazardResults = as<bool>(x["fireHazardResults"]);
  
  ctl.fireHazard.standardWind = as<double>(x["fireHazardStandardWind"]);
  ctl.fireHazard.standardDFMC = as<double>(x["fireHazardStandardDFMC"]);
  
  
  ctl.weather.windMeasurementHeight = as<double>(x["windMeasurementHeight"]);
  ctl.weather.defaultWindSpeed = as<double>(x["defaultWindSpeed"]);
  ctl.weather.defaultCO2 = as<double>(x["defaultCO2"]);
  ctl.weather.defaultRainfallIntensityPerMonth = as< std::vector<double> >(x["defaultRainfallIntensityPerMonth"]);
  
  ctl.phenology.leafPhenology = as<bool>(x["leafPhenology"]);
  ctl.phenology.unfoldingDD = as<double>(x["unfoldingDD"]);
  
  ctl.fireHazard.lfmcComponent = as<std::string>(x["lfmcComponent"]);

  ctl.commonWB.truncateRootDistribution = as<bool>(x["truncateRootDistribution"]);
  ctl.commonWB.fullRhizosphereOverlapConductivity = as<double>(x["fullRhizosphereOverlapConductivity"]);
  ctl.commonWB.soilFunctions = as<std::string>(x["soilFunctions"]);
  ctl.commonWB.VG_PTF = as<std::string>(x["VG_PTF"]);
  ctl.commonWB.max_nsubsteps_soil = as<int>(x["max_nsubsteps_soil"]);
  ctl.commonWB.infiltrationCorrection = as<double>(x["infiltrationCorrection"]);
  ctl.commonWB.verticalLayerSize = as<double>(x["verticalLayerSize"]);
  ctl.commonWB.bareSoilEvaporation = as<bool>(x["bareSoilEvaporation"]);
  ctl.commonWB.unlimitedSoilWater = as<bool>(x["unlimitedSoilWater"]);
  ctl.commonWB.interceptionMode = as<std::string>(x["interceptionMode"]);
  ctl.commonWB.infiltrationMode = as<std::string>(x["infiltrationMode"]);
  ctl.commonWB.cavitationRecoveryMaximumRate = as<double>(x["cavitationRecoveryMaximumRate"]);
  ctl.commonWB.stemCavitationRecovery = as<std::string>(x["stemCavitationRecovery"]);
  ctl.commonWB.leafCavitationRecovery = as<std::string>(x["leafCavitationRecovery"]);
  ctl.commonWB.segmentedXylemVulnerability = as<bool>(x["segmentedXylemVulnerability"]);

  ctl.basicWB.hydraulicRedistributionFraction = as<double>(x["hydraulicRedistributionFraction"]);
  
  ctl.advancedWB.nsubsteps_canopy = as<int>(x["nsubsteps_canopy"]);
  ctl.advancedWB.ndailysteps = as<int>(x["ndailysteps"]);
  ctl.advancedWB.taper = as<bool>(x["taper"]);
  ctl.advancedWB.multiLayerBalance = as<bool>(x["multiLayerBalance"]);
  ctl.advancedWB.sapFluidityVariation = as<bool>(x["sapFluidityVariation"]);
  ctl.advancedWB.TPhase_gmin = as<double>(x["TPhase_gmin"]);
  ctl.advancedWB.Q10_1_gmin = as<double>(x["Q10_1_gmin"]);
  ctl.advancedWB.Q10_2_gmin = as<double>(x["Q10_2_gmin"]);
  ctl.advancedWB.rootRadialConductance = as<double>(x["rootRadialConductance"]);
  ctl.advancedWB.averageFracRhizosphereResistance = as<double>(x["averageFracRhizosphereResistance"]);
  ctl.advancedWB.thermalCapacityLAI = as<double>(x["thermalCapacityLAI"]);
  ctl.advancedWB.boundaryLayerSize = as<double>(x["boundaryLayerSize"]);
  ctl.advancedWB.sunlitShade = as<bool>(x["sunlitShade"]);

  
  ctl.sperry.leafCavitationEffects = as<bool>(x["leafCavitationEffects"]);
  ctl.sperry.stemCavitationEffects = as<bool>(x["stemCavitationEffects"]);
  List numericParams = x["numericParams"];
  ctl.sperry.maxNsteps = as<int>(numericParams["maxNsteps"]);
  ctl.sperry.ntrial = as<int>(numericParams["ntrial"]);
  ctl.sperry.psiTol = as<double>(numericParams["psiTol"]);
  ctl.sperry.ETol = as<double>(numericParams["ETol"]);
  
  ctl.sureau.stomatalSubmodel = as<std::string>(x["stomatalSubmodel"]);
  ctl.sureau.plantCapacitance = as<bool>(x["plantCapacitance"]);
  ctl.sureau.cavitationFlux = as<bool>(x["cavitationFlux"]);
  ctl.sureau.leafCuticularTranspiration = as<bool>(x["leafCuticularTranspiration"]);
  ctl.sureau.stemCuticularTranspiration = as<bool>(x["stemCuticularTranspiration"]);
  ctl.sureau.C_SApoInit = as<double>(x["C_SApoInit"]);
  ctl.sureau.C_LApoInit = as<double>(x["C_LApoInit"]);
  ctl.sureau.k_SSym = as<double>(x["k_SSym"]);
  ctl.sureau.fractionLeafSymplasm = as<double>(x["fractionLeafSymplasm"]);
  ctl.sureau.gs_NightFrac = as<double>(x["gs_NightFrac"]);
  ctl.sureau.JarvisPAR = as<double>(x["JarvisPAR"]);
  ctl.sureau.fTRBToLeaf = as<double>(x["fTRBToLeaf"]);
  
  ctl.growth.subdailyCarbonBalance = as<bool>(x["subdailyCarbonBalance"]);
  ctl.growth.allowDessication = as<bool>(x["allowDessication"]);
  ctl.growth.allowStarvation = as<bool>(x["allowStarvation"]);
  ctl.growth.sinkLimitation = as<bool>(x["sinkLimitation"]);
  ctl.growth.shrubDynamics = as<bool>(x["shrubDynamics"]);
  ctl.growth.herbDynamics = as<bool>(x["herbDynamics"]);
  ctl.growth.allocationStrategy = as<std::string>(x["allocationStrategy"]);
  ctl.growth.phloemConductanceFactor = as<double>(x["phloemConductanceFactor"]);
  ctl.growth.nonSugarConcentration = as<double>(x["nonSugarConcentration"]);
  
  List eqOsmConcs = x["equilibriumOsmoticConcentration"];
  ctl.growth.equilibriumOsmoticConcentration.leaf = as<double>(eqOsmConcs["leaf"]);
  ctl.growth.equilibriumOsmoticConcentration.sapwood = as<double>(eqOsmConcs["sapwood"]);
  ctl.growth.minimumRelativeStarchForGrowth = as<double>(x["minimumRelativeStarchForGrowth"]);
  List constructionCosts = x["constructionCosts"];
  ctl.growth.constructionCosts.leaf = as<double>(constructionCosts["leaf"]);
  ctl.growth.constructionCosts.sapwood = as<double>(constructionCosts["sapwood"]);
  ctl.growth.constructionCosts.fineroot = as<double>(constructionCosts["fineroot"]);
  List senescenceRates = x["senescenceRates"];
  ctl.growth.senescenceRates.sapwood = as<double>(senescenceRates["sapwood"]);
  ctl.growth.senescenceRates.fineroot = as<double>(senescenceRates["fineroot"]);
  List maximumRelativeGrowthRates = x["maximumRelativeGrowthRates"];
  ctl.growth.maximumRelativeGrowthRates.leaf = as<double>(maximumRelativeGrowthRates["leaf"]);
  ctl.growth.maximumRelativeGrowthRates.cambium = as<double>(maximumRelativeGrowthRates["cambium"]);
  ctl.growth.maximumRelativeGrowthRates.sapwood = as<double>(maximumRelativeGrowthRates["sapwood"]);
  ctl.growth.maximumRelativeGrowthRates.fineroot = as<double>(maximumRelativeGrowthRates["fineroot"]);
  
  ctl.mortality.mortalityMode = as<std::string>(x["mortalityMode"]);
  ctl.mortality.mortalityBaselineRate = as<double>(x["mortalityBaselineRate"]);
  ctl.mortality.mortalityRelativeSugarThreshold = as<double>(x["mortalityRelativeSugarThreshold"]);
  ctl.mortality.mortalityRWCThreshold = as<double>(x["mortalityRWCThreshold"]);
  
  ctl.recruitment.recrTreeDBH = as<double>(x["recrTreeDBH"]);
  ctl.recruitment.recrTreeDensity = as<double>(x["recrTreeDensity"]);
  ctl.recruitment.ingrowthTreeDBH = as<double>(x["ingrowthTreeDBH"]);
  ctl.recruitment.ingrowthTreeDensity = as<double>(x["ingrowthTreeDensity"]);

  NumericVector decompRates = x["decompositionAnnualBaseRates"];  
  ctl.decomposition.decompositionAnnualTurnoverRate = as<double>(x["decompositionAnnualTurnoverRate"]);
  ctl.decomposition.annualBaseRates.SurfaceMetabolic = decompRates["SurfaceMetabolic"];
  ctl.decomposition.annualBaseRates.SoilMetabolic = decompRates["SoilMetabolic"];
  ctl.decomposition.annualBaseRates.Leaves = decompRates["Leaves"];
  ctl.decomposition.annualBaseRates.FineRoots = decompRates["FineRoots"];
  ctl.decomposition.annualBaseRates.Twigs = decompRates["Twigs"];
  ctl.decomposition.annualBaseRates.SmallBranches = decompRates["SmallBranches"];
  ctl.decomposition.annualBaseRates.LargeWood = decompRates["LargeWood"];
  ctl.decomposition.annualBaseRates.CoarseRoots = decompRates["CoarseRoots"];
  ctl.decomposition.annualBaseRates.SurfaceActive = decompRates["SurfaceActive"];
  ctl.decomposition.annualBaseRates.SurfaceSlow = decompRates["SurfaceSlow"];
  ctl.decomposition.annualBaseRates.SoilActive = decompRates["SoilActive"];
  ctl.decomposition.annualBaseRates.SoilSlow = decompRates["SoilSlow"];
  ctl.decomposition.annualBaseRates.SoilPassive = decompRates["SoilPassive"];  
  
  return(ctl);
}

// [[Rcpp::export(.testControlListToStructure)]]
int testControlListToStructure(List x) {
  ControlParameters ctl = controlListToStructure(x);
  // Rcout << ctl.sureau.fTRBToLeaf << "\n";
  return(sizeof(ctl));
}
