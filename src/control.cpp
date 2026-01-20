#include <Rcpp.h>
#include "control.h"
using namespace Rcpp;

// Copies Rcpp control list to ControlParameters
ControlParameters controlListToStructure(List x) {
  ControlParameters ctl;
  ctl.fillMissingRootParams = as<bool>(x["fillMissingRootParams"]);
  ctl.fillMissingSpParams = as<bool>(x["fillMissingSpParams"]);
  ctl.fillMissingWithGenusParams = as<bool>(x["fillMissingWithGenusParams"]);
  ctl.verbose = as<bool>(x["verbose"]);
  ctl.subdailyResults = as<bool>(x["subdailyResults"]);
  ctl.standResults = as<bool>(x["standResults"]);
  ctl.soilResults = as<bool>(x["soilResults"]);
  ctl.soilPoolResults = as<bool>(x["soilPoolResults"]);
  ctl.snowResults = as<bool>(x["snowResults"]);
  ctl.plantResults = as<bool>(x["plantResults"]);
  ctl.labileCarbonBalanceResults = as<bool>(x["labileCarbonBalanceResults"]);
  ctl.plantStructureResults = as<bool>(x["plantStructureResults"]);
  ctl.growthMortalityResults = as<bool>(x["growthMortalityResults"]);
  ctl.decompositionPoolResults = as<bool>(x["decompositionPoolResults"]);
  ctl.leafResults = as<bool>(x["leafResults"]);
  ctl.temperatureResults = as<bool>(x["temperatureResults"]);
  ctl.fireHazardResults = as<bool>(x["fireHazardResults"]);
  ctl.fireHazardStandardWind = as<double>(x["fireHazardStandardWind"]);
  ctl.fireHazardStandardDFMC = as<double>(x["fireHazardStandardDFMC"]);
  
  ctl.transpirationMode = as<std::string>(x["transpirationMode"]);
  ctl.soilDomains = as<std::string>(x["soilDomains"]);
  ctl.rhizosphereOverlap = as<std::string>(x["rhizosphereOverlap"]);
  ctl.truncateRootDistribution = as<bool>(x["truncateRootDistribution"]);
  ctl.fullRhizosphereOverlapConductivity = as<double>(x["fullRhizosphereOverlapConductivity"]);
  ctl.soilFunctions = as<std::string>(x["soilFunctions"]);
  ctl.VG_PTF = as<std::string>(x["VG_PTF"]);
  ctl.ndailysteps = as<int>(x["ndailysteps"]);
  ctl.max_nsubsteps_soil = as<int>(x["max_nsubsteps_soil"]);
  ctl.defaultWindSpeed = as<double>(x["defaultWindSpeed"]);
  ctl.defaultCO2 = as<double>(x["defaultCO2"]);
  ctl.defaultRainfallIntensityPerMonth = as< std::vector<double> >(x["defaultRainfallIntensityPerMonth"]);
  ctl.leafPhenology = as<bool>(x["leafPhenology"]);
  ctl.bareSoilEvaporation = as<bool>(x["bareSoilEvaporation"]);
  ctl.unlimitedSoilWater = as<bool>(x["unlimitedSoilWater"]);
  ctl.interceptionMode = as<std::string>(x["interceptionMode"]);
  ctl.infiltrationMode = as<std::string>(x["infiltrationMode"]);
  ctl.infiltrationCorrection = as<double>(x["infiltrationCorrection"]);
  ctl.unfoldingDD = as<double>(x["unfoldingDD"]);
  ctl.verticalLayerSize = as<double>(x["verticalLayerSize"]);
  ctl.windMeasurementHeight = as<double>(x["windMeasurementHeight"]);
  ctl.segmentedXylemVulnerability = as<bool>(x["segmentedXylemVulnerability"]);
  ctl.stemCavitationRecovery = as<std::string>(x["stemCavitationRecovery"]);
  ctl.leafCavitationRecovery = as<std::string>(x["leafCavitationRecovery"]);
  ctl.lfmcComponent = as<std::string>(x["lfmcComponent"]);
  ctl.hydraulicRedistributionFraction = as<double>(x["hydraulicRedistributionFraction"]);
  ctl.nsubsteps_canopy = as<int>(x["nsubsteps_canopy"]);
  ctl.taper = as<bool>(x["taper"]);
  ctl.multiLayerBalance = as<bool>(x["multiLayerBalance"]);
  ctl.sapFluidityVariation = as<bool>(x["sapFluidityVariation"]);
  ctl.TPhase_gmin = as<double>(x["TPhase_gmin"]);
  ctl.Q10_1_gmin = as<double>(x["Q10_1_gmin"]);
  ctl.Q10_2_gmin = as<double>(x["Q10_2_gmin"]);
  ctl.rootRadialConductance = as<double>(x["rootRadialConductance"]);
  ctl.averageFracRhizosphereResistance = as<double>(x["averageFracRhizosphereResistance"]);
  ctl.thermalCapacityLAI = as<double>(x["thermalCapacityLAI"]);
  ctl.boundaryLayerSize = as<double>(x["boundaryLayerSize"]);
  ctl.cavitationRecoveryMaximumRate = as<double>(x["cavitationRecoveryMaximumRate"]);
  ctl.sunlitShade = as<bool>(x["sunlitShade"]);
  
  // List numericParams = x["numericParams"];
  // ctl.snp.maxNsteps = as<int>(numericParams["maxNsteps"]);
  // ctl.snp.ntrial = as<int>(numericParams["ntrial"]);
  // ctl.snp.psiTol = as<double>(numericParams["psiTol"]);
  // ctl.snp.ETol = as<double>(numericParams["ETol"]);
  
  ctl.leafCavitationEffects = as<bool>(x["leafCavitationEffects"]);
  ctl.stemCavitationEffects = as<bool>(x["leafCavitationEffects"]);

  return(ctl);
}

// [[Rcpp::export(.testControlListToStructure)]]
void testControlListToStructure(List x) {
  ControlParameters ctl = controlListToStructure(x);
  Rcout << ctl.sunlitShade << "\n";
}
