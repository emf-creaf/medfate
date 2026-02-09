#include <RcppArmadillo.h>
#include "biophysicsutils_c.h"
#include "meteoland/utils_c.hpp"
#include "lowlevel_structures_c.h"
#include "forestutils_c.h"
#include "transpiration_advanced_c.h"
#include "lightextinction_advanced_c.h"
#include "tissuemoisture_c.h"
#include "inner_sureau_c.h"
#include "inner_sperry_c.h"
#include "windextinction_c.h"
#include "windKatul_c.h"
using namespace Rcpp;


Rcpp::DataFrame copyPlantAdvancedTranspirationResult_c(const PlantsAdvancedTranspiration_RESULT& plants, ModelInput& x) {
  DataFrame plantsDF = DataFrame::create(
    _["LAI"] = Rcpp::wrap(plants.LAI),
    _["LAIlive"] = Rcpp::wrap(plants.LAIlive),
    _["FPAR"] = Rcpp::wrap(plants.FPAR),
    _["Extraction"] = Rcpp::wrap(plants.Extraction),
    _["Transpiration"] = Rcpp::wrap(plants.Transpiration),
    _["GrossPhotosynthesis"] = Rcpp::wrap(plants.GrossPhotosynthesis),
    _["NetPhotosynthesis"] = Rcpp::wrap(plants.NetPhotosynthesis),
    _["RootPsi"] = Rcpp::wrap(plants.RootPsi),
    _["StemPsi"] = Rcpp::wrap(plants.StemPsi),
    _["StemPLC"] = Rcpp::wrap(plants.StemPLC),
    _["LeafPLC"] = Rcpp::wrap(plants.LeafPLC),
    _["LeafPsiMin"] = Rcpp::wrap(plants.LeafPsiMin),
    _["LeafPsiMax"] = Rcpp::wrap(plants.LeafPsiMax),
    _["dEdP"] = Rcpp::wrap(plants.dEdP),
    _["DDS"] = Rcpp::wrap(plants.DDS),
    _["StemRWC"] = Rcpp::wrap(plants.StemRWC),
    _["LeafRWC"] = Rcpp::wrap(plants.LeafRWC),
    _["LFMC"] = Rcpp::wrap(plants.LFMC),
    _["WaterBalance"] = Rcpp::wrap(plants.WaterBalance)
  );
  plantsDF.attr("row.names") = x.cohorts.CohortCode;
  return(plantsDF);
}

Rcpp::List copyPlantAdvancedTranspirationInstResult_c(const PlantsAdvancedTranspirationInst_RESULT& plants_inst, ModelInput& x) {
  int numCohorts = x.cohorts.CohortCode.size();
  int ntimesteps = x.control.advancedWB.ndailysteps;
  
  NumericMatrix E = copyNumericMatrix_c(plants_inst.E, numCohorts, ntimesteps);
  E.attr("dimnames") = List::create(x.cohorts.CohortCode, Rcpp::seq(1,ntimesteps));
  NumericMatrix Ag = copyNumericMatrix_c(plants_inst.Ag, numCohorts, ntimesteps);
  Ag.attr("dimnames") = List::create(x.cohorts.CohortCode, Rcpp::seq(1,ntimesteps));
  NumericMatrix An = copyNumericMatrix_c(plants_inst.An, numCohorts, ntimesteps);
  An.attr("dimnames") = List::create(x.cohorts.CohortCode, Rcpp::seq(1,ntimesteps));
  NumericMatrix dEdP = copyNumericMatrix_c(plants_inst.dEdP, numCohorts, ntimesteps);
  dEdP.attr("dimnames") = List::create(x.cohorts.CohortCode, Rcpp::seq(1,ntimesteps));
  NumericMatrix LeafPsi = copyNumericMatrix_c(plants_inst.LeafPsi, numCohorts, ntimesteps);
  LeafPsi.attr("dimnames") = List::create(x.cohorts.CohortCode, Rcpp::seq(1,ntimesteps));
  NumericMatrix StemPsi = copyNumericMatrix_c(plants_inst.StemPsi, numCohorts, ntimesteps);
  StemPsi.attr("dimnames") = List::create(x.cohorts.CohortCode, Rcpp::seq(1,ntimesteps));
  NumericMatrix RootPsi = copyNumericMatrix_c(plants_inst.RootPsi, numCohorts, ntimesteps);
  RootPsi.attr("dimnames") = List::create(x.cohorts.CohortCode, Rcpp::seq(1,ntimesteps));
  NumericMatrix LeafSympPsi = copyNumericMatrix_c(plants_inst.LeafSympPsi, numCohorts, ntimesteps);
  LeafSympPsi.attr("dimnames") = List::create(x.cohorts.CohortCode, Rcpp::seq(1,ntimesteps));
  NumericMatrix StemSympPsi = copyNumericMatrix_c(plants_inst.StemSympPsi, numCohorts, ntimesteps);
  StemSympPsi.attr("dimnames") = List::create(x.cohorts.CohortCode, Rcpp::seq(1,ntimesteps));
  NumericMatrix LeafSympRWC = copyNumericMatrix_c(plants_inst.LeafSympRWC, numCohorts, ntimesteps);
  LeafSympRWC.attr("dimnames") = List::create(x.cohorts.CohortCode, Rcpp::seq(1,ntimesteps));
  NumericMatrix StemSympRWC = copyNumericMatrix_c(plants_inst.StemSympRWC, numCohorts, ntimesteps);
  StemSympRWC.attr("dimnames") = List::create(x.cohorts.CohortCode, Rcpp::seq(1,ntimesteps));
  NumericMatrix StemPLC = copyNumericMatrix_c(plants_inst.StemPLC, numCohorts, ntimesteps);
  StemPLC.attr("dimnames") = List::create(x.cohorts.CohortCode, Rcpp::seq(1,ntimesteps));
  NumericMatrix LeafPLC = copyNumericMatrix_c(plants_inst.LeafPLC, numCohorts, ntimesteps);
  LeafPLC.attr("dimnames") = List::create(x.cohorts.CohortCode, Rcpp::seq(1,ntimesteps));
  NumericMatrix LeafRWC = copyNumericMatrix_c(plants_inst.LeafRWC, numCohorts, ntimesteps);
  LeafRWC.attr("dimnames") = List::create(x.cohorts.CohortCode, Rcpp::seq(1,ntimesteps));
  NumericMatrix StemRWC = copyNumericMatrix_c(plants_inst.StemRWC, numCohorts, ntimesteps);
  StemRWC.attr("dimnames") = List::create(x.cohorts.CohortCode, Rcpp::seq(1,ntimesteps));
  NumericMatrix PWB = copyNumericMatrix_c(plants_inst.PWB, numCohorts, ntimesteps);
  PWB.attr("dimnames") = List::create(x.cohorts.CohortCode, Rcpp::seq(1,ntimesteps));
  
  List PlantsInst = List::create(
    _["E"]=E, _["Ag"]=Ag, _["An"]=An,
    _["dEdP"] = dEdP,
    _["RootPsi"] = RootPsi, 
    _["StemPsi"] = StemPsi,
    _["LeafPsi"] = LeafPsi,
    _["StemSympPsi"] = StemSympPsi,
    _["LeafSympPsi"] = LeafSympPsi,
    _["StemPLC"] = StemPLC, 
    _["LeafPLC"] = LeafPLC, 
    _["StemRWC"] = StemRWC,
    _["LeafRWC"] = LeafRWC,
    _["StemSympRWC"] = StemSympRWC,
    _["LeafSympRWC"] = LeafSympRWC,
    _["PWB"] = PWB);
  return(PlantsInst);
}
Rcpp::DataFrame copyLeafAdvancedTranspirationResult_c(const LeafAdvancedTranspiration_RESULT& leaf, ModelInput& x) {
  DataFrame leafDF = DataFrame::create(
    _["LeafPsiMin"] = Rcpp::wrap(leaf.LeafPsiMin),
    _["LeafPsiMax"] = Rcpp::wrap(leaf.LeafPsiMax),
    _["GSWMin"] = Rcpp::wrap(leaf.GSWMin),
    _["GSWMax"] = Rcpp::wrap(leaf.GSWMax),
    _["TempMin"] = Rcpp::wrap(leaf.TempMin),
    _["TempMax"] = Rcpp::wrap(leaf.TempMax)
  );
  return(leafDF);
}

Rcpp::List copyLeafAdvancedTranspirationInstResult_c(const LeafAdvancedTranspirationInst_RESULT& leaf_inst, ModelInput& x) {
  int numCohorts = x.cohorts.CohortCode.size();
  int ntimesteps = x.control.advancedWB.ndailysteps;
  
  NumericMatrix LAI = copyNumericMatrix_c(leaf_inst.LAI, numCohorts, ntimesteps);
  LAI.attr("dimnames") = List::create(x.cohorts.CohortCode, Rcpp::seq(1,ntimesteps));
  NumericMatrix Vmax298 = copyNumericMatrix_c(leaf_inst.Vmax298, numCohorts, ntimesteps);
  Vmax298.attr("dimnames") = List::create(x.cohorts.CohortCode, Rcpp::seq(1,ntimesteps));
  NumericMatrix Jmax298 = copyNumericMatrix_c(leaf_inst.Jmax298, numCohorts, ntimesteps);
  Jmax298.attr("dimnames") = List::create(x.cohorts.CohortCode, Rcpp::seq(1,ntimesteps));
  NumericMatrix Abs_SWR = copyNumericMatrix_c(leaf_inst.Abs_SWR, numCohorts, ntimesteps);
  Abs_SWR.attr("dimnames") = List::create(x.cohorts.CohortCode, Rcpp::seq(1,ntimesteps));
  NumericMatrix Abs_PAR = copyNumericMatrix_c(leaf_inst.Abs_PAR, numCohorts, ntimesteps);
  Abs_PAR.attr("dimnames") = List::create(x.cohorts.CohortCode, Rcpp::seq(1,ntimesteps));
  NumericMatrix Net_LWR = copyNumericMatrix_c(leaf_inst.Net_LWR, numCohorts, ntimesteps);
  Net_LWR.attr("dimnames") = List::create(x.cohorts.CohortCode, Rcpp::seq(1,ntimesteps));
  NumericMatrix Ag = copyNumericMatrix_c(leaf_inst.Ag, numCohorts, ntimesteps);
  Ag.attr("dimnames") = List::create(x.cohorts.CohortCode, Rcpp::seq(1,ntimesteps));
  NumericMatrix An = copyNumericMatrix_c(leaf_inst.An, numCohorts, ntimesteps);
  An.attr("dimnames") = List::create(x.cohorts.CohortCode, Rcpp::seq(1,ntimesteps));
  NumericMatrix Ci = copyNumericMatrix_c(leaf_inst.Ci, numCohorts, ntimesteps);
  Ci.attr("dimnames") = List::create(x.cohorts.CohortCode, Rcpp::seq(1,ntimesteps));
  NumericMatrix E = copyNumericMatrix_c(leaf_inst.E, numCohorts, ntimesteps);
  E.attr("dimnames") = List::create(x.cohorts.CohortCode, Rcpp::seq(1,ntimesteps));
  NumericMatrix Gsw = copyNumericMatrix_c(leaf_inst.Gsw, numCohorts, ntimesteps);
  Gsw.attr("dimnames") = List::create(x.cohorts.CohortCode, Rcpp::seq(1,ntimesteps));
  NumericMatrix VPD = copyNumericMatrix_c(leaf_inst.VPD, numCohorts, ntimesteps);
  VPD.attr("dimnames") = List::create(x.cohorts.CohortCode, Rcpp::seq(1,ntimesteps));
  NumericMatrix Temp = copyNumericMatrix_c(leaf_inst.Temp, numCohorts, ntimesteps);
  Temp.attr("dimnames") = List::create(x.cohorts.CohortCode, Rcpp::seq(1,ntimesteps));
  NumericMatrix Psi = copyNumericMatrix_c(leaf_inst.Psi, numCohorts, ntimesteps);
  Psi.attr("dimnames") = List::create(x.cohorts.CohortCode, Rcpp::seq(1,ntimesteps));
  
  List LeavesInst = List::create(
    _["LAI"] = LAI,
    _["Vmax298"] = Vmax298,
    _["Jmax298"] = Jmax298,
    _["Abs_SWR"]=Abs_SWR,
    _["Abs_PAR"]=Abs_PAR,
    _["Net_LWR"] = Net_LWR,
    _["Ag"] = Ag,
    _["An"] = An,
    _["Ci"] = Ci,
    _["E"] = E,
    _["Gsw"] = Gsw,
    _["VPD"] = VPD,
    _["Temp"] = Temp,
    _["Psi"] = Psi);
  
  return(LeavesInst);
}

Rcpp::List copyEnergyBalanceResult_c(const EnergyBalance_RESULT& EBres, ModelInput& x) {
  int ncanlayers = x.canopy.zlow.size(); //Number of canopy layers
  int ntimesteps = EBres.Ebalcan.size();
  int nlayers = x.soil.getNlayers();
  DataFrame Tinst = DataFrame::create(
    _["SolarHour"] = Rcpp::wrap(EBres.SolarHour),
    _["Tatm"] = Rcpp::wrap(EBres.Tatm),
    _["Tcan"] = Rcpp::wrap(EBres.Tcan)
  );
  DataFrame CEBinst = DataFrame::create(
    _["SolarHour"] = Rcpp::wrap(EBres.SolarHour),
    _["SWRcan"] = Rcpp::wrap(EBres.SWRcan),
    _["LWRcan"] = Rcpp::wrap(EBres.LWRcan),
    _["LEVcan"] = Rcpp::wrap(EBres.LEVcan),
    _["LEFsnow"] = Rcpp::wrap(EBres.LEFsnow),
    _["Hcan"] = Rcpp::wrap(EBres.Hcan),
    _["Ebalcan"] = Rcpp::wrap(EBres.Ebalcan)
  );
  DataFrame SEBinst = DataFrame::create(
    _["SolarHour"] = Rcpp::wrap(EBres.SolarHour),
    _["Hcansoil"] = Rcpp::wrap(EBres.Hcansoil),
    _["LEVsoil"] = Rcpp::wrap(EBres.LEVsoil),
    _["SWRsoil"] = Rcpp::wrap(EBres.SWRsoil),
    _["LEFsnow"] = Rcpp::wrap(EBres.LEFsnow),
    _["LWRsoil"] = Rcpp::wrap(EBres.LWRsoil),
    _["Ebalsoil"] = Rcpp::wrap(EBres.Ebalsoil)
  );
  NumericMatrix Tcan_mat= copyNumericMatrix_c(EBres.TemperatureLayers, ntimesteps, ncanlayers);
  Tcan_mat.attr("dimnames") = List::create(seq(1,ntimesteps), seq(1,ncanlayers));
  NumericMatrix VPcan_mat= copyNumericMatrix_c(EBres.VaporPressureLayers, ntimesteps, ncanlayers);
  VPcan_mat.attr("dimnames") = List::create(seq(1,ntimesteps), seq(1,ncanlayers));
  NumericMatrix Tsoil_mat= copyNumericMatrix_c(EBres.SoilTemperature, ntimesteps, nlayers);
  Tsoil_mat.attr("dimnames") = List::create(seq(1,ntimesteps), seq(1,nlayers));
  List EnergyBalance = List::create(_["Temperature"]=Tinst, 
                                    _["SoilTemperature"] = Tsoil_mat,
                                    _["CanopyEnergyBalance"] = CEBinst, 
                                    _["SoilEnergyBalance"] = SEBinst,
                                    _["TemperatureLayers"] = Tcan_mat, 
                                    _["VaporPressureLayers"] = VPcan_mat);
  return(EnergyBalance);
}

Rcpp::List copyAdvancedTranspirationResult_c(const AdvancedTranspiration_RESULT& ATres, ModelInput& x) {
  const std::string& rhizosphereOverlap = x.control.rhizosphereOverlap;
  bool plantWaterPools = (rhizosphereOverlap!="total");
  int nlayers = x.soil.getNlayers();
  int numCohorts = x.cohorts.CohortCode.size();
  int ntimesteps = x.control.advancedWB.ndailysteps;
  
  const arma::mat& extractionComm = ATres.extraction;
  Rcpp::NumericMatrix Extraction = copyNumericMatrix_c(extractionComm, numCohorts, nlayers); // this is final extraction of each cohort from each layer
  Extraction.attr("dimnames") = Rcpp::List::create(x.cohorts.CohortCode, Rcpp::seq(1,nlayers));
  
  Rcpp::List ExtractionPools(numCohorts);
  const std::vector< arma::mat>& ExtractionPoolsComm = ATres.extractionPools;
  if(plantWaterPools) {
    for(int c=0;c<numCohorts;c++) {
      const arma::mat& extractionPoolsCohComm = ExtractionPoolsComm[c];
      Rcpp::NumericMatrix ExtractionPoolsCohComm_c = copyNumericMatrix_c(extractionPoolsCohComm, numCohorts, nlayers); // this is final extraction of each cohort from each layer
      ExtractionPoolsCohComm_c.attr("dimnames") = Rcpp::List::create(x.cohorts.CohortCode, Rcpp::seq(1,nlayers));
      ExtractionPools[c] = ExtractionPoolsCohComm_c;
    }
    ExtractionPools.attr("names") = x.cohorts.CohortCode;
  }
  
  NumericMatrix extractionInst = copyNumericMatrix_c(ATres.extractionInst, nlayers, ntimesteps);
  extractionInst.attr("dimnames") = List::create(seq(1,nlayers), seq(1,ntimesteps));
  
  NumericMatrix rhizoPsi = copyNumericMatrix_c(ATres.rhizoPsi, numCohorts, nlayers);
  rhizoPsi.attr("dimnames") = List::create(x.cohorts.CohortCode, seq(1,nlayers));
  
  
  NumericVector standVEC = copyStandBasicTranspirationResult_c(ATres.stand);
  
  List EnergyBalance = copyEnergyBalanceResult_c(ATres.energy, x);
  
  List lwrExtinctionList(ntimesteps);
  for(int n=0;n<ntimesteps;n++) {
    lwrExtinctionList[n] = copyLongWaveRadiationResult_c(ATres.lwrExtinction[n]);
  }
  
  List l = List::create(_["cohorts"] = copyCohorts_c(x.cohorts),
                        _["EnergyBalance"] = EnergyBalance,
                        _["Extraction"] = Extraction,
                        _["ExtractionPools"] = ExtractionPools,
                        _["RhizoPsi"] = rhizoPsi,
                        _["Stand"] = standVEC,
                        _["Plants"] = copyPlantAdvancedTranspirationResult_c(ATres.plants, x),
                        _["SunlitLeaves"] = copyLeafAdvancedTranspirationResult_c(ATres.sunlit, x),
                        _["ShadeLeaves"] = copyLeafAdvancedTranspirationResult_c(ATres.shade, x),
                        _["ExtractionInst"] = extractionInst,
                        _["RadiationInputInst"] = copyDirectDiffuseDayResult_c(ATres.directDiffuseDay),
                        _["PlantsInst"] = copyPlantAdvancedTranspirationInstResult_c(ATres.plants_inst, x),
                        _["SunlitLeavesInst"] = copyLeafAdvancedTranspirationInstResult_c(ATres.sunlit_inst, x),
                        _["ShadeLeavesInst"] = copyLeafAdvancedTranspirationInstResult_c(ATres.shade_inst, x),
                        _["LightExtinction"] = copyInstantaneousLightExtinctionAbsortionResult_c(ATres.lightExtinctionAbsortion),
                        _["LWRExtinction"] = lwrExtinctionList,
                        _["CanopyTurbulence"] = copyCanopyTurbulenceResult_c(ATres.canopyTurbulence));
  
  // List supply = copyList(as<List>(atc["SupplyFunctions"]), numCohorts);
  // supply.attr("names") = above.attr("row.names");
  // List outPhotoSunlit = copyList(as<List>(atc["PhotoSunlitFunctions"]), numCohorts);
  // outPhotoSunlit.attr("names") = above.attr("row.names");
  // List outPhotoShade = copyList(as<List>(atc["PhotoShadeFunctions"]), numCohorts);
  // outPhotoShade.attr("names") = above.attr("row.names");
  // List outPMSunlit = copyList(as<List>(atc["PMSunlitFunctions"]), numCohorts);
  // outPMSunlit.attr("names") = above.attr("row.names");
  // List outPMShade = copyList(as<List>(atc["PMShadeFunctions"]), numCohorts);
  // outPMShade.attr("names") = above.attr("row.names");
  // l.push_back(supply, "SupplyFunctions");
  // l.push_back(outPhotoSunlit, "PhotoSunlitFunctions");
  // l.push_back(outPhotoShade, "PhotoShadeFunctions");
  // l.push_back(outPMSunlit, "PMSunlitFunctions");
  // l.push_back(outPMShade, "PMShadeFunctions");
  l.attr("class") = CharacterVector::create("pwb_day","list");
  return(l);
}


void transpirationAdvanced_c(AdvancedTranspiration_RESULT& ATres, AdvancedTranspiration_COMM& ATcomm, ModelInput& x, 
                             const WeatherInputVector& meteovec,
                             const double latitude, double elevation, double slope, double aspect, 
                             const double solarConstant, const double delta,
                             const double canopyEvaporation = 0.0, const double snowMelt = 0.0, const double soilEvaporation = 0.0, const double herbTranspiration = 0.0, 
                             int stepFunctions = -1) {
  
  // Should have internal communication structures for output
  StandBasicTranspiration_RESULT& outputStand = ATres.stand;
  PlantsAdvancedTranspiration_RESULT& outputPlants = ATres.plants;
  PlantsAdvancedTranspirationInst_RESULT& outputPlantsInst = ATres.plants_inst;
  arma::mat& outputExtraction = ATres.extraction;
  LeafAdvancedTranspiration_RESULT& outputSunlit = ATres.sunlit;
  LeafAdvancedTranspiration_RESULT& outputShade = ATres.shade;
  LeafAdvancedTranspirationInst_RESULT& outputSunlitInst = ATres.sunlit_inst;
  LeafAdvancedTranspirationInst_RESULT& outputShadeInst = ATres.shade_inst;
  EnergyBalance_RESULT& outputEnergyBalance = ATres.energy;
  arma::mat& RHOPCohDyn = ATcomm.RHOPCohDyn;
  
  // DataFrame outputTemperatureInst =   as<DataFrame>(outputEnergyBalance["Temperature"]);
  // DataFrame outputCEBinst =  as<DataFrame>(outputEnergyBalance["CanopyEnergyBalance"]);
  // DataFrame outputSEBinst =  as<DataFrame>(outputEnergyBalance["SoilEnergyBalance"]);
  
  // List supply = transpOutput["SupplyFunctions"]; 

  //Control parameters
  std::string& transpirationMode = x.control.transpirationMode;
  int ntimesteps = x.control.advancedWB.ndailysteps;
  int nsubsteps = x.control.advancedWB.nsubsteps_canopy;
  std::string& rhizosphereOverlap = x.control.rhizosphereOverlap;
  double fullRhizosphereOverlapConductivity = x.control.commonWB.fullRhizosphereOverlapConductivity;
  bool plantWaterPools = (rhizosphereOverlap!="total");
  double verticalLayerSize = x.control.commonWB.verticalLayerSize;
  double windMeasurementHeight  = x.control.weather.windMeasurementHeight;
  double defaultWindSpeed = x.control.weather.defaultWindSpeed;
  double thermalCapacityLAI = x.control.advancedWB.thermalCapacityLAI;
  bool multiLayerBalance = x.control.advancedWB.multiLayerBalance;
  std::string& stemCavitationRecovery = x.control.commonWB.stemCavitationRecovery;
  std::string& leafCavitationRecovery = x.control.commonWB.leafCavitationRecovery;
  double cavitationRecoveryMaximumRate = x.control.commonWB.cavitationRecoveryMaximumRate;
  bool sapFluidityVariation = x.control.advancedWB.sapFluidityVariation;
  std::string& lfmcComponent = x.control.fireHazard.lfmcComponent;

  //Meteo input
  double pet = meteovec.pet;
  double rhmax = meteovec.rhmax;
  double rhmin = meteovec.rhmin;
  double tmax = meteovec.tmax;
  double tmin = meteovec.tmin;
  double tminPrev = meteovec.tminPrev;
  double tmaxPrev = meteovec.tmaxPrev;
  double tminNext = meteovec.tminNext;
  double Catm = meteovec.Catm;
  double Patm = meteovec.Patm;
  double prec = meteovec.prec;
  double rad = meteovec.rad;
  double wind = meteovec.wind;

  //Atmospheric pressure (if missing)
  if(std::isnan(Patm)) Patm = atmosphericPressure_c(elevation);

  //Vegetation input
  std::vector<double>& LAIlive = x.above.LAI_live;
  std::vector<double>& LAIphe = x.above.LAI_expanded;
  std::vector<double>& LAIdead = x.above.LAI_dead;
  std::vector<double>& H = x.above.H;
  std::vector<double>& CR = x.above.CR;
  
  int numCohorts = LAIlive.size();

  //Soil input
  Soil& soil = x.soil;
  int nlayers = x.soil.getNlayers();
  // NumericVector widths = soil["widths"];
  // NumericVector Water_FC = waterFC(soil, soilFunctions);
  // NumericVector Theta_FC = thetaFC(soil, soilFunctions);
  // NumericVector Theta_SAT = thetaSAT(soil, soilFunctions);
  // NumericVector VG_n = Rcpp::as<Rcpp::NumericVector>(soil["VG_n"]);
  // NumericVector VG_alpha = Rcpp::as<Rcpp::NumericVector>(soil["VG_alpha"]);
  // NumericVector Tsoil = soil["Temp"]; 
  // NumericVector sand = soil["sand"];
  // NumericVector clay = soil["clay"];
  // NumericVector Ws = soil["W"]; //Access to soil state variable
  // double snowpack = x["snowpack"];

  //Canopy params
  int ncanlayers = x.canopy.zlow.size();
  std::vector<double>& zlow = x.canopy.zlow;
  std::vector<double>& zmid = x.canopy.zmid;
  std::vector<double>& zup = x.canopy.zup;
  std::vector<double>& LAIpx = x.canopy.LAIlive;
  std::vector<double>& LAIpe = x.canopy.LAIexpanded;
  std::vector<double>& LAIpd = x.canopy.LAIdead;
  std::vector<double>& Tair = x.canopy.Tair;
  std::vector<double>& VPair = x.canopy.VPair;
  std::vector<double>& Cair = x.canopy.Cair;
  for(int l=0;l<ncanlayers;l++) { //If canopy layers have missing values, then initialize with Catm
    if(!multiLayerBalance) Cair[l] = Catm;
    else {
      if(NumericVector::is_na(Cair[l])) Cair[l] = Catm;
    }
  }
  if(multiLayerBalance) Cair[ncanlayers-1] = Catm;

  // //Root distribution input
  // DataFrame belowdf = Rcpp::as<Rcpp::DataFrame>(x["below"]);
  // List belowLayers = Rcpp::as<Rcpp::List>(x["belowLayers"]);
  // NumericMatrix V =Rcpp::as<Rcpp::NumericMatrix>(belowLayers["V"]);
  // NumericMatrix VCroot_kmax= Rcpp::as<Rcpp::NumericMatrix>(belowLayers["VCroot_kmax"]);
  // NumericMatrix VGrhizo_kmax= Rcpp::as<Rcpp::NumericMatrix>(belowLayers["VGrhizo_kmax"]);
  // NumericMatrix RhizoPsiMAT = Rcpp::as<Rcpp::NumericMatrix>(belowLayers["RhizoPsi"]);
  // 
  // //Water pools
  // NumericMatrix Wpool = Rcpp::as<Rcpp::NumericMatrix>(belowLayers["Wpool"]);
  // List RHOP;
  // NumericVector poolProportions(numCohorts);
  // if(plantWaterPools) {
  //   RHOP = belowLayers["RHOP"];
  //   poolProportions = belowdf["poolProportions"];
  // }
  // 
  //Base parameters
  // 
  // //Phenology parameters
  // DataFrame paramsPhenology = Rcpp::as<Rcpp::DataFrame>(x["paramsPhenology"]);
  // CharacterVector phenoType = paramsPhenology["PhenologyType"];
  // 
  // //Anatomy parameters
  // DataFrame paramsAnatomy = Rcpp::as<Rcpp::DataFrame>(x["paramsAnatomy"]);
  // NumericVector leafWidth = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["LeafWidth"]);
  // NumericVector Al2As = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["Al2As"]);
  // NumericVector r635 = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["r635"]);
  // 
  // //Transpiration parameters
  // DataFrame paramsTranspiration = Rcpp::as<Rcpp::DataFrame>(x["paramsTranspiration"]);
  // NumericVector Plant_kmax = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["Plant_kmax"]);
  // NumericVector VCstem_kmax = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["VCstem_kmax"]);
  // NumericVector VCleaf_kmax = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["VCleaf_kmax"]);
  // 
  // NumericVector Vmax298 = paramsTranspiration["Vmax298"];
  // NumericVector Jmax298 = paramsTranspiration["Jmax298"];
  // 
  // //Water storage parameters
  // DataFrame paramsWaterStorage = Rcpp::as<Rcpp::DataFrame>(x["paramsWaterStorage"]);
  // NumericVector maxFMC = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["maxFMC"]);
  // NumericVector maxMCstem = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["maxMCstem"]);
  // NumericVector maxMCleaf = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["maxMCleaf"]);
  // NumericVector StemPI0 = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["StemPI0"]);
  // NumericVector StemEPS = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["StemEPS"]);
  // NumericVector StemAF = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["StemAF"]);
  // NumericVector Vsapwood = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["Vsapwood"]); //l·m-2 = mm
  // NumericVector LeafPI0 = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["LeafPI0"]);
  // NumericVector LeafEPS = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["LeafEPS"]);
  // NumericVector LeafAF = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["LeafAF"]);
  // NumericVector Vleaf = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["Vleaf"]); //l·m-2 = mm
  // 
  // //Comunication with outside
  // DataFrame internalPhenology = Rcpp::as<Rcpp::DataFrame>(x["internalPhenology"]);
  // NumericVector phi = Rcpp::as<Rcpp::NumericVector>(internalPhenology["phi"]);
  // DataFrame internalWater = Rcpp::as<Rcpp::DataFrame>(x["internalWater"]);
  // NumericVector LeafPLCVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["LeafPLC"]);
  // NumericVector StemPLCVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["StemPLC"]);
  // NumericVector RootCrownPsiVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["RootCrownPsi"]);
  // NumericVector StemPsiVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["StemPsi"]);
  // NumericVector LeafPsiVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["LeafPsi"]);
  // NumericVector StemSympPsiVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["StemSympPsi"]);
  // NumericVector LeafSympPsiVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["LeafSympPsi"]);
  
  arma::mat& LAImx = x.internalLAIDistribution.live;
  arma::mat& LAIme = x.internalLAIDistribution.expanded;
  arma::mat& LAImd = x.internalLAIDistribution.dead;
  std::vector<double>& PrevLAIexpanded = x.internalLAIDistribution.PrevLAIexpanded;
  std::vector<double>& PrevLAIdead = x.internalLAIDistribution.PrevLAIdead;
  std::vector<double>& PARcohort = x.internalLAIDistribution.PARcohort;

  if(std::isnan(aspect)) aspect = 0.0;
  if(std::isnan(slope)) slope = 0.0;
  double latrad = latitude * (M_PI/180.0);
  double asprad = aspect * (M_PI/180.0);
  double slorad = slope * (M_PI/180.0);

  //Step in seconds
  double tstep = 86400.0/((double) ntimesteps);

  //Daily average water vapor pressure at the atmosphere (kPa)
  double vpatm = averageDailyVapourPressure_c(tmin, tmax, rhmin,rhmax);
  //If canopy VP is missing or not multilayer initiate it to vpatm
  if(std::isnan(VPair[0]) || (!multiLayerBalance)){
    for(int i=0;i<ncanlayers;i++) VPair[i] = vpatm;
  }
  //Daily cloud cover
  double cloudcover = 0.0;
  if(prec >0.0) cloudcover = 1.0;
  bool clearday = (prec==0);

  ////////////////////////////////////////
  // CREATE COMMUNICATION INPUT OBJECT
  ////////////////////////////////////////
  InnerTranspirationInput_COMM input(numCohorts, nlayers, ncanlayers);
  input.Patm = Patm;
  
  ////////////////////////////////////////
  // INITIAL SOIL STATE (from previous step)
  ////////////////////////////////////////
  for(int l = 0; l<nlayers;l++) {
    input.psiSoil[l] = x.soil.getPsi(l);
  }
  if(plantWaterPools){
    //Store overall soil moisture in a backup copy
    double* Wbackup = new double[nlayers];
    for(int l = 0; l<nlayers;l++) Wbackup[l] = soil.getW(l);
    for(int j = 0; j<numCohorts;j++) {
      //Copy values of soil moisture from pool of cohort j to general soil
      for(int l = 0; l<nlayers;l++) {
        soil.setW(l,x.belowLayers.Wpool(j,l)); // this updates psi, theta, ... 
        input.psiSoilM(j,l) = soil.getPsi(l);
        input.KunsatM(j,l) = soil.getConductivity(l, true);
      }
    }
    //Restore soil moisture
    for(int l = 0; l<nlayers;l++) soil.setW(l, Wbackup[l]);
    //Delete backup
    delete[] Wbackup;
  }

  // ////////////////////////////////////////
  // // DEFINE OUTPUT
  // ////////////////////////////////////////
  // //Transpiration and photosynthesis
  // NumericMatrix minPsiRhizo = transpOutput["RhizoPsi"];
  // 
  // NumericVector outputFPAR = outputPlants["FPAR"];
  // arma::mat& SoilExtractCoh = ATres.extraction;
  // NumericVector DDS = outputPlants["DDS"];
  // NumericVector LFMC = outputPlants["LFMC"];
  // NumericVector Eplant = outputPlants["Transpiration"];
  // NumericVector Anplant = outputPlants["GrossPhotosynthesis"];
  // NumericVector Agplant = outputPlants["NetPhotosynthesis"];
  // NumericVector minStemPsi= outputPlants["StemPsi"];
  // NumericVector minRootPsi= outputPlants["RootPsi"];
  // NumericVector minLeafPsi= outputPlants["LeafPsiMin"];
  // NumericVector maxLeafPsi= outputPlants["LeafPsiMax"];
  // NumericVector PLClm = outputPlants["LeafPLC"];
  // NumericVector PLCsm = outputPlants["StemPLC"];
  // NumericVector dEdPm = outputPlants["dEdP"];
  // NumericVector PWB = outputPlants["WaterBalance"];
  // NumericVector RWCsm = outputPlants["StemRWC"];
  // NumericVector RWClm = outputPlants["LeafRWC"];
  // 
  // 
  // NumericVector maxGSW_SL = outputSunlit["GSWMax"];
  // NumericVector maxGSW_SH = outputShade["GSWMax"];
  // NumericVector minGSW_SL = outputSunlit["GSWMin"];
  // NumericVector minGSW_SH = outputShade["GSWMin"];
  // NumericVector maxTemp_SL = outputSunlit["TempMax"];
  // NumericVector maxTemp_SH = outputShade["TempMax"];
  // NumericVector minTemp_SL = outputSunlit["TempMin"];
  // NumericVector minTemp_SH = outputShade["TempMin"];
  // NumericVector maxLeafPsi_SL = outputSunlit["LeafPsiMax"];
  // NumericVector maxLeafPsi_SH = outputShade["LeafPsiMax"];
  // NumericVector minLeafPsi_SL = outputSunlit["LeafPsiMin"];
  // NumericVector minLeafPsi_SH = outputShade["LeafPsiMin"];
  // 
  // NumericMatrix Einst = outputPlantsInst["E"];
  // NumericMatrix Aginst = outputPlantsInst["Ag"];
  // NumericMatrix Aninst = outputPlantsInst["An"];
  // NumericMatrix dEdPInst = outputPlantsInst["dEdP"];
  // NumericMatrix LeafPsiInst = outputPlantsInst["LeafPsi"];
  // NumericMatrix StemPsiInst = outputPlantsInst["StemPsi"];
  // NumericMatrix RootPsiInst = outputPlantsInst["RootPsi"];
  // NumericMatrix LeafSympPsiInst = outputPlantsInst["LeafSympPsi"];
  // NumericMatrix StemSympPsiInst = outputPlantsInst["StemSympPsi"];
  // NumericMatrix StemPLC = outputPlantsInst["StemPLC"];
  // NumericMatrix LeafPLC = outputPlantsInst["LeafPLC"];
  // NumericMatrix LeafRWCInst = outputPlantsInst["LeafRWC"];
  // NumericMatrix StemRWCInst = outputPlantsInst["StemRWC"];
  // NumericMatrix LeafSympRWCInst = outputPlantsInst["LeafSympRWC"];
  // NumericMatrix StemSympRWCInst = outputPlantsInst["StemSympRWC"];
  // NumericMatrix PWBinst = outputPlantsInst["PWB"];
  // 
  arma::mat& LAI_SL = outputSunlitInst.LAI;
  arma::mat& LAI_SH = outputShadeInst.LAI;
  arma::mat& Vmax298_SL = outputSunlitInst.Vmax298;
  arma::mat& Vmax298_SH = outputShadeInst.Vmax298;
  arma::mat& Jmax298_SL = outputSunlitInst.Jmax298;
  arma::mat& Jmax298_SH = outputShadeInst.Jmax298;
  // NumericMatrix SWR_SL = outputSunlitInst["Abs_SWR"];
  // NumericMatrix SWR_SH = outputShadeInst["Abs_SWR"];
  // NumericMatrix PAR_SL = outputSunlitInst["Abs_PAR"];
  // NumericMatrix PAR_SH = outputShadeInst["Abs_PAR"];
  // NumericMatrix LWR_SL = outputSunlitInst["Net_LWR"];
  // NumericMatrix LWR_SH = outputShadeInst["Net_LWR"];
  // NumericMatrix An_SL = outputSunlitInst["An"];
  // NumericMatrix An_SH = outputShadeInst["An"];
  // NumericMatrix Ag_SL = outputSunlitInst["Ag"];
  // NumericMatrix Ag_SH = outputShadeInst["Ag"];
  // NumericMatrix Ci_SL = outputSunlitInst["Ci"];
  // NumericMatrix Ci_SH = outputShadeInst["Ci"];
  // NumericMatrix E_SL = outputSunlitInst["E"];
  // NumericMatrix E_SH = outputShadeInst["E"];
  // NumericMatrix GSW_SL = outputSunlitInst["Gsw"];
  // NumericMatrix GSW_SH = outputShadeInst["Gsw"];
  // NumericMatrix VPD_SL = outputSunlitInst["VPD"];
  // NumericMatrix VPD_SH = outputShadeInst["VPD"];
  // NumericMatrix Temp_SL = outputSunlitInst["Temp"];
  // NumericMatrix Temp_SH = outputShadeInst["Temp"];
  // NumericMatrix Psi_SL = outputSunlitInst["Psi"];
  // NumericMatrix Psi_SH = outputShadeInst["Psi"];
  // 
  // 
  //Reset output data
  std::fill(outputPlants.Transpiration.begin(), outputPlants.Transpiration.end(), 0.0);
  std::fill(outputPlants.GrossPhotosynthesis.begin(), outputPlants.GrossPhotosynthesis.end(), 0.0);
  std::fill(outputPlants.NetPhotosynthesis.begin(), outputPlants.NetPhotosynthesis.end(), 0.0);
  std::fill(outputPlants.WaterBalance.begin(), outputPlants.WaterBalance.end(), 0.0);
  std::fill(outputPlants.Extraction.begin(), outputPlants.Extraction.end(), 0.0);
  std::fill(outputExtraction.begin(), outputExtraction.end(), 0.0);
  std::fill(ATres.extractionInst.begin(), ATres.extractionInst.end(), 0.0);
  std::fill(outputPlantsInst.E.begin(), outputPlantsInst.E.end(), 0.0);
  if(plantWaterPools) {
    for(int c=0;c<numCohorts;c++) {
      std::fill(ATres.extractionPools[c].begin(), ATres.extractionPools[c].end(), 0.0);
    }
  }
  std::fill(Jmax298_SH.begin(), Jmax298_SH.end(), 0.0);
  std::fill(Jmax298_SL.begin(), Jmax298_SL.end(), 0.0);
  std::fill(Vmax298_SH.begin(), Vmax298_SH.end(), 0.0);
  std::fill(Vmax298_SL.begin(), Vmax298_SL.end(), 0.0);
  std::fill(LAI_SH.begin(), LAI_SH.end(), 0.0);
  std::fill(LAI_SL.begin(), LAI_SL.end(), 0.0);
  
  // std::fill(StemPLC.begin(), StemPLC.end(), NA_REAL);
  // std::fill(LeafPLC.begin(), LeafPLC.end(), NA_REAL);
  // std::fill(dEdPInst.begin(), dEdPInst.end(), NA_REAL);
  // std::fill(LeafPsiInst.begin(), LeafPsiInst.end(), NA_REAL);
  // std::fill(StemPsiInst.begin(), StemPsiInst.end(), NA_REAL);
  // std::fill(RootPsiInst.begin(), RootPsiInst.end(), NA_REAL);
  // std::fill(LeafSympPsiInst.begin(), LeafSympPsiInst.end(), NA_REAL);
  // std::fill(StemSympPsiInst.begin(), StemSympPsiInst.end(), NA_REAL);
  // std::fill(Aninst.begin(), Aninst.end(), NA_REAL);
  // std::fill(Aginst.begin(), Aginst.end(), NA_REAL);
  // std::fill(LeafRWCInst.begin(), LeafRWCInst.end(), NA_REAL);
  // std::fill(StemRWCInst.begin(), StemRWCInst.end(), NA_REAL);
  // std::fill(PWBinst.begin(), PWBinst.end(), NA_REAL);
  // std::fill(LeafSympRWCInst.begin(), LeafSympRWCInst.end(), NA_REAL);
  // std::fill(StemSympRWCInst.begin(), StemSympRWCInst.end(), NA_REAL);
  // std::fill(Psi_SH.begin(), Psi_SH.end(), NA_REAL);
  // std::fill(Psi_SL.begin(), Psi_SL.end(), NA_REAL);
  // std::fill(Temp_SL.begin(), Temp_SL.end(), NA_REAL);
  // std::fill(Temp_SH.begin(), Temp_SH.end(), NA_REAL);
  // std::fill(VPD_SL.begin(), VPD_SL.end(), NA_REAL);
  // std::fill(VPD_SH.begin(), VPD_SH.end(), NA_REAL);
  // std::fill(GSW_SL.begin(), GSW_SL.end(), NA_REAL);
  // std::fill(GSW_SH.begin(), GSW_SH.end(), NA_REAL);
  // std::fill(E_SH.begin(), E_SH.end(), NA_REAL);
  // std::fill(E_SL.begin(), E_SL.end(), NA_REAL);
  // std::fill(Ci_SH.begin(), Ci_SH.end(), NA_REAL);
  // std::fill(Ci_SL.begin(), Ci_SL.end(), NA_REAL);
  // std::fill(An_SH.begin(), An_SH.end(), NA_REAL);
  // std::fill(Ag_SH.begin(), Ag_SH.end(), NA_REAL);
  // std::fill(An_SL.begin(), An_SL.end(), NA_REAL);
  // std::fill(Ag_SL.begin(), Ag_SL.end(), NA_REAL);
  // std::fill(LWR_SL.begin(), LWR_SL.end(), NA_REAL);
  // std::fill(PAR_SH.begin(), PAR_SH.end(), NA_REAL);
  // std::fill(LWR_SH.begin(), LWR_SH.end(), NA_REAL);
  // std::fill(PAR_SL.begin(), PAR_SL.end(), NA_REAL);
  // std::fill(SWR_SL.begin(), SWR_SL.end(), NA_REAL);
  // std::fill(SWR_SH.begin(), SWR_SH.end(), NA_REAL);


  ////////////////////////////////////////
  // STEP 1. Estimate stand-level leaf area values and leaf distribution across layers from leaf-level live/expanded area
  ////////////////////////////////////////
  double LAIcell = 0.0, LAIcelldead = 0.0, LAIcelllive = 0.0, LAIcellexpanded = 0.0;
  double sum_abs_exp = 0.0, sum_abs_dead = 0.0;
  double canopyHeight = 100.0; //Minimum canopy height of 1 m
  for(int c=0;c<numCohorts;c++) {
    LAIcell += (LAIphe[c]+LAIdead[c]);
    LAIcelldead += LAIdead[c];
    LAIcelllive += LAIlive[c];
    LAIcellexpanded +=LAIphe[c];
    sum_abs_exp += std::abs(LAIphe[c] - PrevLAIexpanded[c]);
    sum_abs_dead += std::abs(LAIdead[c] - PrevLAIdead[c]);
    if((canopyHeight<H[c]) && ((LAIphe[c]+LAIdead[c])>0.0)) canopyHeight = H[c];
  }

  std::vector<double> lad(ncanlayers,0.0);
  if(numCohorts>0) {
    bool recalc_LAI = false;
    if(std::isnan(PrevLAIexpanded[0]) || std::isnan(PrevLAIdead[0])) {
      recalc_LAI = true;
    } else{
      if(sum_abs_exp>0.001) {
        recalc_LAI = true;
      } else {
        if(sum_abs_dead>0.001) recalc_LAI = true;
      }
    }
    if(recalc_LAI) {
      std::vector<double> z(ncanlayers+1,0.0);
      for(int i=1;i<=ncanlayers;i++) z[i] = z[i-1] + verticalLayerSize;
      for(int i=0; i<numCohorts;i++) {
        PARcohort[i] = availableLight_c(H[i]*(1.0-(1.0-CR[i])/2.0), H, LAIphe, LAIdead, x.paramsInterception.kPAR, CR);
        PrevLAIexpanded[i] = LAIphe[i];
        PrevLAIdead[i] = LAIdead[i];
      }
      //Update LAI distribution if necessary
      updateLAIdistributionVectors_c(LAIme, z, LAIphe, H, CR);
      updateLAIdistributionVectors_c(LAImd, z, LAIdead, H, CR);
      updateLAIdistributionVectors_c(LAImx, z, LAIlive, H, CR);//Maximum leaf expansion
      //Update LAI profile per layer
      for(int i=0;i<ncanlayers;i++) {
        LAIpx[i] = LAIpd[i] = LAIpe[i] = 0.0;
        for(int j=0;j<numCohorts;j++) {
          LAIpx[i] += LAImx(i,j);
          LAIpe[i] += LAIme(i,j);
          LAIpd[i] += LAImd(i,j);
        }
        // Rcout<< i << " " << LAIpx[i] << " " << LAIpe[i] <<" "<< LAImx[i]<<"\n";
      }
    }
    // Add LAImax to leaf area density to have a wind speed profile in deciduous canopies
    for(int i=0;i<ncanlayers;i++) lad[i] = 100.0*((0.9*LAIpe[i] + 0.1*LAIpx[i]) + LAIpd[i])/verticalLayerSize;
    for(int i=0; i<numCohorts;i++) {
      outputPlants.FPAR[i] = PARcohort[i];
    }
  }
  
  
  ////////////////////////////////////////
  // STEP 2. Determine vertical wind speed profile
  ////////////////////////////////////////
  if(std::isnan(wind)) wind = defaultWindSpeed; //set to default if missing
  wind = std::min(10.0, std::max(wind, 0.1)); //Bound between 0.1 m/s (0.36 km/h)  and 10 m/s (36 km/h)
  std::vector<double> dU(ncanlayers, 0.0), uw(ncanlayers, 0.0);
  if(canopyHeight>0.0) {
    windCanopyTurbulence_inner_c(ATres.canopyTurbulence, ATcomm.canopyTurbulenceModel,
                                 zmid, lad,  
                                 canopyHeight,
                                 wind, windMeasurementHeight,
                                 "k-epsilon");
    for(int i=0;i<ncanlayers;i++) input.zWind[i] = ATres.canopyTurbulence.u[i];
    dU = ATres.canopyTurbulence.du;
    uw = ATres.canopyTurbulence.uw;
  }

  ////////////////////////////////////////
  // STEP 3a. Direct and diffuse shorwave radiation for sub-steps
  ////////////////////////////////////////
  directDiffuseDay_c(ATres.directDiffuseDay, solarConstant, latrad, slorad, asprad, delta,
                     rad, clearday);
  
  ////////////////////////////////////////
  // STEP 3b. Above-canopy air temperature and long-wave radiation emission for sub-steps
  ////////////////////////////////////////
  // NumericVector solarHour = outputTemperatureInst["SolarHour"];
  // NumericVector Hcansoil = outputSEBinst["Hcansoil"];
  // NumericVector Ebalsoil = outputSEBinst["Ebalsoil"];
  // NumericVector LEVsoil = outputSEBinst["LEVsoil"];
  // NumericVector abs_SWR_soil = outputSEBinst["SWRsoil"];
  // NumericVector net_LWR_soil = outputSEBinst["LWRsoil"];
  // NumericVector LEFsnow = outputCEBinst["LEFsnow"];
  // NumericVector abs_SWR_can = outputCEBinst["SWRcan"];
  // NumericVector net_LWR_can = outputCEBinst["LWRcan"];
  // NumericVector LEVcan = outputCEBinst["LEVcan"];
  // NumericVector Hcan_heat = outputCEBinst["Hcan"];
  // NumericVector Ebal = outputCEBinst["Ebalcan"];
  // NumericMatrix Tcan_mat = outputEnergyBalance["TemperatureLayers"];
  // NumericMatrix VPcan_mat = outputEnergyBalance["VaporPressureLayers"];
  // NumericMatrix Tsoil_mat = outputEnergyBalance["SoilTemperature"];


  std::vector<double>& Tatm = outputEnergyBalance.Tatm;
  std::vector<double>& Tcan = outputEnergyBalance.Tcan;
  std::vector<double> lwdr(ntimesteps, medfate::NA_DOUBLE);
  //Daylength in seconds (assuming flat area because we want to model air temperature variation)
  double tauday = daylengthseconds_c(latrad,0.0,0.0, delta);
  for(int n=0;n<ntimesteps;n++) {
    outputEnergyBalance.SolarHour[n] = ATres.directDiffuseDay.SolarHour[n];
    //From solar hour (radians) to seconds from sunrise
    double Tsunrise = (ATres.directDiffuseDay.SolarHour[n]*43200.0/M_PI)+ (tauday/2.0) +(tstep/2.0);
    //Calculate instantaneous temperature and light conditions
    Tatm[n] = temperatureDiurnalPattern_c(Tsunrise, tmin, tmax, tminPrev, tmaxPrev, tminNext, tauday);
    //Longwave sky diffuse radiation (W/m2)
    lwdr[n] = skyLongwaveRadiation_c(Tatm[n], vpatm, cloudcover);
  }
  if(std::isnan(Tair[0])) {//If missing initialize canopy profile with atmospheric air temperature
    for(int i=0;i<ncanlayers;i++) Tair[i] = Tatm[0];
  }
  if(std::isnan(x.soil.getTemp(0))) {//If missing initialize soil temperature with atmospheric air temperature
    for(int l=0;l<nlayers; l++) x.soil.setTemp(l, Tatm[0]);
  }
  //Take initial canopy air temperature from previous day
  double numSum = 0.0;
  double denSum = 0.0;
  for(int i=0;i<ncanlayers;i++) {
    numSum +=Tair[i]*LAIpx[i];
    denSum +=LAIpx[i];
  }
  Tcan[0] = numSum/denSum;
  // Rcout <<ncanlayers << " "<< numSum<< " " <<  denSum << " "<< Tcan[0] << "\n";
  for(int j=0;j<ncanlayers; j++) {
    outputEnergyBalance.TemperatureLayers(0,j) = Tair[j];
    outputEnergyBalance.VaporPressureLayers(0,j) = VPair[j];
  }
  //Take temperature soil vector
  for(int l=0;l<nlayers;l++) {
    outputEnergyBalance.SoilTemperature(0,l) = x.soil.getTemp(l);
  }

  ////////////////////////////////////////
  // STEP 3c. Short-wave radiation extinction and absortion for sub-steps
  ////////////////////////////////////////
  instantaneousLightExtinctionAbsortion_c(ATres.lightExtinctionAbsortion,
                                          LAIme, LAImd, LAImx,
                                          x.paramsInterception.Beta_p, x.paramsInterception.Beta_q, x.paramsInterception.ClumpingIndex,
                                          x.paramsInterception.alphaSWR, x.paramsInterception.gammaSWR,
                                          ATres.directDiffuseDay, ntimesteps, 0.1);

  std::vector<std::vector<double>>& abs_PAR_SL_COH_list = ATres.lightExtinctionAbsortion.sunshade.PAR_SL;
  std::vector<std::vector<double>>& abs_PAR_SH_COH_list = ATres.lightExtinctionAbsortion.sunshade.PAR_SH;
  std::vector<std::vector<double>>& abs_SWR_SL_COH_list = ATres.lightExtinctionAbsortion.sunshade.SWR_SL;
  std::vector<std::vector<double>>& abs_SWR_SH_COH_list = ATres.lightExtinctionAbsortion.sunshade.SWR_SH;
  std::vector<arma::mat>& abs_SWR_SL_ML_list = ATres.lightExtinctionAbsortion.multilayer.SWR_SL;
  std::vector<arma::mat>& abs_SWR_SH_ML_list = ATres.lightExtinctionAbsortion.multilayer.SWR_SH;
  std::vector<std::vector<double>>& fsunlit_list = ATres.lightExtinctionAbsortion.fsunlit;
    
  //Copy to output data structures
  double sum_abs_SWR_soil = 0.0, sum_abs_SWR_can = 0.0;
  for(int n=0; n<ntimesteps;n++) {
    outputEnergyBalance.SWRcan[n] = ATres.lightExtinctionAbsortion.SWR_can[n];
    sum_abs_SWR_can += ATres.lightExtinctionAbsortion.SWR_can[n];
    outputEnergyBalance.SWRsoil[n] = ATres.lightExtinctionAbsortion.SWR_soil[n];
    sum_abs_SWR_soil += ATres.lightExtinctionAbsortion.SWR_soil[n];
  }

  ////////////////////////////////////////
  //  STEP 4. Hydraulics: determine layers where the plant is connected
  //          and supply functions (Sperry transpiration mode)
  ////////////////////////////////////////

  //Average sap fluidity
  double sapFluidityDay = 1.0;
  if(sapFluidityVariation) sapFluidityDay = 1.0/waterDynamicViscosity_c((tmin+tmax)/2.0);

  
  //Define inner networks
  SureauNetwork* sureauNetworks = new SureauNetwork[numCohorts];
  SperryNetwork* sperryNetworks = new SperryNetwork[numCohorts];
  
  //Hydraulics: Define supply functions
  // List supplyAboveground(numCohorts);
  for(int c=0;c<numCohorts;c++) {
    if(!plantWaterPools) {
      //Determine connected layers (non-zero fine root abundance)
      input.nlayerscon[c] = 0;
      for(int l=0;l<nlayers;l++) {
        if(x.belowLayers.V(c,l)>0.0) {
          input.layerConnected(c,l)= 1;
          input.nlayerscon[c]=input.nlayerscon[c]+1;
        } else {
          input.layerConnected(c,l) = 0;
        }
      }
      // Rcout<<c<<" "<< nlayerscon[c]<<"\n";
      if(input.nlayerscon[c]==0) throw medfate::MedfateInternalError("Plant cohort not connected to any soil layer!");
      
      // Copy values from connected layers
      std::vector<double> Vc(input.nlayerscon[c]);
      std::vector<double> VCroot_kmaxc(input.nlayerscon[c]);
      std::vector<double> VGrhizo_kmaxc(input.nlayerscon[c]);
      std::vector<double> psic(input.nlayerscon[c]);
      std::vector<double> VG_nc(input.nlayerscon[c]);
      std::vector<double> VG_alphac(input.nlayerscon[c]);
      int cnt=0;
      for(int l=0;l<nlayers;l++) {
        if(input.layerConnected(c,l)==1) {
          Vc[cnt] = x.belowLayers.V(c,l);
          VCroot_kmaxc[cnt] = x.belowLayers.VCroot_kmax(c,l);
          VGrhizo_kmaxc[cnt] = x.belowLayers.VGrhizo_kmax(c,l);
          psic[cnt] = x.soil.getPsi(l);
          VG_nc[cnt] = x.soil.getVG_n(l);
          VG_alphac[cnt] = x.soil.getVG_alpha(l);
          cnt++;
        }
      }
      
      //Build supply function networks (Sperry transpiration mode)
      if(transpirationMode=="Sperry") {
        initSperryNetwork_inner_c(sperryNetworks[c], c,
                                  x.internalWater, x.paramsTranspiration, x.paramsWaterStorage,
                                  VCroot_kmaxc, VGrhizo_kmaxc,
                                  psic, VG_nc, VG_alphac,
                                  x.control,
                                  sapFluidityDay);
        //       supply[c] = supplyFunctionNetwork(HN, 0.0, 0.001); 
      } else if(transpirationMode == "Sureau") {
        initSureauNetwork_inner_c(sureauNetworks[c], c, LAIphe,
                                  x.internalWater,
                                  x.paramsAnatomy, x.paramsTranspiration, x.paramsWaterStorage,
                                  VCroot_kmaxc, VGrhizo_kmaxc,
                                  psic, VG_nc, VG_alphac,
                                  x.control, sapFluidityDay);
      }

    } else {
      //Calculate dynamic overlap  
      arma::mat& RHOPcoh = x.belowLayers.RHOP[c];
      for(int l=0;l<nlayers;l++) {
        RHOPCohDyn(c,l) = RHOPcoh(c,l);
        for(int j=0; j<numCohorts;j++) {
          if(j!=c) {
            double overlapFactor = std::min(1.0, input.KunsatM(j,l)/(cmdTOmmolm2sMPa*fullRhizosphereOverlapConductivity));
            RHOPCohDyn(j,l) = RHOPcoh(j,l)*overlapFactor;
            RHOPCohDyn(c,l) = RHOPCohDyn(c,l) + (RHOPcoh(j,l) - RHOPCohDyn(j,l));
          }
        }
      }
      //Determine connected layers (non-zero fine root abundance)
      arma::Mat<uint8_t>& layerConnectedCoh = input.layerConnectedPools[c];
      input.nlayerscon[c] = 0;
      for(int j=0; j<numCohorts;j++) {
        for(int l=0;l<nlayers;l++) {
          if((x.belowLayers.V(c,l)>0.0) && (RHOPCohDyn(j,l)>0.0)) {
            layerConnectedCoh(j,l) = 1;
            input.nlayerscon[c]=input.nlayerscon[c] + 1;
          } else {
            layerConnectedCoh(j,l) = 0;
          }
        }
      }
      if(input.nlayerscon[c]==0) throw medfate::MedfateInternalError("Plant cohort not connected to any soil layer!");

      // Copy values from connected layers
      std::vector<double> Vc(input.nlayerscon[c]);
      std::vector<double> VCroot_kmaxc(input.nlayerscon[c]);
      std::vector<double> VGrhizo_kmaxc(input.nlayerscon[c]);
      std::vector<double> psic(input.nlayerscon[c]);
      std::vector<double> VG_nc(input.nlayerscon[c]);
      std::vector<double> VG_alphac(input.nlayerscon[c]);
      int cnt=0;
      for(int j=0; j<numCohorts;j++) {
        for(int l=0;l<nlayers;l++) {
          if(layerConnectedCoh(j,l)) {
            Vc[cnt] = x.belowLayers.V(c,l)*RHOPCohDyn(j,l);
            VCroot_kmaxc[cnt] = x.belowLayers.VCroot_kmax(c,l)*RHOPCohDyn(j,l);
            VGrhizo_kmaxc[cnt] = x.belowLayers.VGrhizo_kmax(c,l)*RHOPCohDyn(j,l);
            psic[cnt] = input.psiSoilM(j,l);
            VG_nc[cnt] = x.soil.getVG_n(l);
            VG_alphac[cnt] = x.soil.getVG_alpha(l);
            cnt++;
          }
        }
      }
      //Build supply function networks (Sperry transpiration mode)
      if(transpirationMode == "Sperry") {
        initSperryNetwork_inner_c(sperryNetworks[c], c,
                                  x.internalWater, x.paramsTranspiration, x.paramsWaterStorage,
                                  VCroot_kmaxc, VGrhizo_kmaxc,
                                  psic, VG_nc, VG_alphac,
                                  x.control,
                                  sapFluidityDay);
  //       supply[c] = supplyFunctionNetwork(HN, 0.0, 0.001); 
      } else if(transpirationMode == "Sureau") {
        initSureauNetwork_inner_c(sureauNetworks[c], c, LAIphe,
                                  x.internalWater,
                                  x.paramsAnatomy, x.paramsTranspiration, x.paramsWaterStorage,
                                  VCroot_kmaxc, VGrhizo_kmaxc,
                                  psic, VG_nc, VG_alphac,
                                  x.control, sapFluidityDay);
      }
    }
  }

  ////////////////////////////////
  // Create input and output objects to be filled in inner functions
  ////////////////////////////////
  // List innerOutput = List::create(
  //   _["Extraction"] = outputExtraction,
  //   _["ExtractionPools"] = transpOutput["ExtractionPools"],
  //  _["ExtractionInst"] = transpOutput["ExtractionInst"],
  //  _["RhizoPsi"] = minPsiRhizo,
  //   _["Plants"] = outputPlants,
  //  _["SunlitLeaves"] = transpOutput["SunlitLeaves"],
  //  _["ShadeLeaves"] = transpOutput["ShadeLeaves"],
  //  _["PlantsInst"] = transpOutput["PlantsInst"],
  //  _["SunlitLeavesInst"] = transpOutput["SunlitLeavesInst"],
  //  _["ShadeLeavesInst"] = transpOutput["ShadeLeavesInst"],
  //  _["LightExtinction"] = lightExtinctionAbsortion,
  //  _["LWRExtinction"] = lwrExtinctionList)
  // _["SupplyFunctions"] = supply,
  // _["PhotoSunlitFunctions"] = transpOutput["PhotoSunlitFunctions"],
  // _["PhotoShadeFunctions"] = transpOutput["PhotoShadeFunctions"],
  // _["PMSunlitFunctions"] = transpOutput["PMSunlitFunctions"],
  // _["PMShadeFunctions"] = transpOutput["PMShadeFunctions"]);
                                       
                                       
  ////////////////////////////////////////
  // STEP 5. Sub-daily (e.g. hourly) loop
  ////////////////////////////////////////
  for(int n=0;n<ntimesteps;n++) { //Time loop

    // Determine soil evaporation and snow melt for the corresponding step
    double soilEvapStep = outputEnergyBalance.SWRsoil[n]*(soilEvaporation/sum_abs_SWR_soil);
    double snowMeltStep = outputEnergyBalance.SWRsoil[n]*(snowMelt/sum_abs_SWR_soil);
    //Canopy evaporation (mm) in the current step and fraction of dry canopy
    double canEvapStep = canopyEvaporation*(outputEnergyBalance.SWRcan[n]/sum_abs_SWR_can);
    if(sum_abs_SWR_can==0.0) canEvapStep = 0.0;
    input.f_dry = 1.0;
    if(canEvapStep>0.0) {
      input.f_dry = 1.0 - std::min(1.0, canopyEvaporation/pet);
    }
    if(sum_abs_SWR_soil==0.0) { // avoid zero sums
      soilEvapStep = 0.0;
      snowMeltStep = 0.0;
    }
    if(sum_abs_SWR_can==0.0) { // avoid zero sums
      canEvapStep = 0.0;
      input.f_dry = 1.0;
    }

    //Retrieve fraction of sunlit and short-wave radiation absorbed for the current time step
    std::vector<double>& absPAR_SL_COH = abs_PAR_SL_COH_list[n];
    std::vector<double>& absPAR_SH_COH = abs_PAR_SH_COH_list[n];
    std::vector<double>& absSWR_SL_COH = abs_SWR_SL_COH_list[n];
    std::vector<double>& absSWR_SH_COH = abs_SWR_SH_COH_list[n];
    arma::mat& absSWR_SL_ML = abs_SWR_SL_ML_list[n];
    arma::mat& absSWR_SH_ML = abs_SWR_SH_ML_list[n];
    std::vector<double>& fsunlit = fsunlit_list[n];

    //Leaf area and Vmax/Jmax corresponding to sunlit and shade leaves
    for(int c=0;c<numCohorts;c++) {
      // Rcout<<"cohort "<<c<<":\n";
      //Constant properties through time steps
      std::vector<double> Vmax298layer(ncanlayers), Jmax298layer(ncanlayers);
      std::vector<double> SLarealayer(ncanlayers), SHarealayer(ncanlayers);
      double sn =0.0;
      double sumLAIme_c = std::accumulate(LAIme.col(c).begin(), LAIme.col(c).end(), 0.0);
      for(int i=(ncanlayers-1);i>=0.0;i--) {
        //Effect of nitrogen concentration decay through the canopy (Improvement: see 10.5194/bg-7-1833-2010)
        double fn = exp(-0.713*(sn+LAIme(i,c)/2.0)/sumLAIme_c);
        // Rcout<<" l"<<i<<" fsunlit: "<< fsunlit[i]<<" lai: "<< LAIme(i,c)<<" fn: "<< fn <<"\n";
        sn+=LAIme(i,c);
        SLarealayer[i] = LAIme(i,c)*fsunlit[i];
        SHarealayer[i] = LAIme(i,c)*(1.0-fsunlit[i]);
        Vmax298layer[i] = x.paramsTranspiration.Vmax298[c]*fn;
        Jmax298layer[i] = x.paramsTranspiration.Jmax298[c]*fn;
      }
      for(int i=0;i<ncanlayers;i++) {
        LAI_SL(c,n) +=SLarealayer[i];
        LAI_SH(c,n) +=SHarealayer[i];
        Vmax298_SL(c,n) +=Vmax298layer[i]*LAIme(i,c)*fsunlit[i];
        Jmax298_SL(c,n) +=Jmax298layer[i]*LAIme(i,c)*fsunlit[i];
        Vmax298_SH(c,n) +=Vmax298layer[i]*LAIme(i,c)*(1.0-fsunlit[i]);
        Jmax298_SH(c,n) +=Jmax298layer[i]*LAIme(i,c)*(1.0-fsunlit[i]);
      }
    }

    //Determine canopy vertical layer corresponding to cohort canopy, sunlit and shade leaves for each cohort
    for(int c=0;c<numCohorts;c++) {
      double num = 0.0, den = 0.0, numsl=0.0, densl =0.0, numsh = 0.0, densh=0.0;
      for(int i=0;i<ncanlayers;i++) {
        num += LAIme(i,c)*zmid[i];
        den += LAIme(i,c);
        numsl += LAIme(i,c)*zmid[i]*fsunlit[i];
        densl += LAIme(i,c)*fsunlit[i];
        numsh += LAIme(i,c)*zmid[i]*(1.0 - fsunlit[i]);
        densh += LAIme(i,c)*(1.0-fsunlit[i]);
      }
      double hc_sl = numsl/densl;
      double hc_sh = numsh/densh;
      double hc  = num/den;
      for(int i=0;i<ncanlayers;i++) {
        if((hc > zlow[i]) && (hc <=zup[i])) input.iLayerCohort[c] = i;
        if((hc_sl > zlow[i]) && (hc_sl <=zup[i])) input.iLayerSunlit[c] = i;
        if((hc_sh > zlow[i]) && (hc_sh <=zup[i])) input.iLayerShade[c] = i;
      }
    }

  //   List innerInput;
    if(transpirationMode =="Sperry") {
  //     innerInput = List::create(_["Patm"] = Patm,
  //                               _["zWind"] = zWind,
  //                               _["f_dry"] = f_dry,
  //                               _["iPMSunlit"] = iPMSunlit,
  //                               _["iPMShade"] = iPMShade,
  //                               _["nlayerscon"] = nlayerscon,
  //                               _["layerConnected"] = layerConnected,
  //                               _["layerConnectedPools"] = layerConnectedPools,
  //                               _["supply"] = supply);
    } 

    ////////////////////////////////////////
    // STEP 5.1 Long-wave radiation balance
    ////////////////////////////////////////
    longwaveRadiationSHAW_inner_c(ATres.lwrExtinction[n], 
                                  LAIme, LAImd, LAImx,
                                  lwdr[n], soil.getTemp(0), Tair, 0.1);
    outputEnergyBalance.LWRsoil[n] = ATres.lwrExtinction[n].Lnet_ground;
    outputEnergyBalance.LWRcan[n]= ATres.lwrExtinction[n].Lnet_canopy; 
    arma::mat& Lnet_cohort_layer = ATres.lwrExtinction[n].Lnet_cohort_layer;

    ////////////////////////////////////////
    // STEP 5.2 Sunlit/shade leaf energy balance, stomatal conductance and plant hydraulics
    ////////////////////////////////////////
    for(int c=0;c<numCohorts;c++) {
      //default values
      outputPlantsInst.dEdP(c,n) = medfate::NA_DOUBLE;
      outputPlantsInst.E(c,n) = 0.0;
      outputPlantsInst.Ag(c,n) = 0.0;
      outputPlantsInst.An(c,n) = 0.0;
      if(LAIphe[c]>0.0) {
        outputSunlitInst.Abs_PAR(c,n) = absPAR_SL_COH[c];
        outputShadeInst.Abs_PAR(c,n) = absPAR_SH_COH[c];
        outputSunlitInst.Abs_SWR(c,n) = absSWR_SL_COH[c];
        outputShadeInst.Abs_SWR(c,n) = absSWR_SH_COH[c];
        outputSunlitInst.Net_LWR(c,n) = 0.0;
        outputShadeInst.Net_LWR(c,n) = 0.0;
        for(int i=0;i<ncanlayers;i++) {
          outputSunlitInst.Net_LWR(c,n) += Lnet_cohort_layer(i,c)*fsunlit[i];
          outputShadeInst.Net_LWR(c,n) += Lnet_cohort_layer(i,c)*(1.0 - fsunlit[i]);
        }
      }
    }

    if(transpirationMode == "Sperry") {
  //     innerSperry(x, innerInput, innerOutput, n, tstep, 
  //                 verbose, stepFunctions);
    } else if(transpirationMode == "Sureau"){
      innerSureau_c(x, sureauNetworks, input, ATres , n, tstep);
    }

    for(int c=0;c<numCohorts;c++) {
      if(LAIlive[c]>0.0 && (x.internalWater.LeafPLC[c] < 0.999)) {
        //Store (for output) instantaneous leaf, stem and root potential, plc and rwc values
        outputPlantsInst.StemPLC(c,n) = x.internalWater.StemPLC[c];
        outputPlantsInst.LeafPLC(c,n) = x.internalWater.LeafPLC[c];
        outputPlantsInst.StemSympRWC(c,n) = symplasticRelativeWaterContent_c(x.internalWater.StemSympPsi[c], x.paramsWaterStorage.StemPI0[c], x.paramsWaterStorage.StemEPS[c]);
        outputPlantsInst.LeafSympRWC(c,n) = symplasticRelativeWaterContent_c(x.internalWater.LeafSympPsi[c], x.paramsWaterStorage.LeafPI0[c], x.paramsWaterStorage.LeafEPS[c]);
        outputPlantsInst.StemRWC(c,n) = outputPlantsInst.StemSympRWC(c,n)*(1.0 - x.paramsWaterStorage.StemAF[c]) + (1.0 - x.internalWater.StemPLC[c])*x.paramsWaterStorage.StemAF[c];
        outputPlantsInst.LeafRWC(c,n) = outputPlantsInst.LeafSympRWC(c,n)*(1.0 - x.paramsWaterStorage.LeafAF[c]) + (1.0 - x.internalWater.LeafPLC[c])*x.paramsWaterStorage.LeafAF[c];
        outputPlantsInst.StemPsi(c,n) = x.internalWater.StemPsi[c];
        outputPlantsInst.LeafPsi(c,n) = x.internalWater.LeafPsi[c]; //Store instantaneous (average) leaf potential
        outputPlantsInst.RootPsi(c,n) = x.internalWater.RootCrownPsi[c]; //Store instantaneous root crown potential
        outputPlantsInst.LeafSympPsi(c,n) = x.internalWater.LeafSympPsi[c];
        outputPlantsInst.StemSympPsi(c,n) = x.internalWater.StemSympPsi[c];

        if(n==0) {
          outputSunlit.GSWMin[c] = outputSunlitInst.Gsw(c,n);
          outputShade.GSWMin[c] = outputShadeInst.Gsw(c,n);
          outputSunlit.GSWMax[c] = outputSunlitInst.Gsw(c,n);
          outputShade.GSWMax[c] = outputShadeInst.Gsw(c,n);
          outputSunlit.TempMin[c] = outputSunlitInst.Temp(c,n);
          outputShade.TempMin[c] = outputShadeInst.Temp(c,n);
          outputSunlit.TempMax[c] = outputSunlitInst.Temp(c,n);
          outputShade.TempMax[c] = outputShadeInst.Temp(c,n);
          outputSunlit.LeafPsiMin[c] = outputSunlitInst.Psi(c,n);
          outputShade.LeafPsiMin[c] = outputShadeInst.Psi(c,n);
          outputSunlit.LeafPsiMax[c] = outputSunlitInst.Psi(c,n);
          outputShade.LeafPsiMax[c] = outputShadeInst.Psi(c,n);
          outputPlants.LeafPsiMin[c] = outputPlantsInst.LeafPsi(c,n);
          outputPlants.LeafPsiMax[c] = outputPlantsInst.LeafPsi(c,n);
          outputPlants.StemPsi[c] = outputPlantsInst.StemPsi(c,n);
          outputPlants.RootPsi[c] = outputPlantsInst.RootPsi(c,n);
          for(int l=0;l<nlayers;l++) {
            ATres.rhizoPsi(c,l) = x.belowLayers.RhizoPsi(c,l);
          }
        } else {
          outputSunlit.GSWMin[c] = std::min(outputSunlit.GSWMin[c], outputSunlitInst.Gsw(c,n));
          outputShade.GSWMin[c] = std::min(outputShade.GSWMin[c], outputShadeInst.Gsw(c,n));
          outputSunlit.GSWMax[c] = std::max(outputSunlit.GSWMax[c], outputSunlitInst.Gsw(c,n));
          outputShade.GSWMax[c] = std::max(outputShade.GSWMax[c], outputShadeInst.Gsw(c,n));
          outputSunlit.TempMin[c] = std::min(outputSunlit.TempMin[c], outputSunlitInst.Temp(c,n));
          outputShade.TempMin[c] = std::min(outputShade.TempMin[c], outputShadeInst.Temp(c,n));
          outputSunlit.TempMax[c] = std::max(outputSunlit.TempMax[c], outputSunlitInst.Temp(c,n));
          outputShade.TempMax[c] = std::max(outputShade.TempMax[c], outputShadeInst.Temp(c,n));
          outputSunlit.LeafPsiMin[c] = std::min(outputSunlit.LeafPsiMin[c], outputSunlitInst.Psi(c,n));
          outputShade.LeafPsiMin[c] = std::min(outputShade.LeafPsiMin[c], outputShadeInst.Psi(c,n));
          outputSunlit.LeafPsiMax[c] = std::max(outputSunlit.LeafPsiMax[c], outputSunlitInst.Psi(c,n));
          outputShade.LeafPsiMax[c] = std::max(outputShade.LeafPsiMax[c], outputShadeInst.Psi(c,n));
          outputPlants.LeafPsiMin[c] = std::min(outputPlants.LeafPsiMin[c], outputPlantsInst.LeafPsi(c,n));
          outputPlants.LeafPsiMax[c] = std::max(outputPlants.LeafPsiMax[c], outputPlantsInst.LeafPsi(c,n));
          outputPlants.StemPsi[c] = std::min(outputPlants.StemPsi[c], outputPlantsInst.StemPsi(c,n));
          outputPlants.RootPsi[c] = std::min(outputPlants.RootPsi[c], outputPlantsInst.RootPsi(c,n));
          
          for(int l=0;l<nlayers;l++) {
            ATres.rhizoPsi(c,l)  = std::min(ATres.rhizoPsi(c,l) , x.belowLayers.RhizoPsi(c,l));
          }
        }
      } else {
        // Assume constant PLC (so that it can be decreased in the future)
        outputPlantsInst.StemPLC(c,n) = x.internalWater.StemPLC[c];
        outputPlantsInst.LeafPLC(c,n) = x.internalWater.LeafPLC[c];
        outputPlantsInst.StemPsi(c,n) = x.internalWater.StemPsi[c];
        outputPlantsInst.LeafPsi(c,n) = x.internalWater.LeafPsi[c]; //Store instantaneous (average) leaf potential
        outputPlantsInst.RootPsi(c,n) = x.internalWater.RootCrownPsi[c]; //Store instantaneous root crown potential
        outputPlantsInst.LeafSympPsi(c,n) = x.internalWater.LeafSympPsi[c];
        outputPlantsInst.StemSympPsi(c,n) = x.internalWater.StemSympPsi[c];
        outputPlants.LeafPsiMin[c] = x.internalWater.LeafPsi[c];
        outputPlants.LeafPsiMax[c] = x.internalWater.LeafPsi[c];
        outputPlants.StemPsi[c] =  x.internalWater.StemPsi[c];
        outputPlants.RootPsi[c] = x.internalWater.RootCrownPsi[c];
        outputSunlit.LeafPsiMin[c] = x.internalWater.LeafPsi[c];
        outputShade.LeafPsiMin[c] = x.internalWater.LeafPsi[c];
        outputSunlit.LeafPsiMax[c] = x.internalWater.LeafPsi[c];
        outputShade.LeafPsiMax[c] = x.internalWater.LeafPsi[c];
        outputPlantsInst.StemSympRWC(c,n) = symplasticRelativeWaterContent_c(x.internalWater.StemSympPsi[c], x.paramsWaterStorage.StemPI0[c], x.paramsWaterStorage.StemEPS[c]);
        outputPlantsInst.LeafSympRWC(c,n) = symplasticRelativeWaterContent_c(x.internalWater.LeafSympPsi[c], x.paramsWaterStorage.LeafPI0[c], x.paramsWaterStorage.LeafEPS[c]);
        outputPlantsInst.StemRWC(c,n) = outputPlantsInst.StemSympRWC(c,n)*(1.0 - x.paramsWaterStorage.StemAF[c]) + (1.0 - x.internalWater.StemPLC[c])*x.paramsWaterStorage.StemAF[c];
        outputPlantsInst.LeafRWC(c,n) = outputPlantsInst.LeafSympRWC(c,n)*(1.0 - x.paramsWaterStorage.LeafAF[c]) + (1.0 - x.internalWater.LeafPLC[c])*x.paramsWaterStorage.LeafAF[c];
      }
    }

    ////////////////////////////////////////
    // STEP 5.3 Soil and canopy energy balances (single or multiple canopy layers)
    ////////////////////////////////////////

    //Soil latent heat (soil evaporation)
    //Latent heat (snow fusion) as J/m2/s
    if(x.snowpack>0.0) {
      outputEnergyBalance.SWRsoil[n] = 0.0; //Set SWR absorbed by soil to zero (for energy balance) if snow pack is present
      outputEnergyBalance.LWRsoil[n] = 0.0; //Set net LWR to zero
    }
    outputEnergyBalance.LEVsoil[n] = (1e6)*latentHeatVaporisation_c(x.soil.getTemp(0))*soilEvapStep/tstep;
    // Rcout<<n<<" "<<sum_abs_SWR_soil<<" "<<soilEvapStep << " "<<Tsoil[0]<<" " << LEVsoil[n]<<"\n";
    outputEnergyBalance.LEFsnow[n] = (1e6)*(snowMeltStep*0.33355)/tstep; // 0.33355 = latent heat of fusion

    //Herbaceous transpiration (mm) in the current step
    double herbTranspStep = herbTranspiration*(outputEnergyBalance.SWRcan[n]/sum_abs_SWR_can);

    //Canopy convective heat exchange
    double RAcan = aerodynamicResistance_c(canopyHeight,std::max(wind,1.0)); //Aerodynamic resistance to convective heat transfer
    outputEnergyBalance.Hcan[n] = (airDensity_c(Tatm[n],Patm)*Cp_JKG*(Tcan[n]-Tatm[n]))/RAcan;

    if(!multiLayerBalance) {//Canopy balance assuming a single layer
      //Soil-canopy turbulent heat exchange
      double wind2m = windSpeedMassmanExtinction_c(200.0, wind, LAIcell, canopyHeight);
      double RAsoil = aerodynamicResistance_c(200.0, std::max(wind2m,1.0)); //Aerodynamic resistance to convective heat transfer from soil
      if(x.snowpack==0.0) {
        outputEnergyBalance.Hcansoil[n] = (airDensity_c(Tcan[n],Patm)*Cp_JKG*(Tcan[n]-x.soil.getTemp(0)))/RAsoil;
      } else {
        outputEnergyBalance.Hcansoil[n] = (airDensity_c(Tcan[n],Patm)*Cp_JKG*(Tcan[n] - 0.0))/RAsoil; //Assumes a zero degree for soil surface (snow)
      }
      //Latent heat (evaporation + transpiration)
      double sum_Einst_n = 0.0;
      for(int c=0;c<numCohorts;c++) sum_Einst_n += outputPlantsInst.E(c, n);
      double LEwat = (1e6)*latentHeatVaporisation_c(Tcan[n])*(sum_Einst_n + canEvapStep + herbTranspStep)/tstep;
      outputEnergyBalance.LEVcan[n] = LEwat;
      // Rcout<< n <<" " << sum_Einst_n << " " << canEvapStep <<" "<< outputEnergyBalance.LEVcan[n]<<"\n";
      
      //Canopy temperature changes
      outputEnergyBalance.Ebalcan[n] = outputEnergyBalance.SWRcan[n] + outputEnergyBalance.LWRcan[n] - outputEnergyBalance.LEVcan[n] - outputEnergyBalance.LEFsnow[n] - outputEnergyBalance.Hcan[n] - outputEnergyBalance.Hcansoil[n];
      double canopyAirThermalCapacity = airDensity_c(Tcan[n],Patm)*Cp_JKG;
      double canopyThermalCapacity =  canopyAirThermalCapacity + (0.5*(0.8*LAIcelllive + 1.2*LAIcell) + LAIcelldead)*thermalCapacityLAI; //Avoids zero capacity for winter deciduous
      double Tcannext = Tcan[n]+ std::max(-3.0, std::min(3.0, tstep*outputEnergyBalance.Ebalcan[n]/canopyThermalCapacity)); //Avoids changes in temperature that are too fast
      if(n<(ntimesteps-1)) Tcan[n+1] = Tcannext;
      for(int i=0;i<ncanlayers;i++) Tair[i] = Tcannext;

      //Soil energy balance including exchange with canopy
      if(x.snowpack==0.0) {
        outputEnergyBalance.Ebalsoil[n] = outputEnergyBalance.SWRsoil[n] + outputEnergyBalance.LWRsoil[n] + outputEnergyBalance.Hcansoil[n] - outputEnergyBalance.LEVsoil[n]; //Here we use all energy escaping to atmosphere
      } else {
        //Heat conduction between soil and snow
        double Hcond_snow = 0.05*(0.0 - x.soil.getTemp(0));//0.05 Wm-1K-1 for fresh snow
        outputEnergyBalance.Ebalsoil[n] = Hcond_snow - outputEnergyBalance.LWRsoil[n] - outputEnergyBalance.LEVsoil[n];
      }

      //Soil temperature changes
      std::vector<double> Ws(nlayers), Tsoil(nlayers);
      for(int l=0;l<nlayers;l++) {
        Ws[l] = x.soil.getW(l);
        Tsoil[l] = x.soil.getTemp(l);
      }
      temperatureChange_inner_c(ATcomm.SEBcomm, 
                                x.soil.getWidths(), 
                                Tsoil, x.soil.getSand(), x.soil.getClay(), Ws, 
                                x.soil.getThetaSAT(), x.soil.getThetaFC(), 
                                outputEnergyBalance.Ebalsoil[n], tstep);
      for(int l=0;l<nlayers;l++) {
        x.soil.setTemp(l, x.soil.getTemp(l) + std::max(-3.0, std::min(3.0, ATcomm.SEBcomm.tempch[l])));
        if(n<(ntimesteps-1)) outputEnergyBalance.SoilTemperature(n+1,l)= x.soil.getTemp(l);
      }
    } else { //Multilayer canopy balance
      double moistureAtm = 0.622*(vpatm/Patm)*airDensity_c(Tatm[n],Patm);
      double CO2Atm = 0.409*Catm*44.01; //mg/m3

      double tsubstep = tstep/((double) nsubsteps);
      double maxTchange = 3.0/((double) nsubsteps);
      double maxMoistureChange = 0.001/((double)nsubsteps); //=0.16 kPa per step
      double maxCO2Change = 180.0/((double)nsubsteps); //= 10 ppm per step
      double deltaZ = (verticalLayerSize/100.0); //Vertical layer size in m
      std::vector<double>& LWRnet_layer = ATres.lwrExtinction[n].LWR_layer.Lnet;
      outputEnergyBalance.Ebalcan[n] = 0.0;
      outputEnergyBalance.LEVcan[n] = 0.0;
      std::vector<double> Tairnext(ncanlayers), LElayer(ncanlayers), absSWRlayer(ncanlayers), Rnlayer(ncanlayers), Hleaflayer(ncanlayers);
      std::vector<double> layerThermalCapacity(ncanlayers);
      std::vector<double> moistureET(ncanlayers), rho(ncanlayers), moistureLayer(ncanlayers), moistureLayernext(ncanlayers);
      std::vector<double> CO2An(ncanlayers), CO2Layer(ncanlayers), CO2Layernext(ncanlayers);
      for(int i=0;i<ncanlayers;i++) {
        rho[i] = airDensity_c(Tair[i],Patm);
        absSWRlayer[i] = 0.0;
        for(int c; c< numCohorts; c++) absSWRlayer[i] += absSWR_SL_ML(i,c) + absSWR_SH_ML(i,c);
        //Radiation balance
        Rnlayer[i] = absSWRlayer[i] + LWRnet_layer[i];
        // std::vector<double> pLayer = LAIme(i,_)/LAIphe; //Proportion of each cohort LAI in layer i
        //Instantaneous layer transpiration
        //from mmolH2O/m2/s to kgH2O/m2/s
        double ElayerInst = 0.0;
        for(int c; c< numCohorts; c++) ElayerInst += 0.001*0.01802*LAIme(i,c)*(outputSunlitInst.E(c,n)*fsunlit[i] + outputShadeInst.E(c,n)*(1.0-fsunlit[i]));
        //Assumes Layers contribute to evaporation proportionally to their LAI fraction
        double layerEvapInst = (canEvapStep/tstep)*(LAIpe[i]/LAIcellexpanded);
        //Instantaneous herbaceous transpiration (for bottom layer)
        double herbTranspInst = 0.0;
        if(i==0) herbTranspInst = (herbTranspStep/tstep);
        //Estimate instantaneous mgCO2/m2 absorption for the layer, taking into account the proportion of sunlit and shade leaves of each cohort
        //from micro.molCO2/m2/s to mgCO2/m2/s
        double Anlayer = 0.0;
        for(int c; c< numCohorts; c++) Anlayer += (1e-3)*44.01*LAIme(i,c)*(outputSunlitInst.An(c,n)*fsunlit[i] + outputShadeInst.An(c,n)*(1.0-fsunlit[i]));
        // 1000.0*(44.01/12.0)*sum(Aninst(_,n)*pLayer);
        double LEwat = (1e6)*latentHeatVaporisation_c(Tair[i])*(ElayerInst + layerEvapInst+ herbTranspInst);
        LElayer[i] = LEwat; //Energy spent in vaporisation
        if(i==0) LElayer[i] = LElayer[i] - outputEnergyBalance.LEFsnow[n]; //Add latent heat of fusion to first layer
        outputEnergyBalance.LEVcan[n] = LElayer[i];
        // layerThermalCapacity[i] = (0.5*(0.8*LAIcelllive + 1.2*LAIcell) + LAIcelldead)*thermalCapacityLAI/((double) ncanlayers);
        layerThermalCapacity[i] =  (0.5*(0.8*LAIpx[i] + 1.2*LAIpe[i]) + LAIpd[i])*thermalCapacityLAI; //Avoids zero capacity for winter deciduous

        moistureLayer[i] = 0.622*(VPair[i]/Patm)*rho[i]; //kg water vapour/m3
        moistureET[i] = (ElayerInst + layerEvapInst + herbTranspInst)/(deltaZ); //kg water vapour /m3/s

        CO2Layer[i] = 0.409*Cair[i]*44.01; //mg/m3
        CO2An[i] = -1.0*Anlayer/(deltaZ); //mg/m3/s
        // Rcout<<n<< " "<< i<< " - Rn: "<<Rnlayer[i]<<" LE: "<<LElayer[i]<<" Hleaf: "<<Hleaflayer[i]<< " Tini: "<< Tair[i]<<"\n";
        // Rcout<<n<< " "<< i<< " - moistureET: "<<moistureET[i]<<" moistureLayer: "<<moistureLayer[i]<<" CO2An: "<<CO2An[i]<< " CO2Layer: "<< CO2Layer[i]<<"\n";
      }
      //Add soil moisture evaporation
      moistureET[0] += soilEvapStep/(deltaZ*tstep); //kg/m3/s

      for(int s=0;s<nsubsteps;s++) {
        double RAsoil = aerodynamicResistance_c(200.0, std::max(input.zWind[0],1.0)); //Aerodynamic resistance to convective heat transfer from soil
        double Hcansoils = Cp_JKG*airDensity_c(Tair[0],Patm)*(Tair[0] - x.soil.getTemp(0))/RAsoil;
        // double Hcan_heats = (airDensity_c(Tatm[n],Patm)*Cp_JKG*(Tair[ncanlayers-1]-Tatm[n]))/RAcan;
        for(int i=0;i<ncanlayers;i++) {
          double deltaH = 0.0;
          double deltaMoisture = 0.0;
          double deltaCO2 = 0.0;
          // double Hlayers = 0.0;
          //Add turbulent heat flow (positive gradient when temperature is larger above)
          if(i==0) { //Lower layer
            deltaH -= (Cp_JKG*rho[i]*(Tair[i+1] - Tair[i])*uw[i])/(deltaZ*dU[i]);
            deltaH -= Hcansoils;
            deltaMoisture -= ((moistureLayer[i+1] - moistureLayer[i])*uw[i])/(deltaZ*dU[i]);
            deltaCO2 -= ((CO2Layer[i+1] - CO2Layer[i])*uw[i])/(deltaZ*dU[i]);
          } else if((i > 0) && (i<(ncanlayers-1))) { //Intermediate layers
            deltaH -= (Cp_JKG*rho[i]*(Tair[i+1] - Tair[i-1])*uw[i])/(2.0*deltaZ*dU[i]);
            deltaMoisture -= ((moistureLayer[i+1] - moistureLayer[i-1])*uw[i])/(2.0*deltaZ*dU[i]);
            deltaCO2 -= ((CO2Layer[i+1] - CO2Layer[i-1])*uw[i])/(2.0*deltaZ*dU[i]);
          } else if(i==(ncanlayers-1)){ //Upper layer
            deltaH -= (Cp_JKG*rho[i]*(Tatm[n] - Tair[i-1])*uw[i])/(2.0*deltaZ*dU[i]);
            deltaMoisture -= ((moistureAtm - moistureLayer[i-1])*uw[i])/(2.0*deltaZ*dU[i]);
            deltaCO2 -= ((CO2Atm - CO2Layer[i-1])*uw[i])/(2.0*deltaZ*dU[i]);
          }
          Hleaflayer[i] = 0.0;
          for(int c=0;c<numCohorts;c++) {
            double gHa = 0.189*pow(std::max(input.zWind[i],0.1)/(x.paramsAnatomy.LeafWidth[c]*0.0072), 0.5);
            double Hsunlit = 2.0*Cp_Jmol*rho[i]*(outputSunlitInst.Temp(c, n)-Tair[i])*gHa;
            double Hshade = 2.0*Cp_Jmol*rho[i]*(outputShadeInst.Temp(c, n)-Tair[i])*gHa;
            // Rcout<<c<<" " << Hsunlit<< " "<<Hshade<<" \n";
            Hleaflayer[i] +=(Hsunlit*fsunlit[i] + Hshade*(1.0-fsunlit[i]))*LAIme(i,c);
          }
          double EbalLayer = Rnlayer[i] - LElayer[i] + Hleaflayer[i] + deltaH;
          //Instantaneous changes in temperature due to internal energy balance
          double deltaT = EbalLayer/(rho[i]*Cp_JKG + layerThermalCapacity[i]);

          // if(s==0) Rcout<<n<< " "<< i<< " "<< s <<" - Rn: "<<Rnlayer[i]<<" LE: "<<LElayer[i]<<" Hleaf: "<<Hleaflayer[i]<<" H: "<<Hlayers<< " Ebal: "<<EbalLayer<< " LTC: " << rholayer*Cp_JKG + layerThermalCapacity[i]<< " Tini: "<< Tair[i]<< " deltaT: "<<deltaT<<"\n";
          Tairnext[i] = Tair[i] +  std::max(-1.0*maxTchange, std::min(maxTchange, tsubstep*deltaT)); //Avoids changes in temperature that are too fast
          //Changes in water vapour
          moistureLayernext[i] = moistureLayer[i] + std::max(-1.0*maxMoistureChange, std::min(maxMoistureChange, tsubstep*(moistureET[i]+ deltaMoisture)));
          // if(i==0) Rcout<<n<< " "<< i<< " "<< s <<" - moisture: "<<moistureLayer[i]<<" delta: "<<deltaMoisture<<" next: "<< moistureLayernext[i]<<"\n";
          //Changes in CO2
          CO2Layernext[i] = CO2Layer[i] + std::max(-1.0*maxCO2Change, std::min(maxCO2Change, tsubstep*(CO2An[i] + deltaCO2)));
        }
        for(int i=0;i<ncanlayers;i++) {
          Tair[i] = Tairnext[i];
          moistureLayer[i] = moistureLayernext[i];
          CO2Layer[i] = CO2Layernext[i];
        }

        //Soil energy balance including exchange with canopy
        double Ebalsoils = 0.0;
        //Soil energy balance including exchange with canopy
        if(x.snowpack==0.0) {
          Ebalsoils = outputEnergyBalance.SWRsoil[n] + outputEnergyBalance.LWRsoil[n]  - outputEnergyBalance.LEVsoil[n] + Hcansoils; //Here we use all energy escaping to atmosphere
        } else {
          //Heat conduction between soil and snow
          double Hcond_snow = 0.05*(0.0 - x.soil.getTemp(0));//0.05 Wm-1K-1 for fresh snow
          Ebalsoils = Hcond_snow - outputEnergyBalance.LWRsoil[n] - outputEnergyBalance.LEVsoil[n];
        }
        outputEnergyBalance.Ebalsoil[n] +=Ebalsoils;
        outputEnergyBalance.Hcansoil[n] +=Hcansoils;
        //Soil temperature changes
        std::vector<double> Ws(nlayers), Tsoil(nlayers);
        for(int l=0;l<nlayers;l++) {
          Ws[l] = x.soil.getW(l);
          Tsoil[l] = x.soil.getTemp(l);
        }
        temperatureChange_inner_c(ATcomm.SEBcomm, 
                                  x.soil.getWidths(), 
                                  Tsoil, x.soil.getSand(), x.soil.getClay(), Ws, 
                                  x.soil.getThetaSAT(), x.soil.getThetaFC(), 
                                  outputEnergyBalance.Ebalsoil[n], tsubstep);
        for(int l=0;l<nlayers;l++) {
          x.soil.setTemp(l, x.soil.getTemp(l) + std::max(-3.0, std::min(3.0, ATcomm.SEBcomm.tempch[l])));
        }
      }
      outputEnergyBalance.Hcansoil[n] = outputEnergyBalance.Hcansoil[n]/((double) nsubsteps);
      outputEnergyBalance.Ebalsoil[n] = outputEnergyBalance.Ebalsoil[n]/((double) nsubsteps);
      for(int i=0;i<ncanlayers;i++) {
        VPair[i] = ((moistureLayer[i]/rho[i])*Patm)/0.622;
        Cair[i] = CO2Layer[i]/(0.409*44.01);
        // Rcout<< n << " "<<i << " - " << moistureLayer[i]<< " "<< VPair[i]<<"\n";
      }
      // Rcout<< n << " "<<abs_SWR_can[n]<< " = "<< sum(absSWRlayer)<<" = "<< (sum(absSWR_SL_ML) + sum(absSWR_SH_ML))<<" "<<net_LWR_can[n]<< " = "<< sum(LWRnet_layer)<<"\n";
      //Canopy energy balance
      outputEnergyBalance.Ebalcan[n] = outputEnergyBalance.SWRcan[n]+ outputEnergyBalance.LWRcan[n] - outputEnergyBalance.LEVcan[n] - outputEnergyBalance.Hcan[n] - outputEnergyBalance.Hcansoil[n];
      if(n<(ntimesteps-1)) {
        double numSum = 0.0;
        double denSum = 0.0;
        for(int i=0;i<ncanlayers;i++) {
          numSum +=Tair[i]*LAIpx[i];
          denSum +=LAIpx[i];
        }
        Tcan[n+1] = numSum/denSum;
        for(int l=0;l<nlayers;l++) outputEnergyBalance.SoilTemperature(n+1,l)= x.soil.getTemp(l);
      }
    }
    if(n<(ntimesteps-1)) for(int i=0;i<ncanlayers;i++) {
      outputEnergyBalance.TemperatureLayers(n+1,i) = Tair[i];
      outputEnergyBalance.VaporPressureLayers(n+1,i) = VPair[i];
    }
  } //End of timestep loop

  //Delete Sureau Networks
  if(transpirationMode == "Sureau") {
    for(int c=0;c<numCohorts;c++) {
      deleteSureauNetworkPointers_c(sureauNetworks[c]);
    }
  }
  delete[] sureauNetworks;
  delete[] sperryNetworks;
  
  ////////////////////////////////////////
  // STEP 6. Plant drought stress (relative whole-plant conductance), cavitation and live fuel moisture
  ////////////////////////////////////////
  for(int c=0;c<numCohorts;c++) {
    outputPlants.Extraction[c] = 0.0;
    for(int l=0;l<nlayers;l++) {
      outputPlants.Extraction[c] += outputExtraction(c,l);
    }
    outputPlants.StemPLC[c] = 0.0;
    outputPlants.LeafPLC[c] = 0.0;
    outputPlants.StemRWC[c] = 0.0;
    outputPlants.LeafRWC[c] = 0.0;
    outputPlants.dEdP[c] = 0.0;
    for(int n=0;n<ntimesteps;n++) {
      outputPlants.StemPLC[c] += outputPlantsInst.StemPLC(c,n);
      outputPlants.LeafPLC[c] += outputPlantsInst.LeafPLC(c,n);
      outputPlants.StemRWC[c] += outputPlantsInst.StemRWC(c,n);
      outputPlants.LeafRWC[c] += outputPlantsInst.LeafRWC(c,n);
      outputPlants.dEdP[c] += outputPlantsInst.dEdP(c,n);
    }
    outputPlants.StemPLC[c] = outputPlants.StemPLC[c]/((double) ntimesteps);
    outputPlants.LeafPLC[c] = outputPlants.LeafPLC[c]/((double) ntimesteps);
    outputPlants.StemRWC[c] = outputPlants.StemRWC[c]/((double) ntimesteps);
    outputPlants.LeafRWC[c] = outputPlants.LeafRWC[c]/((double) ntimesteps);
    outputPlants.dEdP[c] = outputPlants.dEdP[c]/((double) ntimesteps);
    // The fraction of leaves will decrease due to phenology or processes leading to defoliation
    double fleaf = (1.0/x.paramsAnatomy.r635[c])*(LAIphe[c]/LAIlive[c]);
    if(lfmcComponent=="fine") { //fine fuel moisture
      outputPlants.LFMC[c] = x.paramsWaterStorage.maxMCleaf[c]*outputPlants.LeafRWC[c]*fleaf + x.paramsWaterStorage.maxMCstem[c]*outputPlants.StemRWC[c]*(1.0 - fleaf);
    } else { //"leaf"
      outputPlants.LFMC[c] = x.paramsWaterStorage.maxFMC[c]*outputPlants.LeafRWC[c];
    }

    outputPlants.DDS[c] = (1.0 - (outputPlants.dEdP[c]/(sapFluidityDay*x.paramsTranspiration.Plant_kmax[c])));
    if(x.paramsPhenology.phenoType[c] == "winter-deciduous" || x.paramsPhenology.phenoType[c] == "winter-semideciduous") {
      outputPlants.DDS[c] = x.internalPhenology.phi[c]*outputPlants.DDS[c];
      if(x.internalPhenology.phi[c] == 0.0) outputPlants.DDS[c] = 0.0;
    }

    double SAmax = 10e4/x.paramsAnatomy.Al2As[c]; //cm2·m-2 of leaf area
    double r = cavitationRecoveryMaximumRate*std::max(0.0, (x.internalWater.StemSympPsi[c] + 1.5)/1.5);
    if(stemCavitationRecovery=="rate") {
      x.internalWater.StemPLC[c] = std::max(0.0, x.internalWater.StemPLC[c] - (r/SAmax));
    }
    if(leafCavitationRecovery=="rate") {
      x.internalWater.LeafPLC[c] = std::max(0.0, x.internalWater.LeafPLC[c] - (r/SAmax));
    }
  }
  
  // Copy output stand
  outputStand.LAI = LAIcell;
  outputStand.LAIlive = LAIcelllive;
  outputStand.LAIexpanded = LAIcellexpanded;
  outputStand.LAIdead = LAIcelldead;
  
  // Copy output plants
  for(int c =0;c<numCohorts;c++) {
    outputPlants.LAI[c] = LAIphe[c];
    outputPlants.LAIlive[c] = LAIlive[c];
  }
}
