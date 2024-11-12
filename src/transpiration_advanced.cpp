#define STRICT_R_HEADERS
#include <Rcpp.h>
#include <numeric>
#include "communication_structures.h"
#include "lightextinction_basic.h"
#include "lightextinction_advanced.h"
#include "windextinction.h"
#include "windKatul.h"
#include "hydraulics.h"
#include "hydrology.h"
#include "biophysicsutils.h"
#include "phenology.h"
#include "forestutils.h"
#include "tissuemoisture.h"
#include "carbon.h"
#include "spwb.h"
#include "root.h"
#include "soil.h"
#include "soil_thermodynamics.h"
#include "inner_sperry.h"
#include "inner_sureau.h"
#include <meteoland.h>
using namespace Rcpp;

const double SIGMA_Wm2 = 5.67*1e-8;
const double Cp_JKG = 1013.86; // J * kg^-1 * ºC^-1
const double Cp_Jmol = 29.37152; // J * mol^-1 * ºC^-1

// SCHEDULE - Following steps for one day, given a weather vector:
//
// STEP 1. Estimate stand-level leaf area values and leaf distribution across layers from leaf-level live/expanded area
// STEP 2. Determine vertical wind speed profile
// STEP 3a. Direct and diffuse short-wave radiation for sub-steps
// STEP 3b. Above-canopy air temperature  and long-wave radiation emission for sub-steps
// STEP 3c. Short-wave radiation extinction and absorption for sub-steps
// STEP 4. Hydraulics: determine layers where the plant is connected and supply functions (Sperry mode)
// STEP 5. Sub-daily (e.g. hourly) loop
// STEP 5.1 Long-wave radiation balance
// STEP 5.2 Leaf energy balance, stomatal conductance and plant hydraulics  (Sperry or Sureau inner functions)
// STEP 5.3 Soil and canopy energy balances (single or multiple canopy layers)
// STEP 6. Update plant drought stress (relative whole-plant conductance), cavitation and live fuel moisture
void transpirationAdvanced(List transpOutput, List x, NumericVector meteovec, 
                  double latitude, double elevation, double slope, double aspect, 
                  double solarConstant, double delta,
                  double canopyEvaporation = 0.0, double snowMelt = 0.0, double soilEvaporation = 0.0, double herbTranspiration = 0.0,
                  bool verbose = false, int stepFunctions = NA_INTEGER, 
                  bool modifyInput = true) {
  
  //Will not modify input x 
  if(!modifyInput) {
    x = clone(x);
  }
  // Should have internal communication 
  NumericVector outputStand = as<NumericVector>(transpOutput["Stand"]);
  DataFrame outputPlants = as<DataFrame>(transpOutput["Plants"]);
  NumericMatrix outputExtraction = as<NumericMatrix>(transpOutput["Extraction"]);
  DataFrame outputSunlit = as<DataFrame>(transpOutput["SunlitLeaves"]);
  DataFrame outputShade = as<DataFrame>(transpOutput["ShadeLeaves"]);
  List outputPlantsInst = transpOutput["PlantsInst"];
  List outputSunlitInst = transpOutput["SunlitLeavesInst"];
  List outputShadeInst = transpOutput["ShadeLeavesInst"];
  List outputEnergyBalance = transpOutput["EnergyBalance"];
  DataFrame outputTemperatureInst =   as<DataFrame>(outputEnergyBalance["Temperature"]);
  DataFrame outputCEBinst =  as<DataFrame>(outputEnergyBalance["CanopyEnergyBalance"]);
  DataFrame outputSEBinst =  as<DataFrame>(outputEnergyBalance["SoilEnergyBalance"]);
  List lwrExtinctionList = transpOutput["LWRExtinction"];
  List supply = transpOutput["SupplyFunctions"]; 

  
  //Control parameters
  List control = x["control"];
  String transpirationMode = control["transpirationMode"];
  String soilFunctions = control["soilFunctions"];

  int ntimesteps = control["ndailysteps"];
  int nsubsteps = control["nsubsteps_canopy"];
  String rhizosphereOverlap = control["rhizosphereOverlap"];
  bool plantWaterPools = (rhizosphereOverlap!="total");
  double verticalLayerSize = control["verticalLayerSize"];
  double windMeasurementHeight  = control["windMeasurementHeight"];
  double thermalCapacityLAI = control["thermalCapacityLAI"];
  bool multiLayerBalance = control["multiLayerBalance"];
  double defaultWindSpeed = control["defaultWindSpeed"];
  String stemCavitationRecovery = control["stemCavitationRecovery"];
  String leafCavitationRecovery = control["leafCavitationRecovery"];
  double cavitationRecoveryMaximumRate = control["cavitationRecoveryMaximumRate"];
  bool sapFluidityVariation = control["sapFluidityVariation"];
  String lfmcComponent = control["lfmcComponent"];

  //Meteo input
  double tmin = meteovec["tmin"];
  double tmax = meteovec["tmax"];
  double tminPrev = meteovec["tminPrev"];
  double tmaxPrev = meteovec["tmaxPrev"];
  double tminNext = meteovec["tminNext"];
  double prec = meteovec["prec"];
  double rhmin = meteovec["rhmin"];
  double rhmax = meteovec["rhmax"];
  double rad = meteovec["rad"];
  double wind = meteovec["wind"];
  double Catm = meteovec["Catm"];
  double Patm = meteovec["Patm"];
  double pet = meteovec["pet"];
  
  //Atmospheric pressure (if missing)
  if(NumericVector::is_na(Patm)) Patm = meteoland::utils_atmosphericPressure(elevation);
  
  //Vegetation input
  DataFrame cohorts = Rcpp::as<Rcpp::DataFrame>(x["cohorts"]);
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  NumericVector LAIlive = Rcpp::as<Rcpp::NumericVector>(above["LAI_live"]);
  NumericVector LAI = Rcpp::as<Rcpp::NumericVector>(above["LAI_expanded"]);
  NumericVector LAIdead = Rcpp::as<Rcpp::NumericVector>(above["LAI_dead"]);
  NumericVector H = Rcpp::as<Rcpp::NumericVector>(above["H"]);
  NumericVector CR = Rcpp::as<Rcpp::NumericVector>(above["CR"]);

  int numCohorts = above.nrow();
  
  //Soil input
  DataFrame soil = Rcpp::as<Rcpp::DataFrame>(x["soil"]);
  int nlayers = soil.nrow();
  NumericVector widths = soil["widths"];
  NumericVector Water_FC = waterFC(soil, soilFunctions);
  NumericVector Theta_FC = thetaFC(soil, soilFunctions);
  NumericVector Theta_SAT = thetaSAT(soil, soilFunctions);
  NumericVector VG_n = Rcpp::as<Rcpp::NumericVector>(soil["VG_n"]);
  NumericVector VG_alpha = Rcpp::as<Rcpp::NumericVector>(soil["VG_alpha"]);
  NumericVector Tsoil = soil["Temp"]; 
  NumericVector sand = soil["sand"];
  NumericVector clay = soil["clay"];
  NumericVector Ws = soil["W"]; //Access to soil state variable
  double snowpack = x["snowpack"];

  //Canopy params
  DataFrame canopyParams = Rcpp::as<Rcpp::DataFrame>(x["canopy"]);
  int ncanlayers = canopyParams.nrow();
  NumericVector zlow = canopyParams["zlow"];
  NumericVector zmid = canopyParams["zmid"];
  NumericVector zup = canopyParams["zup"];
  NumericVector LAIpx = canopyParams["LAIlive"];
  NumericVector LAIpe = canopyParams["LAIexpanded"];
  NumericVector LAIpd = canopyParams["LAIdead"];
  NumericVector Tair = canopyParams["Tair"];
  NumericVector VPair = canopyParams["VPair"];
  NumericVector Cair = canopyParams["Cair"];
  for(int l=0;l<ncanlayers;l++) { //If canopy layers have missing values, then initialize with Catm
    if(!multiLayerBalance) Cair[l] = Catm;
    else {
      if(NumericVector::is_na(Cair[l])) Cair[l] = Catm; 
    }
  }
  if(multiLayerBalance) Cair[ncanlayers-1] = Catm;
  
  //Root distribution input
  DataFrame belowdf = Rcpp::as<Rcpp::DataFrame>(x["below"]);
  List belowLayers = Rcpp::as<Rcpp::List>(x["belowLayers"]);
  NumericMatrix V =Rcpp::as<Rcpp::NumericMatrix>(belowLayers["V"]);
  NumericMatrix VCroot_kmax= Rcpp::as<Rcpp::NumericMatrix>(belowLayers["VCroot_kmax"]);
  NumericMatrix VGrhizo_kmax= Rcpp::as<Rcpp::NumericMatrix>(belowLayers["VGrhizo_kmax"]);
  NumericMatrix RhizoPsiMAT = Rcpp::as<Rcpp::NumericMatrix>(belowLayers["RhizoPsi"]);
  
  //Water pools
  NumericMatrix Wpool = Rcpp::as<Rcpp::NumericMatrix>(belowLayers["Wpool"]);
  List RHOP;
  NumericVector poolProportions(numCohorts);
  if(plantWaterPools) {
    RHOP = belowLayers["RHOP"];
    poolProportions = belowdf["poolProportions"];
  }
  
  //Base parameters
  DataFrame paramsInterception = Rcpp::as<Rcpp::DataFrame>(x["paramsInterception"]);
  NumericVector alphaSWR = Rcpp::as<Rcpp::NumericVector>(paramsInterception["alphaSWR"]);
  NumericVector gammaSWR = Rcpp::as<Rcpp::NumericVector>(paramsInterception["gammaSWR"]);
  NumericVector kPAR = Rcpp::as<Rcpp::NumericVector>(paramsInterception["kPAR"]); //Not used in light extinction
  NumericVector Beta_p = Rcpp::as<Rcpp::NumericVector>(paramsInterception["Beta_p"]);
  NumericVector Beta_q = Rcpp::as<Rcpp::NumericVector>(paramsInterception["Beta_q"]);
  NumericVector ClumpingIndex = Rcpp::as<Rcpp::NumericVector>(paramsInterception["ClumpingIndex"]);
  
  //Phenology parameters
  DataFrame paramsPhenology = Rcpp::as<Rcpp::DataFrame>(x["paramsPhenology"]);
  CharacterVector phenoType = paramsPhenology["PhenologyType"];
  
  //Anatomy parameters
  DataFrame paramsAnatomy = Rcpp::as<Rcpp::DataFrame>(x["paramsAnatomy"]);
  NumericVector leafWidth = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["LeafWidth"]);
  NumericVector Al2As = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["Al2As"]);
  NumericVector r635 = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["r635"]);
  
  //Transpiration parameters
  DataFrame paramsTranspiration = Rcpp::as<Rcpp::DataFrame>(x["paramsTranspiration"]);
  NumericVector Plant_kmax = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["Plant_kmax"]);
  NumericVector VCstem_kmax = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["VCstem_kmax"]);
  NumericVector VCleaf_kmax = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["VCleaf_kmax"]);

  NumericVector Vmax298 = paramsTranspiration["Vmax298"];
  NumericVector Jmax298 = paramsTranspiration["Jmax298"];

  //Water storage parameters
  DataFrame paramsWaterStorage = Rcpp::as<Rcpp::DataFrame>(x["paramsWaterStorage"]);
  NumericVector maxFMC = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["maxFMC"]);
  NumericVector StemPI0 = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["StemPI0"]);
  NumericVector StemEPS = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["StemEPS"]);
  NumericVector StemAF = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["StemAF"]);
  NumericVector Vsapwood = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["Vsapwood"]); //l·m-2 = mm
  NumericVector LeafPI0 = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["LeafPI0"]);
  NumericVector LeafEPS = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["LeafEPS"]);
  NumericVector LeafAF = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["LeafAF"]);
  NumericVector Vleaf = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["Vleaf"]); //l·m-2 = mm
  
  //Comunication with outside
  DataFrame internalPhenology = Rcpp::as<Rcpp::DataFrame>(x["internalPhenology"]);
  NumericVector phi = Rcpp::as<Rcpp::NumericVector>(internalPhenology["phi"]);
  DataFrame internalWater = Rcpp::as<Rcpp::DataFrame>(x["internalWater"]);
  NumericVector LeafPLCVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["LeafPLC"]);
  NumericVector StemPLCVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["StemPLC"]);
  NumericVector RootCrownPsiVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["RootCrownPsi"]);
  NumericVector StemPsiVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["StemPsi"]);
  NumericVector LeafPsiVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["LeafPsi"]);
  NumericVector StemSympPsiVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["StemSympPsi"]);
  NumericVector LeafSympPsiVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["LeafSympPsi"]);
  List internalLAIDistribution = x["internalLAIDistribution"];
  NumericMatrix LAImx = internalLAIDistribution["live"];
  NumericMatrix LAIme = internalLAIDistribution["expanded"];
  NumericMatrix LAImd = internalLAIDistribution["dead"];
  NumericVector PrevLAIexpanded = internalLAIDistribution["PrevLAIexpanded"];
  NumericVector PrevLAIdead = internalLAIDistribution["PrevLAIdead"];
  NumericVector PARcohort = internalLAIDistribution["PARcohort"];
  
  if(NumericVector::is_na(aspect)) aspect = 0.0;
  if(NumericVector::is_na(slope)) slope = 0.0;
  double latrad = latitude * (M_PI/180.0);
  double asprad = aspect * (M_PI/180.0);
  double slorad = slope * (M_PI/180.0);
  
  //Step in seconds
  double tstep = 86400.0/((double) ntimesteps);

  
 
  //Daily average water vapor pressure at the atmosphere (kPa)
  double vpatm = meteoland::utils_averageDailyVP(tmin, tmax, rhmin,rhmax);
  //If canopy VP is missing or not multilayer initiate it to vpatm
  if(NumericVector::is_na(VPair[0]) || (!multiLayerBalance)){
    for(int i=0;i<ncanlayers;i++) VPair[i] = vpatm;
  }
  //Daily cloud cover
  double cloudcover = 0.0;
  if(prec >0.0) cloudcover = 1.0;
  bool clearday = (prec==0);
  
  ////////////////////////////////////////
  // INITIAL SOIL STATE (from previous step)
  ////////////////////////////////////////
  NumericVector psiSoil = psi(soil, soilFunctions); //Get soil water potential
  NumericMatrix psiSoilM(numCohorts, nlayers);
  if(plantWaterPools){
    //Copy soil water potentials from pools
    List soil_pool = clone(soil);
    NumericVector Ws_pool = soil_pool["W"];
    for(int j = 0; j<numCohorts;j++) {
      //Copy values of soil moisture from pool of cohort j
      for(int l = 0; l<nlayers;l++) Ws_pool[l] = Wpool(j,l);
      //Calculate soil water potential
      psiSoilM(j,_) = psi(soil_pool, soilFunctions);
    }
  }
  
  ////////////////////////////////////////
  // DEFINE OUTPUT
  ////////////////////////////////////////
  //Transpiration and photosynthesis
  NumericMatrix minPsiRhizo = transpOutput["RhizoPsi"];
  
  NumericVector outputFPAR = outputPlants["FPAR"];
  NumericVector SoilExtractCoh = outputPlants["Extraction"];
  NumericVector DDS = outputPlants["DDS"];
  NumericVector LFMC = outputPlants["LFMC"];
  NumericVector Eplant = outputPlants["Transpiration"];
  NumericVector Anplant = outputPlants["GrossPhotosynthesis"];
  NumericVector Agplant = outputPlants["NetPhotosynthesis"];
  NumericVector minStemPsi= outputPlants["StemPsi"];
  NumericVector minRootPsi= outputPlants["RootPsi"];
  NumericVector minLeafPsi= outputPlants["LeafPsiMin"];
  NumericVector maxLeafPsi= outputPlants["LeafPsiMax"];
  NumericVector PLClm = outputPlants["LeafPLC"];
  NumericVector PLCsm = outputPlants["StemPLC"];
  NumericVector dEdPm = outputPlants["dEdP"];
  NumericVector PWB = outputPlants["WaterBalance"];
  NumericVector RWCsm = outputPlants["StemRWC"];
  NumericVector RWClm = outputPlants["LeafRWC"];
  
  
  NumericVector maxGSW_SL = outputSunlit["GSWMax"];
  NumericVector maxGSW_SH = outputShade["GSWMax"];
  NumericVector minGSW_SL = outputSunlit["GSWMin"];
  NumericVector minGSW_SH = outputShade["GSWMin"];
  NumericVector maxTemp_SL = outputSunlit["TempMax"];
  NumericVector maxTemp_SH = outputShade["TempMax"];
  NumericVector minTemp_SL = outputSunlit["TempMin"];
  NumericVector minTemp_SH = outputShade["TempMin"];
  NumericVector maxLeafPsi_SL = outputSunlit["LeafPsiMax"];
  NumericVector maxLeafPsi_SH = outputShade["LeafPsiMax"];
  NumericVector minLeafPsi_SL = outputSunlit["LeafPsiMin"];
  NumericVector minLeafPsi_SH = outputShade["LeafPsiMin"];
  
  NumericMatrix Einst = outputPlantsInst["E"];
  NumericMatrix Aginst = outputPlantsInst["Ag"];
  NumericMatrix Aninst = outputPlantsInst["An"];
  NumericMatrix dEdPInst = outputPlantsInst["dEdP"];
  NumericMatrix LeafPsiInst = outputPlantsInst["LeafPsi"];
  NumericMatrix StemPsiInst = outputPlantsInst["StemPsi"];
  NumericMatrix RootPsiInst = outputPlantsInst["RootPsi"];
  NumericMatrix LeafSympPsiInst = outputPlantsInst["LeafSympPsi"];
  NumericMatrix StemSympPsiInst = outputPlantsInst["StemSympPsi"];
  NumericMatrix StemPLC = outputPlantsInst["StemPLC"];
  NumericMatrix LeafPLC = outputPlantsInst["LeafPLC"];
  NumericMatrix LeafRWCInst = outputPlantsInst["LeafRWC"];
  NumericMatrix StemRWCInst = outputPlantsInst["StemRWC"];
  NumericMatrix LeafSympRWCInst = outputPlantsInst["LeafSympRWC"];
  NumericMatrix StemSympRWCInst = outputPlantsInst["StemSympRWC"];
  NumericMatrix PWBinst = outputPlantsInst["PWB"];

  NumericMatrix LAI_SL = outputSunlitInst["LAI"];
  NumericMatrix LAI_SH = outputShadeInst["LAI"];
  NumericMatrix Vmax298_SL = outputSunlitInst["Vmax298"];
  NumericMatrix Vmax298_SH = outputShadeInst["Vmax298"];
  NumericMatrix Jmax298_SL = outputSunlitInst["Jmax298"];
  NumericMatrix Jmax298_SH = outputShadeInst["Jmax298"];
  NumericMatrix SWR_SL = outputSunlitInst["Abs_SWR"];
  NumericMatrix SWR_SH = outputShadeInst["Abs_SWR"];
  NumericMatrix PAR_SL = outputSunlitInst["Abs_PAR"];
  NumericMatrix PAR_SH = outputShadeInst["Abs_PAR"];
  NumericMatrix LWR_SL = outputSunlitInst["Net_LWR"];
  NumericMatrix LWR_SH = outputShadeInst["Net_LWR"];
  NumericMatrix An_SL = outputSunlitInst["An"];
  NumericMatrix An_SH = outputShadeInst["An"];
  NumericMatrix Ag_SL = outputSunlitInst["Ag"];
  NumericMatrix Ag_SH = outputShadeInst["Ag"];
  NumericMatrix Ci_SL = outputSunlitInst["Ci"];
  NumericMatrix Ci_SH = outputShadeInst["Ci"];
  NumericMatrix E_SL = outputSunlitInst["E"];
  NumericMatrix E_SH = outputShadeInst["E"];
  NumericMatrix GSW_SL = outputSunlitInst["Gsw"];
  NumericMatrix GSW_SH = outputShadeInst["Gsw"];
  NumericMatrix VPD_SL = outputSunlitInst["VPD"];
  NumericMatrix VPD_SH = outputShadeInst["VPD"];
  NumericMatrix Temp_SL = outputSunlitInst["Temp"];
  NumericMatrix Temp_SH = outputShadeInst["Temp"];
  NumericMatrix Psi_SL = outputSunlitInst["Psi"];
  NumericMatrix Psi_SH = outputShadeInst["Psi"];


  IntegerVector iPMSunlit(numCohorts,0), iPMShade(numCohorts,0); //Initial values set to closed stomata
  
  
  //Reset output data
  std::fill(Eplant.begin(), Eplant.end(), 0.0);
  std::fill(Agplant.begin(), Agplant.end(), 0.0);
  std::fill(Anplant.begin(), Anplant.end(), 0.0);
  std::fill(PWB.begin(), PWB.end(), 0.0);
  std::fill(SoilExtractCoh.begin(), SoilExtractCoh.end(), 0.0);
  std::fill(outputExtraction.begin(), outputExtraction.end(), 0.0);
  NumericMatrix soilLayerExtractInst = transpOutput["ExtractionInst"];
  std::fill(soilLayerExtractInst.begin(), soilLayerExtractInst.end(), 0.0);
  if(plantWaterPools) {
    List ExtractionPools = transpOutput["ExtractionPools"];
    for(int c=0;c<numCohorts;c++) {
      NumericMatrix ExtractionPoolsCoh = Rcpp::as<Rcpp::NumericMatrix>(ExtractionPools[c]);
      std::fill(ExtractionPoolsCoh.begin(), ExtractionPoolsCoh.end(), 0.0);
    }
  }
  std::fill(StemPLC.begin(), StemPLC.end(), NA_REAL);
  std::fill(LeafPLC.begin(), LeafPLC.end(), NA_REAL);
  std::fill(dEdPInst.begin(), dEdPInst.end(), NA_REAL);
  std::fill(LeafPsiInst.begin(), LeafPsiInst.end(), NA_REAL);
  std::fill(StemPsiInst.begin(), StemPsiInst.end(), NA_REAL);
  std::fill(RootPsiInst.begin(), RootPsiInst.end(), NA_REAL);
  std::fill(LeafSympPsiInst.begin(), LeafSympPsiInst.end(), NA_REAL);
  std::fill(StemSympPsiInst.begin(), StemSympPsiInst.end(), NA_REAL);
  std::fill(Einst.begin(), Einst.end(), NA_REAL);
  std::fill(Aninst.begin(), Aninst.end(), NA_REAL);
  std::fill(Aginst.begin(), Aginst.end(), NA_REAL);
  std::fill(LeafRWCInst.begin(), LeafRWCInst.end(), NA_REAL);
  std::fill(StemRWCInst.begin(), StemRWCInst.end(), NA_REAL);
  std::fill(PWBinst.begin(), PWBinst.end(), NA_REAL);
  std::fill(LeafSympRWCInst.begin(), LeafSympRWCInst.end(), NA_REAL);
  std::fill(StemSympRWCInst.begin(), StemSympRWCInst.end(), NA_REAL);
  std::fill(Psi_SH.begin(), Psi_SH.end(), NA_REAL);
  std::fill(Psi_SL.begin(), Psi_SL.end(), NA_REAL);
  std::fill(Temp_SL.begin(), Temp_SL.end(), NA_REAL);
  std::fill(Temp_SH.begin(), Temp_SH.end(), NA_REAL);
  std::fill(VPD_SL.begin(), VPD_SL.end(), NA_REAL);
  std::fill(VPD_SH.begin(), VPD_SH.end(), NA_REAL);
  std::fill(GSW_SL.begin(), GSW_SL.end(), NA_REAL);
  std::fill(GSW_SH.begin(), GSW_SH.end(), NA_REAL);
  std::fill(E_SH.begin(), E_SH.end(), NA_REAL);
  std::fill(E_SL.begin(), E_SL.end(), NA_REAL);
  std::fill(Ci_SH.begin(), Ci_SH.end(), NA_REAL);
  std::fill(Ci_SL.begin(), Ci_SL.end(), NA_REAL);
  std::fill(An_SH.begin(), An_SH.end(), NA_REAL);
  std::fill(Ag_SH.begin(), Ag_SH.end(), NA_REAL);
  std::fill(An_SL.begin(), An_SL.end(), NA_REAL);
  std::fill(Ag_SL.begin(), Ag_SL.end(), NA_REAL);
  std::fill(LWR_SL.begin(), LWR_SL.end(), NA_REAL);
  std::fill(PAR_SH.begin(), PAR_SH.end(), NA_REAL);
  std::fill(LWR_SH.begin(), LWR_SH.end(), NA_REAL);
  std::fill(PAR_SL.begin(), PAR_SL.end(), NA_REAL);
  std::fill(Jmax298_SH.begin(), Jmax298_SH.end(), 0.0);
  std::fill(Jmax298_SL.begin(), Jmax298_SL.end(), 0.0);
  std::fill(Vmax298_SH.begin(), Vmax298_SH.end(), 0.0);
  std::fill(Vmax298_SL.begin(), Vmax298_SL.end(), 0.0);
  std::fill(LAI_SH.begin(), LAI_SH.end(), 0.0);
  std::fill(LAI_SL.begin(), LAI_SL.end(), 0.0);
  std::fill(SWR_SL.begin(), SWR_SL.end(), NA_REAL);
  std::fill(SWR_SH.begin(), SWR_SH.end(), NA_REAL);
  
  
  ////////////////////////////////////////
  // STEP 1. Estimate stand-level leaf area values and leaf distribution across layers from leaf-level live/expanded area
  ////////////////////////////////////////
  double LAIcell = 0.0, LAIcelldead = 0.0, LAIcelllive = 0.0, LAIcellexpanded = 0.0;
  double canopyHeight = 100.0; //Minimum canopy height of 1 m
  for(int c=0;c<numCohorts;c++) {
    LAIcell += (LAI[c]+LAIdead[c]);
    LAIcelldead += LAIdead[c];
    LAIcelllive += LAIlive[c];
    LAIcellexpanded +=LAI[c];
    if((canopyHeight<H[c]) && ((LAI[c]+LAIdead[c])>0.0)) canopyHeight = H[c];
  }
  
  if(numCohorts>0) {
    bool recalc_LAI = false;
    if(NumericVector::is_na(PrevLAIexpanded[0]) || NumericVector::is_na(PrevLAIdead[0])) {
      recalc_LAI = true; 
    } else{
      if(sum(abs(LAI - PrevLAIexpanded))>0.001) {
        recalc_LAI = true; 
      } else {
        if(sum(abs(LAIdead - PrevLAIdead))>0.001) recalc_LAI = true;
      }
    }
    if(recalc_LAI) {
      NumericVector z(ncanlayers+1,0.0);
      for(int i=1;i<=ncanlayers;i++) z[i] = z[i-1] + verticalLayerSize;
      for(int i=0; i<numCohorts;i++) {
        PARcohort[i] = availableLight(H[i]*(1.0-(1.0-CR[i])/2.0), H, LAI, LAIdead, kPAR, CR);
        PrevLAIexpanded[i] = LAI[i];
        PrevLAIdead[i] = LAIdead[i];
      }
      //Update LAI distribution if necessary
      updateLAIdistributionVectors(LAIme, z, LAI, H, CR);
      updateLAIdistributionVectors(LAImd, z, LAIdead, H, CR);
      updateLAIdistributionVectors(LAImx, z, LAIlive, H, CR); //Maximum leaf expansion
      //Update LAI profile per layer
      for(int i=0;i<ncanlayers;i++) {
        LAIpx[i] = sum(LAImx(i,_));
        LAIpe[i] = sum(LAIme(i,_));
        LAIpd[i] = sum(LAImd(i,_));
      }
    }
  }

  NumericVector lad = 100.0*(LAIpe + LAIpd)/verticalLayerSize;
  for(int i=0; i<numCohorts;i++) {
    outputFPAR[i] = PARcohort[i];
  }
  
  ////////////////////////////////////////
  // STEP 2. Determine vertical wind speed profile
  ////////////////////////////////////////
  if(NumericVector::is_na(wind)) wind = defaultWindSpeed; //set to default if missing
  wind = std::min(10.0, std::max(wind, 0.1)); //Bound between 0.1 m/s (0.36 km/h)  and 10 m/s (36 km/h)
  NumericVector zWind(ncanlayers,wind), dU(ncanlayers, 0.0), uw(ncanlayers, 0.0);
  if(canopyHeight>0.0) {
    DataFrame canopyTurbulence = as<DataFrame>(transpOutput["CanopyTurbulence"]);
    windCanopyTurbulence_inner(canopyTurbulence, zmid, lad,  canopyHeight, 
                               wind, windMeasurementHeight);
    zWind = canopyTurbulence["u"]; 
    dU = canopyTurbulence["du"];
    uw = canopyTurbulence["uw"];
  } 
  
  ////////////////////////////////////////
  // STEP 3a. Direct and diffuse shorwave radiation for sub-steps
  ////////////////////////////////////////
  DataFrame ddd = meteoland::radiation_directDiffuseDay(solarConstant, latrad, slorad, asprad, delta,
                                                        rad, clearday, ntimesteps);
  transpOutput["RadiationInputInst"] = ddd;
  
  ////////////////////////////////////////
  // STEP 3b. Above-canopy air temperature and long-wave radiation emission for sub-steps
  ////////////////////////////////////////
  NumericVector solarHour = outputTemperatureInst["SolarHour"];
  NumericVector Tatm = outputTemperatureInst["Tatm"];
  NumericVector Tcan = outputTemperatureInst["Tcan"];
  NumericVector Hcansoil = outputSEBinst["Hcansoil"];
  NumericVector Ebalsoil = outputSEBinst["Ebalsoil"];
  NumericVector LEVsoil = outputSEBinst["LEVsoil"];
  NumericVector abs_SWR_soil = outputSEBinst["SWRsoil"];
  NumericVector net_LWR_soil = outputSEBinst["LWRsoil"];
  NumericVector LEFsnow = outputCEBinst["LEFsnow"];
  NumericVector abs_SWR_can = outputCEBinst["SWRcan"];
  NumericVector net_LWR_can = outputCEBinst["LWRcan"];
  NumericVector LEVcan = outputCEBinst["LEVcan"];
  NumericVector Hcan_heat = outputCEBinst["Hcan"];
  NumericVector Ebal = outputCEBinst["Ebalcan"];
  NumericMatrix Tcan_mat = outputEnergyBalance["TemperatureLayers"];
  NumericMatrix VPcan_mat = outputEnergyBalance["VaporPressureLayers"];
  NumericMatrix Tsoil_mat = outputEnergyBalance["SoilTemperature"];


  NumericVector lwdr(ntimesteps, NA_REAL);
  //Daylength in seconds (assuming flat area because we want to model air temperature variation)
  double tauday = meteoland::radiation_daylengthseconds(latrad,0.0,0.0, delta); 
  NumericVector solarHour_ddd = ddd["SolarHour"]; //in radians
  for(int n=0;n<ntimesteps;n++) {
    solarHour[n] = solarHour_ddd[n];
    //From solar hour (radians) to seconds from sunrise
    double Tsunrise = (solarHour[n]*43200.0/M_PI)+ (tauday/2.0) +(tstep/2.0); 
    //Calculate instantaneous temperature and light conditions
    Tatm[n] = temperatureDiurnalPattern(Tsunrise, tmin, tmax, tminPrev, tmaxPrev, tminNext, tauday);
    //Longwave sky diffuse radiation (W/m2)
    lwdr[n] = meteoland::radiation_skyLongwaveRadiation(Tatm[n], vpatm, cloudcover);
  }
  if(NumericVector::is_na(Tair[0])) {//If missing initialize canopy profile with atmospheric air temperature 
    for(int i=0;i<ncanlayers;i++) Tair[i] = Tatm[0];
  }
  if(NumericVector::is_na(Tsoil[0])) {//If missing initialize soil temperature with atmospheric air temperature 
    for(int l=0;l<nlayers; l++) Tsoil[l] = Tatm[0];
  }
  //Take initial canopy air temperature from previous day
  Tcan[0] = sum(Tair*LAIpx)/sum(LAIpx);
  for(int j=0;j<ncanlayers; j++) {
    Tcan_mat(0,j) = Tair[j];
    VPcan_mat(0,j) = VPair[j];
  }
  //Take temperature soil vector 
  for(int l=0;l<nlayers;l++) {
    Tsoil_mat(0,l) = Tsoil[l]; 
  }
  
  ////////////////////////////////////////
  // STEP 3c. Short-wave radiation extinction and absortion for sub-steps
  ////////////////////////////////////////
  List lightExtinctionAbsortion = instantaneousLightExtinctionAbsortion(LAIme, LAImd, LAImx,
                                                                        Beta_p, Beta_q, ClumpingIndex, 
                                                                        alphaSWR, gammaSWR,
                                                                        ddd, ntimesteps, 0.1);
  transpOutput["LightExtinction"] = lightExtinctionAbsortion;
  List sunshade = lightExtinctionAbsortion["sunshade"];
  List abs_PAR_SL_COH_list = sunshade["PAR_SL"];
  List abs_PAR_SH_COH_list = sunshade["PAR_SH"];
  List abs_SWR_SL_COH_list = sunshade["SWR_SL"];
  List abs_SWR_SH_COH_list = sunshade["SWR_SH"];
  List multilayer = lightExtinctionAbsortion["multilayer"];
  List abs_SWR_SL_ML_list = multilayer["SWR_SL"];
  List abs_SWR_SH_ML_list = multilayer["SWR_SH"];
  List fsunlit_list = lightExtinctionAbsortion["fsunlit"];
  NumericVector abs_SWR_can_LEA = lightExtinctionAbsortion["SWR_can"];
  NumericVector abs_SWR_soil_LEA = lightExtinctionAbsortion["SWR_soil"];
  //Copy to output data structures
  double sum_abs_SWR_soil = 0.0, sum_abs_SWR_can = 0.0;
  for(int n=0; n<ntimesteps;n++) {
    abs_SWR_can[n] = abs_SWR_can_LEA[n];
    abs_SWR_soil[n] = abs_SWR_soil_LEA[n];
    sum_abs_SWR_soil += abs_SWR_soil[n];
    sum_abs_SWR_can += abs_SWR_can[n];
  }

  ////////////////////////////////////////
  //  STEP 4. Hydraulics: determine layers where the plant is connected 
  //          and supply functions (Sperry transpiration mode)
  ////////////////////////////////////////
  IntegerVector nlayerscon(numCohorts,0);
  LogicalMatrix layerConnected(numCohorts, nlayers);
  List layerConnectedPools(numCohorts);



  //Average sap fluidity
  double sapFluidityDay = 1.0;
  if(sapFluidityVariation) sapFluidityDay = 1.0/waterDynamicViscosity((tmin+tmax)/2.0);
  
  //Hydraulics: Define supply functions
  List hydraulicNetwork(numCohorts);
  List supplyAboveground(numCohorts);
  for(int c=0;c<numCohorts;c++) {
    
    if(!plantWaterPools) {
      //Determine connected layers (non-zero fine root abundance)
      nlayerscon[c] = 0;
      for(int l=0;l<nlayers;l++) {
        if(V(c,l)>0.0) {
          layerConnected(c,l)= true;
          nlayerscon[c]=nlayerscon[c]+1;
        } else {
          layerConnected(c,l) = false;
        }
      }
      // Rcout<<c<<" "<< nlayerscon[c]<<"\n";
      if(nlayerscon[c]==0) stop("Plant cohort not connected to any soil layer!");
      
      // Copy values from connected layers
      NumericVector Vc = NumericVector(nlayerscon[c]);
      NumericVector VCroot_kmaxc = NumericVector(nlayerscon[c]);
      NumericVector VGrhizo_kmaxc = NumericVector(nlayerscon[c]);
      NumericVector psic = NumericVector(nlayerscon[c]);
      NumericVector VG_nc = NumericVector(nlayerscon[c]);
      NumericVector VG_alphac= NumericVector(nlayerscon[c]);
      int cnt=0;
      for(int l=0;l<nlayers;l++) {
        if(layerConnected(c,l)) {
          Vc[cnt] = V(c,l);
          VCroot_kmaxc[cnt] = VCroot_kmax(c,l);
          VGrhizo_kmaxc[cnt] = VGrhizo_kmax(c,l);
          psic[cnt] = psiSoil[l];
          VG_nc[cnt] = VG_n[l];
          VG_alphac[cnt] = VG_alpha[l];
          cnt++;
        }
      }
      
      //Build supply function networks (Sperry transpiration mode)
      if(transpirationMode=="Sperry") {
        List HN = initSperryNetwork(c,
                                    internalWater, paramsTranspiration, paramsWaterStorage,
                                    VCroot_kmaxc, VGrhizo_kmaxc,
                                    psic, VG_nc, VG_alphac,
                                    control,
                                    sapFluidityDay);
        hydraulicNetwork[c] = HN;
        supply[c] = supplyFunctionNetwork(HN, 0.0, 0.001); 
      } else if(transpirationMode == "Sureau") {
        hydraulicNetwork[c] = initSureauNetwork(c, LAI,
                                                internalWater, 
                                                paramsAnatomy, paramsTranspiration, paramsWaterStorage,
                                                VCroot_kmaxc, VGrhizo_kmaxc,
                                                psic, VG_nc, VG_alphac,
                                                control, sapFluidityDay);
      }
      
    } else {
      //Determine connected layers (non-zero fine root abundance)
      NumericMatrix RHOPcoh = Rcpp::as<Rcpp::NumericMatrix>(RHOP[c]);
      LogicalMatrix layerConnectedCoh(numCohorts, nlayers);
      NumericMatrix RHOPcohDyn(numCohorts, nlayers);
      nlayerscon[c] = 0;
      for(int l=0;l<nlayers;l++) {
        double overlapFactor = Psi2K(psiSoil[l], -1.0, 4.0);
        RHOPcohDyn(c,l) = RHOPcoh(c,l);
        for(int j=0; j<numCohorts;j++) {
          if(j!=c) {
            RHOPcohDyn(j,l) = RHOPcoh(j,l)*overlapFactor;
            RHOPcohDyn(c,l) = RHOPcohDyn(c,l) + (RHOPcoh(j,l) - RHOPcohDyn(j,l));
          } 
        }
      }
      for(int j=0; j<numCohorts;j++) {
        for(int l=0;l<nlayers;l++) {
          if((V(c,l)>0.0) && (RHOPcohDyn(j,l)>0.0)) {
            layerConnectedCoh(j,l)= true;
            nlayerscon[c]=nlayerscon[c] + 1;
          } else {
            layerConnectedCoh(j,l) = false;
          }
        }
      }
      if(nlayerscon[c]==0) stop("Plant cohort not connected to any soil layer!");
      //Store in list
      layerConnectedPools[c] = layerConnectedCoh; 
      
      // Copy values from connected layers
      NumericVector Vc = NumericVector(nlayerscon[c]);
      NumericVector VCroot_kmaxc = NumericVector(nlayerscon[c]);
      NumericVector VGrhizo_kmaxc = NumericVector(nlayerscon[c]);
      NumericVector psic = NumericVector(nlayerscon[c]);
      NumericVector VG_nc = NumericVector(nlayerscon[c]);
      NumericVector VG_alphac= NumericVector(nlayerscon[c]);
      int cnt=0;
      for(int j=0; j<numCohorts;j++) {
        for(int l=0;l<nlayers;l++) {
          if(layerConnectedCoh(j,l)) {
            Vc[cnt] = V(c,l)*RHOPcohDyn(j,l);
            VCroot_kmaxc[cnt] = VCroot_kmax(c,l)*RHOPcohDyn(j,l);
            VGrhizo_kmaxc[cnt] = VGrhizo_kmax(c,l)*RHOPcohDyn(j,l);
            psic[cnt] = psiSoilM(j,l);
            VG_nc[cnt] = VG_n[l];
            VG_alphac[cnt] = VG_alpha[l];
            cnt++;
          }
        }
      }
      //Build supply function networks (Sperry transpiration mode)
      if(transpirationMode == "Sperry") {
        List HN = initSperryNetwork(c,
                                    internalWater, paramsTranspiration, paramsWaterStorage,
                                    VCroot_kmaxc, VGrhizo_kmaxc,
                                    psic, VG_nc, VG_alphac,
                                    control, 
                                    sapFluidityDay);
        hydraulicNetwork[c] = HN;
        supply[c] = supplyFunctionNetwork(HN, 0.0, 0.001); 
      } else if(transpirationMode == "Sureau") {
        hydraulicNetwork[c] = initSureauNetwork(c, LAI,
                                                internalWater, 
                                                paramsAnatomy, paramsTranspiration, paramsWaterStorage,
                                                VCroot_kmaxc, VGrhizo_kmaxc,
                                                psic, VG_nc, VG_alphac,
                                                control, sapFluidityDay);
      }
    }
  }
  
  
  
  ////////////////////////////////
  // Create input and output objects to be filled in inner functions
  ////////////////////////////////

  List innerOutput = List::create(
                             _["Extraction"] = outputExtraction,
                             _["ExtractionPools"] = transpOutput["ExtractionPools"],
                             _["ExtractionInst"] = transpOutput["ExtractionInst"],
                             _["RhizoPsi"] = minPsiRhizo,
                             _["Plants"] = outputPlants,
                             _["SunlitLeaves"] = transpOutput["SunlitLeaves"],
                             _["ShadeLeaves"] = transpOutput["ShadeLeaves"],
                             _["PlantsInst"] = transpOutput["PlantsInst"],
                             _["SunlitLeavesInst"] = transpOutput["SunlitLeavesInst"],
                             _["ShadeLeavesInst"] = transpOutput["ShadeLeavesInst"],
                             _["LightExtinction"] = lightExtinctionAbsortion,
                             _["LWRExtinction"] = lwrExtinctionList,
                             _["SupplyFunctions"] = supply,
                             _["PhotoSunlitFunctions"] = transpOutput["PhotoSunlitFunctions"],
                             _["PhotoShadeFunctions"] = transpOutput["PhotoShadeFunctions"],
                             _["PMSunlitFunctions"] = transpOutput["PMSunlitFunctions"],
                             _["PMShadeFunctions"] = transpOutput["PMShadeFunctions"]);
  
  
  ////////////////////////////////////////
  // STEP 5. Sub-daily (e.g. hourly) loop
  ////////////////////////////////////////
  for(int n=0;n<ntimesteps;n++) { //Time loop
    
    // Determine soil evaporation and snow melt for the corresponding step
    double soilEvapStep = abs_SWR_soil[n]*(soilEvaporation/sum_abs_SWR_soil);
    double snowMeltStep = abs_SWR_soil[n]*(snowMelt/sum_abs_SWR_soil);
    //Canopy evaporation (mm) in the current step and fraction of dry canopy
    double canEvapStep = canopyEvaporation*(abs_SWR_can[n]/sum_abs_SWR_can);
    double f_dry = 1.0;
    if(canEvapStep>0.0) {
      f_dry = 1.0 - std::min(1.0, canopyEvaporation/pet);
    }
    if(sum_abs_SWR_soil==0.0) { // avoid zero sums
      soilEvapStep = 0.0; 
      snowMeltStep = 0.0;
    }
    if(sum_abs_SWR_can==0.0) { // avoid zero sums
      canEvapStep = 0.0;
      f_dry = 1.0;
    }
    
    //Retrieve fraction of sunlit and short-wave radiation absorbed for the current time step
    NumericVector absPAR_SL_COH = abs_PAR_SL_COH_list[n];
    NumericVector absPAR_SH_COH = abs_PAR_SH_COH_list[n];
    NumericVector absSWR_SL_COH = abs_SWR_SL_COH_list[n];
    NumericVector absSWR_SH_COH = abs_SWR_SH_COH_list[n];
    NumericMatrix absSWR_SL_ML = abs_SWR_SL_ML_list[n];
    NumericMatrix absSWR_SH_ML = abs_SWR_SH_ML_list[n];
    NumericVector fsunlit = fsunlit_list[n];
    
    //Leaf area and Vmax/Jmax corresponding to sunlit and shade leaves
    for(int c=0;c<numCohorts;c++) {
      // Rcout<<"cohort "<<c<<":\n";
      //Constant properties through time steps
      NumericVector Vmax298layer(ncanlayers), Jmax298layer(ncanlayers);
      NumericVector SLarealayer(ncanlayers), SHarealayer(ncanlayers);
      double sn =0.0;
      for(int i=(ncanlayers-1);i>=0.0;i--) {
        //Effect of nitrogen concentration decay through the canopy (Improvement: see 10.5194/bg-7-1833-2010)
        double fn = exp(-0.713*(sn+LAIme(i,c)/2.0)/sum(LAIme(_,c)));
        // Rcout<<" l"<<i<<" fsunlit: "<< fsunlit[i]<<" lai: "<< LAIme(i,c)<<" fn: "<< fn <<"\n";
        sn+=LAIme(i,c);
        SLarealayer[i] = LAIme(i,c)*fsunlit[i];
        SHarealayer[i] = LAIme(i,c)*(1.0-fsunlit[i]);
        Vmax298layer[i] = Vmax298[c]*fn;
        Jmax298layer[i] = Jmax298[c]*fn;
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
    IntegerVector iLayerCohort(numCohorts), iLayerSunlit(numCohorts), iLayerShade(numCohorts);
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
        if((hc > zlow[i]) && (hc <=zup[i])) iLayerCohort[c] = i;
        if((hc_sl > zlow[i]) && (hc_sl <=zup[i])) iLayerSunlit[c] = i;
        if((hc_sh > zlow[i]) && (hc_sh <=zup[i])) iLayerShade[c] = i;
      }
      // Rcout << c << " "<< hc_sl<<" "<< iLayerSunlit[c]<< " "<< hc_sh<<" "<< iLayerShade[c]<<"\n";
    }
    
    List innerInput;
    if(transpirationMode =="Sperry") {
      innerInput = List::create(_["Patm"] = Patm,
                                _["zWind"] = zWind,
                                _["f_dry"] = f_dry,
                                _["iLayerCohort"] = iLayerCohort,
                                _["iLayerSunlit"] = iLayerSunlit,
                                _["iLayerShade"] = iLayerShade,
                                _["iPMSunlit"] = iPMSunlit,
                                _["iPMShade"] = iPMShade,
                                _["nlayerscon"] = nlayerscon,
                                _["layerConnected"] = layerConnected,
                                _["layerConnectedPools"] = layerConnectedPools,
                                _["supply"] = supply);
    } else if(transpirationMode =="Sureau") {
      //To do, create initial plant state
      innerInput = List::create(_["Patm"] = Patm,
                                _["zWind"] = zWind,
                                _["f_dry"] = f_dry,
                                _["iLayerCohort"] = iLayerCohort,
                                _["iLayerSunlit"] = iLayerSunlit,
                                _["iLayerShade"] = iLayerShade,
                                _["nlayerscon"] = nlayerscon,
                                _["layerConnected"] = layerConnected,
                                _["layerConnectedPools"] = layerConnectedPools,
                                _["psiSoil"] = psiSoil,
                                _["psiSoilM"] = psiSoilM,
                                _["networks"] = hydraulicNetwork);
    }
    
    
    ////////////////////////////////////////
    // STEP 5.1 Long-wave radiation balance
    ////////////////////////////////////////
    longwaveRadiationSHAW_inner(lwrExtinctionList[n], LAIme, LAImd, LAImx, 
                                lwdr[n], Tsoil[0], Tair);
    List internalLWR = lwrExtinctionList[n];
    net_LWR_soil[n] = internalLWR["Lnet_ground"];
    net_LWR_can[n]= internalLWR["Lnet_canopy"];
    NumericMatrix Lnet_cohort_layer = internalLWR["Lnet_cohort_layer"];

    ////////////////////////////////////////
    // STEP 5.2 Sunlit/shade leaf energy balance, stomatal conductance and plant hydraulics
    ////////////////////////////////////////
    for(int c=0;c<numCohorts;c++) {
      //default values
      dEdPInst(c,n) = NA_REAL;
      Einst(c,n) = 0.0;
      Aginst(c,n) = 0.0;
      Aninst(c,n) = 0.0;
      if(LAI[c]>0.0) {
        PAR_SL(c,n) = absPAR_SL_COH[c];
        PAR_SH(c,n) = absPAR_SH_COH[c];
        SWR_SL(c,n) = absSWR_SL_COH[c];
        SWR_SH(c,n) = absSWR_SH_COH[c];
        LWR_SL(c,n) = 0.0;
        LWR_SH(c,n) = 0.0;
        for(int i=0;i<ncanlayers;i++) {
          LWR_SL(c,n) += Lnet_cohort_layer(i,c)*fsunlit[i];
          LWR_SH(c,n) += Lnet_cohort_layer(i,c)*(1.0 - fsunlit[i]);
        }
      }
    }

    if(transpirationMode == "Sperry") {
      innerSperry(x, innerInput, innerOutput, n, tstep, 
                  verbose, stepFunctions);
    } else if(transpirationMode == "Sureau"){
      innerSureau(x, innerInput, innerOutput, n, tstep,
                   verbose);
    }
    
    for(int c=0;c<numCohorts;c++) {
      if(LAIlive[c]>0.0 && (LeafPLCVEC[c] < 0.999)) {
        //Store (for output) instantaneous leaf, stem and root potential, plc and rwc values
        StemPLC(c,n) = StemPLCVEC[c];
        LeafPLC(c,n) = LeafPLCVEC[c];
        StemSympRWCInst(c,n) = symplasticRelativeWaterContent(StemSympPsiVEC[c], StemPI0[c], StemEPS[c]);
        LeafSympRWCInst(c,n) = symplasticRelativeWaterContent(LeafSympPsiVEC[c], LeafPI0[c], LeafEPS[c]);
        StemRWCInst(c,n) = StemSympRWCInst(c,n)*(1.0 - StemAF[c]) + (1.0 - StemPLCVEC[c])*StemAF[c];
        LeafRWCInst(c,n) = LeafSympRWCInst(c,n)*(1.0 - LeafAF[c]) + (1.0 - LeafPLCVEC[c])*LeafAF[c];
        StemPsiInst(c,n) = StemPsiVEC[c]; 
        LeafPsiInst(c,n) = LeafPsiVEC[c]; //Store instantaneous (average) leaf potential
        RootPsiInst(c,n) = RootCrownPsiVEC[c]; //Store instantaneous root crown potential
        LeafSympPsiInst(c,n) = LeafSympPsiVEC[c];
        StemSympPsiInst(c,n) = StemSympPsiVEC[c];
        
        if(n==0) {
          minGSW_SL[c] = GSW_SL(c,n);
          minGSW_SH[c] = GSW_SH(c,n);
          maxGSW_SL[c] = GSW_SL(c,n);
          maxGSW_SH[c] = GSW_SH(c,n);
          minTemp_SL[c] = Temp_SL(c,n);
          minTemp_SH[c] = Temp_SH(c,n);
          maxTemp_SL[c] = Temp_SL(c,n);
          maxTemp_SH[c] = Temp_SH(c,n);
          minLeafPsi_SL[c] = Psi_SL(c,n);
          minLeafPsi_SH[c] = Psi_SH(c,n);
          maxLeafPsi_SL[c] = Psi_SL(c,n);
          maxLeafPsi_SH[c] = Psi_SH(c,n);
          minLeafPsi[c] = LeafPsiInst(c,n);
          maxLeafPsi[c] = LeafPsiInst(c,n);
          minStemPsi[c] = StemPsiInst(c,n);
          minRootPsi[c] = RootPsiInst(c,n);
          for(int l=0;l<nlayers;l++) {
            minPsiRhizo(c,l) = RhizoPsiMAT(c,l);
          }
        } else {
          minGSW_SL[c] = std::min(minGSW_SL[c], GSW_SL(c,n));
          minGSW_SH[c] = std::min(minGSW_SH[c], GSW_SH(c,n));
          maxGSW_SL[c] = std::max(maxGSW_SL[c], GSW_SL(c,n));
          maxGSW_SH[c] = std::max(maxGSW_SH[c], GSW_SH(c,n));
          minTemp_SL[c] = std::min(minTemp_SL[c], Temp_SL(c,n));
          minTemp_SH[c] = std::min(minTemp_SH[c], Temp_SH(c,n));
          maxTemp_SL[c] = std::max(maxTemp_SL[c], Temp_SL(c,n));
          maxTemp_SH[c] = std::max(maxTemp_SH[c], Temp_SH(c,n));
          minLeafPsi_SL[c] = std::min(minLeafPsi_SL[c], Psi_SL(c,n));
          minLeafPsi_SH[c] = std::min(minLeafPsi_SH[c], Psi_SH(c,n));
          maxLeafPsi_SL[c] = std::max(maxLeafPsi_SL[c], Psi_SL(c,n));
          maxLeafPsi_SH[c] = std::max(maxLeafPsi_SH[c], Psi_SH(c,n));
          minLeafPsi[c] = std::min(minLeafPsi[c], LeafPsiInst(c,n));
          maxLeafPsi[c] = std::max(maxLeafPsi[c], LeafPsiInst(c,n));
          minStemPsi[c] = std::min(minStemPsi[c], StemPsiInst(c,n));
          minRootPsi[c] = std::min(minRootPsi[c], RootPsiInst(c,n));
          for(int l=0;l<nlayers;l++) {
            minPsiRhizo(c,l) = std::min(minPsiRhizo(c,l), RhizoPsiMAT(c,l));
          }
        }
      } else {
        // Assume constant PLC (so that it can be decreased in the future)
        StemPLC(c,n) = StemPLCVEC[c];
        LeafPLC(c,n) = LeafPLCVEC[c];
        StemPsiInst(c,n) = StemPsiVEC[c]; 
        LeafPsiInst(c,n) = LeafPsiVEC[c]; //Store instantaneous (average) leaf potential
        RootPsiInst(c,n) = RootCrownPsiVEC[c]; //Store instantaneous root crown potential
        LeafSympPsiInst(c,n) = LeafSympPsiVEC[c];
        StemSympPsiInst(c,n) = StemSympPsiVEC[c];
        minLeafPsi[c] = LeafPsiVEC[c];
        maxLeafPsi[c] = LeafPsiVEC[c];
        minStemPsi[c] =  StemPsiVEC[c];
        minRootPsi[c] = RootCrownPsiVEC[c];
        minLeafPsi_SL[c] = LeafPsiVEC[c];
        minLeafPsi_SH[c] = LeafPsiVEC[c];
        maxLeafPsi_SL[c] = LeafPsiVEC[c];
        maxLeafPsi_SH[c] = LeafPsiVEC[c];
        StemSympRWCInst(c,n) = symplasticRelativeWaterContent(StemSympPsiVEC[c], StemPI0[c], StemEPS[c]);
        LeafSympRWCInst(c,n) = symplasticRelativeWaterContent(LeafSympPsiVEC[c], LeafPI0[c], LeafEPS[c]);
        StemRWCInst(c,n) = StemSympRWCInst(c,n)*(1.0 - StemAF[c]) + (1.0 - StemPLCVEC[c])*StemAF[c];
        LeafRWCInst(c,n) = LeafSympRWCInst(c,n)*(1.0 - LeafAF[c]) + (1.0 - LeafPLCVEC[c])*LeafAF[c];
        
      }
    }
    
    ////////////////////////////////////////
    // STEP 5.3 Soil and canopy energy balances (single or multiple canopy layers)
    ////////////////////////////////////////
    
    //Soil latent heat (soil evaporation)
    //Latent heat (snow fusion) as J/m2/s
    if(snowpack>0.0) {
      abs_SWR_soil[n] = 0.0; //Set SWR absorbed by soil to zero (for energy balance) if snow pack is present
      net_LWR_soil[n] = 0.0; //Set net LWR to zero
    } 
    LEVsoil[n] = (1e6)*meteoland::utils_latentHeatVaporisation(Tsoil[0])*soilEvapStep/tstep;
    // Rcout<<n<<" "<<sum_abs_SWR_soil<<" "<<soilEvapStep << " "<<Tsoil[0]<<" " << LEVsoil[n]<<"\n";
    LEFsnow[n] = (1e6)*(snowMeltStep*0.33355)/tstep; // 0.33355 = latent heat of fusion

    //Herbaceous transpiration (mm) in the current step
    double herbTranspStep = herbTranspiration*(abs_SWR_can[n]/sum_abs_SWR_can);
    
    //Canopy convective heat exchange
    double RAcan = aerodynamicResistance(canopyHeight,std::max(wind,1.0)); //Aerodynamic resistance to convective heat transfer
    Hcan_heat[n] = (meteoland::utils_airDensity(Tatm[n],Patm)*Cp_JKG*(Tcan[n]-Tatm[n]))/RAcan;
    
    if(!multiLayerBalance) {//Canopy balance assuming a single layer
      //Soil-canopy turbulent heat exchange
      double wind2m = windSpeedMassmanExtinction(200.0, wind, LAIcell, canopyHeight);
      double RAsoil = aerodynamicResistance(200.0, std::max(wind2m,1.0)); //Aerodynamic resistance to convective heat transfer from soil
      if(snowpack==0.0) {
        Hcansoil[n] = (meteoland::utils_airDensity(Tcan[n],Patm)*Cp_JKG*(Tcan[n]-Tsoil[0]))/RAsoil;
      } else {
        Hcansoil[n] = (meteoland::utils_airDensity(Tcan[n],Patm)*Cp_JKG*(Tcan[n] - 0.0))/RAsoil; //Assumes a zero degree for soil surface (snow)
      } 
      //Latent heat (evaporation + transpiration)
      double sum_Einst_n = 0.0;
      for(int c=0;c<numCohorts;c++) sum_Einst_n +=Einst(c, n);
      double LEwat = (1e6)*meteoland::utils_latentHeatVaporisation(Tcan[n])*(sum_Einst_n + canEvapStep + herbTranspStep)/tstep;
      LEVcan[n] = LEwat; 
      //Canopy temperature changes
      Ebal[n] = abs_SWR_can[n]+ net_LWR_can[n] - LEVcan[n] - LEFsnow[n] - Hcan_heat[n] - Hcansoil[n];
      double canopyAirThermalCapacity = meteoland::utils_airDensity(Tcan[n],Patm)*Cp_JKG;
      double canopyThermalCapacity =  canopyAirThermalCapacity + (0.5*(0.8*LAIcelllive + 1.2*LAIcell) + LAIcelldead)*thermalCapacityLAI; //Avoids zero capacity for winter deciduous
      double Tcannext = Tcan[n]+ std::max(-3.0, std::min(3.0, tstep*Ebal[n]/canopyThermalCapacity)); //Avoids changes in temperature that are too fast
      if(n<(ntimesteps-1)) Tcan[n+1] = Tcannext;
      for(int i=0;i<ncanlayers;i++) Tair[i] = Tcannext;
      
      
      //Soil energy balance including exchange with canopy
      if(snowpack==0.0) {
        Ebalsoil[n] = abs_SWR_soil[n] + net_LWR_soil[n] + Hcansoil[n] - LEVsoil[n]; //Here we use all energy escaping to atmosphere
      } else {
        //Heat conduction between soil and snow
        double Hcond_snow = 0.05*(0.0 - Tsoil[0]);//0.05 Wm-1K-1 for fresh snow
        Ebalsoil[n] = Hcond_snow - net_LWR_soil[n] - LEVsoil[n]; 
      }
      
      //Soil temperature changes
      NumericVector soilTchange = temperatureChange(widths, Tsoil, sand, clay, Ws, Theta_SAT, Theta_FC, Ebalsoil[n], tstep);
      for(int l=0;l<nlayers;l++) {
        Tsoil[l] = Tsoil[l] + std::max(-3.0, std::min(3.0, soilTchange[l])); 
        if(n<(ntimesteps-1)) Tsoil_mat(n+1,l)= Tsoil[l];
      }

    } else { //Multilayer canopy balance
      double moistureAtm = 0.622*(vpatm/Patm)*meteoland::utils_airDensity(Tatm[n],Patm);
      double CO2Atm = 0.409*Catm*44.01; //mg/m3
        
      double tsubstep = tstep/((double) nsubsteps); 
      double maxTchange = 3.0/((double) nsubsteps);
      double maxMoistureChange = 0.001/((double)nsubsteps); //=0.16 kPa per step
      double maxCO2Change = 180.0/((double)nsubsteps); //= 10 ppm per step
      double deltaZ = (verticalLayerSize/100.0); //Vertical layer size in m
      DataFrame LWR_layer = Rcpp::as<Rcpp::DataFrame>(internalLWR["LWR_layer"]);
      NumericVector LWRnet_layer = LWR_layer["Lnet"];
      Ebal[n] = 0.0;
      LEVcan[n] = 0.0;
      NumericVector Tairnext(ncanlayers), LElayer(ncanlayers), absSWRlayer(ncanlayers), Rnlayer(ncanlayers), Hleaflayer(ncanlayers);
      NumericVector layerThermalCapacity(ncanlayers);
      NumericVector moistureET(ncanlayers), rho(ncanlayers), moistureLayer(ncanlayers), moistureLayernext(ncanlayers);
      NumericVector CO2An(ncanlayers), CO2Layer(ncanlayers), CO2Layernext(ncanlayers);
      for(int i=0;i<ncanlayers;i++) {
        rho[i] = meteoland::utils_airDensity(Tair[i],Patm);
        absSWRlayer[i] = sum(absSWR_SL_ML(i,_)) + sum(absSWR_SH_ML(i,_));
        //Radiation balance
        Rnlayer[i] = absSWRlayer[i] + LWRnet_layer[i];
        NumericVector pLayer = LAIme(i,_)/LAI; //Proportion of each cohort LAI in layer i
        //Instantaneous layer transpiration
        //from mmolH2O/m2/s to kgH2O/m2/s
        double ElayerInst = 0.001*0.01802*sum(LAIme(i,_)*(E_SL(_,n)*fsunlit[i] + E_SH(_,n)*(1.0-fsunlit[i])));
        //Assumes Layers contribute to evaporation proportionally to their LAI fraction
        double layerEvapInst = (canEvapStep/tstep)*(LAIpe[i]/LAIcellexpanded);
        //Instantaneous herbaceous transpiration (for bottom layer)
        double herbTranspInst = 0.0;
        if(i==0) herbTranspInst = (herbTranspStep/tstep);
        //Estimate instantaneous mgCO2/m2 absorption for the layer, taking into account the proportion of sunlit and shade leaves of each cohort
        //from micro.molCO2/m2/s to mgCO2/m2/s
        double Anlayer =(1e-3)*44.01*sum(LAIme(i,_)*(An_SL(_,n)*fsunlit[i] + An_SH(_,n)*(1.0-fsunlit[i])));
        // 1000.0*(44.01/12.0)*sum(Aninst(_,n)*pLayer); 
        double LEwat = (1e6)*meteoland::utils_latentHeatVaporisation(Tair[i])*(ElayerInst + layerEvapInst+ herbTranspInst);
        LElayer[i] = LEwat; //Energy spent in vaporisation
        if(i==0) LElayer[i] = LElayer[i] - LEFsnow[n]; //Add latent heat of fusion to first layer
        LEVcan[n] = LElayer[i];
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
        double RAsoil = aerodynamicResistance(200.0, std::max(zWind[0],1.0)); //Aerodynamic resistance to convective heat transfer from soil
        double Hcansoils = Cp_JKG*meteoland::utils_airDensity(Tair[0],Patm)*(Tair[0]-Tsoil[0])/RAsoil;
        // double Hcan_heats = (meteoland::utils_airDensity(Tatm[n],Patm)*Cp_JKG*(Tair[ncanlayers-1]-Tatm[n]))/RAcan;
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
            double gHa = 0.189*pow(std::max(zWind[i],0.1)/(leafWidth[c]*0.0072), 0.5);
            double Hsunlit = 2.0*Cp_Jmol*rho[i]*(Temp_SL(c, n)-Tair[i])*gHa;
            double Hshade = 2.0*Cp_Jmol*rho[i]*(Temp_SH(c, n)-Tair[i])*gHa;
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
        if(snowpack==0.0) {
          Ebalsoils = abs_SWR_soil[n] + net_LWR_soil[n]  - LEVsoil[n] + Hcansoils; //Here we use all energy escaping to atmosphere
        } else {
          //Heat conduction between soil and snow
          double Hcond_snow = 0.05*(0.0 - Tsoil[0]);//0.05 Wm-1K-1 for fresh snow
          Ebalsoils = Hcond_snow - net_LWR_soil[n] - LEVsoil[n]; 
        }
        Ebalsoil[n] +=Ebalsoils;
        Hcansoil[n] +=Hcansoils;
        //Soil temperature changes
        NumericVector soilTchange = temperatureChange(widths, Tsoil, sand, clay, Ws, Theta_SAT, Theta_FC, Ebalsoils, tsubstep);
        for(int l=0;l<nlayers;l++) Tsoil[l] = Tsoil[l] + soilTchange[l];
      }
      Hcansoil[n] = Hcansoil[n]/((double) nsubsteps);
      Ebalsoil[n] = Ebalsoil[n]/((double) nsubsteps);
      for(int i=0;i<ncanlayers;i++) {
        VPair[i] = ((moistureLayer[i]/rho[i])*Patm)/0.622;
        Cair[i] = CO2Layer[i]/(0.409*44.01);
        // Rcout<< n << " "<<i << " - " << moistureLayer[i]<< " "<< VPair[i]<<"\n";
      }
      // stop("kk");
      // Rcout<< n << " "<<abs_SWR_can[n]<< " = "<< sum(absSWRlayer)<<" = "<< (sum(absSWR_SL_ML) + sum(absSWR_SH_ML))<<" "<<net_LWR_can[n]<< " = "<< sum(LWRnet_layer)<<"\n";
      //Canopy energy balance
      Ebal[n] = abs_SWR_can[n]+ net_LWR_can[n] - LEVcan[n] - Hcan_heat[n] - Hcansoil[n];
      if(n<(ntimesteps-1)) {
        Tcan[n+1] = sum(Tair*LAIpx)/sum(LAIpx); 
        Tsoil_mat(n+1,_)= Tsoil;
      }
    }
    if(n<(ntimesteps-1)) for(int i=0;i<ncanlayers;i++) {
      Tcan_mat(n+1,i) = Tair[i];
      VPcan_mat(n+1,i) = VPair[i];
    }
  } //End of timestep loop

  ////////////////////////////////////////
  // STEP 6. Plant drought stress (relative whole-plant conductance), cavitation and live fuel moisture
  ////////////////////////////////////////
  for(int c=0;c<numCohorts;c++) {
    SoilExtractCoh[c] = 0.0;
    for(int l=0;l<nlayers;l++) {
      SoilExtractCoh[c] += outputExtraction(c,l);
    }
    PLCsm[c] = 0.0;
    PLClm[c] = 0.0;
    RWCsm[c] = 0.0;
    RWClm[c] = 0.0;
    dEdPm[c] = 0.0;
    for(int n=0;n<ntimesteps;n++) {
      PLCsm[c] += StemPLC(c,n);
      PLClm[c] += LeafPLC(c,n);
      RWCsm[c] += StemRWCInst(c,n);
      RWClm[c] += LeafRWCInst(c,n);
      dEdPm[c] += dEdPInst(c,n);
    }
    PLCsm[c] = PLCsm[c]/((double) ntimesteps);
    PLClm[c] = PLClm[c]/((double) ntimesteps);
    RWCsm[c] = RWCsm[c]/((double) ntimesteps);
    RWClm[c] = RWClm[c]/((double) ntimesteps);
    dEdPm[c] = dEdPm[c]/((double) ntimesteps);
    if(lfmcComponent=="fine") {
      LFMC[c] = maxFMC[c]*((1.0/r635[c])*RWClm[c]+(1.0 - (1.0/r635[c]))*RWCsm[c]);
    } else { //leaf
      LFMC[c] = maxFMC[c]*RWClm[c];
    }
    DDS[c] = (1.0 - (dEdPm[c]/(sapFluidityDay*Plant_kmax[c])));
    if(phenoType[c] == "winter-deciduous" || phenoType[c] == "winter-semideciduous") {
      DDS[c] = phi[c]*DDS[c];
      if(phi[c] == 0.0) DDS[c] = 0.0;
    }
    
    double SAmax = 10e4/Al2As[c]; //cm2·m-2 of leaf area
    double r = cavitationRecoveryMaximumRate*std::max(0.0, (StemSympPsiVEC[c] + 1.5)/1.5);
    if(stemCavitationRecovery=="rate") {
      StemPLCVEC[c] = std::max(0.0, StemPLCVEC[c] - (r/SAmax));
    }
    if(leafCavitationRecovery=="rate") {
      LeafPLCVEC[c] = std::max(0.0, LeafPLCVEC[c] - (r/SAmax));
    }
  }
  

  // Copy output stand
  outputStand["LAI"] = LAIcell;
  outputStand["LAIlive"] = LAIcelllive;
  outputStand["LAIexpanded"] = LAIcellexpanded;
  outputStand["LAIdead"] = LAIcelldead;
  
  NumericVector outputLAI = outputPlants["LAI"];
  NumericVector outputLAIlive = outputPlants["LAIlive"];
  for(int c =0;c<numCohorts;c++) {
    outputLAI[c] = LAI[c];
    outputLAIlive[c] = LAIlive[c];
  }
}

//' @rdname transp_modes
//' 
//' @param canopyEvaporation Canopy evaporation (from interception) for \code{day} (mm).
//' @param snowMelt Snow melt values  for \code{day} (mm).
//' @param soilEvaporation Bare soil evaporation for \code{day} (mm).
//' @param herbTranspiration Transpiration of herbaceous plants for \code{day} (mm).
//' @param stepFunctions An integer to indicate a simulation step for which photosynthesis and profit maximization functions are desired.
//' 
//' @keywords internal
// [[Rcpp::export("transp_transpirationSperry")]]
List transpirationSperry(List x, DataFrame meteo, int day,
                        double latitude, double elevation, double slope, double aspect,
                        double canopyEvaporation = 0.0, double snowMelt = 0.0, double soilEvaporation = 0.0, double herbTranspiration = 0.0,
                        int stepFunctions = NA_INTEGER, 
                        bool modifyInput = true) {

  List control = x["control"];
  String transpirationMode = control["transpirationMode"];
  if(transpirationMode != "Sperry") stop("Transpiration mode in 'x' must be 'Sperry'");
  
  if(!meteo.containsElementNamed("MinTemperature")) stop("Please include variable 'MinTemperature' in weather input.");
  NumericVector MinTemperature = meteo["MinTemperature"];
  if(!meteo.containsElementNamed("MaxTemperature")) stop("Please include variable 'MaxTemperature' in weather input.");
  NumericVector MaxTemperature = meteo["MaxTemperature"];
  if(!meteo.containsElementNamed("MinRelativeHumidity")) stop("Please include variable 'MinRelativeHumidity' in weather input.");
  NumericVector MinRelativeHumidity = meteo["MinRelativeHumidity"];
  if(!meteo.containsElementNamed("MaxRelativeHumidity")) stop("Please include variable 'MaxRelativeHumidity' in weather input.");
  NumericVector MaxRelativeHumidity = meteo["MaxRelativeHumidity"];
  if(!meteo.containsElementNamed("Radiation")) stop("Please include variable 'Radiation' in weather input.");
  NumericVector Radiation = meteo["Radiation"];
  if(!meteo.containsElementNamed("Precipitation")) stop("Please include variable 'Precipitation' in weather input.");
  NumericVector Precipitation = meteo["Precipitation"];
  NumericVector WindSpeed(Precipitation.length(), NA_REAL);
  if(meteo.containsElementNamed("WindSpeed")) WindSpeed = meteo["WindSpeed"];
  NumericVector CO2(Precipitation.length(), NA_REAL);
  if(meteo.containsElementNamed("CO2")) CO2 = meteo["CO2"];
  NumericVector Patm(MinTemperature.length(), NA_REAL);
  if(meteo.containsElementNamed("Patm")) Patm = meteo["Patm"];

  CharacterVector dateStrings = getWeatherDates(meteo);
  std::string c = as<std::string>(dateStrings[day-1]);
  int J = meteoland::radiation_julianDay(std::atoi(c.substr(0, 4).c_str()),std::atoi(c.substr(5,2).c_str()),std::atoi(c.substr(8,2).c_str()));
  double delta = meteoland::radiation_solarDeclination(J);
  double solarConstant = meteoland::radiation_solarConstant(J);
  double latrad = latitude * (M_PI/180.0);
  if(NumericVector::is_na(aspect)) aspect = 0.0;
  if(NumericVector::is_na(slope)) slope = 0.0;
  double asprad = aspect * (M_PI/180.0);
  double slorad = slope * (M_PI/180.0);
  
  double prec = Precipitation[day-1];
  double rad = Radiation[day-1];
  double tmax = MaxTemperature[day-1];
  double tmin = MinTemperature[day-1];
  double tmaxPrev = tmax;
  double tminPrev = tmin;
  double tminNext = tmin;
  if(day>1) {
    tmaxPrev = MaxTemperature[day-2];
    tminPrev = MinTemperature[day-2];
  }
  if(day<(MaxTemperature.length()-1)) tminNext = MinTemperature[day];
  double rhmax = MaxRelativeHumidity[day-1];
  double rhmin = MinRelativeHumidity[day-1];
  double wind = WindSpeed[day-1];
  double Catm = CO2[day-1];
  if(NumericVector::is_na(Catm)) Catm = control["defaultCO2"];
  
  double pet = meteoland::penman(latrad, elevation, slorad, asprad, J, tmin, tmax, rhmin, rhmax, rad, wind);
  
  NumericVector meteovec = NumericVector::create(
    Named("tmin") = tmin, 
    Named("tmax") = tmax,
    Named("tminPrev") = tminPrev, 
    Named("tmaxPrev") = tmaxPrev, 
    Named("tminNext") = tminNext, 
    Named("prec") = prec,
    Named("rhmin") = rhmin, 
    Named("rhmax") = rhmax, 
    Named("rad") = rad, 
    Named("wind") = wind, 
    Named("Catm") = Catm,
    Named("Patm") = Patm[day-1],
    Named("pet") = pet);
  
  DataFrame cohorts = Rcpp::as<Rcpp::DataFrame>(x["cohorts"]);
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  DataFrame soil = Rcpp::as<Rcpp::DataFrame>(x["soil"]);
  DataFrame canopyParams = Rcpp::as<Rcpp::DataFrame>(x["canopy"]);
  int ncanlayers = canopyParams.nrow(); //Number of canopy layers
  int nlayers = soil.nrow();
  int numCohorts = cohorts.nrow();
  int ntimesteps = control["ndailysteps"];
  List transpOutput = advancedTranspirationCommunicationOutput(numCohorts, nlayers, ncanlayers, ntimesteps);
  
  transpirationAdvanced(transpOutput, x, meteovec,
                        latitude, elevation, slope, aspect,
                        solarConstant, delta,
                        canopyEvaporation, snowMelt, soilEvaporation, herbTranspiration,
                        false, stepFunctions, 
                        modifyInput);
  
  List transpAdvanced = copyAdvancedTranspirationOutput(transpOutput, x);
  
  return(transpAdvanced);
} 

//' @rdname transp_modes
//' @keywords internal
// [[Rcpp::export("transp_transpirationSureau")]]
List transpirationSureau(List x, DataFrame meteo, int day,
                         double latitude, double elevation, double slope, double aspect,
                         double canopyEvaporation = 0.0, double snowMelt = 0.0, double soilEvaporation = 0.0, double herbTranspiration = 0.0,
                         bool modifyInput = true) {
  
  List control = x["control"];
  String transpirationMode = control["transpirationMode"];
  if(transpirationMode != "Sureau") stop("Transpiration mode in 'x' must be 'Sureau'");
  if(!meteo.containsElementNamed("MinTemperature")) stop("Please include variable 'MinTemperature' in weather input.");
  NumericVector MinTemperature = meteo["MinTemperature"];
  if(!meteo.containsElementNamed("MaxTemperature")) stop("Please include variable 'MaxTemperature' in weather input.");
  NumericVector MaxTemperature = meteo["MaxTemperature"];
  if(!meteo.containsElementNamed("MinRelativeHumidity")) stop("Please include variable 'MinRelativeHumidity' in weather input.");
  NumericVector MinRelativeHumidity = meteo["MinRelativeHumidity"];
  if(!meteo.containsElementNamed("MaxRelativeHumidity")) stop("Please include variable 'MaxRelativeHumidity' in weather input.");
  NumericVector MaxRelativeHumidity = meteo["MaxRelativeHumidity"];
  if(!meteo.containsElementNamed("Radiation")) stop("Please include variable 'Radiation' in weather input.");
  NumericVector Radiation = meteo["Radiation"];
  if(!meteo.containsElementNamed("Precipitation")) stop("Please include variable 'Precipitation' in weather input.");
  NumericVector Precipitation = meteo["Precipitation"];
  NumericVector WindSpeed(Precipitation.length(), NA_REAL);
  if(meteo.containsElementNamed("WindSpeed")) WindSpeed = meteo["WindSpeed"];
  NumericVector CO2(Precipitation.length(), NA_REAL);
  if(meteo.containsElementNamed("CO2")) CO2 = meteo["CO2"];
  NumericVector Patm(MinTemperature.length(), NA_REAL);
  if(meteo.containsElementNamed("Patm")) Patm = meteo["Patm"];
  
  CharacterVector dateStrings = getWeatherDates(meteo);
  std::string c = as<std::string>(dateStrings[day-1]);
  int J = meteoland::radiation_julianDay(std::atoi(c.substr(0, 4).c_str()),std::atoi(c.substr(5,2).c_str()),std::atoi(c.substr(8,2).c_str()));
  double delta = meteoland::radiation_solarDeclination(J);
  double solarConstant = meteoland::radiation_solarConstant(J);
  double latrad = latitude * (M_PI/180.0);
  if(NumericVector::is_na(aspect)) aspect = 0.0;
  if(NumericVector::is_na(slope)) slope = 0.0;
  double asprad = aspect * (M_PI/180.0);
  double slorad = slope * (M_PI/180.0);
  
  double prec = Precipitation[day-1];
  double rad = Radiation[day-1];
  double tmax = MaxTemperature[day-1];
  double tmin = MinTemperature[day-1];
  double tmaxPrev = tmax;
  double tminPrev = tmin;
  double tminNext = tmin;
  if(day>1) {
    tmaxPrev = MaxTemperature[day-2];
    tminPrev = MinTemperature[day-2];
  }
  if(day<(MaxTemperature.length()-1)) tminNext = MinTemperature[day];
  double rhmax = MaxRelativeHumidity[day-1];
  double rhmin = MinRelativeHumidity[day-1];
  double wind = WindSpeed[day-1];
  double Catm = CO2[day-1];
  if(NumericVector::is_na(Catm)) Catm = control["defaultCO2"];
  
  double pet = meteoland::penman(latrad, elevation, slorad, asprad, J, tmin, tmax, rhmin, rhmax, rad, wind);
  
  NumericVector meteovec = NumericVector::create(
    Named("tmin") = tmin, 
    Named("tmax") = tmax,
    Named("tminPrev") = tminPrev, 
    Named("tmaxPrev") = tmaxPrev, 
    Named("tminNext") = tminNext, 
    Named("prec") = prec,
    Named("rhmin") = rhmin, 
    Named("rhmax") = rhmax, 
    Named("rad") = rad, 
    Named("wind") = wind, 
    Named("Catm") = Catm,
    Named("Patm") = Patm[day-1],
    Named("pet") = pet);
  
  
  DataFrame cohorts = Rcpp::as<Rcpp::DataFrame>(x["cohorts"]);
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  DataFrame soil = Rcpp::as<Rcpp::DataFrame>(x["soil"]);
  DataFrame canopyParams = Rcpp::as<Rcpp::DataFrame>(x["canopy"]);
  int ncanlayers = canopyParams.nrow(); //Number of canopy layers
  int nlayers = soil.nrow();
  int numCohorts = cohorts.nrow();
  int ntimesteps = control["ndailysteps"];
  List transpOutput = advancedTranspirationCommunicationOutput(numCohorts, nlayers, ncanlayers, ntimesteps);
  transpirationAdvanced(transpOutput, x, meteovec,
                                            latitude, elevation, slope, aspect,
                                            solarConstant, delta,
                                            canopyEvaporation, snowMelt, soilEvaporation, herbTranspiration,
                                            false, NA_INTEGER, 
                                            modifyInput);
  
  List transpAdvanced = copyAdvancedTranspirationOutput(transpOutput, x);
  
  return(transpAdvanced);
} 




