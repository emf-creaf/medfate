// [[Rcpp::interfaces(r,cpp)]]
#define STRICT_R_HEADERS
#include <Rcpp.h>
#include <numeric>
#include <math.h>
#include "communication_structures.h"
#include "lightextinction_basic.h"
#include "lightextinction_advanced.h"
#include "windextinction.h"
#include "hydraulics.h"
#include "hydrology.h"
#include "biophysicsutils.h"
#include "forestutils.h"
#include "photosynthesis.h"
#include "phenology.h"
#include "transpiration.h"
#include "fuelstructure.h"
#include "firebehaviour.h"
#include "tissuemoisture.h"
#include "soil.h"
#include <meteoland.h>
using namespace Rcpp;


void fccsHazard(NumericVector fireHazard, List x, NumericVector meteovec, List transpOutput, double slope) {
  List control = x["control"];
  
  double fireHazardStandardWind = control["fireHazardStandardWind"];
  double fireHazardStandardDFMC = control["fireHazardStandardDFMC"];
  
  DataFrame FCCSprops = Rcpp::as<Rcpp::DataFrame>(x["internalFCCS"]);
  DataFrame Plants = Rcpp::as<Rcpp::DataFrame>(transpOutput["Plants"]);
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  
  
  double fm_dead = NA_REAL;
  if(!NumericVector::is_na(fireHazardStandardDFMC)) {
    fm_dead = fireHazardStandardDFMC;
  } else {
    double tmin = meteovec["tmin"];
    double tmax = meteovec["tmax"];
    double rhmin = meteovec["rhmin"];
    double rhmax = meteovec["rhmax"];
    // Estimate moisture of dead fine fuels (Resco de Dios et al. 2015)
    double vp = meteoland::utils_averageDailyVP(tmin, tmax, rhmin, rhmax);
    double D = std::max(0.0, meteoland::utils_saturationVP(tmax) - vp);
    fm_dead = 5.43 + 52.91*exp(-0.64*D); 
  }
  double wind = meteovec["wind"];
  if(!NumericVector::is_na(fireHazardStandardWind)) wind = fireHazardStandardWind;
  
  //Calculate cohort canopy moisture to the average of canopy live and dead fuels, considering that a fraction of LAI is dead
  //proportionally to stem PLC (Ruffault et al. 2023)
  NumericVector LFMC = Plants["LFMC"];
  NumericVector PLC = Plants["StemPLC"];
  NumericVector canopyFMC = (LFMC*(1.0 - PLC) + fm_dead*PLC);
  
  NumericVector cohHeight = above["H"];
  NumericVector cohCR = above["CR"];
  NumericVector cohLoading = above["Loading"];
  //Correct loading for phenology
  NumericVector LAI_expanded = above["LAI_expanded"];
  NumericVector LAI_live = above["LAI_live"];
  for(int i=0;i<LAI_live.size();i++){
    if(LAI_live[i]>0.0) cohLoading[i] = cohLoading[i]*(LAI_expanded[i]/LAI_live[i]);
    else cohLoading[i] = 0.0;
  }
  
  //Average canopy moisture in the crown and surface layers
  NumericVector ActFMC = FCCSprops["ActFMC"];
  NumericVector w = FCCSprops["w"];
  if(w[0] > 0.0) ActFMC[0] = layerFuelAverageParameter(200.0, 10000.0, canopyFMC, cohLoading, cohHeight, cohCR);
  if(w[1] > 0.0) ActFMC[1] = layerFuelAverageParameter(0.0, 200.0, canopyFMC, cohLoading, cohHeight, cohCR);

  NumericVector MdeadSI = NumericVector::create(fm_dead, fm_dead, fm_dead, fm_dead, fm_dead); 
  NumericVector MliveSI = NumericVector::create(90.0, 90.0, 60.0); //Default values (not actually used if ActFMC is non-missing)
  List fccs = FCCSbehaviour(FCCSprops, MliveSI, MdeadSI, slope, wind); 
  
  //Copy results to fireBehavior vector
  List surfaceFire = fccs["SurfaceFire"];
  List crownFire = fccs["CrownFire"];
  List firePotentials = fccs["FirePotentials"];
  fireHazard["Loading_overstory [kg/m2]"] = w[0];
  fireHazard["Loading_understory [kg/m2]"] = w[1];
  fireHazard["CFMC_overstory [%]"] = ActFMC[0];
  fireHazard["CFMC_understory [%]"] = ActFMC[1];
  fireHazard["DFMC [%]"] = fm_dead;
  fireHazard["ROS_surface [m/min]"] = surfaceFire["ROS [m/min]"];
  fireHazard["I_b_surface [kW/m]"] = surfaceFire["I_b [kW/m]"];
  fireHazard["t_r_surface [s]"] = surfaceFire["t_r [s]"];
  fireHazard["FL_surface [m]"] = surfaceFire["FL [m]"];
  fireHazard["Ic_ratio"] = crownFire["Ic_ratio"];
  fireHazard["ROS_crown [m/min]"] = crownFire["ROS_crown [m/min]"];
  fireHazard["I_b_crown [kW/m]"] = crownFire["I_b_crown [kW/m]"];
  fireHazard["t_r_crown [s]"] = crownFire["t_r_crown [s]"];
  fireHazard["FL_crown [m]"] = crownFire["FL_crown [m]"];
  fireHazard["SFP"] = firePotentials["SFP"];
  fireHazard["CFP"] = firePotentials["CFP"];
}


// Soil water balance with simple hydraulic model
void spwbDay_basic(List internalCommunication, List x, NumericVector meteovec, 
                   double elevation, double slope, double aspect,
                   double runon = 0.0, Nullable<NumericVector> lateralFlows = R_NilValue, double waterTableDepth = NA_REAL, 
                   bool verbose = false) {
  
  //Retrieve communication structures
  List modelOutputComm = internalCommunication["basicSPWBOutput"];
  List transpOutput = internalCommunication["basicTranspirationOutput"];
  List SWBcommunication = internalCommunication["SWBcommunication"];
  
  NumericVector weather = modelOutputComm["weather"];
  weather["tday"] = meteovec["tday"];
  weather["prec"] = meteovec["prec"];
  weather["tmin"] = meteovec["tmin"];
  weather["tmax"] = meteovec["tmax"];
  weather["rhmin"] = meteovec["rhmin"];
  weather["rhmax"] = meteovec["rhmax"];
  weather["rad"] = meteovec["rad"];
  weather["wind"] = meteovec["wind"];
  weather["Catm"] = meteovec["Catm"];
  weather["Patm"] = meteovec["Patm"];
  weather["pet"] = meteovec["pet"];
  weather["rint"] = meteovec["rint"];
  NumericVector topo = modelOutputComm["topography"];
  topo["elevation"] = elevation;
  topo["slope"] = slope;
  topo["aspect"] = aspect;
  
  //Control parameters
  List control = x["control"];
  bool bareSoilEvaporation = control["bareSoilEvaporation"];
  String rhizosphereOverlap = control["rhizosphereOverlap"];
  bool plantWaterPools = (rhizosphereOverlap!="total");
  String soilFunctions = control["soilFunctions"];
  String infiltrationMode = control["infiltrationMode"];
  double infiltrationCorrection = control["infiltrationCorrection"];
  String soilDomains = control["soilDomains"];
  int ndailysteps = control["ndailysteps"];
  int max_nsubsteps_soil = control["max_nsubsteps_soil"];
  
  //Soil parameters
  DataFrame soil = Rcpp::as<Rcpp::DataFrame>(x["soil"]);
  DataFrame belowdf = Rcpp::as<Rcpp::DataFrame>(x["below"]);
  int nlayers = soil.nrow();
  NumericVector Wsoil = soil["W"];
  NumericVector Tsoil = soil["Temp"];
  
  List belowLayers = x["belowLayers"];
  NumericMatrix Wpool = Rcpp::as<Rcpp::NumericMatrix>(belowLayers["Wpool"]);
  
  //Weather input
  double tday = meteovec["tday"];
  double pet = meteovec["pet"]; 
  double prec  = meteovec["prec"];
  double rainfallIntensity  = meteovec["rint"]; 
  double rad = NA_REAL; 
  if(meteovec.containsElementNamed("rad")) rad = meteovec["rad"];
  
  //Set soil temperature to tday
  for(int l=0; l<nlayers; l++) Tsoil[l] = tday;
  
  //Vegetation input
  DataFrame cohorts = Rcpp::as<Rcpp::DataFrame>(x["cohorts"]);
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  NumericVector LAIlive = Rcpp::as<Rcpp::NumericVector>(above["LAI_live"]);
  NumericVector LAIphe = Rcpp::as<Rcpp::NumericVector>(above["LAI_expanded"]);
  NumericVector LAIdead = Rcpp::as<Rcpp::NumericVector>(above["LAI_dead"]);
  NumericVector H = Rcpp::as<Rcpp::NumericVector>(above["H"]);
  NumericVector CR = Rcpp::as<Rcpp::NumericVector>(above["CR"]);
  int numCohorts = LAIphe.size();
  
  
  //Parameters  
  DataFrame paramsInterception = Rcpp::as<Rcpp::DataFrame>(x["paramsInterception"]);
  NumericVector kPAR = Rcpp::as<Rcpp::NumericVector>(paramsInterception["kPAR"]);
  NumericVector gRainIntercept = Rcpp::as<Rcpp::NumericVector>(paramsInterception["g"]);
  
  DataFrame paramsPhenology = Rcpp::as<Rcpp::DataFrame>(x["paramsPhenology"]);
  NumericVector Sgdd = Rcpp::as<Rcpp::NumericVector>(paramsPhenology["Sgdd"]);
  
  
  //Copy clone soil and copy from Wpool to soil pools
  List soilPools(numCohorts);
  if(plantWaterPools) {
    for(int c=0;c<numCohorts;c++) {
      //Clone soil and copy moisture values from x
      List soil_c =  clone(soil);
      NumericVector W_c = soil_c["W"];
      for(int l=0;l<nlayers;l++) {
        W_c[l] = Wpool(c,l); 
      }
      soilPools[c] = soil_c;
    }
  }
  
  //STEP 1 - Update leaf area values according to the phenology of species and recalculate radiation extinction 
  double s = 0.0, LAIcell = 0.0, LAIcelllive = 0.0, LAIcellexpanded = 0.0, LAIcelldead = 0.0, Cm = 0.0;
  for(int c=0;c<numCohorts;c++) {
    s += (kPAR[c]*(LAIphe[c]+LAIdead[c]));
    LAIcell += LAIphe[c]+LAIdead[c];
    LAIcelldead += LAIdead[c];
    LAIcelllive += LAIlive[c];
    LAIcellexpanded +=LAIphe[c];
    Cm += (LAIphe[c]+LAIdead[c])*gRainIntercept[c]; //LAI dead also counts on interception
  }
  //Percentage of irradiance reaching the herb layer
  double LherbSWR = 100.0*exp((-1.0)*s/1.35);
  //Herb layer effects on light extinction and interception
  double herbLAI = x["herbLAI"];
  s += 0.5*herbLAI;
  Cm += herbLAI*1.0;
  LAIcell += herbLAI;
  //Percentage of irradiance reaching the ground
  double LgroundPAR = 100.0*exp((-1.0)*s);
  double LgroundSWR = 100.0*exp((-1.0)*s/1.35);
  
  //STEP 2 - Hidrological inputs (modifies snowpack)
  NumericVector hydroInputs = waterInputs(x, 
                                          prec, rainfallIntensity, 
                                          pet, tday, rad, elevation,
                                          Cm, LgroundPAR, LgroundSWR, 
                                          true);
  double RainfallInput = hydroInputs["NetRain"];
  double Snowmelt = hydroInputs["Snowmelt"];
  
  //STEP 3 - Evaporation from bare soil and herbaceous transpiration
  double snowpack = x["snowpack"];
  NumericVector EherbVec(nlayers,0.0);
  double Esoil = 0.0;
  NumericVector EsoilPools(numCohorts, 0.0);
  NumericMatrix EherbPools(numCohorts, nlayers);
  if(!plantWaterPools) {
    //Evaporation from bare soil if there is no snow (do not yet modify soil)
    if(bareSoilEvaporation) Esoil = soilEvaporation(soil, snowpack, soilFunctions, pet, LgroundSWR, false);
    //Herbaceous transpiration (do not yet modify soil)
    EherbVec = herbaceousTranspiration(pet, LherbSWR, herbLAI, soil, soilFunctions, false);
  } else {
    NumericVector poolProportions = belowdf["poolProportions"];
    for(int c=0;c<numCohorts;c++) {
      //Get soil pool
      List soil_c =  soilPools[c];
      //Evaporation from bare soil_c (if there is no snow), do not modify soil
      if(bareSoilEvaporation) {
        EsoilPools[c] = soilEvaporation(soil_c, snowpack, soilFunctions, pet, LgroundSWR, false);
        Esoil = Esoil + poolProportions[c]*EsoilPools[c]; 
      }
      //Herbaceous transpiration, do not modify soil
      NumericVector EherbVec_c = herbaceousTranspiration(pet, LherbSWR, herbLAI, soil_c, soilFunctions, false);
      //Update average soil evaporation and herbaceous transpiration 
      for(int l=0;l<nlayers;l++) {
        EherbPools(c,l) = EherbVec_c[l];
        EherbVec[l] = EherbVec[l] + poolProportions[c]*EherbVec_c[l]; 
      }
    }
  }
  
  //STEP 4 - Woody plant transpiration  (does not modify soil, only plants)
  transpirationBasic(transpOutput, x, meteovec, elevation, true);
  //Determine hydraulic redistribution and source sink for overall soil
  NumericMatrix soilLayerExtract = Rcpp::as<Rcpp::NumericMatrix>(transpOutput["Extraction"]);
  NumericVector ExtractionVec(nlayers, 0.0);
  NumericVector soilHydraulicInput(nlayers, 0.0); //Water that entered into the layer across all time steps
  NumericVector soilHydraulicOutput(nlayers, 0.0);  //Water that left the layer across all time steps
  for(int l=0;l<nlayers;l++) {
    for(int c=0;c<numCohorts;c++) {
      soilHydraulicInput[l] += (-1.0)*std::min(soilLayerExtract(c,l),0.0);
      soilHydraulicOutput[l] += std::max(soilLayerExtract(c,l),0.0);
      ExtractionVec[l] += soilLayerExtract(c,l);
    }
  }
  
  //STEP 5 - Soil flows
  double DeepDrainage = 0.0;
  double Infiltration = 0.0;
  double InfiltrationExcess = 0.0;
  double SaturationExcess = 0.0;
  double Runoff = 0.0;
  double CapillarityRise = 0.0;
  NumericVector sourceSinkVec(nlayers, 0.0);
  for(int l=0;l<nlayers;l++) {
    sourceSinkVec[l] -= (ExtractionVec[l] + EherbVec[l]);
    if(l ==0) sourceSinkVec[l] -= Esoil;
  }
  
  if(!plantWaterPools) {
    // determine water flows (no mass conservation)
    NumericVector sf = soilWaterBalance_inner(SWBcommunication, soil, soilFunctions,
                                        RainfallInput, rainfallIntensity, Snowmelt, sourceSinkVec, 
                                        runon, lateralFlows, waterTableDepth,
                                        infiltrationMode, infiltrationCorrection, soilDomains, 
                                        ndailysteps, max_nsubsteps_soil, true);
    DeepDrainage = sf["DeepDrainage"];
    Infiltration = sf["Infiltration"];
    Runoff = sf["Runoff"];
    InfiltrationExcess = sf["InfiltrationExcess"];
    SaturationExcess = sf["SaturationExcess"];
    CapillarityRise = sf["CapillarityRise"];
  } else { //Apply soil flows to water pools
    NumericVector poolProportions = belowdf["poolProportions"];
    List ExtractionPools = Rcpp::as<Rcpp::List>(transpOutput["ExtractionPools"]);
    // NumericVector sourceSinkCheck(nlayers, 0.0);
    //Set Wsoil to zero
    for(int l=0;l<nlayers;l++) Wsoil[l] = 0.0;
    NumericMatrix ExtractionPoolMat(numCohorts, nlayers);
    ExtractionPoolMat.fill(0.0);
    for(int c=0;c<numCohorts;c++) {
      //this is used to store extraction of a SINGLE plant cohort from all pools
      NumericMatrix ExtractionPoolsCoh = Rcpp::as<Rcpp::NumericMatrix>(ExtractionPools[c]);
      for(int l=0;l<nlayers;l++) {
        for(int c2=0;c2<numCohorts;c2++) {
          ExtractionPoolMat(c2,l) += ExtractionPoolsCoh(c2,l)/poolProportions[c2];
        }
      }
    }
    for(int c=0;c<numCohorts;c++) {
      List soil_c = soilPools[c];
      NumericVector sourceSinkPoolVec(nlayers, 0.0);
      for(int l=0;l<nlayers;l++) {
        sourceSinkPoolVec[l] -= (ExtractionPoolMat(c,l) + EherbPools(c,l));
        if(l ==0) sourceSinkPoolVec[l] -= EsoilPools[c];
      }
      NumericVector sf_c = soilWaterBalance_inner(SWBcommunication, soil_c, soilFunctions,
                                            RainfallInput, rainfallIntensity, Snowmelt, sourceSinkPoolVec, 
                                            runon, lateralFlows, waterTableDepth,
                                            infiltrationMode, infiltrationCorrection, soilDomains, 
                                            ndailysteps, max_nsubsteps_soil, true);
      double DeepDrainage_c = sf_c["DeepDrainage"];
      double Infiltration_c = sf_c["Infiltration"];
      double InfiltrationExcess_c = sf_c["InfiltrationExcess"];
      double Runoff_c = sf_c["Runoff"];
      double SaturationExcess_c = sf_c["SaturationExcess"];
      double CapillarityRise_c = sf_c["CapillarityRise"];
      DeepDrainage += DeepDrainage_c*poolProportions[c]; 
      Runoff += Runoff_c*poolProportions[c]; 
      Infiltration += Infiltration_c*poolProportions[c]; 
      SaturationExcess += SaturationExcess_c*poolProportions[c]; 
      InfiltrationExcess += InfiltrationExcess_c*poolProportions[c];
      CapillarityRise += CapillarityRise_c*poolProportions[c];
      
      //copy to Wpool and update Wsoil
      NumericVector W_c = soil_c["W"];
      for(int l=0;l<nlayers;l++) {
        Wpool(c,l) = W_c[l];
        Wsoil[l] = Wsoil[l] + W_c[l]*poolProportions[c];
      }
    }
    // for(int l=0; l<nlayers;l++) Rcout<< sourceSinkCheck[l] << " " << sourceSinkVec[l]<<"\n";
  }
  
  //Calculate current soil water potential for output
  NumericVector psiVec = psi(soil, soilFunctions); 
  
  //STEP 6 - Fire hazard
  bool fireHazardResults = control["fireHazardResults"];
  if(fireHazardResults) {
    NumericVector fireHazard = modelOutputComm["FireHazard"];
    fccsHazard(fireHazard, x, meteovec, transpOutput, slope);
  } 
  
  // Arrange output
  NumericVector WaterBalance = modelOutputComm["WaterBalance"];
  WaterBalance["PET"] = pet;
  WaterBalance["Rain"] = hydroInputs["Rain"];
  WaterBalance["Snow"] = hydroInputs["Snow"]; 
  WaterBalance["NetRain"] = hydroInputs["NetRain"];
  WaterBalance["Snowmelt"] = Snowmelt;
  WaterBalance["Runon"] = runon; 
  WaterBalance["Infiltration"] = Infiltration; 
  WaterBalance["InfiltrationExcess"] = InfiltrationExcess;
  WaterBalance["SaturationExcess"] = SaturationExcess;
  WaterBalance["Runoff"] = Runoff; 
  WaterBalance["DeepDrainage"] = DeepDrainage;
  WaterBalance["CapillarityRise"] = CapillarityRise;
  WaterBalance["SoilEvaporation"] = Esoil;
  WaterBalance["HerbTranspiration"] = sum(EherbVec);
  WaterBalance["PlantExtraction"] = sum(ExtractionVec);
  DataFrame outputPlants = Rcpp::as<Rcpp::DataFrame>(transpOutput["Plants"]);
  NumericVector Eplant = Rcpp::as<Rcpp::NumericVector>(outputPlants["Transpiration"]);
  double Transpiration = 0.0;
  for(int c=0;c<numCohorts;c++) Transpiration += Eplant[c];
  WaterBalance["Transpiration"] = Transpiration;
  WaterBalance["HydraulicRedistribution"] = sum(soilHydraulicInput);
  
  NumericVector Stand = modelOutputComm["Stand"];
  Stand["LAI"] = LAIcell;
  Stand["LAIherb"] = herbLAI; 
  Stand["LAIlive"] = LAIcelllive;
  Stand["LAIexpanded"] = LAIcellexpanded;
  Stand["LAIdead"] = LAIcelldead;
  Stand["Cm"] = Cm; 
  Stand["LgroundPAR"] = LgroundPAR; 
  Stand["LgroundSWR"] = LgroundSWR;
  
  
  DataFrame Soil = as<DataFrame>(modelOutputComm["Soil"]);
  NumericVector Psi = Soil["Psi"];
  NumericVector HerbTranspiration = Soil["HerbTranspiration"];
  NumericVector HydraulicInput = Soil["HydraulicInput"];
  NumericVector HydraulicOutput = Soil["HydraulicOutput"];
  NumericVector PlantExtraction = Soil["PlantExtraction"];
  for(int l=0;l<nlayers;l++) {
    Psi[l] = psiVec[l];
    HerbTranspiration[l] = EherbVec[l];
    HydraulicInput[l] = soilHydraulicInput[l];
    HydraulicOutput[l] = soilHydraulicOutput[l];
    PlantExtraction[l] = ExtractionVec[l];
  }
}

// Soil water balance with Sperry or Sureau hydraulic and stomatal conductance models
void spwbDay_advanced(List internalCommunication, List x, NumericVector meteovec, 
                      double latitude, double elevation, double slope, double aspect,
                      double solarConstant, double delta, 
                      double runon = 0.0, Nullable<NumericVector> lateralFlows = R_NilValue, double waterTableDepth = NA_REAL, 
                      bool verbose = false) {
  
  
  //Retrieve communication structures
  List modelOutputComm = internalCommunication["advancedSPWBOutput"];
  List transpOutput = internalCommunication["advancedTranspirationOutput"];
  List SWBcommunication = internalCommunication["SWBcommunication"];
  List SEBcommunication = internalCommunication["SEBcommunication"];
  
  NumericVector weather = modelOutputComm["weather"];
  weather["prec"] = meteovec["prec"];
  weather["tmin"] = meteovec["tmin"];
  weather["tmax"] = meteovec["tmax"];
  weather["tminPrev"] = meteovec["tminPrev"];
  weather["tmaxPrev"] = meteovec["tmaxPrev"];
  weather["tminNext"] = meteovec["tminNext"];
  weather["rhmin"] = meteovec["rhmin"];
  weather["rhmax"] = meteovec["rhmax"];
  weather["rad"] = meteovec["rad"];
  weather["wind"] = meteovec["wind"];
  weather["Catm"] = meteovec["Catm"];
  weather["Patm"] = meteovec["Patm"];
  weather["pet"] = meteovec["pet"];
  weather["rint"] = meteovec["rint"];
  NumericVector topo = modelOutputComm["topography"];
  topo["elevation"] = elevation;
  topo["slope"] = slope;
  topo["aspect"] = aspect;
  
  //Control parameters
  List control = x["control"];
  int ntimesteps = control["ndailysteps"];
  bool bareSoilEvaporation = control["bareSoilEvaporation"];
  String rhizosphereOverlap = control["rhizosphereOverlap"];
  bool plantWaterPools = (rhizosphereOverlap!="total");
  String soilFunctions = control["soilFunctions"];
  String infiltrationMode = control["infiltrationMode"];
  String soilDomains = control["soilDomains"];
  double infiltrationCorrection = control["infiltrationCorrection"];
  int ndailysteps = control["ndailysteps"];
  int max_nsubsteps_soil = control["max_nsubsteps_soil"];
  
  //Soil parameters
  DataFrame soil = Rcpp::as<Rcpp::DataFrame>(x["soil"]);
  DataFrame belowdf = Rcpp::as<Rcpp::DataFrame>(x["below"]);
  int nlayers = Rcpp::as<Rcpp::NumericVector>(soil["widths"]).size();
  
  List belowLayers = x["belowLayers"];
  NumericMatrix Wpool = belowLayers["Wpool"];
  NumericVector Wsoil = soil["W"];
  
  //Meteo input
  double tmin = meteovec["tmin"];
  double tmax = meteovec["tmax"];
  double prec = meteovec["prec"];
  double rad = meteovec["rad"];
  double pet = meteovec["pet"];
  double rainfallIntensity = meteovec["rint"];
  
  //Vegetation input
  DataFrame cohorts = Rcpp::as<Rcpp::DataFrame>(x["cohorts"]);
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  NumericVector LAIlive = Rcpp::as<Rcpp::NumericVector>(above["LAI_live"]);
  NumericVector LAIphe = Rcpp::as<Rcpp::NumericVector>(above["LAI_expanded"]);
  NumericVector LAIdead = Rcpp::as<Rcpp::NumericVector>(above["LAI_dead"]);
  NumericVector H = Rcpp::as<Rcpp::NumericVector>(above["H"]);
  NumericVector CR = Rcpp::as<Rcpp::NumericVector>(above["CR"]);
  int numCohorts = LAIlive.size();
  
  //Base parameters
  DataFrame paramsPhenology = Rcpp::as<Rcpp::DataFrame>(x["paramsPhenology"]);
  NumericVector Sgdd = Rcpp::as<Rcpp::NumericVector>(paramsPhenology["Sgdd"]);
  
  DataFrame paramsInterception = Rcpp::as<Rcpp::DataFrame>(x["paramsInterception"]);
  NumericVector kPAR = Rcpp::as<Rcpp::NumericVector>(paramsInterception["kPAR"]);
  NumericVector gRainIntercept = Rcpp::as<Rcpp::NumericVector>(paramsInterception["g"]);
  
  //Copy clone soil and copy from Wpool to soil pools
  List soilPools(numCohorts);
  if(plantWaterPools) {
    for(int c=0;c<numCohorts;c++) {
      //Clone soil and copy moisture values from x
      List soil_c =  clone(soil);
      NumericVector W_c = soil_c["W"];
      for(int l=0;l<nlayers;l++) W_c[l] = Wpool(c,l);
      soilPools[c] = soil_c;
    }
  }
  
  //STEP 1 - Leaf Phenology: Adjusted leaf area index
  double tday = meteoland::utils_averageDaylightTemperature(tmin, tmax);
  double s = 0.0, LAIcell = 0.0, LAIcelldead = 0.0, LAIcelllive = 0.0,  LAIcellexpanded = 0.0, Cm = 0.0, LAIcellmax = 0.0;
  for(int c=0;c<numCohorts;c++) {
    LAIcell += (LAIphe[c]+LAIdead[c]);
    LAIcelldead += LAIdead[c];
    LAIcellmax += LAIlive[c];
    LAIcelllive += LAIlive[c];
    LAIcellexpanded +=LAIphe[c];
    s += (kPAR[c]*(LAIphe[c]+LAIdead[c]));
    Cm += (LAIphe[c]+LAIdead[c])*gRainIntercept[c]; //LAI dead also counts on interception
  }
  //Percentage of irradiance reaching the herb layer
  double LherbSWR = 100.0*exp((-1.0)*s/1.35);
  //Herb layer effects on light extinction and interception
  double herbLAI = x["herbLAI"];
  s += 0.5*herbLAI;
  Cm += herbLAI*1.0;
  LAIcell +=herbLAI;
  double LgroundPAR = 100.0*exp((-1.0)*s);
  double LgroundSWR = 100.0*exp((-1.0)*s/1.35);
  
  //STEP 2 - Interception, snow pack dynamics and soil water input (modifies snowpack)
  NumericVector hydroInputs = waterInputs(x,
                                          prec, rainfallIntensity, 
                                          pet, tday, rad, elevation,
                                          Cm, LgroundPAR, LgroundSWR, 
                                          true);
  double RainfallInput = hydroInputs["NetRain"];
  double Snowmelt = hydroInputs["Snowmelt"];
  
  //STEP 3 - Evaporation from bare soil and herbaceous transpiration
  double snowpack = x["snowpack"];
  NumericVector EherbVec(nlayers,0.0);
  double Esoil = 0.0;
  NumericVector EsoilPools(numCohorts, 0.0);
  NumericMatrix EherbPools(numCohorts, nlayers);
  if(!plantWaterPools) {
    //Evaporation from bare soil if there is no snow (do not yet modify soil)
    if(bareSoilEvaporation) Esoil = soilEvaporation(soil, snowpack, 
       soilFunctions, pet, LgroundSWR, false);
    //Herbaceous transpiration (do not yet modify soil)
    EherbVec = herbaceousTranspiration(pet, LherbSWR, herbLAI, soil, soilFunctions, false);
  } else {
    NumericVector poolProportions = belowdf["poolProportions"];
    for(int c=0;c<numCohorts;c++) {
      //Get soil pool
      List soil_c =  soilPools[c];
      //Evaporation from bare soil_c (if there is no snow), do not modify soil
      if(bareSoilEvaporation) {
        EsoilPools[c] = soilEvaporation(soil_c, snowpack, soilFunctions, pet, LgroundSWR, false);
        Esoil = Esoil + poolProportions[c]*EsoilPools[c]; 
      }
      //Herbaceous transpiration, do not modify soil
      NumericVector EherbVec_c = herbaceousTranspiration(pet, LherbSWR, herbLAI, soil_c, soilFunctions, false);
      for(int l = 0;l<nlayers;l++) EherbPools(c,l) = EherbVec_c[l];
      //Update average soil evaporation and herbaceous transpiration 
      for(int l=0;l<nlayers;l++) {
        EherbVec[l] = EherbVec[l] + poolProportions[c]*EherbPools(c,l); 
      }
    }
  }
  
  //STEPS 4-8 - Energy balance, transpiration, photosynthesis, uptake 
  transpirationAdvanced(SEBcommunication, transpOutput, x, meteovec, 
                        latitude, elevation, slope, aspect, 
                        solarConstant, delta, 
                        hydroInputs["Interception"], hydroInputs["Snowmelt"], Esoil, sum(EherbVec),
                        verbose, NA_INTEGER, true);
  NumericMatrix soilLayerExtractInst = Rcpp::as<Rcpp::NumericMatrix>(transpOutput["ExtractionInst"]);
  
  NumericVector ExtractionVec(nlayers, 0.0);
  NumericVector soilHydraulicInput(nlayers, 0.0); //Water that entered into the layer across all time steps
  NumericVector soilHydraulicOutput(nlayers, 0.0);  //Water that left the layer across all time steps
  for(int l=0;l<nlayers;l++) {
    for(int n=0;n<ntimesteps;n++) {
      soilHydraulicInput[l] += (-1.0)*std::min(soilLayerExtractInst(l,n),0.0);
      soilHydraulicOutput[l] += std::max(soilLayerExtractInst(l,n),0.0);
    }
    ExtractionVec[l] = sum(soilLayerExtractInst(l,_));
  }
  
  //STEP 9 - SOIL FLOWS
  double DeepDrainage = 0.0;
  double Infiltration = 0.0;
  double Runoff = 0.0;
  double InfiltrationExcess = 0.0;
  double SaturationExcess = 0.0;
  double CapillarityRise = 0.0;
  NumericVector sourceSinkVec(nlayers, 0.0);
  for(int l=0;l<nlayers;l++) {
    sourceSinkVec[l] -= (ExtractionVec[l] + EherbVec[l]);
    if(l ==0) sourceSinkVec[l] -= Esoil;
  }
  if(!plantWaterPools) {
    // determine water flows (no mass conservation)
    NumericVector sf = soilWaterBalance_inner(SWBcommunication, soil, soilFunctions,
                                        RainfallInput, rainfallIntensity, Snowmelt, sourceSinkVec, 
                                        runon, lateralFlows, waterTableDepth,
                                        infiltrationMode, infiltrationCorrection, soilDomains, 
                                        ndailysteps, max_nsubsteps_soil, true);
    DeepDrainage = sf["DeepDrainage"];
    Infiltration = sf["Infiltration"];
    Runoff = sf["Runoff"];
    InfiltrationExcess = sf["InfiltrationExcess"];
    SaturationExcess = sf["SaturationExcess"];
    CapillarityRise = sf["CapillarityRise"]; 
  } else {
    NumericVector poolProportions = belowdf["poolProportions"];
    List ExtractionPools = Rcpp::as<Rcpp::List>(transpOutput["ExtractionPools"]);
    //Set Wsoil to zero
    for(int l=0;l<nlayers;l++) Wsoil[l] = 0.0;
    NumericMatrix ExtractionPoolMat(numCohorts, nlayers);
    ExtractionPoolMat.fill(0.0);
    for(int c=0;c<numCohorts;c++) {
      //this is used to store extraction of a SINGLE plant cohort from all pools
      NumericMatrix ExtractionPoolsCoh = Rcpp::as<Rcpp::NumericMatrix>(ExtractionPools[c]);
      for(int l=0;l<nlayers;l++) {
        for(int c2=0;c2<numCohorts;c2++) {
          ExtractionPoolMat(c2,l) += ExtractionPoolsCoh(c2,l)/poolProportions[c2];
        }
      }
    }
    for(int c=0;c<numCohorts;c++) {
      List soil_c = soilPools[c];
      NumericVector sourceSinkPoolVec(nlayers, 0.0);
      for(int l=0;l<nlayers;l++) {
        sourceSinkPoolVec[l] -= (ExtractionPoolMat(c,l) + EherbPools(c,l));
        if(l ==0) sourceSinkPoolVec[l] -= EsoilPools[c];
      }
      NumericVector sf_c = soilWaterBalance_inner(SWBcommunication, soil_c, soilFunctions,
                                            RainfallInput, rainfallIntensity, Snowmelt, sourceSinkPoolVec, 
                                            runon, lateralFlows, waterTableDepth,
                                            infiltrationMode, infiltrationCorrection, soilDomains, 
                                            ndailysteps, max_nsubsteps_soil, true);
      double DeepDrainage_c = sf_c["DeepDrainage"];
      double Infiltration_c = sf_c["Infiltration"];
      double Runoff_c = sf_c["Runoff"];
      double InfiltrationExcess_c = sf_c["InfiltrationExcess"];
      double SaturationExcess_c = sf_c["SaturationExcess"];
      double CapillarityRise_c = sf_c["CapillarityRise"]; 
      DeepDrainage += DeepDrainage_c*poolProportions[c]; 
      Runoff += Runoff_c*poolProportions[c]; 
      Infiltration += Infiltration_c*poolProportions[c]; 
      InfiltrationExcess += InfiltrationExcess_c*poolProportions[c]; 
      SaturationExcess += SaturationExcess_c*poolProportions[c]; 
      CapillarityRise += CapillarityRise_c*poolProportions[c]; 
      
      //copy to Wpool and update Wsoil
      NumericVector W_c = soil_c["W"];
      for(int l=0;l<nlayers;l++) {
        Wpool(c,l) = W_c[l];
        Wsoil[l] = Wsoil[l] + W_c[l]*poolProportions[c];
      }
    }
  }
  //Calculate current soil water potential for output
  NumericVector psiVec = psi(soil, soilFunctions); 
  
  //STEP 11 - Fire hazard
  bool fireHazardResults = control["fireHazardResults"];
  if(fireHazardResults) {
    NumericVector fireHazard = modelOutputComm["FireHazard"];
    fccsHazard(fireHazard, x, meteovec, transpOutput, slope);
  }
  
  // Arrange output
  NumericVector WaterBalance = modelOutputComm["WaterBalance"];
  WaterBalance["PET"] = pet;
  WaterBalance["Rain"] = hydroInputs["Rain"];
  WaterBalance["Snow"] = hydroInputs["Snow"]; 
  WaterBalance["NetRain"] = hydroInputs["NetRain"];
  WaterBalance["Snowmelt"] = Snowmelt;
  WaterBalance["Runon"] = runon; 
  WaterBalance["Infiltration"] = Infiltration; 
  WaterBalance["InfiltrationExcess"] = InfiltrationExcess;
  WaterBalance["SaturationExcess"] = SaturationExcess;
  WaterBalance["Runoff"] = Runoff; 
  WaterBalance["DeepDrainage"] = DeepDrainage;
  WaterBalance["CapillarityRise"] = CapillarityRise;
  WaterBalance["SoilEvaporation"] = Esoil;
  WaterBalance["HerbTranspiration"] = sum(EherbVec);
  WaterBalance["PlantExtraction"] = sum(ExtractionVec);
  DataFrame outputPlants = Rcpp::as<Rcpp::DataFrame>(transpOutput["Plants"]);
  NumericVector Eplant = outputPlants["Transpiration"];
  WaterBalance["Transpiration"] = sum(Eplant);
  WaterBalance["HydraulicRedistribution"] = sum(soilHydraulicInput);
  
  NumericVector Stand = modelOutputComm["Stand"];
  Stand["LAI"] = LAIcell;
  Stand["LAIherb"] = herbLAI; 
  Stand["LAIlive"] = LAIcelllive;
  Stand["LAIexpanded"] = LAIcellexpanded;
  Stand["LAIdead"] = LAIcelldead;
  Stand["Cm"] = Cm; 
  Stand["LgroundPAR"] = LgroundPAR; 
  Stand["LgroundSWR"] = LgroundSWR;
  
  DataFrame Soil = as<DataFrame>(modelOutputComm["Soil"]);
  NumericVector Psi = Soil["Psi"];
  NumericVector HerbTranspiration = Soil["HerbTranspiration"];
  NumericVector HydraulicInput = Soil["HydraulicInput"];
  NumericVector HydraulicOutput = Soil["HydraulicOutput"];
  NumericVector PlantExtraction = Soil["PlantExtraction"];
  for(int l=0;l<nlayers;l++) {
    Psi[l] = psiVec[l];
    HerbTranspiration[l] = EherbVec[l];
    HydraulicInput[l] = soilHydraulicInput[l];
    HydraulicOutput[l] = soilHydraulicOutput[l];
    PlantExtraction[l] = ExtractionVec[l];
  }
}



//' @rdname communication
//' @keywords internal
// [[Rcpp::export("spwb_day_inner")]]
void spwbDay_inner(List internalCommunication, List x, CharacterVector date, NumericVector meteovec, 
                   double latitude, double elevation, double slope = NA_REAL, double aspect = NA_REAL,  
                   double runon = 0.0, Nullable<NumericVector> lateralFlows = R_NilValue, double waterTableDepth = NA_REAL,
                   bool modifyInput = true) {
  double tmin = meteovec["MinTemperature"];
  double tmax = meteovec["MaxTemperature"];
  if(tmin > tmax) {
    warning("tmin > tmax. Swapping values.");
    double swap = tmin;
    tmin = tmax;
    tmax = swap;
  }
  double rhmin = meteovec["MinRelativeHumidity"];
  double rhmax = meteovec["MaxRelativeHumidity"];
  if(NumericVector::is_na(rhmax)) {
    warning("Maximum relative humidity assumed 100");
    rhmax = 100.0;
  }
  if(NumericVector::is_na(rhmin)) {
    warning("Minimum relative humidity estimated from temperature range");
    double vp_tmin = meteoland::utils_saturationVP(tmin);
    double vp_tmax = meteoland::utils_saturationVP(tmax);
    rhmin = std::min(rhmax, 100.0*(vp_tmin/vp_tmax));
  }
  if(rhmin > rhmax) {
    warning("rhmin > rhmax. Swapping values.");
    double swap = rhmin;
    rhmin = rhmax;
    rhmax = swap;
  }
  double rad = meteovec["Radiation"];
  double prec = meteovec["Precipitation"];
  double wind = NA_REAL;
  if(meteovec.containsElementNamed("WindSpeed")) wind = meteovec["WindSpeed"];
  double Catm = NA_REAL; 
  if(meteovec.containsElementNamed("CO2")) Catm = meteovec["CO2"];
  double Patm = NA_REAL; 
  if(meteovec.containsElementNamed("Patm")) Patm = meteovec["Patm"];
  double Rint = NA_REAL; 
  if(meteovec.containsElementNamed("RainfallIntensity")) Rint = meteovec["RainfallIntensity"];
  //Control parameters
  List control = x["control"];
  bool verbose = control["verbose"];
  
  bool leafPhenology = control["leafPhenology"];
  String transpirationMode = control["transpirationMode"];
  if(NumericVector::is_na(Catm)) Catm = control["defaultCO2"];
  
  //Will not modify input x 
  if(!modifyInput) {
    x = clone(x);
  }
  
  //Soul parameters
  DataFrame soil = Rcpp::as<Rcpp::DataFrame>(x["soil"]);
  
  std::string c = as<std::string>(date[0]);
  int month = std::atoi(c.substr(5,2).c_str());
  int J = meteoland::radiation_julianDay(std::atoi(c.substr(0, 4).c_str()),std::atoi(c.substr(5,2).c_str()),std::atoi(c.substr(8,2).c_str()));
  double delta = meteoland::radiation_solarDeclination(J);
  double solarConstant = meteoland::radiation_solarConstant(J);
  double latrad = latitude * (M_PI/180.0);
  if(NumericVector::is_na(aspect)) aspect = 0.0;
  if(NumericVector::is_na(slope)) slope = 0.0;
  double asprad = aspect * (M_PI/180.0);
  double slorad = slope * (M_PI/180.0);
  double photoperiod = meteoland::radiation_daylength(latrad, 0.0, 0.0, delta);
  double tday = meteoland::utils_averageDaylightTemperature(tmin, tmax);
  if(NumericVector::is_na(rad)) {
    warning("Estimating solar radiation");
    double vpa = meteoland::utils_averageDailyVP(tmin, tmax, rhmin, rhmax);
    rad = meteoland::radiation_solarRadiation(solarConstant, latrad, elevation,
                                              slorad, asprad, delta, tmax -tmin, tmax-tmin,
                                              vpa, prec);
  }
  double pet = meteoland::penman(latrad, elevation, slorad, asprad, J, tmin, tmax, rhmin, rhmax, rad, wind);
  
  //Derive doy from date  
  int J0101 = meteoland::radiation_julianDay(std::atoi(c.substr(0, 4).c_str()),1,1);
  int doy = J - J0101+1;
  
  if(NumericVector::is_na(wind)) wind = control["defaultWindSpeed"]; 
  if(wind<0.1) wind = 0.1; //Minimum windspeed abovecanopy
  
  NumericVector defaultRainfallIntensityPerMonth = control["defaultRainfallIntensityPerMonth"];
  if(NumericVector::is_na(Rint)) Rint = rainfallIntensity(month, prec, defaultRainfallIntensityPerMonth);
  
  //Update phenology
  if(leafPhenology) {
    updatePhenology(x, doy, photoperiod, tday);
    updateLeaves(x, wind, false);
  }
  
  List modelOutput;
  if(transpirationMode=="Granier") {
    NumericVector meteovec_bas = NumericVector::create(
      Named("tday") = tday, 
      Named("prec") = prec,
      Named("tmin") = tmin, 
      Named("tmax") = tmax,
      Named("rhmin") = rhmin, 
      Named("rhmax") = rhmax, 
      Named("rad") = rad, 
      Named("wind") = wind, 
      Named("Catm") = Catm,
      Named("Patm") = Patm,
      Named("pet") = pet,
      Named("rint") = Rint);
    spwbDay_basic(internalCommunication, x, meteovec_bas,
                  elevation, slope, aspect, 
                  runon, lateralFlows, waterTableDepth, 
                  verbose);
  } else {
    NumericVector meteovec_adv = NumericVector::create(
      Named("tmin") = tmin, 
      Named("tmax") = tmax,
      Named("tminPrev") = tmin, 
      Named("tmaxPrev") = tmax, 
      Named("tminNext") = tmin, 
      Named("prec") = prec,
      Named("rhmin") = rhmin, 
      Named("rhmax") = rhmax, 
      Named("rad") = rad, 
      Named("wind") = wind, 
      Named("Catm") = Catm,
      Named("Patm") = Patm,
      Named("pet") = pet,
      Named("rint") = Rint);
    spwbDay_advanced(internalCommunication, x, meteovec_adv,
                     latitude, elevation, slope, aspect,
                     solarConstant, delta, 
                     runon, lateralFlows, waterTableDepth, 
                     verbose);
  }
}

//' Single-day simulation
//'
//' Function \code{spwb_day} performs water balance for a single day and \code{growth_day} 
//' performs water and carbon balance for a single day.
//' 
//' @param x An object of class \code{\link{spwbInput}} or \code{\link{growthInput}}.
//' @param date Date as string "yyyy-mm-dd".
//' @param meteovec A named numerical vector with weather data. See variable names in parameter \code{meteo} of \code{\link{spwb}}.
//' @param latitude Latitude (in degrees).
//' @param elevation,slope,aspect Elevation above sea level (in m), slope (in degrees) and aspect (in degrees from North). 
//' @param runon Surface water amount running on the target area from upslope (in mm).
//' @param lateralFlows Lateral source/sink terms for each soil layer (interflow/to from adjacent locations) as mm/day.
//' @param waterTableDepth Water table depth (in mm). When not missing, capillarity rise will be allowed if lower than total soil depth.
//' @param modifyInput Boolean flag to indicate that the input \code{x} object is allowed to be modified during the simulation.
//' 
//' @details
//' The simulation functions allow using three different sub-models of transpiration and photosynthesis:
//' \itemize{
//'   \item{The sub-model corresponding to 'Granier' transpiration mode is illustrated by function \code{\link{transp_transpirationGranier}} and was described in De Caceres et al. (2015),
//'   and implements an approach originally described in Granier et al. (1999).} 
//'   \item{The sub-model corresponding to 'Sperry' transpiration mode is illustrated by function \code{\link{transp_transpirationSperry}} and was described in De Caceres et al. (2021), and
//'   implements a modelling approach originally described in Sperry et al. (2017).}  
//'   \item{The sub-model corresponding to 'Sureau' transpiration mode is illustrated by function \code{\link{transp_transpirationSureau}} and was described for model SurEau-Ecos v2.0 in Ruffault et al. (2022).} 
//' }
//' 
//' Simulations using the 'Sperry' or 'Sureau' transpiration mode are computationally much more expensive than 'Granier'.
//' 
//' @return
//' Function \code{spwb_day()} returns a list of class \code{spwb_day} with the 
//' following elements:
//' \itemize{
//'   \item{\code{"cohorts"}: A data frame with cohort information, copied from \code{\link{spwbInput}}.}
//'   \item{\code{"topography"}: Vector with elevation, slope and aspect given as input.} 
//'   \item{\code{"weather"}: A vector with the input weather.}
//'   \item{\code{"WaterBalance"}: A vector of water balance components (rain, snow, net rain, infiltration, ...) for the simulated day, equivalent to one row of 'WaterBalance' object given in \code{\link{spwb}}.}
//'   \item{\code{"Soil"}: A data frame with results for each soil layer:
//'     \itemize{
//'       \item{\code{"Psi"}: Soil water potential (in MPa) at the end of the day.}
//'       \item{\code{"HerbTranspiration"}: Water extracted by herbaceous plants from each soil layer (in mm).}
//'       \item{\code{"HydraulicInput"}: Water entering each soil layer from other layers, transported via plant roots (in mm).}
//'       \item{\code{"HydraulicOutput"}: Water leaving each soil layer (going to other layers or the transpiration stream) (in mm).}
//'       \item{\code{"PlantExtraction"}: Water extracted by woody plants from each soil layer (in mm).}
//'     }
//'   }
//'   \item{\code{"Stand"}: A named vector with with stand values for the simulated day, equivalent to one row of 'Stand' object returned by \code{\link{spwb}}.}
//'   \item{\code{"Plants"}: A data frame of results for each plant cohort (see \code{\link{transp_transpirationGranier}} or \code{\link{transp_transpirationSperry}}).}
//' }
//' The following items are only returned when \code{transpirationMode = "Sperry"} or  \code{transpirationMode = "Sureau"}:
//' \itemize{
//'   \item{\code{"EnergyBalance"}: Energy balance of the stand (see \code{\link{transp_transpirationSperry}}).}
//'   \item{\code{"RhizoPsi"}: Minimum water potential (in MPa) inside roots, after crossing rhizosphere, per cohort and soil layer.}
//'   \item{\code{"SunlitLeaves"} and \code{"ShadeLeaves"}: For each leaf type, a data frame with values of LAI, Vmax298 and Jmax298 for leaves of this type in each plant cohort.}
//'   \item{\code{"ExtractionInst"}: Water extracted by each plant cohort during each time step.}
//'   \item{\code{"PlantsInst"}: A list with instantaneous (per time step) results for each plant cohort (see \code{\link{transp_transpirationSperry}}).}
//'   \item{\code{"LightExtinction"}: A list of information regarding radiation balance through the canopy, as returned by function \code{\link{light_instantaneousLightExtinctionAbsortion}}.}
//'   \item{\code{"CanopyTurbulence"}: Canopy turbulence (see \code{\link{wind_canopyTurbulence}}).}
//' }
//'   
//' @references
//' De \enc{Cáceres}{Caceres} M, \enc{Martínez}{Martinez}-Vilalta J, Coll L, Llorens P, Casals P, Poyatos R, Pausas JG, Brotons L. (2015) Coupling a water balance model with forest inventory data to predict drought stress: the role of forest structural changes vs. climate changes. Agricultural and Forest Meteorology 213: 77-90 (doi:10.1016/j.agrformet.2015.06.012).
//' 
//' De \enc{Cáceres}{Caceres} M, Mencuccini M, Martin-StPaul N, Limousin JM, Coll L, Poyatos R, Cabon A, Granda V, Forner A, Valladares F, \enc{Martínez}{Martinez}-Vilalta J (2021) Unravelling the effect of species mixing on water use and drought stress in holm oak forests: a modelling approach. Agricultural and Forest Meteorology 296 (doi:10.1016/j.agrformet.2020.108233).
//' 
//' Granier A, \enc{Bréda}{Breda} N, Biron P, Villette S (1999) A lumped water balance model to evaluate duration and intensity of drought constraints in forest stands. Ecol Modell 116:269–283. https://doi.org/10.1016/S0304-3800(98)00205-1.
//' 
//' Ruffault J, Pimont F, Cochard H, Dupuy JL, Martin-StPaul N (2022) 
//' SurEau-Ecos v2.0: a trait-based plant hydraulics model for simulations of plant water status and drought-induced mortality at the ecosystem level.
//' Geoscientific Model Development 15, 5593-5626 (doi:10.5194/gmd-15-5593-2022).
//' 
//' Sperry, J. S., M. D. Venturas, W. R. L. Anderegg, M. Mencuccini, D. S. Mackay, Y. Wang, and D. M. Love. 2017. Predicting stomatal responses to the environment from the optimization of photosynthetic gain and hydraulic cost. Plant Cell and Environment 40, 816-830 (doi: 10.1111/pce.12852).
//' 
//' @author
//' \itemize{
//'   \item{Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF}
//'   \item{Nicolas Martin-StPaul, URFM-INRAE}
//' }
//' 
//' @seealso
//' \code{\link{spwbInput}}, \code{\link{spwb}},  \code{\link{plot.spwb_day}},  
//' \code{\link{growthInput}}, \code{\link{growth}},  \code{\link{plot.growth_day}}  
//' 
//' @examples
//' #Load example daily meteorological data
//' data(examplemeteo)
//' 
//' #Load example plot plant data
//' data(exampleforest)
//' 
//' #Default species parameterization
//' data(SpParamsMED)
//' 
//' #Define soil parameters
//' examplesoil <- defaultSoilParams(4)
//' 
//' # Day to be simulated
//' d <- 100
//' meteovec <- unlist(examplemeteo[d,-1])
//' date <- as.character(examplemeteo$dates[d])
//' 
//' #Simulate water balance one day only (Granier mode)
//' control <- defaultControl("Granier")
//' x1 <- spwbInput(exampleforest,examplesoil, SpParamsMED, control)
//' sd1 <- spwb_day(x1, date, meteovec,  
//'                 latitude = 41.82592, elevation = 100, slope=0, aspect=0) 
//' 
//' #Simulate water balance for one day only (Sperry mode)
//' control <- defaultControl("Sperry")
//' x2 <- spwbInput(exampleforest, examplesoil, SpParamsMED, control)
//' sd2 <-spwb_day(x2, date, meteovec,
//'               latitude = 41.82592, elevation = 100, slope=0, aspect=0)
//' 
//' #Simulate water balance for one day only (Sureau mode)
//' control <- defaultControl("Sureau")
//' x3 <- spwbInput(exampleforest, examplesoil, SpParamsMED, control)
//' sd3 <-spwb_day(x3, date, meteovec,
//'               latitude = 41.82592, elevation = 100, slope=0, aspect=0)
//' 
//' 
//' #Simulate water and carbon balance for one day only (Granier mode)
//' control <- defaultControl("Granier")
//' x4  <- growthInput(exampleforest,examplesoil, SpParamsMED, control)
//' sd4 <- growth_day(x4, date, meteovec,
//'                 latitude = 41.82592, elevation = 100, slope=0, aspect=0)
//' 
//' #Simulate water and carbon balance for one day only (Sperry mode)
//' control <- defaultControl("Sperry")
//' x5  <- growthInput(exampleforest,examplesoil, SpParamsMED, control)
//' sd5 <- growth_day(x5, date, meteovec,
//'                 latitude = 41.82592, elevation = 100, slope=0, aspect=0)
//' 
//' #Simulate water and carbon balance for one day only (Sureau mode)
//' control <- defaultControl("Sureau")
//' x6  <- growthInput(exampleforest,examplesoil, SpParamsMED, control)
//' sd6 <- growth_day(x6, date, meteovec,
//'                 latitude = 41.82592, elevation = 100, slope=0, aspect=0)
//' 
//' @name spwb_day
// [[Rcpp::export("spwb_day")]]
List spwbDay(List x, CharacterVector date, NumericVector meteovec, 
              double latitude, double elevation, double slope = NA_REAL, double aspect = NA_REAL,  
              double runon = 0.0, Nullable<NumericVector> lateralFlows = R_NilValue, double waterTableDepth = NA_REAL,
              bool modifyInput = true) {
   
   //Instance communication structures
   List internalCommunication = instanceCommunicationStructures(x, "spwb");
   
   spwbDay_inner(internalCommunication, x, date, meteovec,
                 latitude, elevation, slope, aspect,
                 runon, lateralFlows, waterTableDepth,
                 modifyInput);
   
   List modelOutput = copyModelOutput(internalCommunication, x, "spwb");
   return(modelOutput);
 }


