#include <RcppArmadillo.h>
#include <numeric>
#include <math.h>
#include "biophysicsutils_c.h"
#include "communication_structures_c.h"
#include "lightextinction_basic_c.h"
#include "windextinction_c.h"
#include "firebehaviour_c.h"
#include "hydraulics_c.h"
#include "hydrology_c.h"
#include "forestutils_c.h"
#include "modelInput_c.h"
#include "phenology.h"
#include "transpiration.h"
#include "fuelstructure.h"
#include "firebehaviour.h"
#include "spwb_day_basic_c.h"
#include "tissuemoisture.h"
#include "soil.h"
#include <meteoland.h>

// Soil water balance with simple hydraulic model
void spwbDay_basic_c(BasicSPWB_RESULT& BSPWBres, BasicSPWB_COMM& BSPWB_comm, ModelInput& x, 
                     const WeatherInputVector& meteovec, 
                     const double elevation, const double slope, const double aspect,
                     const double runon, 
                     const std::vector<double>& lateralFlows, const double waterTableDepth) {
  
  //Retrieve communication structures
  SoilWaterBalance_COMM& SWBcomm = BSPWB_comm.SWBcomm;
  BasicTranspiration_COMM& BTcomm = BSPWB_comm.BTcomm;
   
   
   //Meteo input
   double pet = meteovec.pet;
   double rhmax = meteovec.rhmax;
   double rhmin = meteovec.rhmin;
   double tday = meteovec.tday;
   double prec = meteovec.prec;
   double tmax = meteovec.tmax;
   double tmin = meteovec.tmin;
   double Catm = meteovec.Catm;
   double Patm = meteovec.Patm;
   double rad = meteovec.rad;
   double wind = meteovec.wind;
   double rainfallIntensity = meteovec.rint;

   //Store topography for output
   BSPWBres.topo.elevation = elevation;
   BSPWBres.topo.slope = slope;
   BSPWBres.topo.aspect = aspect;
   
  //Control parameters
  bool bareSoilEvaporation = x.control.commonWB.bareSoilEvaporation;
  std::string& rhizosphereOverlap = x.control.rhizosphereOverlap;
  bool plantWaterPools = (rhizosphereOverlap!="total");
  std::string& infiltrationMode = x.control.commonWB.infiltrationMode;
  double infiltrationCorrection = x.control.commonWB.infiltrationCorrection;
  std::string& soilDomains = x.control.soilDomains;
  int ndailysteps = x.control.advancedWB.ndailysteps; // MAYBE CHANGE TO COMMONWB parameters
  int max_nsubsteps_soil = x.control.commonWB.max_nsubsteps_soil;

  //Soil
  Soil& soil = x.soil;
  int nlayers = soil.getNlayers();
  
  //Set soil temperature to tday
  for(int l=0; l<nlayers; l++) soil.setTemp(l, tday);
  
  //Water pools
  arma::mat& Wpool = x.belowLayers.Wpool;
  
  // Canopy
  int ncanlayers = x.canopy.zlow.size();
  
  //Vegetation input
  std::vector<double>& LAIlive = x.above.LAI_live;
  std::vector<double>& LAIphe = x.above.LAI_expanded;
  std::vector<double>& LAIdead = x.above.LAI_dead;
  std::vector<double>& H = x.above.H;
  std::vector<double>& CR = x.above.CR;
  int numCohorts = LAIphe.size();
  
  
  // //Soil parameters
  // DataFrame belowdf = Rcpp::as<Rcpp::DataFrame>(x["below"]);
  // 
  // List belowLayers = x["belowLayers"];
  // 
  // 
  // 
  // //Vegetation input
  // DataFrame cohorts = Rcpp::as<Rcpp::DataFrame>(x["cohorts"]);
  // DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  // 
  // 
  
  //Parameters
  const std::vector<double>& kPAR = x.paramsInterception.kPAR;
  const std::vector<double>& gRainIntercept = x.paramsInterception.g;
  
  const std::vector<double>& Sgdd = x.paramsPhenology.Sgdd;
  
  // //Copy clone soil and copy from Wpool to soil pools
  // List soilPools(numCohorts);
  // if(plantWaterPools) {
  //   for(int c=0;c<numCohorts;c++) {
  //     //Clone soil and copy moisture values from x
  //     List soil_c =  clone(soil);
  //     NumericVector W_c = soil_c["W"];
  //     for(int l=0;l<nlayers;l++) {
  //       W_c[l] = Wpool(c,l); 
  //     }
  //     soilPools[c] = soil_c;
  //   }
  // }
  

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
  s += 0.5*x.herbLAI;
  Cm += x.herbLAI*1.0;
  LAIcell += x.herbLAI;
  
  //Percentage of irradiance reaching the ground
  double LgroundPAR = 100.0*exp((-1.0)*s);
  double LgroundSWR = 100.0*exp((-1.0)*s/1.35);

  //STEP 2 - Hidrological inputs (modifies snowpack)
  waterInputs_c(BSPWB_comm.waterInputs, x,
                prec, rainfallIntensity,
                pet, tday, rad, elevation,
                Cm, LgroundPAR, LgroundSWR,
                true);
  double RainfallInput = BSPWB_comm.waterInputs.netrain;
  double Snowmelt = BSPWB_comm.waterInputs.melt;

  //STEP 3 - Evaporation from bare soil and herbaceous transpiration
  double snowpack = x.snowpack;
  std::vector<double> EherbVec(nlayers,0.0);
  double Esoil = 0.0;
  // NumericVector EsoilPools(numCohorts, 0.0);
  // NumericMatrix EherbPools(numCohorts, nlayers);
  if(!plantWaterPools) {
    //Evaporation from bare soil if there is no snow (do not yet modify soil)
    if(bareSoilEvaporation) Esoil = soilEvaporation_c(soil, snowpack, pet, LgroundSWR, false);
    //Herbaceous transpiration (do not yet modify soil)
    // herbaceousTranspiration_c(EherbVec, soil, pet, LherbSWR, herbLAI,  false);
  } else {
  //   NumericVector poolProportions = belowdf["poolProportions"];
  //   for(int c=0;c<numCohorts;c++) {
  //     //Get soil pool
  //     List soil_c =  soilPools[c];
  //     //Evaporation from bare soil_c (if there is no snow), do not modify soil
  //     if(bareSoilEvaporation) {
  //       EsoilPools[c] = soilEvaporation(soil_c, snowpack, soilFunctions, pet, LgroundSWR, false);
  //       Esoil = Esoil + poolProportions[c]*EsoilPools[c]; 
  //     }
  //     //Herbaceous transpiration, do not modify soil
  //     NumericVector EherbVec_c = herbaceousTranspiration(pet, LherbSWR, herbLAI, soil_c, soilFunctions, false);
  //     //Update average soil evaporation and herbaceous transpiration 
  //     for(int l=0;l<nlayers;l++) {
  //       EherbPools(c,l) = EherbVec_c[l];
  //       EherbVec[l] = EherbVec[l] + poolProportions[c]*EherbVec_c[l]; 
  //     }
  //   }
  }
  // 
  // //STEP 4 - Woody plant transpiration  (does not modify soil, only plants)
  // transpirationBasic(transpOutput, x, meteovec, elevation, true);
  // //Determine hydraulic redistribution and source sink for overall soil
  // NumericMatrix soilLayerExtract = Rcpp::as<Rcpp::NumericMatrix>(transpOutput["Extraction"]);
  // NumericVector ExtractionVec(nlayers, 0.0);
  // NumericVector soilHydraulicInput(nlayers, 0.0); //Water that entered into the layer across all time steps
  // NumericVector soilHydraulicOutput(nlayers, 0.0);  //Water that left the layer across all time steps
  // for(int l=0;l<nlayers;l++) {
  //   for(int c=0;c<numCohorts;c++) {
  //     soilHydraulicInput[l] += (-1.0)*std::min(soilLayerExtract(c,l),0.0);
  //     soilHydraulicOutput[l] += std::max(soilLayerExtract(c,l),0.0);
  //     ExtractionVec[l] += soilLayerExtract(c,l);
  //   }
  // }
  // 
  // //STEP 5 - Soil flows
  // double DeepDrainage = 0.0;
  // double Infiltration = 0.0;
  // double InfiltrationExcess = 0.0;
  // double SaturationExcess = 0.0;
  // double Runoff = 0.0;
  // double CapillarityRise = 0.0;
  // NumericVector sourceSinkVec(nlayers, 0.0);
  // for(int l=0;l<nlayers;l++) {
  //   sourceSinkVec[l] -= (ExtractionVec[l] + EherbVec[l]);
  //   if(l ==0) sourceSinkVec[l] -= Esoil;
  // }
  // 
  // if(soilDomains != "none") {
  //   if(!plantWaterPools) {
  //     // determine water flows (no mass conservation)
  //     NumericVector sf = soilWaterBalance_inner(SWBcommunication, soil, soilFunctions,
  //                                               RainfallInput, rainfallIntensity, Snowmelt, sourceSinkVec, 
  //                                               runon, lateralFlows, waterTableDepth,
  //                                               infiltrationMode, infiltrationCorrection, soilDomains, 
  //                                               ndailysteps, max_nsubsteps_soil, true);
  //     DeepDrainage = sf["DeepDrainage"];
  //     Infiltration = sf["Infiltration"];
  //     Runoff = sf["Runoff"];
  //     InfiltrationExcess = sf["InfiltrationExcess"];
  //     SaturationExcess = sf["SaturationExcess"];
  //     CapillarityRise = sf["CapillarityRise"];
  //   } else { //Apply soil flows to water pools
  //     NumericVector poolProportions = belowdf["poolProportions"];
  //     List ExtractionPools = Rcpp::as<Rcpp::List>(transpOutput["ExtractionPools"]);
  //     // NumericVector sourceSinkCheck(nlayers, 0.0);
  //     //Set Wsoil to zero
  //     for(int l=0;l<nlayers;l++) Wsoil[l] = 0.0;
  //     NumericMatrix ExtractionPoolMat(numCohorts, nlayers);
  //     ExtractionPoolMat.fill(0.0);
  //     for(int c=0;c<numCohorts;c++) {
  //       //this is used to store extraction of a SINGLE plant cohort from all pools
  //       NumericMatrix ExtractionPoolsCoh = Rcpp::as<Rcpp::NumericMatrix>(ExtractionPools[c]);
  //       for(int l=0;l<nlayers;l++) {
  //         for(int c2=0;c2<numCohorts;c2++) {
  //           ExtractionPoolMat(c2,l) += ExtractionPoolsCoh(c2,l)/poolProportions[c2];
  //         }
  //       }
  //     }
  //     for(int c=0;c<numCohorts;c++) {
  //       List soil_c = soilPools[c];
  //       NumericVector sourceSinkPoolVec(nlayers, 0.0);
  //       for(int l=0;l<nlayers;l++) {
  //         sourceSinkPoolVec[l] -= (ExtractionPoolMat(c,l) + EherbPools(c,l));
  //         if(l ==0) sourceSinkPoolVec[l] -= EsoilPools[c];
  //       }
  //       NumericVector sf_c = soilWaterBalance_inner(SWBcommunication, soil_c, soilFunctions,
  //                                                   RainfallInput, rainfallIntensity, Snowmelt, sourceSinkPoolVec, 
  //                                                   runon, lateralFlows, waterTableDepth,
  //                                                   infiltrationMode, infiltrationCorrection, soilDomains, 
  //                                                   ndailysteps, max_nsubsteps_soil, true);
  //       double DeepDrainage_c = sf_c["DeepDrainage"];
  //       double Infiltration_c = sf_c["Infiltration"];
  //       double InfiltrationExcess_c = sf_c["InfiltrationExcess"];
  //       double Runoff_c = sf_c["Runoff"];
  //       double SaturationExcess_c = sf_c["SaturationExcess"];
  //       double CapillarityRise_c = sf_c["CapillarityRise"];
  //       DeepDrainage += DeepDrainage_c*poolProportions[c]; 
  //       Runoff += Runoff_c*poolProportions[c]; 
  //       Infiltration += Infiltration_c*poolProportions[c]; 
  //       SaturationExcess += SaturationExcess_c*poolProportions[c]; 
  //       InfiltrationExcess += InfiltrationExcess_c*poolProportions[c];
  //       CapillarityRise += CapillarityRise_c*poolProportions[c];
  //       
  //       //copy to Wpool and update Wsoil
  //       NumericVector W_c = soil_c["W"];
  //       for(int l=0;l<nlayers;l++) {
  //         Wpool(c,l) = W_c[l];
  //         Wsoil[l] = Wsoil[l] + W_c[l]*poolProportions[c];
  //       }
  //     }
  //     // for(int l=0; l<nlayers;l++) Rcout<< sourceSinkCheck[l] << " " << sourceSinkVec[l]<<"\n";
  //   }
  // }
  // 
  // //Calculate current soil water potential for output
  // NumericVector psiVec = psi(soil, soilFunctions); 
  // 
  // //STEP 6 - Fire hazard
  // bool fireHazardResults = control["fireHazardResults"];
  // if(fireHazardResults) {
  //   NumericVector fireHazard = modelOutputComm["FireHazard"];
  //   fccsHazard(fireHazard, x, meteovec, transpOutput, slope);
  // } 
  // 
  // // Arrange output
  // NumericVector WaterBalance = modelOutputComm["WaterBalance"];
  // WaterBalance["PET"] = pet;
  // WaterBalance["Rain"] = hydroInputs["Rain"];
  // WaterBalance["Snow"] = hydroInputs["Snow"]; 
  // WaterBalance["NetRain"] = hydroInputs["NetRain"];
  // WaterBalance["Snowmelt"] = Snowmelt;
  // WaterBalance["Runon"] = runon; 
  // WaterBalance["Infiltration"] = Infiltration; 
  // WaterBalance["InfiltrationExcess"] = InfiltrationExcess;
  // WaterBalance["SaturationExcess"] = SaturationExcess;
  // WaterBalance["Runoff"] = Runoff; 
  // WaterBalance["DeepDrainage"] = DeepDrainage;
  // WaterBalance["CapillarityRise"] = CapillarityRise;
  // WaterBalance["SoilEvaporation"] = Esoil;
  // WaterBalance["HerbTranspiration"] = sum(EherbVec);
  // WaterBalance["PlantExtraction"] = sum(ExtractionVec);
  // DataFrame outputPlants = Rcpp::as<Rcpp::DataFrame>(transpOutput["Plants"]);
  // NumericVector Eplant = Rcpp::as<Rcpp::NumericVector>(outputPlants["Transpiration"]);
  // double Transpiration = 0.0;
  // for(int c=0;c<numCohorts;c++) Transpiration += Eplant[c];
  // WaterBalance["Transpiration"] = Transpiration;
  // WaterBalance["HydraulicRedistribution"] = sum(soilHydraulicInput);
  // 
  // NumericVector Stand = modelOutputComm["Stand"];
  // Stand["LAI"] = LAIcell;
  // Stand["LAIherb"] = herbLAI; 
  // Stand["LAIlive"] = LAIcelllive;
  // Stand["LAIexpanded"] = LAIcellexpanded;
  // Stand["LAIdead"] = LAIcelldead;
  // Stand["Cm"] = Cm; 
  // Stand["LgroundPAR"] = LgroundPAR; 
  // Stand["LgroundSWR"] = LgroundSWR;
  // 
  // 
  // DataFrame Soil = as<DataFrame>(modelOutputComm["Soil"]);
  // NumericVector Psi = Soil["Psi"];
  // NumericVector HerbTranspiration = Soil["HerbTranspiration"];
  // NumericVector HydraulicInput = Soil["HydraulicInput"];
  // NumericVector HydraulicOutput = Soil["HydraulicOutput"];
  // NumericVector PlantExtraction = Soil["PlantExtraction"];
  // for(int l=0;l<nlayers;l++) {
  //   Psi[l] = psiVec[l];
  //   HerbTranspiration[l] = EherbVec[l];
  //   HydraulicInput[l] = soilHydraulicInput[l];
  //   HydraulicOutput[l] = soilHydraulicOutput[l];
  //   PlantExtraction[l] = ExtractionVec[l];
  // }
}

Rcpp::List copyBasicSPWBResult_c(const BasicSPWB_RESULT& BSPWBres, ModelInput& x) {
  Rcpp::List l = Rcpp::List::create(_["cohorts"] = copyCohorts_c(x.cohorts),
                                    _["topography"] = copyTopo_c(BSPWBres.topo),
                                    _["weather"] = copyWeather_c(BSPWBres.meteovec),
                                    _["WaterBalance"] = copyWaterBalanceResult_c(BSPWBres.WaterBalance));
  if(x.control.results.soilResults) {
    l.push_back(copySoilResult_c(BSPWBres.Soil), "Soil");
  }
  if(x.control.results.standResults) {
    l.push_back(copyStandResult_c(BSPWBres.Stand), "Stand");
  }
  if(x.control.results.plantResults) {
    l.push_back(copyPlantBasicTranspirationResult_c(BSPWBres.Plants, x), "Plants");
  }
  if(x.control.results.fireHazardResults) {
    l.push_back(copyFCCSResult_c(BSPWBres.fccs), "FireHazard");
  }
  l.attr("class") = CharacterVector::create("spwb_day","list");
  return(l);
}