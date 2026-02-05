#include <RcppArmadillo.h>
#include <numeric>
#include <math.h>
#include "biophysicsutils_c.h"
#include "lowlevel_structures_c.h"
#include "lightextinction_basic_c.h"
#include "windextinction_c.h"
#include "firebehaviour_c.h"
#include "hydraulics_c.h"
#include "hydrology_c.h"
#include "forestutils_c.h"
#include "modelInput_c.h"
#include "phenology_c.h"
#include "transpiration_basic_c.h"
#include "transpiration_advanced_c.h"
#include "fuelstructure.h"
#include "firebehaviour.h"
#include "spwb_day_c.h"
#include "tissuemoisture.h"
#include "soil.h"
#include "root_c.h"
#include <meteoland.h>


// Soil water balance with basic hydraulic model
void spwbDay_basic_c(BasicSPWB_RESULT& BSPWBres, BasicSPWB_COMM& BSPWB_comm, ModelInput& x, 
                     const WeatherInputVector& meteovec, 
                     const double elevation, const double slope, const double aspect,
                     const double runon, 
                     const std::vector<double>& lateralFlows, const double waterTableDepth) {
  
  //Retrieve communication structures
  BasicTranspiration_COMM& BTcomm = BSPWB_comm.BTcomm;
  BasicTranspiration_RESULT& BTres = BSPWBres.BTres;
  SoilWaterBalance_COMM& SWBcomm = BSPWB_comm.SWBcomm;
  SoilWaterBalance_RESULT& SWBres = BSPWBres.SWBres;
  
   
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
   
   //Store weather for output
   BSPWBres.meteovec.pet = pet;
   BSPWBres.meteovec.rhmax = rhmax;
   BSPWBres.meteovec.rhmin = rhmin;
   BSPWBres.meteovec.tmax = tmax;
   BSPWBres.meteovec.tmin = tmin;
   BSPWBres.meteovec.tday = tday;
   BSPWBres.meteovec.prec = prec;
   BSPWBres.meteovec.wind = wind;
   BSPWBres.meteovec.rad = rad;
   BSPWBres.meteovec.Catm = Catm;
   BSPWBres.meteovec.Patm = Patm;
   BSPWBres.meteovec.rint = rainfallIntensity;
   
   
  //Control parameters
  bool bareSoilEvaporation = x.control.commonWB.bareSoilEvaporation;
  std::string& rhizosphereOverlap = x.control.rhizosphereOverlap;
  bool plantWaterPools = (rhizosphereOverlap!="total");
  std::string& infiltrationMode = x.control.commonWB.infiltrationMode;
  double infiltrationCorrection = x.control.commonWB.infiltrationCorrection;
  std::string& soilDomains = x.control.soilDomains;
  int ndailysteps = x.control.advancedWB.ndailysteps; // MAYBE CHANGE TO COMMONWB parameters
  int max_nsubsteps_soil = x.control.commonWB.max_nsubsteps_soil;
  bool fireHazardResults = x.control.results.fireHazardResults;

  //Soil
  Soil& soil = x.soil;
  int nlayers = soil.getNlayers();
  
  
  //Set soil temperature to tday
  for(int l=0; l<nlayers; l++) soil.setTemp(l, tday);
  
  //Water pools
  arma::mat& Wpool = x.belowLayers.Wpool;
  std::vector<double>& poolProportions = x.below.poolProportions;
  
  //Vegetation input
  std::vector<double>& LAIlive = x.above.LAI_live;
  std::vector<double>& LAIphe = x.above.LAI_expanded;
  std::vector<double>& LAIdead = x.above.LAI_dead;
  int numCohorts = LAIphe.size();
  
  //Parameters
  const std::vector<double>& kPAR = x.paramsInterception.kPAR;
  const std::vector<double>& gRainIntercept = x.paramsInterception.g;
  
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
  
  double Esoil = 0.0;
  std::vector<double> EsoilPools(numCohorts,0.0);
  std::vector<double> EherbVec(nlayers,0.0);
  std::vector<double> EherbVec_c(nlayers,0.0);
  arma::mat EherbPools(numCohorts, nlayers);
  EherbPools.fill(0.0);
  std::vector<double> V(nlayers);
  ldrRS_one_c(V, 50, 500, NA_REAL, soil.getWidths());
  if(!plantWaterPools) {
    //Evaporation from bare soil if there is no snow (do not yet modify soil)
    if(bareSoilEvaporation) Esoil = soilEvaporation_c(soil, snowpack, pet, LgroundSWR, false);
    //Herbaceous transpiration (do not yet modify soil)
    herbaceousTranspiration_c(EherbVec, soil, pet, LherbSWR, x.herbLAI,  V, false);
  } else {
    //Store overall soil moisture in a backup copy
    double* Wbackup = new double[nlayers];
    for(int l = 0; l<nlayers;l++) Wbackup[l] = soil.getW(l);
    
    for(int c=0;c<numCohorts;c++) {
      for(int l = 0; l<nlayers;l++) soil.setW(l,Wpool(c,l)); // this updates psi, theta, ... 
      
      //Evaporation from bare  (if there is no snow), do not modify soil
      if(bareSoilEvaporation) {
        EsoilPools[c] = soilEvaporation_c(soil, snowpack, pet, LgroundSWR, false);
        Esoil = Esoil + poolProportions[c]*EsoilPools[c];
      }
      //Herbaceous transpiration, do not modify soil
      herbaceousTranspiration_c(EherbVec_c, soil, pet, LherbSWR, x.herbLAI,  V, false);
      //Update average soil evaporation and herbaceous transpiration
      for(int l=0;l<nlayers;l++) {
        EherbPools(c,l) = EherbVec_c[l];
        EherbVec[l] = EherbVec[l] + poolProportions[c]*EherbVec_c[l];
      }
    }
    //Restore soil moisture
    for(int l = 0; l<nlayers;l++) soil.setW(l, Wbackup[l]);
    delete[] Wbackup;
  }
  
  //STEP 4 - Woody plant transpiration  (does not modify soil, only plants)
  transpirationBasic_c(BTres, BTcomm, x, meteovec, elevation);

  //Determine hydraulic redistribution and source sink for overall soil
  for(int l=0;l<nlayers;l++) {
    BSPWBres.Soil.HerbTranspiration[l] = 0.0;
    BSPWBres.Soil.HydraulicInput[l] = 0.0;
    BSPWBres.Soil.HydraulicOutput[l] = 0.0;
    BSPWBres.Soil.PlantExtraction[l] = 0.0;
    for(int c=0;c<numCohorts;c++) {
      BSPWBres.Soil.HydraulicInput[l] += (-1.0)*std::min(BTres.extraction(c,l),0.0);
      BSPWBres.Soil.HydraulicOutput[l] += std::max(BTres.extraction(c,l),0.0);
      BSPWBres.Soil.PlantExtraction[l] += BTres.extraction(c,l);
    }
  }

  //STEP 5 - Soil flows
  double DeepDrainage = 0.0;
  double Infiltration = 0.0;
  double InfiltrationExcess = 0.0;
  double SaturationExcess = 0.0;
  double Runoff = 0.0;
  double CapillarityRise = 0.0;
  std::vector<double> sourceSinkVec(nlayers, 0.0);

  if(soilDomains != "none") {
    if(!plantWaterPools) {
      // determine water flows (no mass conservation)
      for(int l=0;l<nlayers;l++) {
        sourceSinkVec[l] -= (BSPWBres.Soil.PlantExtraction[l] + EherbVec[l]);
        if(l ==0) sourceSinkVec[l] -= Esoil;
      }
      soilWaterBalance_inner_c(SWBres, SWBcomm, soil,
                               RainfallInput, rainfallIntensity, Snowmelt, sourceSinkVec,
                               runon, lateralFlows, waterTableDepth,
                               infiltrationMode, infiltrationCorrection,
                               soilDomains,
                               ndailysteps, max_nsubsteps_soil);
      DeepDrainage = SWBres.deepDrainage_mm;
      Infiltration = SWBres.infiltration_mm;
      Runoff = SWBres.runoff_mm;
      InfiltrationExcess = SWBres.infiltrationExcess_mm;
      SaturationExcess = SWBres.saturationExcess_mm;
      CapillarityRise = SWBres.capillarityRise_mm;
    } else { //Apply soil flows to water pools
      arma::mat ExtractionPoolMat(numCohorts, nlayers);
      ExtractionPoolMat.fill(0.0);
      for(int c=0;c<numCohorts;c++) {
        //this is used to store extraction of a SINGLE plant cohort from all pools
        arma::mat& ExtractionPoolsCoh = BTres.extractionPools[c];
        for(int l=0;l<nlayers;l++) {
          for(int c2=0;c2<numCohorts;c2++) {
            ExtractionPoolMat(c2,l) += ExtractionPoolsCoh(c2,l)/poolProportions[c2];
          }
        }
      }
      //Create vector to store averaged soil moisture
      //Set Wsoil to zero
      double* Wsoil = new double[nlayers];
      for(int l = 0; l<nlayers;l++) Wsoil[l] = 0.0;
      for(int c=0;c<numCohorts;c++) {
        for(int l = 0; l<nlayers;l++) soil.setW(l,Wpool(c,l)); // this updates psi, theta, ... 
        for(int l=0;l<nlayers;l++) {
          sourceSinkVec[l] -= (ExtractionPoolMat(c,l) + EherbPools(c,l));
          if(l ==0) sourceSinkVec[l] -= EsoilPools[c];
        }
        soilWaterBalance_inner_c(SWBres, SWBcomm, soil,
                                 RainfallInput, rainfallIntensity, Snowmelt, sourceSinkVec,
                                 runon, lateralFlows, waterTableDepth,
                                 infiltrationMode, infiltrationCorrection,
                                 soilDomains,
                                 ndailysteps, max_nsubsteps_soil);

        DeepDrainage +=  SWBres.deepDrainage_mm*poolProportions[c];
        Runoff += SWBres.runoff_mm*poolProportions[c];
        Infiltration += SWBres.infiltration_mm*poolProportions[c];
        SaturationExcess += SWBres.saturationExcess_mm*poolProportions[c];
        InfiltrationExcess += SWBres.infiltrationExcess_mm*poolProportions[c];
        CapillarityRise += SWBres.capillarityRise_mm*poolProportions[c];

        //copy to Wpool and update averaged soil moisture (Wsoil)
        for(int l=0;l<nlayers;l++) {
          Wpool(c,l) = soil.getW(l);
          Wsoil[l] += soil.getW(l)*poolProportions[c];
        }
      }

      //Copy averaged soil moisture
      for(int l = 0; l<nlayers;l++) soil.setW(l, Wsoil[l]);
      delete[] Wsoil;
    }
  }
  //STEP 6 - Fire hazard
  if(fireHazardResults) {
    fccsHazard_c(BSPWBres.fccsbeh, BSPWBres.fccs, x, 
                 meteovec, 
                 BTres.plants.LFMC,
                 BTres.plants.StemPLC, 
                 slope);
  }

  // Arrange output
  BSPWBres.WaterBalance.PET = pet;
  BSPWBres.WaterBalance.Rain = BSPWB_comm.waterInputs.rain;
  BSPWBres.WaterBalance.Snow = BSPWB_comm.waterInputs.snow;
  BSPWBres.WaterBalance.NetRain = BSPWB_comm.waterInputs.netrain;
  BSPWBres.WaterBalance.Snowmelt = Snowmelt;
  BSPWBres.WaterBalance.Runon = runon;
  BSPWBres.WaterBalance.Infiltration = Infiltration;
  BSPWBres.WaterBalance.InfiltrationExcess = InfiltrationExcess;
  BSPWBres.WaterBalance.SaturationExcess = SaturationExcess;
  BSPWBres.WaterBalance.Runoff = Runoff;
  BSPWBres.WaterBalance.DeepDrainage = DeepDrainage;
  BSPWBres.WaterBalance.CapillarityRise = CapillarityRise;
  BSPWBres.WaterBalance.SoilEvaporation = Esoil;
  BSPWBres.WaterBalance.HerbTranspiration = std::accumulate(EherbVec.begin(), EherbVec.end(), 0.0);
  BSPWBres.WaterBalance.PlantExtraction = std::accumulate(BSPWBres.Soil.PlantExtraction.begin(), BSPWBres.Soil.PlantExtraction.end(), 0.0);
  BSPWBres.WaterBalance.Transpiration = std::accumulate(BTres.plants.Transpiration.begin(), BTres.plants.Transpiration.end(), 0.0);
  BSPWBres.WaterBalance.HydraulicRedistribution =std::accumulate(BSPWBres.Soil.HydraulicInput.begin(), BSPWBres.Soil.HydraulicInput.end(), 0.0);
  
  BSPWBres.Stand.LAI = LAIcell;
  BSPWBres.Stand.LAIherb = x.herbLAI; 
  BSPWBres.Stand.LAIlive = LAIcelllive;
  BSPWBres.Stand.LAIexpanded = LAIcellexpanded;
  BSPWBres.Stand.LAIdead = LAIcelldead;
  BSPWBres.Stand.Cm = Cm;
  BSPWBres.Stand.LgroundPAR = LgroundPAR;
  BSPWBres.Stand.LgroundSWR = LgroundSWR;
  
  //Copy final soil state to output
  for(int l=0;l<nlayers;l++) {
    BSPWBres.Soil.Psi[l] = soil.getPsi(l);
  }
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
    l.push_back(copyPlantBasicTranspirationResult_c(BSPWBres.BTres.plants, x), "Plants");
  }
  if(x.control.results.fireHazardResults) {
    l.push_back(copyFCCSResult_c(BSPWBres.fccs), "FireHazard");
  }
  l.attr("class") = CharacterVector::create("spwb_day","list");
  return(l);
}


// Soil water balance with advanced hydraulic model
void spwbDay_advanced_c(AdvancedSPWB_RESULT& ASPWBres, AdvancedSPWB_COMM& ASPWB_comm, ModelInput& x, 
                     const WeatherInputVector& meteovec, 
                     const double latitude, const double elevation, const double slope, const double aspect,
                     const double solarConstant, const double delta, 
                     const double runon, 
                     const std::vector<double>& lateralFlows, const double waterTableDepth) {
  
  //Retrieve communication structures
  AdvancedTranspiration_COMM& ATcomm = ASPWB_comm.ATcomm;
  AdvancedTranspiration_RESULT& ATres = ASPWBres.ATres;
  SoilWaterBalance_COMM& SWBcomm = ASPWB_comm.SWBcomm;
  SoilWaterBalance_RESULT& SWBres = ASPWBres.SWBres;
  
  
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
  ASPWBres.topo.elevation = elevation;
  ASPWBres.topo.slope = slope;
  ASPWBres.topo.aspect = aspect;
  
  //Store weather for output
  ASPWBres.meteovec.pet = pet;
  ASPWBres.meteovec.rhmax = rhmax;
  ASPWBres.meteovec.rhmin = rhmin;
  ASPWBres.meteovec.tmax = tmax;
  ASPWBres.meteovec.tmin = tmin;
  ASPWBres.meteovec.tday = tday;
  ASPWBres.meteovec.prec = prec;
  ASPWBres.meteovec.wind = wind;
  ASPWBres.meteovec.rad = rad;
  ASPWBres.meteovec.Catm = Catm;
  ASPWBres.meteovec.Patm = Patm;
  ASPWBres.meteovec.rint = rainfallIntensity;
  
  
  //Control parameters
  bool bareSoilEvaporation = x.control.commonWB.bareSoilEvaporation;
  std::string& rhizosphereOverlap = x.control.rhizosphereOverlap;
  bool plantWaterPools = (rhizosphereOverlap!="total");
  std::string& infiltrationMode = x.control.commonWB.infiltrationMode;
  double infiltrationCorrection = x.control.commonWB.infiltrationCorrection;
  std::string& soilDomains = x.control.soilDomains;
  int ndailysteps = x.control.advancedWB.ndailysteps; // MAYBE CHANGE TO COMMONWB parameters
  int max_nsubsteps_soil = x.control.commonWB.max_nsubsteps_soil;
  bool fireHazardResults = x.control.results.fireHazardResults;
  
  //Soil
  Soil& soil = x.soil;
  int nlayers = soil.getNlayers();
  
  
  //Set soil temperature to tday
  for(int l=0; l<nlayers; l++) soil.setTemp(l, tday);
  
  //Water pools
  arma::mat& Wpool = x.belowLayers.Wpool;
  std::vector<double>& poolProportions = x.below.poolProportions;
  
  //Vegetation input
  std::vector<double>& LAIlive = x.above.LAI_live;
  std::vector<double>& LAIphe = x.above.LAI_expanded;
  std::vector<double>& LAIdead = x.above.LAI_dead;
  int numCohorts = LAIphe.size();
  
  //Parameters
  const std::vector<double>& kPAR = x.paramsInterception.kPAR;
  const std::vector<double>& gRainIntercept = x.paramsInterception.g;
  
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
  
  
  //STEP 2 - Interception, snow pack dynamics and soil water input (modifies snowpack)
  waterInputs_c(ASPWB_comm.waterInputs, x,
                prec, rainfallIntensity,
                pet, tday, rad, elevation,
                Cm, LgroundPAR, LgroundSWR,
                true);
  double RainfallInput = ASPWB_comm.waterInputs.netrain;
  double Snowmelt = ASPWB_comm.waterInputs.melt;
  
  //STEP 3 - Evaporation from bare soil and herbaceous transpiration
  double snowpack = x.snowpack;
  
  double Esoil = 0.0;
  std::vector<double> EsoilPools(numCohorts,0.0);
  std::vector<double> EherbVec(nlayers,0.0);
  std::vector<double> EherbVec_c(nlayers,0.0);
  arma::mat EherbPools(numCohorts, nlayers);
  EherbPools.fill(0.0);
  std::vector<double> V(nlayers);
  ldrRS_one_c(V, 50, 500, NA_REAL, soil.getWidths());
  if(!plantWaterPools) {
    //Evaporation from bare soil if there is no snow (do not yet modify soil)
    if(bareSoilEvaporation) Esoil = soilEvaporation_c(soil, snowpack, pet, LgroundSWR, false);
    //Herbaceous transpiration (do not yet modify soil)
    herbaceousTranspiration_c(EherbVec, soil, pet, LherbSWR, x.herbLAI,  V, false);
  } else {
    //Store overall soil moisture in a backup copy
    double* Wbackup = new double[nlayers];
    for(int l = 0; l<nlayers;l++) Wbackup[l] = soil.getW(l);
    
    for(int c=0;c<numCohorts;c++) {
      for(int l = 0; l<nlayers;l++) soil.setW(l,Wpool(c,l)); // this updates psi, theta, ... 
      
      //Evaporation from bare  (if there is no snow), do not modify soil
      if(bareSoilEvaporation) {
        EsoilPools[c] = soilEvaporation_c(soil, snowpack, pet, LgroundSWR, false);
        Esoil = Esoil + poolProportions[c]*EsoilPools[c];
      }
      //Herbaceous transpiration, do not modify soil
      herbaceousTranspiration_c(EherbVec_c, soil, pet, LherbSWR, x.herbLAI,  V, false);
      //Update average soil evaporation and herbaceous transpiration
      for(int l=0;l<nlayers;l++) {
        EherbPools(c,l) = EherbVec_c[l];
        EherbVec[l] = EherbVec[l] + poolProportions[c]*EherbVec_c[l];
      }
    }
    //Restore soil moisture
    for(int l = 0; l<nlayers;l++) soil.setW(l, Wbackup[l]);
    delete[] Wbackup;
  }
  
  //STEPS 4-8 - Energy balance, transpiration, photosynthesis, uptake 
  transpirationAdvanced_c(ATres, ATcomm, x, meteovec, 
                          latitude, elevation, slope, aspect, 
                          solarConstant, delta, 
                          ASPWB_comm.waterInputs.interception, ASPWB_comm.waterInputs.melt, Esoil, 
                          std::accumulate(EherbVec.begin(), EherbVec.end(), 0.0),
                          -1);
  
  for(int l=0;l<nlayers;l++) {
    ASPWBres.Soil.HydraulicInput[l] = 0.0;
    ASPWBres.Soil.HydraulicOutput[l] = 0.0;
    ASPWBres.Soil.PlantExtraction[l] = 0.0;
    for(int n=0;n<ndailysteps;n++) {
      ASPWBres.Soil.HydraulicInput[l] += (-1.0)*std::min(ATres.extractionInst(l,n),0.0);
      ASPWBres.Soil.HydraulicOutput[l] += std::max(ATres.extractionInst(l,n),0.0);
      ASPWBres.Soil.PlantExtraction[l] += ATres.extractionInst(l,n);
    }
  }
  
  //STEP 9 - SOIL FLOWS
  double DeepDrainage = 0.0;
  double Infiltration = 0.0;
  double InfiltrationExcess = 0.0;
  double SaturationExcess = 0.0;
  double Runoff = 0.0;
  double CapillarityRise = 0.0;
  std::vector<double> sourceSinkVec(nlayers, 0.0);
  
  if(soilDomains != "none") {
    if(!plantWaterPools) {
      // determine water flows (no mass conservation)
      for(int l=0;l<nlayers;l++) {
        sourceSinkVec[l] -= (ASPWBres.Soil.PlantExtraction[l] + EherbVec[l]);
        if(l ==0) sourceSinkVec[l] -= Esoil;
      }
      soilWaterBalance_inner_c(SWBres, SWBcomm, soil,
                               RainfallInput, rainfallIntensity, Snowmelt, sourceSinkVec,
                               runon, lateralFlows, waterTableDepth,
                               infiltrationMode, infiltrationCorrection,
                               soilDomains,
                               ndailysteps, max_nsubsteps_soil);
      DeepDrainage = SWBres.deepDrainage_mm;
      Infiltration = SWBres.infiltration_mm;
      Runoff = SWBres.runoff_mm;
      InfiltrationExcess = SWBres.infiltrationExcess_mm;
      SaturationExcess = SWBres.saturationExcess_mm;
      CapillarityRise = SWBres.capillarityRise_mm;
    } else { //Apply soil flows to water pools
      arma::mat ExtractionPoolMat(numCohorts, nlayers);
      ExtractionPoolMat.fill(0.0);
      for(int c=0;c<numCohorts;c++) {
        //this is used to store extraction of a SINGLE plant cohort from all pools
        arma::mat& ExtractionPoolsCoh = ATres.extractionPools[c];
        for(int l=0;l<nlayers;l++) {
          for(int c2=0;c2<numCohorts;c2++) {
            ExtractionPoolMat(c2,l) += ExtractionPoolsCoh(c2,l)/poolProportions[c2];
          }
        }
      }
      //Create vector to store averaged soil moisture
      //Set Wsoil to zero
      double* Wsoil = new double[nlayers];
      for(int l = 0; l<nlayers;l++) Wsoil[l] = 0.0;
      for(int c=0;c<numCohorts;c++) {
        for(int l = 0; l<nlayers;l++) soil.setW(l,Wpool(c,l)); // this updates psi, theta, ... 
        for(int l=0;l<nlayers;l++) {
          sourceSinkVec[l] -= (ExtractionPoolMat(c,l) + EherbPools(c,l));
          if(l ==0) sourceSinkVec[l] -= EsoilPools[c];
        }
        soilWaterBalance_inner_c(SWBres, SWBcomm, soil,
                                 RainfallInput, rainfallIntensity, Snowmelt, sourceSinkVec,
                                 runon, lateralFlows, waterTableDepth,
                                 infiltrationMode, infiltrationCorrection,
                                 soilDomains,
                                 ndailysteps, max_nsubsteps_soil);
        
        DeepDrainage +=  SWBres.deepDrainage_mm*poolProportions[c];
        Runoff += SWBres.runoff_mm*poolProportions[c];
        Infiltration += SWBres.infiltration_mm*poolProportions[c];
        SaturationExcess += SWBres.saturationExcess_mm*poolProportions[c];
        InfiltrationExcess += SWBres.infiltrationExcess_mm*poolProportions[c];
        CapillarityRise += SWBres.capillarityRise_mm*poolProportions[c];
        
        //copy to Wpool and update averaged soil moisture (Wsoil)
        for(int l=0;l<nlayers;l++) {
          Wpool(c,l) = soil.getW(l);
          Wsoil[l] += soil.getW(l)*poolProportions[c];
        }
      }
      
      //Copy averaged soil moisture
      for(int l = 0; l<nlayers;l++) soil.setW(l, Wsoil[l]);
      delete[] Wsoil;
    }
  }
  //STEP 11 - Fire hazard
  if(fireHazardResults) {
    fccsHazard_c(ASPWBres.fccsbeh, ASPWBres.fccs, x, 
                 meteovec, 
                 ATres.plants.LFMC,
                 ATres.plants.StemPLC, 
                 slope);
  }
  
  // Arrange output
  ASPWBres.WaterBalance.PET = pet;
  ASPWBres.WaterBalance.Rain = ASPWB_comm.waterInputs.rain;
  ASPWBres.WaterBalance.Snow = ASPWB_comm.waterInputs.snow;
  ASPWBres.WaterBalance.NetRain = ASPWB_comm.waterInputs.netrain;
  ASPWBres.WaterBalance.Snowmelt = Snowmelt;
  ASPWBres.WaterBalance.Runon = runon;
  ASPWBres.WaterBalance.Infiltration = Infiltration;
  ASPWBres.WaterBalance.InfiltrationExcess = InfiltrationExcess;
  ASPWBres.WaterBalance.SaturationExcess = SaturationExcess;
  ASPWBres.WaterBalance.Runoff = Runoff;
  ASPWBres.WaterBalance.DeepDrainage = DeepDrainage;
  ASPWBres.WaterBalance.CapillarityRise = CapillarityRise;
  ASPWBres.WaterBalance.SoilEvaporation = Esoil;
  ASPWBres.WaterBalance.HerbTranspiration = std::accumulate(EherbVec.begin(), EherbVec.end(), 0.0);
  ASPWBres.WaterBalance.PlantExtraction = std::accumulate(ASPWBres.Soil.PlantExtraction.begin(), ASPWBres.Soil.PlantExtraction.end(), 0.0);
  ASPWBres.WaterBalance.Transpiration = std::accumulate(ATres.plants.Transpiration.begin(), ATres.plants.Transpiration.end(), 0.0);
  ASPWBres.WaterBalance.HydraulicRedistribution =std::accumulate(ASPWBres.Soil.HydraulicInput.begin(), ASPWBres.Soil.HydraulicInput.end(), 0.0);
  
  ASPWBres.Stand.LAI = LAIcell;
  ASPWBres.Stand.LAIherb = x.herbLAI; 
  ASPWBres.Stand.LAIlive = LAIcelllive;
  ASPWBres.Stand.LAIexpanded = LAIcellexpanded;
  ASPWBres.Stand.LAIdead = LAIcelldead;
  ASPWBres.Stand.Cm = Cm;
  ASPWBres.Stand.LgroundPAR = LgroundPAR;
  ASPWBres.Stand.LgroundSWR = LgroundSWR;
  
  //Copy final soil state to output
  for(int l=0;l<nlayers;l++) {
    ASPWBres.Soil.Psi[l] = soil.getPsi(l);
  }
}

Rcpp::List copyAdvancedSPWBResult_c(const AdvancedSPWB_RESULT& ASPWBres, ModelInput& x) {
  
  int nlayers = x.soil.getNlayers();
  int numCohorts = x.cohorts.CohortCode.size();
  int ntimesteps = x.control.advancedWB.ndailysteps;
  
  Rcpp::List l = Rcpp::List::create(_["cohorts"] = copyCohorts_c(x.cohorts),
                                    _["topography"] = copyTopo_c(ASPWBres.topo),
                                    _["weather"] = copyWeather_c(ASPWBres.meteovec),
                                    _["WaterBalance"] = copyWaterBalanceResult_c(ASPWBres.WaterBalance));
  if(x.control.results.soilResults) {
    l.push_back(copySoilResult_c(ASPWBres.Soil), "Soil");
  }
  if(x.control.results.standResults) {
    l.push_back(copyStandResult_c(ASPWBres.Stand), "Stand");
  }
  if(x.control.results.plantResults) {
    l.push_back(copyPlantAdvancedTranspirationResult_c(ASPWBres.ATres.plants, x), "Plants");
  }
  if(x.control.results.leafResults) {
    l.push_back(copyLeafAdvancedTranspirationResult_c(ASPWBres.ATres.sunlit, x), "SunlitLeaves");
    l.push_back(copyLeafAdvancedTranspirationResult_c(ASPWBres.ATres.shade, x), "ShadeLeaves");
  }
  
  NumericMatrix rhizoPsi = copyNumericMatrix_c(ASPWBres.ATres.rhizoPsi, numCohorts, nlayers);
  rhizoPsi.attr("dimnames") = List::create(x.cohorts.CohortCode, seq(1,nlayers));
  l.push_back(rhizoPsi, "RhizoPsi");
  
  if(x.control.results.subdailyResults) {
    NumericMatrix extractionInst = copyNumericMatrix_c(ASPWBres.ATres.extractionInst, nlayers, ntimesteps);
    extractionInst.attr("dimnames") = List::create(seq(1,nlayers), seq(1,ntimesteps));
    l.push_back(extractionInst, "ExtractionInst");
    l.push_back(copyPlantAdvancedTranspirationInstResult_c(ASPWBres.ATres.plants_inst, x), "PlantsInst");
    l.push_back(copyDirectDiffuseDayResult_c(ASPWBres.ATres.directDiffuseDay), "RadiationInputInst");
    l.push_back(copyLeafAdvancedTranspirationInstResult_c(ASPWBres.ATres.sunlit_inst, x), "SunlitLeavesInst");
    l.push_back(copyLeafAdvancedTranspirationInstResult_c(ASPWBres.ATres.shade_inst, x), "ShadeLeavesInst");
  }
  List lwrExtinctionList(ntimesteps);
  for(int n=0;n<ntimesteps;n++) {
    lwrExtinctionList[n] = copyLongWaveRadiationResult_c(ASPWBres.ATres.lwrExtinction[n]);
  }
  l.push_back(copyInstantaneousLightExtinctionAbsortionResult_c(ASPWBres.ATres.lightExtinctionAbsortion), "LightExtinction");
  l.push_back(lwrExtinctionList, "LWRExtinction");
  l.push_back(copyCanopyTurbulenceResult_c(ASPWBres.ATres.canopyTurbulence), "CanopyTurbulence");
  if(x.control.results.fireHazardResults) {
    l.push_back(copyFCCSResult_c(ASPWBres.fccs), "FireHazard");
  }
  l.attr("class") = CharacterVector::create("spwb_day","list");
  return(l);
}

void spwbDay_inner_c(SPWBCommunicationStructures& SPWBcomm, ModelInput& x, 
                     std::string date,
                     WeatherInputVector meteovec, 
                     double latitude, double elevation, double slope, double aspect,
                     const double runon, 
                     const std::vector<double>& lateralFlows, const double waterTableDepth) {
  
  if(std::isnan(meteovec.prec)) throw medfate::MedfateInternalError("Missing precipitation value");
  if(std::isnan(meteovec.tmin)) throw medfate::MedfateInternalError("Missing minimum temperature value");
  if(std::isnan(meteovec.tmax)) throw medfate::MedfateInternalError("Missing maximum temperature value");
  if(meteovec.tmin > meteovec.tmax) {
    double swap = meteovec.tmin;
    meteovec.tmin = meteovec.tmax;
    meteovec.tmax = swap;
  }
  if(std::isnan(meteovec.rhmax)) {
    meteovec.rhmax = 100.0;
  }
  if(std::isnan(meteovec.rhmin)) {
    double vp_tmin = meteoland::utils_saturationVP(meteovec.tmin);
    double vp_tmax = meteoland::utils_saturationVP(meteovec.tmax);
    meteovec.rhmin = std::min(meteovec.rhmax, 100.0*(vp_tmin/vp_tmax));
  }
  if(meteovec.rhmin > meteovec.rhmax) {
    // warning("rhmin > rhmax. Swapping values.");
    double swap = meteovec.rhmin;
    meteovec.rhmin = meteovec.rhmax;
    meteovec.rhmax = swap;
  }
  if(std::isnan(meteovec.wind)) meteovec.wind = x.control.weather.defaultWindSpeed; 
  if(meteovec.wind<0.1) meteovec.wind = 0.1; //Minimum windspeed abovecanopy

  if(std::isnan(meteovec.Catm)) meteovec.Catm = x.control.weather.defaultCO2;
  
  int month = std::atoi(date.substr(5,2).c_str());
  int J = julianDay_c(std::atoi(date.substr(0, 4).c_str()),std::atoi(date.substr(5,2).c_str()),std::atoi(date.substr(8,2).c_str()));
  double delta = solarDeclination_c(J);
  double solarConstant = solarConstant_c(J);
  double latrad = latitude * (M_PI/180.0);
  if(std::isnan(aspect)) aspect = 0.0;
  if(std::isnan(slope)) slope = 0.0;
  double asprad = aspect * (M_PI/180.0);
  double slorad = slope * (M_PI/180.0);
  double photoperiod = daylength_c(latrad, 0.0, 0.0, delta);
  if(std::isnan(meteovec.tday)) {
    meteovec.tday = meteoland::utils_averageDaylightTemperature(meteovec.tmin, meteovec.tmax);
  }
  if(std::isnan(meteovec.rad)) {
    // warning("Estimating solar radiation");
    double vpa = meteoland::utils_averageDailyVP(meteovec.tmin, meteovec.tmax, meteovec.rhmin, meteovec.rhmax);
    meteovec.rad = RDay_c(solarConstant, latrad, elevation,
                          slorad, asprad, delta, meteovec.tmax -meteovec.tmin, meteovec.tmax-meteovec.tmin,
                          vpa, meteovec.prec);
  }
  if(std::isnan(meteovec.pet)) {
    meteovec.pet = meteoland::penman(latrad, elevation, slorad, asprad, J, 
                                     meteovec.tmin, meteovec.tmax, meteovec.rhmin, meteovec.rhmax, meteovec.rad, meteovec.wind);
  }
  
  //Derive doy from date  
  int J0101 = julianDay_c(std::atoi(date.substr(0, 4).c_str()),1,1);
  int doy = J - J0101+1;
  
  std::vector<double> defaultRainfallIntensityPerMonth = x.control.weather.defaultRainfallIntensityPerMonth;
  if(std::isnan(meteovec.rint)) meteovec.rint = rainfallIntensity_c(month, meteovec.prec, defaultRainfallIntensityPerMonth);
  
  bool leafPhenology = x.control.phenology.leafPhenology;
  
  //Update phenology
  if(leafPhenology) {
    updatePhenology_c(x, doy, photoperiod, meteovec.tday);
    updateLeaves_c(x, meteovec.wind, false);
  }
  
  meteovec.tminPrev = meteovec.tmin;
  meteovec.tmaxPrev = meteovec.tmax;
  meteovec.tminNext = meteovec.tmin;

  if(x.control.transpirationMode=="Granier") {
    spwbDay_basic_c(SPWBcomm.BSPWBres, SPWBcomm.BSPWBcomm, x, 
                    meteovec, 
                    elevation, slope, aspect,
                    runon, 
                    lateralFlows, waterTableDepth);
  } else {
    spwbDay_advanced_c(SPWBcomm.ASPWBres, SPWBcomm.ASPWBcomm, x, 
                       meteovec, 
                       latitude, elevation, slope, aspect,
                       solarConstant, delta,
                       runon, 
                       lateralFlows, waterTableDepth);
  }
}
