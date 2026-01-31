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
#include "transpiration_basic_c.h"
#include "fuelstructure.h"
#include "firebehaviour.h"
#include "spwb_day_basic_c.h"
#include "tissuemoisture.h"
#include "soil.h"
#include "root_c.h"
#include <meteoland.h>


// Soil water balance with simple hydraulic model
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