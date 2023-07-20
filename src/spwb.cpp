// [[Rcpp::interfaces(r,cpp)]]
#define STRICT_R_HEADERS
#include <Rcpp.h>
#include <numeric>
#include <math.h>
#include "lightextinction.h"
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

CharacterVector getWeatherDates(DataFrame meteo){
  CharacterVector dateStrings;
  if(meteo.containsElementNamed("dates")){
    RObject vector = Rcpp::as<Rcpp::RObject>(meteo["dates"]);
    if(is<DateVector>(vector)) {
      DateVector dateVector = Rcpp::as<Rcpp::DateVector>(vector);
      CharacterVector dS(dateVector.size(), NA_STRING);
      for(int i=0;i< dateVector.size();i++) {
        Date d = dateVector[i];
        dS[i] = d.format("%Y-%m-%d");
      }
      dateStrings = dS;
    } else if(is<DatetimeVector>(vector)) {
      DatetimeVector datetimeVector = Rcpp::as<Rcpp::DatetimeVector>(vector);
      CharacterVector dS(datetimeVector.size(), NA_STRING);
      for(int i=0;i< datetimeVector.size();i++) {
        Datetime dt = datetimeVector[i];
        Date d(dt.getYear(), dt.getMonth(), dt.getDay());
        dS[i] = d.format("%Y-%m-%d");
      }
      dateStrings = dS;
    } else if(is<StringVector>(vector)) {
      dateStrings = Rcpp::as<Rcpp::StringVector>(vector);
    }
  } else {
    dateStrings = meteo.attr("row.names"); 
  }
  return(dateStrings);
}

NumericVector fccsHazard(List x, NumericVector meteovec, List transp, double slope) {
  List control = x["control"];
  
  double fireHazardStandardWind = control["fireHazardStandardWind"];
  double fireHazardStandardDFMC = control["fireHazardStandardDFMC"];
  
  DataFrame FCCSprops = Rcpp::as<Rcpp::DataFrame>(x["internalFCCS"]);
  DataFrame Plants = Rcpp::as<Rcpp::DataFrame>(transp["Plants"]);
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  
  double tmin = meteovec["tmin"];
  double tmax = meteovec["tmax"];
  double rhmin = meteovec["rhmin"];
  double rhmax = meteovec["rhmax"];
  double wind = meteovec["wind"];
  

  // Estimate moisture of dead fine fuels (Resco de Dios et al. 2015)
  double vp = meteoland::utils_averageDailyVP(tmin, tmax, rhmin, rhmax);
  double D = std::max(0.0, meteoland::utils_saturationVP(tmax) - vp);
  double fm_dead = 5.43 + 52.91*exp(-0.64*D); 
  
  //Calculate cohort canopy moisture to the average of canopy live and dead fuels, considering that a fraction of LAI is dead
  //proportionally to stem PLC (Ruffault et al. 2023)
  NumericVector LFMC = Plants["LFMC"];
  NumericVector PLC = Plants["StemPLC"];
  NumericVector deadFMC(LFMC.size(), fm_dead);
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
  if(w[1] > 0.0) ActFMC[1] = layerFuelAverageParameter(200.0, 10000.0, canopyFMC, cohLoading, cohHeight, cohCR);
  else ActFMC[1] = NA_REAL;
  if(w[2] > 0.0) ActFMC[2] = layerFuelAverageParameter(0.0, 200.0, canopyFMC, cohLoading, cohHeight, cohCR);
  else ActFMC[1] = NA_REAL;
      
  NumericVector MdeadSI = NumericVector::create(fm_dead, fm_dead, fm_dead, fm_dead, fm_dead); 
  if(!NumericVector::is_na(fireHazardStandardDFMC)) {
    MdeadSI = NumericVector::create(fireHazardStandardDFMC, fireHazardStandardDFMC, 
                                    fireHazardStandardDFMC, fireHazardStandardDFMC, fireHazardStandardDFMC); 
  }
  NumericVector MliveSI = NumericVector::create(90.0, 90.0, 60.0); //Default values (not actually used)
  List fccs;
  if(!NumericVector::is_na(fireHazardStandardWind)) {
    fccs = FCCSbehaviour(FCCSprops, MliveSI, MdeadSI, slope, fireHazardStandardWind); 
  } else {
    fccs = FCCSbehaviour(FCCSprops, MliveSI, MdeadSI, slope, wind); 
  }
  List surfaceFire = fccs["SurfaceFire"];
  List crownFire = fccs["CrownFire"];
  List firePotentials = fccs["FirePotentials"];
  NumericVector fireHazard = NumericVector::create(
    _["DFMC [%]"] = fm_dead,
    _["CFMC_understory [%]"] = ActFMC[1],
    _["CFMC_overstory [%]"] = ActFMC[2],
    _["ROS_surface [m/min]"] = surfaceFire["ROS [m/min]"],
    _["I_b_surface [kW/m]"] = surfaceFire["I_b [kW/m]"],
    _["t_r_surface [s]"] = surfaceFire["t_r [s]"],
    _["FL_surface [m]"] = surfaceFire["FL [m]"],
    _["Ic_ratio"] = crownFire["Ic_ratio"],
    _["ROS_crown [m/min]"] = crownFire["ROS_crown [m/min]"],
    _["I_b_crown [kW/m]"] = crownFire["I_b_crown [kW/m]"],
    _["t_r_crown [s]"] = crownFire["t_r_crown [s]"],
    _["FL_crown [m]"] = crownFire["FL_crown [m]"],
    _["SFP"] = firePotentials["SFP"],
    _["CFP"] = firePotentials["CFP"]
  );
  return(fireHazard);
}

// Soil water balance with simple hydraulic model
List spwbDay_basic(List x, NumericVector meteovec, 
              double elevation, double slope, double aspect,
              double runon=0.0, bool verbose = false) {

  //Control parameters
  List control = x["control"];
  bool snowpack = control["snowpack"];
  bool rockyLayerDrainage = control["rockyLayerDrainage"];
  String rhizosphereOverlap = control["rhizosphereOverlap"];
  bool plantWaterPools = (rhizosphereOverlap!="total");
  String soilFunctions = control["soilFunctions"];

  //Soil parameters
  List soil = x["soil"];
  int nlayers = Rcpp::as<Rcpp::NumericVector>(soil["dVec"]).size();
  
  List belowLayers = x["belowLayers"];
  NumericMatrix Wpool = Rcpp::as<Rcpp::NumericMatrix>(belowLayers["Wpool"]);
  NumericVector Wsoil = soil["W"];
  
  //Weather input
  double tday = meteovec["tday"];
  double pet = meteovec["pet"]; 
  double prec  = meteovec["prec"];
  double er  = meteovec["er"]; 
  double rad = NA_REAL; 
  if(meteovec.containsElementNamed("rad")) rad = meteovec["rad"];
    
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
  
  //Snow pack dynamics and hydrology input
  NumericVector hydroInputs = soilWaterInputs(soil, soilFunctions, prec, er, tday, rad, elevation,
                                             Cm, LgroundPAR, LgroundSWR, 
                                             runon,
                                             snowpack, true);
  
  NumericVector infilPerc, EsoilVec;
  double Eherb=0.0;
  
  if(!plantWaterPools) {
    //Soil infiltration and percolation
    infilPerc = soilInfiltrationPercolation(soil, soilFunctions, 
                                            hydroInputs["Input"],
                                            rockyLayerDrainage, true);
    //Evaporation from bare soil (if there is no snow)
    EsoilVec = soilEvaporation(soil, soilFunctions, pet, LgroundSWR, true);

    //Herbaceous transpiration
    Eherb = herbaceousTranspiration(pet, LherbSWR, herbLAI, soil, soilFunctions, true);
    //Copy soil status to x
    for(int c=0;c<numCohorts;c++) for(int l=0;l<nlayers;l++) Wpool(c,l) = Wsoil[l];
  } else {
    DataFrame belowdf = Rcpp::as<Rcpp::DataFrame>(x["below"]);
    NumericVector poolProportions = belowdf["poolProportions"];
    
    //Reset soil moisture
    for(int l=0;l<nlayers;l++) Wsoil[l] = 0.0;
    
    //Initialize result vectors
    infilPerc = NumericVector::create(_["Infiltration"] = 0.0, 
                                      _["Runoff"] = 0.0, 
                                      _["DeepDrainage"] = 0.0);
    EsoilVec = NumericVector(nlayers,0.0);
    for(int c=0;c<numCohorts;c++) {

      //Clone soil and copy moisture values from x
      List soil_c =  clone(soil);
      NumericVector W_c = soil_c["W"];
      for(int l=0;l<nlayers;l++) W_c[l] = Wpool(c,l);
      
      //Soil_c infiltration and percolation
      NumericVector infilPerc_c = soilInfiltrationPercolation(soil_c, soilFunctions, 
                                              hydroInputs["Input"],
                                              rockyLayerDrainage, true);
      //Evaporation from bare soil_c (if there is no snow)
      NumericVector EsoilVec_c = soilEvaporation(soil_c, soilFunctions, pet, LgroundSWR, true);
      //Herbaceous transpiration
      double Eherb_c = herbaceousTranspiration(pet, LherbSWR, herbLAI, soil_c, soilFunctions, true);
      
      //Copy result vectors
      infilPerc["Infiltration"] = infilPerc["Infiltration"] + poolProportions[c]*infilPerc_c["Infiltration"];
      infilPerc["Runoff"] = infilPerc["Runoff"] + poolProportions[c]*infilPerc_c["Runoff"];
      infilPerc["DeepDrainage"] = infilPerc["DeepDrainage"] + poolProportions[c]*infilPerc_c["DeepDrainage"];
      for(int l=0;l<nlayers;l++) EsoilVec[l] = EsoilVec[l] + poolProportions[c]*EsoilVec_c[l];
      Eherb = Eherb + Eherb_c;
      
      // Copy soil_c status back to x
      for(int l=0;l<nlayers;l++) {
        Wpool(c,l) = W_c[l];
        Wsoil[l] = Wsoil[l] + poolProportions[c]*Wpool(c,l); //weighted average for soil moisture
      }
    }
  }

  
  //Woody plant transpiration  
  // Rcout<<"hola";
  List transp = transpirationBasic(x, meteovec, elevation, true);
  // Rcout<<"hola2";
  NumericMatrix soilLayerExtract = Rcpp::as<Rcpp::NumericMatrix>(transp["Extraction"]);
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
  
  //Fire hazard
  bool fireHazardResults = control["fireHazardResults"];
  NumericVector fireHazard(1,NA_REAL);
  if(fireHazardResults) fireHazard = fccsHazard(x, meteovec, transp, slope);
  
  
  // Arrange output
  DataFrame Plants = Rcpp::as<Rcpp::DataFrame>(transp["Plants"]);
  NumericVector Eplant = Rcpp::as<Rcpp::NumericVector>(Plants["Transpiration"]);

  NumericVector psiVec = psi(soil, soilFunctions); //Calculate current soil water potential for output
  
  NumericVector topo = NumericVector::create(elevation, slope, aspect);
  topo.attr("names") = CharacterVector::create("elevation", "slope", "aspect");
  
  NumericVector DB = NumericVector::create(_["PET"] = pet, 
                                           _["Rain"] = hydroInputs["Rain"], _["Snow"] = hydroInputs["Snow"], 
                                           _["NetRain"] = hydroInputs["NetRain"], _["Snowmelt"] = hydroInputs["Snowmelt"],
                                           _["Runon"] = hydroInputs["Runon"], 
                                           _["Infiltration"] = infilPerc["Infiltration"], _["Runoff"] = infilPerc["Runoff"], _["DeepDrainage"] = infilPerc["DeepDrainage"],
                                           _["SoilEvaporation"] = sum(EsoilVec), _["HerbTranspiration"] = Eherb,
                                           _["PlantExtraction"] = sum(ExtractionVec), _["Transpiration"] = sum(Eplant),
                                           _["HydraulicRedistribution"] = sum(soilHydraulicInput));
  
  NumericVector Stand = NumericVector::create(_["LAI"] = LAIcell, _["LAIherb"] = herbLAI, 
                                              _["LAIlive"] = LAIcelllive,  _["LAIexpanded"] = LAIcellexpanded, _["LAIdead"] = LAIcelldead,
                                              _["Cm"] = Cm, _["LgroundPAR"] = LgroundPAR, _["LgroundSWR"] = LgroundSWR);
  DataFrame SB = DataFrame::create(_["SoilEvaporation"] = EsoilVec, 
                                   _["HydraulicInput"] = soilHydraulicInput, 
                                   _["HydraulicOutput"] = soilHydraulicOutput, 
                                   _["PlantExtraction"] = ExtractionVec, 
                                   _["psi"] = psiVec);
  List l = List::create(_["cohorts"] = clone(cohorts),
                        _["topography"] = topo,
                        _["weather"] = clone(meteovec),
                        _["WaterBalance"] = DB, 
                        _["Soil"] = SB,
                        _["Stand"] = Stand,
                        _["Plants"] = Plants);
  if(control["fireHazardResults"]) l.push_back(fireHazard, "FireHazard");
  l.attr("class") = CharacterVector::create("spwb_day","list");
  return(l);
}



// Soil water balance with Sperry or Cochard hydraulic and stomatal conductance models
List spwbDay_advanced(List x, NumericVector meteovec, 
             double latitude, double elevation, double slope, double aspect,
             double solarConstant, double delta, 
             double runon=0.0, bool verbose = false) {
  
  //Control parameters
  List control = x["control"];
  bool rockyLayerDrainage = control["rockyLayerDrainage"];
  bool snowpack = control["snowpack"];
  String rhizosphereOverlap = control["rhizosphereOverlap"];
  bool plantWaterPools = (rhizosphereOverlap!="total");
  String soilFunctions = control["soilFunctions"];
  int ntimesteps = control["ndailysteps"];
  String transpirationMode = control["transpirationMode"];

  //Soil parameters
  List soil = x["soil"];
  int nlayers = Rcpp::as<Rcpp::NumericVector>(soil["dVec"]).size();

  List belowLayers = x["belowLayers"];
  NumericMatrix Wpool = belowLayers["Wpool"];
  NumericVector Wsoil = soil["W"];
  
  //Meteo input
  double tmin = meteovec["tmin"];
  double tmax = meteovec["tmax"];
  double prec = meteovec["prec"];
  double rad = meteovec["rad"];
  double pet = meteovec["pet"];
  double er = meteovec["er"];
  
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

  //1. Leaf Phenology: Adjusted leaf area index
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
  
  //A.1 - Snow pack dynamics and soil water input
  NumericVector hydroInputs = soilWaterInputs(soil, soilFunctions, prec, er, tday, rad, elevation,
                                              Cm, LgroundPAR, LgroundSWR, 
                                              runon,
                                              snowpack, true);
  
  NumericVector infilPerc, EsoilVec;
  double Eherb = 0.0;
  if(!plantWaterPools) {
    //A.2 - Soil infiltration and percolation
    infilPerc = soilInfiltrationPercolation(soil, soilFunctions, 
                                            hydroInputs["Input"],
                                            rockyLayerDrainage, true);
    //B.1 - Evaporation from bare soil if there is no snow
    EsoilVec = soilEvaporation(soil, soilFunctions, pet, LgroundSWR, true);
    //Herbaceous transpiration
    Eherb = herbaceousTranspiration(pet, LherbSWR, herbLAI, soil, soilFunctions, true);
    
    //Copy soil status to x
    for(int c=0;c<numCohorts;c++) for(int l=0;l<nlayers;l++) Wpool(c,l) = Wsoil[l];
  } else {
    DataFrame belowdf = Rcpp::as<Rcpp::DataFrame>(x["below"]);
    NumericVector poolProportions = belowdf["poolProportions"];
    //Reset soil moisture
    for(int l=0;l<nlayers;l++) Wsoil[l] = 0.0;
    
    //Initialize result vectors
    infilPerc = NumericVector::create(_["Infiltration"] = 0.0, 
                                      _["Runoff"] = 0.0, 
                                      _["DeepDrainage"] = 0.0);
    EsoilVec = NumericVector(nlayers,0.0);
    for(int c=0;c<numCohorts;c++) {

      //Clone soil and copy moisture values from x
      List soil_c =  clone(soil);
      NumericVector W_c = soil_c["W"];
      for(int l=0;l<nlayers;l++) W_c[l] = Wpool(c,l);
      
      //Soil_c infiltration and percolation
      NumericVector infilPerc_c = soilInfiltrationPercolation(soil_c, soilFunctions, 
                                                              hydroInputs["Input"],
                                                              rockyLayerDrainage, true);
      //Evaporation from bare soil_c (if there is no snow)
      NumericVector EsoilVec_c = soilEvaporation(soil_c, soilFunctions, pet, LgroundSWR, true);
      //Herbaceous transpiration
      double Eherb_c = herbaceousTranspiration(pet, LherbSWR, herbLAI, soil_c, soilFunctions, true);
      //Copy result vectors
      infilPerc["Infiltration"] = infilPerc["Infiltration"] + poolProportions[c]*infilPerc_c["Infiltration"];
      infilPerc["Runoff"] = infilPerc["Runoff"] + poolProportions[c]*infilPerc_c["Runoff"];
      infilPerc["DeepDrainage"] = infilPerc["DeepDrainage"] + poolProportions[c]*infilPerc_c["DeepDrainage"];
      for(int l=0;l<nlayers;l++) EsoilVec[l] = EsoilVec[l] + poolProportions[c]*EsoilVec_c[l];
      Eherb = Eherb + Eherb_c;
      // Copy soil_c status back to x
      for(int l=0;l<nlayers;l++) {
        Wpool(c,l) = W_c[l];
        Wsoil[l] = Wsoil[l] + poolProportions[c]*Wpool(c,l); //weighted average for soil moisture
      }
    }
  }

  //B.2 - Canopy transpiration  
  List transp = transpirationAdvanced(x, meteovec, 
                                    latitude, elevation, slope, aspect, 
                                    solarConstant, delta, 
                                    hydroInputs["Interception"], hydroInputs["Snowmelt"], sum(EsoilVec), Eherb,
                                    verbose, NA_INTEGER, true);

  
  NumericMatrix soilLayerExtractInst = Rcpp::as<Rcpp::NumericMatrix>(transp["ExtractionInst"]);
  NumericMatrix RhizoPsi = Rcpp::as<Rcpp::NumericMatrix>(transp["RhizoPsi"]);
  DataFrame Plants = Rcpp::as<Rcpp::DataFrame>(transp["Plants"]);
  NumericVector Eplant = Plants["Transpiration"];
  
  List PlantsInst = Rcpp::as<Rcpp::List>(transp["PlantsInst"]);
  List EnergyBalance = Rcpp::as<Rcpp::List>(transp["EnergyBalance"]);
  
  //B.3 - Substract evaporated and extracted water from soil moisture 
  NumericVector EplantVec(nlayers, 0.0);
  NumericVector soilHydraulicInput(nlayers, 0.0); //Water that entered into the layer across all time steps
  NumericVector soilHydraulicOutput(nlayers, 0.0);  //Water that left the layer across all time steps
  for(int l=0;l<nlayers;l++) {
    for(int n=0;n<ntimesteps;n++) {
      soilHydraulicInput[l] += (-1.0)*std::min(soilLayerExtractInst(l,n),0.0);
      soilHydraulicOutput[l] += std::max(soilLayerExtractInst(l,n),0.0);
    }
    EplantVec[l] = sum(soilLayerExtractInst(l,_));
  }
  NumericVector psiVec = psi(soil, soilFunctions); //Calculate current soil water potential for output
  
  
  //Fire hazard
  bool fireHazardResults = control["fireHazardResults"];
  NumericVector fireHazard(1,NA_REAL);
  if(fireHazardResults) fireHazard = fccsHazard(x, meteovec, transp, slope);

  //Arrange output
  NumericVector topo = NumericVector::create(elevation, slope, aspect);
  topo.attr("names") = CharacterVector::create("elevation", "slope", "aspect");
  
  NumericVector DB = NumericVector::create(_["PET"] = pet,
                                           _["Rain"] = hydroInputs["Rain"],_["Snow"] = hydroInputs["Snow"],_["NetRain"] = hydroInputs["NetRain"], _["Snowmelt"] = hydroInputs["Snowmelt"],
                                           _["Runon"] = hydroInputs["Runon"], 
                                           _["Infiltration"] = infilPerc["Infiltration"], _["Runoff"] = infilPerc["Runoff"], _["DeepDrainage"] = infilPerc["DeepDrainage"],
                                           _["SoilEvaporation"] = sum(EsoilVec), _["HerbTranspiration"] = Eherb,
                                           _["PlantExtraction"] = sum(EplantVec), _["Transpiration"] = sum(Eplant),
                                           _["HydraulicRedistribution"] = sum(soilHydraulicInput));
  
  NumericVector Stand = NumericVector::create(_["LAI"] = LAIcell, _["LAIherb"] = herbLAI, 
                                              _["LAIlive"] = LAIcelllive, _["LAIexpanded"] = LAIcellexpanded, _["LAIdead"] = LAIcelldead,
                                              _["Cm"] = Cm, _["LgroundPAR"] = LgroundPAR, _["LgroundSWR"] = LgroundSWR);
  
  DataFrame SB = DataFrame::create(_["SoilEvaporation"] = EsoilVec, 
                                   _["HydraulicInput"] = soilHydraulicInput, 
                                   _["HydraulicOutput"] = soilHydraulicOutput, 
                                   _["PlantExtraction"] = EplantVec, 
                                   _["psi"] = psiVec);
  
  List l = List::create(_["cohorts"] = clone(cohorts),
                        _["topography"] = topo,
                        _["weather"] = clone(meteovec),
                        _["WaterBalance"] = DB, 
                        _["EnergyBalance"] = EnergyBalance,
                        _["Soil"] = SB, 
                        _["Stand"] = Stand, 
                        _["Plants"] = Plants,
                        _["RhizoPsi"] = RhizoPsi,
                        _["SunlitLeaves"] = transp["SunlitLeaves"],
                        _["ShadeLeaves"] = transp["ShadeLeaves"],
                        _["ExtractionInst"] = soilLayerExtractInst,
                        _["PlantsInst"] = PlantsInst,
                        _["SunlitLeavesInst"] = transp["SunlitLeavesInst"],
                        _["ShadeLeavesInst"] = transp["ShadeLeavesInst"],
                        _["LightExtinction"] = transp["LightExtinction"],
                        _["LWRExtinction"] = transp["LWRExtinction"],
                        _["CanopyTurbulence"] = transp["CanopyTurbulence"]);
  if(control["fireHazardResults"]) l.push_back(fireHazard, "FireHazard");
  l.attr("class") = CharacterVector::create("spwb_day","list");
  return(l);
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
//' @param modifyInput Boolean flag to indicate that the input \code{x} object is allowed to be modified during the simulation.
//' 
//' @details
//' The simulation functions allow using three different sub-models of transpiration and photosynthesis:
//' \itemize{
//'   \item{The sub-model corresponding to 'Granier' transpiration mode is illustrated by function \code{\link{transp_transpirationGranier}} and was described in De Caceres et al. (2015),
//'   and implements an approach originally described in Granier et al. (1999).} 
//'   \item{The sub-model corresponding to 'Sperry' transpiration mode is illustrated by function \code{\link{transp_transpirationSperry}} and was described in De Caceres et al. (2021), and
//'   implements a modelling approach originally described in Sperry et al. (2017).}  
//'   \item{The sub-model corresponding to 'Cochard' transpiration mode is illustrated by function \code{\link{transp_transpirationCochard}} and was described for model SurEau-Ecos v2.0 in Ruffault et al. (2022).} 
//' }
//' Simulations using the 'Sperry' or 'Cochard' transpiration mode are computationally much more expensive than 'Granier'.
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
//'       \item{\code{"SoilEvaporation"}: Water evaporated from the soil surface (in mm).}
//'       \item{\code{"HydraulicInput"}: Water entering each soil layer from other layers, transported via plant hydraulic network (in mm) (only for \code{transpirationMode = "Sperry"}).}
//'       \item{\code{"HydraulicOutput"}: Water leaving each soil layer (going to other layers or the transpiration stream) (in mm) (only for \code{transpirationMode = "Sperry"}).}
//'       \item{\code{"PlantExtraction"}: Water extracted by plants from each soil layer (in mm).}
//'       \item{\code{"psi"}: Soil water potential (in MPa).}
//'     }
//'   }
//'   \item{\code{"Stand"}: A named vector with with stand values for the simulated day, equivalent to one row of 'Stand' object returned by \code{\link{spwb}}.}
//'   \item{\code{"Plants"}: A data frame of results for each plant cohort (see \code{\link{transp_transpirationGranier}} or \code{\link{transp_transpirationSperry}}).}
//' }
//' The following items are only returned when \code{transpirationMode = "Sperry"} or  \code{transpirationMode = "Cochard"}:
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
//' data(exampleforestMED)
//' 
//' #Default species parameterization
//' data(SpParamsMED)
//' 
//' #Initialize control parameters
//' control <- defaultControl("Granier")
//' 
//' # Day to be simulated
//' d <- 100
//' meteovec <- unlist(examplemeteo[d,])
//' date <- rownames(examplemeteo)[d]
//' 
//' #Simulate water balance one day only (Granier mode)
//' examplesoil <- soil(defaultSoilParams(4))
//' x1 <- forest2spwbInput(exampleforestMED,examplesoil, SpParamsMED, control)
//' sd1 <- spwb_day(x1, date, meteovec,  
//'                 latitude = 41.82592, elevation = 100, slope=0, aspect=0) 
//' 
//' #Simulate water balance for one day only (Sperry mode)
//' control <- defaultControl("Sperry")
//' x2 <- forest2spwbInput(exampleforestMED, examplesoil, SpParamsMED, control)
//' sd2 <-spwb_day(x2, date, meteovec,
//'               latitude = 41.82592, elevation = 100, slope=0, aspect=0)
//' 
//' #Plot plant transpiration (see function 'plot.swb.day()')
//' plot(sd2)
//' 
//' #Simulate water balance for one day only (Cochard mode)
//' control <- defaultControl("Cochard")
//' x3 <- forest2spwbInput(exampleforestMED, examplesoil, SpParamsMED, control)
//' sd3 <-spwb_day(x3, date, meteovec,
//'               latitude = 41.82592, elevation = 100, slope=0, aspect=0)
//' 
//' 
//' #Simulate water and carbon balance for one day only (Granier mode)
//' control <- defaultControl("Granier")
//' x4  <- forest2growthInput(exampleforestMED,examplesoil, SpParamsMED, control)
//' sd4 <- growth_day(x4, date, meteovec,
//'                 latitude = 41.82592, elevation = 100, slope=0, aspect=0)
//' 
//' @name spwb_day
// [[Rcpp::export("spwb_day")]]
List spwbDay(List x, CharacterVector date, NumericVector meteovec, 
            double latitude, double elevation, double slope, double aspect,  
            double runon=0.0, 
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
  List soil = x["soil"];
  
  std::string c = as<std::string>(date[0]);
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
  
  //Update phenology
  if(leafPhenology) {
    updatePhenology(x, doy, photoperiod, tday);
    updateLeaves(x, wind, false);
  }
  
  
  double er = erFactor(doy, pet, prec);
  List s;
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
      Named("er") = er);
    s = spwbDay_basic(x, meteovec_bas,
                 elevation, slope, aspect, 
                 runon, verbose);
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
      Named("er") = er);
    s = spwbDay_advanced(x, meteovec_adv,
                 latitude, elevation, slope, aspect,
                 solarConstant, delta, 
                 runon, verbose);
  }
  // Rcout<<"hola4\n";
  return(s);
}

  

IntegerVector order_vector(NumericVector x) {
  if (is_true(any(duplicated(x)))) {
    Rf_warning("There are duplicates in 'x'; order not guaranteed to match that of R's base::order");
  }
  NumericVector sorted = clone(x).sort();
  return match(sorted, x);
}



void checkspwbInput(List x,  String transpirationMode, String soilFunctions) {
  
  List soil  = x["soil"];
  if(!x.containsElementNamed("above")) stop("above missing in spwbInput");
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  if(!above.containsElementNamed("LAI_live")) stop("LAI_live missing in spwbInput$above");
  if(!above.containsElementNamed("CR")) stop("CR missing in spwbInput$above");
  if(!above.containsElementNamed("H")) stop("H missing in spwbInput$above");
  
  if(!x.containsElementNamed("belowLayers")) stop("belowLayers missing in spwbInput");
  List belowLayers = Rcpp::as<Rcpp::List>(x["belowLayers"]);
  if(!belowLayers.containsElementNamed("V")) stop("V missing in spwbInput$belowLayers");
  if(transpirationMode=="Sperry"){
    if(!belowLayers.containsElementNamed("VGrhizo_kmax")) stop("VGrhizo_kmax missing in spwbInput$belowLayers");
    if(!belowLayers.containsElementNamed("VCroot_kmax")) stop("VCroot_kmax missing in spwbInput$belowLayers");
  }  
  
  if(!x.containsElementNamed("paramsPhenology")) stop("paramsPhenology missing in spwbInput");
  DataFrame paramsPhenology = Rcpp::as<Rcpp::DataFrame>(x["paramsPhenology"]);
  if(!paramsPhenology.containsElementNamed("Sgdd")) stop("Sgdd missing in spwbInput$paramsPhenology");
  if(!x.containsElementNamed("paramsInterception")) stop("paramsInterception missing in spwbInput");
  DataFrame paramsInterception = Rcpp::as<Rcpp::DataFrame>(x["paramsInterception"]);
  if(!paramsInterception.containsElementNamed("kPAR")) stop("kPAR missing in spwbInput$paramsInterception");
  if(!paramsInterception.containsElementNamed("g")) stop("g missing in spwbInput$paramsInterception");
  
  if(!x.containsElementNamed("paramsTranspiration")) stop("paramsTranspiration missing in spwbInput");
  DataFrame paramsTranspiration = Rcpp::as<Rcpp::DataFrame>(x["paramsTranspiration"]);
  if(transpirationMode=="Granier") {
    if(!paramsTranspiration.containsElementNamed("Psi_Extract")) stop("Psi_Extract missing in spwbInput$paramsTranspiration");
    if(!paramsTranspiration.containsElementNamed("WUE")) stop("WUE missing in spwbInput$paramsTranspiration");
  } else if(transpirationMode=="Sperry") {
    if(!paramsTranspiration.containsElementNamed("VCstem_kmax")) stop("VCstem_kmax missing in spwbInput$paramsTranspiration");
    if(!paramsTranspiration.containsElementNamed("VCstem_c")) stop("VCstem_c missing in spwbInput$paramsTranspiration");
    if(!paramsTranspiration.containsElementNamed("VCstem_d")) stop("VCstem_d missing in spwbInput$paramsTranspiration");
    if(!paramsTranspiration.containsElementNamed("VCroot_c")) stop("VCroot_c missing in spwbInput$paramsTranspiration");
    if(!paramsTranspiration.containsElementNamed("VCroot_d")) stop("VCroot_d missing in spwbInput$paramsTranspiration");
    if(!paramsTranspiration.containsElementNamed("Gswmax")) stop("Gswmax missing in spwbInput$paramsTranspiration");
    if(!paramsTranspiration.containsElementNamed("Vmax298")) stop("Vmax298 missing in spwbInput$paramsTranspiration");
    if(!paramsTranspiration.containsElementNamed("Jmax298")) stop("Jmax298 missing in spwbInput$paramsTranspiration");
  }
  if(transpirationMode=="Sperry") {
    if(!soil.containsElementNamed("VG_n")) stop("VG_n missing in soil");
    if(!soil.containsElementNamed("VG_alpha")) stop("VG_alpha missing in soil");
  }
  if(!soil.containsElementNamed("W")) stop("W missing in soil");
  if(!soil.containsElementNamed("dVec")) stop("dVec missing in soil");
  if(!soil.containsElementNamed("macro")) stop("macro missing in soil");
  if(soilFunctions=="SX") {
    if(!soil.containsElementNamed("clay")) stop("clay missing in soil");
    if(!soil.containsElementNamed("sand")) stop("sand missing in soil");
  }
  if(soilFunctions=="VG") {
    if(!soil.containsElementNamed("VG_n")) stop("VG_n missing in soil");
    if(!soil.containsElementNamed("VG_alpha")) stop("VG_alpha missing in soil");
    if(!soil.containsElementNamed("VG_theta_res")) stop("VG_theta_res missing in soil");
    if(!soil.containsElementNamed("VG_theta_sat")) stop("VG_theta_sat missing in soil");
  }
}

DataFrame defineWaterBalanceDailyOutput(CharacterVector dateStrings, NumericVector PET, String transpirationMode) {
  int numDays = dateStrings.length();
  
  NumericVector Precipitation(numDays), Evapotranspiration(numDays);
  NumericVector Runoff(numDays),Rain(numDays),Snow(numDays);
  NumericVector Snowmelt(numDays),NetRain(numDays);
  NumericVector Interception(numDays),Infiltration(numDays),DeepDrainage(numDays);
  NumericVector SoilEvaporation(numDays), HerbTranspiration(numDays), Transpiration(numDays),PlantExtraction(numDays);
  NumericVector HydraulicRedistribution(numDays, 0.0);
  
  DataFrame DWB = DataFrame::create(_["PET"]=PET, 
                          _["Precipitation"] = Precipitation, _["Rain"] = Rain, _["Snow"] = Snow, 
                          _["NetRain"]=NetRain, _["Snowmelt"] = Snowmelt, _["Infiltration"]=Infiltration, _["Runoff"]=Runoff, _["DeepDrainage"]=DeepDrainage, 
                            _["Evapotranspiration"]=Evapotranspiration,_["Interception"] = Interception, 
                            _["SoilEvaporation"]=SoilEvaporation, _["HerbTranspiration"] = HerbTranspiration,
                            _["PlantExtraction"] = PlantExtraction, _["Transpiration"]=Transpiration, 
                            _["HydraulicRedistribution"] = HydraulicRedistribution);
  DWB.attr("row.names") = dateStrings;
  return(DWB);
}
DataFrame defineSoilWaterBalanceDailyOutput(CharacterVector dateStrings, List soil, String transpirationMode) {
  int numDays = dateStrings.length();
  NumericVector W = soil["W"];
  int nlayers = W.length();
  
  NumericMatrix Eplantdays(numDays, nlayers);
  NumericMatrix HydrIndays(numDays, nlayers);
  NumericMatrix Wdays(numDays, nlayers); //Soil moisture content in relation to field capacity
  NumericMatrix psidays(numDays, nlayers);
  NumericMatrix MLdays(numDays, nlayers);
  NumericVector WaterTable(numDays, NA_REAL);
  NumericVector MLTot(numDays, 0.0);
  NumericVector SWE(numDays, 0.0);
  
  NumericVector Wini = soil["W"];
  Wdays(0,_) = clone(Wini);
  
  DataFrame SWB = DataFrame::create(_["W"]=Wdays, _["ML"]=MLdays,_["MLTot"]=MLTot,
                                    _["WTD"] = WaterTable,
                                    _["SWE"] = SWE, 
                                    _["PlantExt"]=Eplantdays,
                                    _["HydraulicInput"] = HydrIndays,
                                    _["psi"]=psidays); 
  SWB.attr("row.names") = dateStrings;

  return(SWB);  
}

DataFrame defineEnergyBalanceDailyOutput(CharacterVector dateStrings) {
  int numDays = dateStrings.length();
  NumericVector SWRcan(numDays, NA_REAL);
  NumericVector LWRcan(numDays, NA_REAL);
  NumericVector LEcan_heat(numDays, NA_REAL);
  NumericVector LEsoil_heat(numDays, NA_REAL);
  NumericVector Hcan_heat(numDays, NA_REAL);
  NumericVector Ebalcan(numDays, NA_REAL);
  NumericVector SWRsoil(numDays, NA_REAL);
  NumericVector LWRsoil(numDays, NA_REAL);
  NumericVector Ebalsoil(numDays, NA_REAL);
  NumericVector Hcansoil(numDays, NA_REAL);

  DataFrame DEB = DataFrame::create(_["SWRcan"] = SWRcan, _["LWRcan"] = LWRcan,
                                    _["LEcan"] = LEcan_heat, _["Hcan"] = Hcan_heat, _["Ebalcan"] = Ebalcan, 
                                    _["Hcansoil"] = Hcansoil, _["SWRsoil"] = SWRsoil, _["LWRsoil"] = LWRsoil, 
                                    _["LEsoil"] = LEsoil_heat, _["Ebalsoil"] = Ebalsoil);  
  DEB.attr("row.names") = dateStrings;
  return(DEB);
}
DataFrame defineTemperatureDailyOutput(CharacterVector dateStrings) {
  int numDays = dateStrings.length();
  
  NumericVector Tatm_mean(numDays, NA_REAL);
  NumericVector Tatm_min(numDays, NA_REAL);
  NumericVector Tatm_max(numDays, NA_REAL);
  NumericVector Tcan_mean(numDays, NA_REAL);
  NumericVector Tcan_min(numDays, NA_REAL);
  NumericVector Tcan_max(numDays, NA_REAL);
  NumericVector Tsoil_mean(numDays, NA_REAL);
  NumericVector Tsoil_min(numDays, NA_REAL);
  NumericVector Tsoil_max(numDays, NA_REAL);
  DataFrame DT = DataFrame::create(_["Tatm_mean"] = Tatm_mean, _["Tatm_min"] = Tatm_min, _["Tatm_max"] = Tatm_max,
                                   _["Tcan_mean"] = Tcan_mean, _["Tcan_min"] = Tcan_min, _["Tcan_max"] = Tcan_max,
                                     _["Tsoil_mean"] = Tsoil_mean, _["Tsoil_min"] = Tsoil_min, _["Tsoil_max"] = Tsoil_max);
  DT.attr("row.names") = dateStrings;
  return(DT);
}
NumericMatrix defineTemperatureLayersDailyOutput(CharacterVector dateStrings, DataFrame canopy) {
  int numDays = dateStrings.length();
  int ncanlayers = canopy.nrow();
  NumericMatrix DLT(numDays, ncanlayers);
  DLT.attr("dimnames") = List::create(dateStrings, seq(1,ncanlayers));
  return(DLT);
}
List defineSunlitShadeLeavesDailyOutput(CharacterVector dateStrings, DataFrame above) {
  int numDays = dateStrings.length();
  int numCohorts = above.nrow();
  NumericMatrix LeafPsiMin(numDays, numCohorts);
  NumericMatrix LeafPsiMax(numDays, numCohorts);
  NumericMatrix LeafGSWMin(numDays, numCohorts);
  NumericMatrix LeafGSWMax(numDays, numCohorts);
  NumericMatrix TempMin(numDays, numCohorts);
  NumericMatrix TempMax(numDays, numCohorts);
  LeafPsiMin.attr("dimnames") = List::create(dateStrings, above.attr("row.names")) ;
  LeafPsiMax.attr("dimnames") = List::create(dateStrings, above.attr("row.names")) ;
  LeafGSWMin.attr("dimnames") = List::create(dateStrings, above.attr("row.names")) ;
  LeafGSWMax.attr("dimnames") = List::create(dateStrings, above.attr("row.names")) ;
  TempMin.attr("dimnames") = List::create(dateStrings, above.attr("row.names")) ;
  TempMax.attr("dimnames") = List::create(dateStrings, above.attr("row.names")) ;
  List shade = List::create(Named("LeafPsiMin") = LeafPsiMin, 
                            Named("LeafPsiMax") = LeafPsiMax,
                            Named("TempMin") = TempMin, 
                            Named("TempMax") = TempMax,
                            Named("GSWMin") = LeafGSWMin,
                            Named("GSWMax") = LeafGSWMax);
  return(shade);
}

List definePlantWaterDailyOutput(CharacterVector dateStrings, DataFrame above, List soil, List control) {
  
  String transpirationMode = control["transpirationMode"];
  int numDays = dateStrings.length();
  NumericVector W = soil["W"];
  int nlayers = W.length();
  int numCohorts = above.nrow();
 
  NumericMatrix PlantStress(numDays, numCohorts);
  NumericMatrix PlantTranspiration(numDays, numCohorts);
  NumericMatrix PlantLAI(numDays, numCohorts);
  NumericMatrix PlantLAIlive(numDays, numCohorts);
  NumericMatrix StemPLC(numDays, numCohorts);
  NumericMatrix StemRWC(numDays, numCohorts), LeafRWC(numDays, numCohorts), LFMC(numDays, numCohorts);
  NumericMatrix PlantWaterBalance(numDays, numCohorts);
  
  PlantTranspiration.attr("dimnames") = List::create(dateStrings, above.attr("row.names"));
  PlantStress.attr("dimnames") = List::create(dateStrings, above.attr("row.names")) ;
  
  PlantLAI.attr("dimnames") = List::create(dateStrings, above.attr("row.names")) ;
  PlantLAIlive.attr("dimnames") = List::create(dateStrings, above.attr("row.names")) ;
  PlantTranspiration.attr("dimnames") = List::create(dateStrings, above.attr("row.names"));
  PlantStress.attr("dimnames") = List::create(dateStrings, above.attr("row.names")) ;
  StemPLC.attr("dimnames") = List::create(dateStrings, above.attr("row.names")) ;
  StemRWC.attr("dimnames") = List::create(dateStrings, above.attr("row.names")) ;
  LeafRWC.attr("dimnames") = List::create(dateStrings, above.attr("row.names")) ;
  LFMC.attr("dimnames") = List::create(dateStrings, above.attr("row.names")) ;
  PlantWaterBalance.attr("dimnames") = List::create(dateStrings, above.attr("row.names")) ;
  
  List plants;
  if(transpirationMode=="Granier") {
    NumericMatrix PlantPsi(numDays, numCohorts);
    NumericMatrix PlantGrossPhotosynthesis(numDays, numCohorts);
    NumericMatrix PlantAbsSWRFraction(numDays, numCohorts);
    NumericMatrix PlantFPAR(numDays, numCohorts);
    PlantFPAR.attr("dimnames") = List::create(dateStrings, above.attr("row.names")) ;
    PlantAbsSWRFraction.attr("dimnames") = List::create(dateStrings, above.attr("row.names")) ;
    PlantGrossPhotosynthesis.attr("dimnames") = List::create(dateStrings, above.attr("row.names")) ;
    PlantPsi.attr("dimnames") = List::create(dateStrings, above.attr("row.names")) ;
    plants = List::create(Named("LAI") = PlantLAI,
                          Named("LAIlive") = PlantLAIlive,
                          Named("FPAR") = PlantFPAR,
                          Named("AbsorbedSWRFraction") = PlantAbsSWRFraction,
                          Named("Transpiration") = PlantTranspiration,
                               Named("GrossPhotosynthesis") = PlantGrossPhotosynthesis,
                               Named("PlantPsi") = PlantPsi, 
                               Named("StemPLC") = StemPLC,
                               Named("PlantWaterBalance") = PlantWaterBalance,
                               Named("LeafRWC") = LeafRWC, 
                               Named("StemRWC") = StemRWC, 
                               Named("LFMC") = LFMC,
                               Named("PlantStress") = PlantStress);
  } else {
    NumericMatrix dEdP(numDays, numCohorts);
    NumericMatrix LeafPsiMin(numDays, numCohorts);
    NumericMatrix LeafPsiMax(numDays, numCohorts);
    NumericMatrix StemPsi(numDays, numCohorts);
    NumericMatrix RootPsi(numDays, numCohorts);
    List RhizoPsi(numCohorts);
    for(int c=0;c<numCohorts;c++) {
      NumericMatrix nm = NumericMatrix(numDays, nlayers);
      nm.attr("dimnames") = List::create(dateStrings, seq(1,nlayers)) ;
      RhizoPsi[c] = nm;
    }
    RhizoPsi.attr("names") = above.attr("row.names");
    
    NumericMatrix PlantNetPhotosynthesis(numDays, numCohorts);
    NumericMatrix PlantGrossPhotosynthesis(numDays, numCohorts);
    NumericMatrix PlantAbsSWR(numDays, numCohorts);
    NumericMatrix PlantNetLWR(numDays, numCohorts);
    dEdP.attr("dimnames") = List::create(dateStrings, above.attr("row.names")) ;
    LeafPsiMin.attr("dimnames") = List::create(dateStrings, above.attr("row.names")) ;
    LeafPsiMax.attr("dimnames") = List::create(dateStrings, above.attr("row.names")) ;
    StemPsi.attr("dimnames") = List::create(dateStrings, above.attr("row.names")) ;
    RootPsi.attr("dimnames") = List::create(dateStrings, above.attr("row.names")) ;
    PlantGrossPhotosynthesis.attr("dimnames") = List::create(dateStrings, above.attr("row.names")) ;
    PlantNetPhotosynthesis.attr("dimnames") = List::create(dateStrings, above.attr("row.names")) ;
    PlantAbsSWR.attr("dimnames") = List::create(dateStrings, above.attr("row.names")) ;
    PlantNetLWR.attr("dimnames") = List::create(dateStrings, above.attr("row.names")) ;
    
    plants = List::create(Named("LAI") = PlantLAI,
                          Named("LAIlive") = PlantLAIlive,
                          Named("AbsorbedSWR") = PlantAbsSWR,
                               Named("NetLWR") = PlantNetLWR,
                               Named("Transpiration") = PlantTranspiration,
                               Named("GrossPhotosynthesis") = PlantGrossPhotosynthesis,
                               Named("NetPhotosynthesis") = PlantNetPhotosynthesis,
                               Named("dEdP") = dEdP, 
                               Named("PlantWaterBalance") = PlantWaterBalance,
                               Named("LeafPsiMin") = LeafPsiMin, 
                               Named("LeafPsiMax") = LeafPsiMax, 
                               Named("LeafRWC") = LeafRWC, 
                               Named("StemRWC") = StemRWC, 
                               Named("StemPsi") = StemPsi, 
                               Named("StemPLC") = StemPLC, 
                               Named("RootPsi") = RootPsi, 
                               Named("RhizoPsi") = RhizoPsi, 
                               Named("LFMC") = LFMC,
                               Named("PlantStress") = PlantStress);
    
  }
  return(plants);
}
DataFrame defineFireHazardOutput(CharacterVector dateStrings){
  int numDays = dateStrings.length();
  
  NumericVector DFMC(numDays), CFMC_understory(numDays), CFMC_overstory(numDays);
  NumericVector ROS_surface(numDays), I_b_surface(numDays), t_r_surface(numDays),  Ic_ratio(numDays), FL_surface(numDays);
  NumericVector ROS_crown(numDays), I_b_crown(numDays), t_r_crown(numDays), FL_crown(numDays);
  NumericVector SFP(numDays), CFP(numDays);
  DataFrame df = DataFrame::create(_["DFMC"] = DFMC,
                                   _["CFMC_understory"] = CFMC_understory,
                                   _["CFMC_overstory"] = CFMC_overstory,
                                   _["ROS_surface"] = ROS_surface,
                                   _["I_b_surface"] = I_b_surface,
                                   _["t_r_surface"] = t_r_surface,
                                   _["FL_surface"] = FL_surface,
                                   _["Ic_ratio"] = Ic_ratio,
                                   _["ROS_crown"] = ROS_crown,
                                   _["I_b_crown"] = I_b_crown,
                                   _["t_r_crown"] = t_r_crown,
                                   _["FL_crown"] = FL_crown,
                                   _["SFP"] = SFP,
                                   _["CFP"] = CFP);
  
  df.attr("row.names") = dateStrings;
  return(df);
}
void fillWaterBalanceDailyOutput(DataFrame DWB, List sDay, int iday, String transpirationMode) {
  List db = sDay["WaterBalance"];
  NumericVector Precipitation = DWB["Precipitation"];
  NumericVector DeepDrainage = DWB["DeepDrainage"];
  NumericVector Infiltration = DWB["Infiltration"];
  NumericVector Runoff = DWB["Runoff"];
  NumericVector Rain = DWB["Rain"];
  NumericVector Snow = DWB["Snow"];
  NumericVector Snowmelt = DWB["Snowmelt"];
  NumericVector NetRain = DWB["NetRain"];
  NumericVector PlantExtraction = DWB["PlantExtraction"];
  NumericVector Transpiration = DWB["Transpiration"];
  NumericVector SoilEvaporation = DWB["SoilEvaporation"];
  NumericVector HerbTranspiration = DWB["HerbTranspiration"];
  NumericVector Interception = DWB["Interception"];
  NumericVector Evapotranspiration = DWB["Evapotranspiration"];
  DeepDrainage[iday] = db["DeepDrainage"];
  Infiltration[iday] = db["Infiltration"];
  Runoff[iday] = db["Runoff"];
  Rain[iday] = db["Rain"];
  Snow[iday] = db["Snow"];
  Precipitation[iday] = Rain[iday]+Snow[iday];
  Snowmelt[iday] = db["Snowmelt"];
  NetRain[iday] = db["NetRain"];
  PlantExtraction[iday] = db["PlantExtraction"];
  NumericVector HydraulicRedistribution = DWB["HydraulicRedistribution"];
  HydraulicRedistribution[iday] = db["HydraulicRedistribution"];
  Transpiration[iday] = db["Transpiration"];
  SoilEvaporation[iday] = db["SoilEvaporation"];
  HerbTranspiration[iday] = db["HerbTranspiration"];
  Interception[iday] = Rain[iday]-NetRain[iday];
  Evapotranspiration[iday] = Transpiration[iday] + SoilEvaporation[iday] + HerbTranspiration[iday] + Interception[iday];
}
void fillSoilWaterBalanceDailyOutput(DataFrame SWB, List soil, List sDay, 
                                     int iday, int numDays, String transpirationMode,
                                     String soilFunctions) {
  NumericVector W = soil["W"];
  int nlayers = W.length();
  NumericVector Water_FC = waterFC(soil, soilFunctions);
  
  List sb = sDay["Soil"];
  List Plants = sDay["Plants"];
  NumericVector psi = sb["psi"];
  NumericVector EplantVec = sb["PlantExtraction"];
  
  NumericVector MLTot = as<Rcpp::NumericVector>(SWB["MLTot"]);
  NumericVector WaterTable = as<Rcpp::NumericVector>(SWB["WTD"]);
  NumericVector SWE = as<Rcpp::NumericVector>(SWB["SWE"]);
  
  for(int l=0; l<nlayers; l++) {
    String wS = "W.";
    wS += (l+1);
    String mlS = "ML.";
    mlS += (l+1);
    String psiS = "psi.";
    psiS += (l+1);
    String peS = "PlantExt.";
    peS += (l+1);
    NumericVector Wdays = as<Rcpp::NumericVector>(SWB[wS]);
    NumericVector MLdays = as<Rcpp::NumericVector>(SWB[mlS]);
    NumericVector psidays = as<Rcpp::NumericVector>(SWB[psiS]);
    NumericVector Eplantdays = as<Rcpp::NumericVector>(SWB[peS]);
    psidays[iday] = psi[l];
    Eplantdays[iday] = EplantVec[l];
    if(iday<(numDays-1)) Wdays[iday+1] = W[l];
    MLdays[iday] = Wdays[iday]*Water_FC[l]; 
    MLTot[iday] = MLTot[iday] + MLdays[iday];
  }
  NumericVector HydrInVec = sb["HydraulicInput"];
  for(int l=0; l<nlayers; l++) {
    String hiS = "HydraulicInput.";
    hiS += (l+1);
    NumericVector HydrIndays = as<Rcpp::NumericVector>(SWB[hiS]);
    HydrIndays[iday] = HydrInVec[l];
  }
  SWE[iday] = soil["SWE"];
  WaterTable[iday] = waterTableDepth(soil, soilFunctions);
}

void fillEnergyBalanceTemperatureDailyOutput(DataFrame DEB, DataFrame DT, NumericMatrix DLT, List sDay, 
                                             int iday, bool multiLayerBalance) {
  List EB = Rcpp::as<Rcpp::List>(sDay["EnergyBalance"]);
  DataFrame Tinst = Rcpp::as<Rcpp::DataFrame>(EB["Temperature"]); 
  DataFrame CEBinst = Rcpp::as<Rcpp::DataFrame>(EB["CanopyEnergyBalance"]); 
  DataFrame SEBinst = Rcpp::as<Rcpp::DataFrame>(EB["SoilEnergyBalance"]); 
  NumericVector Tatm = Rcpp::as<Rcpp::NumericVector>(Tinst["Tatm"]);
  NumericVector Tcan = Rcpp::as<Rcpp::NumericVector>(Tinst["Tcan"]);
  NumericVector Tsoil = Rcpp::as<Rcpp::NumericVector>(Tinst["Tsoil.1"]);
  int ntimesteps = Tcan.length();
  int ncanlayers = DLT.ncol();
  double tstep = 86400.0/((double) ntimesteps);
  
  NumericVector SWRcan = DEB["SWRcan"];
  NumericVector LWRcan = DEB["LWRcan"];
  NumericVector LEcan_heat = DEB["LEcan"];
  NumericVector Hcan_heat = DEB["Hcan"];
  NumericVector Ebalcan = DEB["Ebalcan"];
  SWRcan[iday] = 0.000001*sum(Rcpp::as<Rcpp::NumericVector>(CEBinst["SWRcan"]))*tstep;
  LWRcan[iday] = 0.000001*sum(Rcpp::as<Rcpp::NumericVector>(CEBinst["LWRcan"]))*tstep;
  LEcan_heat[iday] = 0.000001*sum(Rcpp::as<Rcpp::NumericVector>(CEBinst["LEcan"]))*tstep;
  Hcan_heat[iday] = 0.000001*sum(Rcpp::as<Rcpp::NumericVector>(CEBinst["Hcan"]))*tstep;
  Ebalcan[iday] = 0.000001*sum(Rcpp::as<Rcpp::NumericVector>(CEBinst["Ebalcan"]))*tstep;
  NumericVector SWRsoil = DEB["SWRsoil"];
  NumericVector LWRsoil = DEB["LWRsoil"];
  NumericVector LEsoil_heat = DEB["LEsoil"];
  NumericVector Hcansoil = DEB["Hcansoil"];
  NumericVector Ebalsoil = DEB["Ebalsoil"];
  SWRsoil[iday] = 0.000001*sum(Rcpp::as<Rcpp::NumericVector>(SEBinst["SWRsoil"]))*tstep;
  LWRsoil[iday] = 0.000001*sum(Rcpp::as<Rcpp::NumericVector>(SEBinst["LWRsoil"]))*tstep;
  LEsoil_heat[iday] = 0.000001*sum(Rcpp::as<Rcpp::NumericVector>(SEBinst["LEsoil"]))*tstep;
  Hcansoil[iday] = 0.000001*sum(Rcpp::as<Rcpp::NumericVector>(SEBinst["Hcansoil"]))*tstep;
  Ebalsoil[iday] = 0.000001*sum(Rcpp::as<Rcpp::NumericVector>(SEBinst["Ebalsoil"]))*tstep;

  NumericVector Tatm_min = DT["Tatm_min"];
  NumericVector Tatm_max = DT["Tatm_max"];
  NumericVector Tatm_mean = DT["Tatm_mean"];
  NumericVector Tcan_min = DT["Tcan_min"];
  NumericVector Tcan_max = DT["Tcan_max"];
  NumericVector Tcan_mean = DT["Tcan_mean"];
  NumericVector Tsoil_min = DT["Tsoil_min"];
  NumericVector Tsoil_max = DT["Tsoil_max"];
  NumericVector Tsoil_mean = DT["Tsoil_mean"];
  Tatm_min[iday] = min(Tatm);
  Tatm_max[iday] = max(Tatm);
  Tatm_mean[iday] = sum(Tatm)/((double) ntimesteps);
  Tcan_min[iday] = min(Tcan);
  Tcan_max[iday] = max(Tcan);
  Tcan_mean[iday] = sum(Tcan)/((double) ntimesteps);
  Tsoil_min[iday] = min(Tsoil);
  Tsoil_max[iday] = max(Tsoil);
  Tsoil_mean[iday] = sum(Tsoil)/((double) ntimesteps);
  if(multiLayerBalance) {
    NumericMatrix LTinst = Rcpp::as<Rcpp::NumericMatrix>(EB["TemperatureLayers"]); 
    for(int l=0;l<ncanlayers;l++) DLT(iday, l) = sum(LTinst(_,l))/((double) ntimesteps);
  }
}
void fillPlantWaterDailyOutput(List x, List sunlit, List shade, List sDay, int iday, String transpirationMode) {
  List Plants = sDay["Plants"];
  
  NumericMatrix PlantStress= Rcpp::as<Rcpp::NumericMatrix>(x["PlantStress"]);
  NumericMatrix PlantTranspiration= Rcpp::as<Rcpp::NumericMatrix>(x["Transpiration"]);
  NumericMatrix PlantLAI= Rcpp::as<Rcpp::NumericMatrix>(x["LAI"]);
  NumericMatrix PlantLAIlive= Rcpp::as<Rcpp::NumericMatrix>(x["LAIlive"]);
  NumericMatrix StemPLC= Rcpp::as<Rcpp::NumericMatrix>(x["StemPLC"]);
  NumericMatrix StemRWC= Rcpp::as<Rcpp::NumericMatrix>(x["StemRWC"]);
  NumericMatrix LeafRWC= Rcpp::as<Rcpp::NumericMatrix>(x["LeafRWC"]);
  NumericMatrix LFMC= Rcpp::as<Rcpp::NumericMatrix>(x["LFMC"]);
  NumericMatrix PlantWaterBalance= Rcpp::as<Rcpp::NumericMatrix>(x["PlantWaterBalance"]);
  
  int numCohorts = PlantLAI.ncol();
  
  PlantTranspiration(iday,_) =  Rcpp::as<Rcpp::NumericVector>(Plants["Transpiration"]);
  PlantStress(iday,_) = Rcpp::as<Rcpp::NumericVector>(Plants["DDS"]);
  PlantLAI(iday,_) = Rcpp::as<Rcpp::NumericVector>(Plants["LAI"]);
  PlantLAIlive(iday,_) = Rcpp::as<Rcpp::NumericVector>(Plants["LAIlive"]);
  StemPLC(iday,_) = Rcpp::as<Rcpp::NumericVector>(Plants["StemPLC"]); 
  StemRWC(iday,_) = as<Rcpp::NumericVector>(Plants["StemRWC"]);
  LeafRWC(iday,_) = as<Rcpp::NumericVector>(Plants["LeafRWC"]); 
  LFMC(iday,_) = as<Rcpp::NumericVector>(Plants["LFMC"]); 
  PlantWaterBalance(iday,_) = Rcpp::as<Rcpp::NumericVector>(Plants["WaterBalance"]); 
  
  
  if(transpirationMode=="Granier") {
    NumericMatrix PlantGrossPhotosynthesis= Rcpp::as<Rcpp::NumericMatrix>(x["GrossPhotosynthesis"]);
    NumericMatrix PlantPsi = Rcpp::as<Rcpp::NumericMatrix>(x["PlantPsi"]);
    NumericMatrix PlantAbsSWRFraction= Rcpp::as<Rcpp::NumericMatrix>(x["AbsorbedSWRFraction"]);
    NumericMatrix PlantFPAR= Rcpp::as<Rcpp::NumericMatrix>(x["FPAR"]);
    
    PlantPsi(iday,_) =  Rcpp::as<Rcpp::NumericVector>(Plants["PlantPsi"]);  
    PlantGrossPhotosynthesis(iday,_) =  Rcpp::as<Rcpp::NumericVector>(Plants["GrossPhotosynthesis"]);  
    PlantAbsSWRFraction(iday,_) =  Rcpp::as<Rcpp::NumericVector>(Plants["AbsorbedSWRFraction"]);  
    PlantFPAR(iday,_) =  Rcpp::as<Rcpp::NumericVector>(Plants["FPAR"]);  
  } else {
    NumericMatrix RhizoPsiStep = Rcpp::as<Rcpp::NumericMatrix>(sDay["RhizoPsi"]);
    List PlantsInst = sDay["PlantsInst"];
    
    NumericMatrix dEdP = Rcpp::as<Rcpp::NumericMatrix>(x["dEdP"]);
    NumericMatrix LeafPsiMin = Rcpp::as<Rcpp::NumericMatrix>(x["LeafPsiMin"]);
    NumericMatrix LeafPsiMax = Rcpp::as<Rcpp::NumericMatrix>(x["LeafPsiMax"]);
    NumericMatrix StemPsi= Rcpp::as<Rcpp::NumericMatrix>(x["StemPsi"]);
    NumericMatrix RootPsi= Rcpp::as<Rcpp::NumericMatrix>(x["RootPsi"]);
    NumericMatrix LeafPsiMin_SL = Rcpp::as<Rcpp::NumericMatrix>(sunlit["LeafPsiMin"]);
    NumericMatrix LeafPsiMax_SL = Rcpp::as<Rcpp::NumericMatrix>(sunlit["LeafPsiMax"]);
    NumericMatrix LeafGSWMin_SL = Rcpp::as<Rcpp::NumericMatrix>(sunlit["GSWMin"]);
    NumericMatrix LeafGSWMax_SL = Rcpp::as<Rcpp::NumericMatrix>(sunlit["GSWMax"]);
    NumericMatrix LeafTempMin_SL = Rcpp::as<Rcpp::NumericMatrix>(sunlit["TempMin"]);
    NumericMatrix LeafTempMax_SL = Rcpp::as<Rcpp::NumericMatrix>(sunlit["TempMax"]);
    NumericMatrix LeafPsiMin_SH = Rcpp::as<Rcpp::NumericMatrix>(shade["LeafPsiMin"]);
    NumericMatrix LeafPsiMax_SH = Rcpp::as<Rcpp::NumericMatrix>(shade["LeafPsiMax"]);
    NumericMatrix LeafGSWMin_SH = Rcpp::as<Rcpp::NumericMatrix>(shade["GSWMin"]);
    NumericMatrix LeafGSWMax_SH = Rcpp::as<Rcpp::NumericMatrix>(shade["GSWMax"]);
    NumericMatrix LeafTempMin_SH = Rcpp::as<Rcpp::NumericMatrix>(shade["TempMin"]);
    NumericMatrix LeafTempMax_SH = Rcpp::as<Rcpp::NumericMatrix>(shade["TempMax"]);    
    List RhizoPsi = x["RhizoPsi"];
    NumericMatrix PlantNetPhotosynthesis= Rcpp::as<Rcpp::NumericMatrix>(x["NetPhotosynthesis"]);
    NumericMatrix PlantGrossPhotosynthesis= Rcpp::as<Rcpp::NumericMatrix>(x["GrossPhotosynthesis"]);
    NumericMatrix PlantAbsSWR= Rcpp::as<Rcpp::NumericMatrix>(x["AbsorbedSWR"]);
    NumericMatrix PlantNetLWR= Rcpp::as<Rcpp::NumericMatrix>(x["NetLWR"]);

    List SunlitLeavesInst = sDay["SunlitLeavesInst"]; 
    List ShadeLeavesInst = sDay["ShadeLeavesInst"]; 

    NumericMatrix SWR_SL = Rcpp::as<Rcpp::NumericMatrix>(SunlitLeavesInst["Abs_SWR"]);
    NumericMatrix SWR_SH = Rcpp::as<Rcpp::NumericMatrix>(ShadeLeavesInst["Abs_SWR"]);
    NumericMatrix LWR_SL = Rcpp::as<Rcpp::NumericMatrix>(SunlitLeavesInst["Net_LWR"]);
    NumericMatrix LWR_SH = Rcpp::as<Rcpp::NumericMatrix>(ShadeLeavesInst["Net_LWR"]);
    
    int ntimesteps = LWR_SH.ncol();
    double tstep = 86400.0/((double) ntimesteps);
    
    for(int j=0;j<numCohorts;j++) {
      for(int n=0;n<ntimesteps;n++){
        PlantAbsSWR(iday,j) += 0.000001*(SWR_SL(j,n)+SWR_SH(j,n))*tstep;
        PlantNetLWR(iday,j) += 0.000001*(LWR_SL(j,n)+LWR_SH(j,n))*tstep;
      }
    }
    
    List SunlitLeaves = sDay["SunlitLeaves"]; 
    List ShadeLeaves = sDay["ShadeLeaves"]; 
    LeafGSWMin_SL(iday,_) = Rcpp::as<Rcpp::NumericVector>(SunlitLeaves["GSWMin"]);
    LeafGSWMin_SH(iday,_) = Rcpp::as<Rcpp::NumericVector>(ShadeLeaves["GSWMin"]);
    LeafGSWMax_SL(iday,_) = Rcpp::as<Rcpp::NumericVector>(SunlitLeaves["GSWMax"]);
    LeafGSWMax_SH(iday,_) = Rcpp::as<Rcpp::NumericVector>(ShadeLeaves["GSWMax"]);
    LeafPsiMin_SL(iday,_) = Rcpp::as<Rcpp::NumericVector>(SunlitLeaves["LeafPsiMin"]);
    LeafPsiMax_SL(iday,_) = Rcpp::as<Rcpp::NumericVector>(SunlitLeaves["LeafPsiMax"]);
    LeafPsiMin_SH(iday,_) = Rcpp::as<Rcpp::NumericVector>(ShadeLeaves["LeafPsiMin"]);
    LeafPsiMax_SH(iday,_) = Rcpp::as<Rcpp::NumericVector>(ShadeLeaves["LeafPsiMax"]);
    LeafTempMin_SL(iday,_) = Rcpp::as<Rcpp::NumericVector>(SunlitLeaves["TempMin"]);
    LeafTempMax_SL(iday,_) = Rcpp::as<Rcpp::NumericVector>(SunlitLeaves["TempMax"]);
    LeafTempMin_SH(iday,_) = Rcpp::as<Rcpp::NumericVector>(ShadeLeaves["TempMin"]);
    LeafTempMax_SH(iday,_) = Rcpp::as<Rcpp::NumericVector>(ShadeLeaves["TempMax"]);
    
    PlantNetPhotosynthesis(iday,_) = Rcpp::as<Rcpp::NumericVector>(Plants["NetPhotosynthesis"]);
    PlantGrossPhotosynthesis(iday,_) = Rcpp::as<Rcpp::NumericVector>(Plants["GrossPhotosynthesis"]);
    LeafPsiMin(iday,_) = Rcpp::as<Rcpp::NumericVector>(Plants["LeafPsiMin"]);
    LeafPsiMax(iday,_) = Rcpp::as<Rcpp::NumericVector>(Plants["LeafPsiMax"]);
    RootPsi(iday,_) = Rcpp::as<Rcpp::NumericVector>(Plants["RootPsi"]); 
    StemPsi(iday,_) = Rcpp::as<Rcpp::NumericVector>(Plants["StemPsi"]); 
    
    dEdP(iday,_) = Rcpp::as<Rcpp::NumericVector>(Plants["dEdP"]); 
    for(int c=0;c<numCohorts;c++) {
      NumericMatrix nm = Rcpp::as<Rcpp::NumericMatrix>(RhizoPsi[c]);
      nm(iday,_) =  RhizoPsiStep(c,_);
    }

  }
}
void fillFireHazardOutput(DataFrame fireHazard, List sDay, int iday) {
  NumericVector fhd = sDay["FireHazard"];
  NumericVector DFMC = fireHazard["DFMC"];
  NumericVector CFMC_understory = fireHazard["CFMC_understory"];
  NumericVector CFMC_overstory = fireHazard["CFMC_overstory"];
  NumericVector ROS_surface = fireHazard["ROS_surface"];
  NumericVector I_b_surface = fireHazard["I_b_surface"];
  NumericVector t_r_surface = fireHazard["t_r_surface"];
  NumericVector FL_surface = fireHazard["FL_surface"];
  NumericVector Ic_ratio = fireHazard["Ic_ratio"];
  NumericVector ROS_crown = fireHazard["ROS_crown"];
  NumericVector I_b_crown = fireHazard["I_b_crown"];
  NumericVector t_r_crown = fireHazard["t_r_crown"];
  NumericVector FL_crown = fireHazard["FL_crown"];
  NumericVector SFP = fireHazard["SFP"];
  NumericVector CFP = fireHazard["CFP"];
  DFMC[iday] = fhd["DFMC [%]"];
  CFMC_understory[iday] = fhd["CFMC_understory [%]"];
  CFMC_overstory[iday] = fhd["CFMC_overstory [%]"];
  ROS_surface[iday] = fhd["ROS_surface [m/min]"];
  I_b_surface[iday] = fhd["I_b_surface [kW/m]"];
  t_r_surface[iday] = fhd["t_r_surface [s]"];
  FL_surface[iday] = fhd["FL_surface [m]"];
  Ic_ratio[iday] = fhd["Ic_ratio"];
  ROS_crown[iday] = fhd["ROS_crown [m/min]"];
  I_b_crown[iday] = fhd["I_b_crown [kW/m]"];
  t_r_crown[iday] = fhd["t_r_crown [s]"];
  FL_crown[iday] = fhd["FL_crown [m]"];
  SFP[iday] = fhd["SFP"];
  CFP[iday] = fhd["CFP"];
}
void printWaterBalanceResult(DataFrame DWB, List plantDWOL, List x,
                             NumericVector initialPlantContent, NumericVector initialSoilContent, double initialSnowContent,
                             String transpirationMode) {
  List soil = x["soil"];
  List control = x["control"];
  String soilFunctions = control["soilFunctions"];
  
  NumericVector finalPlantContent = plantWaterContent(x);
  NumericVector finalSoilContent = water(soil, soilFunctions);
  double finalSnowContent = soil["SWE"];
  Rcout<<"Final plant water content (mm): "<< sum(finalPlantContent)<<"\n";
  Rcout<<"Final soil water content (mm): "<< sum(finalSoilContent)<<"\n";
  Rcout<<"Final snowpack content (mm): "<< finalSnowContent<<"\n";
  
  NumericVector Precipitation = DWB["Precipitation"];
  NumericVector DeepDrainage = DWB["DeepDrainage"];
  NumericVector Infiltration = DWB["Infiltration"];
  NumericVector Runoff = DWB["Runoff"];
  NumericVector Rain = DWB["Rain"];
  NumericVector Snow = DWB["Snow"];
  NumericVector Snowmelt = DWB["Snowmelt"];
  NumericVector NetRain = DWB["NetRain"];
  NumericVector PlantExtraction = DWB["PlantExtraction"];
  NumericVector Transpiration = DWB["Transpiration"];
  NumericVector SoilEvaporation = DWB["SoilEvaporation"];
  NumericVector HerbTranspiration = DWB["HerbTranspiration"];
  NumericVector Interception = DWB["Interception"];
  NumericVector Evapotranspiration = DWB["Evapotranspiration"];
  
  NumericMatrix PlantWaterBalance = Rcpp::as<Rcpp::NumericMatrix>(plantDWOL["PlantWaterBalance"]);
  
  double Precipitationsum = sum(Precipitation);
  double Rainfallsum = sum(Rain);
  double NetRainsum = sum(NetRain);
  double Interceptionsum = sum(Interception);
  double SoilEvaporationsum = sum(SoilEvaporation);
  double Runoffsum  = sum(Runoff);
  double Infiltrationsum  = sum(Infiltration);
  double DeepDrainagesum = sum(DeepDrainage);
  double Transpirationsum = sum(Transpiration);
  double Snowmeltsum = sum(Snowmelt);
  double Snowsum = sum(Snow);
  double HerbTranspirationsum = sum(HerbTranspiration);
  
  double soil_wb = (Rainfallsum - Interceptionsum) + Snowmeltsum - Runoffsum - DeepDrainagesum - SoilEvaporationsum - HerbTranspirationsum - sum(PlantExtraction);
  double snowpack_wb = Snowsum - Snowmeltsum;
  Rcout<<"Change in plant water content (mm): "<< sum(finalPlantContent) - sum(initialPlantContent)<<"\n";
  Rcout<<"Plant water balance result (mm): "<< sum(PlantWaterBalance)<<"\n";
  Rcout<<"Change in soil water content (mm): "<< sum(finalSoilContent) - sum(initialSoilContent)<<"\n";
  Rcout<<"Soil water balance result (mm): "<< soil_wb<<"\n";
  Rcout<<"Change in snowpack water content (mm): "<< finalSnowContent - initialSnowContent<<"\n";
  Rcout<<"Snowpack water balance result (mm): "<< snowpack_wb<<"\n";
  Rcout<<"Water balance components:\n";
  Rcout<<"  Precipitation (mm) "  <<round(Precipitationsum) <<"\n";
  Rcout<<"  Rain (mm) "  <<round(Rainfallsum) <<" Snow (mm) "  <<round(Snowsum) <<"\n";
  Rcout<<"  Interception (mm) " << round(Interceptionsum)  <<" Net rainfall (mm) " << round(NetRainsum) <<"\n";
  Rcout<<"  Infiltration (mm) " << round(Infiltrationsum)  <<
    " Runoff (mm) " << round(Runoffsum) <<
      " Deep drainage (mm) "  << round(DeepDrainagesum)  <<"\n";
  Rcout<<"  Soil evaporation (mm) " << round(SoilEvaporationsum);
  Rcout<<"  Herbaceous transpiration (mm) " << round(HerbTranspirationsum);
  Rcout<<" Woody plant transpiration (mm) "  <<round(Transpirationsum) <<"\n";
  Rcout<<"  Plant extraction from soil (mm) " << round(sum(PlantExtraction));
  Rcout<<"  Plant water balance (mm) " << round(sum(PlantWaterBalance));
  NumericVector HydraulicRedistribution = DWB["HydraulicRedistribution"];
  Rcout<<" Hydraulic redistribution (mm) " << round(sum(HydraulicRedistribution));
  Rcout <<"\n";
}

//' Soil-plant water balance
//' 
//' Function \code{spwb()} is a water balance model that determines changes in soil moisture, 
//' soil water potentials, plant transpiration and drought stress at daily steps for a given forest stand 
//' during a period specified in the input climatic data. Function \code{pwb()} performs plant water balance 
//' only (i.e. soil moisture dynamics is an input) at daily steps for a given forest stand 
//' during a period specified in the input climatic data. On both simulation functions plant transpiration 
//' and photosynthesis processes are conducted with different level of detail depending on the transpiration mode.
//' 
//' @param x An object of class \code{\link{spwbInput}}.
//' @param meteo A data frame with daily meteorological data series. 
//' Row names of the data frame should correspond to date strings with format "yyyy-mm-dd" (see \code{\link{Date}}). Alternatively,
//' a column \code{dates} can contain \code{\link{Date}} or \code{\link{POSIXct}} classes.
//' The following columns are required and cannot have missing values:
//'   \itemize{
//'     \item{\code{MinTemperature}: Minimum temperature (in degrees Celsius).}
//'     \item{\code{MaxTemperature}: Maximum temperature (in degrees Celsius).}
//'     \item{\code{Precipitation}: Precipitation (in mm).}
//'   }
//' The following columns are required but can contain missing values (NOTE: missing values will raise warnings):
//'   \itemize{
//'     \item{\code{MinRelativeHumidity}: Minimum relative humidity (in percent).}
//'     \item{\code{MaxRelativeHumidity}: Maximum relative humidity (in percent).}
//'     \item{\code{Radiation}: Solar radiation (in MJ/m2/day).}
//'   }
//' The following columns are optional:
//'   \itemize{
//'     \item{\code{WindSpeed}: Above-canopy wind speed (in m/s). This column may not exist, or can be left with \code{NA} values. In both cases simulations will assume a constant value specified in \code{\link{defaultControl}}.}
//'     \item{\code{CO2}: Atmospheric (above-canopy) CO2 concentration (in ppm). This column may not exist, or can be left with \code{NA} values. In both cases simulations will assume a constant value specified in \code{\link{defaultControl}}.}
//'     \item{\code{Patm}: Atmospheric pressure (in kPa). This column may not exist, or can be left with \code{NA} values. In both cases, a value is estimated from elevation.}
//'   }
//' @param latitude Latitude (in degrees).
//' @param elevation,slope,aspect Elevation above sea level (in m), slope (in degrees) and aspect (in degrees from North).
//' @param CO2ByYear A named numeric vector with years as names and atmospheric CO2 concentration (in ppm) as values. Used to specify annual changes in CO2 concentration along the simulation (as an alternative to specifying daily values in \code{meteo}).
//' 
//' @details 
//' The simulation functions allow using three different sub-models of transpiration and photosynthesis:
//' \itemize{
//'   \item{The sub-model corresponding to 'Granier' transpiration mode is illustrated by function \code{\link{transp_transpirationGranier}} and was described in De Caceres et al. (2015),
//'   and implements an approach originally described in Granier et al. (1999).} 
//'   \item{The sub-model corresponding to 'Sperry' transpiration mode is illustrated by function \code{\link{transp_transpirationSperry}} and was described in De Caceres et al. (2021), and
//'   implements a modelling approach originally described in Sperry et al. (2017).}  
//'   \item{The sub-model corresponding to 'Cochard' transpiration mode is illustrated by function \code{\link{transp_transpirationCochard}} and was described for model SurEau-Ecos v2.0 in Ruffault et al. (2022).} 
//' }
//' Simulations using the 'Sperry' or 'Cochard' transpiration mode are computationally much more expensive than 'Granier'.
//' 
//' @return
//' Function \code{spwb} returns a list of class 'spwb' whereas function \code{pwb} returns a list of class 'pwb'. 
//' There are many elements in common in these lists, so they are listed here together:
//' \itemize{
//'   \item{\code{"latitude"}: Latitude (in degrees) given as input.} 
//'   \item{\code{"topography"}: Vector with elevation, slope and aspect given as input.} 
//'   \item{\code{"weather"}: A copy of the input weather data frame.}
//'   \item{\code{"spwbInput"}: An copy of the object \code{x} of class \code{\link{spwbInput}} given as input.}
//'   \item{\code{"spwbOutput"}: An copy of the final state of the object \code{x} of class \code{\link{spwbInput}}.}
//'   \item{\code{"WaterBalance"}: A data frame where different variables (in columns) are given for each simulated day (in rows):}
//'   \itemize{
//'     \item{\code{"PET"}: Potential evapotranspiration (in mm).}
//'     \item{\code{"Precipitation"}: Input precipitation (in mm).}
//'     \item{\code{"Rain"}: Precipitation as rain (in mm).}
//'     \item{\code{"Snow"}: Precipitation as snow (in mm).}
//'     \item{\code{"NetRain"}: Net rain, after accounting for interception (in mm).}
//'     \item{\code{"Infiltration"}: The amount of water infiltrating into the soil (in mm).}
//'     \item{\code{"Runoff"}: The amount of water exported via surface runoff (in mm).}
//'     \item{\code{"DeepDrainage"}: The amount of water exported via deep drainage (in mm).}
//'     \item{\code{"Evapotranspiration"}: Evapotranspiration (in mm).}
//'     \item{\code{"SoilEvaporation"}: Bare soil evaporation (in mm).}
//'     \item{\code{"HerbTranspiration"}: Transpiration due to the herbaceous layer (in mm).}
//'     \item{\code{"PlantExtraction"}: Amount of water extracted from soil by woody plants (in mm).}
//'     \item{\code{"Transpiration"}: Woody plant transpiration (in mm).}
//'     \item{\code{"HydraulicRedistribution"}: Water redistributed among soil layers, transported through the plant hydraulic network.}
//'   }
//'   \item{\code{"EnergyBalance"}: A data frame with the daily values of energy balance components for the soil and the canopy (only for \code{transpirationMode = "Sperry"} or \code{transpirationMode = "Cochard"}).}
//'   \item{\code{"Temperature"}: A data frame with the daily values of minimum/mean/maximum temperatures for the atmosphere (input), canopy and soil (only for \code{transpirationMode = "Sperry"} or \code{transpirationMode = "Cochard"}).}
//'   \item{\code{"Soil"}: A data frame where different variables (in columns) are given for each simulated day (in rows):}
//'   \itemize{
//'     \item{\code{"W.1"}, \code{...}, \code{"W.k"}: Relative soil moisture content (relative to field capacity) in each soil layer.}
//'     \item{\code{"ML.1"}, \code{...}, \code{"ML.k"}: Soil water volume in each soil layer (in L/m2).}
//'     \item{\code{"MLTot"}: Total soil water volume (in L/m2).}
//'     \item{\code{"SWE"}: Snow water equivalent (mm) of the snow pack.}
//'     \item{\code{"PlantExt.1"}, \code{...}, \code{"PlantExt.k"}: Plant extraction from each soil layer (in mm).}
//'     \item{\code{"HydraulicInput.1"}, \code{...}, \code{"HydraulicInput.k"}: Water that entered the layer coming from other layers and transported via the plant hydraulic network (in mm).}
//'     \item{\code{"psi.1"}, \code{...}, \code{"psi.k"}: Soil water potential in each soil layer (in MPa).}
//'   }
//'   \item{\code{"Stand"}: A data frame where different variables (in columns) are given for each simulated day (in rows):}
//'   \itemize{
//'     \item{\code{"LAI"}: LAI of the stand (including the herbaceous layer and live + dead leaves of woody plants) (in m2/m2).}
//'     \item{\code{"LAIherb"}: LAI of the herbaceous layer (in m2/m2).}
//'     \item{\code{"LAIlive"}: LAI of the woody plants assuming all leaves are unfolded (in m2/m2).}
//'     \item{\code{"LAIexpanded"}: LAI of the woody plants with leaves actually unfolded (in m2/m2).}
//'     \item{\code{"LAIdead"}: LAI of the woody plants corresponding to dead leaves (in m2/m2).}
//'     \item{\code{"Cm"}: Water retention capacity of the canopy (in mm) (accounting for leaf phenology).}
//'     \item{\code{"LgroundPAR"}: The percentage of PAR that reaches the ground (accounting for leaf phenology).}
//'     \item{\code{"LgroundSWR"}: The percentage of SWR that reaches the ground (accounting for leaf phenology).}
//'   }
//'   \item{\code{"Plants"}: A list of daily results for plant cohorts (see below).}
//'   \item{\code{"subdaily"}: A list of objects of class \code{\link{spwb_day}}, one per day simulated (only if required in \code{control} parameters, see \code{\link{defaultControl}}).}
//' }
//' 
//' When \code{transpirationMode = "Granier"}, element \code{"Plants"} is a list with the following subelements:
//'   \itemize{
//'     \item{\code{"LAI"}: A data frame with the daily leaf area index for each plant cohort.}
//'     \item{\code{"LAIlive"}: A data frame with the daily leaf area index for each plant cohort, assuming all leaves are unfolded (in m2/m2).}
//'     \item{\code{"FPAR"}: A data frame with the fraction of PAR at the canopy level of each plant cohort. }
//'     \item{\code{"AbsorbedSWRFraction"}: A data frame with the fraction of SWR absorbed by each plant cohort. }
//'     \item{\code{"Transpiration"}: A data frame with the amount of daily transpiration (in mm) for each plant cohort.}
//'     \item{\code{"GrossPhotosynthesis"}: A data frame with the amount of daily gross photosynthesis (in g C·m-2) for each plant cohort. }
//'     \item{\code{"PlantPsi"}: A data frame with the average daily water potential of each plant (in MPa).}
//'     \item{\code{"StemPLC"}: A data frame with the average daily proportion of stem conductance loss of each plant ([0-1]).}
//'     \item{\code{"PlantWaterBalance"}: A data frame with the daily balance between transpiration and soil water extraction for each plant cohort. }
//'     \item{\code{"LeafRWC"}: A data frame with the average daily leaf relative water content of each plant (in percent).}
//'     \item{\code{"StemRWC"}: A data frame with the average daily stem relative water content of each plant (in percent). }
//'     \item{\code{"LFMC"}: A data frame with the daily live fuel moisture content (in percent of dry weight).}
//'     \item{\code{"PlantStress"}: A data frame with the amount of daily stress [0-1] suffered by each plant cohort (relative whole-plant conductance).}
//'   }
//' If \code{transpirationMode="Sperry"} or \code{transpirationMode="Cochard"}, element \code{"Plants"} is a list with the following subelements:
//'   \itemize{
//'     \item{\code{"LAI"}: A data frame with the daily leaf area index for each plant cohort.}
//'     \item{\code{"AbsorbedSWR"}: A data frame with the daily SWR absorbed by each plant cohort.}
//'     \item{\code{"NetLWR"}: A data frame with the daily net LWR by each plant cohort.}
//'     \item{\code{"Transpiration"}: A data frame with the amount of daily transpiration (in mm) for each plant cohorts.}
//'     \item{\code{"GrossPhotosynthesis"}: A data frame with the amount of daily gross photosynthesis (in g C·m-2) for each plant cohort. }
//'     \item{\code{"NetPhotosynthesis"}: A data frame with the amount of daily net photosynthesis (in g C·m-2) for each plant cohort. }
//'     \item{\code{"dEdP"}: A data frame with mean daily values of soil-plant conductance (derivative of the supply function) for each plant cohort.}
//'     \item{\code{"PlantWaterBalance"}: A data frame with the daily balance between transpiration and soil water extraction for each plant cohort. }
//'     \item{\code{"SunlitLeaves"} and \code{"ShadeLeaves"}: A list with daily results for sunlit and shade leaves:
//'       \itemize{
//'         \item{\code{"PsiMin"}: A data frame with the minimum (midday) daily sunlit or shade leaf water potential (in MPa). }
//'         \item{\code{"PsiMax"}: A data frame with the maximum (predawn) daily sunlit or shade leaf water potential (in MPa). }
//'       }
//'     }
//'     \item{\code{"LeafPsiMin"}: A data frame with the minimum (midday) daily (average) leaf water potential of each plant (in MPa).}
//'     \item{\code{"LeafPsiMax"}: A data frame with the maximum (predawn) daily (average) leaf water potential of each plant (in MPa).}
//'     \item{\code{"LeafRWC"}: A data frame with the average daily leaf relative water content of each plant (in percent).}
//'     \item{\code{"StemRWC"}: A data frame with the average daily stem relative water content of each plant (in percent). }
//'     \item{\code{"LFMC"}: A data frame with the daily live fuel moisture content (in percent of dry weight).}
//'     \item{\code{"StemPsi"}: A data frame with the minimum daily stem water potential of each plant (in MPa). }
//'     \item{\code{"StemPLC"}: A data frame with the average daily proportion of stem conductance loss of each plant ([0-1]).}
//'     \item{\code{"RootPsi"}: A data frame with the minimum daily root water potential of each plant (in MPa). }
//'     \item{\code{"RhizoPsi"}: A list of data frames (one per plant cohort) with the minimum daily root water potential of each plant (in MPa).}
//'     \item{\code{"PlantStress"}: A data frame with the amount of daily stress [0-1] suffered by each plant cohort (relative whole-plant conductance).}
//'   }
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
//' \code{\link{spwbInput}}, \code{\link{spwb_day}}, \code{\link{plot.spwb}}, 
//' \code{\link{forest}}
//' 
//' @examples
//' #Load example daily meteorological data
//' data(examplemeteo)
//' 
//' #Load example plot plant data
//' data(exampleforestMED)
//' 
//' #Default species parameterization
//' data(SpParamsMED)
//' 
//' #Initialize soil with default soil params (4 layers)
//' examplesoil <- soil(defaultSoilParams(4))
//' 
//' #Initialize control parameters
//' control <- defaultControl("Granier")
//' 
//' #Initialize input
//' x1 <- forest2spwbInput(exampleforestMED,examplesoil, SpParamsMED, control)
//' 
//' #Call simulation function
//' S1 <- spwb(x1, examplemeteo, latitude = 41.82592, elevation = 100)
//' 
//' #Plot results
//' plot(S1)
//' 
//' #Monthly summary (averages) of soil water balance
//' summary(S1, freq="months",FUN=mean, output="Soil")
//'                   
//' \donttest{
//' #Switch to 'Sperry' transpiration mode
//' control <- defaultControl("Sperry")
//' 
//' #Initialize input
//' x2 <- forest2spwbInput(exampleforestMED,examplesoil, SpParamsMED, control)
//' 
//' #Call simulation function
//' S2 <- spwb(x2, examplemeteo, latitude = 41.82592, elevation = 100)
//' 
//' #Switch to 'Cochard' transpiration mode
//' control <- defaultControl("Cochard")
//' 
//' #Initialize input
//' x3 <- forest2spwbInput(exampleforestMED,examplesoil, SpParamsMED, control)
//' 
//' #Call simulation function
//' S3 <- spwb(x3, examplemeteo, latitude = 41.82592, elevation = 100)
//' }
//'                 
//' @name spwb
// [[Rcpp::export("spwb")]]
List spwb(List x, DataFrame meteo, double latitude, double elevation = NA_REAL, double slope = NA_REAL, double aspect = NA_REAL,
          NumericVector CO2ByYear = NumericVector(0)) {
  List control = x["control"];
  String transpirationMode = control["transpirationMode"];
  String soilFunctions = control["soilFunctions"];
  String cavitationRefill = control["cavitationRefill"];
  bool verbose = control["verbose"];
  
  bool leafPhenology = control["leafPhenology"];
  bool unlimitedSoilWater = control["unlimitedSoilWater"];
  bool multiLayerBalance = control["multiLayerBalance"];
  checkspwbInput(x,transpirationMode, soilFunctions);
  
  //Store input
  List spwbInput = x; // Store initial object
  x = clone(x); //Ensure a copy will be modified
  
  List soil = x["soil"];
  
  

  //Meteorological input    
  NumericVector MinTemperature, MaxTemperature;
  NumericVector MinRelativeHumidity, MaxRelativeHumidity;
  NumericVector Radiation;
  
  if(NumericVector::is_na(latitude)) stop("Value for 'latitude' should not be missing.");
  double latrad = latitude * (M_PI/180.0);
  if(NumericVector::is_na(aspect)) aspect = 0.0;
  if(NumericVector::is_na(slope)) slope = 0.0;
  double asprad = aspect * (M_PI/180.0);
  double slorad = slope * (M_PI/180.0);
  
  
  if(!meteo.containsElementNamed("Precipitation")) stop("Please include variable 'Precipitation' in weather input.");
  NumericVector Precipitation = meteo["Precipitation"];
  int numDays = Precipitation.size();
  NumericVector WindSpeed(numDays, NA_REAL);
  if(meteo.containsElementNamed("WindSpeed")) WindSpeed = meteo["WindSpeed"];
  
  NumericVector PET(numDays, NA_REAL);
  
  
  if(NumericVector::is_na(elevation)) stop("Value for 'elevation' should not be missing.");
  
  if(!meteo.containsElementNamed("MinTemperature")) stop("Please include variable 'MinTemperature' in weather input.");
  MinTemperature = meteo["MinTemperature"];
  if(!meteo.containsElementNamed("MaxTemperature")) stop("Please include variable 'MaxTemperature' in weather input.");
  MaxTemperature = meteo["MaxTemperature"];
  if(!meteo.containsElementNamed("MinRelativeHumidity")) stop("Please include variable 'MinRelativeHumidity' in weather input.");
  MinRelativeHumidity = meteo["MinRelativeHumidity"];
  if(!meteo.containsElementNamed("MaxRelativeHumidity")) stop("Please include variable 'MaxRelativeHumidity' in weather input.");
  MaxRelativeHumidity = meteo["MaxRelativeHumidity"];
  if(!meteo.containsElementNamed("Radiation")) stop("Please include variable 'Radiation' in weather input.");
  Radiation = meteo["Radiation"];
  
  if(any(is_na(Precipitation))) stop("Missing values in 'Precipitation'");
  if(any(is_na(MinTemperature))) stop("Missing values in 'MinTemperature'");
  if(any(is_na(MaxTemperature))) stop("Missing values in 'MaxTemperature'");
  if(any(is_na(MinRelativeHumidity))) warning("Missing values in 'MinRelativeHumidity' were estimated from temperature range");
  if(any(is_na(MaxRelativeHumidity))) warning("Missing values in 'MaxRelativeHumidity' were assumed to be 100");
  if(any(is_na(Radiation))) warning("Missing values in 'Radiation' were estimated");
  
  NumericVector CO2(Precipitation.length(), NA_REAL);
  if(meteo.containsElementNamed("CO2")) {
    CO2 = meteo["CO2"];
    if(verbose) {
      Rcout<<"CO2 taken from input column 'CO2'\n";
    }
    if(any(is_na(CO2))) stop("Missing values in 'CO2'");
  }
  NumericVector Patm(Precipitation.length(), NA_REAL);
  if(meteo.containsElementNamed("Patm")) {
    Patm = meteo["Patm"];
    if(verbose) {
      Rcout<<"Patm taken from input column 'Patm'\n";
    }
  }
  
  IntegerVector DOY, JulianDay;
  NumericVector Photoperiod;
  bool doy_input = false, photoperiod_input = false, julianday_input = false;
  if(meteo.containsElementNamed("DOY")) {
    DOY = meteo["DOY"];
    doy_input = true;
    if(verbose) {
      Rcout<<"DOY taken from input column 'DOY'\n";
    }
  }
  if(meteo.containsElementNamed("Photoperiod")) {
    Photoperiod = meteo["Photoperiod"];
    photoperiod_input = true;
    if(verbose) {
      Rcout<<"Photoperiod taken from input column 'Photoperiod'\n";
    }
  }
  if(meteo.containsElementNamed("JulianDay")) {
    JulianDay = meteo["JulianDay"];
    julianday_input = true;
    if(verbose) {
      Rcout<<"Julian day taken from input column 'JulianDay'\n";
    }
  }
  
  // Dates
  CharacterVector dateStrings = getWeatherDates(meteo);
  if(!doy_input) DOY = date2doy(dateStrings);
  if(!photoperiod_input) Photoperiod = date2photoperiod(dateStrings, latrad);
  

  //Canopy scalars
  DataFrame canopy = Rcpp::as<Rcpp::DataFrame>(x["canopy"]);
  
  //Plant input
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);

  //Detailed subday results
  List subdailyRes(numDays);
  
  //Stand output variables
  NumericVector LAI(numDays), LAIherb(numDays), LAIexpanded(numDays),LAIlive(numDays),LAIdead(numDays);
  NumericVector Cm(numDays);
  NumericVector LgroundPAR(numDays);
  NumericVector LgroundSWR(numDays);
  
  
  //Water balance output variables
  DataFrame DWB = defineWaterBalanceDailyOutput(dateStrings, PET, transpirationMode);
  DataFrame SWB = defineSoilWaterBalanceDailyOutput(dateStrings, soil, transpirationMode);
  
  //EnergyBalance output variables
  DataFrame DEB = defineEnergyBalanceDailyOutput(dateStrings);
  DataFrame DT = defineTemperatureDailyOutput(dateStrings);
  NumericMatrix DLT;
  if(transpirationMode!="Granier") DLT =  defineTemperatureLayersDailyOutput(dateStrings, canopy);
  
  //Plant output variables
  List sunlitDO = defineSunlitShadeLeavesDailyOutput(dateStrings, above);
  List shadeDO = defineSunlitShadeLeavesDailyOutput(dateStrings, above);
  List plantDWOL = definePlantWaterDailyOutput(dateStrings, above, soil, control);

  //Fire hazard output variables
  DataFrame fireHazard;
  if(control["fireHazardResults"]) fireHazard = defineFireHazardOutput(dateStrings);

  NumericVector initialSoilContent = water(soil, soilFunctions);
  NumericVector initialPlantContent = plantWaterContent(x);
  double initialSnowContent = soil["SWE"];
  if(verbose) {
    Rcout<<"Initial plant water content (mm): "<< sum(initialPlantContent)<<"\n";
    Rcout<<"Initial soil water content (mm): "<< sum(initialSoilContent)<<"\n";
    Rcout<<"Initial snowpack content (mm): "<< initialSnowContent<<"\n";
  }
  
  bool error_occurence = false;
  if(verbose) Rcout << "Performing daily simulations\n";
  NumericVector Eplanttot(numDays,0.0);
  List s;
  std::string yearString;
  for(int i=0;(i<numDays) && (!error_occurence);i++) {
     std::string c = as<std::string>(dateStrings[i]);
     yearString = c.substr(0, 4);
     if(verbose) {
        if(DOY[i]==1 || i==0) {
          Rcout<<"\n [Year "<< yearString << "]:";
        } 
        else if(i%10 == 0) Rcout<<".";//<<i;
      } 
      
      double wind = WindSpeed[i];
      if(NumericVector::is_na(wind)) wind = control["defaultWindSpeed"]; //Default 1 m/s -> 10% of fall every day
      if(wind<0.1) wind = 0.1; //Minimum windspeed abovecanopy
      
      
      double Catm = CO2[i];
      //If missing, use
      if(NumericVector::is_na(Catm)) {
        if(CO2ByYear.attr("names") != R_NilValue) Catm = CO2ByYear[yearString];
      }
      //If still missing, use default control value
      if(NumericVector::is_na(Catm)) {
        Catm = control["defaultCO2"];
      }
      
      //If DOY == 1 reset PLC (Growth assumed)
      if(cavitationRefill=="annual") {
        if(DOY[i]==1) {
          DataFrame internalWater = Rcpp::as<Rcpp::DataFrame>(x["internalWater"]);
          NumericVector StemPLC = Rcpp::as<Rcpp::NumericVector>(internalWater["StemPLC"]);
          for(int j=0;j<StemPLC.length();j++) StemPLC[j] = 0.0;
        }
      }

      if(unlimitedSoilWater) {
        NumericVector W = soil["W"];
        for(int h=0;h<W.size();h++) W[h] = 1.0;
      }
      
      
      //Julian day from either input column or date
      int J = NA_INTEGER;
      if(julianday_input) J = JulianDay[i];
      if(IntegerVector::is_na(J)){
        std::string c = as<std::string>(dateStrings[i]);
        J = meteoland::radiation_julianDay(std::atoi(c.substr(0, 4).c_str()),std::atoi(c.substr(5,2).c_str()),std::atoi(c.substr(8,2).c_str())); 
      }
      double delta = meteoland::radiation_solarDeclination(J);
      double solarConstant = meteoland::radiation_solarConstant(J);

      double tmin = MinTemperature[i];
      double tmax = MaxTemperature[i];
      double tday = meteoland::utils_averageDaylightTemperature(tmin, tmax);
      double rhmin = MinRelativeHumidity[i];
      double rhmax = MaxRelativeHumidity[i];
      double prec = Precipitation[i];
      double rad = Radiation[i];
      if(tmin > tmax) {
        warning("tmin > tmax. Swapping values.");
        double swap = tmin;
        tmin = tmax;
        tmax = swap;
      }
      if(NumericVector::is_na(rhmax)) {
        rhmax = 100.0;
      }
      if(NumericVector::is_na(rhmin)) {
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
      if(NumericVector::is_na(rad)) {
        double vpa = meteoland::utils_averageDailyVP(tmin, tmax, rhmin, rhmax);
        rad = meteoland::radiation_solarRadiation(solarConstant, latrad, elevation,
                                                  slorad, asprad, delta, tmax -tmin, tmax-tmin,
                                                  vpa, prec);
      }
      PET[i] = meteoland::penman(latrad, elevation, slorad, asprad, J, 
                                 tmin, tmax, rhmin, rhmax, rad, wind);
      
      //1. Phenology and leaf fall
      if(leafPhenology) {
        updatePhenology(x, DOY[i], Photoperiod[i], tday);
        updateLeaves(x, wind, false);
      }

      //2. Water balance and photosynthesis
      if(transpirationMode=="Granier") {
        NumericVector meteovec = NumericVector::create(
          Named("tday") = tday, Named("tmax") = tmax, Named("tmin") = tmin,
          Named("prec") = prec, Named("rhmin") = rhmin, Named("rhmax") = rhmax,
          Named("rad") = rad, 
          Named("wind") = wind, 
          Named("Catm") = Catm,
          Named("Patm") = Patm[i],
          Named("pet") = PET[i],
          Named("er") = erFactor(DOY[i], PET[i], prec));
        try{
          s = spwbDay_basic(x, meteovec, 
                            elevation, slope, aspect, 
                            0.0, verbose); //No Runon in simulations for a single cell
        } catch(std::exception& ex) {
          Rcerr<< "c++ error: "<< ex.what() <<"\n";
          error_occurence = true;
        }
      } else {
        // int ntimesteps = control["ndailysteps"];
        double tmaxPrev = tmax;
        double tminPrev = tmin;
        double tminNext = tmin;
        if(i>0) {
          tmaxPrev = MaxTemperature[i-1];
          tminPrev = MinTemperature[i-1];
        }
        if(i<(numDays-1)) tminNext = MinTemperature[i+1]; 
        NumericVector meteovec = NumericVector::create(
          Named("tmin") = tmin, 
          Named("tmax") = tmax,
          Named("tminPrev") = tminPrev, 
          Named("tmaxPrev") = tmaxPrev, 
          Named("tminNext") = tminNext, 
          Named("prec") = Precipitation[i],
          Named("rhmin") = rhmin, 
          Named("rhmax") = rhmax, 
          Named("rad") = rad, 
          Named("wind") = wind, 
          Named("Catm") = Catm,
          Named("Patm") = Patm[i],
          Named("pet") = PET[i],
          Named("er") = erFactor(DOY[i], PET[i], Precipitation[i]));
          try{
            s = spwbDay_advanced(x, meteovec, 
                         latitude, elevation, slope, aspect,
                         solarConstant, delta, 
                         0.0, verbose); 
          } catch(std::exception& ex) {
            Rcerr<< "c++ error: "<< ex.what() <<"\n";
            error_occurence = true;
          }
        
        fillEnergyBalanceTemperatureDailyOutput(DEB,DT,DLT, s,i, multiLayerBalance);
      }
      
      //Update plant daily water output
      fillPlantWaterDailyOutput(plantDWOL, sunlitDO, shadeDO, s, i, transpirationMode);
      fillWaterBalanceDailyOutput(DWB, s,i, transpirationMode);
      fillSoilWaterBalanceDailyOutput(SWB, soil, s,
                                      i, numDays, transpirationMode, soilFunctions);
      if(control["fireHazardResults"]) fillFireHazardOutput(fireHazard, s, i);
      
      List stand = s["Stand"];
      LgroundPAR[i] = stand["LgroundPAR"];
      LgroundSWR[i] = stand["LgroundSWR"];
      LAI[i] = stand["LAI"];
      LAIherb[i] = stand["LAIherb"];
      LAIexpanded[i] = stand["LAIexpanded"];
      LAIlive[i] = stand["LAIlive"];
      LAIdead[i] = stand["LAIdead"];
      Cm[i] = stand["Cm"];
      

      if(control["subdailyResults"]) {
        subdailyRes[i] = clone(s);
      }
  }
  if(verbose) Rcout << "\n\n";
  
  if(verbose) {
    printWaterBalanceResult(DWB, plantDWOL, x,
                            initialPlantContent, initialSoilContent, initialSnowContent,
                            transpirationMode);
    if(error_occurence) {
      Rcout<< " ERROR: Calculations stopped because of numerical error: Revise parameters\n";
    }
  }


  subdailyRes.attr("names") = dateStrings ;
  
  NumericVector topo = NumericVector::create(elevation, slope, aspect);
  topo.attr("names") = CharacterVector::create("elevation", "slope", "aspect");
  
  DataFrame Stand = DataFrame::create(_["LAI"]=LAI, _["LAIherb"]=LAIherb, 
                                      _["LAIlive"]=LAIlive, _["LAIexpanded"] = LAIexpanded, _["LAIdead"] = LAIdead,  
                                      _["Cm"]=Cm, _["LgroundPAR"] = LgroundPAR, _["LgroundSWR"] = LgroundSWR);
  Stand.attr("row.names") = dateStrings;
  
  List l;
  if(transpirationMode=="Granier") {
    l = List::create(Named("latitude") = latitude,
                     Named("topography") = topo,
                     Named("weather") = clone(meteo),
                     Named("spwbInput") = spwbInput,
                     Named("spwbOutput") = clone(x),
                     Named("WaterBalance")=DWB);
    if(control["soilResults"]) l.push_back(SWB, "Soil");
    if(control["standResults"]) l.push_back(Stand, "Stand");
    if(control["plantResults"]) l.push_back(plantDWOL, "Plants");
    if(control["fireHazardResults"]) l.push_back(fireHazard, "FireHazard");
  } else {
    l = List::create(Named("latitude") = latitude,
                     Named("topography") = topo,
                     Named("weather") = clone(meteo),
                     Named("spwbInput") = spwbInput,
                     Named("spwbOutput") = clone(x),
                     Named("WaterBalance")=DWB, 
                     Named("EnergyBalance") = DEB);
    if(control["temperatureResults"]) {
      l.push_back(DT, "Temperature");
      if(multiLayerBalance) l.push_back(DLT,"TemperatureLayers");
    }
    if(control["soilResults"]) l.push_back(SWB, "Soil");
    if(control["standResults"]) l.push_back(Stand, "Stand");
    if(control["plantResults"]) l.push_back(plantDWOL, "Plants");
    if(control["leafResults"]) {
      l.push_back(sunlitDO, "SunlitLeaves");
      l.push_back(shadeDO, "ShadeLeaves");
    }
    if(control["fireHazardResults"]) l.push_back(fireHazard, "FireHazard");
  }
  if(control["subdailyResults"]) l.push_back(subdailyRes,"subdaily");
  l.attr("class") = CharacterVector::create("spwb","list");
  return(l);
}


//' @rdname spwb
//' 
//' @param W A matrix with the same number of rows as \code{meteo} and as many columns as soil layers, containing the soil moisture of each layer as proportion of field capacity.
//' @param canopyEvaporation A vector of daily canopy evaporation (from interception) values (mm). The length should match the number of rows in \code{meteo}.
//' @param snowMelt A vector of daily snow melt values (mm). The length should match the number of rows in \code{meteo}.
//' @param soilEvaporation A vector of daily bare soil evaporation values (mm). The length should match the number of rows in \code{meteo}.
//' @param herbTranspiration A vector of daily herbaceous transpiration values (mm). The length should match the number of rows in \code{meteo}.
//' 
// [[Rcpp::export("pwb")]]
List pwb(List x, DataFrame meteo, NumericMatrix W,
         double latitude, double elevation = NA_REAL, double slope = NA_REAL, double aspect = NA_REAL, 
         NumericVector canopyEvaporation = NumericVector(0), 
         NumericVector snowMelt = NumericVector(0), 
         NumericVector soilEvaporation = NumericVector(0),
         NumericVector herbTranspiration = NumericVector(0),
         NumericVector CO2ByYear = NumericVector(0)) {
  List control = x["control"];
  String transpirationMode = control["transpirationMode"];
  String soilFunctions = control["soilFunctions"];
  String cavitationRefill = control["cavitationRefill"];
  bool verbose = control["verbose"];
  bool subdailyResults = control["subdailyResults"];
  bool leafPhenology = control["leafPhenology"];
  bool multiLayerBalance = control["multiLayerBalance"];
  
  //Store input
  List spwbInput = x; // Store initial object
  x = clone(x); //Ensure a copy will be modified
  
  
  List soil = x["soil"];
  
  if(NumericVector::is_na(latitude)) stop("Value for 'latitude' should not be missing.");
  double latrad = latitude * (M_PI/180.0);
  if(NumericVector::is_na(aspect)) aspect = 0.0;
  if(NumericVector::is_na(slope)) slope = 0.0;
  double asprad = aspect * (M_PI/180.0);
  double slorad = slope * (M_PI/180.0);
  
    //Meteorological input    
  NumericVector MinTemperature, MaxTemperature;
  NumericVector MinRelativeHumidity, MaxRelativeHumidity;
  NumericVector Radiation;
  if(!meteo.containsElementNamed("Precipitation")) stop("Please include variable 'Precipitation' in weather input.");
  NumericVector Precipitation = meteo["Precipitation"];
  int numDays = Precipitation.size();
  if(NumericVector::is_na(elevation)) stop("Value for 'elevation' should not be missing.");
  if(!meteo.containsElementNamed("MinTemperature")) stop("Please include variable 'MinTemperature' in weather input.");
  MinTemperature = meteo["MinTemperature"];
  if(!meteo.containsElementNamed("MaxTemperature")) stop("Please include variable 'MaxTemperature' in weather input.");
  MaxTemperature = meteo["MaxTemperature"];
  if(!meteo.containsElementNamed("MinRelativeHumidity")) stop("Please include variable 'MinRelativeHumidity' in weather input.");
  MinRelativeHumidity = meteo["MinRelativeHumidity"];
  if(!meteo.containsElementNamed("MaxRelativeHumidity")) stop("Please include variable 'MaxRelativeHumidity' in weather input.");
  MaxRelativeHumidity = meteo["MaxRelativeHumidity"];
  if(!meteo.containsElementNamed("Radiation")) stop("Please include variable 'Radiation' in weather input.");
  Radiation = meteo["Radiation"];
  
  if(any(is_na(Precipitation))) stop("Missing values in 'Precipitation'");
  if(any(is_na(MinTemperature))) stop("Missing values in 'MinTemperature'");
  if(any(is_na(MaxTemperature))) stop("Missing values in 'MaxTemperature'");
  if(any(is_na(MinRelativeHumidity))) warning("Missing values in 'MinRelativeHumidity' were estimated from temperature range");
  if(any(is_na(MaxRelativeHumidity))) warning("Missing values in 'MaxRelativeHumidity' were assumed to be 100");
  if(any(is_na(Radiation))) warning("Missing values in 'Radiation' were estimated");
  
  NumericVector WindSpeed(numDays, NA_REAL);
  if(meteo.containsElementNamed("WindSpeed")) WindSpeed = meteo["WindSpeed"];
  
  NumericVector PET = NumericVector(numDays,NA_REAL);
  
  NumericVector CO2(Precipitation.length(), NA_REAL);
  if(meteo.containsElementNamed("CO2")) {
    CO2 = meteo["CO2"];
    if(verbose) {
      Rcout<<"CO2 taken from input column 'CO2'\n";
    }
    if(any(is_na(CO2))) stop("Missing values in 'CO2'");
  }
  NumericVector Patm(Precipitation.length(), NA_REAL);
  if(meteo.containsElementNamed("Patm")) {
    Patm = meteo["Patm"];
    if(verbose) {
      Rcout<<"Patm taken from input column 'Patm'\n";
    }
  }
  IntegerVector DOY, JulianDay;
  NumericVector Photoperiod;
  bool doy_input = false, photoperiod_input = false, julianday_input = false;
  if(meteo.containsElementNamed("DOY")) {
    DOY = meteo["DOY"];
    doy_input = true;
    if(verbose) {
      Rcout<<"DOY taken from input column 'DOY'\n";
    }
  }
  if(meteo.containsElementNamed("Photoperiod")) {
    Photoperiod = meteo["Photoperiod"];
    photoperiod_input = true;
    if(verbose) {
      Rcout<<"Photoperiod taken from input column 'Photoperiod'\n";
    }
  }
  if(meteo.containsElementNamed("JulianDay")) {
    JulianDay = meteo["JulianDay"];
    julianday_input = true;
    if(verbose) {
      Rcout<<"Julian day taken from input column 'JulianDay'\n";
    }
  }
  

  if(canopyEvaporation.length()==0) {
    canopyEvaporation = NumericVector(numDays,0.0);
  }
  if(snowMelt.length()==0) {
    snowMelt = NumericVector(numDays,0.0);
  }
  if(soilEvaporation.length()==0) {
    soilEvaporation = NumericVector(numDays,0.0);
  }
  
  // Dates
  CharacterVector dateStrings = getWeatherDates(meteo);
  if(!doy_input) DOY = date2doy(dateStrings);
  if(!photoperiod_input) Photoperiod = date2photoperiod(dateStrings, latrad);
  
  
  //Canopy scalars
  DataFrame canopy = Rcpp::as<Rcpp::DataFrame>(x["canopy"]);
  
  //Plant input
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  int numCohorts = above.nrow();
  
  
  //Soil input
  NumericVector Water_FC = waterFC(soil, soilFunctions);
  int nlayers = Water_FC.size();
  
  //Detailed subday results
  List subdailyRes(numDays);
  
  //Transpiration output variables
  NumericVector Transpiration(numDays);
  NumericVector PlantExtraction(numDays);
  NumericVector HydraulicRedistribution(numDays, 0.0);
  NumericMatrix HydrIndays(numDays, nlayers);
  
  //EnergyBalance output variables
  DataFrame DEB = defineEnergyBalanceDailyOutput(dateStrings);
  DataFrame DT = defineTemperatureDailyOutput(dateStrings);
  NumericMatrix DLT;
  if(transpirationMode=="Sperry") DLT =  defineTemperatureLayersDailyOutput(dateStrings, canopy);
  
  //Stand output variables
  NumericVector LAI(numDays),LAIherb(numDays), LAIlive(numDays),LAIexpanded(numDays),LAIdead(numDays);


  //Soil output variables
  NumericMatrix Wdays(numDays, nlayers); //Soil moisture content in relation to field capacity
  NumericMatrix psidays(numDays, nlayers);
  NumericMatrix Eplantdays(numDays, nlayers);
  
  //Plant output variables
  List sunlitDO = defineSunlitShadeLeavesDailyOutput(dateStrings, above);
  List shadeDO = defineSunlitShadeLeavesDailyOutput(dateStrings, above);
  List plantDWOL = definePlantWaterDailyOutput(dateStrings, above, soil, control);
  NumericVector EplantCohTot(numCohorts, 0.0);

  //Fire hazard output variables
  DataFrame fireHazard;
  if(control["fireHazardResults"]) fireHazard = defineFireHazardOutput(dateStrings);
  
  bool error_occurence = false;
  if(verbose) Rcout << "Performing daily simulations ";
  NumericVector Eplanttot(numDays,0.0);
  List s;
  std::string yearString;
  for(int i=0;i<numDays;i++) {
    std::string c = as<std::string>(dateStrings[i]);
    yearString = c.substr(0, 4);
    if(verbose) {
      if(DOY[i]==1 || i==0) {
        Rcout<<"\n Year "<< yearString<< ":";
      } 
      else if(i%10 == 0) Rcout<<".";//<<i;
    } 
    double wind = WindSpeed[i];
    if(NumericVector::is_na(wind)) wind = control["defaultWindSpeed"]; //Default 1 m/s -> 10% of fall every day
    if(wind<0.1) wind = 0.1; //Minimum windspeed abovecanopy
    
    
    double Catm = CO2[i];
    //If missing, use
    if(NumericVector::is_na(Catm)) {
      if(CO2ByYear.attr("names") != R_NilValue) Catm = CO2ByYear[yearString];
    }
    //If still missing, use default control value
    if(NumericVector::is_na(Catm)) {
      Catm = control["defaultCO2"];
    }
    
    
    //Julian day from either input column or date
    int J = NA_INTEGER;
    if(julianday_input) J = JulianDay[i];
    if(IntegerVector::is_na(J)){
      std::string c = as<std::string>(dateStrings[i]);
      J = meteoland::radiation_julianDay(std::atoi(c.substr(0, 4).c_str()),std::atoi(c.substr(5,2).c_str()),std::atoi(c.substr(8,2).c_str())); 
    }
    double delta = meteoland::radiation_solarDeclination(J);
    double solarConstant = meteoland::radiation_solarConstant(J);
    
    double tmin = MinTemperature[i];
    double tmax = MaxTemperature[i];
    double tday = meteoland::utils_averageDaylightTemperature(tmin, tmax);
    double rhmin = MinRelativeHumidity[i];
    double rhmax = MaxRelativeHumidity[i];
    double rad = Radiation[i];
    double prec = Precipitation[i];
    if(tmin > tmax) {
      warning("tmin > tmax. Swapping values.");
      double swap = tmin;
      tmin = tmax;
      tmax = swap;
    }
    if(NumericVector::is_na(rhmax)) {
      rhmax = 100.0;
    }
    if(NumericVector::is_na(rhmin)) {
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
    if(NumericVector::is_na(rad)) {
      double vpa = meteoland::utils_averageDailyVP(tmin, tmax, rhmin, rhmax);
      rad = meteoland::radiation_solarRadiation(solarConstant, latrad, elevation,
                                                slorad, asprad, delta, tmax -tmin, tmax-tmin,
                                                vpa, prec);
    }
    PET[i] = meteoland::penman(latrad, elevation, slorad, asprad, J, 
                               tmin, tmax, rhmin, rhmax, rad, wind);

    //0. Soil moisture
    soil["W"] = W(i,_);
    Wdays(i,_) = W(i,_);
    psidays(i,_) = psi(soil, soilFunctions); //Get soil water potential
      
    //If DOY == 1 reset PLC (Growth assumed)
    if(cavitationRefill=="annual") {
        if(DOY[i]==1) {
          DataFrame internalWater = Rcpp::as<Rcpp::DataFrame>(x["internalWater"]);
          if(transpirationMode=="Granier") {
            NumericVector PLC = Rcpp::as<Rcpp::NumericVector>(internalWater["PLC"]);
            for(int j=0;j<PLC.length();j++) PLC[j] = 0.0;
          } else {
            NumericVector StemPLC = Rcpp::as<Rcpp::NumericVector>(internalWater["StemPLC"]);
            for(int j=0;j<StemPLC.length();j++) StemPLC[j] = 0.0;
          }
        }
    }
      
      
    //1. Phenology and leaf fall
    if(leafPhenology) {
      updatePhenology(x, DOY[i], Photoperiod[i], tday);
      updateLeaves(x, wind, false);
    }
      
    
    int ntimesteps = control["ndailysteps"];
   
    //2. transpiration and photosynthesis
    if(transpirationMode=="Granier") {
      NumericVector meteovec = NumericVector::create(
        Named("tday") = tday, Named("tmax") = tmax, Named("tmin") = tmin,
        Named("prec") = prec, Named("rhmin") = rhmin, Named("rhmax") = rhmax,
        Named("rad") = rad, 
        Named("wind") = wind, 
        Named("Catm") = Catm,
        Named("Patm") = Patm[i],
        Named("pet") = PET[i]);
      try{
        s = transpirationBasic(x, meteovec, 
                               elevation, true);
      } catch(std::exception& ex) {
        Rcerr<< "c++ error: "<< ex.what() <<"\n";
        error_occurence = true;
      }
    } else {
      double tmaxPrev = tmax;
      double tminPrev = tmin;
      double tminNext = tmin;
      if(i>0) {
        tmaxPrev = MaxTemperature[i-1];
        tminPrev = MinTemperature[i-1];
      }
      if(i<(numDays-1)) tminNext = MinTemperature[i+1]; 
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
        Named("Patm") = Patm[i]);
      try{
        s = transpirationAdvanced(x, meteovec, 
                                latitude, elevation, slope, aspect,
                                solarConstant, delta,
                                canopyEvaporation[i], snowMelt[i], soilEvaporation[i], herbTranspiration[i],
                                verbose, NA_INTEGER, 
                                true);
      } catch(std::exception& ex) {
        Rcerr<< "c++ error: "<< ex.what() <<"\n";
        error_occurence = true;
      }
      fillEnergyBalanceTemperatureDailyOutput(DEB,DT,DLT,s,i, multiLayerBalance);
    }
    
    //Update plant daily water output
    fillPlantWaterDailyOutput(plantDWOL, sunlitDO, shadeDO, s, i, transpirationMode);
    if(control["fireHazardResults"]) fillFireHazardOutput(fireHazard, s, i);
    
    List Plants = s["Plants"];
    NumericVector EplantCoh = Plants["Transpiration"];
    NumericMatrix SoilWaterExtract = s["Extraction"];
    for(int l=0;l<nlayers;l++) {
      Eplantdays(i,l) = sum(SoilWaterExtract(_,l));
    }

    PlantExtraction[i] = sum(SoilWaterExtract);
    Transpiration[i] = sum(EplantCoh);
    NumericVector HydrInVec(nlayers, 0.0);

    if(transpirationMode=="Sperry")  {
      NumericMatrix soilLayerExtractInst = s["ExtractionInst"];
      for(int l=0;l<nlayers;l++) {
        for(int n=0;n<ntimesteps;n++) {
          HydrInVec[l] += (-1.0)*std::min(soilLayerExtractInst(l,n),0.0);
        }
      }
      HydraulicRedistribution[i] = sum(HydrInVec);
      HydrIndays(i,_) = HydrInVec;
    } 
    List stand = s["Stand"];
    LAI[i] = stand["LAI"];
    LAIherb[i] = stand["LAIherb"];
    LAIlive[i] = stand["LAIlive"];
    LAIexpanded[i] = stand["LAIexpanded"];
    LAIdead[i] = stand["LAIdead"];
        
    EplantCohTot = EplantCohTot + EplantCoh;
    Eplanttot[i] = sum(EplantCoh);
    
    if(subdailyResults) {
      subdailyRes[i] = clone(s);
    }
  }
  if(verbose) Rcout << "done\n";
  
  if(verbose) {
    double Transpirationsum = sum(Transpiration);
    
    Rcout<<"Transpiration (mm) "  <<round(Transpirationsum);
    Rcout<<" Plant extraction from soil (mm) " << round(sum(PlantExtraction));
    Rcout<<" Hydraulic redistribution (mm) " << round(sum(HydraulicRedistribution)) <<"\n";
    if(error_occurence) {
      Rcout<< " ERROR: Calculations stopped because of numerical error: Revise parameters\n";
    }
  }

  
  DataFrame SWB;
  if(transpirationMode=="Granier") {
    SWB = DataFrame::create(_["W"]=Wdays, _["PlantExt"]=Eplantdays, _["psi"]=psidays); 
  } else {
    SWB = DataFrame::create(_["W"]=Wdays, _["PlantExt"]=Eplantdays,
                            _["HydraulicInput"] = HydrIndays,
                            _["psi"]=psidays); 
  }
  SWB.attr("row.names") = dateStrings;
  DataFrame Stand = DataFrame::create(_["LAI"]=LAI,_["LAIherb"]=LAIherb, 
                                      _["LAIlive"]=LAIlive,_["LAIexpanded"]=LAIexpanded, _["LAIdead"] = LAIdead);
  Stand.attr("row.names") = dateStrings;
  
  DataFrame DWB = DataFrame::create(_["PlantExtraction"] = PlantExtraction, _["Transpiration"]=Transpiration, 
                                    _["HydraulicRedistribution"] = HydraulicRedistribution);
  DWB.attr("row.names") = dateStrings;
  
  subdailyRes.attr("names") = dateStrings;
  
  NumericVector topo = NumericVector::create(elevation, slope, aspect);
  topo.attr("names") = CharacterVector::create("elevation", "slope", "aspect");

  List l;
  if(transpirationMode=="Granier") {
    l = List::create(Named("latitude") = latitude,
                     Named("topography") = topo,
                     Named("weather") = clone(meteo),
                     Named("spwbInput") = spwbInput,
                     Named("spwbOutput") = clone(x),
                     Named("WaterBalance")=DWB);
    if(control["soilResults"]) l.push_back(SWB, "Soil");
    if(control["standResults"]) l.push_back(Stand, "Stand");
    if(control["plantResults"]) l.push_back(plantDWOL, "Plants");
    if(control["fireHazardResults"]) l.push_back(fireHazard, "FireHazard");
  } else {
    l = List::create(Named("latitude") = latitude,
                     Named("topography") = topo,
                     Named("weather") = clone(meteo),
                     Named("spwbInput") = spwbInput,
                     Named("spwbOutput") = clone(x),
                     Named("WaterBalance")=DWB, 
                     Named("EnergyBalance") = DEB);
    if(control["temperatureResults"]) {
      l.push_back(DT, "Temperature");
      if(multiLayerBalance) l.push_back(DLT,"TemperatureLayers");
    }
    if(control["soilResults"]) l.push_back(SWB, "Soil");
    if(control["standResults"]) l.push_back(Stand, "Stand");
    if(control["plantResults"]) l.push_back(plantDWOL, "Plants");
    if(control["leafResults"]) {
      l.push_back(sunlitDO, "SunlitLeaves");
      l.push_back(shadeDO, "ShadeLeaves");
    }
    if(control["fireHazardResults"]) l.push_back(fireHazard, "FireHazard");
  }
  if(control["subdailyResults"]) l.push_back(subdailyRes,"subdaily");
  l.attr("class") = CharacterVector::create("pwb","list");
  return(l);                    
}
