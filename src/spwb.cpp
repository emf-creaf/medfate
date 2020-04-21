// [[Rcpp::interfaces(r,cpp)]]

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
#include "soil.h"
#include <Rcpp.h>
#include <meteoland.h>
using namespace Rcpp;


// Soil water balance with simple hydraulic model
List spwbDay1(List x, List soil, double tday, double pet, double prec, double er, double runon=0.0, 
             double rad = NA_REAL, double elevation = NA_REAL, bool verbose = false) {

  //Control parameters
  List control = x["control"];
  bool snowpack = control["snowpack"];
  bool drainage = control["drainage"];
  bool plantWaterPools = control["plantWaterPools"];
  String soilFunctions = control["soilFunctions"];

  //Number of soil layers
  int nlayers = Rcpp::as<Rcpp::NumericVector>(soil["dVec"]).size();
  
  List below = x["below"];
  NumericMatrix Wpool = below["Wpool"];
  NumericVector Wsoil = soil["W"];

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
  
  //Determine whether leaves are out (phenology) and the adjusted Leaf area
  NumericVector Phe(numCohorts,0.0);
  double s = 0.0, LAIcell = 0.0, LAIcelllive, LAIcelldead = 0.0, Cm = 0.0;
  for(int c=0;c<numCohorts;c++) {
    if(LAIlive[c]>0) Phe[c]=LAIphe[c]/LAIlive[c]; //Phenological status
    else Phe[c]=0.0;
    s += (kPAR[c]*(LAIphe[c]+LAIdead[c]));
    LAIcell += LAIphe[c]+LAIdead[c];
    LAIcelldead += LAIdead[c];
    LAIcelllive += LAIlive[c];
    Cm += (LAIphe[c]+LAIdead[c])*gRainIntercept[c]; //LAI dead also counts on interception
  }
  double LgroundPAR = exp((-1.0)*s);
  double LgroundSWR = exp((-1.0)*s/1.35);
  
  //Snow pack dynamics and hydrology input
  NumericVector hydroInputs = soilWaterInputs(soil, soilFunctions, prec, er, tday, rad, elevation,
                                             Cm, LgroundPAR, LgroundSWR, 
                                             runon,
                                             snowpack, true);
  
  NumericVector infilPerc, EsoilVec;
  NumericVector EplantVec(nlayers, 0.0);
  
  if(!plantWaterPools) {
    //Soil infiltration and percolation
    infilPerc = soilInfiltrationPercolation(soil, soilFunctions, 
                                            hydroInputs["Input"],
                                            drainage, true);
    //Evaporation from bare soil (if there is no snow)
    EsoilVec = soilEvaporation(soil, soilFunctions, pet, LgroundSWR, true);
    
    //Copy soil status to x
    for(int c=0;c<numCohorts;c++) for(int l=0;l<nlayers;l++) Wpool(c,l) = Wsoil[l];
  } else {
    //Reset soil moisture
    for(int l=0;l<nlayers;l++) Wsoil[l] = 0.0;
    
    //Initialize result vectors
    infilPerc = NumericVector::create(_["Infiltration"] = 0.0, 
                                      _["Runoff"] = 0.0, 
                                      _["DeepDrainage"] = 0.0);
    EsoilVec = NumericVector(nlayers,0.0);
    for(int c=0;c<numCohorts;c++) {
      double f_soil_c = LAIlive[c]/LAIcelllive;
      
      //Clone soil and copy moisture values from x
      List soil_c =  clone(soil);
      NumericVector W_c = soil_c["W"];
      for(int l=0;l<nlayers;l++) W_c[l] = Wpool(c,l);
      
      //Soil_c infiltration and percolation
      NumericVector infilPerc_c = soilInfiltrationPercolation(soil_c, soilFunctions, 
                                              hydroInputs["Input"],
                                              drainage, true);
      //Evaporation from bare soil_c (if there is no snow)
      NumericVector EsoilVec_c = soilEvaporation(soil_c, soilFunctions, pet, LgroundSWR, true);
      //Copy result vectors
      infilPerc["Infiltration"] = infilPerc["Infiltration"] + f_soil_c*infilPerc_c["Infiltration"];
      infilPerc["Runoff"] = infilPerc["Runoff"] + f_soil_c*infilPerc_c["Runoff"];
      infilPerc["DeepDrainage"] = infilPerc["DeepDrainage"] + f_soil_c*infilPerc_c["DeepDrainage"];
      for(int l=0;l<nlayers;l++) EsoilVec[l] = EsoilVec[l] + f_soil_c*EsoilVec_c[l];
      // Copy soil_c status back to x
      for(int l=0;l<nlayers;l++) {
        Wpool(c,l) = W_c[l];
        Wsoil[l] = Wsoil[l] + f_soil_c*Wpool(c,l); //weighted average for soil moisture
      }
    }
  }
  
  //Canopy transpiration  
  // Rcout<<"hola";
  List transp = transpirationGranier(x, soil, tday, pet, true, true);
  // Rcout<<"hola2";
  NumericMatrix EplantCoh = Rcpp::as<Rcpp::NumericMatrix>(transp["Extraction"]);
  for(int l=0;l<nlayers;l++) EplantVec[l] = sum(EplantCoh(_,l));
  DataFrame Plants = Rcpp::as<Rcpp::DataFrame>(transp["Plants"]);

  NumericVector psiVec = psi(soil, soilFunctions); //Calculate current soil water potential for output
  
  NumericVector DB = NumericVector::create(_["PET"] = pet, 
                                           _["Rain"] = hydroInputs["Rain"], _["Snow"] = hydroInputs["Snow"], 
                                           _["NetRain"] = hydroInputs["NetRain"], _["Snowmelt"] = hydroInputs["Snowmelt"],
                                           _["Runon"] = hydroInputs["Runon"], 
                                           _["Infiltration"] = infilPerc["Infiltration"], _["Runoff"] = infilPerc["Runoff"], _["DeepDrainage"] = infilPerc["DeepDrainage"],
                                           _["SoilEvaporation"] = sum(EsoilVec), _["PlantExtraction"] = sum(EplantVec), _["Transpiration"] = sum(EplantVec));
  
  NumericVector Stand = NumericVector::create(_["LAIcell"] = LAIcell, _["LAIcelldead"] = LAIcelldead, 
                                           _["Cm"] = Cm, _["LgroundPAR"] = LgroundPAR, _["LgroundSWR"] = LgroundSWR);
  DataFrame SB = DataFrame::create(_["SoilEvaporation"] = EsoilVec, 
                                   _["PlantExtraction"] = EplantVec, 
                                   _["psi"] = psiVec);
  List l = List::create(_["cohorts"] = clone(cohorts),
                        _["WaterBalance"] = DB, 
                        _["Soil"] = SB,
                        _["Stand"] = Stand,
                        _["Plants"] = Plants);
  l.attr("class") = CharacterVector::create("spwb_day","list");
  return(l);
}



// Soil water balance with Sperry hydraulic and stomatal conductance models
List spwbDay2(List x, List soil, double tmin, double tmax, double tminPrev, double tmaxPrev, double tminNext, double rhmin, double rhmax, double rad, double wind, 
             double latitude, double elevation, double slope, double aspect,
             double solarConstant, double delta, 
             double prec, double pet, double er, double runon=0.0, bool verbose = false) {
  
  //Control parameters
  List control = x["control"];
  bool drainage = control["drainage"];
  bool snowpack = control["snowpack"];
  bool plantWaterPools = control["plantWaterPools"];
  String soilFunctions = control["soilFunctions"];
  int ntimesteps = control["ndailysteps"];

  //Number of soil layers
  int nlayers = Rcpp::as<Rcpp::NumericVector>(soil["dVec"]).size();

  List below = x["below"];
  NumericMatrix Wpool = below["Wpool"];
  NumericVector Wsoil = soil["W"];
  
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
  NumericVector Phe(numCohorts);
  double s = 0.0, LAIcell = 0.0, LAIcelldead = 0.0, LAIcelllive = 0.0, Cm = 0.0, LAIcellmax = 0.0;
  for(int c=0;c<numCohorts;c++) {
    Phe[c]=LAIphe[c]/LAIlive[c]; //Phenological status
    LAIcell += (LAIphe[c]+LAIdead[c]);
    LAIcelldead += LAIdead[c];
    LAIcellmax += LAIlive[c];
    LAIcelllive += LAIlive[c];
    s += (kPAR[c]*(LAIphe[c]+LAIdead[c]));
    Cm += (LAIphe[c]+LAIdead[c])*gRainIntercept[c]; //LAI dead also counts on interception
  }
  double LgroundPAR = exp((-1.0)*s);
  double LgroundSWR = exp((-1.0)*s/1.35);
  
  //A.1 - Snow pack dynamics and soil water input
  NumericVector hydroInputs = soilWaterInputs(soil, soilFunctions, prec, er, tday, rad, elevation,
                                              Cm, LgroundPAR, LgroundSWR, 
                                              runon,
                                              snowpack, true);
  
  NumericVector infilPerc, EsoilVec;
  if(!plantWaterPools) {
    //A.2 - Soil infiltration and percolation
    infilPerc = soilInfiltrationPercolation(soil, soilFunctions, 
                                            hydroInputs["Input"],
                                                       drainage, true);
    //B.1 - Evaporation from bare soil if there is no snow
    EsoilVec = soilEvaporation(soil, soilFunctions, pet, LgroundSWR, true);
    
    //Copy soil status to x
    for(int c=0;c<numCohorts;c++) for(int l=0;l<nlayers;l++) Wpool(c,l) = Wsoil[l];
  } else {
    //Reset soil moisture
    for(int l=0;l<nlayers;l++) Wsoil[l] = 0.0;
    
    //Initialize result vectors
    infilPerc = NumericVector::create(_["Infiltration"] = 0.0, 
                                      _["Runoff"] = 0.0, 
                                      _["DeepDrainage"] = 0.0);
    EsoilVec = NumericVector(nlayers,0.0);
    for(int c=0;c<numCohorts;c++) {
      double f_soil_c = LAIlive[c]/LAIcelllive;
      
      //Clone soil and copy moisture values from x
      List soil_c =  clone(soil);
      NumericVector W_c = soil_c["W"];
      for(int l=0;l<nlayers;l++) W_c[l] = Wpool(c,l);
      
      //Soil_c infiltration and percolation
      NumericVector infilPerc_c = soilInfiltrationPercolation(soil_c, soilFunctions, 
                                                              hydroInputs["Input"],
                                                                         drainage, true);
      //Evaporation from bare soil_c (if there is no snow)
      NumericVector EsoilVec_c = soilEvaporation(soil_c, soilFunctions, pet, LgroundSWR, true);
      //Copy result vectors
      infilPerc["Infiltration"] = infilPerc["Infiltration"] + f_soil_c*infilPerc_c["Infiltration"];
      infilPerc["Runoff"] = infilPerc["Runoff"] + f_soil_c*infilPerc_c["Runoff"];
      infilPerc["DeepDrainage"] = infilPerc["DeepDrainage"] + f_soil_c*infilPerc_c["DeepDrainage"];
      for(int l=0;l<nlayers;l++) EsoilVec[l] = EsoilVec[l] + f_soil_c*EsoilVec_c[l];
      // Copy soil_c status back to x
      for(int l=0;l<nlayers;l++) {
        Wpool(c,l) = W_c[l];
        Wsoil[l] = Wsoil[l] + f_soil_c*Wpool(c,l); //weighted average for soil moisture
      }
    }
  }

  //B.2 - Canopy transpiration  
  List transp = transpirationSperry(x, soil,tmin, tmax, tminPrev, tmaxPrev, tminNext, 
                                    rhmin, rhmax, rad, wind, 
                                    latitude, elevation, slope, aspect, 
                                    solarConstant, delta, prec, 
                                    hydroInputs["Interception"], hydroInputs["Snowmelt"], sum(EsoilVec),
                                    verbose, NA_INTEGER, true, true);

  
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
  
  NumericVector DB = NumericVector::create(_["PET"] = pet,
                                           _["Rain"] = hydroInputs["Rain"],_["Snow"] = hydroInputs["Snow"],_["NetRain"] = hydroInputs["NetRain"], _["Snowmelt"] = hydroInputs["Snowmelt"],
                                           _["Runon"] = hydroInputs["Runon"], 
                                           _["Infiltration"] = infilPerc["Infiltration"], _["Runoff"] = infilPerc["Runoff"], _["DeepDrainage"] = infilPerc["DeepDrainage"],
                                           _["SoilEvaporation"] = sum(EsoilVec), _["PlantExtraction"] = sum(EplantVec), _["Transpiration"] = sum(Eplant),
                                           _["HydraulicRedistribution"] = sum(soilHydraulicInput));
  
  NumericVector Stand = NumericVector::create(_["LAIcell"] = LAIcell, _["LAIcelldead"] = LAIcelldead, _["Cm"] = Cm, 
                                              _["LgroundPAR"] = LgroundPAR, _["LgroundSWR"] = LgroundSWR);
  
  DataFrame SB = DataFrame::create(_["SoilEvaporation"] = EsoilVec, 
                                   _["HydraulicInput"] = soilHydraulicInput, 
                                   _["HydraulicOutput"] = soilHydraulicOutput, 
                                   _["PlantExtraction"] = EplantVec, 
                                   _["psi"] = psiVec);
  
  List l = List::create(_["cohorts"] = clone(cohorts),
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
                        _["LightExtinction"] = transp["LightExtinction"],
                        _["WindExtinction"] = transp["WindExtinction"]);
  l.attr("class") = CharacterVector::create("spwb_day","list");
  return(l);
}

// [[Rcpp::export("spwb_day")]]
List spwbDay(List x, List soil, CharacterVector date, double tmin, double tmax, double rhmin, double rhmax, double rad, double wind, 
            double latitude, double elevation, double slope, double aspect,  
            double prec, double runon=0.0) {
  //Control parameters
  List control = x["control"];
  bool verbose = control["verbose"];
  bool leafPhenology = control["leafPhenology"];
  String transpirationMode = control["transpirationMode"];
  std::string c = as<std::string>(date[0]);
  int J = meteoland::radiation_julianDay(std::atoi(c.substr(0, 4).c_str()),std::atoi(c.substr(5,2).c_str()),std::atoi(c.substr(8,2).c_str()));
  double delta = meteoland::radiation_solarDeclination(J);
  double solarConstant = meteoland::radiation_solarConstant(J);
  double tday = meteoland::utils_averageDaylightTemperature(tmin, tmax);
  double latrad = latitude * (PI/180.0);
  double asprad = aspect * (PI/180.0);
  double slorad = slope * (PI/180.0);
  double photoperiod = meteoland::radiation_daylength(latrad, 0.0, 0.0, delta);
  double pet = meteoland::penman(latrad, elevation, slorad, asprad, J, tmin, tmax, rhmin, rhmax, rad, wind);

  //Derive doy from date  
  int J0101 = meteoland::radiation_julianDay(std::atoi(c.substr(0, 4).c_str()),1,1);
  int doy = J - J0101+1;
  if(NumericVector::is_na(wind)) wind = control["defaultWindSpeed"]; 
  if(wind<0.1) wind = 0.1; //Minimum windspeed abovecanopy
  
  //Update phenology
  if(leafPhenology) updateLeaves(x, doy, photoperiod, tday, wind);
  
  double er = erFactor(doy, pet, prec);
  List s;
  if(transpirationMode=="Granier") {
    s = spwbDay1(x,soil, tday, pet, prec, er, runon, rad, elevation, verbose);
  } else {
    s = spwbDay2(x,soil, tmin, tmax, tmin, tmax, tmin, rhmin, rhmax, rad, wind, 
                 latitude, elevation, slope, aspect,
                 solarConstant, delta, prec, pet, er, runon, verbose);
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



void checkspwbInput(List x, List soil, String transpirationMode, String soilFunctions) {
  if(!x.containsElementNamed("above")) stop("above missing in spwbInput");
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  if(!above.containsElementNamed("LAI_live")) stop("LAI_live missing in spwbInput$above");
  if(!above.containsElementNamed("CR")) stop("CR missing in spwbInput$above");
  if(!above.containsElementNamed("H")) stop("H missing in spwbInput$above");
  
  if(!x.containsElementNamed("below")) stop("below missing in spwbInput");
  List below = Rcpp::as<Rcpp::List>(x["below"]);
  if(!below.containsElementNamed("V")) stop("V missing in spwbInput$below");
  if(transpirationMode=="Sperry"){
    if(!below.containsElementNamed("VGrhizo_kmax")) stop("VGrhizo_kmax missing in spwbInput$below");
    if(!below.containsElementNamed("VCroot_kmax")) stop("VCroot_kmax missing in spwbInput$below");
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
    if(!paramsTranspiration.containsElementNamed("pRootDisc")) stop("pRootDisc missing in spwbInput$paramsTranspiration");
    if(!paramsTranspiration.containsElementNamed("Psi_Extract")) stop("Psi_Extract missing in spwbInput$paramsTranspiration");
    if(!paramsTranspiration.containsElementNamed("WUE")) stop("WUE missing in spwbInput$paramsTranspiration");
  } else if(transpirationMode=="Sperry") {
    if(!paramsTranspiration.containsElementNamed("VCstem_kmax")) stop("VCstem_kmax missing in spwbInput$paramsTranspiration");
    if(!paramsTranspiration.containsElementNamed("VCstem_c")) stop("VCstem_c missing in spwbInput$paramsTranspiration");
    if(!paramsTranspiration.containsElementNamed("VCstem_d")) stop("VCstem_d missing in spwbInput$paramsTranspiration");
    if(!paramsTranspiration.containsElementNamed("VCroot_c")) stop("VCroot_c missing in spwbInput$paramsTranspiration");
    if(!paramsTranspiration.containsElementNamed("VCroot_d")) stop("VCroot_d missing in spwbInput$paramsTranspiration");
    if(!paramsTranspiration.containsElementNamed("Gwmax")) stop("Gwmax missing in spwbInput$paramsTranspiration");
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


// [[Rcpp::export("spwb")]]
List spwb(List x, List soil, DataFrame meteo, double latitude, double elevation = NA_REAL, double slope = NA_REAL, double aspect = NA_REAL) {
  List control = x["control"];
  String transpirationMode = control["transpirationMode"];
  String soilFunctions = control["soilFunctions"];
  String cavitationRefill = control["cavitationRefill"];
  
  bool verbose = control["verbose"];
  bool subdailyResults = control["subdailyResults"];
  bool leafPhenology = control["leafPhenology"];
  bool unlimitedSoilWater = control["unlimitedSoilWater"];
  checkspwbInput(x, soil, transpirationMode, soilFunctions);
  
  //Store input
  List spwbInput = clone(x);
  List soilInput = clone(soil);
    
  //Meteorological input    
  NumericVector MinTemperature, MaxTemperature;
  NumericVector MinRelativeHumidity, MaxRelativeHumidity;
  NumericVector Radiation;
  
  if(NumericVector::is_na(latitude)) stop("Value for 'latitude' should not be missing.");
  double latrad = latitude * (PI/180.0);
  
  if(!meteo.containsElementNamed("Precipitation")) stop("Please include variable 'Precipitation' in weather input.");
  NumericVector Precipitation = meteo["Precipitation"];
  int numDays = Precipitation.size();
  if(!meteo.containsElementNamed("MeanTemperature")) stop("Please include variable 'MeanTemperature' in weather input.");
  NumericVector MeanTemperature = meteo["MeanTemperature"];
  NumericVector WindSpeed(numDays, NA_REAL);
  if(meteo.containsElementNamed("WindSpeed")) WindSpeed = meteo["WindSpeed"];
  NumericVector PET = NumericVector(numDays,0.0);
  if(transpirationMode=="Granier") {
    if(!meteo.containsElementNamed("PET")) stop("Please include variable 'PET' in weather input.");
    PET = meteo["PET"];
    if(control["snowpack"]) {
      if(!meteo.containsElementNamed("Radiation")) stop("If 'snowpack = TRUE', variable 'Radiation' must be provided.");
      else Radiation = meteo["Radiation"];
    }
  } else if(transpirationMode=="Sperry") {
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
  }
  CharacterVector dateStrings = meteo.attr("row.names");
  
  IntegerVector DOY = date2doy(dateStrings);
  IntegerVector Photoperiod = date2photoperiod(dateStrings, latrad);
  
  //Canpopy parameters
  List canopyParams = x["canopy"];
  

  //Plant input
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  int numCohorts = above.nrow();
  

  //Soil input
  NumericVector Water_FC = waterFC(soil, soilFunctions);
  int nlayers = Water_FC.size();
  
  //Detailed subday results
  List subdailyRes(numDays);
  

  //Water balance output variables
  NumericVector LAIcell(numDays),LAIcelldead(numDays);
  NumericVector Cm(numDays);
  NumericVector LgroundPAR(numDays);
  NumericVector LgroundSWR(numDays);
  NumericVector Runoff(numDays);
  NumericVector Rain(numDays);
  NumericVector Snow(numDays);
  NumericVector Snowmelt(numDays);
  NumericVector NetRain(numDays);
  NumericVector Interception(numDays);
  NumericVector Infiltration(numDays);
  NumericVector DeepDrainage(numDays);
  NumericVector SoilEvaporation(numDays);
  NumericVector Transpiration(numDays);
  NumericVector PlantExtraction(numDays);
  NumericVector HydraulicRedistribution(numDays, 0.0);
  
  NumericMatrix Eplantdays(numDays, nlayers);
  NumericMatrix HydrIndays(numDays, nlayers);
  NumericMatrix Wdays(numDays, nlayers); //Soil moisture content in relation to field capacity
  NumericMatrix psidays(numDays, nlayers);
  NumericMatrix MLdays(numDays, nlayers);
  NumericVector WaterTable(numDays, NA_REAL);
  NumericVector MLTot(numDays, 0.0);
  NumericVector SWE(numDays, 0.0);
  
  //EnergyBalance output variables
  NumericVector Tatm_mean(numDays, NA_REAL);
  NumericVector Tatm_min(numDays, NA_REAL);
  NumericVector Tatm_max(numDays, NA_REAL);
  NumericVector Tcan_mean(numDays, NA_REAL);
  NumericVector Tcan_min(numDays, NA_REAL);
  NumericVector Tcan_max(numDays, NA_REAL);
  NumericVector Tsoil_mean(numDays, NA_REAL);
  NumericVector Tsoil_min(numDays, NA_REAL);
  NumericVector Tsoil_max(numDays, NA_REAL);
  NumericVector SWRcanin(numDays, NA_REAL);
  NumericVector LWRcanin(numDays, NA_REAL);
  NumericVector LWRcanout(numDays, NA_REAL);
  NumericVector LEcan_heat(numDays, NA_REAL);
  NumericVector LEsoil_heat(numDays, NA_REAL);
  NumericVector Hcan_heat(numDays, NA_REAL);
  NumericVector Ebalcan(numDays, NA_REAL);
  NumericVector SWRsoilin(numDays, NA_REAL);
  NumericVector LWRsoilin(numDays, NA_REAL);
  NumericVector LWRsoilout(numDays, NA_REAL);
  NumericVector Ebalsoil(numDays, NA_REAL);
  NumericVector LWRsoilcan(numDays, NA_REAL);
  NumericVector Hcansoil(numDays, NA_REAL);
  NumericVector RAcan(numDays, NA_REAL);
  NumericVector RAsoil(numDays, NA_REAL);
  
  //Plant output variables
  NumericMatrix PlantPsi(numDays, numCohorts);
  NumericMatrix dEdP(numDays, numCohorts);
  NumericMatrix LeafPsiMin(numDays, numCohorts);
  NumericMatrix LeafPsiMax(numDays, numCohorts);
  NumericMatrix LeafPsiMin_SL(numDays, numCohorts);
  NumericMatrix LeafPsiMax_SL(numDays, numCohorts);
  NumericMatrix LeafPsiMin_SH(numDays, numCohorts);
  NumericMatrix LeafPsiMax_SH(numDays, numCohorts);
  NumericMatrix StemPsi(numDays, numCohorts);
  NumericMatrix RootPsi(numDays, numCohorts);
  NumericMatrix StemPLC(numDays, numCohorts);
  NumericMatrix PlantWaterBalance(numDays, numCohorts);
  List RhizoPsi(numCohorts);
  for(int c=0;c<numCohorts;c++) {
    NumericMatrix nm = NumericMatrix(numDays, nlayers);
    nm.attr("dimnames") = List::create(meteo.attr("row.names"), seq(1,nlayers)) ;
    RhizoPsi[c] = nm;
  }
  RhizoPsi.attr("names") = above.attr("row.names");
  
  NumericMatrix PlantStress(numDays, numCohorts);
  NumericMatrix StemRWC(numDays, numCohorts);
  NumericMatrix LeafRWC(numDays, numCohorts);
  NumericMatrix PlantTranspiration(numDays, numCohorts);
  NumericMatrix PlantPhotosynthesis(numDays, numCohorts);
  NumericMatrix PlantNetPhotosynthesis(numDays, numCohorts);
  NumericMatrix PlantGrossPhotosynthesis(numDays, numCohorts);
  NumericVector EplantCohTot(numCohorts, 0.0);
  NumericMatrix PlantAbsSWRFraction(numDays, numCohorts);
  NumericMatrix PlantAbsSWR(numDays, numCohorts);
  NumericMatrix PlantAbsLWR(numDays, numCohorts);
  NumericMatrix PlantLAI(numDays, numCohorts);
  
  
  NumericVector Wini = soil["W"];
  Wdays(0,_) = Wini;
  NumericVector initialContent = water(soil, soilFunctions);
  double initialSnowContent = soil["SWE"];
  if(verbose) {
    Rcout<<"Initial soil water content (mm): "<< sum(initialContent)<<"\n";
    Rcout<<"Initial snowpack content (mm): "<< initialSnowContent<<"\n";
  }

  if(verbose) Rcout << "Performing daily simulations ";
  NumericVector Eplanttot(numDays,0.0);
  List s;
  for(int i=0;i<numDays;i++) {
      if(verbose & (i%10 == 0)) Rcout<<".";//<<i;
      
      double wind = WindSpeed[i];
      if(NumericVector::is_na(wind)) wind = control["defaultWindSpeed"]; //Default 1 m/s -> 10% of fall every day
      if(wind<0.1) wind = 0.1; //Minimum windspeed abovecanopy
      
      //If DOY == 1 reset PLC (Growth assumed)
      if(cavitationRefill=="annual") {
        if(DOY[i]==1) {
          DataFrame internalWater = Rcpp::as<Rcpp::DataFrame>(x["internalWater"]);
          if(transpirationMode=="Granier") {
            NumericVector PLC = Rcpp::as<Rcpp::NumericVector>(internalWater["PLC"]);
            for(int j=0;j<PLC.length();j++) PLC[j] = 0.0;
          } else {
            NumericVector StemPLC = Rcpp::as<Rcpp::NumericVector>(internalWater["PLCstem"]);
            for(int j=0;j<StemPLC.length();j++) StemPLC[j] = 0.0;
          }
        }
      }

      if(unlimitedSoilWater) {
        NumericVector W = soil["W"];
        for(int h=0;h<W.size();h++) W[h] = 1.0;
      }
      
      //1. Phenology and leaf fall
      if(leafPhenology) updateLeaves(x, DOY[i], Photoperiod[i], MeanTemperature[i], wind);
      
      //2. Water balance and photosynthesis
      if(transpirationMode=="Granier") {
        double er = erFactor(DOY[i], PET[i], Precipitation[i]);
        s = spwbDay1(x, soil, MeanTemperature[i], PET[i], Precipitation[i], er, 0.0, 
                     Radiation[i], elevation, verbose); //No Runon in simulations for a single cell
      } else if(transpirationMode=="Sperry") {
        int ntimesteps = control["ndailysteps"];
        double tstep = 86400.0/((double) ntimesteps);
        std::string c = as<std::string>(dateStrings[i]);
        int J = meteoland::radiation_julianDay(std::atoi(c.substr(0, 4).c_str()),std::atoi(c.substr(5,2).c_str()),std::atoi(c.substr(8,2).c_str()));
        double delta = meteoland::radiation_solarDeclination(J);
        double solarConstant = meteoland::radiation_solarConstant(J);
        if(NumericVector::is_na(aspect)) aspect = 0.0;
        if(NumericVector::is_na(slope)) slope = 0.0;
        double asprad = aspect * (PI/180.0);
        double slorad = slope * (PI/180.0);
        double tmin = MinTemperature[i];
        double tmax = MaxTemperature[i];
        double tmaxPrev = tmax;
        double tminPrev = tmin;
        double tminNext = tmin;
        if(i>0) {
          tmaxPrev = MaxTemperature[i-1];
          tminPrev = MinTemperature[i-1];
        }
        if(i<(numDays-1)) tminNext = MinTemperature[i+1]; 
        double rhmin = MinRelativeHumidity[i];
        double rhmax = MaxRelativeHumidity[i];
        double rad = Radiation[i];
        PET[i] = meteoland::penman(latrad, elevation, slorad, asprad, J, tmin, tmax, rhmin, rhmax, rad, wind);
        double er = erFactor(DOY[i], PET[i], Precipitation[i]);
        s = spwbDay2(x, soil, tmin, tmax, tminPrev, tmaxPrev, tminNext,
                     rhmin, rhmax, rad, wind, 
                     latitude, elevation, slope, aspect,
                     solarConstant, delta, Precipitation[i], PET[i], 
                     er, 0.0, verbose);
        List Plants = Rcpp::as<Rcpp::List>(s["Plants"]);
        List PlantsInst = Rcpp::as<Rcpp::List>(s["PlantsInst"]);
        List ShadeLeaves = Rcpp::as<Rcpp::List>(PlantsInst["ShadeLeaves"]);
        List SunlitLeaves = Rcpp::as<Rcpp::List>(PlantsInst["SunlitLeaves"]);
        
        NumericMatrix SWR_SL = Rcpp::as<Rcpp::NumericMatrix>(SunlitLeaves["Abs_SWR"]);
        NumericMatrix SWR_SH = Rcpp::as<Rcpp::NumericMatrix>(ShadeLeaves["Abs_SWR"]);
        NumericMatrix LWR_SL = Rcpp::as<Rcpp::NumericMatrix>(SunlitLeaves["Abs_LWR"]);
        NumericMatrix LWR_SH = Rcpp::as<Rcpp::NumericMatrix>(ShadeLeaves["Abs_LWR"]);
        for(int j=0;j<numCohorts;j++) {
          for(int n=0;n<ntimesteps;n++){
            PlantAbsSWR(i,j) += 0.000001*(SWR_SL(j,n)+SWR_SH(j,n))*tstep;
            PlantAbsLWR(i,j) += 0.000001*(LWR_SL(j,n)+LWR_SH(j,n))*tstep;
          }
        }
        List EB = Rcpp::as<Rcpp::List>(s["EnergyBalance"]);
        DataFrame Tinst = Rcpp::as<Rcpp::DataFrame>(EB["Temperature"]); 
        DataFrame CEBinst = Rcpp::as<Rcpp::DataFrame>(EB["CanopyEnergyBalance"]); 
        DataFrame SEBinst = Rcpp::as<Rcpp::DataFrame>(EB["SoilEnergyBalance"]); 
        NumericVector Tcan = Rcpp::as<Rcpp::NumericVector>(Tinst["Tcan"]);
        NumericVector Tsoil = Rcpp::as<Rcpp::NumericVector>(Tinst["Tsoil.1"]);

        SWRcanin[i] = 0.000001*sum(Rcpp::as<Rcpp::NumericVector>(CEBinst["SWRcanin"]))*tstep;
        LWRcanin[i] = 0.000001*sum(Rcpp::as<Rcpp::NumericVector>(CEBinst["LWRcanin"]))*tstep;
        LWRcanout[i] = 0.000001*sum(Rcpp::as<Rcpp::NumericVector>(CEBinst["LWRcanout"]))*tstep;
        LEcan_heat[i] = 0.000001*sum(Rcpp::as<Rcpp::NumericVector>(CEBinst["LEcan"]))*tstep;
        Hcan_heat[i] = 0.000001*sum(Rcpp::as<Rcpp::NumericVector>(CEBinst["Hcan"]))*tstep;
        Ebalcan[i] = 0.000001*sum(Rcpp::as<Rcpp::NumericVector>(CEBinst["Ebalcan"]))*tstep;
        LWRsoilcan[i] = 0.000001*sum(Rcpp::as<Rcpp::NumericVector>(CEBinst["LWRsoilcan"]))*tstep;
        RAcan[i] = sum(Rcpp::as<Rcpp::NumericVector>(CEBinst["RAcan"]))/((double) ntimesteps);
        SWRsoilin[i] = 0.000001*sum(Rcpp::as<Rcpp::NumericVector>(SEBinst["SWRsoilin"]))*tstep;
        LWRsoilin[i] = 0.000001*sum(Rcpp::as<Rcpp::NumericVector>(SEBinst["LWRsoilin"]))*tstep;
        LWRsoilout[i] = 0.000001*sum(Rcpp::as<Rcpp::NumericVector>(SEBinst["LWRsoilout"]))*tstep;
        LEsoil_heat[i] = 0.000001*sum(Rcpp::as<Rcpp::NumericVector>(SEBinst["LEsoil"]))*tstep;
        Hcansoil[i] = 0.000001*sum(Rcpp::as<Rcpp::NumericVector>(SEBinst["Hcansoil"]))*tstep;
        Ebalsoil[i] = 0.000001*sum(Rcpp::as<Rcpp::NumericVector>(SEBinst["Ebalsoil"]))*tstep;
        RAsoil[i] = sum(Rcpp::as<Rcpp::NumericVector>(SEBinst["RAsoil"]))/((double) ntimesteps);
        
        Tatm_min[i] = MinTemperature[i];
        Tatm_max[i] = MaxTemperature[i];
        Tatm_mean[i] = MeanTemperature[i];
        Tcan_min[i] = min(Tcan);
        Tcan_max[i] = max(Tcan);
        Tcan_mean[i] = sum(Tcan)/((double) ntimesteps);
        Tsoil_min[i] = min(Tsoil);
        Tsoil_max[i] = max(Tsoil);
        Tsoil_mean[i] = sum(Tsoil)/((double) ntimesteps);
      }
      List stand = s["Stand"];
      LgroundPAR[i] = stand["LgroundPAR"];
      LgroundSWR[i] = stand["LgroundSWR"];
      LAIcell[i] = stand["LAIcell"];
      LAIcelldead[i] = stand["LAIcelldead"];
      Cm[i] = stand["Cm"];
      
      List db = s["WaterBalance"];
      DeepDrainage[i] = db["DeepDrainage"];
      Infiltration[i] = db["Infiltration"];
      Runoff[i] = db["Runoff"];
      Rain[i] = db["Rain"];
      Snow[i] = db["Snow"];
      Snowmelt[i] = db["Snowmelt"];
      NetRain[i] = db["NetRain"];
      PlantExtraction[i] = db["PlantExtraction"];
      if(transpirationMode=="Sperry")  {
        HydraulicRedistribution[i] = db["HydraulicRedistribution"];
      }
      Transpiration[i] = db["Transpiration"];
      SoilEvaporation[i] = db["SoilEvaporation"];
      Interception[i] = Rain[i]-NetRain[i];
      
      List sb = s["Soil"];
      NumericVector psi = sb["psi"];
      psidays(i,_) = psi;
      NumericVector EplantVec = sb["PlantExtraction"];
      SWE[i] = soil["SWE"];
        
      List Plants = s["Plants"];
      PlantLAI(i,_) = Rcpp::as<Rcpp::NumericVector>(Plants["LAI"]);
      NumericVector EplantCoh = Plants["Transpiration"];
      Eplantdays(i,_) = EplantVec;
      PlantTranspiration(i,_) = EplantCoh;
      PlantStress(i,_) = Rcpp::as<Rcpp::NumericVector>(Plants["DDS"]);
      if(transpirationMode=="Sperry") {
        PlantGrossPhotosynthesis(i,_) = Rcpp::as<Rcpp::NumericVector>(Plants["GrossPhotosynthesis"]);
        PlantNetPhotosynthesis(i,_) = Rcpp::as<Rcpp::NumericVector>(Plants["NetPhotosynthesis"]);
        NumericVector HydrInVec = sb["HydraulicInput"];
        HydrIndays(i,_) = HydrInVec;
        LeafPsiMin(i,_) = Rcpp::as<Rcpp::NumericVector>(Plants["LeafPsiMin"]);
        LeafPsiMax(i,_) = Rcpp::as<Rcpp::NumericVector>(Plants["LeafPsiMax"]);
        LeafPsiMin_SL(i,_) = Rcpp::as<Rcpp::NumericVector>(Plants["LeafPsiMin_SL"]);
        LeafPsiMax_SL(i,_) = Rcpp::as<Rcpp::NumericVector>(Plants["LeafPsiMax_SL"]);
        LeafPsiMin_SH(i,_) = Rcpp::as<Rcpp::NumericVector>(Plants["LeafPsiMin_SH"]);
        LeafPsiMax_SH(i,_) = Rcpp::as<Rcpp::NumericVector>(Plants["LeafPsiMax_SH"]);
        RootPsi(i,_) = Rcpp::as<Rcpp::NumericVector>(Plants["RootPsi"]); 
        StemPsi(i,_) = Rcpp::as<Rcpp::NumericVector>(Plants["StemPsi"]); 
        StemPLC(i,_) = Rcpp::as<Rcpp::NumericVector>(Plants["StemPLC"]); 
        PlantWaterBalance(i,_) = Rcpp::as<Rcpp::NumericVector>(Plants["WaterBalance"]); 
        dEdP(i,_) = Rcpp::as<Rcpp::NumericVector>(Plants["dEdP"]); 
        NumericMatrix RhizoPsiStep = Rcpp::as<Rcpp::NumericMatrix>(s["RhizoPsi"]);
        for(int c=0;c<numCohorts;c++) {
          NumericMatrix nm = Rcpp::as<Rcpp::NumericMatrix>(RhizoPsi[c]);
          nm(i,_) =  RhizoPsiStep(c,_);
        }
      } else {
        PlantPhotosynthesis(i,_) = Rcpp::as<Rcpp::NumericVector>(Plants["Photosynthesis"]);
        PlantAbsSWRFraction(i,_) = Rcpp::as<Rcpp::NumericVector>(Plants["AbsorbedSWRFraction"]);
        PlantPsi(i,_) = Rcpp::as<Rcpp::NumericVector>(Plants["psi"]);
      }
      EplantCohTot = EplantCohTot + EplantCoh;
      Eplanttot[i] = sum(EplantCoh);
      if(transpirationMode=="Sperry"){
        StemRWC(i,_) = as<Rcpp::NumericVector>(Plants["StemRWC"]);
        LeafRWC(i,_) = as<Rcpp::NumericVector>(Plants["LeafRWC"]); 
      }
      
      if(subdailyResults) {
        subdailyRes[i] = clone(s);
      }
      if(i<(numDays-1)) Wdays(i+1,_) = as<Rcpp::NumericVector>(soil["W"]);
      WaterTable[i] = waterTableDepth(soil, soilFunctions);
  }
  if(verbose) Rcout << "done.\n";
  
  for(int l=0;l<nlayers;l++) {
    MLdays(_,l) = Wdays(_,l)*Water_FC[l]; 
    MLTot = MLTot + MLdays(_,l);
  }
  NumericVector Evapotranspiration = Transpiration+SoilEvaporation + Interception;
  
  if(verbose) {
    NumericVector finalContent = water(soil, soilFunctions);
    double finalSnowContent = soil["SWE"];
    Rcout<<"Final soil water content (mm): "<< sum(finalContent)<<"\n";
    Rcout<<"Final snowpack content (mm): "<< finalSnowContent<<"\n";
    
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
    
    double soil_wb = (Rainfallsum - Interceptionsum) + Snowmeltsum - Runoffsum - DeepDrainagesum - SoilEvaporationsum - sum(PlantExtraction);
    double snowpack_wb = Snowsum - Snowmeltsum;
    Rcout<<"Change in soil water content (mm): "<< sum(finalContent) - sum(initialContent)<<"\n";
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
    Rcout<<" Transpiration (mm) "  <<round(Transpirationsum) <<"\n";
    if(transpirationMode =="Sperry") {
      Rcout<<"  Plant extraction from soil (mm) " << round(sum(PlantExtraction));
      Rcout<<"  Plant water balance (mm) " << round(sum(PlantWaterBalance));
      Rcout<<" Hydraulic redistribution (mm) " << round(sum(HydraulicRedistribution)) <<"\n";
    }
  }

   DataFrame SWB;
   if(transpirationMode=="Granier") {
     SWB = DataFrame::create(_["W"]=Wdays, _["ML"]=MLdays,_["MLTot"]=MLTot,
                             _["WTD"] = WaterTable,
                             _["SWE"] = SWE, 
                             _["PlantExt"]=Eplantdays,
                             _["psi"]=psidays); 
   } else {
     SWB = DataFrame::create(_["W"]=Wdays, _["ML"]=MLdays,_["MLTot"]=MLTot,
                             _["WTD"] = WaterTable,
                             _["SWE"] = SWE, 
                             _["PlantExt"]=Eplantdays,
                             _["HydraulicInput"] = HydrIndays,
                             _["psi"]=psidays); 
   }
   SWB.attr("row.names") = meteo.attr("row.names") ;
   DataFrame Stand = DataFrame::create(_["LAIcell"]=LAIcell, _["LAIcelldead"] = LAIcelldead,  _["Cm"]=Cm, 
                               _["LgroundPAR"] = LgroundPAR, _["LgroundSWR"] = LgroundSWR);
   Stand.attr("row.names") = meteo.attr("row.names") ;
   DataFrame DWB;
   if(transpirationMode=="Granier") {
     DWB = DataFrame::create(_["PET"]=PET, 
                             _["Precipitation"] = Precipitation, _["Rain"] = Rain, _["Snow"] = Snow, 
                             _["NetRain"]=NetRain, _["Snowmelt"] = Snowmelt, _["Infiltration"]=Infiltration, _["Runoff"]=Runoff, _["DeepDrainage"]=DeepDrainage, 
                             _["Evapotranspiration"]=Evapotranspiration,_["Interception"] = Interception, _["SoilEvaporation"]=SoilEvaporation,
                             _["PlantExtraction"] = PlantExtraction, _["Transpiration"]=Transpiration);
   } else {
     DWB = DataFrame::create(_["PET"]=PET, 
                             _["Precipitation"] = Precipitation, _["Rain"] = Rain, _["Snow"] = Snow, 
                             _["NetRain"]=NetRain, _["Snowmelt"] = Snowmelt, _["Infiltration"]=Infiltration, _["Runoff"]=Runoff, _["DeepDrainage"]=DeepDrainage, 
                             _["Evapotranspiration"]=Evapotranspiration,_["Interception"] = Interception, _["SoilEvaporation"]=SoilEvaporation,
                             _["PlantExtraction"] = PlantExtraction, _["Transpiration"]=Transpiration, 
                             _["HydraulicRedistribution"] = HydraulicRedistribution);
   }
   DWB.attr("row.names") = meteo.attr("row.names") ;

   if(transpirationMode=="Granier") PlantAbsSWRFraction.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")); 
  PlantTranspiration.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names"));
  PlantStress.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  StemPLC.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  StemRWC.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  LeafRWC.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  PlantPsi.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  dEdP.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  LeafPsiMin.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  LeafPsiMax.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  LeafPsiMin_SL.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  LeafPsiMax_SL.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  LeafPsiMin_SH.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  LeafPsiMax_SH.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  StemPsi.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  PlantWaterBalance.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  RootPsi.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  PlantPhotosynthesis.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  PlantGrossPhotosynthesis.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  PlantNetPhotosynthesis.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  PlantAbsSWR.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  PlantAbsLWR.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  PlantLAI.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  
  DataFrame DEB = DataFrame::create(_["SWRcanin"] = SWRcanin, _["LWRcanin"] = LWRcanin, _["LWRcanout"] = LWRcanout,
                                    _["LEcan"] = LEcan_heat, _["Hcan"] = Hcan_heat, _["LWRsoilcan"] = LWRsoilcan, _["Ebalcan"] = Ebalcan, 
                                    _["Hcansoil"] = Hcansoil, _["SWRsoilin"] = SWRsoilin, _["LWRsoilin"] = LWRsoilin, _["LWRsoilout"] = LWRsoilout,
                                    _["LEsoil"] = LEsoil_heat, _["Ebalsoil"] = Ebalsoil, _["RAcan"] = RAcan, _["RAsoil"] = RAsoil);  
  DataFrame DT = DataFrame::create(_["Tatm_mean"] = Tatm_mean, _["Tatm_min"] = Tatm_min, _["Tatm_max"] = Tatm_max,
                                   _["Tcan_mean"] = Tcan_mean, _["Tcan_min"] = Tcan_min, _["Tcan_max"] = Tcan_max,
                                   _["Tsoil_mean"] = Tsoil_mean, _["Tsoil_min"] = Tsoil_min, _["Tsoil_max"] = Tsoil_max);
  DEB.attr("row.names") = meteo.attr("row.names") ;
  DT.attr("row.names") = meteo.attr("row.names") ;
  subdailyRes.attr("names") = meteo.attr("row.names") ;
  
  NumericVector topo = NumericVector::create(elevation, slope, aspect);
  topo.attr("names") = CharacterVector::create("elevation", "slope", "aspect");
  List l;
  if(transpirationMode=="Granier") {
    l = List::create(Named("latitude") = latitude,
                     Named("topography") = topo,
                     Named("spwbInput") = spwbInput,
                     Named("soilInput") = soilInput,
                     Named("WaterBalance")=DWB, 
                     Named("Soil")=SWB,
                     Named("Stand")=Stand, 
                     Named("PlantLAI") = PlantLAI,
                     Named("PlantAbsorbedSWRFraction") = PlantAbsSWRFraction,
                     Named("PlantTranspiration") = PlantTranspiration,
                     Named("PlantPhotosynthesis") = PlantPhotosynthesis,
                     Named("PlantPsi") = PlantPsi, 
                     Named("PlantStress") = PlantStress,
                     Named("subdaily") =  subdailyRes);
  } else {
    CharacterVector ln = CharacterVector(31);
    l = List(31);
    l[0] = latitude;
    ln[0] = "latitude";
    l[1] = topo;
    ln[1] = "topography";
    l[2] = spwbInput;
    ln[2] = "spwbInput";
    l[3] = soilInput;
    ln[3] = "soilInput";
    l[4] = DWB;
    ln[4] = "WaterBalance";
    l[5] = SWB;
    ln[5] = "Soil";
    l[6] = Stand;
    ln[6] = "Stand";
    l[7] = DEB;
    ln[7] = "EnergyBalance";
    l[8] = DT;
    ln[8] = "Temperature";
    l[9] = PlantLAI;
    ln[9] = "PlantLAI";
    l[10] = PlantAbsSWR;
    ln[10] = "PlantAbsorbedSWR";
    l[11] = PlantAbsLWR;
    ln[11] = "PlantAbsorbedLWR";
    l[12] = PlantTranspiration;
    ln[12] = "PlantTranspiration";
    l[13] = PlantGrossPhotosynthesis;
    ln[13] = "PlantGrossPhotosynthesis";
    l[14] = PlantNetPhotosynthesis;
    ln[14] = "PlantNetPhotosynthesis";
    l[15] = dEdP;
    ln[15] = "dEdP";
    l[16] = LeafPsiMin;
    ln[16] = "LeafPsiMin";
    l[17] = LeafPsiMax;
    ln[17] = "LeafPsiMax";
    l[18] = LeafPsiMin_SL;
    ln[18] = "LeafPsiMin_SL";
    l[19] = LeafPsiMax_SL;
    ln[19] = "LeafPsiMax_SL";
    l[20] = LeafPsiMin_SH;
    ln[20] = "LeafPsiMin_SH";
    l[21] = LeafPsiMax_SH;
    ln[21] = "LeafPsiMax_SH";
    l[22] = LeafRWC;
    ln[22] = "LeafRWC";
    l[23] = StemPsi;
    ln[23] = "StemPsi";
    l[24] = StemPLC;
    ln[24] = "StemPLC";
    l[25] = StemRWC;
    ln[25] = "StemRWC";
    l[26] = RootPsi;
    ln[26] = "RootPsi";
    l[27] = RhizoPsi;
    ln[27] = "RhizoPsi";
    l[28] = PlantStress;
    ln[28] = "PlantStress";
    l[29] = PlantWaterBalance;
    ln[29] = "PlantWaterBalance";
    l[30] = subdailyRes;
    ln[30] = "subdaily";
    l.attr("names") = ln;
  }
  l.attr("class") = CharacterVector::create("spwb","list");
  return(l);
}


// [[Rcpp::export("pwb")]]
List pwb(List x, List soil, DataFrame meteo, NumericMatrix W,
            double latitude, double elevation = NA_REAL, double slope = NA_REAL, double aspect = NA_REAL, 
            NumericVector canopyEvaporation = NumericVector(0), 
            NumericVector snowMelt = NumericVector(0), 
            NumericVector soilEvaporation = NumericVector(0)) {
  List control = x["control"];
  String transpirationMode = control["transpirationMode"];
  String soilFunctions = control["soilFunctions"];
  String cavitationRefill = control["cavitationRefill"];

  bool verbose = control["verbose"];
  bool subdailyResults = control["subdailyResults"];
  bool leafPhenology = control["leafPhenology"];
  
  
  
  //Store input
  List spwbInput = clone(x);
  List soilInput = clone(soil);
  
  if(NumericVector::is_na(latitude)) stop("Value for 'latitude' should not be missing.");
  double latrad = latitude * (PI/180.0);

    //Meteorological input    
  NumericVector MinTemperature, MaxTemperature;
  NumericVector MinRelativeHumidity, MaxRelativeHumidity;
  NumericVector Radiation;
  if(!meteo.containsElementNamed("Precipitation")) stop("Please include variable 'Precipitation' in weather input.");
  NumericVector Precipitation = meteo["Precipitation"];
  int numDays = Precipitation.size();
  if(!meteo.containsElementNamed("MeanTemperature")) stop("Please include variable 'MeanTemperature' in weather input.");
  NumericVector MeanTemperature = meteo["MeanTemperature"];
  NumericVector WindSpeed(numDays, NA_REAL);
  if(meteo.containsElementNamed("WindSpeed")) WindSpeed = meteo["WindSpeed"];
  NumericVector PET = NumericVector(numDays,0.0);
  if(transpirationMode=="Granier") {
    if(!meteo.containsElementNamed("PET")) stop("Please include variable 'PET' in weather input.");
    PET = meteo["PET"];
  } else if(transpirationMode=="Sperry") {
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
  CharacterVector dateStrings = meteo.attr("row.names");
  
  IntegerVector DOY = date2doy(dateStrings);
  IntegerVector Photoperiod = date2photoperiod(dateStrings, latrad);
  
  //Canpopy parameters
  List canopyParams = x["canopy"];
  
  
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
  NumericVector Tatm_mean(numDays, NA_REAL);
  NumericVector Tatm_min(numDays, NA_REAL);
  NumericVector Tatm_max(numDays, NA_REAL);
  NumericVector Tcan_mean(numDays, NA_REAL);
  NumericVector Tcan_min(numDays, NA_REAL);
  NumericVector Tcan_max(numDays, NA_REAL);
  NumericVector Tsoil_mean(numDays, NA_REAL);
  NumericVector Tsoil_min(numDays, NA_REAL);
  NumericVector Tsoil_max(numDays, NA_REAL);
  NumericVector SWRcanin(numDays, NA_REAL);
  NumericVector LWRcanin(numDays, NA_REAL);
  NumericVector LWRcanout(numDays, NA_REAL);
  NumericVector LEcan_heat(numDays, NA_REAL);
  NumericVector LEsoil_heat(numDays, NA_REAL);
  NumericVector Hcan_heat(numDays, NA_REAL);
  NumericVector Ebalcan(numDays, NA_REAL);
  NumericVector SWRsoilin(numDays, NA_REAL);
  NumericVector LWRsoilin(numDays, NA_REAL);
  NumericVector LWRsoilout(numDays, NA_REAL);
  NumericVector Ebalsoil(numDays, NA_REAL);
  NumericVector LWRsoilcan(numDays, NA_REAL);
  NumericVector Hcansoil(numDays, NA_REAL);
  NumericVector RAcan(numDays, NA_REAL);
  NumericVector RAsoil(numDays, NA_REAL);
  

  //Soil output variables
  NumericMatrix Wdays(numDays, nlayers); //Soil moisture content in relation to field capacity
  NumericMatrix psidays(numDays, nlayers);
  NumericMatrix Eplantdays(numDays, nlayers);
  
  //Plant output variables
  NumericMatrix PlantPsi(numDays, numCohorts);
  NumericMatrix dEdP(numDays, numCohorts);
  NumericMatrix LeafPsiMin(numDays, numCohorts);
  NumericMatrix LeafPsiMax(numDays, numCohorts);
  NumericMatrix LeafPsiMin_SL(numDays, numCohorts);
  NumericMatrix LeafPsiMax_SL(numDays, numCohorts);
  NumericMatrix LeafPsiMin_SH(numDays, numCohorts);
  NumericMatrix LeafPsiMax_SH(numDays, numCohorts);
  NumericMatrix StemPsi(numDays, numCohorts);
  NumericMatrix RootPsi(numDays, numCohorts);
  NumericMatrix StemPLC(numDays, numCohorts);
  NumericMatrix PlantWaterBalance(numDays, numCohorts);
  List RhizoPsi(numCohorts);
  for(int c=0;c<numCohorts;c++) {
    NumericMatrix nm = NumericMatrix(numDays, nlayers);
    nm.attr("dimnames") = List::create(meteo.attr("row.names"), seq(1,nlayers)) ;
    RhizoPsi[c] = nm;
  }
  RhizoPsi.attr("names") = above.attr("row.names");
  
  NumericMatrix PlantStress(numDays, numCohorts);
  NumericMatrix StemRWC(numDays, numCohorts);
  NumericMatrix LeafRWC(numDays, numCohorts);
  NumericMatrix PlantTranspiration(numDays, numCohorts);
  NumericMatrix PlantPhotosynthesis(numDays, numCohorts);
  NumericMatrix PlantGrossPhotosynthesis(numDays, numCohorts);
  NumericMatrix PlantNetPhotosynthesis(numDays, numCohorts);
  NumericVector EplantCohTot(numCohorts, 0.0);
  NumericMatrix PlantAbsSWRFraction(numDays, numCohorts);
  NumericMatrix PlantAbsSWR(numDays, numCohorts);
  NumericMatrix PlantAbsLWR(numDays, numCohorts);
  NumericMatrix PlantLAI(numDays, numCohorts);
  
  

  if(verbose) Rcout << "Performing daily simulations ";
  NumericVector Eplanttot(numDays,0.0);
  List s;
  for(int i=0;i<numDays;i++) {
    if(verbose & (i%10 == 0)) Rcout<<".";//<<i;
    double wind = WindSpeed[i];
    if(NumericVector::is_na(wind)) wind = control["defaultWindSpeed"]; //Default 1 m/s -> 10% of fall every day
    if(wind<0.1) wind = 0.1; //Minimum windspeed abovecanopy
    
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
            NumericVector StemPLC = Rcpp::as<Rcpp::NumericVector>(internalWater["PLCstem"]);
            for(int j=0;j<StemPLC.length();j++) StemPLC[j] = 0.0;
          }
        }
      }
      
      
      //1. Phenology and leaf fall
      if(leafPhenology) updateLeaves(x, DOY[i], Photoperiod[i], MeanTemperature[i], wind);
      
    
    int ntimesteps = control["ndailysteps"];
    double tstep = 86400.0/((double) ntimesteps);
    
    //2. transpiration and photosynthesis
    if(transpirationMode=="Granier") {
      s = transpirationGranier(x, soil, MeanTemperature[i], PET[i], true, true);
    } else if(transpirationMode=="Sperry") {
      std::string c = as<std::string>(dateStrings[i]);
      int J = meteoland::radiation_julianDay(std::atoi(c.substr(0, 4).c_str()),std::atoi(c.substr(5,2).c_str()),std::atoi(c.substr(8,2).c_str()));
      double delta = meteoland::radiation_solarDeclination(J);
      double solarConstant = meteoland::radiation_solarConstant(J);
      double tmin = MinTemperature[i];
      double tmax = MaxTemperature[i];
      double tmaxPrev = tmax;
      double tminPrev = tmin;
      double tminNext = tmin;
      if(i>0) {
        tmaxPrev = MaxTemperature[i-1];
        tminPrev = MinTemperature[i-1];
      }
      if(i<(numDays-1)) tminNext = MinTemperature[i+1]; 
      double rhmin = MinRelativeHumidity[i];
      double rhmax = MaxRelativeHumidity[i];
      double rad = Radiation[i];
      double prec = Precipitation[i];
      
      s = transpirationSperry(x, soil, tmin, tmax, tminPrev, tmaxPrev, tminNext, 
                              rhmin, rhmax, rad, wind, 
                              latitude, elevation, slope, aspect,
                              solarConstant, delta, prec,
                              canopyEvaporation[i], snowMelt[i], soilEvaporation[i],
                              verbose, NA_INTEGER, true, true);
      
    }
    List Plants = s["Plants"];
    PlantLAI(i,_) = Rcpp::as<Rcpp::NumericVector>(Plants["LAI"]);
    NumericVector EplantCoh = Plants["Transpiration"];
    NumericMatrix SoilWaterExtract = s["Extraction"];
    for(int l=0;l<nlayers;l++) {
      Eplantdays(i,l) = sum(SoilWaterExtract(_,l));
    }

    PlantExtraction[i] = sum(SoilWaterExtract);
    Transpiration[i] = sum(EplantCoh);
    NumericVector HydrInVec(nlayers, 0.0);
    PlantTranspiration(i,_) = EplantCoh;
    PlantStress(i,_) = Rcpp::as<Rcpp::NumericVector>(Plants["DDS"]);
    
    if(transpirationMode=="Granier") {
      PlantPhotosynthesis(i,_) = Rcpp::as<Rcpp::NumericVector>(Plants["Photosynthesis"]);
      PlantAbsSWRFraction(i,_) = Rcpp::as<Rcpp::NumericVector>(Plants["AbsorbedSWRFraction"]);
      PlantPsi(i,_) = Rcpp::as<Rcpp::NumericVector>(Plants["psi"]);
    }
    else if(transpirationMode=="Sperry")  {
      PlantGrossPhotosynthesis(i,_) = Rcpp::as<Rcpp::NumericVector>(Plants["GrossPhotosynthesis"]);
      PlantNetPhotosynthesis(i,_) = Rcpp::as<Rcpp::NumericVector>(Plants["NetPhotosynthesis"]);
      NumericMatrix soilLayerExtractInst = s["ExtractionInst"];
      for(int l=0;l<nlayers;l++) {
        for(int n=0;n<ntimesteps;n++) {
          HydrInVec[l] += (-1.0)*std::min(soilLayerExtractInst(l,n),0.0);
        }
      }
      HydraulicRedistribution[i] = sum(HydrInVec);
      
      List PlantsInst = Rcpp::as<Rcpp::List>(s["PlantsInst"]);
      List ShadeLeaves = Rcpp::as<Rcpp::List>(PlantsInst["ShadeLeaves"]);
      List SunlitLeaves = Rcpp::as<Rcpp::List>(PlantsInst["SunlitLeaves"]);
      
      NumericMatrix SWR_SL = Rcpp::as<Rcpp::NumericMatrix>(SunlitLeaves["Abs_SWR"]);
      NumericMatrix SWR_SH = Rcpp::as<Rcpp::NumericMatrix>(ShadeLeaves["Abs_SWR"]);
      NumericMatrix LWR_SL = Rcpp::as<Rcpp::NumericMatrix>(SunlitLeaves["Abs_LWR"]);
      NumericMatrix LWR_SH = Rcpp::as<Rcpp::NumericMatrix>(ShadeLeaves["Abs_LWR"]);
      for(int j=0;j<numCohorts;j++) {
        for(int n=0;n<ntimesteps;n++){
          PlantAbsSWR(i,j) += 0.000001*(SWR_SL(j,n)+SWR_SH(j,n))*tstep;
          PlantAbsLWR(i,j) += 0.000001*(LWR_SL(j,n)+LWR_SH(j,n))*tstep;
        }
      }
      List EB = Rcpp::as<Rcpp::List>(s["EnergyBalance"]);
      DataFrame Tinst = Rcpp::as<Rcpp::DataFrame>(EB["Temperature"]); 
      DataFrame CEBinst = Rcpp::as<Rcpp::DataFrame>(EB["CanopyEnergyBalance"]); 
      DataFrame SEBinst = Rcpp::as<Rcpp::DataFrame>(EB["SoilEnergyBalance"]); 
      NumericVector Tcan = Rcpp::as<Rcpp::NumericVector>(Tinst["Tcan"]);
      NumericVector Tsoil = Rcpp::as<Rcpp::NumericVector>(Tinst["Tsoil.1"]);
      
      SWRcanin[i] = 0.000001*sum(Rcpp::as<Rcpp::NumericVector>(CEBinst["SWRcanin"]))*tstep;
      LWRcanin[i] = 0.000001*sum(Rcpp::as<Rcpp::NumericVector>(CEBinst["LWRcanin"]))*tstep;
      LWRcanout[i] = 0.000001*sum(Rcpp::as<Rcpp::NumericVector>(CEBinst["LWRcanout"]))*tstep;
      LEcan_heat[i] = 0.000001*sum(Rcpp::as<Rcpp::NumericVector>(CEBinst["LEcan"]))*tstep;
      Hcan_heat[i] = 0.000001*sum(Rcpp::as<Rcpp::NumericVector>(CEBinst["Hcan"]))*tstep;
      Ebalcan[i] = 0.000001*sum(Rcpp::as<Rcpp::NumericVector>(CEBinst["Ebalcan"]))*tstep;
      LWRsoilcan[i] = 0.000001*sum(Rcpp::as<Rcpp::NumericVector>(CEBinst["LWRsoilcan"]))*tstep;
      RAcan[i] = sum(Rcpp::as<Rcpp::NumericVector>(CEBinst["RAcan"]))/((double) ntimesteps);
      SWRsoilin[i] = 0.000001*sum(Rcpp::as<Rcpp::NumericVector>(SEBinst["SWRsoilin"]))*tstep;
      LWRsoilin[i] = 0.000001*sum(Rcpp::as<Rcpp::NumericVector>(SEBinst["LWRsoilin"]))*tstep;
      LWRsoilout[i] = 0.000001*sum(Rcpp::as<Rcpp::NumericVector>(SEBinst["LWRsoilout"]))*tstep;
      LEsoil_heat[i] = 0.000001*sum(Rcpp::as<Rcpp::NumericVector>(SEBinst["LEsoil"]))*tstep;
      Hcansoil[i] = 0.000001*sum(Rcpp::as<Rcpp::NumericVector>(SEBinst["Hcansoil"]))*tstep;
      Ebalsoil[i] = 0.000001*sum(Rcpp::as<Rcpp::NumericVector>(SEBinst["Ebalsoil"]))*tstep;
      RAsoil[i] = sum(Rcpp::as<Rcpp::NumericVector>(SEBinst["RAsoil"]))/((double) ntimesteps);
      
      Tatm_min[i] = MinTemperature[i];
      Tatm_max[i] = MaxTemperature[i];
      Tatm_mean[i] = MeanTemperature[i];
      Tcan_min[i] = min(Tcan);
      Tcan_max[i] = max(Tcan);
      Tcan_mean[i] = sum(Tcan)/((double) ntimesteps);
      Tsoil_min[i] = min(Tsoil);
      Tsoil_max[i] = max(Tsoil);
      Tsoil_mean[i] = sum(Tsoil)/((double) ntimesteps);
      HydrIndays(i,_) = HydrInVec;
      LeafPsiMin(i,_) = Rcpp::as<Rcpp::NumericVector>(Plants["LeafPsiMin"]);
      LeafPsiMax(i,_) = Rcpp::as<Rcpp::NumericVector>(Plants["LeafPsiMax"]);
      LeafPsiMin_SL(i,_) = Rcpp::as<Rcpp::NumericVector>(Plants["LeafPsiMin_SL"]);
      LeafPsiMax_SL(i,_) = Rcpp::as<Rcpp::NumericVector>(Plants["LeafPsiMax_SL"]);
      LeafPsiMin_SH(i,_) = Rcpp::as<Rcpp::NumericVector>(Plants["LeafPsiMin_SH"]);
      LeafPsiMax_SH(i,_) = Rcpp::as<Rcpp::NumericVector>(Plants["LeafPsiMax_SH"]);
      RootPsi(i,_) = Rcpp::as<Rcpp::NumericVector>(Plants["RootPsi"]); 
      StemPsi(i,_) = Rcpp::as<Rcpp::NumericVector>(Plants["StemPsi"]); 
      StemPLC(i,_) = Rcpp::as<Rcpp::NumericVector>(Plants["StemPLC"]); 
      PlantWaterBalance(i,_) = Rcpp::as<Rcpp::NumericVector>(Plants["WaterBalance"]); 
      dEdP(i,_) = Rcpp::as<Rcpp::NumericVector>(Plants["dEdP"]); 
      NumericMatrix RhizoPsiStep = Rcpp::as<Rcpp::NumericMatrix>(s["RhizoPsi"]);
      for(int c=0;c<numCohorts;c++) {
        NumericMatrix nm = Rcpp::as<Rcpp::NumericMatrix>(RhizoPsi[c]);
        nm(i,_) =  RhizoPsiStep(c,_);
      }
      StemRWC(i,_) = as<Rcpp::NumericVector>(Plants["StemRWC"]);
      LeafRWC(i,_) = as<Rcpp::NumericVector>(Plants["LeafRWC"]); 
    } 
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
    if(transpirationMode =="Sperry") {
      Rcout<<" Plant extraction from soil (mm) " << round(sum(PlantExtraction));
      Rcout<<" Hydraulic redistribution (mm) " << round(sum(HydraulicRedistribution)) <<"\n";
    } else {
      Rcout <<"\n";
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
  SWB.attr("row.names") = meteo.attr("row.names") ;
  
  DataFrame DWB;
  if(transpirationMode=="Granier") {
    DWB = DataFrame::create(_["PlantExtraction"] = PlantExtraction, _["Transpiration"]=Transpiration);
  } else {
    DWB = DataFrame::create(_["PlantExtraction"] = PlantExtraction, _["Transpiration"]=Transpiration, 
                            _["HydraulicRedistribution"] = HydraulicRedistribution);
  }
  DWB.attr("row.names") = meteo.attr("row.names") ;
  
  if(transpirationMode=="Granier") PlantAbsSWRFraction.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")); 
  PlantTranspiration.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names"));
  PlantStress.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  StemPLC.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  StemRWC.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  LeafRWC.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  PlantPsi.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  dEdP.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  LeafPsiMin.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  LeafPsiMax.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  LeafPsiMin_SL.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  LeafPsiMax_SL.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  LeafPsiMin_SH.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  LeafPsiMax_SH.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  StemPsi.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  PlantWaterBalance.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  RootPsi.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  PlantPhotosynthesis.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  PlantGrossPhotosynthesis.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  PlantNetPhotosynthesis.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  PlantAbsSWR.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  PlantAbsLWR.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  PlantLAI.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  
  DataFrame DEB = DataFrame::create(_["SWRcanin"] = SWRcanin, _["LWRcanin"] = LWRcanin, _["LWRcanout"] = LWRcanout,
                                    _["LEcan"] = LEcan_heat, _["Hcan"] = Hcan_heat, _["LWRsoilcan"] = LWRsoilcan, _["Ebalcan"] = Ebalcan, 
                                      _["Hcansoil"] = Hcansoil, _["SWRsoilin"] = SWRsoilin, _["LWRsoilin"] = LWRsoilin, _["LWRsoilout"] = LWRsoilout,
                                        _["LEsoil"] = LEsoil_heat, _["Ebalsoil"] = Ebalsoil, _["RAcan"] = RAcan, _["RAsoil"] = RAsoil);  
  DataFrame DT = DataFrame::create(_["Tatm_mean"] = Tatm_mean, _["Tatm_min"] = Tatm_min, _["Tatm_max"] = Tatm_max,
                                   _["Tcan_mean"] = Tcan_mean, _["Tcan_min"] = Tcan_min, _["Tcan_max"] = Tcan_max,
                                     _["Tsoil_mean"] = Tsoil_mean, _["Tsoil_min"] = Tsoil_min, _["Tsoil_max"] = Tsoil_max);
  DEB.attr("row.names") = meteo.attr("row.names") ;
  DT.attr("row.names") = meteo.attr("row.names") ;
  subdailyRes.attr("names") = meteo.attr("row.names") ;
  
  NumericVector topo = NumericVector::create(elevation, slope, aspect);
  topo.attr("names") = CharacterVector::create("elevation", "slope", "aspect");

  List l;
  if(transpirationMode=="Granier") {
    l = List::create(Named("latitude") = latitude,
                     Named("topography") = topo,
                     Named("spwbInput") = spwbInput,
                     Named("soilInput") = soilInput,
                     Named("WaterBalance")=DWB, 
                     Named("Soil") = SWB,
                     Named("PlantLAI") = PlantLAI,
                     Named("PlantAbsorbedSWRFraction") = PlantAbsSWRFraction,
                     Named("PlantTranspiration") = PlantTranspiration,
                     Named("PlantPhotosynthesis") = PlantPhotosynthesis,
                     Named("PlantPsi") = PlantPsi, 
                     Named("PlantStress") = PlantStress,
                     Named("subdaily") =  subdailyRes);
  } else {
    CharacterVector ln = CharacterVector(30);
    l = List(30);
    l[0] = latitude;
    ln[0] = "latitude";
    l[1] = topo;
    ln[1] = "topography";
    l[2] = spwbInput;
    ln[2] = "spwbInput";
    l[3] = soilInput;
    ln[3] = "soilInput";
    l[4] = DWB;
    ln[4] = "WaterBalance";
    l[5] = SWB;
    ln[5] = "Soil";
    l[6] = DEB;
    ln[6] = "EnergyBalance";
    l[7] = DT;
    ln[7] = "Temperature";
    l[8] = PlantLAI;
    ln[8] = "PlantLAI";
    l[9] = PlantAbsSWR;
    ln[9] = "PlantAbsorbedSWR";
    l[10] = PlantAbsLWR;
    ln[10] = "PlantAbsorbedLWR";
    l[11] = PlantTranspiration;
    ln[11] = "PlantTranspiration";
    l[12] = PlantGrossPhotosynthesis;
    ln[12] = "PlantGrossPhotosynthesis";
    l[13] = PlantNetPhotosynthesis;
    ln[13] = "PlantNetPhotosynthesis";
    l[14] = dEdP;
    ln[14] = "dEdP";
    l[15] = LeafPsiMin;
    ln[15] = "LeafPsiMin";
    l[16] = LeafPsiMax;
    ln[16] = "LeafPsiMax";
    l[17] = LeafPsiMin_SL;
    ln[17] = "LeafPsiMin_SL";
    l[18] = LeafPsiMax_SL;
    ln[18] = "LeafPsiMax_SL";
    l[19] = LeafPsiMin_SH;
    ln[19] = "LeafPsiMin_SH";
    l[20] = LeafPsiMax_SH;
    ln[20] = "LeafPsiMax_SH";
    l[21] = LeafRWC;
    ln[21] = "LeafRWC";
    l[22] = StemPsi;
    ln[22] = "StemPsi";
    l[23] = StemPLC;
    ln[23] = "StemPLC";
    l[24] = StemRWC;
    ln[24] = "StemRWC";
    l[25] = RootPsi;
    ln[25] = "RootPsi";
    l[26] = RhizoPsi;
    ln[26] = "RhizoPsi";
    l[27] = PlantStress;
    ln[27] = "PlantStress";
    l[28] = PlantWaterBalance;
    ln[28] = "PlantWaterBalance";
    l[29] = subdailyRes;
    ln[29] = "subdaily";
    l.attr("names") = ln;
  }
  l.attr("class") = CharacterVector::create("pwb","list");
  return(l);                    
}