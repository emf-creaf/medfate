// [[Rcpp::interfaces(r,cpp)]]
#define STRICT_R_HEADERS
#include <Rcpp.h>
#include <numeric>
#include <math.h>
#include "spwb_day.h"
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

// [[Rcpp::export(".getWeatherDates")]]
CharacterVector getWeatherDates(DataFrame meteo){
  CharacterVector dateStrings;
  String dateColumnName = NA_STRING;
  bool is_date_column = false;
  if(meteo.containsElementNamed("dates")) {
    dateColumnName = "dates";
    is_date_column = true;
  } else if(meteo.containsElementNamed("date")) {
    dateColumnName = "date";
    is_date_column = true;
  } else if(meteo.containsElementNamed("Dates")) {
    dateColumnName = "Dates";
    is_date_column = true;
  } else if(meteo.containsElementNamed("Date")) {
    dateColumnName = "Date";
    is_date_column = true;
  }
  if(is_date_column){
    RObject vector = Rcpp::as<Rcpp::RObject>(meteo[dateColumnName]);
    if(is<DateVector>(vector)) {
      DateVector dateVector = Rcpp::as<Rcpp::DateVector>(vector);
      CharacterVector dS(dateVector.size(), NA_STRING);
      for(int i=0;i< dateVector.size();i++) {
        Date d = dateVector[i];
        dS[i] = d.format("%Y-%m-%d");
      }
      dateStrings = dS;
    } else if(vector.inherits("Date")) {
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
    } else if(vector.inherits("POSIXct")) {
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
    } else {
      stop("Could not parse date column.");
    }
  } else {
    dateStrings = meteo.attr("row.names"); 
  }
  return(dateStrings);
}



IntegerVector order_vector(NumericVector x) {
  if (is_true(any(duplicated(x)))) {
    Rf_warning("There are duplicates in 'x'; order not guaranteed to match that of R's base::order");
  }
  NumericVector sorted = clone(x).sort();
  return match(sorted, x);
}



void checkspwbInput(List x,  String transpirationMode, String soilFunctions) {
  
  DataFrame soil = Rcpp::as<Rcpp::DataFrame>(x["soil"]);
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
  if(!soil.containsElementNamed("widths")) stop("widths missing in soil");
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

DataFrame defineStandDailyOutput(CharacterVector dateStrings) {
  int numDays = dateStrings.length();
  NumericVector LAI(numDays), LAIherb(numDays), LAIexpanded(numDays),LAIlive(numDays),LAIdead(numDays);
  NumericVector Cm(numDays);
  NumericVector LgroundPAR(numDays);
  NumericVector LgroundSWR(numDays);
  DataFrame Stand = DataFrame::create(_["LAI"]=LAI, _["LAIherb"]=LAIherb, 
                                      _["LAIlive"]=LAIlive, _["LAIexpanded"] = LAIexpanded, _["LAIdead"] = LAIdead,  
                                      _["Cm"]=Cm, _["LgroundPAR"] = LgroundPAR, _["LgroundSWR"] = LgroundSWR);
  Stand.attr("row.names") = dateStrings;
  return(Stand);
}

DataFrame defineWaterBalanceDailyOutput(CharacterVector dateStrings, String transpirationMode) {
  int numDays = dateStrings.length();
  
  NumericVector PET(numDays), Precipitation(numDays), Evapotranspiration(numDays);
  NumericVector Runoff(numDays),Rain(numDays),Snow(numDays);
  NumericVector Snowmelt(numDays),NetRain(numDays);
  NumericVector Interception(numDays),Infiltration(numDays), InfiltrationExcess(numDays), DeepDrainage(numDays), SaturationExcess(numDays), CapillarityRise(numDays);
  NumericVector SoilEvaporation(numDays), HerbTranspiration(numDays), Transpiration(numDays),PlantExtraction(numDays);
  NumericVector HydraulicRedistribution(numDays, 0.0);
  
  DataFrame DWB = DataFrame::create(_["PET"]=PET, 
                          _["Precipitation"] = Precipitation, _["Rain"] = Rain, _["Snow"] = Snow, 
                          _["NetRain"]=NetRain, _["Snowmelt"] = Snowmelt, _["Infiltration"]=Infiltration,
                          _["InfiltrationExcess"] = InfiltrationExcess, _["SaturationExcess"] = SaturationExcess, _["Runoff"]=Runoff, 
                          _["DeepDrainage"]=DeepDrainage, _["CapillarityRise"] = CapillarityRise,
                          _["Evapotranspiration"]=Evapotranspiration,_["Interception"] = Interception, 
                          _["SoilEvaporation"]=SoilEvaporation, _["HerbTranspiration"] = HerbTranspiration,
                          _["PlantExtraction"] = PlantExtraction, _["Transpiration"]=Transpiration, 
                          _["HydraulicRedistribution"] = HydraulicRedistribution);
  DWB.attr("row.names") = dateStrings;
  return(DWB);
}
DataFrame defineSnowDailyOutput(CharacterVector dateStrings) {
  int numDays = dateStrings.length();
  NumericVector SWE(numDays, 0.0);
  DataFrame Snow = DataFrame::create(_["SWE"] = SWE);
  Snow.attr("row.names") = dateStrings;
  return(Snow);
}
List defineSoilDailyOutput(CharacterVector dateStrings, DataFrame soil, bool includePlants = true) {
  int numDays = dateStrings.length();
  NumericVector W = soil["W"];
  int nlayers = W.length();
  
  CharacterVector layerNames(nlayers + 1);
  for(int l=0; l<nlayers; l++) {
    String s = "";
    s += (l+1);
    layerNames[l] = s;
  }
  layerNames[nlayers] = "Overall";
  
  NumericMatrix SWCdays(numDays, nlayers+1); //Soil moisture content in percent volume
  SWCdays.attr("dimnames") = List::create(dateStrings, layerNames);
  NumericMatrix RWCdays(numDays, nlayers+1); //Soil moisture content in relation to field capacity
  RWCdays.attr("dimnames") = List::create(dateStrings, layerNames);
  NumericMatrix REWdays(numDays, nlayers+1); //Soil moisture content in relation to extractable water
  REWdays.attr("dimnames") = List::create(dateStrings, layerNames);
  NumericMatrix Psidays(numDays, nlayers+1);
  Psidays.attr("dimnames") = List::create(dateStrings, layerNames);
  NumericMatrix MLdays(numDays, nlayers+1);
  MLdays.attr("dimnames") = List::create(dateStrings, layerNames);
  
  List Soil = List::create(_["SWC"]=SWCdays, 
                           _["RWC"]=RWCdays, 
                           _["REW"]=REWdays, 
                           _["ML"]=MLdays,
                           _["Psi"]=Psidays);
  if(includePlants) {
    NumericMatrix Eplantdays(numDays, nlayers+1);
    Eplantdays.attr("dimnames") = List::create(dateStrings, layerNames);
    NumericMatrix HydrIndays(numDays, nlayers+1);
    HydrIndays.attr("dimnames") = List::create(dateStrings, layerNames);
    Soil.push_back(Eplantdays , "PlantExt");
    Soil.push_back(HydrIndays , "HydraulicInput");
  }
  return(Soil);  
}
List defineSoilPoolDailyOutput(CharacterVector dateStrings, DataFrame above, DataFrame soil) {
  int numDays = dateStrings.length();
  NumericVector W = soil["W"];
  int nlayers = W.length();
  int numCohorts = above.nrow();
  
  List REWpool(numCohorts);
  for(int c=0;c<numCohorts;c++) {
    NumericMatrix nm = NumericMatrix(numDays, nlayers);
    nm.attr("dimnames") = List::create(dateStrings, seq(1,nlayers)) ;
    REWpool[c] = nm;
  }
  REWpool.attr("names") = above.attr("row.names");

  List Psipool(numCohorts);
  for(int c=0;c<numCohorts;c++) {
    NumericMatrix nm = NumericMatrix(numDays, nlayers);
    nm.attr("dimnames") = List::create(dateStrings, seq(1,nlayers)) ;
    Psipool[c] = nm;
  }
  Psipool.attr("names") = above.attr("row.names");
  List soilPool = List::create(_["REW"] = REWpool,
                               _["Psi"] = Psipool);
  return(soilPool);
}
DataFrame defineEnergyBalanceDailyOutput(CharacterVector dateStrings) {
  int numDays = dateStrings.length();
  NumericVector SWRcan(numDays, NA_REAL);
  NumericVector LWRcan(numDays, NA_REAL);
  NumericVector LEVcan(numDays, NA_REAL);
  NumericVector LEVsoil(numDays, NA_REAL);
  NumericVector LEFsnow(numDays, NA_REAL);
  NumericVector Hcan_heat(numDays, NA_REAL);
  NumericVector Ebalcan(numDays, NA_REAL);
  NumericVector SWRsoil(numDays, NA_REAL);
  NumericVector LWRsoil(numDays, NA_REAL);
  NumericVector Ebalsoil(numDays, NA_REAL);
  NumericVector Hcansoil(numDays, NA_REAL);

  DataFrame DEB = DataFrame::create(_["SWRcan"] = SWRcan, _["LWRcan"] = LWRcan,
                                    _["LEVcan"] = LEVcan, _["LEFsnow"] = LEFsnow, _["Hcan"] = Hcan_heat, _["Ebalcan"] = Ebalcan, 
                                    _["Hcansoil"] = Hcansoil, _["SWRsoil"] = SWRsoil, _["LWRsoil"] = LWRsoil, 
                                    _["LEVsoil"] = LEVsoil, _["Ebalsoil"] = Ebalsoil);  
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

List definePlantWaterDailyOutput(CharacterVector dateStrings, DataFrame above, DataFrame soil, List control) {
  
  String transpirationMode = control["transpirationMode"];
  int numDays = dateStrings.length();
  NumericVector W = soil["W"];
  int nlayers = W.length();
  int numCohorts = above.nrow();
 
  NumericMatrix PlantStress(numDays, numCohorts);
  NumericMatrix PlantTranspiration(numDays, numCohorts);
  NumericMatrix PlantLAI(numDays, numCohorts);
  NumericMatrix PlantLAIlive(numDays, numCohorts);
  NumericMatrix LeafPLC(numDays, numCohorts);
  NumericMatrix StemPLC(numDays, numCohorts);
  NumericMatrix StemRWC(numDays, numCohorts), LeafRWC(numDays, numCohorts), LFMC(numDays, numCohorts);
  NumericMatrix PlantWaterBalance(numDays, numCohorts);
  NumericMatrix PlantFPAR(numDays, numCohorts);
  
  PlantFPAR.attr("dimnames") = List::create(dateStrings, above.attr("row.names")) ;
  PlantTranspiration.attr("dimnames") = List::create(dateStrings, above.attr("row.names"));
  PlantStress.attr("dimnames") = List::create(dateStrings, above.attr("row.names")) ;
  
  PlantLAI.attr("dimnames") = List::create(dateStrings, above.attr("row.names")) ;
  PlantLAIlive.attr("dimnames") = List::create(dateStrings, above.attr("row.names")) ;
  PlantTranspiration.attr("dimnames") = List::create(dateStrings, above.attr("row.names"));
  PlantStress.attr("dimnames") = List::create(dateStrings, above.attr("row.names")) ;
  LeafPLC.attr("dimnames") = List::create(dateStrings, above.attr("row.names")) ;
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
                               Named("LeafPLC") = LeafPLC,
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
                          Named("FPAR") = PlantFPAR,
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
                               Named("LeafPLC") = LeafPLC,
                               Named("StemPLC") = StemPLC, 
                               Named("RootPsi") = RootPsi, 
                               Named("RhizoPsi") = RhizoPsi, 
                               Named("LFMC") = LFMC);
    plants.push_back(PlantStress, "PlantStress");
    
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

// [[Rcpp::export(".defineSPWBDailyOutput")]]
List defineSPWBDailyOutput(double latitude, double elevation, double slope, double aspect, 
                           CharacterVector dateStrings, List x) {
  
  NumericVector topo = NumericVector::create(elevation, slope, aspect);
  topo.attr("names") = CharacterVector::create("elevation", "slope", "aspect");
  
  List spwbInput = clone(x);
  List control = x["control"];
  DataFrame soil = Rcpp::as<Rcpp::DataFrame>(x["soil"]);
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  String transpirationMode = control["transpirationMode"];
  String rhizosphereOverlap = control["rhizosphereOverlap"];
  bool plantWaterPools = (rhizosphereOverlap!="total");
  DataFrame DWB = defineWaterBalanceDailyOutput(dateStrings, transpirationMode);
  List Soil = defineSoilDailyOutput(dateStrings, soil, true);
  List soilPool = defineSoilPoolDailyOutput(dateStrings, above, soil);
  DataFrame Snow = defineSnowDailyOutput(dateStrings);
  List sunlitDO = defineSunlitShadeLeavesDailyOutput(dateStrings, above);
  List shadeDO = defineSunlitShadeLeavesDailyOutput(dateStrings, above);
  List plantDWOL = definePlantWaterDailyOutput(dateStrings, above, soil, control);
  DataFrame Stand = defineStandDailyOutput(dateStrings);
  //Detailed subdaily results
  int numDays = dateStrings.size();
  List subdailyRes(numDays);
  subdailyRes.attr("names") = dateStrings;
  
  List l;
  if(transpirationMode=="Granier") {
    l = List::create(Named("latitude") = latitude,
                     Named("topography") = topo,
                     Named("weather") = NA_REAL,
                     Named("spwbInput") = spwbInput,
                     Named("spwbOutput") = x,
                     Named("WaterBalance")= DWB);
    if(control["soilResults"]) l.push_back(Soil, "Soil");
    if(control["soilPoolResults"] && plantWaterPools) l.push_back(soilPool, "SoilPools");
    if(control["snowResults"]) l.push_back(Snow, "Snow");
    if(control["standResults"]) l.push_back(Stand, "Stand");
    if(control["plantResults"]) l.push_back(plantDWOL, "Plants");
    if(control["fireHazardResults"]) {
      DataFrame fireHazard = defineFireHazardOutput(dateStrings);
      l.push_back(fireHazard, "FireHazard");
    }
  } else {
    DataFrame canopy = Rcpp::as<Rcpp::DataFrame>(x["canopy"]);
    DataFrame DEB = defineEnergyBalanceDailyOutput(dateStrings);
    DataFrame DT = defineTemperatureDailyOutput(dateStrings);
    NumericMatrix DLT =  defineTemperatureLayersDailyOutput(dateStrings, canopy);
    l = List::create(Named("latitude") = latitude,
                     Named("topography") = topo,
                     Named("weather") = NA_REAL,
                     Named("spwbInput") = spwbInput,
                     Named("spwbOutput") = x,
                     Named("WaterBalance")= DWB,
                     Named("EnergyBalance") = DEB);
      if(control["temperatureResults"]) {
        l.push_back(DT, "Temperature");
        if(control["multiLayerBalance"]) l.push_back(DLT,"TemperatureLayers");
      }
    if(control["soilResults"]) l.push_back(Soil, "Soil");
    if(control["soilPoolResults"] && plantWaterPools) l.push_back(soilPool, "SoilPools");
    if(control["snowResults"]) l.push_back(Snow, "Snow");
    if(control["standResults"]) l.push_back(Stand, "Stand");
    if(control["plantResults"]) l.push_back(plantDWOL, "Plants");
    if(control["leafResults"]) {
      l.push_back(sunlitDO, "SunlitLeaves");
      l.push_back(shadeDO, "ShadeLeaves");
    }
    if(control["fireHazardResults"]) {
      DataFrame fireHazard = defineFireHazardOutput(dateStrings);
      l.push_back(fireHazard, "FireHazard");
    }
  }
  if(control["subdailyResults"]) l.push_back(subdailyRes,"subdaily");
  l.attr("class") = CharacterVector::create("spwb","list");
  return(l);
}

void fillWaterBalanceDailyOutput(DataFrame DWB, List sDay, int iday, String transpirationMode) {
  List db = sDay["WaterBalance"];
  NumericVector PET = DWB["PET"];
  NumericVector Precipitation = DWB["Precipitation"];
  NumericVector DeepDrainage = DWB["DeepDrainage"];
  NumericVector Infiltration = DWB["Infiltration"];
  NumericVector InfiltrationExcess = DWB["InfiltrationExcess"];
  NumericVector SaturationExcess = DWB["SaturationExcess"];
  NumericVector CapillarityRise = DWB["CapillarityRise"];
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
  InfiltrationExcess[iday] = db["InfiltrationExcess"];
  SaturationExcess[iday] = db["SaturationExcess"];
  CapillarityRise[iday] = db["CapillarityRise"];
  Runoff[iday] = db["Runoff"];
  Rain[iday] = db["Rain"];
  Snow[iday] = db["Snow"];
  Precipitation[iday] = Rain[iday]+Snow[iday];
  PET[iday] = db["PET"];
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

void fillStandDailyOutput(DataFrame Stand, List sDay, int iday) {
  List stand = sDay["Stand"];
  NumericVector LgroundPAR = Stand["LgroundPAR"];
  NumericVector LgroundSWR = Stand["LgroundSWR"];
  NumericVector LAI = Stand["LAI"];
  NumericVector LAIherb = Stand["LAIherb"];
  NumericVector LAIexpanded = Stand["LAIexpanded"];
  NumericVector LAIlive = Stand["LAIlive"];
  NumericVector LAIdead = Stand["LAIdead"];
  NumericVector Cm = Stand["Cm"];
  
  LgroundPAR[iday] = stand["LgroundPAR"];
  LgroundSWR[iday] = stand["LgroundSWR"];
  LAI[iday] = stand["LAI"];
  LAIherb[iday] = stand["LAIherb"];
  LAIexpanded[iday] = stand["LAIexpanded"];
  LAIlive[iday] = stand["LAIlive"];
  LAIdead[iday] = stand["LAIdead"];
  Cm[iday] = stand["Cm"];
}

void fillSoilDailyOutput(List SWB, DataFrame soil, List sDay, 
                         int iday, int numDays, String soilFunctions,
                         bool includePlants = true) {
  int nlayers = soil.nrow();
  NumericVector W = soil["W"];
  NumericVector Water_SAT = waterSAT(soil, soilFunctions);
  NumericVector Theta_FC = thetaFC(soil, soilFunctions);
  NumericVector Water_FC = waterFC(soil, soilFunctions);
  NumericVector Water_min = waterPsi(soil, -5.0, soilFunctions);
  NumericVector Water_ext = waterExtractable(soil, soilFunctions, -5.0);
  
  List sb = sDay["Soil"];
  NumericVector psi = sb["Psi"];
  
  NumericMatrix SWCdays = as<Rcpp::NumericMatrix>(SWB["SWC"]);
  NumericMatrix MLdays = as<Rcpp::NumericMatrix>(SWB["ML"]);
  NumericMatrix RWCdays = as<Rcpp::NumericMatrix>(SWB["RWC"]);
  NumericMatrix Psidays = as<Rcpp::NumericMatrix>(SWB["Psi"]);
  NumericMatrix REWdays = as<Rcpp::NumericMatrix>(SWB["REW"]);

  for(int l=0; l<nlayers; l++) {
    Psidays(iday,l) = psi[l];
    SWCdays(iday,l) = W[l]*Theta_FC[l];
    RWCdays(iday,l) = W[l];
    MLdays(iday,l) = RWCdays(iday,l)*Water_FC[l]; 
    REWdays(iday,l) = (MLdays(iday,l)-Water_min[l])/Water_ext[l];
    MLdays(iday,nlayers) = MLdays(iday,nlayers) + MLdays(iday,l);
    Psidays(iday,nlayers) += psi[l]*Water_SAT[l];
  }
  SWCdays(iday,nlayers) = sum((W*Theta_FC)*Water_SAT)/sum(Water_SAT);
  RWCdays(iday,nlayers) = MLdays(iday,nlayers)/sum(Water_FC);
  REWdays(iday,nlayers) = (MLdays(iday,nlayers) - sum(Water_min))/sum(Water_ext);
  Psidays(iday,nlayers) = Psidays(iday,nlayers)/sum(Water_SAT);

  if(includePlants) {
    NumericMatrix Eplantdays = as<Rcpp::NumericMatrix>(SWB["PlantExt"]);
    NumericMatrix HydrIndays = as<Rcpp::NumericMatrix>(SWB["HydraulicInput"]);
    NumericVector HydrInVec = sb["HydraulicInput"];
    List Plants = sDay["Plants"];
    NumericVector EplantVec = sb["PlantExtraction"];
    Eplantdays(iday,nlayers) = 0.0;
    HydrIndays(iday,nlayers) = 0.0;
    for(int l=0; l<nlayers; l++) {
      HydrIndays(iday,l) = HydrInVec[l];
      Eplantdays(iday,l) = EplantVec[l];
      Eplantdays(iday,nlayers) += EplantVec[l];
      HydrIndays(iday,nlayers) += HydrIndays[l];
    }
  }
}
void fillSoilPoolDailyOutput(List soilPools, DataFrame soil, NumericMatrix Wpool, int iday, String soilFunctions) {
  List REWpool = soilPools["REW"];
  List Psipool = soilPools["Psi"];
  int numCohorts = Psipool.size();
  int nlayers = soil.nrow();
  for(int c=0;c<numCohorts;c++) {
    NumericMatrix rewM = Rcpp::as<Rcpp::NumericMatrix>(REWpool[c]);
    rewM(iday,_) = Wpool(c,_);
    NumericMatrix psiM = Rcpp::as<Rcpp::NumericMatrix>(Psipool[c]);
    DataFrame soil_pool = clone(soil);
    NumericVector soil_pool_W = soil_pool["W"];
    for(int l=0;l<nlayers;l++) soil_pool_W[l] = Wpool(c,l);
    psiM(iday,_) = psi(soil_pool, soilFunctions);
  }
}
void fillSnowDailyOutput(DataFrame Snow, List x, int iday) {
  NumericVector SWE = Snow["SWE"];
  SWE[iday] = x["snowpack"];
}
void fillEnergyBalanceDailyOutput(DataFrame DEB, List sDay, int iday, int ntimesteps) {
  List EB = Rcpp::as<Rcpp::List>(sDay["EnergyBalance"]);
  DataFrame Tinst = Rcpp::as<Rcpp::DataFrame>(EB["Temperature"]); 
  DataFrame CEBinst = Rcpp::as<Rcpp::DataFrame>(EB["CanopyEnergyBalance"]); 
  DataFrame SEBinst = Rcpp::as<Rcpp::DataFrame>(EB["SoilEnergyBalance"]); 
  double tstep = 86400.0/((double) ntimesteps);
  
  NumericVector SWRcan = DEB["SWRcan"];
  NumericVector LWRcan = DEB["LWRcan"];
  NumericVector LEVcan = DEB["LEVcan"];
  NumericVector LEFsnow = DEB["LEFsnow"];
  NumericVector Hcan_heat = DEB["Hcan"];
  NumericVector Ebalcan = DEB["Ebalcan"];
  SWRcan[iday] = 0.000001*sum(Rcpp::as<Rcpp::NumericVector>(CEBinst["SWRcan"]))*tstep;
  LWRcan[iday] = 0.000001*sum(Rcpp::as<Rcpp::NumericVector>(CEBinst["LWRcan"]))*tstep;
  LEVcan[iday] = 0.000001*sum(Rcpp::as<Rcpp::NumericVector>(CEBinst["LEVcan"]))*tstep;
  LEFsnow[iday] = 0.000001*sum(Rcpp::as<Rcpp::NumericVector>(CEBinst["LEFsnow"]))*tstep;
  Hcan_heat[iday] = 0.000001*sum(Rcpp::as<Rcpp::NumericVector>(CEBinst["Hcan"]))*tstep;
  Ebalcan[iday] = 0.000001*sum(Rcpp::as<Rcpp::NumericVector>(CEBinst["Ebalcan"]))*tstep;
  NumericVector SWRsoil = DEB["SWRsoil"];
  NumericVector LWRsoil = DEB["LWRsoil"];
  NumericVector LEVsoil = DEB["LEVsoil"];
  NumericVector Hcansoil = DEB["Hcansoil"];
  NumericVector Ebalsoil = DEB["Ebalsoil"];
  SWRsoil[iday] = 0.000001*sum(Rcpp::as<Rcpp::NumericVector>(SEBinst["SWRsoil"]))*tstep;
  LWRsoil[iday] = 0.000001*sum(Rcpp::as<Rcpp::NumericVector>(SEBinst["LWRsoil"]))*tstep;
  LEVsoil[iday] = 0.000001*sum(Rcpp::as<Rcpp::NumericVector>(SEBinst["LEVsoil"]))*tstep;
  Hcansoil[iday] = 0.000001*sum(Rcpp::as<Rcpp::NumericVector>(SEBinst["Hcansoil"]))*tstep;
  Ebalsoil[iday] = 0.000001*sum(Rcpp::as<Rcpp::NumericVector>(SEBinst["Ebalsoil"]))*tstep;
}

void fillTemperatureDailyOutput(DataFrame DT, List sDay, int iday, int ntimesteps) {
  List EB = Rcpp::as<Rcpp::List>(sDay["EnergyBalance"]);
  DataFrame Tinst = Rcpp::as<Rcpp::DataFrame>(EB["Temperature"]); 
  NumericVector Tatm = Rcpp::as<Rcpp::NumericVector>(Tinst["Tatm"]);
  NumericVector Tcan = Rcpp::as<Rcpp::NumericVector>(Tinst["Tcan"]);
  NumericMatrix TsoilMat = Rcpp::as<Rcpp::NumericMatrix>(EB["SoilTemperature"]); 
  NumericVector Tsoil_0(ntimesteps);
  for(int n = 0;n<ntimesteps; n++) Tsoil_0[n] = TsoilMat(n,0);

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
  Tsoil_min[iday] = min(Tsoil_0);
  Tsoil_max[iday] = max(Tsoil_0);
  Tsoil_mean[iday] = sum(Tsoil_0)/((double) ntimesteps);
}

void fillTemperatureLayersDailyOutput(NumericMatrix DLT, List sDay, int iday, int ncanlayers, int ntimesteps) {
  List EB = Rcpp::as<Rcpp::List>(sDay["EnergyBalance"]);
  DataFrame Tinst = Rcpp::as<Rcpp::DataFrame>(EB["Temperature"]); 
  NumericVector Tcan = Rcpp::as<Rcpp::NumericVector>(Tinst["Tcan"]);

  NumericMatrix LTinst = Rcpp::as<Rcpp::NumericMatrix>(EB["TemperatureLayers"]); 
  for(int l=0;l<ncanlayers;l++) {
    DLT(iday, l) = 0.0;
    for(int n=0;n<ntimesteps;n++) {
       DLT(iday, l) += LTinst(n,l);
    }
    DLT(iday, l) = DLT(iday, l)/((double) ntimesteps);
  }
}

void fillPlantWaterDailyOutput(List x, List sDay, int iday, String transpirationMode) {
  List Plants = sDay["Plants"];
  
  NumericMatrix PlantStress= Rcpp::as<Rcpp::NumericMatrix>(x["PlantStress"]);
  NumericMatrix PlantTranspiration= Rcpp::as<Rcpp::NumericMatrix>(x["Transpiration"]);
  NumericMatrix PlantLAI= Rcpp::as<Rcpp::NumericMatrix>(x["LAI"]);
  NumericMatrix PlantLAIlive= Rcpp::as<Rcpp::NumericMatrix>(x["LAIlive"]);
  NumericMatrix LeafPLC= Rcpp::as<Rcpp::NumericMatrix>(x["LeafPLC"]);
  NumericMatrix StemPLC= Rcpp::as<Rcpp::NumericMatrix>(x["StemPLC"]);
  NumericMatrix StemRWC= Rcpp::as<Rcpp::NumericMatrix>(x["StemRWC"]);
  NumericMatrix LeafRWC= Rcpp::as<Rcpp::NumericMatrix>(x["LeafRWC"]);
  NumericMatrix LFMC= Rcpp::as<Rcpp::NumericMatrix>(x["LFMC"]);
  NumericMatrix PlantWaterBalance= Rcpp::as<Rcpp::NumericMatrix>(x["PlantWaterBalance"]);
  NumericMatrix PlantFPAR= Rcpp::as<Rcpp::NumericMatrix>(x["FPAR"]);
  
  int numCohorts = PlantLAI.ncol();
  
  PlantTranspiration(iday,_) =  Rcpp::as<Rcpp::NumericVector>(Plants["Transpiration"]);
  PlantStress(iday,_) = Rcpp::as<Rcpp::NumericVector>(Plants["DDS"]);
  PlantLAI(iday,_) = Rcpp::as<Rcpp::NumericVector>(Plants["LAI"]);
  PlantLAIlive(iday,_) = Rcpp::as<Rcpp::NumericVector>(Plants["LAIlive"]);
  LeafPLC(iday,_) = Rcpp::as<Rcpp::NumericVector>(Plants["LeafPLC"]); 
  StemPLC(iday,_) = Rcpp::as<Rcpp::NumericVector>(Plants["StemPLC"]); 
  StemRWC(iday,_) = as<Rcpp::NumericVector>(Plants["StemRWC"]);
  LeafRWC(iday,_) = as<Rcpp::NumericVector>(Plants["LeafRWC"]); 
  LFMC(iday,_) = as<Rcpp::NumericVector>(Plants["LFMC"]); 
  PlantWaterBalance(iday,_) = Rcpp::as<Rcpp::NumericVector>(Plants["WaterBalance"]); 
  PlantFPAR(iday,_) =  Rcpp::as<Rcpp::NumericVector>(Plants["FPAR"]);  
  
  
  if(transpirationMode=="Granier") {
    NumericMatrix PlantGrossPhotosynthesis= Rcpp::as<Rcpp::NumericMatrix>(x["GrossPhotosynthesis"]);
    NumericMatrix PlantPsi = Rcpp::as<Rcpp::NumericMatrix>(x["PlantPsi"]);
    NumericMatrix PlantAbsSWRFraction= Rcpp::as<Rcpp::NumericMatrix>(x["AbsorbedSWRFraction"]);
    
    PlantPsi(iday,_) =  Rcpp::as<Rcpp::NumericVector>(Plants["PlantPsi"]);  
    PlantGrossPhotosynthesis(iday,_) =  Rcpp::as<Rcpp::NumericVector>(Plants["GrossPhotosynthesis"]);  
    PlantAbsSWRFraction(iday,_) =  Rcpp::as<Rcpp::NumericVector>(Plants["AbsorbedSWRFraction"]);  
  } else {
    NumericMatrix RhizoPsiStep = Rcpp::as<Rcpp::NumericMatrix>(sDay["RhizoPsi"]);
    List PlantsInst = sDay["PlantsInst"];
    
    NumericMatrix dEdP = Rcpp::as<Rcpp::NumericMatrix>(x["dEdP"]);
    NumericMatrix LeafPsiMin = Rcpp::as<Rcpp::NumericMatrix>(x["LeafPsiMin"]);
    NumericMatrix LeafPsiMax = Rcpp::as<Rcpp::NumericMatrix>(x["LeafPsiMax"]);
    NumericMatrix StemPsi= Rcpp::as<Rcpp::NumericMatrix>(x["StemPsi"]);
    NumericMatrix RootPsi= Rcpp::as<Rcpp::NumericMatrix>(x["RootPsi"]);

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

void fillSunlitShadeLeavesDailyOutput(List sunlit, List shade, List sDay, int iday, int numCohorts) {

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

  DataFrame SunlitLeaves = as<DataFrame>(sDay["SunlitLeaves"]);
  DataFrame ShadeLeaves = as<DataFrame>(sDay["ShadeLeaves"]); 
  NumericVector GSWMin_SL = Rcpp::as<Rcpp::NumericVector>(SunlitLeaves["GSWMin"]);
  NumericVector GSWMin_SH = Rcpp::as<Rcpp::NumericVector>(ShadeLeaves["GSWMin"]);
  NumericVector GSWMax_SL = Rcpp::as<Rcpp::NumericVector>(SunlitLeaves["GSWMax"]);
  NumericVector GSWMax_SH = Rcpp::as<Rcpp::NumericVector>(ShadeLeaves["GSWMax"]);
  NumericVector PsiMin_SL = Rcpp::as<Rcpp::NumericVector>(SunlitLeaves["LeafPsiMin"]);
  NumericVector PsiMin_SH =  Rcpp::as<Rcpp::NumericVector>(SunlitLeaves["LeafPsiMin"]);
  NumericVector PsiMax_SL = Rcpp::as<Rcpp::NumericVector>(ShadeLeaves["LeafPsiMax"]);
  NumericVector PsiMax_SH = Rcpp::as<Rcpp::NumericVector>(ShadeLeaves["LeafPsiMax"]);
  NumericVector TempMin_SL = Rcpp::as<Rcpp::NumericVector>(SunlitLeaves["TempMin"]);
  NumericVector TempMin_SH = Rcpp::as<Rcpp::NumericVector>(ShadeLeaves["TempMin"]);
  NumericVector TempMax_SL = Rcpp::as<Rcpp::NumericVector>(SunlitLeaves["TempMax"]);
  NumericVector TempMax_SH = Rcpp::as<Rcpp::NumericVector>(ShadeLeaves["TempMax"]);
  
  for(int i=0;i<numCohorts;i++) {
    LeafGSWMin_SL(iday,i) = GSWMin_SL[i];
    LeafGSWMin_SH(iday,i) = GSWMin_SH[i];
    LeafGSWMax_SL(iday,i) = GSWMax_SL[i];
    LeafGSWMax_SH(iday,i) = GSWMax_SH[i];
    LeafPsiMin_SL(iday,i) = PsiMin_SL[i];
    LeafPsiMin_SH(iday,i) = PsiMin_SH[i];
    LeafPsiMax_SL(iday,i) = PsiMax_SL[i];
    LeafPsiMax_SH(iday,i) = PsiMax_SH[i];
    LeafTempMin_SL(iday,i) = TempMin_SL[i];
    LeafTempMin_SH(iday,i) = TempMin_SH[i];
    LeafTempMax_SL(iday,i) = TempMax_SL[i];
    LeafTempMax_SH(iday,i) = TempMax_SH[i];
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

// [[Rcpp::export(".fillSPWBDailyOutput")]]
void fillSPWBDailyOutput(List l, List x, List sDay, int iday) {
  
  List control = x["control"];
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  int numCohorts = above.nrow();
  int ntimesteps = control["ndailysteps"];
  DataFrame canopyParams = Rcpp::as<Rcpp::DataFrame>(x["canopy"]);
  DataFrame soil = Rcpp::as<Rcpp::DataFrame>(x["soil"]);
  int ncanlayers = canopyParams.nrow(); //Number of canopy layers
  String transpirationMode = control["transpirationMode"];
  String rhizosphereOverlap = control["rhizosphereOverlap"];
  String soilFunctions = control["soilFunctions"];
  bool plantWaterPools = (rhizosphereOverlap!="total");
  
  DataFrame DWB = Rcpp::as<Rcpp::DataFrame>(l["WaterBalance"]);
  int numDays = DWB.nrow();
  fillWaterBalanceDailyOutput(DWB, sDay, iday, transpirationMode);
  
  if(control["soilResults"]) {
    List Soil = Rcpp::as<Rcpp::List>(l["Soil"]);
    fillSoilDailyOutput(Soil, soil, sDay, 
                                    iday, numDays, soilFunctions,
                                    true);
  }
  if(control["soilPoolResults"] && plantWaterPools) {
    List soilPools = Rcpp::as<Rcpp::List>(l["SoilPools"]);
    List belowLayers = Rcpp::as<Rcpp::List>(x["belowLayers"]);
    NumericMatrix Wpool = Rcpp::as<Rcpp::NumericMatrix>(belowLayers["Wpool"]);
    fillSoilPoolDailyOutput(soilPools, soil, Wpool, iday, soilFunctions);
  }
  if(control["snowResults"]) {
    DataFrame Snow = Rcpp::as<Rcpp::DataFrame>(l["Snow"]);
    fillSnowDailyOutput(Snow, x, iday);
  }
  if(control["standResults"]) {
    DataFrame Stand = Rcpp::as<Rcpp::DataFrame>(l["Stand"]);
    fillStandDailyOutput(Stand, sDay,iday); 
  }
  if(control["plantResults"]) {
    List plantDWOL = l["Plants"];
    fillPlantWaterDailyOutput(plantDWOL, sDay, iday, transpirationMode); 
    if(transpirationMode!= "Granier") {
      if(control["leafResults"]) {
        List sunlitDO = l["SunlitLeaves"];
        List shadeDO = l["ShadeLeaves"];
        fillSunlitShadeLeavesDailyOutput(sunlitDO, shadeDO, sDay, iday, numCohorts); 
      }
    } 
  }
  if(transpirationMode!= "Granier") {
    List DEB = l["EnergyBalance"];
    fillEnergyBalanceDailyOutput(DEB,sDay, iday, ntimesteps);
    
    if(control["temperatureResults"]) {
      List DT = l["Temperature"];
      fillTemperatureDailyOutput(DT,sDay, iday, ntimesteps);
      if(control["multiLayerBalance"]) {
        NumericMatrix DLT = l["TemperatureLayers"];
        fillTemperatureLayersDailyOutput(DLT,sDay, iday, ncanlayers, ntimesteps);
      }
    }
  } 
  if(control["fireHazardResults"]) {
    DataFrame fireHazard = Rcpp::as<Rcpp::DataFrame>(l["FireHazard"]);
    fillFireHazardOutput(fireHazard, sDay, iday);
  }
  
  if(control["subdailyResults"]) {
    List subdailyRes = Rcpp::as<Rcpp::List>(l["subdaily"]);
    subdailyRes[iday] = copyAdvancedSPWBOutput(sDay, x); //Clones subdaily results because they are communication structures
  }
}

void printWaterBalanceResult(List outputList, List x,
                             NumericVector initialPlantContent, NumericVector initialSoilContent, double initialSnowContent,
                             String transpirationMode) {
  DataFrame soil = Rcpp::as<Rcpp::DataFrame>(x["soil"]);
  List control = x["control"];
  String soilFunctions = control["soilFunctions"];
  bool plantResults = control["plantResults"];

  DataFrame DWB = Rcpp::as<Rcpp::DataFrame>(outputList["WaterBalance"]);
  List plantDWOL;
  if(plantResults)  plantDWOL = outputList["Plants"];
  
  NumericVector finalPlantContent = plantWaterContent(x);
  NumericVector finalSoilContent = water(soil, soilFunctions);
  double finalSnowContent = x["snowpack"];
  if(plantResults) Rcout<<"Final plant water content (mm): "<< sum(finalPlantContent)<<"\n";
  Rcout<<"Final soil water content (mm): "<< sum(finalSoilContent)<<"\n";
  Rcout<<"Final snowpack content (mm): "<< finalSnowContent<<"\n";
  
  NumericVector Precipitation = DWB["Precipitation"];
  NumericVector DeepDrainage = DWB["DeepDrainage"];
  NumericVector Infiltration = DWB["Infiltration"];
  NumericVector InfiltrationExcess = DWB["InfiltrationExcess"];
  NumericVector SaturationExcess = DWB["SaturationExcess"];
  NumericVector CapillarityRise = DWB["CapillarityRise"];
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
  
  NumericMatrix PlantWaterBalance;
  if(plantResults) PlantWaterBalance = Rcpp::as<Rcpp::NumericMatrix>(plantDWOL["PlantWaterBalance"]);
  
  double Precipitationsum = sum(Precipitation);
  double Rainfallsum = sum(Rain);
  double NetRainsum = sum(NetRain);
  double Interceptionsum = sum(Interception);
  double SoilEvaporationsum = sum(SoilEvaporation);
  double Runoffsum  = sum(Runoff);
  double Infiltrationsum  = sum(Infiltration);
  double InfiltrationExcesssum = sum(InfiltrationExcess);
  double SaturationExcesssum  = sum(SaturationExcess);
  double CapillarityRisesum = sum(CapillarityRise);
  double DeepDrainagesum = sum(DeepDrainage);
  double Transpirationsum = sum(Transpiration);
  double Snowmeltsum = sum(Snowmelt);
  double Snowsum = sum(Snow);
  double HerbTranspirationsum = sum(HerbTranspiration);
  
  double soil_wb = Infiltrationsum + CapillarityRisesum - SaturationExcesssum - DeepDrainagesum - SoilEvaporationsum - HerbTranspirationsum - sum(PlantExtraction);
  double snowpack_wb = Snowsum - Snowmeltsum;
  if(plantResults) {
    Rcout<<"Change in plant water content (mm): "<< sum(finalPlantContent) - sum(initialPlantContent)<<"\n";
    Rcout<<"Plant water balance result (mm): "<< sum(PlantWaterBalance)<<"\n"; 
  }
  Rcout<<"Change in soil water content (mm): "<< sum(finalSoilContent) - sum(initialSoilContent)<<"\n";
  Rcout<<"Soil water balance result (mm): "<< soil_wb<<"\n";
  Rcout<<"Change in snowpack water content (mm): "<< finalSnowContent - initialSnowContent<<"\n";
  Rcout<<"Snowpack water balance result (mm): "<< snowpack_wb<<"\n";
  Rcout<<"Water balance components:\n";
  Rcout<<"  Precipitation (mm) "  <<round(Precipitationsum) << " Rain (mm) "  <<round(Rainfallsum) <<" Snow (mm) "  <<round(Snowsum) <<"\n";
  Rcout<<"  Interception (mm) " << round(Interceptionsum)  <<" Net rainfall (mm) " << round(NetRainsum) <<"\n";
  Rcout<<"  Infiltration (mm) " << round(Infiltrationsum)  << " Infiltration excess (mm) " << round(InfiltrationExcesssum) << " Saturation excess (mm) " << round(SaturationExcesssum) << " Capillarity rise (mm) " << round(CapillarityRisesum)  <<"\n";
  Rcout<<"  Soil evaporation (mm) " << round(SoilEvaporationsum);
  Rcout<<"  Herbaceous transpiration (mm) " << round(HerbTranspirationsum);
  Rcout<<" Woody plant transpiration (mm) "  <<round(Transpirationsum) <<"\n";
  Rcout<<"  Plant extraction from soil (mm) " << round(sum(PlantExtraction));
  if(plantResults) Rcout<<"  Plant water balance (mm) " << round(sum(PlantWaterBalance));
  NumericVector HydraulicRedistribution = DWB["HydraulicRedistribution"];
  Rcout<<" Hydraulic redistribution (mm) " << round(sum(HydraulicRedistribution)) <<"\n";
  Rcout<<"  Runoff (mm) " << round(Runoffsum) << " Deep drainage (mm) "  << round(DeepDrainagesum)  <<"\n";
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
//' a column called \code{"dates"} or \code{"Dates"} can contain \code{\link{Date}} or \code{\link{POSIXct}} classes.
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
//' @param waterTableDepth Water table depth (in mm). When not missing, capillarity rise will be allowed if lower than total soil depth.
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
//' Simulations using the 'Sperry' or 'Sureau' transpiration mode are computationally much more expensive than 'Granier'.
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
//'     \item{\code{"Rain"}: Precipitation as rainfall (in mm).}
//'     \item{\code{"Snow"}: Precipitation as snowfall (in mm).}
//'     \item{\code{"NetRain"}: Net rain, after accounting for interception (in mm).}
//'     \item{\code{"Infiltration"}: The amount of water infiltrating into the soil (in mm).}
//'     \item{\code{"InfiltrationExcess"}: Excess infiltration in the topmost layer leading to an increase in runoff (in mm).}
//'     \item{\code{"SaturationExcess"}: Excess saturation in the topmost layer leading to an increase in runoff (in mm).}
//'     \item{\code{"CapillarityRise"}: Water entering the soil via capillarity rise (mm) from the water table, if \code{waterTableDepth} is supplied.}
//'     \item{\code{"Runoff"}: The amount of water exported via surface runoff (in mm).}
//'     \item{\code{"DeepDrainage"}: The amount of water exported via deep drainage (in mm).}
//'     \item{\code{"Evapotranspiration"}: Evapotranspiration (in mm).}
//'     \item{\code{"SoilEvaporation"}: Bare soil evaporation (in mm).}
//'     \item{\code{"HerbTranspiration"}: Transpiration due to the herbaceous layer (in mm).}
//'     \item{\code{"PlantExtraction"}: Amount of water extracted from soil by woody plants (in mm).}
//'     \item{\code{"Transpiration"}: Woody plant transpiration (in mm).}
//'     \item{\code{"HydraulicRedistribution"}: Water redistributed among soil layers, transported through the plant hydraulic network.}
//'   }
//'   \item{\code{"EnergyBalance"}: A data frame with the daily values of energy balance components for the soil and the canopy (only for \code{transpirationMode = "Sperry"} or \code{transpirationMode = "Sureau"}).}
//'   \item{\code{"Temperature"}: A data frame with the daily values of minimum/mean/maximum temperatures for the atmosphere (input), canopy and soil (only for \code{transpirationMode = "Sperry"} or \code{transpirationMode = "Sureau"}).}
//'   \item{\code{"Soil"}: A list with the following subelements:}
//'   \itemize{
//'     \item{\code{"SWC"}: Soil water content (percent of soil volume) in each soil layer (and overall).}
//'     \item{\code{"RWC"}: Relative soil moisture content (relative to field capacity) in each soil layer (and overall).}
//'     \item{\code{"REW"}: Relative extractable water (min. psi = -5 MPa) in each soil layer (and overall).}
//'     \item{\code{"ML"}: Soil water volume in each soil layer (in L/m2) (and overall).}
//'     \item{\code{"Psi"}: Soil water potential in each soil layer (in MPa) (and overall).}
//'     \item{\code{"PlantExt"}: Plant extraction from each soil layer (in mm) (and overall).}
//'     \item{\code{"HydraulicInput"}: Water that entered the layer coming from other layers and transported via the plant hydraulic network (in mm) (and overall).}
//'   }
//'   \item{\code{"Snow"}: A data frame where the following variable (in columns) is given for each simulated day (in rows):}
//'   \itemize{
//'     \item{\code{"SWE"}: Snow water equivalent (mm) of the snow pack.}
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
//'     \item{\code{"LeafPLC"}: A data frame with the average daily proportion of leaf conductance loss of each plant (\[0-1\]).}
//'     \item{\code{"StemPLC"}: A data frame with the average daily proportion of stem conductance loss of each plant (\[0-1\]).}
//'     \item{\code{"PlantWaterBalance"}: A data frame with the daily balance between transpiration and soil water extraction for each plant cohort. }
//'     \item{\code{"LeafRWC"}: A data frame with the average daily leaf relative water content of each plant (in percent).}
//'     \item{\code{"StemRWC"}: A data frame with the average daily stem relative water content of each plant (in percent). }
//'     \item{\code{"LFMC"}: A data frame with the daily live fuel moisture content (in percent of dry weight).}
//'     \item{\code{"PlantStress"}: A data frame with the amount of daily stress \[0-1\] suffered by each plant cohort (relative whole-plant conductance).}
//'   }
//' If \code{transpirationMode="Sperry"} or \code{transpirationMode="Sureau"}, element \code{"Plants"} is a list with the following subelements:
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
//'     \item{\code{"LeafPLC"}: A data frame with the average daily proportion of leaf conductance loss of each plant (\[0-1\]).}
//'     \item{\code{"StemPLC"}: A data frame with the average daily proportion of stem conductance loss of each plant (\[0-1\]).}
//'     \item{\code{"RootPsi"}: A data frame with the minimum daily root water potential of each plant (in MPa). }
//'     \item{\code{"RhizoPsi"}: A list of data frames (one per plant cohort) with the minimum daily root water potential of each plant (in MPa).}
//'     \item{\code{"PlantStress"}: A data frame with the amount of daily stress \[0-1\] suffered by each plant cohort (relative whole-plant conductance).}
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
//' \code{\link{extract}}, \code{\link{summary.spwb}},  \code{\link{forest}}, \code{\link{aspwb}}
//' 
//' @examples
//' \donttest{
//' #Load example daily meteorological data
//' data(examplemeteo)
//' 
//' #Load example plot plant data
//' data(exampleforest)
//' 
//' #Default species parameterization
//' data(SpParamsMED)
//' 
//' #Define soil with default soil params (4 layers)
//' examplesoil <- defaultSoilParams(4)
//' 
//' #Initialize control parameters
//' control <- defaultControl("Granier")
//' 
//' #Initialize input
//' x1 <- spwbInput(exampleforest,examplesoil, SpParamsMED, control)
//' 
//' #Call simulation function
//' S1 <- spwb(x1, examplemeteo, latitude = 41.82592, elevation = 100)
//' 
//' #Switch to 'Sperry' transpiration mode
//' control <- defaultControl("Sperry")
//' 
//' #Initialize input
//' x2 <- spwbInput(exampleforest,examplesoil, SpParamsMED, control)
//' 
//' #Call simulation function
//' S2 <- spwb(x2, examplemeteo, latitude = 41.82592, elevation = 100)
//' 
//' #Switch to 'Sureau' transpiration mode
//' control <- defaultControl("Sureau")
//' 
//' #Initialize input
//' x3 <- spwbInput(exampleforest,examplesoil, SpParamsMED, control)
//' 
//' #Call simulation function
//' S3 <- spwb(x3, examplemeteo, latitude = 41.82592, elevation = 100)
//' }
//'                 
//' @name spwb
// [[Rcpp::export("spwb")]]
List spwb(List x, DataFrame meteo, 
          double latitude, double elevation, double slope = NA_REAL, double aspect = NA_REAL,
          NumericVector CO2ByYear = NumericVector(0), double waterTableDepth = NA_REAL) {
  
  //Clone input
  x = clone(x);
  
  List control = x["control"];
  String transpirationMode = control["transpirationMode"];
  String soilFunctions = control["soilFunctions"];
  String stemCavitationRecovery = control["stemCavitationRecovery"];
  String leafCavitationRecovery = control["leafCavitationRecovery"];
  bool verbose = control["verbose"];
  NumericVector defaultRainfallIntensityPerMonth = control["defaultRainfallIntensityPerMonth"];
  
  bool leafPhenology = control["leafPhenology"];
  bool unlimitedSoilWater = control["unlimitedSoilWater"];
  checkspwbInput(x,transpirationMode, soilFunctions);


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
  NumericVector RainfallIntensity(Precipitation.length(), NA_REAL);
  if(meteo.containsElementNamed("RainfallIntensity")) {
    RainfallIntensity = meteo["RainfallIntensity"];
    if(verbose) {
      Rcout<<"Rainfall intensity taken from input column 'RainfallIntensity'\n";
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
  
  //Soil
  DataFrame soil = Rcpp::as<Rcpp::DataFrame>(x["soil"]);
  
  //Define output list
  List outputList = defineSPWBDailyOutput(latitude, elevation, slope, aspect,
                                          dateStrings, x);
  outputList["weather"] = clone(meteo);

  //Initial soil status
  NumericVector initialSoilContent = water(soil, soilFunctions);
  NumericVector initialPlantContent = plantWaterContent(x);
  double initialSnowContent = x["snowpack"];
  if(verbose) {
    Rcout<<"Initial plant water content (mm): "<< sum(initialPlantContent)<<"\n";
    Rcout<<"Initial soil water content (mm): "<< sum(initialSoilContent)<<"\n";
    Rcout<<"Initial snowpack content (mm): "<< initialSnowContent<<"\n";
  }
  
  //Instance communication structures
  List internalCommunication = instanceCommunicationStructures(x, "spwb");

  
  bool error_occurence = false;
  if(verbose) Rcout << "Performing daily simulations\n";

  std::string yearString;
  for(int i=0;(i<numDays) && (!error_occurence);i++) {
     std::string c = as<std::string>(dateStrings[i]);
     yearString = c.substr(0, 4);
     if(verbose) {
        if(DOY[i]==1 || i==0) {
          Rcout<<"\n [Year "<< yearString << "]:";
        } 
        else if(i%30 == 0) Rcout<<".";//<<i;
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
      
      double Rint = RainfallIntensity[i];
      if(NumericVector::is_na(Rint)) {
        int month = std::atoi(c.substr(5,2).c_str());
        Rint = rainfallIntensity(month, Precipitation[i], defaultRainfallIntensityPerMonth);
      }
      
      //If DOY == 1 reset PLC (Growth assumed)
      if(stemCavitationRecovery=="annual") {
        if(DOY[i]==1) {
          DataFrame internalWater = Rcpp::as<Rcpp::DataFrame>(x["internalWater"]);
          NumericVector StemPLC = Rcpp::as<Rcpp::NumericVector>(internalWater["StemPLC"]);
          for(int j=0;j<StemPLC.length();j++) StemPLC[j] = 0.0;
          if(transpirationMode =="Sureau") {
            NumericVector StemPsi = Rcpp::as<Rcpp::NumericVector>(internalWater["StemPsi"]);
            NumericVector StemSympPsi = Rcpp::as<Rcpp::NumericVector>(internalWater["StemSympPsi"]);
            for(int j=0;j<StemPsi.length();j++) {
              StemPsi[j] = -0.033;
              StemSympPsi[j] = -0.033;
            } 
          }
        }
      }
      if(leafCavitationRecovery=="annual") {
        if(DOY[i]==1) {
          DataFrame internalWater = Rcpp::as<Rcpp::DataFrame>(x["internalWater"]);
          NumericVector LeafPLC = Rcpp::as<Rcpp::NumericVector>(internalWater["LeafPLC"]);
          for(int j=0;j<LeafPLC.length();j++) LeafPLC[j] = 0.0;
          if(transpirationMode =="Sureau") {
            NumericVector LeafPsi = Rcpp::as<Rcpp::NumericVector>(internalWater["LeafPsi"]);
            NumericVector LeafSympPsi = Rcpp::as<Rcpp::NumericVector>(internalWater["LeafSympPsi"]);
            for(int j=0;j<LeafPsi.length();j++) {
              LeafPsi[j] = -0.033;
              LeafSympPsi[j] = -0.033;
            }
          }
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
          Named("rint") = Rint);
        try{
          spwbDay_basic(internalCommunication, x, meteovec, 
                        elevation, slope, aspect, 
                        0.0, R_NilValue, waterTableDepth, 
                        verbose);
          fillSPWBDailyOutput(outputList, x, internalCommunication["basicSPWBOutput"],i);
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
          Named("rint") = Rint);
        try{
          spwbDay_advanced(internalCommunication, x, meteovec, 
                               latitude, elevation, slope, aspect,
                               solarConstant, delta, 
                               0.0, R_NilValue, waterTableDepth, 
                               verbose); 
          fillSPWBDailyOutput(outputList, x, internalCommunication["advancedSPWBOutput"],i);
        } catch(std::exception& ex) {
          Rcerr<< "c++ error: "<< ex.what() <<"\n";
          error_occurence = true;
        }
      }
  }
  if(verbose) Rcout << "\n\n";
  
  if(verbose) {
    printWaterBalanceResult(outputList, x,
                            initialPlantContent, initialSoilContent, initialSnowContent,
                            transpirationMode);
    if(error_occurence) {
      Rcout<< " ERROR: Calculations stopped because of numerical error: Revise parameters\n";
    }
  }
  
  return(outputList);
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
         double latitude, double elevation, double slope = NA_REAL, double aspect = NA_REAL, 
         NumericVector canopyEvaporation = NumericVector(0), 
         NumericVector snowMelt = NumericVector(0), 
         NumericVector soilEvaporation = NumericVector(0),
         NumericVector herbTranspiration = NumericVector(0),
         NumericVector CO2ByYear = NumericVector(0)) {
  
  //Clone input
  x = clone(x);

  
  List control = x["control"];
  String transpirationMode = control["transpirationMode"];
  String soilFunctions = control["soilFunctions"];
  String stemCavitationRecovery = control["stemCavitationRecovery"];
  String leafCavitationRecovery = control["leafCavitationRecovery"];
  bool verbose = control["verbose"];
  bool subdailyResults = control["subdailyResults"];
  bool leafPhenology = control["leafPhenology"];
  bool multiLayerBalance = control["multiLayerBalance"];
  
  //Store input
  List spwbInput = x; // Store initial object
  x = clone(x); //Ensure a copy will be modified
  
  DataFrame canopyParams = Rcpp::as<Rcpp::DataFrame>(x["canopy"]);
  int ncanlayers = canopyParams.nrow(); //Number of canopy layers
  DataFrame soil = Rcpp::as<Rcpp::DataFrame>(x["soil"]);
  
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
  if(transpirationMode!="Granier") DLT =  defineTemperatureLayersDailyOutput(dateStrings, canopy);
  
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

  //Fire hazard output variables
  DataFrame fireHazard;
  if(control["fireHazardResults"]) fireHazard = defineFireHazardOutput(dateStrings);
    
    
  List internalCommunication = instanceCommunicationStructures(x, "spwb");
  List SEBcommunication = internalCommunication["SEBcommunication"];
  List transpOutput;
  if(transpirationMode == "Granier") {
    transpOutput  = internalCommunication["basicTranspirationOutput"];
  } else {
    transpOutput  = internalCommunication["advancedTranspirationOutput"];
  }
  bool error_occurence = false;
  if(verbose) Rcout << "Performing daily simulations ";
  List s;
  std::string yearString;
  for(int i=0;i<numDays;i++) {
    std::string c = as<std::string>(dateStrings[i]);
    yearString = c.substr(0, 4);
    if(verbose) {
      if(DOY[i]==1 || i==0) {
        Rcout<<"\n Year "<< yearString<< ":";
      } 
      else if(i%30 == 0) Rcout<<".";//<<i;
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
    if(stemCavitationRecovery=="annual") {
        if(DOY[i]==1) {
          DataFrame internalWater = Rcpp::as<Rcpp::DataFrame>(x["internalWater"]);
          NumericVector StemPLC = Rcpp::as<Rcpp::NumericVector>(internalWater["StemPLC"]);
          for(int j=0;j<StemPLC.length();j++) StemPLC[j] = 0.0;
          if(transpirationMode =="Sureau") {
            NumericVector StemPsi = Rcpp::as<Rcpp::NumericVector>(internalWater["StemPsi"]);
            NumericVector StemSympPsi = Rcpp::as<Rcpp::NumericVector>(internalWater["StemSympPsi"]);
            for(int j=0;j<StemPsi.length();j++) {
              StemPsi[j] = -0.033;
              StemSympPsi[j] = -0.033;
            } 
          }
        }
    }
    if(leafCavitationRecovery=="annual") {
      if(DOY[i]==1) {
        DataFrame internalWater = Rcpp::as<Rcpp::DataFrame>(x["internalWater"]);
        NumericVector LeafPLC = Rcpp::as<Rcpp::NumericVector>(internalWater["LeafPLC"]);
        for(int j=0;j<LeafPLC.length();j++) LeafPLC[j] = 0.0;
        if(transpirationMode =="Sureau") {
          NumericVector LeafPsi = Rcpp::as<Rcpp::NumericVector>(internalWater["LeafPsi"]);
          NumericVector LeafSympPsi = Rcpp::as<Rcpp::NumericVector>(internalWater["LeafSympPsi"]);
          for(int j=0;j<LeafPsi.length();j++) {
            LeafPsi[j] = -0.033;
            LeafSympPsi[j] = -0.033;
          }
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
        transpirationBasic(transpOutput, x, meteovec, elevation, true);
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
        transpirationAdvanced(SEBcommunication, transpOutput, x, meteovec, 
                                latitude, elevation, slope, aspect,
                                solarConstant, delta,
                                canopyEvaporation[i], snowMelt[i], soilEvaporation[i], herbTranspiration[i],
                                verbose, NA_INTEGER, 
                                false);
      } catch(std::exception& ex) {
        Rcerr<< "c++ error: "<< ex.what() <<"\n";
        error_occurence = true;
      }
      fillEnergyBalanceDailyOutput(DEB,transpOutput, i, ntimesteps);
      if(control["temperatureResults"]) {
        fillTemperatureDailyOutput(DT,transpOutput, i, ntimesteps);
        if(multiLayerBalance) {
          fillTemperatureLayersDailyOutput(DLT,transpOutput, i, ncanlayers, ntimesteps);
        }
      }
    }
    
    //Update plant daily water output
    fillPlantWaterDailyOutput(plantDWOL, transpOutput, i, transpirationMode);
    if(transpirationMode!="Granier") fillSunlitShadeLeavesDailyOutput(sunlitDO, shadeDO, transpOutput, i, numCohorts);
    if(control["fireHazardResults"]) fillFireHazardOutput(fireHazard, transpOutput, i);
    
    List Plants = transpOutput["Plants"];
    NumericVector EplantCoh = Plants["Transpiration"];
    NumericMatrix SoilWaterExtract = transpOutput["Extraction"];
    for(int l=0;l<nlayers;l++) {
      Eplantdays(i,l) = sum(SoilWaterExtract(_,l));
    }

    PlantExtraction[i] = sum(SoilWaterExtract);
    Transpiration[i] = sum(EplantCoh);
    NumericVector HydrInVec(nlayers, 0.0);

    if(transpirationMode=="Sperry")  {
      NumericMatrix soilLayerExtractInst = transpOutput["ExtractionInst"];
      for(int l=0;l<nlayers;l++) {
        for(int n=0;n<ntimesteps;n++) {
          HydrInVec[l] += (-1.0)*std::min(soilLayerExtractInst(l,n),0.0);
        }
      }
      HydraulicRedistribution[i] = sum(HydrInVec);
      HydrIndays(i,_) = HydrInVec;
    } 
    NumericVector stand = transpOutput["Stand"];
    LAI[i] = stand["LAI"];
    LAIlive[i] = stand["LAIlive"];
    LAIexpanded[i] = stand["LAIexpanded"];
    LAIdead[i] = stand["LAIdead"];
        

    if(subdailyResults) {
      subdailyRes[i] = copyAdvancedTranspirationOutput(transpOutput, x);
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
  DataFrame Stand = DataFrame::create(_["LAI"]=LAI,
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
