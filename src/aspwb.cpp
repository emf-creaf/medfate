// [[Rcpp::interfaces(r,cpp)]]
#define STRICT_R_HEADERS
#include <meteoland.h>
#include "soil.h"
#include "root.h"
#include "spwb.h"
#include "biophysicsutils.h"
#include "hydrology.h"
using namespace Rcpp;
using namespace meteoland;

NumericVector agricultureWaterInputs(List x, 
                                     double prec, double tday, double rad, double elevation,
                                     double LgroundSWR, 
                                     bool modifyInput = true) {

  double swe = x["snowpack"];
  
  //Snow pack dynamics
  double snow = 0.0, rain=0.0;
  double melt = 0.0;
  //Turn rain into snow and add it into the snow pack
  if(tday < 0.0) { 
    snow = prec; 
    swe = swe + snow;
  } else {
    rain = prec;
  }
  //Apply snow melting
  if(swe > 0.0) {
    melt = std::min(swe, snowMelt(tday, rad, LgroundSWR, elevation));
    // Rcout<<" swe: "<< swe<<" temp: "<<ten<< " rad: "<< ren << " melt : "<< melt<<"\n";
    swe = swe-melt;
  }
  
  //Hydrologic input
  double NetRain = 0.0, Interception = 0.0;
  if(rain>0.0)  {
    NetRain = rain - Interception; 
  }
  if(modifyInput) {
    x["snowpack"] = swe;
  }
  NumericVector WI = NumericVector::create(_["Rain"] = rain, _["Snow"] = snow,
                                           _["Interception"] = Interception,
                                           _["NetRain"] = NetRain, 
                                           _["Snowmelt"] = melt);
  return(WI);
}


//' @rdname aspwb
//' @keywords internal
// [[Rcpp::export("aspwbInput")]]
List aspwbInput(double crop_factor, List control, DataFrame soil) {
  
  String VG_PTF = control["VG_PTF"]; 
  DataFrame soil_out;
  if(soil.inherits("soil")) {
    soil_out = clone(soil);
  } else {
    soil_out = soilInit(soil, VG_PTF);
  }
  
  List input = List::create(_["control"] = clone(control),
                            _["crop_factor"] = crop_factor, 
                            _["snowpack"] = 0.0,
                            _["soil"] = soil_out);
  input.attr("class") = CharacterVector::create("aspwbInput","list");
  return(input);
}

List aspwb_day_internal(List x, NumericVector meteovec, 
              double elevation, double slope, double aspect, 
              double runon =  0.0, Nullable<NumericVector> lateralFlows = R_NilValue, double waterTableDepth = NA_REAL,
              bool verbose = false) {
  
  double crop_factor = x["crop_factor"];
  List control = x["control"];
  DataFrame soil = Rcpp::as<Rcpp::DataFrame>(x["soil"]);
  String soilFunctions = control["soilFunctions"];
  String infiltrationMode = control["infiltrationMode"];
  double infiltrationCorrection = control["infiltrationCorrection"];
  String soilDomains = control["soilDomains"];
  int ndailysteps = control["ndailysteps"];
  int max_nsubsteps_soil = control["max_nsubsteps_soil"];
  
  int nlayers = Rcpp::as<Rcpp::NumericVector>(soil["widths"]).size();
  
  //Weather input
  double tday = meteovec["tday"];
  double pet = meteovec["pet"]; 
  double prec  = meteovec["prec"];
  double rad = meteovec["rad"];
  double rainfallIntensity = meteovec["rint"];
  
  // Assume SWR is reduced with crop factor
  double LgroundSWR = 100.0 * (1.0 - crop_factor);
  
  //Snow pack dynamics and hydrology input (update snowpack)
  NumericVector hydroInputs = agricultureWaterInputs(x,
                                                     prec, tday, rad, elevation,
                                                     LgroundSWR,
                                                     true);
  double NetRain = hydroInputs["NetRain"];
  double Snowmelt = hydroInputs["Snowmelt"];  
  

  NumericVector widths = soil["widths"];
  NumericVector psiVec = psi(soil, soilFunctions); 
  double snowpack = x["snowpack"];
  //Evaporation from bare soil (if there is no snow), do not update soil yet
  double Esoil = soilEvaporation(soil, snowpack, 
                                 soilFunctions, pet, LgroundSWR, false);
  
  //Define plant net extraction 
  NumericVector ExtractionVec(nlayers, 0.0);
  // Transpiration is the product of PET and CROP FACTOR. HOWEVER, it is reduced with 
  double transp_max = pet*crop_factor; 
  //Calculate current soil water potential for transpiration
  NumericVector V = ldrRS_one(50, 500, widths);
  for(int l=0;l<nlayers;l++) {
    ExtractionVec[l] = V[l]*transp_max * exp(-0.6931472*pow(std::abs(psiVec[l]/(-2.0)),3.0)); //Reduce transpiration when soil is dry 
  }
  
  //Define source/sink with soil evaporation, herb transpiration and woody plant transpiration
  NumericVector sourceSinkVec(nlayers, 0.0);
  for(int l=0;l<nlayers;l++) {
    sourceSinkVec[l] -= ExtractionVec[l];
    if(l ==0) sourceSinkVec[l] -= Esoil;
  }
  
  //Determine water flows, returning deep drainage
  NumericVector sf = soilWaterBalance(soil, soilFunctions,
                                      NetRain, rainfallIntensity, Snowmelt, sourceSinkVec, 
                                      runon, lateralFlows, waterTableDepth,
                                      infiltrationMode, infiltrationCorrection, soilDomains, 
                                      ndailysteps, max_nsubsteps_soil, true);
  double Infiltration = sf["Infiltration"];
  double DeepDrainage = sf["DeepDrainage"];
  double Runoff = sf["Runoff"];
  double InfiltrationExcess = sf["InfiltrationExcess"];
  double SaturationExcess = sf["SaturationExcess"];
  double CapillarityRise = sf["CapillarityRise"];
  
  //Recalculate current soil water potential for output
  psiVec = psi(soil, soilFunctions); 
  
  NumericVector DB = NumericVector::create(_["PET"] = pet, 
                                           _["Rain"] = hydroInputs["Rain"], _["Snow"] = hydroInputs["Snow"], 
                                           _["NetRain"] = NetRain, _["Snowmelt"] = Snowmelt,
                                           _["Runon"] = runon, 
                                           _["Infiltration"] = Infiltration, _["InfiltrationExcess"] = InfiltrationExcess, _["SaturationExcess"] = SaturationExcess,
                                           _["Runoff"] = Runoff, _["DeepDrainage"] = DeepDrainage, _["CapillarityRise"] = CapillarityRise,
                                           _["SoilEvaporation"] = Esoil, _["Transpiration"] = sum(ExtractionVec));
  
  DataFrame SB = DataFrame::create(_["Psi"] = psiVec,
                                   _["PlantExtraction"] = ExtractionVec);
  List l = List::create(_["WaterBalance"] = DB, 
                        _["Soil"] = SB);
  l.attr("class") = CharacterVector::create("aspwb_day","list");
  return(l);
}

//' @rdname aspwb
//' @keywords internal
// [[Rcpp::export("aspwb_day")]]
List aspwb_day(List x, CharacterVector date, NumericVector meteovec, 
               double latitude, double elevation, double slope = NA_REAL, double aspect = NA_REAL,
               double runon =  0.0, Nullable<NumericVector> lateralFlows = R_NilValue, double waterTableDepth = NA_REAL,
               bool modifyInput = true) {
  
  double tmin = meteovec["MinTemperature"];
  double tmax = meteovec["MaxTemperature"];
  double rhmin = meteovec["MinRelativeHumidity"];
  double rhmax = meteovec["MaxRelativeHumidity"];
  double rad = meteovec["Radiation"];
  double prec = meteovec["Precipitation"];
  double wind = NA_REAL;
  if(meteovec.containsElementNamed("WindSpeed")) wind = meteovec["WindSpeed"];
  double Rint = NA_REAL; 
  if(meteovec.containsElementNamed("RainfallIntensity")) Rint = meteovec["RainfallIntensity"];
  
  List control = x["control"];
  bool verbose = control["verbose"];
  
  
  std::string c = as<std::string>(date[0]);
  int month = std::atoi(c.substr(5,2).c_str());
  int J = meteoland::radiation_julianDay(std::atoi(c.substr(0, 4).c_str()),std::atoi(c.substr(5,2).c_str()),std::atoi(c.substr(8,2).c_str()));
  double latrad = latitude * (M_PI/180.0);
  if(NumericVector::is_na(aspect)) aspect = 0.0;
  if(NumericVector::is_na(slope)) slope = 0.0;
  double asprad = aspect * (M_PI/180.0);
  double slorad = slope * (M_PI/180.0);
  double tday = meteoland::utils_averageDaylightTemperature(tmin, tmax);
  double pet = meteoland::penman(latrad, elevation, slorad, asprad, J, tmin, tmax, rhmin, rhmax, rad, wind);
  
  NumericVector defaultRainfallIntensityPerMonth = control["defaultRainfallIntensityPerMonth"];
  if(NumericVector::is_na(Rint)) Rint = rainfallIntensity(month, prec, defaultRainfallIntensityPerMonth);
  
  //Will not modify input x 
  if(!modifyInput) {
    x = clone(x);
  }
  
  
  NumericVector meteovec_inner = NumericVector::create(
    Named("tday") = tday, 
    Named("prec") = prec,
    Named("rad") = rad, 
    Named("pet") = pet,
    Named("rint") = Rint);
  
  return(aspwb_day_internal(x, meteovec_inner,
                            elevation, slope, aspect, 
                            runon, lateralFlows, waterTableDepth,
                            verbose));
}




DataFrame defineAgricultureWaterBalanceDailyOutput(CharacterVector dateStrings) {
  int numDays = dateStrings.length();
  
  NumericVector PET(numDays), Precipitation(numDays), Evapotranspiration(numDays);
  NumericVector Runoff(numDays),Rain(numDays),Snow(numDays);
  NumericVector Snowmelt(numDays),NetRain(numDays);
  NumericVector Infiltration(numDays),InfiltrationExcess(numDays),SaturationExcess(numDays),DeepDrainage(numDays), CapillarityRise(numDays);
  NumericVector SoilEvaporation(numDays),Transpiration(numDays);
  
  DataFrame DWB = DataFrame::create(_["PET"]=PET, 
                                    _["Precipitation"] = Precipitation, _["Rain"] = Rain, 
                                    _["Snow"] = Snow, 
                                    _["Snowmelt"] = Snowmelt, 
                                    _["Infiltration"]=Infiltration, 
                                    _["InfiltrationExcess"]=InfiltrationExcess, 
                                    _["SaturationExcess"]=SaturationExcess, 
                                    _["CapillarityRise"] = CapillarityRise,
                                    _["Runoff"]=Runoff, 
                                    _["DeepDrainage"]=DeepDrainage, 
                                    _["Evapotranspiration"]=Evapotranspiration,
                                    _["SoilEvaporation"]=SoilEvaporation,
                                    _["Transpiration"]=Transpiration);
  DWB.attr("row.names") = dateStrings;
  return(DWB);
}

// [[Rcpp::export(".defineASPWBDailyOutput")]]
List defineASPWBDailyOutput(double latitude, double elevation, double slope, double aspect, 
                           CharacterVector dateStrings, List x) {
  
  NumericVector topo = NumericVector::create(elevation, slope, aspect);
  topo.attr("names") = CharacterVector::create("elevation", "slope", "aspect");
  
  List aspwbInput = clone(x);
  List control = x["control"];
  DataFrame soil = Rcpp::as<Rcpp::DataFrame>(x["soil"]);
  
  //Water balance output variables
  DataFrame DWB = defineAgricultureWaterBalanceDailyOutput(dateStrings);
  List Soil = defineSoilDailyOutput(dateStrings, soil, false);
  DataFrame Snow = defineSnowDailyOutput(dateStrings);
  
  List l;
  l = List::create(Named("latitude") = latitude,
                   Named("topography") = topo,
                   Named("weather") = NA_REAL,
                   Named("aspwbInput") = aspwbInput,
                   Named("aspwbOutput") = x,
                   Named("WaterBalance")=DWB);
  if(control["soilResults"]) l.push_back(Soil, "Soil");
  if(control["snowResults"]) l.push_back(Snow, "Snow");
  l.attr("class") = CharacterVector::create("aspwb","list");
  return(l);
}

void fillAgricultureWaterBalanceDailyOutput(DataFrame DWB, List sDay, int iday) {
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
  NumericVector Transpiration = DWB["Transpiration"];
  NumericVector SoilEvaporation = DWB["SoilEvaporation"];
  NumericVector Evapotranspiration = DWB["Evapotranspiration"];
  DeepDrainage[iday] = db["DeepDrainage"];
  Infiltration[iday] = db["Infiltration"];
  InfiltrationExcess[iday] = db["InfiltrationExcess"];
  SaturationExcess[iday] = db["SaturationExcess"];
  CapillarityRise[iday] = db["CapillarityRise"];
  Runoff[iday] = db["Runoff"];
  Rain[iday] = db["Rain"];
  Snow[iday] = db["Snow"];
  PET[iday] = db["PET"];
  Precipitation[iday] = Rain[iday]+Snow[iday];
  Snowmelt[iday] = db["Snowmelt"];
  Transpiration[iday] = db["Transpiration"];
  SoilEvaporation[iday] = db["SoilEvaporation"];
  Evapotranspiration[iday] = Transpiration[iday]+SoilEvaporation[iday];
}

// [[Rcpp::export(".fillASPWBDailyOutput")]]
void fillASPWBDailyOutput(List l, List x, List sDay, int iday) {
  
  DataFrame soil = Rcpp::as<Rcpp::DataFrame>(x["soil"]);
  List control = x["control"];
  String transpirationMode = control["transpirationMode"];
  
  DataFrame DWB = Rcpp::as<Rcpp::DataFrame>(l["WaterBalance"]);
  int numDays = DWB.nrow();
  fillAgricultureWaterBalanceDailyOutput(DWB, sDay, iday);
  
  if(control["soilResults"]) {
    String soilFunctions = control["soilFunctions"];
    List Soil = Rcpp::as<Rcpp::List>(l["Soil"]);
    fillSoilDailyOutput(Soil, soil, sDay, 
                        iday, numDays, soilFunctions,
                        false);
  }
  if(control["snowResults"]) {
    DataFrame Snow = Rcpp::as<Rcpp::DataFrame>(l["Snow"]);
    fillSnowDailyOutput(Snow, x, iday);
  }
}

void printAgricultureWaterBalanceResult(DataFrame DWB,
                                        List x,
                                        NumericVector initialContent, double initialSnowContent) {
  
  List control = x["control"];
  DataFrame soil = Rcpp::as<Rcpp::DataFrame>(x["soil"]);
  String soilFunctions = control["soilFunctions"];
  NumericVector finalContent = water(soil, soilFunctions);
  double finalSnowContent = x["snowpack"];
  Rcout<<"Final soil water content (mm): "<< sum(finalContent)<<"\n";
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
  NumericVector Transpiration = DWB["Transpiration"];
  NumericVector SoilEvaporation = DWB["SoilEvaporation"];
  NumericVector Evapotranspiration = DWB["Evapotranspiration"];
  
  double Precipitationsum = sum(Precipitation);
  double Rainfallsum = sum(Rain);
  double SoilEvaporationsum = sum(SoilEvaporation);
  double Runoffsum  = sum(Runoff);
  double Infiltrationsum  = sum(Infiltration);
  double SaturationExcesssum  = sum(SaturationExcess);
  double InfiltrationExcesssum = sum(InfiltrationExcess);
  double CapillarityRisesum  = sum(CapillarityRise);
  double DeepDrainagesum = sum(DeepDrainage);
  double Transpirationsum = sum(Transpiration);
  double Snowmeltsum = sum(Snowmelt);
  double Snowsum = sum(Snow);
  
  double soil_wb = Infiltrationsum + CapillarityRisesum - SaturationExcesssum - DeepDrainagesum - SoilEvaporationsum - Transpirationsum;
  double snowpack_wb = Snowsum - Snowmeltsum;
  Rcout<<"Change in soil water content (mm): "<< sum(finalContent) - sum(initialContent)<<"\n";
  Rcout<<"Soil water balance result (mm): "<< soil_wb<<"\n";
  Rcout<<"Change in snowpack water content (mm): "<< finalSnowContent - initialSnowContent<<"\n";
  Rcout<<"Snowpack water balance result (mm): "<< snowpack_wb<<"\n";
  Rcout<<"Water balance components:\n";
  Rcout<<"  Precipitation (mm) "  <<round(Precipitationsum) <<"\n";
  Rcout<<"  Rain (mm) "  <<round(Rainfallsum) <<" Snow (mm) "  <<round(Snowsum) <<"\n";
  Rcout<<"  Infiltration (mm) " << round(Infiltrationsum)  << " Infiltration excess (mm) " << round(InfiltrationExcesssum) <<" Saturation excess (mm) " << round(SaturationExcesssum)<<" Capillarity rise (mm) " << round(CapillarityRisesum) << "\n";
  Rcout<<"  Soil evaporation (mm) " << round(SoilEvaporationsum) <<" Transpiration (mm) "  <<round(Transpirationsum) <<"\n";
  Rcout<<"  Runoff (mm) " << round(Runoffsum) << " Deep drainage (mm) "  << round(DeepDrainagesum)  <<"\n";
}

//' Simulation in agricultural areas
//'
//' Function \code{aspwb_day} performs water balance for a single day in an agriculture location.
//' Function \code{aspwb} performs water balance for multiple days in an agriculture location.
//' 
//' @param crop_factor Agriculture crop factor.
//' @param soil An object of class \code{\link{data.frame}} or \code{\link{soil}}.
//' @param control A list with default control parameters (see \code{\link{defaultControl}}).
//' @param x An object of class \code{\link{aspwbInput}}.
//' @param meteo A data frame with daily meteorological data series (see \code{\link{spwb}}). 
//' @param date Date as string "yyyy-mm-dd".
//' @param meteovec A named numerical vector with weather data. See variable names in parameter \code{meteo} of \code{\link{spwb}}.
//' @param latitude Latitude (in degrees).
//' @param elevation,slope,aspect Elevation above sea level (in m), slope (in degrees) and aspect (in degrees from North). 
//' @param runon Surface water amount running on the target area from upslope (in mm).
//' @param lateralFlows Lateral source/sink terms for each soil layer (interflow/to from adjacent locations) as mm/day.
//' @param waterTableDepth Water table depth (in mm). When not missing, capillarity rise will be allowed if lower than total soil depth.
//' @param modifyInput Boolean flag to indicate that the input \code{x} object is allowed to be modified during the simulation.
//' @param waterTableDepth Water table depth (in mm). When not missing, capillarity rise will be allowed if lower than total soil depth.
//'   
//' @author
//' Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
//' 
//' @seealso
//' \code{\link{spwbInput}}, \code{\link{spwb}}  
//' 
//' @examples
//' 
//' control <- defaultControl()
//' examplesoil <- defaultSoilParams(4)
//' 
//' x <- aspwbInput(0.75, control, examplesoil)
//' 
//' # Day to be simulated
//' d <- 100
//' meteovec <- unlist(examplemeteo[d,-1])
//' date <- as.character(examplemeteo$dates[d])
//' 
//' #Call simulation function for a single days
//' sd <- aspwb_day(x, date, meteovec,  
//'                latitude = 41.82592, elevation = 100) 
//' 
//' #Call simulation function for multiple days
//' S <- aspwb(x, examplemeteo, latitude = 41.82592, elevation = 100)
//' 
//' @name aspwb
//' @keywords internal
// [[Rcpp::export("aspwb")]]
List aspwb(List x, DataFrame meteo, double latitude, 
           double elevation, double slope = NA_REAL, double aspect = NA_REAL, 
           double waterTableDepth = NA_REAL) {
  List control = x["control"];
  String soilFunctions = control["soilFunctions"];
  bool verbose = control["verbose"];
  bool unlimitedSoilWater = control["unlimitedSoilWater"];
  NumericVector defaultRainfallIntensityPerMonth = control["defaultRainfallIntensityPerMonth"];
  
  //Clone input
  x = clone(x);
  
  DataFrame soil = Rcpp::as<Rcpp::DataFrame>(x["soil"]);
  

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
  if(any(is_na(MinRelativeHumidity))) stop("Missing values in 'MinRelativeHumidity'");
  if(any(is_na(MaxRelativeHumidity))) stop("Missing values in 'MaxRelativeHumidity'");
  if(any(is_na(Radiation))) stop("Missing values in 'Radiation'");
  
  NumericVector RainfallIntensity(Precipitation.length(), NA_REAL);
  if(meteo.containsElementNamed("RainfallIntensity")) {
    RainfallIntensity = meteo["RainfallIntensity"];
    if(verbose) {
      Rcout<<"Rainfall intensity taken from input column 'RainfallIntensity'\n";
    }
  }
  
  IntegerVector DOY, JulianDay;
  bool doy_input = false;
  if(meteo.containsElementNamed("DOY")) {
    DOY = meteo["DOY"];
    doy_input = true;
    if(verbose) {
      Rcout<<"DOY taken from input column 'DOY'\n";
    }
  }
  
  bool julianday_input = false;
  if(meteo.containsElementNamed("JulianDay")) {
    JulianDay = meteo["JulianDay"];
    julianday_input = true;
    if(verbose) {
      Rcout<<"Julian day taken from input column 'JulianDay'\n";
    }
  }
  
  CharacterVector dateStrings = getWeatherDates(meteo);
  if(!doy_input) DOY = date2doy(dateStrings);

  //Detailed subday results
  List subdailyRes(numDays);
  

  
  //Define output list
  List outputList = defineASPWBDailyOutput(latitude, elevation, slope, aspect,
                                          dateStrings, x);
  outputList["weather"] = clone(meteo);
  
  NumericVector initialContent = water(soil, soilFunctions);
  double initialSnowContent = x["snowpack"];
  if(verbose) {
    Rcout<<"Initial soil water content (mm): "<< sum(initialContent)<<"\n";
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
    
    
    double Rint = RainfallIntensity[i];
    if(NumericVector::is_na(Rint)) {
      int month = std::atoi(c.substr(5,2).c_str());
      Rint = rainfallIntensity(month, Precipitation[i], defaultRainfallIntensityPerMonth);
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

    double tmin = MinTemperature[i];
    double tmax = MaxTemperature[i];
    double tday = meteoland::utils_averageDaylightTemperature(tmin, tmax);
    double rhmin = MinRelativeHumidity[i];
    double rhmax = MaxRelativeHumidity[i];
    double rad = Radiation[i];
    
    PET[i] = meteoland::penman(latrad, elevation, slorad, asprad, J, 
                               tmin, tmax, rhmin, rhmax, rad, wind);
    

    //2. Water balance and photosynthesis
    NumericVector meteovec = NumericVector::create(
      Named("tday") = tday, 
      Named("prec") = Precipitation[i], 
      Named("rad") = rad, 
      Named("pet") = PET[i],
      Named("rint") = Rint);
    try{
      s = aspwb_day_internal(x, meteovec, 
                             elevation, slope, aspect,
                             0.0, R_NilValue, waterTableDepth, 
                             verbose);
    } catch(std::exception& ex) {
      Rcerr<< "c++ error: "<< ex.what() <<"\n";
      error_occurence = true;
    }
    
    //Fill output list      
    fillASPWBDailyOutput(outputList, x, s,i);

    if(control["subdailyResults"]) {
      subdailyRes[i] = clone(s);
    }
  }
  if(verbose) Rcout << "\n\n";
  
  if(verbose) {
    List DWB = outputList["WaterBalance"];
    printAgricultureWaterBalanceResult(DWB, x,
                                       initialContent, initialSnowContent);
    if(error_occurence) {
      Rcout<< " ERROR: Calculations stopped because of numerical error: Revise parameters\n";
    }
  }
  
  return(outputList);
}
