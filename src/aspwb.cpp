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

NumericVector agricultureSoilWaterInputs(List soil, 
                                         double prec, double tday, double rad, double elevation,
                                         double LgroundSWR, 
                                         double runon = 0.0,
                                         bool snowpack = true, bool modifySoil = true) {
  //Soil input
  double swe = soil["SWE"]; //snow pack
  
  //Snow pack dynamics
  double snow = 0.0, rain=0.0;
  double melt = 0.0;
  if(snowpack) {
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
  } else {
    rain = prec;
  }
  
  //Hydrologic input
  double NetRain = 0.0, Interception = 0.0;
  if(rain>0.0)  {
    NetRain = rain - Interception; 
  }
  if(modifySoil) {
    soil["SWE"] = swe;
  }
  NumericVector WI = NumericVector::create(_["Rain"] = rain, _["Snow"] = snow,
                                           _["Interception"] = Interception,
                                           _["NetRain"] = NetRain, 
                                           _["Snowmelt"] = melt,
                                           _["Runon"] = runon,
                                           _["RainfallInput"] = runon+NetRain,
                                           _["TotalInput"] = runon+melt+NetRain);
  return(WI);
}


//' @rdname aspwb
// [[Rcpp::export("aspwbInput")]]
List aspwbInput(double crop_factor, List control, List soil) {
  List input = List::create(_["crop_factor"] = crop_factor, 
                            _["control"] = clone(control),
                            _["soil"] = clone(soil));
  input.attr("class") = CharacterVector::create("aspwbInput","list");
  return(input);
}

List aspwb_day_internal(List x, NumericVector meteovec, 
              double elevation, double slope, double aspect, 
              double runon=0.0, bool verbose=false) {
  
  double crop_factor = x["crop_factor"];
  List control = x["control"];
  List soil = x["soil"];
  bool snowpack = control["snowpack"];
  String soilFunctions = control["soilFunctions"];
  
  int nlayers = Rcpp::as<Rcpp::NumericVector>(soil["dVec"]).size();
  
  //Weather input
  double tday = meteovec["tday"];
  double pet = meteovec["pet"]; 
  double prec  = meteovec["prec"];
  double rad = meteovec["rad"];

  // Assume SWR is reduced with crop factor
  double LgroundSWR = 100.0 * (1.0 - crop_factor);
  
  //Snow pack dynamics and hydrology input (update SWE)
  NumericVector hydroInputs = agricultureSoilWaterInputs(soil, prec,
                                                         tday, rad, elevation,
                                                         LgroundSWR, runon,
                                                         snowpack, true);

  //Soil infiltration and percolation (update W)
  NumericVector dVec = soil["dVec"];
  NumericVector macro = soil["macro"];
  NumericVector Water_FC = waterFC(soil, soilFunctions);
  NumericVector psiVec = psi(soil, soilFunctions); 
  
  double RainfallInput = hydroInputs["RainfallInput"];
  double Snowmelt = hydroInputs["Snowmelt"];
  double Infiltration = infiltrationAmount(RainfallInput, Water_FC[0]);
  double Runoff = RainfallInput - Infiltration;
  //Decide infiltration repartition among layers
  NumericVector Ivec = infiltrationRepartition(Infiltration, dVec, macro);
  
  //Evaporation from bare soil (if there is no snow), do not update soil yet
  double Esoil = soilEvaporation(soil, soilFunctions, pet, LgroundSWR, false);
  
  //Define plant net extraction 
  NumericVector ExtractionVec(nlayers, 0.0);
  // Transpiration is the product of PET and CROP FACTOR. HOWEVER, it is reduced with 
  double transp_max = pet*crop_factor; 
  //Calculate current soil water potential for transpiration
  NumericVector V = ldrRS_one(50, 500, dVec);
  for(int l=0;l<nlayers;l++) {
    ExtractionVec[l] = V[l]*transp_max * exp(-0.6931472*pow(std::abs(psiVec[l]/(-2.0)),3.0)); //Reduce transpiration when soil is dry 
  }
  
  //Modify soil with soil evaporation, herb transpiration and woody plant transpiration
  //and determine water flows, returning deep drainage
  NumericVector sourceSinkVec(nlayers, 0.0);
  for(int l=0;l<nlayers;l++) {
    sourceSinkVec[l] += Ivec[l] - ExtractionVec[l];
    if(l ==0) sourceSinkVec[l] += Snowmelt - Esoil;
  }
  double DeepDrainage = soilFlows(soil, sourceSinkVec, 24, true);
    
  //Recalculate current soil water potential for output
  psiVec = psi(soil, soilFunctions); 
  
  NumericVector DB = NumericVector::create(_["PET"] = pet, 
                                           _["Rain"] = hydroInputs["Rain"], _["Snow"] = hydroInputs["Snow"], 
                                           _["NetRain"] = hydroInputs["NetRain"], _["Snowmelt"] = Snowmelt,
                                           _["Runon"] = hydroInputs["Runon"], 
                                           _["Infiltration"] = Infiltration, _["Runoff"] = Runoff, 
                                           _["DeepDrainage"] = DeepDrainage,
                                           _["SoilEvaporation"] = Esoil, _["Transpiration"] = sum(ExtractionVec));
  
  DataFrame SB = DataFrame::create(_["PlantExtraction"] = ExtractionVec, 
                                   _["psi"] = psiVec);
  List l = List::create(_["WaterBalance"] = DB, 
                        _["Soil"] = SB);
  l.attr("class") = CharacterVector::create("aspwb_day","list");
  return(l);
}

//' @rdname aspwb
// [[Rcpp::export("aspwb_day")]]
List aspwb_day(List x, CharacterVector date, NumericVector meteovec, 
               double latitude, double elevation, double slope, double aspect,
               double runon=0.0, bool modifyInput = true) {
  
  double tmin = meteovec["MinTemperature"];
  double tmax = meteovec["MaxTemperature"];
  double rhmin = meteovec["MinRelativeHumidity"];
  double rhmax = meteovec["MaxRelativeHumidity"];
  double rad = meteovec["Radiation"];
  double prec = meteovec["Precipitation"];
  double wind = NA_REAL;
  if(meteovec.containsElementNamed("WindSpeed")) wind = meteovec["WindSpeed"];

  List control = x["control"];
  bool verbose = control["verbose"];
  
  
  std::string c = as<std::string>(date[0]);
  int J = meteoland::radiation_julianDay(std::atoi(c.substr(0, 4).c_str()),std::atoi(c.substr(5,2).c_str()),std::atoi(c.substr(8,2).c_str()));
  double latrad = latitude * (M_PI/180.0);
  if(NumericVector::is_na(aspect)) aspect = 0.0;
  if(NumericVector::is_na(slope)) slope = 0.0;
  double asprad = aspect * (M_PI/180.0);
  double slorad = slope * (M_PI/180.0);
  double tday = meteoland::utils_averageDaylightTemperature(tmin, tmax);
  double pet = meteoland::penman(latrad, elevation, slorad, asprad, J, tmin, tmax, rhmin, rhmax, rad, wind);
  
  
  //Will not modify input x 
  if(!modifyInput) {
    x = clone(x);
  }
  
  NumericVector meteovec_inner = NumericVector::create(
    Named("tday") = tday, 
    Named("prec") = prec,
    Named("rad") = rad, 
    Named("pet") = pet);
  
  return(aspwb_day_internal(x, meteovec_inner,
                  elevation, slope, aspect, 
                  runon, verbose));
}


DataFrame defineAgricultureWaterBalanceDailyOutput(CharacterVector dateStrings, NumericVector PET) {
  int numDays = dateStrings.length();
  
  NumericVector Precipitation(numDays), Evapotranspiration(numDays);
  NumericVector Runoff(numDays),Rain(numDays),Snow(numDays);
  NumericVector Snowmelt(numDays),NetRain(numDays);
  NumericVector Infiltration(numDays),DeepDrainage(numDays);
  NumericVector SoilEvaporation(numDays),Transpiration(numDays);
  
  DataFrame DWB = DataFrame::create(_["PET"]=PET, 
                                    _["Precipitation"] = Precipitation, _["Rain"] = Rain, 
                                    _["Snow"] = Snow, 
                                    _["Snowmelt"] = Snowmelt, 
                                    _["Infiltration"]=Infiltration, 
                                    _["Runoff"]=Runoff, 
                                    _["DeepDrainage"]=DeepDrainage, 
                                    _["Evapotranspiration"]=Evapotranspiration,
                                    _["SoilEvaporation"]=SoilEvaporation,
                                    _["Transpiration"]=Transpiration);
  DWB.attr("row.names") = dateStrings;
  return(DWB);
}

DataFrame defineAgricultureSoilWaterBalanceDailyOutput(CharacterVector dateStrings, List soil) {
  int numDays = dateStrings.length();
  NumericVector W = soil["W"];
  int nlayers = W.length();
  
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
                                    _["psi"]=psidays); 
  SWB.attr("row.names") = dateStrings;
  return(SWB);  
}

void fillAgricultureWaterBalanceDailyOutput(DataFrame DWB, List sDay, int iday) {
  List db = sDay["WaterBalance"];
  NumericVector Precipitation = DWB["Precipitation"];
  NumericVector DeepDrainage = DWB["DeepDrainage"];
  NumericVector Infiltration = DWB["Infiltration"];
  NumericVector Runoff = DWB["Runoff"];
  NumericVector Rain = DWB["Rain"];
  NumericVector Snow = DWB["Snow"];
  NumericVector Snowmelt = DWB["Snowmelt"];
  NumericVector Transpiration = DWB["Transpiration"];
  NumericVector SoilEvaporation = DWB["SoilEvaporation"];
  NumericVector Evapotranspiration = DWB["Evapotranspiration"];
  DeepDrainage[iday] = db["DeepDrainage"];
  Infiltration[iday] = db["Infiltration"];
  Runoff[iday] = db["Runoff"];
  Rain[iday] = db["Rain"];
  Snow[iday] = db["Snow"];
  Precipitation[iday] = Rain[iday]+Snow[iday];
  Snowmelt[iday] = db["Snowmelt"];
  Transpiration[iday] = db["Transpiration"];
  SoilEvaporation[iday] = db["SoilEvaporation"];
  Evapotranspiration[iday] = Transpiration[iday]+SoilEvaporation[iday];
}
void fillAgricultureSoilWaterBalanceDailyOutput(DataFrame SWB, List soil, List sDay, 
                                     int iday, int numDays, String soilFunctions) {
  NumericVector W = soil["W"];
  int nlayers = W.length();
  NumericVector Water_FC = waterFC(soil, soilFunctions);
  
  List sb = sDay["Soil"];
  NumericVector psi = sb["psi"];

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
    NumericVector Wdays = as<Rcpp::NumericVector>(SWB[wS]);
    NumericVector MLdays = as<Rcpp::NumericVector>(SWB[mlS]);
    NumericVector psidays = as<Rcpp::NumericVector>(SWB[psiS]);
    psidays[iday] = psi[l];
    if(iday<(numDays-1)) Wdays[iday+1] = W[l];
    MLdays[iday] = Wdays[iday]*Water_FC[l]; 
    MLTot[iday] = MLTot[iday] + MLdays[iday];
  }
  SWE[iday] = soil["SWE"];
  WaterTable[iday] = waterTableDepth(soil, soilFunctions);
}

void printAgricultureWaterBalanceResult(DataFrame DWB,
                             List soil, String soilFunctions,
                             NumericVector initialContent, double initialSnowContent) {
  
  NumericVector finalContent = water(soil, soilFunctions);
  double finalSnowContent = soil["SWE"];
  Rcout<<"Final soil water content (mm): "<< sum(finalContent)<<"\n";
  Rcout<<"Final snowpack content (mm): "<< finalSnowContent<<"\n";
  
  NumericVector Precipitation = DWB["Precipitation"];
  NumericVector DeepDrainage = DWB["DeepDrainage"];
  NumericVector Infiltration = DWB["Infiltration"];
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
  double DeepDrainagesum = sum(DeepDrainage);
  double Transpirationsum = sum(Transpiration);
  double Snowmeltsum = sum(Snowmelt);
  double Snowsum = sum(Snow);
  
  double soil_wb = Rainfallsum + Snowmeltsum - Runoffsum - DeepDrainagesum - SoilEvaporationsum - Transpirationsum;
  double snowpack_wb = Snowsum - Snowmeltsum;
  Rcout<<"Change in soil water content (mm): "<< sum(finalContent) - sum(initialContent)<<"\n";
  Rcout<<"Soil water balance result (mm): "<< soil_wb<<"\n";
  Rcout<<"Change in snowpack water content (mm): "<< finalSnowContent - initialSnowContent<<"\n";
  Rcout<<"Snowpack water balance result (mm): "<< snowpack_wb<<"\n";
  Rcout<<"Water balance components:\n";
  Rcout<<"  Precipitation (mm) "  <<round(Precipitationsum) <<"\n";
  Rcout<<"  Rain (mm) "  <<round(Rainfallsum) <<" Snow (mm) "  <<round(Snowsum) <<"\n";
  Rcout<<"  Infiltration (mm) " << round(Infiltrationsum)  <<
    " Runoff (mm) " << round(Runoffsum) <<
      " Deep drainage (mm) "  << round(DeepDrainagesum)  <<"\n";
  Rcout<<"  Soil evaporation (mm) " << round(SoilEvaporationsum);
  Rcout<<" Transpiration (mm) "  <<round(Transpirationsum) <<"\n";
  Rcout <<"\n";
}

//' Simulation in agricultural areas
//'
//' Function \code{aspwb_day} performs water balance for a single day in an agriculture location.
//' 
//' @param crop_factor Agriculture crop factor.
//' @param soil An object of class \code{\link{soil}}.
//' @param control A list with default control parameters (see \code{\link{defaultControl}}).
//' @param x An object of class \code{\link{aspwbInput}}.
//' @param meteo A data frame with daily meteorological data series (see \code{\link{spwb}}). 
//' @param date Date as string "yyyy-mm-dd".
//' @param meteovec A named numerical vector with weather data. See variable names in parameter \code{meteo} of \code{\link{spwb}}.
//' @param latitude Latitude (in degrees).
//' @param elevation,slope,aspect Elevation above sea level (in m), slope (in degrees) and aspect (in degrees from North). 
//' @param runon Surface water amount running on the target area from upslope (in mm).
//' @param modifyInput Boolean flag to indicate that the input \code{x} object is allowed to be modified during the simulation.
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
//' examplesoil <- soil(defaultSoilParams())
//' x <- aspwbInput(0.75, control, examplesoil)
//' 
//' #Call simulation function
//' S <- aspwb(x, examplemeteo, latitude = 41.82592, elevation = 100)
//' 
//' @name aspwb
// [[Rcpp::export("aspwb")]]
List apwb(List x, DataFrame meteo, double latitude, double elevation = NA_REAL, double slope = NA_REAL, double aspect = NA_REAL) {
  List control = x["control"];
  String soilFunctions = control["soilFunctions"];
  bool verbose = control["verbose"];
  bool unlimitedSoilWater = control["unlimitedSoilWater"];
  
  //Store input
  List aspwbInput = x; // Store initial object
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
  if(any(is_na(MinRelativeHumidity))) stop("Missing values in 'MinRelativeHumidity'");
  if(any(is_na(MaxRelativeHumidity))) stop("Missing values in 'MaxRelativeHumidity'");
  if(any(is_na(Radiation))) stop("Missing values in 'Radiation'");
  
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
  

  //Water balance output variables
  DataFrame DWB = defineAgricultureWaterBalanceDailyOutput(dateStrings, PET);
  DataFrame SWB = defineAgricultureSoilWaterBalanceDailyOutput(dateStrings, soil);
  
  NumericVector initialContent = water(soil, soilFunctions);
  double initialSnowContent = soil["SWE"];
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
      Named("pet") = PET[i]);
    try{
      s = aspwb_day_internal(x, meteovec, 
                             elevation, 
                             0.0, verbose); //No Runon in simulations for a single cell
    } catch(std::exception& ex) {
      Rcerr<< "c++ error: "<< ex.what() <<"\n";
      error_occurence = true;
    }
    
    //Update plant daily water output
    fillAgricultureWaterBalanceDailyOutput(DWB, s,i);
    fillAgricultureSoilWaterBalanceDailyOutput(SWB, soil, s,
                                    i, numDays, soilFunctions);

    if(control["subdailyResults"]) {
      subdailyRes[i] = clone(s);
    }
  }
  if(verbose) Rcout << "\n\n";
  
  if(verbose) {
    printAgricultureWaterBalanceResult(DWB, soil, soilFunctions,
                                       initialContent, initialSnowContent);
    if(error_occurence) {
      Rcout<< " ERROR: Calculations stopped because of numerical error: Revise parameters\n";
    }
  }
  
  
  subdailyRes.attr("names") = dateStrings;
  
  NumericVector topo = NumericVector::create(elevation, slope, aspect);
  topo.attr("names") = CharacterVector::create("elevation", "slope", "aspect");
  

  List l;
  l = List::create(Named("latitude") = latitude,
                   Named("topography") = topo,
                   Named("weather") = clone(meteo),
                   Named("aspwbInput") = aspwbInput,
                   Named("aspwbOutput") = clone(x),
                   Named("WaterBalance")=DWB);
  if(control["soilResults"]) l.push_back(SWB, "Soil");
  if(control["subdailyResults"]) l.push_back(subdailyRes,"subdaily");
  l.attr("class") = CharacterVector::create("aspwb","list");
  return(l);
}
