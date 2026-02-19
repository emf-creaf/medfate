#define STRICT_R_HEADERS
#include <RcppArmadillo.h>
#include <numeric>
#include "biophysicsutils_c.h"
#include "lightextinction_basic.h"
#include "lightextinction_advanced.h"
#include "windextinction_c.h"
#include "windKatul.h"
#include "hydraulics.h"
#include "hydrology.h"
#include "phenology.h"
#include "forestutils.h"
#include "tissuemoisture_c.h"
#include "carbon.h"
#include "carbon_c.h"
#include "spwb.h"
#include "root.h"
#include "soil.h"
#include "soil_c.h"
#include "soil_thermodynamics.h"
#include "inner_sperry.h"
#include "inner_sureau.h"
#include "inner_sureau_c.h"
#include "transpiration_advanced_c.h"
#include "meteoland/utils_c.hpp"
#include "meteoland/radiation_c.hpp"
#include "meteoland/pet_c.hpp"
#include <meteoland.h>
using namespace Rcpp;

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
  int J = julianDay_c(std::atoi(c.substr(0, 4).c_str()),std::atoi(c.substr(5,2).c_str()),std::atoi(c.substr(8,2).c_str()));
  double delta = solarDeclination_c(J);
  double solarConstant = solarConstant_c(J);
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
  
  double pet = PenmanPET_c(latrad, elevation, slorad, asprad, J, tmin, tmax, rhmin, rhmax, rad, wind);
  
  
  WeatherInputVector meteovec;
  meteovec.tmax = tmax;
  meteovec.tmin = tmin;
  meteovec.tminPrev = tminPrev;
  meteovec.tmaxPrev = tmaxPrev;
  meteovec.tminNext = tminNext;
  meteovec.rhmin = rhmin; 
  meteovec.rhmax = rhmax;
  meteovec.pet = pet;
  meteovec.wind = wind;
  meteovec.rad = rad;
  meteovec.prec = prec;
  meteovec.Catm = Catm;
  meteovec.Patm = Patm[day-1];
  
  ModelInput x_c = ModelInput(x);
  int nlayers = x_c.soil.getNlayers();
  int numCohorts = x_c.cohorts.SpeciesIndex.size();
  int ncanlayers = x_c.canopy.zlow.size();
  int ntimesteps = x_c.control.advancedWB.ndailysteps;
  
  
  AdvancedTranspiration_RESULT ATres(numCohorts, nlayers, ncanlayers, ntimesteps);
  AdvancedTranspiration_COMM ATcomm(numCohorts, nlayers, ncanlayers, ntimesteps);
  
  transpirationAdvanced_c(ATres, ATcomm, x_c, 
                          meteovec,
                          latitude, elevation, slope, aspect, 
                          solarConstant, delta,
                          canopyEvaporation, snowMelt, soilEvaporation, herbTranspiration, 
                          stepFunctions);
  
  List transpAdvanced = copyAdvancedTranspirationResult_c(ATres, x_c);
  
  if(modifyInput) {
    x_c.copyStateToList(x);
  }
  
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
  int J = julianDay_c(std::atoi(c.substr(0, 4).c_str()),std::atoi(c.substr(5,2).c_str()),std::atoi(c.substr(8,2).c_str()));
  double delta = solarDeclination_c(J);
  double solarConstant = solarConstant_c(J);
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
  
  double pet = PenmanPET_c(latrad, elevation, slorad, asprad, J, tmin, tmax, rhmin, rhmax, rad, wind);
  
  WeatherInputVector meteovec;
  meteovec.tmax = tmax;
  meteovec.tmin = tmin;
  meteovec.tminPrev = tminPrev;
  meteovec.tmaxPrev = tmaxPrev;
  meteovec.tminNext = tminNext;
  meteovec.rhmin = rhmin; 
  meteovec.rhmax = rhmax;
  meteovec.pet = pet;
  meteovec.wind = wind;
  meteovec.rad = rad;
  meteovec.prec = prec;
  meteovec.Catm = Catm;
  meteovec.Patm = Patm[day-1];
  
  ModelInput x_c = ModelInput(x);
  int nlayers = x_c.soil.getNlayers();
  int numCohorts = x_c.cohorts.SpeciesIndex.size();
  int ncanlayers = x_c.canopy.zlow.size();
  int ntimesteps = x_c.control.advancedWB.ndailysteps;
  

  AdvancedTranspiration_RESULT ATres(numCohorts, nlayers, ncanlayers, ntimesteps);
  AdvancedTranspiration_COMM ATcomm(numCohorts, nlayers, ncanlayers, ntimesteps);
  
  transpirationAdvanced_c(ATres, ATcomm, x_c, 
                          meteovec,
                          latitude, elevation, slope, aspect, 
                          solarConstant, delta,
                          canopyEvaporation, snowMelt, soilEvaporation, herbTranspiration, -1);
    
  List transpAdvanced = copyAdvancedTranspirationResult_c(ATres, x_c);

  if(modifyInput) {
    // Modify all state variables of input object from structure
    x_c.copyStateToList(x);
  }  
  return(transpAdvanced);
} 




