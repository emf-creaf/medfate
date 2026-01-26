#define STRICT_R_HEADERS
#include <RcppArmadillo.h>
#include <numeric>
#include <math.h>
#include <meteoland.h>
#include "biophysicsutils.h"
#include "biophysicsutils_c.h"

using namespace Rcpp;

/**
 * Transforms dates (yyyy-mm-dd) into day of the year (DOY)
 */
IntegerVector date2doy(CharacterVector dateStrings) {
  IntegerVector doy(dateStrings.size());
  //Derive doy from date  
  for(int i=0;i<dateStrings.size();i++) {
    std::string c = as<std::string>(dateStrings[i]);
    int J = meteoland::radiation_julianDay(std::atoi(c.substr(0, 4).c_str()),std::atoi(c.substr(5,2).c_str()),std::atoi(c.substr(8,2).c_str()));
    int J0101 = meteoland::radiation_julianDay(std::atoi(c.substr(0, 4).c_str()),1,1);
    doy[i] = J - J0101+1;
  }
  return(doy);
}

NumericVector date2photoperiod(CharacterVector dateStrings, double latitude) {
  NumericVector photoperiod(dateStrings.size());
  //Derive photoperiod from date and latitude
  for(int i=0;i<dateStrings.size();i++) {
    std::string c = as<std::string>(dateStrings[i]);
    int J = meteoland::radiation_julianDay(std::atoi(c.substr(0, 4).c_str()),std::atoi(c.substr(5,2).c_str()),std::atoi(c.substr(8,2).c_str()));
    double delta = meteoland::radiation_solarDeclination(J);
    photoperiod[i] = meteoland::radiation_daylength(latitude, 0.0, 0.0, delta);
  }
  return(photoperiod);
}


