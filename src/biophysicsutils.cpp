#include <Rcpp.h>
#include <numeric>
#include <math.h>
#include <meteoland.h>

using namespace Rcpp;

const double Cp_Jmol = 29.37152; // J * mol^-1 * ºC^-1
const double SIGMA_W = 5.67*1e-8; //Stefan-Boltzmann constant W * K^-4 * m^-2

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


/**
 * Returns the proportion of daily radiation corresponding to the input time 
 * 
 * t - time of the day (in seconds from sunrise)
 * daylength - duration of the day (in seconds)
 * 
 * B. Y. H. Liu and R. C. Jordan, “The interrelationship and characteristic distribution of direct, diffuse and total solar radiation,” 
 * Solar Energy, vol. 4, no. 3, pp. 1–19, 1960. 
 */
// [[Rcpp::export("biophysics.radiationDiurnalPattern")]]
double radiationDiurnalPattern(double t, double daylength) {
  double ws = (daylength/3600.0)*(PI/24.0); //sunrise
  double w = ws - (t/daylength)*(ws*2.0);
  double prop = ((PI/24.0)*(cos(w)-cos(ws)))/(sin(ws)-ws*cos(ws));
  return(prop/3600.0);
}
/**
 * Calculated diurnal pattern of temperature assuming a sinusoidal pattern with T = tmin at sunrise
 * and T = (tmin+tmax)/2 at sunset. From sunset to sunrise follows a linear trend 
 * 
 * t - time of the day (in seconds from sunrise)
 * daylength - duration of the day (in seconds)
 *
 * McMurtrie, R. E., D. A. Rook, and F. M. Kelliher. 1990. 
 * Modelling the yield of Pinus radiata on a site limited by water and nitrogen. 
 * Forest Ecology and Management 30:381–413.
 */
// [[Rcpp::export("biophysics.temperatureDiurnalPattern")]]
double temperatureDiurnalPattern(double t, double tmin, double tmax, double daylength) {
  double temp;
  if((t<0) | (t>daylength)) {
    if(t<0) t = t + 86400.0;
    t = t - daylength;
    double tfin = 86400.0-daylength;
    temp = (0.5*(tmax-tmin)*(1.0-(t/tfin)) + tmin*(t/tfin));
  } else {
    temp = 0.5*(tmin+tmax-(tmax-tmin)*cos(1.5*PI*t/daylength));
  }
  return(temp);
}




/**
 * Calculates leaf temperature
 *   Campbell & Norman 1998 (eqns. 14.1 & 14.3)
 * 
 *  airTemperature - Air temperature (in ºC)
 *  absRad - Absorbed long- and short-wave radiation (in W*m^-2)
 *  E - Transpiration flow (in mmol H20 * m^-2 * s^-1) one sided leaf area basis
 *  leafWidth - Leaf width (here in cm)
 *  u - wind speed above the leaf boundary layer (in m/s)
 */
// [[Rcpp::export("biophysics.leafTemperature")]]
double leafTemperature(double absRad, double airTemperature, double u, double E,  double leafWidth = 1.0) {
  double lambda = meteoland::utils_latentHeatVaporisationMol(airTemperature);
  u = std::max(u, 0.1);//Force minimum wind speed to avoid excessive heating
  double gHa = 0.189*pow(u/(leafWidth*0.0072), 0.5);
  double gr = 4.0*0.97*SIGMA_W*pow(273.16+airTemperature,3.0)/Cp_Jmol;
  double deltaTemp = (absRad- (0.97*SIGMA_W*pow(273.16+airTemperature,4.0)) - (lambda*(E/2000.0)))/(Cp_Jmol*(gr+gHa));
  return(airTemperature+deltaTemp);
}


/**
 * Converts irradiance units (W*m-2) to quantum flux (micromol * m-2 * s-1), 
 * defined as the number of photons (in micromol) per second and unit area
 * 
 *  I - Irradiance (in W*m-2)
 *  lambda - wavelength (in nm)
 */
double irradianceToPhotonFlux(double I, double lambda = 546.6507) {
  return(I*lambda*0.836*1e-2);
}

