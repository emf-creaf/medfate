#include <Rcpp.h>
#include <numeric>
#include <math.h>
#include <meteoland.h>

using namespace Rcpp;

const double Cp_Jmol = 29.37152; // J * mol^-1 * ºC^-1
const double SIGMA_W = 5.67*pow(10,-8.0); //Stefan-Boltzmann constant W * K^-4 * m^-2


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
 * and T = (tmin+tmax)/2 at sunset.
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
  double temp = 0.5*(tmin+tmax-(tmax-tmin)*cos(1.5*PI*t/daylength));
  return(temp);
}




/**
 * Calculates leaf temperature
 *   Campbell & Norman 1998 (eqns. 14.1 & 14.3)
 * 
 *  airTemperature - Air temperature (in ºC)
 *  absRad - Absorbed long- and short-wave radiation (in W*m^-2)
 *  E - Transpiration flow (in mmol H20 * m^-2 * s^-1) one sided leaf area basis
 *  leafWidth - Leaf width (in m)
 *  u - wind speed above the leaf boundary layer (in m/s)
 */
// [[Rcpp::export("biophysics.leafTemperature")]]
double leafTemperature(double absRad, double airTemperature, double u, double E,  double leafWidth = 0.01) {
  double lambda = meteoland::utils_latentHeatVaporisationMol(airTemperature);
  u = std::max(u, 0.1);//Force minimum wind speed to avoid excessive heating
  double gHa = 0.189*pow(u/(leafWidth*0.72), 0.5);
  double gr = 4.0*0.97*SIGMA_W*pow(273.16+airTemperature,3.0)/Cp_Jmol;
  double deltaTemp = (absRad- (0.97*SIGMA_W*pow(273.16+airTemperature,4.0)) - (lambda*(E/2000.0)))/(Cp_Jmol*(gr+gHa));
  return(airTemperature+deltaTemp);
}


/**
 * Calculates leaf relative water content from leaf water potential
 * 
 *  Bartlett, M. K., C. Scoffoni, and L. Sack. 2012. 
 *  The determinants of leaf turgor loss point and prediction of drought tolerance of species and biomes: a global meta-analysis. 
 *  Ecology letters 15:393–405.
 *  
 *  psi - Leaf water potential (MPa)
 *  pi0 - Full turgor osmotic potential (MPa)
 *  epsilon - bulk modulus elasticity (MPa)
 *  af - Apoplastic fraction (percentage)
 *  
 *  Returns Leaf RWC as percentage of maximum hydration (including apoplastic fraction)
 */
double leafRelativeWaterContent(double psi, double pi0, double epsilon, double af) {
  double psi_tl = (pi0*epsilon)/(pi0+epsilon);
  double rwc = 0;
  if(psi< psi_tl) {
    rwc = (-std::abs(pi0))/psi;
  } else {
    double c = std::abs(pi0);
    double b = psi+epsilon - c;
    double a = -epsilon;
    rwc = ((-b)-sqrt(pow(b,2.0)-4.0*a*c))/(2.0*a);
  }
  return((100.0-af)*rwc+af);
}


/**
 * Converts irradiance units (W*m-2) to quantum flux (micromol * m-2 * s-1), 
 * defined as the number of photons (in micromol) per second and unit area
 * 
 *  I - Irradiance (in W*m-2)
 *  lambda - wavelength (in nm)
 */
double irradianceToPhotonFlux(double I, double lambda = 546.6507) {
  return(I*lambda*0.836*pow(10.0,-2.0));
}

