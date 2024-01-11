#define STRICT_R_HEADERS
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

//' Physical and biophysical utility functions
//' 
//' Utility functions for the calculation of biophysical variables. 
//'
//' @param t Time of the day (in seconds).
//' @param daylength Day length (in seconds).
//' 
//' @details 
//' Functions \code{biophysics_leafTemperature} and \code{biophysics_leafTemperature2} calculate leaf temperature according to energy balance equation given in Campbell and Norman (1988). 
//' 
//' Function \code{biophysics_radiationDiurnalPattern} follows the equations given in Liu and Jordan (1960). 
//' 
//' Function \code{biophysics_temperatureDiurnalPattern} determines diurnal temperature pattern assuming a sinusoidal pattern with T = Tmin at sunrise and T = (Tmin+Tmax)/2 at sunset and a linear change in temperature between sunset and Tmin of the day after (McMurtrie et al. 1990). 
//' 
//' Function \code{biophysics_waterDynamicViscosity} calculates water dynamic viscosity following the Vogel (1921) equation.
//' 
//' @return
//' Values returned for each function are:
//' \itemize{
//'   \item{\code{biophysics_leafTemperature} and \code{biophysics_leafTemperature2}: leaf temperature (in ºC)}
//'   \item{\code{biophysics_leafVapourPressure}: leaf vapour pressure (in kPa)} 
//'   \item{\code{biophysics_radiationDiurnalPattern}: the proportion of daily radiation corresponding to the input time in seconds after sunrise.} 
//'   \item{\code{biophysics_temperatureDiurnalPattern}: diurnal pattern of temperature.}
//'   \item{\code{biophysics_waterDynamicViscosity}: Water dynamic viscosity relative to 20ºC.} 
//' }
//' 
//' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
//' 
//' @references
//' Campbell, G. S., and J. M. Norman. 1998. An introduction to environmental biophysics: 2nd edition. (eqns. 14.1 & 14.3)
//'   
//' B. Y. H. Liu and R. C. Jordan, “The interrelationship and characteristic distribution of direct, diffuse and total solar radiation,” Solar Energy, vol. 4, no. 3, pp. 1–19, 1960. 
//' 
//' McMurtrie, R. E., D. A. Rook, and F. M. Kelliher. 1990. Modelling the yield of Pinus radiata on a site limited by water and nitrogen. Forest Ecology and Management 30:381–413.
//' 
//' H. Vogel, "Das Temperaturabhangigkeitsgesetz der Viskositat von Flussigkeiten", Physikalische Zeitschrift, vol. 22, pp. 645–646, 1921.
//' 
//' @seealso \code{\link{spwb}}
//' 
//' @name biophysics
// [[Rcpp::export("biophysics_radiationDiurnalPattern")]]
double radiationDiurnalPattern(double t, double daylength) {
  double ws = (daylength/3600.0)*(M_PI/24.0); //sunrise
  double w = ws - (t/daylength)*(ws*2.0);
  double prop = ((M_PI/24.0)*(cos(w)-cos(ws)))/(sin(ws)-ws*cos(ws));
  return(prop/3600.0);
}

/**
 * Calculated diurnal pattern of temperature assuming a sinusoidal pattern with T = tmin at sunrise
 * and T = (tmin+tmax)/2 at sunset. From sunset to sunrise follows a linear trend 
 * 
 * t - time of the day (in seconds from sunrise)
 * daylength - duration of the day (in seconds)
 *
 */
//' @rdname biophysics
//' @param tmin,tmax Minimum and maximum daily temperature (ºC).
//' @param tminPrev,tmaxPrev,tminNext Maximum and minimum daily temperatures of the previous and following day (ºC).
// [[Rcpp::export("biophysics_temperatureDiurnalPattern")]]
double temperatureDiurnalPattern(double t, double tmin, double tmax, 
                                 double tminPrev, double tmaxPrev, double tminNext, double daylength) {
  double temp;
  if((t<0.0) | (t>daylength)) {
    double tfin = 86400.0-daylength;
    if(t<0.0) {
      t = t + 86400.0 - daylength;
      temp = (0.5*(tmaxPrev+tminPrev)*(1.0-(t/tfin)) + tmin*(t/tfin));
      // Rcout<< t << " dl "<< daylength <<" ("<< tminPrev<< ", "<< tmaxPrev<<") ("<< tmin<< ", "<<tmax<<") to ("<<tminNext<<",) ";
      // Rcout << " from prev tfin "<< tfin<< " ratio "<< t/tfin;
    } else {
      t = t - daylength;
      temp = (0.5*(tmax+tmin)*(1.0-(t/tfin)) + tminNext*(t/tfin));
      // Rcout<< t << " dl "<< daylength <<" ("<< tminPrev<< ", "<< tmaxPrev<<") ("<< tmin<< ", "<<tmax<<") to ("<<tminNext<<",) ";
      // Rcout << " to next tfin "<< tfin<< " ratio "<< t/tfin;
    }
    // Rcout<<" "<< temp <<"\n";
  } else {
    double ct = cos(1.5*M_PI*t/daylength);
    temp = 0.5*(tmin+tmax-(tmax-tmin)*ct);
    // Rcout<<" t "<< t << " dl "<< daylength << " ct "<<ct <<" "<< 0.5*(tmax+tmin)<<" "<<temp <<"\n";
  }
  return(temp);
}

//' @rdname biophysics
//' 
//' @param u Wind speed above the leaf boundary layer (in m/s).
//' @param airTemperature Air temperature (in ºC).
//' @param absRad Absorbed long- and short-wave radiation (in W·m-2).
//' @param E Transpiration flow (in mmol H20·m-2·s-1) per one sided leaf area basis.
//' @param leafWidth Leaf width (in cm).
//' 
// [[Rcpp::export("biophysics_leafTemperature")]]
double leafTemperature(double absRad, double airTemperature, double u, double E,  double leafWidth = 1.0) {
  double lambda = meteoland::utils_latentHeatVaporisationMol(airTemperature);
  u = std::max(u, 0.1);//Force minimum wind speed to avoid excessive heating
  double gHa = 0.189*pow(u/(leafWidth*0.0072), 0.5);
  double gr = 4.0*0.97*SIGMA_W*pow(273.16+airTemperature,3.0)/Cp_Jmol;
  double deltaTemp = (absRad- (0.97*SIGMA_W*pow(273.16+airTemperature,4.0)) - (lambda*(E/2000.0)))/(Cp_Jmol*(gr+gHa));
  return(airTemperature+deltaTemp);
}

//' @rdname biophysics
//' @param SWRabs Absorbed short-wave radiation (in W·m-2).
//' @param LWRnet Net long-wave radiation balance (in W·m-2).
// [[Rcpp::export("biophysics_leafTemperature2")]]
double leafTemperature2(double SWRabs, double LWRnet, double airTemperature, double u, double E,  double leafWidth = 1.0) {
  if(NumericVector::is_na(SWRabs)) SWRabs = 0.0;
  if(NumericVector::is_na(LWRnet)) LWRnet = 0.0;
  double lambda = meteoland::utils_latentHeatVaporisationMol(airTemperature);
  u = std::max(u, 0.1);//Force minimum wind speed to avoid excessive heating
  double gHa = 0.189*pow(u/(leafWidth*0.0072), 0.5);
  double gr = 4.0*0.97*SIGMA_W*pow(273.16+airTemperature,3.0)/Cp_Jmol;
  double deltaTemp = (SWRabs + LWRnet - (lambda*(E/2000.0)))/(Cp_Jmol*(gr+gHa));
  return(airTemperature+deltaTemp);
}

/*
 *  returns leaf vapour pressure in kPa
 */
//' @rdname biophysics
//' @param leafTemp Leaf temperature (ºC).
//' @param leafPsi Leaf water potential (MPa).
// [[Rcpp::export("biophysics_leafVapourPressure")]]
double leafVapourPressure(double leafTemp,  double leafPsi) {
  double vpsl = meteoland::utils_saturationVP(std::max(0.0,leafTemp));
  double vpl = vpsl*exp((2.17*leafPsi)/(leafTemp+273.15));
  return(vpl);
}

/**
 * Converts irradiance units (W*m-2) to quantum flux (micromol * m-2 * s-1), 
 * defined as the number of photons (in micromol) per second and unit area
 */
//' @rdname biophysics
//' @param I Irradiance (in W*m-2).
//' @param lambda Wavelength (in nm).
// [[Rcpp::export("biophysics_irradianceToPhotonFlux")]]
double irradianceToPhotonFlux(double I, double lambda = 546.6507) {
  return(I*lambda*0.836*1e-2);
}

/**
 * Vogel equation for liquid dynamic viscosity (= 1 for 20ºC) (= 1/fluidity)
 * 
 * Hervé Cochard. A new mechanism for tree mortality due to drought and heatwaves. BioRxiv, 2019,
 *  pp.531632. ff10.1101/531632ff. ffhal-02273372f
 *  empirical equation for fluidity = 1/viscosity = 1.012e-4*pow(T,2.0) + 2.042e-2*T + 5.518e-1
 *  where T in degrees C
 */
//' @rdname biophysics
//' @param temp Temperature (ºC).
// [[Rcpp::export("biophysics_waterDynamicViscosity")]]
double waterDynamicViscosity(double temp) {
  return(exp(-3.7188+(578.919/(-137.546+ temp + 273.15))));
}
