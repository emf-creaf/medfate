#include <numeric>
#include <math.h>
#include "radiation_c.h"
#include "biophysicsutils_c.h"

double airDensity_c(double temperature, double Patm) {
  return((Patm/(1.01*(temperature+273.16)*0.287)));
}

double saturationVapourPressure_c(double temperature) {
  return(0.61078 * exp(17.269 * temperature/(temperature + 237.3)));
}
double saturationVaporPressureCurveSlope_c(double temperature) {
  return(4098.0 * (0.6108 * exp((17.27 * temperature)/(temperature + 237.3)))/pow(temperature + 237.3,2.0));
}

double averageDailyVapourPressure_c(double Tmin, double Tmax, double RHmin, double RHmax) {
  double vs_Tmax = saturationVapourPressure_c(Tmax);
  double vs_Tmin = saturationVapourPressure_c(Tmin);
  return((vs_Tmin * (RHmax/100.0) + vs_Tmax * (RHmin/100.0))/2.0);
}
double averageDaylightTemperature_c(double Tmin, double Tmax) {
  return(0.606*Tmax + 0.394*Tmin);
}

double atmosphericPressure_c(double elevation) {
  return(101.32500*pow(1.0-2.2569*pow(10.0,-5)*elevation,5.2353));
}

double latentHeatVaporisation_c(double temperature) {
  return(2.5023 - 0.00243054*temperature);
}

double latentHeatVaporisationMol_c(double temperature) {
  return(latentHeatVaporisation_c(temperature)*pow(10.0,6.0)*0.018);
}

double psychrometricConstant_c(double temperature, double Patm) {
  return((0.00163*Patm)/latentHeatVaporisation_c(temperature));
}
double PenmanPET_c(double latrad, double elevation, double slorad, double asprad, int J,
                   double Tmin, double Tmax, double RHmin, double RHmax, double R_s,
                   double u, double z, double z0,
                   double alpha, std::string windfun) {
  double Tday = (Tmax + Tmin)/2.0;
  double RHmean = (RHmax + RHmin)/2.0;
  //Atmospheric pressure kPa
  double Patm = atmosphericPressure_c(elevation);
  //Slope of the saturation vapor pressure curve kPa.Celsius^-1
  double delta = saturationVaporPressureCurveSlope_c(Tday);
  //Latent heat of vaporisation MJ.kg^-1
  double lambda = latentHeatVaporisation_c(Tday);
  //Psychrometric constant kPa.Celsius^-1
  double gamma = psychrometricConstant_c(Tday, Patm);
  //Solar declination
  double delta2 =  solarDeclination_c(J);
  //Solar constant
  double Gsc = solarConstant_c(J);
  // double d_r2 = 1.0 + 0.033 * cos(2.0 * PI/365.0 * ((double)J));
  // double w_s = acos(-tan(latitude) * tan(delta2));
  // double N = 24.0/PI*w_s; (N not used)
  
  double PET = 0.0;
  if(!std::isnan(u)) {
    double u2 = u * log(2.0/z0)/log(z/z0);
    double vs_Tmax = saturationVapourPressure_c(Tmax);
    double vs_Tmin = saturationVapourPressure_c(Tmin);
    double vas = (vs_Tmax + vs_Tmin)/2.0;
    double vpa = (vs_Tmin * (RHmax/100.0) + vs_Tmax * (RHmin/100.0))/2.0;
    double R_n = netRadiation_c(Gsc, latrad, elevation, slorad, asprad, delta2,
                                vpa, Tmin, Tmax, R_s,
                                alpha);//Net radiation
    double Ea = (vas - vpa); //Saturation vapor deficit
    if(windfun == "1956") {
      Ea = Ea*(1.313 + 1.381 * u2);
    } else if(windfun == "1948") {
      Ea = Ea*(2.626 + 1.381*u2);
    }
    PET = delta/(delta + gamma) * (R_n/lambda) + gamma/(delta + gamma) * Ea;
  } else {
    //Equation by Valiantzas (2006, eq 33) for situations where wind is not available
    //Valiantzas JD (2006) Simplified versions for the Penman evaporation equation using routine weather data. Journal of Hydrology 331, 690–702. doi:10.1016/j.jhydrol.2006.06.012.
    double R_a = RpotDay_c(Gsc, latrad,  slorad, asprad, delta2); //Extraterrestrial (potential) radiation
    double R_ratio = R_s/R_a;
    if(R_a==0.0) R_ratio = 0.0;
    R_ratio = std::min(std::max(R_ratio, 0.0), 1.0);
    double wf = 0.09;
    if(windfun == "1956") {
      wf = 0.06;
    } else if(windfun == "1948") {
      wf = 0.09;
    }
    PET = 0.047 * R_s * sqrt(Tday + 9.5) - 2.4 * pow(R_ratio,2.0) + wf * (Tday + 20.0) * (1.0 - RHmean/100.0);
  }
  if(PET<0.0) PET = 0.0;
  return(PET);
}


//' Physical and biophysical utility functions
//' 
//' Internal utility functions for the calculation of biophysical variables. 
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
//' @keywords internal
// [[Rcpp::export("biophysics_radiationDiurnalPattern")]]
double radiationDiurnalPattern_c(double t, double daylength) {
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
//' @keywords internal
// [[Rcpp::export("biophysics_temperatureDiurnalPattern")]]
double temperatureDiurnalPattern_c(double t, double tmin, double tmax, 
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

/**
 * Leaf temperature from energy balance
 * 
 * Campbell, G. S., and J. M. Norman. 1998. An introduction to environmental biophysics: 2nd edition. (eqns. 14.1 & 14.3)
 */
//' @rdname biophysics
//' 
//' @param u Wind speed above the leaf boundary layer (in m/s).
//' @param airTemperature Air temperature (in ºC).
//' @param absRad Absorbed long- and short-wave radiation (in W·m-2).
//' @param E Transpiration flow (in mmol H20·m-2·s-1) per one sided leaf area basis.
//' @param leafWidth Leaf width (in cm).
//' 
//' @keywords internal
// [[Rcpp::export("biophysics_leafTemperature")]]
double leafTemperature_c(double absRad, double airTemperature, double u, double E,  double leafWidth = 1.0) {
  double lambda = latentHeatVaporisationMol_c(airTemperature);
  u = std::max(u, 0.1);//Force minimum wind speed to avoid excessive heating
  double gHa = 0.189*pow(u/(leafWidth*0.0072), 0.5);
  double gr = 4.0*0.97*SIGMA_Wm2*pow(273.16+airTemperature,3.0)/Cp_Jmol;
  double deltaTemp = (absRad- (0.97*SIGMA_Wm2*pow(273.16+airTemperature,4.0)) - (lambda*(E/2000.0)))/(Cp_Jmol*(gr+gHa));
  return(airTemperature+deltaTemp);
}

//' @rdname biophysics
//' @param SWRabs Absorbed short-wave radiation (in W·m-2).
//' @param LWRnet Net long-wave radiation balance (in W·m-2).
//' @keywords internal
// [[Rcpp::export("biophysics_leafTemperature2")]]
double leafTemperature2_c(double SWRabs, double LWRnet, double airTemperature, double u, double E,  double leafWidth = 1.0) {
  if(std::isnan(SWRabs)) SWRabs = 0.0;
  if(std::isnan(LWRnet)) LWRnet = 0.0;
  double lambda = latentHeatVaporisationMol_c(airTemperature);
  u = std::max(u, 0.1);//Force minimum wind speed to avoid excessive heating
  double gHa = 0.189*pow(u/(leafWidth*0.0072), 0.5);
  double gr = 4.0*0.97*SIGMA_Wm2*pow(273.16+airTemperature,3.0)/Cp_Jmol;
  double deltaTemp = (SWRabs + LWRnet - (lambda*(E/2000.0)))/(Cp_Jmol*(gr+gHa));
  return(airTemperature+deltaTemp);
}

/*
 *  returns leaf vapour pressure in kPa
 */
//' @rdname biophysics
//' @param leafTemp Leaf temperature (ºC).
//' @param leafPsi Leaf water potential (MPa).
//' @keywords internal
// [[Rcpp::export("biophysics_leafVapourPressure")]]
double leafVapourPressure_c(double leafTemp,  double leafPsi) {
  double vpsl = saturationVapourPressure_c(std::max(0.0,leafTemp));
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
//' @keywords internal
// [[Rcpp::export("biophysics_irradianceToPhotonFlux")]]
double irradianceToPhotonFlux_c(double I, double lambda = 546.6507) {
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
//' @keywords internal
// [[Rcpp::export("biophysics_waterDynamicViscosity")]]
double waterDynamicViscosity_c(double temp) {
  return(exp(-3.7188+(578.919/(-137.546+ temp + 273.15))));
}

