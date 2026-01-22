#include <math.h>
#include "fireseverity_c.h"

double erfInv_c(double x){
  if(abs(x) < 6e-3) return(x);
  double tt1, tt2, lnx, sgn;
  sgn = (x < 0.0) ? -1.0 : 1.0;
  
  x = (1.0 - x)*(1.0 + x);        // x = 1 - x*x;
  lnx = logf(x);
  
  tt1 = 2.0/(3.141592*0.147) + 0.5 * lnx;
  tt2 = 1.0/(0.147) * lnx;
  
  return(sgn*sqrtf(-tt1 + sqrtf(tt1*tt1 - tt2)));
}

//' Fire severity functions
//' 
//' Functions to estimate fire effects on foliage, buds and cambium, based on the model
//' by Michaletz & Johnson (2008)
//' 
//' @param Ib_surf Surface fireline intensity (kW/m).
//' @param t_res fire residence time (seconds).
//' @param T_air Air temperature (degrees Celsius).
//' @param rho_air Air density (kg/m3).
//' @param rho_bark Bark density (kg/m3).
//' @param fmc_bark Bark moisture content (% dry weight).
//' @param z height (m).
//' @param SLA Specific leaf area (m2/kg).
//' @param h Heat transfer coefficient
//' @param c Specific heat capacity
//' @param thermal_factor Tissue thermal factor.
//' @param bark_diffusivity Bark thermal diffusivity (m2/s).
//' 
//' @return 
//' \itemize{
//'   \item{Function \code{fire_plumeTemperature} returns the plume temperature at a given height.}
//'   \item{Function \code{fire_barkThermalDiffusivity} returns the bark thermal diffusivity given a bark moisture value.}
//'   \item{Function \code{fire_radialBoleNecrosis} returns the depth of radial bole necrosis in cm.}
//'   \item{Function \code{fire_leafThermalFactor} returns the thermal factor of leaves as a function of specific leaf area.}
//'   \item{Function \code{fire_necrosisCriticalTemperature} returns the (plume) temperature yielding necrosis for a given residence time and tissue thermal factor.}
//'   \item{Function \code{fire_necrosisHeight} returns the height (in m) of necrosis for tissues with given thermal factor.}
//' }
//' 
//' @references
//' 
//'   Michaletz, S.T., and Johnson, E.A. 2006. A heat transfer model of crown scorch in forest fires. Can. J. For. Res. 36: 2839–2851. doi:10.1139/X06-158.
//' 
//'   Michaletz ST, Johnson EA. 2008. A biophysical process model of tree mortality in surface fires. Canadian Journal of Forest Research 38: 2013–2029.
//'   
//' @name fire_severity
//' @keywords internal
// [[Rcpp::export("fire_plumeTemperature")]]
double plumeTemperature_c(double Ib_surf, double z, double T_air = 25.0, double rho_air = 1.169) {
  double C = 2.6; //(Yuan and Cox 1996)
  double c_air = 1.007; //J·kg-1·ºC-1
  return(std::min(900.0, C*(1.0/z)*pow((T_air + 273.15)/9.8,1.0/3.0)*pow(Ib_surf/(c_air*rho_air), 2.0/3.0)+ T_air));
}
 
//' @rdname fire_severity
//' @keywords internal
// [[Rcpp::export("fire_barkThermalDiffusivity")]]
double barkThermalDiffusivity_c(double fmc_bark, double rho_bark = 500.0, double T_air = 25.0) {
  double W = fmc_bark/100.0; // From percent to fraction
  double c_bark = 1105.315 + 4.857*T_air + W*4180.0 + 348.342; // Specific heat capacity of bark, J·kg-1·ºC-1, Martin 1963 
  // Rcout<< c_bark<<"\n";
  double rho_bark_moist = (W + pow(W, 2.0))*rho_bark; // Moisture density kg/m-3
  // Rcout<< rho_bark_moist<<"\n";
  double k_bark = 1.0e-4*(2.104*rho_bark + 5.544*rho_bark_moist + 3.266*T_air - 166.216); //Thermal conductivity, Martin 1963
  // Rcout<< k_bark<<"\n";
  double alpha_bark = k_bark /(rho_bark*c_bark); //Thermal diffusivity
  return(alpha_bark); 
}

//' @rdname fire_severity
//' @keywords internal
// [[Rcpp::export("fire_radialBoleNecrosis")]]
double radialBoleNecrosis_c(double Ib_surf, double t_res, double bark_diffusivity,
                          double T_air = 25.0, double rho_air = 1.169) {
  double T_plume = plumeTemperature_c(Ib_surf, 0.1, T_air, rho_air);
  double theta = std::max(0.0, (temperatureNecrosis - T_plume)/(T_air - T_plume));
  double xn = 2.0*pow(bark_diffusivity*t_res, 0.5)*erfInv_c(theta);
  return(xn*100.0); // from m to cm
}

//' @rdname fire_severity
//' @keywords internal
// [[Rcpp::export("fire_leafThermalFactor")]]
double leafThermalFactor_c(double SLA, double h = 130.0, double c = 2500.0) {
   return(SLA*(h/c));
}

//' @rdname fire_severity
//' @keywords internal
// [[Rcpp::export("fire_necrosisCriticalTemperature")]]
double necrosisCriticalTemperature_c(double t_res, double thermal_factor, double T_air = 25.0) {
   double theta = exp(-1.0*thermal_factor*t_res);
   double T_c = (temperatureNecrosis - (theta*T_air))/(1.0 - theta);
   return(T_c);
}

//' @rdname fire_severity
//' @keywords internal
// [[Rcpp::export("fire_necrosisHeight")]]
double necrosisHeight_c(double Ib_surf, double t_res, double thermal_factor, 
                      double T_air = 25.0, double rho_air = 1.169) {
  double T_c = necrosisCriticalTemperature_c(t_res, thermal_factor, T_air);
  double C = 2.6; //(Yuan and Cox 1996)
  double c_air = 1.007; //J·kg-1·ºC-1
  double z_necrosis = C*(1.0/(T_c - T_air))*pow((T_air + 273.15)/9.8,1.0/3.0)*pow(Ib_surf/(c_air*rho_air), 2.0/3.0);
  return(z_necrosis);
}