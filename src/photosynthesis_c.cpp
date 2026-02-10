#include <RcppArmadillo.h>
#include <numeric>
#include <math.h>
#include "biophysicsutils_c.h"
#include "meteoland/utils_c.hpp"
#include "photosynthesis_c.h"


/**
 * Species-independent photosynthesis terms depenent of leaf temperature
 *  From:
 *  Bernacchi, C. J., E. L. Singsaas, C. Pimentel, A. R. Portis, and S. P. Long. 2001. Improved temperature response functions for models of Rubisco-limited photosynthesis. 
 *  Plant, Cell and Environment 24:253–259.
 *  
 *  Tleaf - Leaf temperature (ºC)
 *  Oi - Oxigen concentration (mmol*mol-1)
 */
//Compensation point (micromol * mol-1)
//' Photosynthesis submodel functions
//' 
//' Set of functions used in the calculation of photosynthesis
//' 
//' @param Tleaf Leaf temperature (in ºC).
//' @param Oi Oxigen concentration (mmol*mol-1).
//' @param Vmax298,Vmax298SL,Vmax298SH Maximum Rubisco carboxylation rate per leaf area at 298ºK (i.e. 25 ºC) (micromol*s-1*m-2) (for each canopy layer in the case of \code{photo_multilayerPhotosynthesisFunction}). 'SH' stands for shade leaves, whereas 'SL' stands for sunlit leaves.
//' @param Jmax298,Jmax298SL,Jmax298SH Maximum electron transport rate per leaf area at 298ºK (i.e. 25 ºC) (micromol*s-1*m-2) (for each canopy layer in the case of \code{photo_multilayerPhotosynthesisFunction}). 'SH' stands for shade leaves, whereas 'SL' stands for sunlit leaves.
//' @param Q Active photon flux density (micromol * s-1 * m-2).
//' @param Ci CO2 internal concentration (micromol * mol-1).
//' @param GT CO2 saturation point corrected by temperature (micromol * mol-1).
//' @param Jmax Maximum electron transport rate per leaf area (micromol*s-1*m-2).
//' @param Km Km = Kc*(1.0+(Oi/Ko)) - Michaelis-Menten term corrected by temperature (in micromol * mol-1).
//' @param Vmax Maximum Rubisco carboxylation rate per leaf area (micromol*s-1*m-2).
//' @param Catm CO2 air concentration (micromol * mol-1).
//' @param Gc CO2 leaf (stomatal) conductance (mol * s-1 * m-2).
//' @param E Transpiration flow rate per leaf area (mmol*s-1*m-2).
//' @param psiLeaf Leaf water potential (MPa).
//' @param Patm Atmospheric air pressure (in kPa).
//' @param Tair Air temperature (in ºC).
//' @param vpa Vapour pressure deficit (in kPa).
//' @param u Wind speed above the leaf boundary (in m/s) (for each canopy layer in the case of \code{photo_multilayerPhotosynthesisFunction}).
//' @param absRad Absorbed long- and short-wave radiation (in W*m^-2).
//' @param SWRabs Absorbed short-wave radiation (in W·m-2).
//' @param LWRnet Net long-wave radiation balance (in W·m-2).
//' @param leafWidth Leaf width (in cm).
//' @param refLeafArea Leaf reference area.
//' @param verbose Boolean flag to indicate console output.
//' @param SLarea,SHarea Leaf area index of sunlit/shade leaves (for each canopy layer in the case of \code{photo_multilayerPhotosynthesisFunction}).
//' @param absRadSL,absRadSH Instantaneous absorbed radiation (W·m-2) per unit of sunlit/shade leaf area (for each canopy layer in the case of \code{photo_multilayerPhotosynthesisFunction}).
//' @param QSL,QSH Active photon flux density (micromol * s-1 * m-2) per unit of sunlit/shade leaf area (for each canopy layer in the case of \code{photo_multilayerPhotosynthesisFunction}).
//' 
//' @details Details of the photosynthesis submodel are given in the medfate book
//' 
//' @return
//' Values returned for each function are:
//' \itemize{
//'   \item{\code{photo_GammaTemp}: CO2 compensation concentration (micromol * mol-1).}
//'   \item{\code{photo_KmTemp}: Michaelis-Menten coefficients of Rubisco for Carbon (micromol * mol-1) and Oxigen (mmol * mol-1).}
//'   \item{\code{photo_VmaxTemp}: Temperature correction of Vmax298.}
//'   \item{\code{photo_JmaxTemp}: Temperature correction of Jmax298.}
//'   \item{\code{photo_electronLimitedPhotosynthesis}: Electron-limited photosynthesis (micromol*s-1*m-2) following Farquhar et al. (1980).}
//'   \item{\code{photo_rubiscoLimitedPhotosynthesis}: Rubisco-limited photosynthesis (micromol*s-1*m-2) following Farquhar et al. (1980).}
//'   \item{\code{photo_photosynthesis}: Calculates gross photosynthesis (micromol*s-1*m-2) following (Farquhar et al. (1980) and Collatz et al (1991).}
//'   \item{\code{photo_leafPhotosynthesisFunction}: Returns a data frame with the following columns:
//'     \itemize{
//'       \item{\code{LeafTemperature}: Leaf temperature (ºC).}
//'       \item{\code{LeafVPD}: Leaf vapor pressure deficit (kPa).}
//'       \item{\code{LeafCi}: Internal CO2 concentration (micromol * mol-1).}
//'       \item{\code{Gsw}: Leaf stomatal conductance to water vapor (mol * s-1 * m-2).}
//'       \item{\code{GrossPhotosynthesis}: Gross photosynthesis (micromol*s-1*m-2).}
//'       \item{\code{NetPhotosynthesis}: Net photosynthesis, after discounting autotrophic respiration (micromol*s-1*m-2).}
//'     }
//'   }
//'   \item{\code{photo_sunshadePhotosynthesisFunction}: Returns a data frame with the following columns:
//'     \itemize{
//'       \item{\code{GrossPhotosynthesis}: Gross photosynthesis (micromol*s-1*m-2).}
//'       \item{\code{NetPhotosynthesis}: Net photosynthesis, after discounting autotrophic respiration (micromol*s-1*m-2).}
//'       \item{\code{LeafCiSL}: Sunlit leaf internal CO2 concentration (micromol * mol-1).}
//'       \item{\code{LeafCiSH}: Shade leaf internal CO2 concentration (micromol * mol-1).}
//'       \item{\code{LeafTempSL}: Sunlit leaf temperature (ºC).}
//'       \item{\code{LeafTempSH}: Shade leaf temperature (ºC).}
//'       \item{\code{LeafVPDSL}: Sunlit leaf vapor pressure deficit (kPa).}
//'       \item{\code{LeafVPDSH}: Shade leaf vapor pressure deficit (kPa).}
//'     }
//'   }
//'   \item{\code{photo_multilayerPhotosynthesisFunction}: Return a data frame with the following columns:
//'     \itemize{
//'       \item{\code{GrossPhotosynthesis}: Gross photosynthesis (micromol*s-1*m-2).}
//'       \item{\code{NetPhotosynthesis}: Net photosynthesis, after discounting autotrophic respiration (micromol*s-1*m-2).}
//'     }
//'   }
//' }
//' 
//' @references
//' Bernacchi, C. J., E. L. Singsaas, C. Pimentel, A. R. Portis, and S. P. Long. 2001. Improved temperature response functions for models of Rubisco-limited photosynthesis. Plant, Cell and Environment 24:253–259.
//' 
//' Collatz, G. J., J. T. Ball, C. Grivet, and J. A. Berry. 1991. Physiological and environmental regulation of stomatal conductance, photosynthesis and transpiration: a model that includes a laminar boundary layer. Agricultural and Forest Meteorology 54:107–136.
//' 
//' Farquhar, G. D., S. von Caemmerer, and J. A. Berry. 1980. A biochemical model of photosynthetic CO2 assimilation in leaves of C3 species. Planta 149:78–90.
//' 
//' Leuning, R. 2002. Temperature dependence of two parameters in a photosynthesis model. Plant, Cell and Environment 25:1205–1210.
//' 
//' Sperry, J. S., M. D. Venturas, W. R. L. Anderegg, M. Mencuccini, D. S. Mackay, Y. Wang, and D. M. Love. 2016. Predicting stomatal responses to the environment from the optimization of photosynthetic gain and hydraulic cost. Plant Cell and Environment.
//' 
//' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
//' 
//' @seealso
//' \code{\link{hydraulics_supplyFunctionNetwork}}, \code{\link{biophysics_leafTemperature}}, \code{\link{spwb}}
//' 
//' @name photo
// [[Rcpp::export("photo_GammaTemp")]]
double gammaTemp_c(double Tleaf) {return(42.75*exp((37830*(Tleaf-25.0))/(298.0*R_gas*(Tleaf+273))));} 



//Michaelis-Menten coefficients of Rubisco for Carbon (micromol * mol-1)
double KcTemp_c(double Tleaf) {return(404.9*exp((79430*(Tleaf-25.0))/(298.0*R_gas*(Tleaf+273))));}

//Michaelis-Menten coefficients of Rubisco for Oxigen (mmol * mol-1)
double KoTemp_c(double Tleaf) {return(278.4*exp((36380*(Tleaf-25.0))/(298.0*R_gas*(Tleaf+273))));}

//' @rdname photo
// [[Rcpp::export("photo_KmTemp")]]
double KmTemp_c(double Tleaf, double Oi) {
   double Kc = KcTemp_c(Tleaf);
   double Ko = KoTemp_c(Tleaf);  
   return(Kc*(1.0+(Oi/Ko)));
}


/**
 * Temperature correction of Vmax298
 * 
 * From:
 *    Leuning, R. 2002. Temperature dependence of two parameters in a photosynthesis model. 
 *    Plant, Cell and Environment 25:1205–1210.
 *    
 * Eq.1 with parameters from Table 2
 * 
 *  Tleaf - Leaf temperature (ºC)
 *  Vmax298 - maximum carboxylation rate at 298ºK (ie. 25 ºC) (micromol*s-1*m-2)
 */
//' @rdname photo
// [[Rcpp::export("photo_VmaxTemp")]]
double VmaxTemp_c(double Vmax298, double Tleaf) {
  if(std::isnan(Vmax298)) Vmax298 = 0.0;
  double Ha = 73637.0; //Energy of activation J * mol-1
  double Hd = 149252.0; //Energy of deactivation J * mol-1
  double Sv = 486.0;  //Entropy term J * mol-1 * K-1
  double C = 1.0+exp((Sv*298.2-Hd)/(R_gas*298.2));
  return(Vmax298*(C*exp((Ha/(R_gas*298.2))*(1.0-298.2/(Tleaf+273.2))))/(1.0+exp((Sv*Tleaf-Hd)/(R_gas*(Tleaf+273.2)))));
}


/**
 * Temperature correction of Jmax298
 * 
 * From:
 *    Leuning, R. 2002. Temperature dependence of two parameters in a photosynthesis model. 
 *    Plant, Cell and Environment 25:1205–1210.
 *    
 * Eq.1 with parameters from Table 2
 * 
 *  Tleaf - Leaf temperature (ºC)
 *  Jmax298 - maximum electron transport rate at 298ºK (ie. 25 ºC) (micromol*s-1*m-2)
 */
//' @rdname photo
// [[Rcpp::export("photo_JmaxTemp")]]
double JmaxTemp_c(double Jmax298, double Tleaf) {
  double Ha = 50300.0; //Energy of activation J * mol-1
  double Hd = 152044.0; //Energy of deactivation J * mol-1
  double Sv = 495.0;  //Entropy term J * mol-1 * K-1
  double C = 1.0+exp((Sv*298.2-Hd)/(R_gas*298.2));
  return(Jmax298*(C*exp((Ha/(R_gas*298.2))*(1.0-298.2/(Tleaf+273.2))))/(1.0+exp((Sv*Tleaf-Hd)/(R_gas*(Tleaf+273.2)))));
}


// Returns canopy conductance in mol
double gCrown_c(double u, double gCrown0){
  u=  std::max(0.1, u); //# to avoid very high conductance values 
  return(gCrown0*pow(u,0.6));
}

// Returns leaf boundary layer conductance in mol
double gLeafBoundary_c(double u, double leafWidth, double gBound0){
  return(gBound0*pow(u/(leafWidth*0.0072), 0.5));
}


/**
 * Calculates electron-limited photosynthesis (Farquhar et al. 1980)
 * 
 * Ci - CO2 internal concentration (micromol * mol-1)
 * GT - CO2 saturation point corrected by temperature (micromol * mol-1)
 * Q - Active photon flux density (micromol * s-1 * m-2)
 * Jmax - maximum electron transport rate per leaf area (micromol*s-1*m-2)
 * 
 * return units: micromol*s-1*m-2
 */
//' @rdname photo
//' @keywords internal
// [[Rcpp::export("photo_electronLimitedPhotosynthesis")]]
double electronLimitedPhotosynthesis_c(double Q, double Ci, double GT, double Jmax) {
  double J = ((quantumYield*Q+Jmax)-sqrt(pow(quantumYield*Q+Jmax, 2.0) - 4.0*lightResponseCurvature*quantumYield*Q*Jmax))/(2.0*lightResponseCurvature);
  return((J/4.0)*((Ci-GT)/(Ci+2.0*GT)));
}
double electronLimitedPhotosynthesisDerivative_c(double Q, double Ci, double GT, double Jmax){
  double J = ((quantumYield*Q+Jmax)-sqrt(pow(quantumYield*Q+Jmax, 2.0) - 4.0*lightResponseCurvature*quantumYield*Q*Jmax))/(2.0*lightResponseCurvature);
  return((J/4.0)*((3.0*GT)/pow(Ci+2.0*GT,2.0)));
}

/**
 * Calculates rubisco-limited photosynthesis (Farquhar et al. 1980)
 * 
 * Ci - CO2 internal concentration (micromol * mol-1)
 * GT - CO2 saturation point corrected by temperature (micromol * mol-1)
 * Km = Kc*(1.0+(Oi/Ko)) - Michaelis-Menten term corrected by temperature (in micromol * mol-1)
 * Vmax - maximum Rubisco carboxylation rate per leaf area (micromol*s-1*m-2)
 * 
 * return units: micromol*s-1*m-2
 */
//' @rdname photo
//' @keywords internal
// [[Rcpp::export("photo_rubiscoLimitedPhotosynthesis")]]
double rubiscoLimitedPhotosynthesis_c(double Ci, double GT, double Km, double Vmax) {
  return(Vmax *(Ci-GT)/(Ci+Km));
}
double rubiscoLimitedPhotosynthesisDerivative_c(double Ci, double GT, double Km, double Vmax) {
  return(Vmax *(Km+GT)/pow(Ci+Km,2.0));
}

/**
 * Calculates photosynthesis (Farquhar et al. 1980/Collatz et al 1991)
 * 
 * Ci - CO2 internal concentration (micromol * mol-1)
 * GT - CO2 saturation point corrected by temperature (micromol * mol-1)
 * Km = Kc*(1.0+(Oi/Ko)) - Michaelis-Menten term corrected by temperature (in micromol * mol-1)
 * Q - Active photon flux density (micromol * s-1 * m-2)
 * Jmax - maximum electron transport rate per leaf area (micromol*s-1*m-2) corrected  by temperature
 * Vmax - maximum Rubisco carboxylation rate per leaf area (micromol*s-1*m-2) corrected  by temperature
 * 
 * return units: micromol*s-1*m-2
 */
double photosynthesis_Ci_c(double Q, double Ci, double GT, double Km, double Vmax, double Jmax) {
  double Je = electronLimitedPhotosynthesis_c(Q, Ci, GT, Jmax);
  double Jc = rubiscoLimitedPhotosynthesis_c(Ci, GT, Km, Vmax);
  return(std::max(0.0,(Je+Jc-sqrt(pow(Je+Jc,2.0)-4.0*0.98*Je*Jc))/(2.0*0.98)));
}



// Auxiliary functions for Newton-Raphson
double f_c(double x, double Q, double Ca, double Gc, double GT, double Km, double Vmax, double Jmax) {
  return(photosynthesis_Ci_c(Q,x, GT, Km, Vmax, Jmax)-(Gc*(Ca-x)));
}
double fder_c(double x, double Q, double Ca, double Gc, double GT, double Km, double Vmax, double Jmax) {
  double Je = electronLimitedPhotosynthesis_c(Q, x, GT, Jmax);
  double dJe = electronLimitedPhotosynthesisDerivative_c(Q, x, GT, Jmax);
  double Jc = rubiscoLimitedPhotosynthesis_c(x, GT, Km, Vmax);
  double dJc = rubiscoLimitedPhotosynthesisDerivative_c(x, GT, Km, Vmax);
  double dA1 = (1.0/(2.0*0.98))*(dJe+dJc-(0.5*pow(pow(Je+Jc,2.0)-4.0*0.98*Je*Jc,-0.5)*(2.0*Je*dJe+2.0*Jc*dJc+(2.0-4.0*0.98)*(dJe*Jc + dJc*Je))));
  double dA2 = -Gc;
  return(dA1-dA2);
}

void leafphotosynthesis_inner_c(Photo& photo, 
                                double Q, double Catm, double Gc, double Tleaf, double Vmax298, double Jmax298) {
  //Corrections per leaf temperature
  double GT = gammaTemp_c(Tleaf);
  double Km = KmTemp_c(Tleaf, O2_conc);
  double Vmax = VmaxTemp_c(Vmax298, Tleaf);
  double Jmax = JmaxTemp_c(Jmax298, Tleaf);
  double x,x1,e,fx,fx1;
  x1 = 0.0;//initial guess
  e = 0.001; // accuracy in micromol * mol-1
  int cnt = 0;
  int mxiter = 100;
  do {
    x=x1; /*make x equal to the last calculated value of                             x1*/
    fx=f_c(x, Q, Catm, Gc, GT, Km, Vmax, Jmax);            //simplifying f(x)to fx
    fx1=fder_c(x, Q, Catm, Gc, GT, Km, Vmax, Jmax);            //simplifying fprime(x) to fx1
    x1=x-(fx/fx1);/*calculate x{1} from x, fx and fx1*/ 
    cnt++;
  } while ((std::abs(x1-x)>=e) && (cnt < mxiter));
  photo.Ci = x1;
  photo.A = photosynthesis_Ci_c(Q,x1,GT,Km,Vmax,Jmax);
}

double third_cubic_root(double p, double q, double r) {
  // Terms of the solution
  double Q = (pow(p, 2.0) - 3.0*q)/9.0;
  double R = (2.0*pow(p, 3.0) - 9.0*p*q + 27.0*r)/54.0;
  // bool sol = pow(R,2.0) < pow(Q, 3.0);
  double theta = std::acos(R/pow(Q, 3.0/2.0));
  //Third root of the cubic equation (Numerical Recipes in C, Press et al. 1989)
  double x3 = -2.0*pow(Q, 0.5)*cos((theta - 2.0*M_PI)/3.0) - (p/3.0);
  return(x3);
}

void photosynthesisBaldocchi_inner_c(BaldocchiPhoto &photoOut,
                                     double Q, 
                                     double Catm, 
                                     double Tleaf, 
                                     double u,
                                     double Vmax298, 
                                     double Jmax298, 
                                     double leafWidth,
                                     double Gsw_AC_slope,
                                     double Gsw_AC_intercept) {
  double Vmax = VmaxTemp_c(Vmax298, Tleaf);
  double Jmax = JmaxTemp_c(Jmax298, Tleaf);
  //Dark respiration
  double Rd = 0.015*Vmax;
  //Boundary layer conductance (UNITS in mol water per s-1 m-2!!!)
  double Gbound = gLeafBoundary_c(u, leafWidth); // mol boundary layer conductance
  double Gbc = Gbound/1.6;
  
  //Translate stomatal model to carbon units
  double b_prime = Gsw_AC_intercept/1.6;
  double mrh = Gsw_AC_slope/1.6;
  
  //Compensation point
  double Gamma_comp = gammaTemp_c(Tleaf);
  
  //Coefficients when Rubisco limits (c)
  double a_c = Vmax;
  double b_c = KmTemp_c(Tleaf, O2_conc);  
  double d_c = Gamma_comp;
  double e_c = 1.0;
  
  //Coefficients when electron transport limits (j)  
  double J = ((quantumYield*Q+Jmax)-sqrt(pow(quantumYield*Q+Jmax, 2.0) - 4.0*lightResponseCurvature*quantumYield*Q*Jmax))/(2.0*lightResponseCurvature);  
  double a_j = J;
  double b_j = 8.0*Gamma_comp;
  double d_j = Gamma_comp;
  double e_j = 4.0;
  
  //Solving for An when Rubisco limits (c) and when electron transport limits (j)
  double alpha = 1.0 + (b_prime/Gbc) - mrh;
  double beta = Catm*(Gbc*mrh - 2.0*b_prime - Gbc);
  double gamma = pow(Catm, 2.0)*b_prime*Gbc;
  double theta = Gbc*mrh - b_prime;
  double p_c = (e_c*beta + b_c*theta - a_c*alpha + e_c*alpha*Rd)/(e_c*alpha);
  double p_j = (e_j*beta + b_j*theta - a_j*alpha + e_j*alpha*Rd)/(e_j*alpha);
  double q_c = ((e_c*gamma) + (b_c*gamma/Catm) - (a_c*beta) + (a_c*d_c*theta) + (e_c*Rd*beta) + (Rd*b_c*theta))/(e_c*alpha);
  double q_j = ((e_j*gamma) + (b_j*gamma/Catm) - (a_j*beta) + (a_j*d_j*theta) + (e_j*Rd*beta) + (Rd*b_j*theta))/(e_j*alpha);
  double r_c = ((-1.0*a_c*gamma) + (a_c*d_c*gamma/Catm) + (e_c*Rd*gamma) + (Rd*b_c*gamma/Catm))/(e_c*alpha);
  double r_j = ((-1.0*a_j*gamma) + (a_j*d_j*gamma/Catm) + (e_j*Rd*gamma) + (Rd*b_j*gamma/Catm))/(e_j*alpha);
  double An_c = third_cubic_root(p_c, q_c, r_c);
  double An_j = third_cubic_root(p_j, q_j, r_j);
  //Take the minimum as An
  double An = std::min(An_c, An_j);
  //Stomatal CO2
  double Cs = Catm - (An/Gbc);
  //Stomatal conductance (for carbon)
  double Gsc = (mrh*An/Cs) + b_prime;
  //Stomatal conductance (for water)
  double Gsw = Gsc*1.6;
  //Internal CO2
  double Ci = Cs - (An/Gsc);
  //Gross photosynthesis
  double Ag = An + Rd;
  photoOut.Gsw = Gsw;
  photoOut.Cs = Cs;
  photoOut.Ci = Ci;
  photoOut.An = An;
  photoOut.Ag = Ag;
}


void leafPhotosynthesisOneFunction2_c(PhotoFunction& photo, 
                                      double E, double psiLeaf, double Catm, double Patm, double Tair, double vpa, double u, 
                                      double SWRabs, double LWRnet, double Q, double Vmax298, double Jmax298, 
                                      double leafWidth, double refLeafArea) {
  double leafTemp, leafVPD;
  double Gwdiff, Gbound;
  leafTemp = leafTemperature2_c(SWRabs/refLeafArea, LWRnet/refLeafArea, Tair, u, E, leafWidth);
  leafVPD = std::max(0.0,leafVapourPressure_c(leafTemp, psiLeaf) - vpa);
  // Separates diffusive conductance into stomatal and boundary layer conductance
  Gwdiff = Patm*(E/1000.0)/leafVPD; //Transform flow from mmol to mol
  Gbound = gLeafBoundary_c(u, leafWidth); // mol boundary layer conductance
  Gwdiff = std::min(Gwdiff, Gbound); //Diffusive conductance cannot be smaller than the boundary layer resistances
  Photo LP; 
  leafphotosynthesis_inner_c(LP, Q/refLeafArea, Catm, Gwdiff/1.6, std::max(0.0,leafTemp), Vmax298/refLeafArea, Jmax298/refLeafArea);
  photo.Gsw  = std::abs(1.0/((1.0/Gwdiff) - (1.0/Gbound))); //Determine stomatal conductance after accounting for leaf boundary conductance
  photo.Ci = LP.Ci;
  photo.GrossPhotosynthesis = LP.A;
  photo.NetPhotosynthesis = LP.A - 0.015*VmaxTemp_c(Vmax298/refLeafArea, leafTemp);
  photo.LeafTemperature = leafTemp;
  photo.LeafVPD = leafVPD;
}

