#include <Rcpp.h>
#include <numeric>
#include <math.h>
#include "biophysicsutils.h"
#include "struct_photosynthesis.h"
#include <meteoland.h>
using namespace Rcpp;

const double R_gas = 8.314; //(J/mol/ºK) Universal gas constant
const double O2_conc = 209.0; //mmol*mol-1 (Collatz et al. 2001)
const double quantumYield = 0.3; //mol photon * mol-1 electron
const double lightResponseCurvature = 0.9;




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
double gammaTemp(double Tleaf) {return(42.75*exp((37830*(Tleaf-25.0))/(298.0*R_gas*(Tleaf+273))));} 

//Michaelis-Menten coefficients of Rubisco for Carbon (micromol * mol-1)
double KcTemp(double Tleaf) {return(404.9*exp((79430*(Tleaf-25.0))/(298.0*R_gas*(Tleaf+273))));}

//Michaelis-Menten coefficients of Rubisco for Oxigen (mmol * mol-1)
double KoTemp(double Tleaf) {return(278.4*exp((36380*(Tleaf-25.0))/(298.0*R_gas*(Tleaf+273))));}

//' @rdname photo
// [[Rcpp::export("photo_KmTemp")]]
double KmTemp(double Tleaf, double Oi = 209.0) {
  double Kc = KcTemp(Tleaf);
  double Ko = KoTemp(Tleaf);  
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
double VmaxTemp(double Vmax298, double Tleaf) {
  if(NumericVector::is_na(Vmax298)) Vmax298 = 0.0;
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
double JmaxTemp(double Jmax298, double Tleaf) {
  double Ha = 50300.0; //Energy of activation J * mol-1
  double Hd = 152044.0; //Energy of deactivation J * mol-1
  double Sv = 495.0;  //Entropy term J * mol-1 * K-1
  double C = 1.0+exp((Sv*298.2-Hd)/(R_gas*298.2));
  return(Jmax298*(C*exp((Ha/(R_gas*298.2))*(1.0-298.2/(Tleaf+273.2))))/(1.0+exp((Sv*Tleaf-Hd)/(R_gas*(Tleaf+273.2)))));
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
double electronLimitedPhotosynthesis(double Q, double Ci, double GT, double Jmax) {
  double J = ((quantumYield*Q+Jmax)-sqrt(pow(quantumYield*Q+Jmax, 2.0) - 4.0*lightResponseCurvature*quantumYield*Q*Jmax))/(2.0*lightResponseCurvature);
  return((J/4.0)*((Ci-GT)/(Ci+2.0*GT)));
}
double electronLimitedPhotosynthesisDerivative(double Q, double Ci, double GT, double Jmax){
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
double rubiscoLimitedPhotosynthesis(double Ci, double GT, double Km, double Vmax) {
  return(Vmax *(Ci-GT)/(Ci+Km));
}
double rubiscoLimitedPhotosynthesisDerivative(double Ci, double GT, double Km, double Vmax) {
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
double photosynthesis_Ci(double Q, double Ci, double GT, double Km, double Vmax, double Jmax) {
  double Je = electronLimitedPhotosynthesis(Q, Ci, GT, Jmax);
  double Jc = rubiscoLimitedPhotosynthesis(Ci, GT, Km, Vmax);
  return(std::max(0.0,(Je+Jc-sqrt(pow(Je+Jc,2.0)-4.0*0.98*Je*Jc))/(2.0*0.98)));
}

// Auxiliary functions for Newton-Raphson
double f(double x, double Q, double Ca, double Gc, double GT, double Km, double Vmax, double Jmax) {
  return(photosynthesis_Ci(Q,x, GT, Km, Vmax, Jmax)-(Gc*(Ca-x)));
}
double fder(double x, double Q, double Ca, double Gc, double GT, double Km, double Vmax, double Jmax) {
  double Je = electronLimitedPhotosynthesis(Q, x, GT, Jmax);
  double dJe = electronLimitedPhotosynthesisDerivative(Q, x, GT, Jmax);
  double Jc = rubiscoLimitedPhotosynthesis(x, GT, Km, Vmax);
  double dJc = rubiscoLimitedPhotosynthesisDerivative(x, GT, Km, Vmax);
  double dA1 = (1.0/(2.0*0.98))*(dJe+dJc-(0.5*pow(pow(Je+Jc,2.0)-4.0*0.98*Je*Jc,-0.5)*(2.0*Je*dJe+2.0*Jc*dJc+(2.0-4.0*0.98)*(dJe*Jc + dJc*Je))));
  double dA2 = -Gc;
  return(dA1-dA2);
}

/**
 * Calculates photosynthesis (Farquhar et al. 1980/Collatz et al 1991)
 * 
 * Q - Active photon flux density (micromol * s-1 * m-2)
 * Ca - CO2 air concentration (micromol * mol-1)
 * Gc - CO2 stomatal conductance (mol * s-1 * m-2)
 * Tleaf - Leaf temperature (ºC)
 * Jmax298 - maximum electron transport rate per leaf area at 298 ºK (i.e. 25 ºC) (micromol*s-1*m-2) 
 * Vmax298 - maximum Rubisco carboxylation rate per leaf area at 298 ºK (i.e. 25 ºC) (micromol*s-1*m-2) 
 * 
 * return units: micromol*s-1*m-2
 */
//' @rdname photo
//' @keywords internal
// [[Rcpp::export("photo_photosynthesis")]]
NumericVector leafphotosynthesis(double Q, double Catm, double Gc, double Tleaf, double Vmax298, double Jmax298, bool verbose=false) {
  //Corrections per leaf temperature
  double GT = gammaTemp(Tleaf);
  double Km = KmTemp(Tleaf, O2_conc);
  double Vmax = VmaxTemp(Vmax298, Tleaf);
  double Jmax = JmaxTemp(Jmax298, Tleaf);
  double x,x1,e,fx,fx1;
  x1 = 0.0;//initial guess
  e = 0.001; // accuracy in micromol * mol-1
  int cnt = 0;
  int mxiter = 100;
  if(verbose) Rcout <<"x{i}"<<"    "<<"x{i+1}"<<"        "<<"|x{i+1}-x{i}|\n";                
  do {
    x=x1; /*make x equal to the last calculated value of                             x1*/
    fx=f(x, Q, Catm, Gc, GT, Km, Vmax, Jmax);            //simplifying f(x)to fx
    fx1=fder(x, Q, Catm, Gc, GT, Km, Vmax, Jmax);            //simplifying fprime(x) to fx1
    x1=x-(fx/fx1);/*calculate x{1} from x, fx and fx1*/ 
    cnt++;
    if(verbose) Rcout<<x<<"     "<<x1<<"           "<<std::abs(x1-x)<<"\n";        
  } while ((std::abs(x1-x)>=e) && (cnt < mxiter));
  double A = photosynthesis_Ci(Q,x1,GT,Km,Vmax,Jmax);
  NumericVector res = NumericVector::create(x1, A);
  res.attr("names") = CharacterVector::create("Ci", "A");
  return(res);
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


// Returns canopy conductance in mol
double gCrown(double u, double gCrown0 = 0.150){
  u=  std::max(0.1, u); //# to avoid very high conductance values 
  return(gCrown0*pow(u,0.6));
}

// Returns leaf boundary layer conductance in mol
double gLeafBoundary(double u, double leafWidth, double gBound0 = 0.397){
  return(gBound0*pow(u/(leafWidth*0.0072), 0.5));
}

void photosynthesisBaldocchi_inner(BaldocchiPhoto &photoOut,
                                   double Q, 
                                   double Catm, 
                                   double Tleaf, 
                                   double u,
                                   double Vmax298, 
                                   double Jmax298, 
                                   double leafWidth,
                                   double Gsw_AC_slope,
                                   double Gsw_AC_intercept) {
  double Vmax = VmaxTemp(Vmax298, Tleaf);
  double Jmax = JmaxTemp(Jmax298, Tleaf);
  //Dark respiration
  double Rd = 0.015*Vmax;
  //Boundary layer conductance (UNITS in mol water per s-1 m-2!!!)
  double Gbound = gLeafBoundary(u, leafWidth); // mol boundary layer conductance
  double Gbc = Gbound/1.6;
  
  //Translate stomatal model to carbon units
  double b_prime = Gsw_AC_intercept/1.6;
  double mrh = Gsw_AC_slope/1.6;
  
  //Compensation point
  double Gamma_comp = gammaTemp(Tleaf);
  
  //Coefficients when Rubisco limits (c)
  double a_c = Vmax;
  double b_c = KmTemp(Tleaf, O2_conc);  
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

// From Baldocchi D (1994). An analytical solution for the coupled leaf photosynthesis and stomatal conductance models. Tree Physiology 14: 1069-1079 
//' @rdname photo
//' @param Gsw_AC_slope Slope of the An/C vs Gsw relationship 
//' @param Gsw_AC_intercept Intercept of the An/C vs Gsw relationship 
//' @keywords internal
// [[Rcpp::export("photo_photosynthesisBaldocchi")]]
NumericVector photosynthesisBaldocchi(double Q, 
                                      double Catm, 
                                      double Tleaf, 
                                      double u,
                                      double Vmax298, 
                                      double Jmax298, 
                                      double leafWidth,
                                      double Gsw_AC_slope,
                                      double Gsw_AC_intercept) {
  BaldocchiPhoto photoOut;
  photosynthesisBaldocchi_inner(photoOut, Q, Catm, Tleaf, u,Vmax298,Jmax298, leafWidth,Gsw_AC_slope,Gsw_AC_intercept);
  NumericVector res = {photoOut.Gsw, 
                       photoOut.Cs,
                       photoOut.Ci,
                       photoOut.An,
                       photoOut.Ag};
  res.attr("names") = CharacterVector::create("Gsw", "Cs" ,"Ci", "An", "Ag");
  return(res);
}

//' @rdname photo
//' @keywords internal
// [[Rcpp::export("photo_leafPhotosynthesisFunction")]]
DataFrame leafPhotosynthesisFunction(NumericVector E, NumericVector psiLeaf, double Catm, double Patm, double Tair, double vpa, double u, 
                             double absRad, double Q, double Vmax298, double Jmax298, 
                             double leafWidth = 1.0, double refLeafArea = 1.0, bool verbose = false) {
  int nsteps = E.size();
  NumericVector leafTemp(nsteps);
  NumericVector leafVPD(nsteps);
  NumericVector Gsw(nsteps), Ci(nsteps);
  NumericVector Ag(nsteps), An(nsteps);
  double Gwdiff, Gbound;
  for(int i=0;i<nsteps;i++){
    leafTemp[i] = leafTemperature(absRad/refLeafArea, Tair, u, E[i], leafWidth);
    leafVPD[i] = std::max(0.0,leafVapourPressure(leafTemp[i], psiLeaf[i]) - vpa);
    // Separates diffusive conductance into stomatal and boundary layer conductance
    Gwdiff = Patm*(E[i]/1000.0)/leafVPD[i]; //Transform flow from mmol to mol
    Gbound = gLeafBoundary(u, leafWidth); // mol boundary layer conductance
    Gwdiff = std::min(Gwdiff, Gbound); //Diffusive resistance cannot be smaller than the boundary layer resistance
    Gsw[i]  = std::abs(1.0/((1.0/Gwdiff) - (1.0/Gbound))); //Determine stomatal conductance after accounting for leaf boundary conductance
    NumericVector LP = leafphotosynthesis(Q/refLeafArea, Catm, Gwdiff/1.6, std::max(0.0,leafTemp[i]), Vmax298/refLeafArea, Jmax298/refLeafArea);
    Ci[i] = LP[0];
    Ag[i] = LP[1];
    An[i] = Ag[i] - 0.015*VmaxTemp(Vmax298/refLeafArea, leafTemp[i]);
  }
  return(DataFrame::create(Named("LeafTemperature") = leafTemp,
                      Named("LeafVPD") = leafVPD,
                      Named("Gsw") = Gsw,
                      Named("Ci") = Ci,
                      Named("GrossPhotosynthesis") = Ag,
                      Named("NetPhotosynthesis") = An));
}


NumericVector leafPhotosynthesisOneFunction2(double E, double psiLeaf, double Catm, double Patm, double Tair, double vpa, double u, 
                                             double SWRabs, double LWRnet, double Q, double Vmax298, double Jmax298, 
                                             double leafWidth = 1.0, double refLeafArea = 1.0, bool verbose = false) {
  double leafTemp, leafVPD, Gsw, Ci;
  double Ag, An;
  double Gwdiff, Gbound;
  leafTemp = leafTemperature2(SWRabs/refLeafArea, LWRnet/refLeafArea, Tair, u, E, leafWidth);
  leafVPD = std::max(0.0,leafVapourPressure(leafTemp, psiLeaf) - vpa);
  // Separates diffusive conductance into stomatal and boundary layer conductance
  Gwdiff = Patm*(E/1000.0)/leafVPD; //Transform flow from mmol to mol
  Gbound = gLeafBoundary(u, leafWidth); // mol boundary layer conductance
  Gwdiff = std::min(Gwdiff, Gbound); //Diffusive conductance cannot be smaller than the boundary layer resistances
  Gsw  = std::abs(1.0/((1.0/Gwdiff) - (1.0/Gbound))); //Determine stomatal conductance after accounting for leaf boundary conductance
  NumericVector LP = leafphotosynthesis(Q/refLeafArea, Catm, Gwdiff/1.6, std::max(0.0,leafTemp), Vmax298/refLeafArea, Jmax298/refLeafArea);
  Ci = LP[0];
  Ag = LP[1];
  An = Ag - 0.015*VmaxTemp(Vmax298/refLeafArea, leafTemp);
  return(NumericVector::create(Named("LeafTemperature") = leafTemp,
                               Named("LeafVPD") = leafVPD,
                               Named("Gsw") = Gsw,
                               Named("Ci") = Ci,
                               Named("GrossPhotosynthesis") = Ag,
                               Named("NetPhotosynthesis") = An));
}


//' @rdname photo
//' @keywords internal
// [[Rcpp::export("photo_leafPhotosynthesisFunction2")]]
DataFrame leafPhotosynthesisFunction2(NumericVector E, NumericVector psiLeaf, double Catm, double Patm, double Tair, double vpa, double u, 
                                     double SWRabs, double LWRnet, double Q, double Vmax298, double Jmax298, 
                                     double leafWidth = 1.0, double refLeafArea = 1.0, bool verbose = false) {
  int nsteps = E.size();
  NumericVector leafTemp(nsteps);
  NumericVector leafVPD(nsteps);
  NumericVector Gsw(nsteps), Ci(nsteps);
  NumericVector Ag(nsteps), An(nsteps);
  for(int i=0;i<nsteps;i++){
    NumericVector lpf = leafPhotosynthesisOneFunction2(E[i], psiLeaf[i], Catm, Patm, Tair, vpa, u,
                                                       SWRabs, LWRnet, Q, Vmax298, Jmax298,
                                                       leafWidth, refLeafArea, verbose);
    leafTemp[i] = lpf["LeafTemperature"];
    leafVPD[i] = lpf["LeafVPD"];
    Gsw[i] = lpf["Gsw"];
    Ci[i] = lpf["Ci"];
    Ag[i] = lpf["GrossPhotosynthesis"];
    An[i] = lpf["NetPhotosynthesis"];
  }
  return(DataFrame::create(Named("LeafTemperature") = leafTemp,
                           Named("LeafVPD") = leafVPD,
                           Named("Gsw") = Gsw,
                           Named("Ci") = Ci,
                           Named("GrossPhotosynthesis") = Ag,
                           Named("NetPhotosynthesis") = An));
}


/**
 * Calculates gross/net canopy photosynthesis function, considering a multilayer canopy 
 * and sunlit/shade leaves.
 * (Farquhar et al. 1980/Collatz et al 1991/De Pury and Farquhar)
 * 
 * supplyFunction - Hydraulic supply function
 * 
 * Catm - CO2 air concentration (micromol * mol-1)
 * Patm - Air pressure (kPa)
 * Tair - Air temperature (ºC) - changes through the day and from one day to the other
 * vpa - Air actual vapour pressure (kPa)
 * 
 * SLarea, SHarea - leaf area index of sunlit/shade leaves for each canopy layer
 * u - Wind speed (m/s) for each canopy layer
 * absRadSL, absRadSH - instantaneous absorbed radiation (W·m-2) per unit of sunlit/shade leaf area, for each canopy layer
 * QSL, QSH - Active photon flux density (micromol * s-1 * m-2) per unit of sunlit/shade leaf area, for each canopy layer
 * Vmax298 - maximum Rubisco carboxylation rate per leaf area at 298 ºK (i.e. 25 ºC) (micromol*s-1*m-2), for each canopy layer
 * Jmax298 - maximum electron transport rate per leaf area at 298 ºK (i.e. 25 ºC) (micromol*s-1*m-2), for each canopy layer
 * 
 * return units: micromol*s-1*m-2
 */
//' @rdname photo
//' @keywords internal
// [[Rcpp::export("photo_sunshadePhotosynthesisFunction")]]
DataFrame sunshadePhotosynthesisFunction(NumericVector E, NumericVector psiLeaf, double Catm, double Patm, double Tair, double vpa, 
                                  double SLarea, double SHarea,
                                  double u, double absRadSL, double absRadSH,
                                  double QSL, double QSH, 
                                  double Vmax298SL, double Vmax298SH, 
                                  double Jmax298SL, double Jmax298SH, 
                                  double leafWidth = 1.0, bool verbose = false) {
  int nsteps = E.size();
  NumericVector Ag(nsteps,0.0), An(nsteps,0.0);
  NumericVector leafCiSL(nsteps,0.0), leafCiSH(nsteps,0.0);
  NumericVector leafTSL(nsteps,0.0), leafTSH(nsteps,0.0);
  NumericVector leafVPDSL(nsteps,0.0), leafVPDSH(nsteps,0.0);
  // Rcout<<"ws "<<u<<" tair "<< Tair<< " SLarea "<< SLarea << " SHarea "<< SHarea<< " absRadSL"<< absRadSL<< " absRadSH "<< absRadSH<< " QSL "<<QSL<<" QSH "<<QSH<<"\n";
  double leafT, Gwdiff, Gbound, Agj, Anj;
  for(int i=0;i<nsteps;i++){
    //Sunlit leaves
    Ag[i]=0.0;
    An[i]=0.0;
    //From rad per ground area to rad per leaf area
    leafT = leafTemperature(absRadSL/SLarea, Tair, u, E[i], leafWidth);
    leafTSL[i]= leafT;
    leafVPDSL[i] = std::max(0.0,leafVapourPressure(leafT, psiLeaf[i]) - vpa);
    // Assumes infinite boundary layer conductance (well-coupling) so that stomatal conductance equals diffusive conductance
    // Gw = Patm*(E[i]/1000.0)/leafVPDSL[i];
    // Gw = Gw*SLarea; //From Gw per leaf area to Gw per ground area
    // Separates diffusive conductance into stomatal and boundary layer conductance
    Gwdiff = Patm*(E[i]/1000.0)/leafVPDSL[i]; //Transform flow from mmol to mol
    Gbound = gLeafBoundary(u, leafWidth); // mol (boundary layer conductance)
    Gwdiff = std::min(Gwdiff, Gbound); //Diffusive conductance cannot be lower than boundary layer conductance
    Gwdiff = Gwdiff*SLarea; //From Gwdiff per leaf area to Gwdiff per ground area
    if(QSL>0.0) {
      NumericVector LP = leafphotosynthesis(QSL, Catm, Gwdiff/1.6, leafT, Vmax298SL, Jmax298SL);//Call photosynthesis with aggregated values
      leafCiSL[i] = LP[0];
      Agj = LP[1];
      Anj = Agj - 0.015*VmaxTemp(Vmax298SL, leafT);
      Ag[i]+=Agj;
      An[i]+=Anj;
    }
    //SHADE leaves
    //From rad per ground area to rad per leaf area
    leafT = leafTemperature(absRadSH/SHarea, Tair, u, E[i], leafWidth);
    leafTSH[i]= leafT;
    leafVPDSH[i] = std::max(0.0,leafVapourPressure(leafT, psiLeaf[i]) - vpa);
    // leafVPDSH[i] = std::max(0.0,meteoland::utils_saturationVP(leafT)-vpa);
    // Assumes infinite boundary layer conductance (well-coupling) so that stomatal conductance equals diffusive conductance
    // Gw = Patm*(E[i]/1000.0)/leafVPDSH[i];
    // Gw = Gw*SHarea; //From Gw per leaf area to Gw per ground area
    // Separates diffusive conductance into stomatal and boundary layer conductance
    Gwdiff = Patm*(E[i]/1000.0)/leafVPDSH[i]; //Transform flow from mmol to mol
    Gbound = gLeafBoundary(u, leafWidth); // mol boundary layer conductance
    Gwdiff = std::min(Gwdiff, Gbound); //Diffusive conductance cannot be lower than boundary layer conductance
    // Gsw  = std::abs(1.0/((1.0/Gwdiff) - (1.0/Gbw))); //Determine stomatal conductance after accounting for leaf boundary conductance
    Gwdiff = Gwdiff*SLarea; //From Gwdiff per leaf area to Gwdiff per ground area
    if(QSH>0.0) {
      NumericVector LP = leafphotosynthesis(QSH, Catm, Gwdiff/1.6, leafT, Vmax298SH, Jmax298SH); //Call photosynthesis with aggregated values
      leafCiSH[i] = LP[0];
      Agj = LP[1];
      Anj = Agj - 0.015*VmaxTemp(Vmax298SH, leafT);
      Ag[i]+=Agj;
      An[i]+=Anj;
    }
    
  }
  return(DataFrame::create(Named("GrossPhotosynthesis") = Ag,
                      Named("NetPhotosynthesis") = An,
                      Named("LeafCiSL") = leafCiSL,
                      Named("LeafCiSH") = leafCiSH,
                      Named("LeafTempSL") = leafTSL,
                      Named("LeafTempSH") = leafTSH,
                      Named("LeafVPDSL") = leafVPDSL,
                      Named("LeafVPDSH") = leafVPDSH));
}

//' @rdname photo
//' @keywords internal
// [[Rcpp::export("photo_multilayerPhotosynthesisFunction")]]
DataFrame multilayerPhotosynthesisFunction(NumericVector E, NumericVector psiLeaf, 
                                           double Catm, double Patm, double Tair, double vpa, 
                                  NumericVector SLarea, NumericVector SHarea,
                                  NumericVector u, NumericVector absRadSL, NumericVector absRadSH,
                                  NumericVector QSL, NumericVector QSH, 
                                  NumericVector Vmax298, NumericVector Jmax298, 
                                  double leafWidth = 1.0, bool verbose = false) {
  int nsteps = E.size();
  int nlayers = SLarea.size();
  NumericVector Ag(nsteps,0.0), An(nsteps,0.0);
  double leafT,leafVPD, Gwdiff, Gbound, Agj, Anj;
  for(int i=0;i<nsteps;i++){
    Ag[i]=0.0;
    An[i]=0.0;
    for(int j=0;j<nlayers;j++) {
      //Sunlit leaves
      leafT = leafTemperature(absRadSL[j], Tair, u[j], E[i], leafWidth);
      leafVPD = std::max(0.0,leafVapourPressure(leafT, psiLeaf[i]) - vpa);
      // leafVPD = std::max(0.0,meteoland::utils_saturationVP(leafT)-vpa);
      // Assumes infinite boundary layer conductance (well-coupling) so that stomatal conductance equals diffusive conductance
      // Gw = Patm*(E[i]/1000.0)/leafVPD;
      // Separates diffusive conductance into stomatal and boundary layer conductance
      Gwdiff = Patm*(E[i]/1000.0)/leafVPD; //Transform flow from mmol to mol
      Gbound = gLeafBoundary(u[j], leafWidth); // mol boundary layer conductance
      Gwdiff = std::min(Gwdiff, Gbound); //Diffusive conductance cannot be lower than boundary layer conductance
      if(QSL[j]>0.0) {
        NumericVector LP = leafphotosynthesis(QSL[j], Catm, Gwdiff/1.6, leafT, Vmax298[j], Jmax298[j]);
        Agj = LP[1];
        Anj = Agj - 0.015*VmaxTemp(Vmax298[j], leafT);
        //From A per leaf area to A per ground area
        Ag[i]+=Agj*SLarea[j];
        An[i]+=Anj*SLarea[j];
      }
      //SHADE leaves
      leafT = leafTemperature(absRadSH[j], Tair, u[j], E[i], leafWidth);
      leafVPD = std::max(0.0,leafVapourPressure(leafT, psiLeaf[i]) - vpa);
      // leafVPD = std::max(0.0,meteoland::utils_saturationVP(leafT)-vpa);
      // Assumes infinite boundary layer conductance (well-coupling) so that stomatal conductance equals diffusive conductance
      // Gw = Patm*(E[i]/1000.0)/leafVPD;
      // Separates diffusive conductance into stomatal and boundary layer conductance
      Gwdiff = Patm*(E[i]/1000.0)/leafVPD; //Transform flow from mmol to mol
      Gbound = gLeafBoundary(u[j], leafWidth); // mol boundary layer conductance
      Gwdiff = std::min(Gwdiff, Gbound); //Diffusive conductance cannot be lower than boundary layer conductance
      if(QSH[j]>0.0) {
        NumericVector LP = leafphotosynthesis(QSH[j], Catm, Gwdiff/1.6, leafT, Vmax298[j], Jmax298[j]);
        Agj = LP[1];
        Anj = Agj - 0.015*VmaxTemp(Vmax298[j], leafT);
        Ag[i]+=Agj*SHarea[j];
        An[i]+=Anj*SHarea[j];
      }
    }
  }
  return(DataFrame::create(Named("GrossPhotosynthesis") = Ag,
                      Named("NetPhotosynthesis") = An));
}



 