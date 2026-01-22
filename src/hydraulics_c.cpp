#include "hydraulics_c.h"
#include <cmath>



//' Hydraulic conductance functions
//' 
//' Set of functions used in the calculation of soil and plant hydraulic conductance.
//'
//' @param psi A scalar (or a vector, depending on the function) with water potential (in MPa).
//' @param K Whole-plant relative conductance (0-1).
//' @param psi_extract Soil water potential (in MPa) corresponding to 50% whole-plant relative transpiration.
//' @param exp_extract Exponent of the whole-plant relative transpiration Weibull function.
//' @param v Proportion of fine roots within each soil layer.
//' @param krhizomax Maximum rhizosphere hydraulic conductance (defined as flow per leaf surface unit and per pressure drop).
//' @param kxylemmax Maximum xylem hydraulic conductance (defined as flow per leaf surface unit and per pressure drop).
//' @param c,d Parameters of the Weibull function (generic xylem vulnerability curve).
//' @param n,alpha Parameters of the Van Genuchten function (rhizosphere vulnerability curve).
//' @param kxylem Xylem hydraulic conductance (defined as flow per surface unit and per pressure drop).
//' @param pCrit Proportion of maximum conductance considered critical for hydraulic functioning.
//' @param psi50,psi88,psi12 Water potentials (in MPa) corresponding to 50%, 88% and 12% percent conductance loss.
//' @param temp Temperature (in degrees Celsius).
//' 
//' @details Details of plant hydraulic models are given the medfate book. 
//' Function \code{hydraulics_vulnerabilityCurvePlot} draws a plot of the vulnerability curves for the given \code{soil} object and network properties of each plant cohort in \code{x}.
//' 
//' @return
//' Values returned for each function are:
//' \itemize{
//'   \item{\code{hydraulics_psi2K}: Whole-plant relative conductance (0-1).}
//'   \item{\code{hydraulics_K2Psi}: Soil water potential (in MPa) corresponding to the given whole-plant relative conductance value (inverse of \code{hydraulics_psi2K()}).}
//'   \item{\code{hydraulics_averagePsi}: The average water potential (in MPa) across soil layers.}
//'   \item{\code{hydraulics_vanGenuchtenConductance}: Rhizosphere conductance corresponding to an input water potential (soil vulnerability curve).}
//'   \item{\code{hydraulics_xylemConductance}: Xylem conductance (flow rate per pressure drop) corresponding to an input water potential (plant vulnerability curve).}
//'   \item{\code{hydraulics_xylemPsi}: Xylem water potential (in MPa) corresponding to an input xylem conductance (flow rate per pressure drop).}
//'   \item{\code{hydraulics_psi2Weibull}: Parameters of the Weibull vulnerability curve that goes through the supplied psi50 and psi88 values.}
//' }
//' 
//' @references
//' Sperry, J. S., F. R. Adler, G. S. Campbell, and J. P. Comstock. 1998. Limitation of plant water use by rhizosphere and xylem conductance: results from a model. Plant, Cell and Environment 21:347–359.
//' 
//' Sperry, J. S., and D. M. Love. 2015. What plant hydraulics can tell us about responses to climate-change droughts. New Phytologist 207:14–27.
//' 
//' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
//' 
//' @seealso
//' \code{\link{hydraulics_supplyFunctionPlot}}, \code{\link{hydraulics_maximumStemHydraulicConductance}}, \code{\link{spwb}}, \code{\link{soil}}
//' 
//' @examples
//' 
//' #Manual display of vulnerability curve
//' kstemmax = 4 # in mmol·m-2·s-1·MPa-1
//' stemc = 3 
//' stemd = -4 # in MPa
//' psiVec = seq(-0.1, -7.0, by =-0.01)
//' kstem = unlist(lapply(psiVec, hydraulics_xylemConductance, kstemmax, stemc, stemd))
//' plot(-psiVec, kstem, type="l",ylab="Xylem conductance (mmol·m-2·s-1·MPa-1)", 
//'      xlab="Canopy pressure (-MPa)", lwd=1.5,ylim=c(0,kstemmax))
//' 
//' #Load example dataset
//' data(exampleforest)
//' 
//' #Default species parameterization
//' data(SpParamsMED)
//' 
//' #Initialize soil with default soil params (4 layers)
//' examplesoil <- defaultSoilParams(4)
//' 
//' #Initialize control parameters
//' control <- defaultControl("Granier")
//' 
//' #Switch to 'Sperry' transpiration mode
//' control <- defaultControl("Sperry")
//' 
//' #Initialize input
//' x <- spwbInput(exampleforest,examplesoil, SpParamsMED, control)
//' 
//' #Leaf vulnerability curves
//' hydraulics_vulnerabilityCurvePlot(x, type="leaf")
//' 
//' #Stem vulnerability curves
//' hydraulics_vulnerabilityCurvePlot(x, type="stem")
//'              
//' @name hydraulics_conductancefunctions
//' @keywords internal
// [[Rcpp::export("hydraulics_psi2K")]]
double Psi2K(double psi, double psi_extract, double exp_extract = 3.0) {
  return(exp(-0.6931472*pow(std::abs(psi/psi_extract),exp_extract)));
}

/**
 * Inverse of the whole-plant conductance function. Used to obtain the 'average' soil water
 * potential perceived by each plant cohort.
 */
//' @rdname hydraulics_conductancefunctions
//' @keywords internal
// [[Rcpp::export("hydraulics_K2Psi")]]
double K2Psi(double K, double psi_extract, double exp_extract = 3.0) {
  double psi = psi_extract*pow(log(K)/(-0.6931472),1.0/exp_extract);
  if(psi>0.0) psi = -psi; //Usually psi_extr is a positive number
  if(psi< -40.0) psi = -40.0; //Minimum value
  return(psi);
}

/**
 *  Returns minimum leaf conductance (in mmolH2O·m-2·s-1·MPa-1) as a function of leaf temperature
 *  according to a two-Q10 model.
 */
double gmin(double leafTemperature, double gmin_20, 
            double TPhase, double Q10_1, double Q10_2) {
  double gmin;
  if (leafTemperature<= TPhase) {
    gmin = gmin_20 * pow(Q10_1,(leafTemperature - 20.0) / 10.0);
  } else if (leafTemperature > TPhase) {
    gmin = gmin_20 * pow(Q10_1, (TPhase - 20.0) / 10.0) * pow(Q10_2, (leafTemperature- TPhase) / 10.0);
  }
  return(gmin);
}

//' @rdname hydraulics_conductancefunctions
// [[Rcpp::export("hydraulics_xylemConductance")]]
double xylemConductance(double psi, double kxylemmax, double c, double d) {
  if(psi>=0.0) return(kxylemmax);
  return(kxylemmax*exp(-pow(psi/d,c)));
}



//' @rdname hydraulics_conductancefunctions
// [[Rcpp::export("hydraulics_xylemConductanceSigmoid")]]
double xylemConductanceSigmoid(double psi, double kxylemmax, double P50, double slope){
  if(psi>=0.0) return(kxylemmax);
  return (kxylemmax * (1.0 - 1.0/(1.0 + exp(slope / 25.0 * (psi - P50)))));
}

//' @rdname hydraulics_conductancefunctions
//' @keywords internal
// [[Rcpp::export("hydraulics_xylemPsi")]]
double xylemPsi(double kxylem, double kxylemmax, double c, double d) {
  double psi = d*pow(-log(kxylem/kxylemmax),1.0/c);
  if(psi< -40.0) psi = -40.0; //Minimum value
  return(psi);
}

//' @rdname hydraulics_conductancefunctions
//' @keywords internal
// [[Rcpp::export("hydraulics_psiCrit")]]
double psiCrit(double c, double d, double pCrit = 0.001) {
  return(d * pow(-log(pCrit), 1.0/c));
}

/**
 * Van genuchten-mualem conductance equation (m = 1 - 1/n; l = 0.5)
 */
//' @rdname hydraulics_conductancefunctions
//' @keywords internal
// [[Rcpp::export("hydraulics_vanGenuchtenConductance")]]
double vanGenuchtenConductance(double psi, double krhizomax, double n, double alpha) {
  double v = 1.0/(pow(alpha*std::abs(psi),n)+1.0);
  return(krhizomax*pow(v,(n-1.0)/(2.0*n))*pow(pow((1.0-v),(n-1.0)/n)-1.0,2.0));
}
