#include "hydraulics_c.h"
#include "medfate.h"
#include "incgamma_c.h"
#include "biophysicsutils_c.h"
#include "numerical_solving_c.h"
#include "meteoland/utils_c.hpp"
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
double Psi2K_c(double psi, double psi_extract, double exp_extract = 3.0) {
  return(std::exp(-0.6931472*std::pow(std::abs(psi/psi_extract),exp_extract)));
}

/**
 * Inverse of the whole-plant conductance function. Used to obtain the 'average' soil water
 * potential perceived by each plant cohort.
 */
//' @rdname hydraulics_conductancefunctions
//' @keywords internal
// [[Rcpp::export("hydraulics_K2Psi")]]
double K2Psi_c(double K, double psi_extract, double exp_extract = 3.0) {
  double psi = psi_extract*pow(log(K)/(-0.6931472),1.0/exp_extract);
  if(psi>0.0) psi = -psi; //Usually psi_extr is a positive number
  if(psi< -40.0) psi = -40.0; //Minimum value
  return(psi);
}

/**
 *  Returns minimum leaf conductance (in mmolH2O·m-2·s-1·MPa-1) as a function of leaf temperature
 *  according to a two-Q10 model.
 */
double gmin_c(double leafTemperature, double gmin_20, 
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
double xylemConductance_c(double psi, double kxylemmax, double c, double d) {
  if(psi>=0.0) return(kxylemmax);
  return(kxylemmax*exp(-pow(psi/d,c)));
}



//' @rdname hydraulics_conductancefunctions
// [[Rcpp::export("hydraulics_xylemConductanceSigmoid")]]
double xylemConductanceSigmoid_c(double psi, double kxylemmax, double P50, double slope){
  if(psi>=0.0) return(kxylemmax);
  return (kxylemmax * (1.0 - 1.0/(1.0 + exp(slope / 25.0 * (psi - P50)))));
}

//' @rdname hydraulics_conductancefunctions
//' @keywords internal
// [[Rcpp::export("hydraulics_xylemPsi")]]
double xylemPsi_c(double kxylem, double kxylemmax, double c, double d) {
  double psi = d*pow(-log(kxylem/kxylemmax),1.0/c);
  if(psi< -40.0) psi = -40.0; //Minimum value
  return(psi);
}

//' @rdname hydraulics_conductancefunctions
//' @keywords internal
// [[Rcpp::export("hydraulics_psiCrit")]]
double psiCrit_c(double c, double d, double pCrit = 0.001) {
  return(d * pow(-log(pCrit), 1.0/c));
}

/**
 * Van genuchten-mualem conductance equation (m = 1 - 1/n; l = 0.5)
 */
//' @rdname hydraulics_conductancefunctions
//' @keywords internal
// [[Rcpp::export("hydraulics_vanGenuchtenConductance")]]
double vanGenuchtenConductance_c(double psi, double krhizomax, double n, double alpha) {
  double v = 1.0/(pow(alpha*std::abs(psi),n)+1.0);
  return(krhizomax*pow(v,(n-1.0)/(2.0*n))*pow(pow((1.0-v),(n-1.0)/n)-1.0,2.0));
}


//' @rdname hydraulics_conductancefunctions
// [[Rcpp::export("hydraulics_correctConductanceForViscosity")]]
double correctConductanceForViscosity_c(double kxylem, double temp) {
  return(kxylem/waterDynamicViscosity_c(temp));
}


double averagePsi_c(const std::vector<double>& psi, const std::vector<double>& v, double exp_extract, double psi_extract) {
  int nlayers = psi.size();
  double K = 0.0, sumKv = 0.0;
  for(int l=0;l<nlayers;l++) {
    K= exp(-0.6931472*pow(std::abs(psi[l]/psi_extract),exp_extract));
    sumKv += K*v[l];
  }
  double psires =  psi_extract*pow(log(sumKv)/(-0.6931472),1.0/exp_extract);
  psires = std::max(psires, -40.0); //Limits plant water potential to -40 MPa
  return(psires);
}

double averagePsiPool_c(const arma::mat& Psi, const arma::mat& RHOPcohV, double exp_extract, double psi_extract) {
  int nlayers = Psi.n_cols;
  int numCohorts = Psi.n_rows;
  double K = 0.0, sumKv = 0.0;
  for(int j =0;j<numCohorts;j++) {
    for(int l=0;l<nlayers;l++) {
      K = exp(-0.6931472*pow(std::abs(Psi(j,l)/psi_extract),exp_extract)); 
      sumKv += K*RHOPcohV(j,l);
    } 
  }
  double psires =  psi_extract*pow(log(sumKv)/(-0.6931472),1.0/exp_extract);
  psires = std::max(psires, -40.0); //Limits plant water potential to -40 MPa
  return(psires);
}


// [[Rcpp::export(".Egamma")]]
double Egamma_c(double psi, double kxylemmax, double c, double d, double psiCav) {
  if(psi>0.0) return(-Egamma_c(-psi, kxylemmax,c,d,0.0));
  else if(psi==0.0) return(0.0);
  double h = 1.0/c;
  double z = pow(psi/d,c);
  std::vector<double> pq = incgam_c(h,z);
  double g = tgamma(h)*pq[0]; //Upper incomplete gamma, without the normalizing factor
  double E = kxylemmax*(-d/c)*g;
  if(psiCav<0.0) { //Decrease E from 0 to psiCav (avoid recursiveness!)
    if(psiCav < psi) {
      E = xylemConductance_c(psiCav,kxylemmax,c,d)*(-psi); //square integral
    } else {
      std::vector<double> pq = incgam_c(h,pow(psiCav/d,c));
      double Epsimin = kxylemmax*(-d/c)*tgamma(h)*pq[0];
      E = E - Epsimin + xylemConductance_c(psiCav,kxylemmax,c,d)*(-psiCav); //Remove part of the integral corresponding to psimin and add square integral
    }
  }
  return(E);
}

// [[Rcpp::export(".Egammainv")]]
double Egammainv_c(double Eg, double kxylemmax, double c, double d, double psiCav) {
  if(psiCav<0.0) {
    double Eq = xylemConductance_c(psiCav,kxylemmax,c,d)*(-psiCav);
    if(Eg > Eq) {
      double Ec = Egamma_c(psiCav, kxylemmax, c, d) - Eq;
      Eg = Eg + Ec; 
    } else {
      return(-1.0*(Eg/xylemConductance_c(psiCav,kxylemmax,c,d)));
    }
  }
  double h = 1.0/c;
  double g = (-c/d)*(Eg/kxylemmax);
  double p = g/tgamma(h);
  double q = 1.0 - p;//Upper incomplete gamma, without the normalizing factor
  double x = invincgam_c(h,p,q);
  double psi = d*pow(x, 1.0/c);
  return(psi);
}

/*
 * Integral of the xylem vulnerability curve
 */
//' Hydraulic supply functions
//' 
//' Set of functions used in the implementation of hydraulic supply functions (Sperry and Love 2015).
//'
//' @param v Proportion of fine roots within each soil layer.
//' @param krhizomax Maximum rhizosphere hydraulic conductance (defined as flow per leaf surface unit and per pressure drop).
//' @param kxylemmax Maximum xylem hydraulic conductance (defined as flow per leaf surface unit and per pressure drop).
//' @param kleafmax Maximum leaf hydraulic conductance (defined as flow per leaf surface unit and per pressure drop).
//' @param kstemmax Maximum stem xylem hydraulic conductance (defined as flow per leaf surface unit and per pressure drop).
//' @param E Flow per surface unit.
//' @param Emax Maximum flow per surface unit.
//' @param Erootcrown Flow per surface unit at the root crown.
//' @param psiDownstream Water potential upstream (in MPa).
//' @param psiUpstream Water potential upstream (in MPa). In a one-component model corresponds to soil potential. In a two-component model corresponds to the potential inside the roots.
//' @param psiCav Minimum water potential (in MPa) experienced (for irreversible cavitation).
//' @param minFlow Minimum flow in supply function.
//' @param psiPlant Plant water potential (in MPa).
//' @param hydraulicNetwork List with the hydraulic characteristics of nodes in the hydraulic network.
//' @param psiSoil Soil water potential (in MPa). A scalar or a vector depending on the function.
//' @param psiRhizo Soil water potential (in MPa) in the rhizosphere (root surface).
//' @param psiRootCrown Soil water potential (in MPa) at the root crown.
//' @param psiStep Water potential precision (in MPa).
//' @param psiIni Vector of initial water potential values (in MPa).
//' @param psiMax Minimum (maximum in absolute value) water potential to be considered (in MPa).
//' @param pCrit Critical water potential (in MPa).
//' @param dE Increment of flow per surface unit.
//' @param c,d Parameters of the Weibull function (generic xylem vulnerability curve).
//' @param stemc,stemd Parameters of the Weibull function for stems (stem xylem vulnerability curve).
//' @param leafc,leafd Parameters of the Weibull function for leaves (leaf vulnerability curve).
//' @param n,alpha,l Parameters of the Van Genuchten function (rhizosphere vulnerability curve).
//' @param allowNegativeFlux A boolean to indicate whether negative flux (i.e. from plant to soil) is allowed.
//' @param maxNsteps Maximum number of steps in the construction of supply functions.
//' 
//' @details 
//' Function \code{hydraulics_supplyFunctionPlot} draws a plot of the supply function for the given \code{soil} object and network properties of each plant cohort in \code{x}. Function \code{hydraulics_vulnerabilityCurvePlot} draws a plot of the vulnerability curves for the given \code{soil} object and network properties of each plant cohort in \code{x}.
//' 
//' @return
//' Values returned for each function are:
//' \itemize{
//'   \item{\code{hydraulics_E2psiXylem}: The plant (leaf) water potential (in MPa) corresponding to the input flow, according to the xylem supply function and given an upstream (soil or root) water potential.}
//'   \item{\code{hydraulics_E2psiVanGenuchten}: The root water potential (in MPa) corresponding to the input flow, according to the rhizosphere supply function and given a soil water potential.}
//'   \item{\code{hydraulics_E2psiTwoElements}: The plant (leaf) water potential (in MPa) corresponding to the input flow, according to the rhizosphere and plant supply functions and given an input soil water potential.}
//'   \item{\code{hydraulics_E2psiNetwork}: The rhizosphere, root crown and plant (leaf water potential (in MPa) corresponding to the input flow, according to the vulnerability curves of rhizosphere, root and stem elements in a network.}
//'   \item{\code{hydraulics_Ecrit}: The critical flow according to the xylem supply function and given an input soil water potential.}
//'   \item{\code{hydraulics_EVanGenuchten}: The flow (integral of the vulnerability curve) according to the rhizosphere supply function and given an input drop in water potential (soil and rhizosphere).}
//'   \item{\code{hydraulics_EXylem}: The flow (integral of the vulnerability curve) according to the xylem supply function and given an input drop in water potential (rhizosphere and plant).}
//'   \item{\code{hydraulics_supplyFunctionOneXylem}, \code{hydraulics_supplyFunctionTwoElements} and
//'     \code{hydraulics_supplyFunctionNetwork}: A list with different numeric vectors with information of the two-element supply function:
//'     \itemize{
//'       \item{\code{E}: Flow values (supply values).}
//'       \item{\code{FittedE}: Fitted flow values (for \code{hydraulics_supplyFunctionTwoElements}).}
//'       \item{\code{Elayers}: Flow values across the roots of each soil layer (only for \code{hydraulics_supplyFunctionNetwork}).}
//'       \item{\code{PsiRhizo}: Water potential values at the root surface (only for \code{hydraulics_supplyFunctionNetwork}).}
//'       \item{\code{PsiRoot}: Water potential values inside the root crown (not for \code{hydraulics_supplyFunctionOneXylem}).}
//'       \item{\code{PsiPlant}: Water potential values at the canopy (leaf).}
//'       \item{\code{dEdP}: Derivatives of the supply function.}
//'     }
//'   }
//'   \item{\code{hydraulics_supplyFunctionPlot}: If \code{draw = FALSE} a list with the result of calling \code{hydraulics_supplyFunctionNetwork} for each cohort. }
//'   \item{\code{hydraulics_regulatedPsiXylem}: Plant water potential after regulation (one-element loss function) given an input water potential.}
//'   \item{\code{hydraulics_regulatedPsiTwoElements}: Plant water potential after regulation (two-element loss function) given an input soil water potential.}
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
//' \code{\link{hydraulics_psi2K}}, \code{\link{hydraulics_maximumStemHydraulicConductance}}, \code{\link{spwb}}, \code{\link{soil}}
//' 
//' @examples
//' kstemmax = 4 # in mmol·m-2·s-1·MPa-1
//' stemc = 3 
//' stemd = -4 # in MPa
//' psiVec = seq(-0.1, -7.0, by =-0.01)
//' 
//' #Vulnerability curve
//' kstem = unlist(lapply(psiVec, hydraulics_xylemConductance, kstemmax, stemc, stemd))
//' plot(-psiVec, kstem, type="l",ylab="Xylem conductance (mmol·m-2·s-1·MPa-1)", 
//'      xlab="Canopy pressure (-MPa)", lwd=1.5,ylim=c(0,kstemmax))
//' 
//' @name hydraulics_supplyfunctions
//' @keywords internal
// [[Rcpp::export("hydraulics_EXylem")]]
double EXylem_c(double psiPlant, double psiUpstream, 
                double kxylemmax, double c, double d, 
                bool allowNegativeFlux, double psiCav) {
  if((psiPlant > psiUpstream) && !allowNegativeFlux) throw std::range_error("Downstream potential larger (less negative) than upstream potential");
  return(Egamma_c(psiPlant, kxylemmax, c, d, psiCav)-Egamma_c(psiUpstream, kxylemmax, c,d, psiCav));
}

//' @rdname hydraulics_supplyfunctions
// [[Rcpp::export("hydraulics_E2psiXylem")]]
double E2psiXylem_c(double E, double psiUpstream, 
                    double kxylemmax, double c, double d, double psiCav) {
  if(E==0) return(psiUpstream);
  double Eg = E + Egamma_c(psiUpstream, kxylemmax, c,d, psiCav);
  return(Egammainv_c(Eg, kxylemmax, c, d, psiCav));
}


//' @rdname hydraulics_supplyfunctions
//' @keywords internal
// [[Rcpp::export("hydraulics_E2psiXylemUp")]]
double E2psiXylemUp_c(double E, double psiDownstream, 
                      double kxylemmax, double c, double d, double psiCav) {
  if(E==0) return(psiDownstream);
  double Eg = Egamma_c(psiDownstream, kxylemmax, c,d, psiCav) - E;
  return(Egammainv_c(Eg, kxylemmax, c, d, psiCav));
}


/**
 * Analytical approximation to the integral of van genuchten model
 * Van Lier QDJ, Neto DD, Metselaar K (2009) Modeling of transpiration reduction in van genuchten-mualem type soils. 
 * Water Resour Res 45:1–9. doi: 10.1029/2008WR006938
 */
//' @rdname hydraulics_supplyfunctions
// [[Rcpp::export("hydraulics_EVanGenuchten")]]
double EVanGenuchten_c(double psiRhizo, double psiSoil, double krhizomax, 
                       double n, double alpha, double l) {
  double m = 1.0 - (1.0/n);
  double thetaR = pow(1.0+pow(alpha*std::abs(psiRhizo),n),-1.0*m);
  double thetaS = pow(1.0+pow(alpha*std::abs(psiSoil),n),-1.0*m);
  double a1 = (1.0/m) + l + 1.0;
  double a2 = (2.0/m) + l + 1.0;
  double a3 = (3.0/m) + l + 1.0;
  double phi = m*(l+1.0);
  double B1 = ((1.0+phi)*(2.0 + m))/(3.0*(2.0+phi));
  double B2 = ((2.0+phi)*(3.0 + m))/(4.0*(3.0+phi));
  double B3 = ((1.0+phi)*(2.0 - m))/(3.0*(2.0+phi));
  double B4 = ((2.0+phi)*(3.0 - m))/(4.0*(3.0+phi));
  
  double GammaR1 = (2.0*m*pow(thetaR, a1));
  double GammaS1 = (2.0*m*pow(thetaS, a1));
  double Gamma2 = ((1.0 + m)*B1 - (1.0-m)*B3);
  double Gamma3 = ((1.0 + m)*B1*B2 - (1.0-m)*B3*B4);
  double GammaR = GammaR1 + Gamma2*pow(thetaR, a2) + Gamma3*pow(thetaR, a3);
  double GammaS = GammaS1 + Gamma2*pow(thetaS, a2) + Gamma3*pow(thetaS, a3);
  double E = ((m*(1.0-m)*krhizomax)/(2.0*alpha*(phi+1)))*(GammaR - GammaS);
  return(-E);
}


//' @rdname hydraulics_supplyfunctions
//' @keywords internal
// [[Rcpp::export("hydraulics_ECrit")]]
double ECrit_c(double psiUpstream, double kxylemmax, double c, double d, double pCrit) {
  return(EXylem_c(psiCrit_c(c,d, pCrit), psiUpstream, kxylemmax, c, d));
}


//' @rdname hydraulics_supplyfunctions
//' @keywords internal
// [[Rcpp::export("hydraulics_E2psiVanGenuchten")]]
double E2psiVanGenuchten_c(double E, double psiSoil, double krhizomax, double n, double alpha, 
                           double psiStep, double psiMax) {
  if(E<0.0) throw medfate::MedfateInternalError("E has to be positive");
  if(E==0) return(psiSoil);
  double psi = psiSoil;
  double psiPrev = psi;
  double vgPrev = vanGenuchtenConductance_c(psi, krhizomax, n, alpha);
  double vg = vgPrev;
  double Eg = 0.0;
  while(Eg<E) {
    psiPrev = psi;
    vgPrev = vg;
    psi = psi + psiStep;
    vg = vanGenuchtenConductance_c(psi, krhizomax, n, alpha);
    Eg = Eg + ((vg+vgPrev)/2.0)*std::abs(psiStep);
    if(psi<psiMax) return(NA_REAL);
  }
  return(psiPrev);
}

//' @rdname hydraulics_supplyfunctions
//' @keywords internal
// [[Rcpp::export("hydraulics_E2psiTwoElements")]]
double E2psiTwoElements_c(double E, double psiSoil, double krhizomax, double kxylemmax, double n, double alpha, double c, double d, double psiCav,
                          double psiStep, double psiMax) {
  if(E<0.0) throw medfate::MedfateInternalError("E has to be positive");
  if(E==0) return(psiSoil);
  double psiRoot = E2psiVanGenuchten_c(E, psiSoil, krhizomax, n, alpha, psiStep, psiMax);
  if(std::isnan(psiRoot)) return(medfate::NA_DOUBLE);
  return(E2psiXylem_c(E, psiRoot, kxylemmax, c, d, psiCav));
}


//' Hydraulic-related defoliation
//' 
//' Functions to calculate the proportion of crown defoliation due to hydraulic disconnection.
//'
//' @param psiLeaf Leaf water potential (in MPa).
//' @param c,d Parameters of the Weibull function.
//' @param P50,slope Parameters of the Sigmoid function.
//' @param PLC_crit Critical leaf PLC corresponding to defoliation
//' @param P50_cv Coefficient of variation (in percent) of leaf P50, to describe the
//' variability in hydraulic vulnerability across crown leaves.
//' 
//' @details The functions assume that crowns are made of a population of leaves whose
//' hydraulic vulnerability (i.e. the water potential corresponding to 50% loss of conductance) 
//' follows a Gaussian distribution centered on the input P50 and with a known coefficient of variation (\code{P50_cv}).
//' The slope parameter (or the c exponent in the case of a Weibull function) is considered constant.
//' Leaves are hydraulically disconnected, and shedded, when their embolism rate exceeds a critical value (\code{PLC_crit}).
//' 
//' @return The proportion of crown defoliation.
//' 
//' @author 
//' Hervé Cochard, INRAE 
//' 
//' Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
//' 
//' @seealso
//' \code{\link{hydraulics_conductancefunctions}}
//' 
//' @name hydraulics_defoliation
//' @keywords internal
// [[Rcpp::export("hydraulics_proportionDefoliationSigmoid")]]
double proportionDefoliationSigmoid_c(double psiLeaf, double P50, double slope, 
                                     double PLC_crit = 0.88, double P50_cv = 10.0) {
   double P50_crit = psiLeaf - (log((1.0 - PLC_crit)/PLC_crit)/(slope/25));
   double pnorm =  normal_cdf(P50_crit, P50, std::abs((P50_cv/100.0)*P50));
   double PDEF = (1.0-pnorm);
   return(PDEF);
}

//' @name hydraulics_defoliation
//' @keywords internal
// [[Rcpp::export("hydraulics_proportionDefoliationWeibull")]]
double proportionDefoliationWeibull_c(double psiLeaf, double c, double d, 
                                       double PLC_crit = 0.88, double P50_cv = 10.0) {
   double d_crit = psiLeaf/pow(-1.0*log(1.0 - PLC_crit), 1.0/c);
   double P50 = xylemPsi_c(0.5,1.0, c, d);
   double P50_crit = xylemPsi_c(0.5,1.0, c, d_crit);
   double pnorm =  normal_cdf(P50_crit, P50, std::abs((P50_cv/100.0)*P50));
   double PDEF = (1.0-pnorm);
   return(PDEF);
 }



/*
 * Parametrization of rhizosphere conductance
 */
double rhizosphereResistancePercent_c(double psiSoil, 
                                    double krhizomax, double n, double alpha,
                                    double krootmax, double rootc, double rootd,
                                    double kstemmax, double stemc, double stemd,
                                    double kleafmax, double leafc, double leafd) {
  double krhizo = vanGenuchtenConductance_c(psiSoil, krhizomax, n, alpha);
  double kroot = xylemConductance_c(psiSoil, krootmax, rootc, rootd);
  double kstem = xylemConductance_c(psiSoil, kstemmax, stemc, stemd);
  double kleaf = xylemConductance_c(psiSoil, kleafmax, leafc, leafd);
  return(100.0*(1.0/krhizo)/((1.0/kroot)+(1.0/kstem)+(1.0/kleaf)+(1.0/krhizo)));
}


//' @rdname hydraulics_scalingconductance
//' @keywords internal
// [[Rcpp::export("hydraulics_averageRhizosphereResistancePercent")]]
double averageRhizosphereResistancePercent_c(double krhizomax, double n, double alpha,
                                             double krootmax, double rootc, double rootd,
                                             double kstemmax, double stemc, double stemd, 
                                             double kleafmax, double leafc, double leafd,
                                             double psiStep = -0.01){
  double psiC = psiCrit_c(stemc, stemd, 0.001);
  double cnt = 0.0;
  double sum = 0.0;
  for(double psi=0.0; psi>psiC;psi += psiStep) {
    sum +=rhizosphereResistancePercent_c(psi, krhizomax, n,alpha,krootmax, rootc, rootd,
                                         kstemmax, stemc,stemd,
                                         kleafmax, leafc,leafd);
    cnt+=1.0;
  }
  return(sum/cnt);
}


//' @rdname hydraulics_scalingconductance
//' @keywords internal
// [[Rcpp::export("hydraulics_findRhizosphereMaximumConductance")]]
double findRhizosphereMaximumConductance_c(double averageResistancePercent, double n, double alpha,
                                         double krootmax, double rootc, double rootd,
                                         double kstemmax, double stemc, double stemd,
                                         double kleafmax, double leafc, double leafd,
                                         double initialValue) {
  double step = 1.0;
  double fTol = 0.1;
  
  // Rcout<<exp(initialValue)<<"\n";
  double krhizomaxlog = initialValue;
  int nsteps = 0;
  int max_nsteps = 100;
  double f = averageRhizosphereResistancePercent_c(exp(krhizomaxlog), n,alpha,krootmax, rootc, rootd,
                                                 kstemmax, stemc,stemd,
                                                 kleafmax, leafc,leafd);
  while((std::abs(f-averageResistancePercent)>fTol) && (nsteps < max_nsteps)) {
    // Rcout<< nsteps<< " "<<exp(krhizomaxlog) << " "<< f << " "<< averageResistancePercent<< " "<< step<<"\n";
    if(f>averageResistancePercent) {
      if(step < 0) step = -step/2.0;
    } else {
      if(step > 0) step = -step/2.0;
    }
    krhizomaxlog += step; 
    f = averageRhizosphereResistancePercent_c(exp(krhizomaxlog), n,alpha,krootmax, rootc, rootd,
                                            kstemmax, stemc,stemd,
                                            kleafmax, leafc,leafd);
    nsteps++;
  }
  // Rcout<< nsteps<< " "<<exp(krhizomaxlog) << " "<< f << " "<< averageResistancePercent<< " "<< step<<"\n";
  return(exp(krhizomaxlog));
}
