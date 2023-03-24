#include "Rcpp.h"
#include "root.h"
#include "biophysicsutils.h"
#include "tissuemoisture.h"
#include "incgamma.h"
#include <math.h>


using namespace Rcpp;
using namespace std;

double const maxPsi = -0.000001;
double const cmhead2MPa = 0.00009804139; //Constant to transform cm head to MPa

//' Hydraulic confuctance functions
//' 
//' Set of functions used in the calculation of soil and plant hydraulic conductance.
//'
//' @param psi A scalar (or a vector, depending on the function) with water potential (in MPa).
//' @param K Whole-plant relative conductance (0-1).
//' @param psi_extract Soil water potential (in MPa) corresponding to 50\% whole-plant relative transpiration.
//' @param exp_extract Exponent of the whole-plant relative transpiration Weibull function.
//' @param v Proportion of fine roots within each soil layer.
//' @param krhizomax Maximum rhizosphere hydraulic conductance (defined as flow per leaf surface unit and per pressure drop).
//' @param kxylemmax Maximum xylem hydraulic conductance (defined as flow per leaf surface unit and per pressure drop).
//' @param c,d Parameters of the Weibull function (generic xylem vulnerability curve).
//' @param n,alpha Parameters of the Van Genuchten function (rhizosphere vulnerability curve).
//' @param kxylem Xylem hydraulic conductance (defined as flow per surface unit and per pressure drop).
//' @param pCrit Proportion of maximum conductance considered critical for hydraulic functioning.
//' @param psi50,psi88,psi12 Water potentials (in MPa) corresponding to 50\%, 88\% and 12\% percent conductance loss.
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
//' data(exampleforestMED)
//' 
//' #Default species parameterization
//' data(SpParamsMED)
//' 
//' #Initialize soil with default soil params (2 layers)
//' examplesoil = soil(defaultSoilParams(2)) 
//' 
//' #Initialize control parameters
//' control = defaultControl("Granier")
//' 
//' #Switch to 'Sperry' transpiration mode
//' control = defaultControl("Sperry")
//' 
//' #Initialize input
//' x = forest2spwbInput(exampleforestMED,examplesoil, SpParamsMED, control)
//' 
//' #Leaf vulnerability curves
//' hydraulics_vulnerabilityCurvePlot(x, type="leaf")
//' 
//' #Stem vulnerability curves
//' hydraulics_vulnerabilityCurvePlot(x, type="stem")
//'              
//' @name hydraulics_conductancefunctions
// [[Rcpp::export("hydraulics_psi2K")]]
double Psi2K(double psi, double psi_extract, double exp_extract = 3.0) {
  return(exp(-0.6931472*pow(std::abs(psi/psi_extract),exp_extract)));
}
NumericVector Psi2K(double psi, NumericVector psi_extract, double exp_extract = 3.0) {
  int n = psi_extract.size();
  NumericVector k(n);
  for(int i=0; i<n; i++) {
    k[i] = Psi2K(psi,psi_extract[i],exp_extract);
  }
  return k;
}

/**
 * Inverse of the whole-plant conductance function. Used to obtain the 'average' soil water
 * potential perceived by each plant cohort.
 */
//' @rdname hydraulics_conductancefunctions
// [[Rcpp::export("hydraulics_K2Psi")]]
double K2Psi(double K, double psi_extract, double exp_extract = 3.0) {
  double psi = psi_extract*pow(log(K)/(-0.6931472),1.0/exp_extract);
  if(psi>0.0) psi = -psi; //Usually psi_extr is a positive number
  return psi;
}
NumericVector K2Psi(NumericVector K, NumericVector psi_extract, double exp_extract = 3.0) {
  int n = psi_extract.size();
  NumericVector psi(n);
  for(int i=0; i<n; i++) {
    psi[i] = K2Psi(K[i], psi_extract[i], exp_extract);
  }
  return psi;
}

//' @rdname hydraulics_conductancefunctions
// [[Rcpp::export("hydraulics_averagePsi")]]
double averagePsi(NumericVector psi, NumericVector v, double exp_extract, double psi_extract) {
  int nlayers = psi.size();
  NumericVector K(nlayers);
  for(int l=0;l<nlayers;l++) K[l]= exp(-0.6931472*pow(std::abs(psi[l]/psi_extract),exp_extract));
  double psires =  psi_extract*pow(log(sum(K*v))/(-0.6931472),1.0/exp_extract);
  psires = std::max(psires, -40.0); //Limits plant water potential to -40 MPa
  return(psires);
}
double averagePsiPool(NumericMatrix Psi, NumericMatrix RHOPcohV, double exp_extract, double psi_extract) {
  int nlayers = Psi.ncol();
  int numCohorts = Psi.nrow();
  NumericMatrix K(numCohorts, nlayers);
  for(int j =0;j<numCohorts;j++) for(int l=0;l<nlayers;l++) K(j,l)= exp(-0.6931472*pow(std::abs(Psi(j,l)/psi_extract),exp_extract));
  double psires =  psi_extract*pow(log(sum(K*RHOPcohV))/(-0.6931472),1.0/exp_extract);
  psires = std::max(psires, -40.0); //Limits plant water potential to -40 MPa
  return(psires);
}

//' @rdname hydraulics_conductancefunctions
// [[Rcpp::export("hydraulics_xylemConductance")]]
double xylemConductance(double psi, double kxylemmax, double c, double d) {
  if(psi>=0.0) return(kxylemmax);
  return(kxylemmax*exp(-pow(psi/d,c)));
}

//' @rdname hydraulics_conductancefunctions
// [[Rcpp::export("hydraulics_xylemPsi")]]
double xylemPsi(double kxylem, double kxylemmax, double c, double d) {
  return(d*pow(-log(kxylem/kxylemmax),1.0/c));
}

//' @rdname hydraulics_conductancefunctions
// [[Rcpp::export("hydraulics_psiCrit")]]
double psiCrit(double c, double d, double pCrit = 0.001) {
  return(d * pow(-log(pCrit), 1.0/c));
}

/**
 * Van genuchten-mualem conductance equation (m = 1 - 1/n; l = 0.5)
 */
//' @rdname hydraulics_conductancefunctions
// [[Rcpp::export("hydraulics_vanGenuchtenConductance")]]
double vanGenuchtenConductance(double psi, double krhizomax, double n, double alpha) {
  double v = 1.0/(pow(alpha*std::abs(psi),n)+1.0);
  return(krhizomax*pow(v,(n-1.0)/(2.0*n))*pow(pow((1.0-v),(n-1.0)/n)-1.0,2.0));
}

//' @rdname hydraulics_conductancefunctions
// [[Rcpp::export("hydraulics_correctConductanceForViscosity")]]
double correctConductanceForViscosity(double kxylem, double temp) {
  return(kxylem/waterDynamicViscosity(temp));
}

//' @rdname hydraulics_conductancefunctions
// [[Rcpp::export("hydraulics_psi2Weibull")]]
NumericVector psi2Weibull(double psi50, double psi88 = NA_REAL, double psi12 = NA_REAL) {
  if(NumericVector::is_na(psi88) && NumericVector::is_na(psi12)) stop("Either 'psi88' or 'psi12' has to be non-missing");
  double a, psiRatio;
  if(!NumericVector::is_na(psi88)) {
    psiRatio = psi50/psi88;
    a = 0.3269156; // log(0.5)/log(0.12);
  } else {
    psiRatio = psi50/psi12;
    a = 5.422271; // log(0.5)/log(0.88);
  }
  double c = log(a)/log(psiRatio);
  double d = psi50/pow(0.6931472,1.0/c);
  NumericVector par = NumericVector::create(c,d);
  par.attr("names") = CharacterVector::create("c","d");
  return(par);
}


// [[Rcpp::export(".Egamma")]]
double Egamma(double psi, double kxylemmax, double c, double d, double psiCav = 0.0) {
  if(psi>0.0) return(-Egamma(-psi, kxylemmax,c,d,0.0));
  else if(psi==0.0) return(0.0);
  double h = 1.0/c;
  double z = pow(psi/d,c);
  NumericVector pq = incgam(h,z);
  double g = tgamma(h)*pq[0]; //Upper incomplete gamma, without the normalizing factor
  double E = kxylemmax*(-d/c)*g;
  if(psiCav<0.0) { //Decrease E from 0 to psiCav (avoid recursiveness!)
    if(psiCav < psi) {
      E = xylemConductance(psiCav,kxylemmax,c,d)*(-psi); //square integral
    } else {
      NumericVector pq = incgam(h,pow(psiCav/d,c));
      double Epsimin = kxylemmax*(-d/c)*tgamma(h)*pq[0];
      E = E - Epsimin + xylemConductance(psiCav,kxylemmax,c,d)*(-psiCav); //Remove part of the integral corresponding to psimin and add square integral
    }
  }
  return(E);
}

// [[Rcpp::export(".Egammainv")]]
double Egammainv(double Eg, double kxylemmax, double c, double d, double psiCav = 0.0) {
  if(psiCav<0.0) {
    double Eq = xylemConductance(psiCav,kxylemmax,c,d)*(-psiCav);
    if(Eg > Eq) {
      double Ec = Egamma(psiCav, kxylemmax, c, d) - Eq;
      Eg = Eg + Ec; 
    } else {
      return(-1.0*(Eg/xylemConductance(psiCav,kxylemmax,c,d)));
    }
  }
  double h = 1.0/c;
  double g = (-c/d)*(Eg/kxylemmax);
  double p = g/tgamma(h);
  double q = 1.0 - p;//Upper incomplete gamma, without the normalizing factor
  double x = invincgam(h,p,q);
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
//' @param psi Water potential (in MPa).
//' @param psiPrev Water potential (in MPa) in the previous time step.
//' @param psiDownstream Water potential upstream (in MPa).
//' @param psiUpstream Water potential upstream (in MPa). In a one-component model corresponds to soil potential. In a two-component model corresponds to the potential inside the roots.
//' @param psiCav Minimum water potential (in MPa) experienced (for irreversible cavitation).
//' @param minFlow Minimum flow in supply function.
//' @param psiPlant Plant water potential (in MPa).
//' @param hydraulicNetwork List with the hydraulic characteristics of nodes in the hydraulic network.
//' @param psiFineRoot Water potential (in MPa) inside fine roots.
//' @param psiSoil Soil water potential (in MPa). A scalar or a vector depending on the function.
//' @param psiRhizo Soil water potential (in MPa) in the rhizosphere (root surface).
//' @param psiRootCrown Soil water potential (in MPa) at the root crown.
//' @param psiStep Water potential precision (in MPa).
//' @param psiTol Precision for water potential estimates (in MPa).
//' @param psiIni Vector of initial water potential values (in MPa).
//' @param psiMax Minimum (maximum in absolute value) water potential to be considered (in MPa).
//' @param pCrit Critical water potential (in MPa).
//' @param PLCprev Previous proportion of loss conductance [0-1].
//' @param V Capacity of the compartment per leaf area (in L/m2).
//' @param fapo Apoplastic fraction (proportion) in the segment.
//' @param pi0 Full turgor osmotic potential (MPa).
//' @param eps Bulk modulus of elasticity (MPa).
//' @param dE Increment of flow per surface unit.
//' @param ETol Precision for water flow per surface unit.
//' @param c,d Parameters of the Weibull function (generic xylem vulnerability curve).
//' @param stemc,stemd Parameters of the Weibull function for stems (stem xylem vulnerability curve).
//' @param leafc,leafd Parameters of the Weibull function for leaves (leaf vulnerability curve).
//' @param n,alpha,l Parameters of the Van Genuchten function (rhizosphere vulnerability curve).
//' @param allowNegativeFlux A boolean to indicate wether negative flux (i.e. from plant to soil) is allowed.
//' @param maxNsteps Maximum number of steps in the construction of supply functions.
//' @param ntrial Maximum number of steps in Newton-Raphson optimization.
//' @param timestep Time step in seconds.
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
// [[Rcpp::export("hydraulics_EXylem")]]
double EXylem(double psiPlant, double psiUpstream, 
              double kxylemmax, double c, double d, 
              bool allowNegativeFlux = true, double psiCav = 0.0) {
  if((psiPlant > psiUpstream) && !allowNegativeFlux) throw std::range_error("Downstream potential larger (less negative) than upstream potential");
  return(Egamma(psiPlant, kxylemmax, c, d, psiCav)-Egamma(psiUpstream, kxylemmax, c,d, psiCav));
}

//' @rdname hydraulics_supplyfunctions
// [[Rcpp::export("hydraulics_E2psiXylem")]]
double E2psiXylem(double E, double psiUpstream, 
                  double kxylemmax, double c, double d, double psiCav = 0.0) {
  if(E==0) return(psiUpstream);
  double Eg = E + Egamma(psiUpstream, kxylemmax, c,d, psiCav);
  return(Egammainv(Eg, kxylemmax, c, d, psiCav));
}

//' @rdname hydraulics_supplyfunctions
// [[Rcpp::export("hydraulics_E2psiXylemUp")]]
double E2psiXylemUp(double E, double psiDownstream, 
                  double kxylemmax, double c, double d, double psiCav = 0.0) {
  if(E==0) return(psiDownstream);
  double Eg = Egamma(psiDownstream, kxylemmax, c,d, psiCav) - E;
  return(Egammainv(Eg, kxylemmax, c, d, psiCav));
}

/**
 * Analytical approximation to the integral of van genuchten model
 * Van Lier QDJ, Neto DD, Metselaar K (2009) Modeling of transpiration reduction in van genuchten-mualem type soils. 
 * Water Resour Res 45:1–9. doi: 10.1029/2008WR006938
 */
//' @rdname hydraulics_supplyfunctions
// [[Rcpp::export("hydraulics_EVanGenuchten")]]
double EVanGenuchten(double psiRhizo, double psiSoil, double krhizomax, 
                     double n, double alpha, double l = 0.5) {
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
// Numerical integral
// double EVanGenuchten(double psiRhizo, double psiSoil, double krhizomax, double n, double alpha, double psiStep = -0.001, double psiTol = 0.0001, bool allowNegativeFlux = true) {
//   if((psiRhizo>psiSoil) && !allowNegativeFlux) ::Rf_error("Downstream potential larger (less negative) than upstream potential");
//   bool reverse = false;
//   if(psiRhizo>psiSoil) reverse = true;
//   if(reverse) {
//     double tmp = psiSoil;
//     psiSoil = psiRhizo;
//     psiRhizo = tmp;
//   }
//   double psi = psiSoil;
//   double vg = vanGenuchtenConductance(psi, krhizomax, n, alpha);
//   double E = 0.0, vgPrev = vg;
//   psiStep = std::max(psiStep, (psiRhizo-psiSoil)/10.0); //Check that step is not too large
//   do {
//     psi = psi + psiStep;
//     if(psi>psiRhizo) {
//       vgPrev = vg;
//       vg = vanGenuchtenConductance(psi, krhizomax, n, alpha);
//       E += ((vg+vgPrev)/2.0)*std::abs(psiStep);
//     } else {
//       psi = psi - psiStep; //retrocedeix
//       psiStep = psiStep/2.0; //canvia pas
//     }
//   } while (std::abs(psi-psiRhizo)>psiTol);
//   if(reverse) E = -E;
//   return(E);
// }

//' @rdname hydraulics_supplyfunctions
// [[Rcpp::export("hydraulics_ECrit")]]
double ECrit(double psiUpstream, double kxylemmax, double c, double d, double pCrit = 0.001) {
  return(EXylem(psiCrit(c,d, pCrit), psiUpstream, kxylemmax, c, d));
}



/*
 * Calculates drop in water potential along a segment with a given flow
 */
// double E2psiXylem(double E, double psiUpstream, double kxylemmax, double c, double d, double psiCav = 0.0, 
//                    double psiStep = -0.0001, double psiMax = -10.0) {
//   // if(E<0.0) stop("E has to be positive");
//   if(E==0) return(psiUpstream);
//   double psi = psiUpstream;
//   double k = xylemConductance(std::min(psi, psiCav),kxylemmax, c, d);
//   double Eg = 0.0;
//   double psiPrev = psi;
//   double kprev = k;
//   if(E<0.0) psiStep = -psiStep;
//   while(std::abs(Eg)<std::abs(E)) {
//     psiPrev = psi;
//     kprev = k;
//     psi = psi + psiStep;
//     k = xylemConductance(std::min(psi, psiCav),kxylemmax, c, d);
//     Eg = Eg + (-1.0*psiStep)*((kprev+k)/2.0);
//     if(psi<psiMax) return(NA_REAL);
//     if(NumericVector::is_na(Eg)) return(NA_REAL);
//   }
//   return(std::min(0.0,(psiPrev+psi)/2.0));
// }

//' @rdname hydraulics_supplyfunctions
// [[Rcpp::export("hydraulics_ECapacitance")]]
double ECapacitance(double psi, double psiPrev, double PLCprev,
                    double V, double fapo, double c, double d, 
                    double pi0, double eps,
                    double timestep) {
  double m3tommol = 55555556.0;
  double RWCprev = tissueRelativeWaterContent(psiPrev, pi0, eps,
                                              psiPrev, c, d,
                                              fapo);
  double RWC = tissueRelativeWaterContent(psi, pi0, eps,
                                          psi, c, d,
                                          fapo);
  return(((m3tommol*V)/timestep)*(RWCprev-RWC));
}


// double E2psiXylemCap(double E, double psiUpstream,  
//                      double psiPrev, double PLCprev,
//                      double kxylemmax, double c, double d, 
//                      double V, double fapo, double pi0, double eps,
//                      double timestep,
//                      double psiStep = -0.0001, double ETol = 0.0001) {
//   double psiPLC = apoplasticWaterPotential(1.0-PLCprev, c, d);
//   //Matrix flux potential upstream
//   double Egup  = Egamma(psiUpstream, kxylemmax, c,d, psiPLC);
//   double Cefprev = (V/timestep)*tissueRelativeWaterContent(psiPrev, pi0, eps,
//                                                 psiPrev, c, d,
//                                                fapo);
//   double Ect = E+Egup-Cefprev;
//   Rcout<<Egup<<" "<<Cefprev<<" "<< Ect<<"\n";
//   
//   double psi = psiUpstream;
//   double Cef = (V/timestep)*tissueRelativeWaterContent(psi, pi0, eps, psi, c, d,fapo);
//   if((Ect+Cef)<0.0) psiStep = -psiStep;
//   double Egd = Egamma(psi, kxylemmax, c,d, psiPLC); 
//   // Rcout<<Cef<<" "<<Egd<<"\n";
//   // for(int i =0;i<10000;i++) {
//   while(std::abs(Ect-Egd+Cef)>ETol) {
//      if(((Ect-Egd+Cef)<0.0) && (psiStep<0.0)) psiStep = -0.5*psiStep;
//      else if(((Ect-Egd+Cef)>0.0) && (psiStep>0.0)) psiStep = -0.5*psiStep;
//      
//      // Rcout<<psi<<" "<<Cef<<" "<<Egd<<" "<<Ect-Egd+Cef<<"\n";
//      psi = psi + psiStep;
//      Egd = Egamma(psi, kxylemmax, c,d, psiPLC);
//      Cef = (V/timestep)*tissueRelativeWaterContent(psi, pi0, eps, psi, c, d,fapo);
//   }
//   return(psi);
// }


//' @rdname hydraulics_supplyfunctions
// [[Rcpp::export("hydraulics_E2psiVanGenuchten")]]
double E2psiVanGenuchten(double E, double psiSoil, double krhizomax, double n, double alpha, 
                         double psiStep = -0.0001, double psiMax = -10.0) {
  if(E<0.0) stop("E has to be positive");
  if(E==0) return(psiSoil);
  double psi = psiSoil;
  double psiPrev = psi;
  double vgPrev = vanGenuchtenConductance(psi, krhizomax, n, alpha);
  double vg = vgPrev;
  double Eg = 0.0;
  while(Eg<E) {
    psiPrev = psi;
    vgPrev = vg;
    psi = psi + psiStep;
    vg = vanGenuchtenConductance(psi, krhizomax, n, alpha);
    Eg = Eg + ((vg+vgPrev)/2.0)*std::abs(psiStep);
    if(psi<psiMax) return(NA_REAL);
  }
  return(psiPrev);
}

//' @rdname hydraulics_supplyfunctions
// [[Rcpp::export("hydraulics_E2psiTwoElements")]]
double E2psiTwoElements(double E, double psiSoil, double krhizomax, double kxylemmax, double n, double alpha, double c, double d, double psiCav = 0.0,
                        double psiStep = -0.0001, double psiMax = -10.0) {
  if(E<0.0) stop("E has to be positive");
  if(E==0) return(psiSoil);
  double psiRoot = E2psiVanGenuchten(E, psiSoil, krhizomax, n, alpha, psiStep, psiMax);
  if(NumericVector::is_na(psiRoot)) return(NA_REAL);
  return(E2psiXylem(E, psiRoot, kxylemmax, c, d, psiCav));
}


double ludcmp(NumericMatrix a, int n, IntegerVector indx) {
    const double TINY=1.0e-20;
    int i,imax=0,j,k;
    double big,dum,sum,temp;
    NumericVector vv(n);
    double d = 1.0;
    for(i=0;i<n;i++) {
      big = 0.0;
      for(j=0;j<n;j++) if((temp=std::abs(a(i,j)))>big) big=temp;
      if(big==0.0) throw std::range_error("Singular matrix in routine ludcmp");
      vv[i] = 1.0/big; //Save the scaling
    }
    //Loop over columns of Crout's method
    for(j=0;j<n;j++){
      for(i=0;i<j;i++) {
        sum=a(i,j);
        for(k=0;k<i;k++) sum-=a(i,k)*a(k,j);
        a(i,j)=sum;
      }
      big=0.0;
      for(i=j;i<n;i++) {
        sum=a(i,j);
        for(k=0;k<j;k++) sum-=a(i,k)*a(k,j);
        a(i,j)=sum;
        if((dum=vv[i]*std::abs(sum))>=big) {
          big=dum;
          imax=i;
        }
      }
      if(j!=imax) {
        for(k=0;k<n;k++){
          dum=a(imax,k);
          a(imax,k) = a(j,k);
          a(j,k) = dum;
        }
        d=-d;
        vv[imax] = vv[j];
      }
      indx[j] = imax;
      if(a(j,j)==0.0) a(j,j) = TINY;
      if(j!=n) {
        dum=1.0/(a(j,j));
        for(i=j+1;i<n;i++) a(i,j)*=dum;
      }
    }
    return(d);
}
void lubksb(NumericMatrix a, int n, IntegerVector indx, NumericVector b) {
  int i,ii=-1,ip,j;
  double sum;
  for(i = 0;i<n;i++){
    ip=indx[i];
    sum=b[ip];
    b[ip]=b[i];
    if(ii>=0) for(j=ii;j<=i-1;j++) sum-=a(i,j)*b[j];
    else if(sum) ii=i;
    b[i] = sum;
  }
  for(i=(n-1);i>=0;i--) {
    sum=b[i];
    for(j=i+1;j<n;j++) sum-=a(i,j)*b[j];
    b[i] = sum/a(i,i);
  }
}

//' @rdname hydraulics_supplyfunctions
// [[Rcpp::export("hydraulics_E2psiBelowground")]]
List E2psiBelowground(double E, List hydraulicNetwork, 
                  NumericVector psiIni = NumericVector::create(0),
                  int ntrial = 10, double psiTol = 0.0001, double ETol = 0.0001) {
  NumericVector psiSoil = hydraulicNetwork["psisoil"]; 
  NumericVector krhizomax = hydraulicNetwork["krhizomax"]; 
  NumericVector nsoil = hydraulicNetwork["nsoil"]; 
  NumericVector alphasoil = hydraulicNetwork["alphasoil"]; 
  NumericVector krootmax = hydraulicNetwork["krootmax"];
  double rootc = hydraulicNetwork["rootc"];
  double rootd  = hydraulicNetwork["rootd"]; 
  
  int nlayers = psiSoil.length();
  //Initialize
  NumericVector x(nlayers+1);
  if(psiIni.size()==(nlayers+1)){
    for(int l=0;l<(nlayers+1);l++) {
      x[l] = psiIni[l];
    }
  } else{
    double minPsi = -0.00001;
    for(int l=0;l<nlayers;l++) {
      x[l] =psiSoil[l];
      minPsi = std::min(minPsi, psiSoil[l]);
      // Rcout<<"("<<x[l]<<") ";
    }
    x[nlayers] = minPsi;
    // Rcout<<"("<<x[nlayers]<<")\n";
    
  }
  
  //Flow across root xylem and rhizosphere elements
  NumericVector Eroot(nlayers), Erhizo(nlayers);
  
  //Newton-Raphson algorithm
  NumericVector p(nlayers+1), fvec(nlayers+1);
  IntegerVector indx(nlayers+1);
  NumericMatrix fjac(nlayers+1,nlayers+1);
  double Esum = 0.0;
  for(int k=0;k<ntrial;k++) {
    // Rcout<<"trial "<<k<<"\n";
    //Calculate steady-state flow functions
    Esum = 0.0;
    bool stop = false;
    for(int l=0;l<nlayers;l++) {
      Eroot[l] = EXylem(x[nlayers], x[l], krootmax[l], rootc, rootd, true, 0.0);
      // Rcout<<"("<<Eroot[l]<<"\n";
      Erhizo[l] = EVanGenuchten(x[l], psiSoil[l], krhizomax[l], nsoil[l], alphasoil[l]);
      fvec[l] = Erhizo[l] - Eroot[l];
      // Rcout<<" Erhizo"<<l<<": "<< Erhizo[l]<<" Eroot"<<l<<": "<<Eroot[l]<<" fvec: "<<fvec[l]<<"\n";
      // Rcout<<"der psi_l "<<d_psi_l<<"der Eroot psi_l "<<d_Eroot_psi_l<<"  der psi_root "<<d_psi_root<<"\n";
      Esum +=Eroot[l];
    }
    fvec[nlayers] = Esum-E;
    // Rcout<<"fvec_nlayers: "<<fvec[nlayers]<<"\n";
    //Fill Jacobian
    for(int l1=0;l1<nlayers;l1++) { //funcio
      for(int l2=0;l2<nlayers;l2++) { //derivada
        if(l1==l2) {
          fjac(l1,l2) = -vanGenuchtenConductance(x[l2],krhizomax[l2], nsoil[l2], alphasoil[l2])-xylemConductance(x[l2], krootmax[l2], rootc, rootd);  
        }
        else fjac(l1,l2) = 0.0;
      }
    }
    fjac(nlayers,nlayers) = 0.0;
    for(int l=0;l<nlayers;l++) { 
      fjac(l,nlayers) = xylemConductance(x[nlayers], krootmax[l], rootc, rootd); //funcio l derivada psi_rootcrown
      fjac(nlayers,l) = xylemConductance(x[l], krootmax[l], rootc, rootd);//funcio nlayers derivada psi_l
      // funcio nlayers derivada psi_rootcrown
      fjac(nlayers,nlayers) +=-xylemConductance(x[nlayers], krootmax[l], rootc, rootd);
    }
    // for(int l1=0;l1<=nlayers;l1++) { //funcio
    //   for(int l2=0;l2<=nlayers;l2++) { //derivada
    //     Rcout<<fjac(l1,l2)<<" ";
    //   }
    //   Rcout<<"\n";
    // }
    //Check function convergence
    double errf = 0.0;
    for(int fi=0;fi<=nlayers;fi++) errf += std::abs(fvec[fi]);
    if(errf<=ETol) break;
    //Right-hand side of linear equations
    for(int fi=0;fi<=nlayers;fi++) p[fi] = -fvec[fi];
    //Solve linear equations using LU decomposition
    ludcmp(fjac,nlayers+1,indx);
    lubksb(fjac,nlayers+1,indx,p);
    //Check root convergence
    double errx = 0.0;
    for(int fi=0;fi<=nlayers;fi++) {
      errx +=std::abs(p[fi]);
      x[fi]+=p[fi];
      x[fi] = std::min(0.0, x[fi]);
      if(x[fi]<-40.0) {
        x[fi] = NA_REAL;
        stop = true;
      }
      // Rcout<<"("<<x[fi]<<") ";
    }
    // Rcout<<"\n";
    if(errx<=psiTol) break;
    else if(k==(ntrial-1)) { //Last trial and no convergence
      for(int fi=0;fi<=nlayers;fi++) x[fi] = NA_REAL;
      // Rcout<<"LC";
      stop = true;
    }
    if(stop) break;
  }
  
  //Initialize and copy output
  NumericVector psiRhizo(nlayers);
  for(int l=0;l<nlayers;l++) {
    psiRhizo[l] = x[l];
  }
  double psiRootCrown = x[nlayers];
  //Calculate final flows
  Esum = 0.0;
  for(int l=0;l<(nlayers-1);l++) {
    Erhizo[l] = EVanGenuchten(x[l], psiSoil[l], krhizomax[l], nsoil[l], alphasoil[l]);
    Esum += Erhizo[l];
  }
  Erhizo[nlayers-1] = E - Esum; //Define as difference to match input
  return(List::create(Named("E") = E, Named("ERhizo")=Erhizo, Named("psiRhizo") = psiRhizo, Named("psiRootCrown") = psiRootCrown, Named("x") = x));
} 


//' @rdname hydraulics_supplyfunctions
// [[Rcpp::export("hydraulics_E2psiAboveground")]]
List E2psiAboveground(double E, double psiRootCrown, 
                      List hydraulicNetwork) {
  double kstemmax = hydraulicNetwork["kstemmax"];
  double stemc = hydraulicNetwork["stemc"];
  double stemd = hydraulicNetwork["stemd"];
  double kleafmax = hydraulicNetwork["kleafmax"];
  double leafc = hydraulicNetwork["leafc"];
  double leafd = hydraulicNetwork["leafd"];
  NumericVector PLCstem = hydraulicNetwork["PLCstem"];
  
  int nStemSegments = PLCstem.size();
  NumericVector psiStem(nStemSegments, 0.0), psiPLCStem(nStemSegments, 0.0);
  double kxsegmax = kstemmax*((double) nStemSegments);
  double psiUp = psiRootCrown;
  for(int i=0;i<nStemSegments; i++) {
    psiPLCStem[i]=  apoplasticWaterPotential(1.0-PLCstem[i], stemc, stemd);
    psiStem[i] = E2psiXylem(E, psiUp, kxsegmax, stemc, stemd, psiPLCStem[i]); //Apliquem la fatiga per cavitacio a la caiguda de potencial a la tija 
    psiUp = psiStem[i];
  }  
  double psiLeaf = E2psiXylem(E, psiStem[nStemSegments-1], kleafmax, leafc, leafd, 0.0); 
  double kterm = xylemConductance(psiLeaf, kleafmax, leafc, leafd);
  return(List::create( Named("E")=E, Named("psiStem") =psiStem,Named("psiLeaf") =psiLeaf,Named("kterm") = kterm));
}


List E2psiAbovegroundCapacitance(double E, double psiRootCrown, 
                         NumericVector psiStemPrev, NumericVector PLCstem,
                         double psiLeafPrev, 
                      double kstemmax, double stemc, double stemd,
                      double kleafmax, double leafc, double leafd,
                      double Vsapwood, double stemfapo, double stempi0, double stemeps,
                      double Vleaf, double leaffapo, double leafpi0, double leafeps,
                      double tstep = 3600.0) {
  int nStemSegments = PLCstem.size();
  
  double Vsegmax = Vsapwood/((double) nStemSegments);
  NumericVector psiStem(nStemSegments, 0.0), psiPLCStem(nStemSegments, 0.0);
  NumericVector EStem(nStemSegments, 0.0);
  double kxsegmax = kstemmax*((double) nStemSegments);
  double psiUp = psiRootCrown;
  double EUp = E;
  for(int i=0;i<nStemSegments; i++) {
    double psiCav = apoplasticWaterPotential(1.0-PLCstem[i], stemc, stemd);
    psiStem[i] = E2psiXylem(EUp, psiUp, 
                            kxsegmax, stemc, stemd, psiCav);
    double Ecap = ECapacitance(psiStem[i],psiStemPrev[i], PLCstem[i],
                               Vsegmax, stemfapo,stemc, stemd,
                               stempi0, stemeps,
                               tstep);                                   
    EStem[i] = EUp+Ecap;
    EUp = EStem[i];
    psiUp = psiStem[i];
  }  
  
  double psiLeaf = E2psiXylem(EUp, psiUp, 
                              kleafmax, leafc, leafd, 0.0);
  double Ecap = ECapacitance(psiLeaf,psiLeafPrev, 0.0,
                             Vleaf, leaffapo,leafc, leafd,
                             leafpi0, leafeps,
                             tstep);             
  double ELeaf = EUp + Ecap;
  double kterm = xylemConductance(psiLeaf, kleafmax, leafc, leafd);
  return(List::create( Named("EStem")=EStem, Named("psiStem") =psiStem,
                       Named("psiLeaf") =psiLeaf,
                       Named("ELeaf") = ELeaf,
                       Named("kterm") = kterm));
}

// List E2psiAbovegroundCapacitance(double E, double psiRootCrown,      
//                       double EPrev, double psiRootCrownPrev,
//                       NumericVector psiStemPrev, NumericVector PLCstemPrev, NumericVector RWCsympstemPrev, 
//                       double psiLeafPrev, double RWCsympleafPrev,
//                       double kstemmax, double stemc, double stemd,
//                       double kleafmax, double leafc, double leafd,
//                       double Vsapwood, double stemfapo, double stempi0, double stemeps,
//                       double Vleaf, double leaffapo, double leafpi0, double leafeps,
//                       double klat,
//                       double tstep = 3600.0, int nSubSteps = 1000) {
// 
// 
//   
//   //Make copy of initial vectors
//   NumericVector psiStem = clone(psiStemPrev);
//   NumericVector PLCstem = clone(PLCstemPrev);
//   NumericVector RWCsympstem = clone(RWCsympstemPrev);
//   double psiLeaf = psiLeafPrev;
//   double RWCsympleaf = RWCsympleafPrev;
//   
//   int n = PLCstem.size();
//   double kxsegmax = kstemmax*((double) n);
//   // double kstoseg = ksto*((double) n);
//   double Vsegmax = Vsapwood/((double) n);
//   double m3tommol = 55555556.0;
//   
//   //Calculate initial apoplastic volumes and water potentials
//   NumericVector VStem(n, NA_REAL);
//   NumericVector psiPLCStem(n, NA_REAL);
//   NumericVector psiStorageStem(n, NA_REAL);
//   double psiLeafSymp = symplasticWaterPotential(RWCsympleaf, leafpi0, leafeps);
//   for(int i=0;i<n;i++) {
//     VStem[i] = Vsegmax*stemfapo*apoplasticRelativeWaterContent(psiStem[i], stemc, stemd);
//     psiPLCStem[i]=  apoplasticWaterPotential(1.0-PLCstem[i], stemc, stemd);
//     psiStorageStem[i] = symplasticWaterPotential(RWCsympstem[i], stempi0, stemeps);
//     // Rcout<< VStem[i] << " "<<psiStem[i] <<" "<<psiPLCStem[i] << " "<<psiStorageStem[i] <<"\n";
//   }
// 
//   double Es = EPrev;
//   double deltaE = (E-EPrev)/((double) (nSubSteps-1));
//   double Efin = 0.0;
//   double tstepsub = tstep/((double) nSubSteps);
//   double psiRootS = psiRootCrownPrev;
//   double deltaPsi = (psiRootCrown-psiRootCrownPrev)/((double) (nSubSteps-1));
//   for(int s=0;s<nSubSteps;s++) {
//     
//     NumericVector FlatStem(n, NA_REAL);
//     // NumericVector Fversym1(n, 0.0), Fversym2(n, 0.0);
//     
//     //SYMPLASTIC FLOWS
//     //Calculate lateral symplastic flow (positive when symplasm has less negative WP)
//     double FlatLeaf = klat*(psiLeafSymp - psiLeaf);
//     //Calculate vertical and lateral symplastic flows among stem segments
//     for(int i=0;i<n;i++) {
//       //Lateral flow
//       FlatStem[i] =  klat*(psiStorageStem[i]-psiStem[i]);
//       //Towards above
//       // if(i<(n-1)) Fversym2[i] = kstoseg*(psiStorageStem[i] - psiStorageStem[i+1]); 
//       // else Fversym2[i] = 0.0;
//       //From below
//       // if(i>0) Fversym1[i] = kstoseg*(psiStorageStem[i-1]-psiStorageStem[i]);
//       // else Fversym1[i] = 0.0;
//       // Rcout<< "Flow "<< i<< " "<<FlatStem[i] << " "<<Fverapo1[i] << " "<<Fverapo2[i] <<" "<<Fversym1[i] << " "<<Fversym2[i] <<"\n";
//     }
//     
//     
//     //Update stem compartments
//     double psiUp = psiRootS;
//     double Ein = Es;
//     for(int i=0;i<n;i++) {
//       //Store previous WP and calculate new one
//       double psiPrev = psiStem[i];
//       psiStem[i] = E2psiXylem(Ein, psiUp, kxsegmax, stemc,stemd, psiPLCStem[i]);
//       if(Vsegmax>0.0) {
//         //Flow from change in volume in the stem apoplasm 
//         double Fapostem = (m3tommol/tstep)*Vsegmax*stemfapo*(apoplasticRelativeWaterContent(psiPrev, stemc, stemd)-apoplasticRelativeWaterContent(psiStem[i], stemc, stemd)); 
//         //Update flow for next segment
//         Ein += (FlatStem[i] + Fapostem);
//         //Symplastic compartment
//         RWCsympstem[i] = (Vsegmax*(1.0-stemfapo)*RWCsympstem[i] - (tstepsub/m3tommol)*(FlatStem[i]))/(Vsegmax*(1.0-stemfapo));
//         // RWCsympstem[i] = (Vsegmax*(1.0-stemfapo)*RWCsympstem[i] + (tstepsub/m3tommol)*(Fversym1[i] - Fversym2[i] - FlatStem[i]))/(Vsegmax*(1.0-stemfapo));
//         psiStorageStem[i] = symplasticWaterPotential(RWCsympstem[i], stempi0, stemeps);
//         // Rcout<< "psi "<<psiStorageStem[i]<<"\n";
//       }
//       //Store water potential for next segment
//       psiUp = psiStem[i];
//     }
//     
//     
//     // Rcout<< "Vertical flow stem to leaf "<<FverStemLeaf<<"\n";
//     // // Rcout<< "Vertical flow leaf to atm "<<FverLeafAtm<<"\n";
//     // Rcout<< "Lateral flow to leaf symplasm "<<FlatLeaf<<"\n";
//     
//     psiLeafPrev = psiLeaf;
//     psiLeaf = E2psiXylem(Ein, psiUp, kleafmax, leafc, leafd,0.0); //apoplasticWaterPotential(Vapoleaf/(Vleaf*leaffapo), leafc, leafd);
//     if(Vleaf>0.0) {
//       //Flow from change in volume in the leaf apoplasm 
//       double Fapoleaf = (m3tommol/tstep)*Vleaf*leaffapo*(apoplasticRelativeWaterContent(psiLeafPrev, leafc, leafd) - apoplasticRelativeWaterContent(psiLeaf, leafc, leafd));
//       //Water balance leaf symplasm 
//       RWCsympleaf = (Vleaf*(1.0-leaffapo)*RWCsympleaf - (tstepsub/m3tommol)*FlatLeaf)/(Vleaf*(1.0-leaffapo));
//       psiLeafSymp = symplasticWaterPotential(RWCsympleaf, leafpi0, leafeps);
//       // Rcout<< "psiLeafSymp"<<psiLeafSymp<<"\n";
//       //Water balance leaf apoplasm 
//       // Rcout<< "new psiLeaf"<<psiLeaf<<"\n";
//       Ein += (FlatLeaf + Fapoleaf);
//     }
//     Efin +=Ein;
//     psiRootS += deltaPsi;
//     Es +=deltaE;
//   }
//   // Rcout<<(Es-deltaE)<<" "<<(psiRootS-deltaPsi)<<"\n";
//   //Update PLC
//   if(Vsegmax>0.0) {
//     for(int i=0;i<n;i++) {
//       PLCstem[i] = std::max(PLCstem[i], 1.0 - apoplasticRelativeWaterContent(psiStem[i], stemc, stemd));
//     }
//   }
//   
//   //Difference between flow with and without capacitance effects
//   double Edif = (Efin/((double)nSubSteps)) - ((EPrev + E)/2.0);
//   // Rcout<<E<< " "<<Edif<<"\n";
//   // newRWCsympleaf = std::min(1.0, newRWCsympleaf);
//   double kterm = xylemConductance(psiLeaf, kleafmax, leafc, leafd);
//   return(List::create( _["E"] = E+Edif,
//                        _["Edif"] = Edif,
//                        _["psiLeaf"] = psiLeaf,
//                        _["psiStem"] = psiStem, 
//                        _["PLCstem"] = PLCstem, 
//                        _["RWCsympstem"] = RWCsympstem,
//                        _["RWCsympleaf"] = RWCsympleaf,
//                        _["kterm"] = kterm));
// }



List E2psiAbovegroundCapacitanceDisconnected(double E,                           
                      NumericVector psiStemPrev, NumericVector PLCstem, NumericVector RWCsympstemPrev, 
                      double psiLeafPrev, double RWCsympleafPrev,
                      double kstemmax, double stemc, double stemd,
                      double kleafmax, double leafc, double leafd,
                      double Vsapwood, double stemfapo, double stempi0, double stemeps,
                      double Vleaf, double leaffapo, double leafpi0, double leafeps,
                      double klat,
                      double tstep = 3600.0) {
  
  int n = PLCstem.size();
  double kxsegmax = kstemmax*((double) n);
  // double kstoseg = ksto*((double) n);
  double Vsegmax = Vsapwood/((double) n);
  double m3tommol = 55555556.0;
  


  //Make copy of initial vectors
  NumericVector RWCsympstem = clone(RWCsympstemPrev);
  NumericVector psiStem = clone(psiStemPrev);
  double psiLeaf = psiLeafPrev;
  double RWCsympleaf = RWCsympleafPrev;
  
  //Calculate initial apoplastic volumes and water potentials
  NumericVector VStem(n, NA_REAL), VStemMax(n, NA_REAL);
  NumericVector psiPLCStem(n, NA_REAL),psiStorageStem(n, NA_REAL);
  double Vapoleafini = Vleaf*leaffapo*apoplasticRelativeWaterContent(psiLeaf, leafc, leafd);
  // Rcout<< "Initial leaf volume "<<Vapoleafini<<"\n";
  double psiLeafSymp = symplasticWaterPotential(RWCsympleaf, leafpi0, leafeps);
  for(int i=0;i<n;i++) {
    VStemMax[i] = Vsegmax*stemfapo; //Maximum apoplastic volume 
    VStem[i]  = Vsegmax*stemfapo*apoplasticRelativeWaterContent(psiStem[i], stemc, stemd);//current apoplastic volume
    psiPLCStem[i] = apoplasticWaterPotential(1.0-PLCstem[i], stemc, stemd);
    psiStorageStem[i] = symplasticWaterPotential(RWCsympstem[i], stempi0, stemeps);
    // Rcout<< "VStem "<<VStem[i] << " VstemMax "<<VStemMax[i]<<" psiPLC "<<psiPLCStem[i] <<"\n";
  }
  double Vleafssub = Vapoleafini;
    
    
  //Substeps of 1 second
  double tstepsub = 1.0;
  double nSubSteps = tstep; //(tstep/((secEmpty*0.001)));

  NumericVector Esub((int)nSubSteps); //Flow every substep
  for(int s=0;s<((int)nSubSteps);s++) {
    Esub[s] = E;
    //Calculate vertical flow from stem to leaf (positive when leaf has more negative WP)
    double FverStemLeaf = xylemConductance(psiLeaf, kleafmax, leafc, leafd)*(psiStem[n-1]-psiLeaf);
    // Rcout<< "Vertical flow stem to leaf "<<FverStemLeaf<<"\n";
    //Calculate lateral symplastic flow (positive when symplasm has less negative WP)
    double FlatLeaf = klat*(psiLeafSymp - psiLeaf);
    // Rcout<< "Lateral flow to leaf symplasm "<<FlatLeaf<<"\n";
    
    //Calculate vertical and lateral flows among stem segments
    NumericVector FlatStem(n, NA_REAL);
    NumericVector Fverapo1(n, 0.0), Fverapo2(n, 0.0);//, Fversym1(n, 0.0), Fversym2(n, 0.0);
    for(int i=0;i<n;i++) {
      //Lateral flow
      FlatStem[i] =  klat*(psiStorageStem[i]-psiStem[i]);
      
      //Towards above
      if(i<(n-1)) {
        Fverapo2[i] = xylemConductance(std::min(psiStem[i], psiPLCStem[i]), kxsegmax, stemc, stemd)*(psiStem[i]-psiStem[i+1]);
        // Fversym2[i] = kstoseg*(psiStorageStem[i] - psiStorageStem[i+1]); 
      } else {
        Fverapo2[i] = FverStemLeaf;
        // Fversym2[i] = 0.0;
      }
      //From below
      if(i>0) {
        Fverapo1[i] = xylemConductance(std::min(psiStem[i-1], psiPLCStem[i-1]), kxsegmax, stemc, stemd)*(psiStem[i-1]-psiStem[i]);
        // Fversym1[i] = kstoseg*(psiStorageStem[i-1]-psiStorageStem[i]);
      } else {
        Fverapo1[i] = 0.0;
        // Fversym1[i] = 0.0;
      }
      // Rcout<< FlatStem[i] << " "<<Fverapo1[i] << " "<<Fverapo2[i] <<" "<<Fversym1[i] << " "<<Fversym2[i] <<"\n";
    }
    
    //Water balance leaf symplasm (volume and water potential)
    RWCsympleaf = RWCsympleaf - ((tstepsub/m3tommol)*FlatLeaf)/(Vleaf*(1.0-leaffapo));
    psiLeafSymp = symplasticWaterPotential(RWCsympleaf, leafpi0, leafeps);

    //Water balance leaf apoplasm (volume and water potential)
    Vleafssub = Vleafssub + (tstepsub/m3tommol)*(FverStemLeaf + FlatLeaf - Esub[s]);
    if(Vleafssub > (Vleaf*leaffapo)) { //Add output flow to close balance
      double dif = Vleafssub - (Vleaf*leaffapo);
      // Rcout<< "Excess "<<dif<<"\n";
      Vleafssub = Vleaf*leaffapo;
      Esub[s] = Esub[s] + dif*(m3tommol/tstepsub);
    }
    psiLeaf = apoplasticWaterPotential(Vleafssub/(Vleaf*leaffapo), leafc, leafd);
    if(NumericVector::is_na(psiLeaf)) psiLeaf = -40.0;
    
    //Update stem compartments
    for(int i=0;i<n;i++) {
      //Apoplastic compartment
      VStem[i] = VStem[i] + (tstepsub/m3tommol)*(Fverapo1[i] - Fverapo2[i] + FlatStem[i]);
      VStem[i] = std::max(0.0,std::min(VStem[i],VStemMax[i]));
      psiStem[i]=  apoplasticWaterPotential((VStem[i]/VStemMax[i]), stemc, stemd);
      if(NumericVector::is_na(psiStem[i])) psiStem[i] = -40.0;
      //Symplastic compartment
      RWCsympstem[i] = RWCsympstem[i] - ((tstepsub/m3tommol)*FlatStem[i])/(Vsegmax*(1.0-stemfapo));
      psiStorageStem[i] = symplasticWaterPotential(RWCsympstem[i], stempi0, stemeps);
    }
  }
  //Final flow
  double Efin = sum(Esub)/((double)nSubSteps); 
  
  // Rcout<<" E "<< E << " Efin "<< Efin << " psiLeaf "<<psiLeaf << "  RWCsympleaf "<< RWCsympleaf;
  // Rcout<< " VStem "<<VStem[0] << " psiStem "<<psiStem[0]<< " RWCsympstem "<< RWCsympstem[0]<<"\n";
  
  return(List::create( _["E"] = E,
                       _["Efin"] = Efin,
                       _["psiStem"] = psiStem,
                       _["psiLeaf"] = psiLeaf,
                       _["kleaf"] = xylemConductance(psiLeaf, kleafmax, leafc, leafd),
                       _["RWCsympstem"] = RWCsympstem,
                       _["RWCsympleaf"] = RWCsympleaf));
}


//' @rdname hydraulics_supplyfunctions
// [[Rcpp::export("hydraulics_E2psiFineRootLeaf")]]
List E2psiFineRootLeaf(double E, double psiFineRoot, 
                       List hydraulicNetwork) {
  
  NumericVector krootmax = hydraulicNetwork["krootmax"];
  double rootc = hydraulicNetwork["rootc"];
  double rootd  = hydraulicNetwork["rootd"]; 
  double kstemmax = hydraulicNetwork["kstemmax"];
  double stemc = hydraulicNetwork["stemc"];
  double stemd = hydraulicNetwork["stemd"];
  double kleafmax = hydraulicNetwork["kleafmax"];
  double leafc = hydraulicNetwork["leafc"];
  double leafd = hydraulicNetwork["leafd"];
  NumericVector PLCstem = hydraulicNetwork["PLCstem"];
  
  double kxsegmax = kstemmax*2.0; //Assume two stem segments
  double psiPLCStem = apoplasticWaterPotential(1.0 - PLCstem[0], stemc, stemd);
  double psiRootCrown = E2psiXylem(E, psiFineRoot, krootmax[0], rootc, rootd);
  double psiStem1 = E2psiXylem(E, psiRootCrown, kxsegmax, stemc, stemd, psiPLCStem); //Apliquem la fatiga per cavitacio a la caiguda de potencial a la tija 
  double psiStem2 = E2psiXylem(E, psiStem1, kxsegmax, stemc, stemd, psiPLCStem); //Apliquem la fatiga per cavitacio a la caiguda de potencial a la tija 
  double psiLeaf = E2psiXylem(E, psiStem2, kleafmax, leafc, leafd, 0.0); 
  double kterm = xylemConductance(psiLeaf, kleafmax, leafc, leafd);
  return(List::create( Named("E")=E, 
                       Named("psiRootCrown") =psiRootCrown,
                       Named("psiStem1") =psiStem1,
                       Named("psiStem2") =psiStem2,
                       Named("psiLeaf") =psiLeaf,Named("kterm") = kterm));
}

//' @rdname hydraulics_supplyfunctions
// [[Rcpp::export("hydraulics_E2psiNetworkStem1")]]
List E2psiNetworkStem1(double E, List hydraulicNetwork,
                  NumericVector psiIni = NumericVector::create(0),
                  int ntrial = 10, 
                  double psiTol = 0.0001, double ETol = 0.0001) {
  double kstemmax = hydraulicNetwork["kstemmax"];
  double stemc = hydraulicNetwork["stemc"];
  double stemd = hydraulicNetwork["stemd"];
  NumericVector PLCstem = hydraulicNetwork["PLCstem"];
  
  List E2psiRS = E2psiBelowground(E, hydraulicNetwork, 
                                  psiIni,
                                  ntrial, psiTol, ETol);
  double psiRootCrown  = E2psiRS["psiRootCrown"];
  NumericVector psiRhizo = E2psiRS["psiRhizo"];  
  NumericVector ERhizo = E2psiRS["ERhizo"];  
  
  double psiStem = NA_REAL;  
  if(!NumericVector::is_na(psiRootCrown)) {
    double kxsegmax = kstemmax*2.0;
    double psiPLCStem =  apoplasticWaterPotential(1.0-PLCstem[0], stemc, stemd);
    psiStem = E2psiXylem(E, psiRootCrown, kxsegmax, stemc, stemd, psiPLCStem); //Apliquem la fatiga per cavitacio a la caiguda de potencial a la tija 
  } 
  return(List::create(Named("E") = E, Named("ERhizo")=ERhizo, Named("psiRhizo") = psiRhizo, 
                      Named("psiRootCrown") = psiRootCrown, Named("psiStem1") = psiStem, 
                      Named("x") = E2psiRS["x"]));
} 

//' @rdname hydraulics_supplyfunctions
// [[Rcpp::export("hydraulics_E2psiNetwork")]]
List E2psiNetwork(double E, List hydraulicNetwork,
                  NumericVector psiIni = NumericVector::create(0),
                  int ntrial = 10, 
                  double psiTol = 0.0001, double ETol = 0.0001) {
  NumericVector PLCstem = hydraulicNetwork["PLCstem"];
  
  int nStemSegments = PLCstem.size();
  List E2psiRS = E2psiBelowground(E, hydraulicNetwork, 
                                 psiIni,
                                 ntrial, psiTol, ETol);
  
  //Copy output
  double psiRootCrown  = E2psiRS["psiRootCrown"];
  NumericVector psiRhizo = E2psiRS["psiRhizo"];  
  NumericVector ERhizo = E2psiRS["ERhizo"];  


  NumericVector psiStem(nStemSegments, NA_REAL);
  double psiLeaf = NA_REAL;
  double kterm = NA_REAL;
  if(!NumericVector::is_na(psiRootCrown)) {
    List E2psiAG = E2psiAboveground(E, psiRootCrown, hydraulicNetwork);
    kterm = E2psiAG["kterm"];
    psiLeaf = E2psiAG["psiLeaf"];
    psiStem = E2psiAG["psiStem"];
  } 
  return(List::create(Named("E") = E, Named("ERhizo")=ERhizo, Named("psiRhizo") = psiRhizo, 
                      Named("psiRootCrown") = psiRootCrown, Named("psiStem") = psiStem, Named("psiLeaf") = psiLeaf, 
                      Named("kterm") = kterm, Named("x") = E2psiRS["x"]));
} 

List E2psiNetworkCapacitance(double E, NumericVector psiSoil, 
                             NumericVector psiStemPrev, NumericVector PLCstem,
                             double psiLeafPrev, 
                             NumericVector krhizomax, NumericVector nsoil, NumericVector alphasoil,
                             NumericVector krootmax, double rootc, double rootd, 
                             double kstemmax, double stemc, double stemd,
                             double kleafmax, double leafc, double leafd,
                             double Vsapwood, double stemfapo, double stempi0, double stemeps,
                             double Vleaf, double leaffapo, double leafpi0, double leafeps,
                             double tstep = 3600.0,
                             NumericVector psiIni = NumericVector::create(0),
                             int ntrial = 10, double psiTol = 0.0001, double ETol = 0.0001) {
  
  List hydraulicNetwork = List::create(_["psisoil"] = psiSoil,
                                       _["krhizomax"] = krhizomax,_["nsoil"] = nsoil,_["alphasoil"] = alphasoil,
                                       _["krootmax"] = krootmax, _["rootc"] = rootc, _["rootd"] = rootd,
                                       _["kstemmax"] = kstemmax, _["stemc"] = stemc, _["stemd"] = stemd,
                                       _["kleafmax"] = kleafmax, _["leafc"] = leafc, _["leafd"] = leafd,
                                       _["PLCstem"] = PLCstem);
  
  int nStemSegments = PLCstem.size();
  List E2psiRS = E2psiBelowground(E, hydraulicNetwork, 
                                  psiIni,
                                  ntrial, psiTol, ETol);
  
  //Copy output
  double psiRootCrown  = E2psiRS["psiRootCrown"];
  NumericVector psiRhizo = E2psiRS["psiRhizo"];  
  NumericVector ERhizo = E2psiRS["ERhizo"];  
  
  
  NumericVector psiStem(nStemSegments, NA_REAL);
  NumericVector EStem(nStemSegments, NA_REAL);
  double psiLeaf = NA_REAL;
  double kterm = NA_REAL;
  double ELeaf = NA_REAL;
  if(!NumericVector::is_na(psiRootCrown)) {
    List E2psiAG = E2psiAbovegroundCapacitance(E, psiRootCrown, 
                                               psiStemPrev, PLCstem,
                                               psiLeafPrev,
                                               kstemmax, stemc, stemd,
                                               kleafmax, leafc, leafd,
                                               Vsapwood, stemfapo, stempi0, stemeps,
                                               Vleaf, leaffapo, leafpi0, leafeps,
                                               tstep);
    kterm = E2psiAG["kterm"];
    psiLeaf = E2psiAG["psiLeaf"];
    psiStem = E2psiAG["psiStem"];
    EStem = E2psiAG["EStem"];
    ELeaf = E2psiAG["ELeaf"];
  } 
  return(List::create(Named("EStem") = EStem, Named("ELeaf") = ELeaf,
                      Named("ERhizo")=ERhizo, Named("psiRhizo") = psiRhizo, 
                      Named("psiRootCrown") = psiRootCrown, Named("psiStem") = psiStem, Named("psiLeaf") = psiLeaf, 
                      Named("kterm") = kterm, Named("x") = E2psiRS["x"]));
} 

//' @rdname hydraulics_supplyfunctions
// [[Rcpp::export("hydraulics_supplyFunctionOneXylem")]]
List supplyFunctionOneXylem(NumericVector psiSoil, NumericVector v,
                            double kstemmax, double stemc, double stemd, double psiCav = 0.0,
                            int maxNsteps=200, double dE=0.01) {
  int nlayers = psiSoil.size();
  NumericVector supplyE(maxNsteps);
  NumericVector supplydEdp(maxNsteps);
  NumericMatrix supplyERhizo(maxNsteps,nlayers);
  NumericVector supplyPsi(maxNsteps);
  
  supplyE[0] = 0;
  for(int l=0;l<nlayers;l++) supplyERhizo[l] = 0.0;
  supplyPsi[0] = averagePsi(psiSoil, v, stemc, stemd);
  NumericVector Psilayers(nlayers);
  //Calculate initial slope
  for(int l=0;l<nlayers;l++) {
    Psilayers[l] = E2psiXylem(dE, psiSoil[l], 
                              kstemmax, stemc,stemd, psiCav);
    // Rcout<<Psilayers[l]<<" ";
  }
  // Rcout<<"\n";
  double psiI = averagePsi(Psilayers, v, stemc, stemd);
  double maxdEdp = dE/std::abs(psiI-supplyPsi[0]);
  // Rcout<<maxdEdp<<"\n";
  
  int nsteps = 1;
  dE = maxdEdp*0.1;
  for(int i=1;i<maxNsteps;i++) {
    supplyE[i] = supplyE[i-1]+dE;
    // Rcout<<supplyE[i]<<" ";
    for(int l=0;l<nlayers;l++) {
      Psilayers[l] = E2psiXylem(supplyE[i], psiSoil[l],
                                kstemmax, stemc,stemd, psiCav);
      // Rcout<<Psilayers[l]<<" ";
    }
    // Rcout<<"\n";
    supplyPsi[i] = averagePsi(Psilayers, v, stemc, stemd);
    for(int l=0;l<nlayers;l++) {
      supplyERhizo(i,l) = supplyE[i]*v[l];
    }
    
    if(!NumericVector::is_na(supplyPsi[i])) {
      if(i==1) {
        supplydEdp[0] = (supplyE[1]-supplyE[0])/std::abs(supplyPsi[1]-supplyPsi[0]);
      } else {
        double d1 = (supplyE[i-1]-supplyE[i-2])/std::abs(supplyPsi[i-1]-supplyPsi[i-2]);
        double d2 = (supplyE[i]-supplyE[i-1])/std::abs(supplyPsi[i]-supplyPsi[i-1]);
        supplydEdp[i-1] = (d1+d2)/2.0;
      }
      dE = supplydEdp[i-1]*0.1;
      nsteps++;
      if(supplydEdp[i-1]<0.01*maxdEdp) break;
    } else {
      break;
    }
  }
  //Calculate last dEdP
  if(nsteps>1) supplydEdp[nsteps-1] = (supplyE[nsteps-1]-supplyE[nsteps-2])/std::abs(supplyPsi[nsteps-1]-supplyPsi[nsteps-2]);
  //Copy values tp nsteps
  NumericVector supplyEDef(nsteps);
  NumericVector supplydEdpDef(nsteps);
  NumericMatrix supplyERhizoDef(nsteps,nlayers);
  NumericVector supplyPsiPlant(nsteps);
  for(int i=0;i<nsteps;i++) {
    supplyEDef[i] = supplyE[i];
    supplydEdpDef[i] = supplydEdp[i];
    supplyPsiPlant[i] = supplyPsi[i];
    for(int l=0;l<nlayers;l++) {
      supplyERhizoDef(i,l) = supplyERhizo(i,l);
    }
  }
  return(List::create(Named("E") = supplyEDef,
                      Named("ERhizo") = supplyERhizoDef,
                      Named("PsiPlant")=supplyPsiPlant,
                      Named("dEdP")=supplydEdpDef));
  
}

//' @rdname hydraulics_supplyfunctions
// [[Rcpp::export("hydraulics_supplyFunctionTwoElements")]]
List supplyFunctionTwoElements(double Emax, double psiSoil, double krhizomax, double kxylemmax, double n, double alpha, double c, double d, 
                               double psiCav = 0.0, 
                               double dE = 0.1, double psiMax = -10.0) {
  dE = std::min(dE,Emax/5.0);
  int maxNsteps = round(Emax/dE)+1;
  NumericVector supplyE(maxNsteps);
  NumericVector supplydEdp(maxNsteps);
  NumericVector supplyFittedE(maxNsteps);
  NumericVector supplyPsiRoot(maxNsteps);
  NumericVector supplyPsiPlant(maxNsteps);
  double Eg1 = 0.0;
  double Eg2 = 0.0;
  double psiStep1 = -0.1;
  double psiStep2 = -0.1;
  double psiRoot = psiSoil;
  double psiPlant = psiSoil;
  double vgPrev = vanGenuchtenConductance(psiRoot, krhizomax, n, alpha);
  double vg = 0.0;
  double wPrev = xylemConductance(std::min(psiCav,psiPlant), kxylemmax, c, d); //conductance can decrease if psiCav < psiPlant
  double w = 0.0;
  double incr = 0.0;
  supplyPsiRoot[0] = psiSoil;
  supplyPsiPlant[0] = psiSoil;
  supplyE[0] = 0.0;
  double psiPrec = -0.000001;
  for(int i=1;i<maxNsteps;i++) {
    supplyE[i] = supplyE[i-1]+dE;
    psiStep1 = -0.01;
    psiRoot = supplyPsiRoot[i-1];
    vgPrev = vanGenuchtenConductance(psiRoot, krhizomax, n, alpha);
    while((psiStep1<psiPrec) && (psiRoot>psiMax))  {
      vg = vanGenuchtenConductance(psiRoot+psiStep1, krhizomax, n, alpha);
      incr = ((vg+vgPrev)/2.0)*std::abs(psiStep1);
      if((Eg1+incr)>supplyE[i]) {
        psiStep1 = psiStep1*0.5;
      } else {
        psiRoot = psiRoot + psiStep1;
        Eg1 = Eg1+incr;
        vgPrev = vg;
      }
    }
    supplyPsiRoot[i] = psiRoot;
    if(supplyPsiRoot[i]<psiMax) supplyPsiRoot[i] = psiMax;
    
    psiStep2 = -0.01;
    Eg2 = 0.0;
    psiPlant = psiRoot;
    wPrev = xylemConductance(std::min(psiCav,psiPlant), kxylemmax, c, d);
    while((psiStep2<psiPrec) && (psiPlant>psiMax))  {
      w = xylemConductance(std::min(psiCav,psiPlant+psiStep2), kxylemmax, c, d);
      incr = ((w+wPrev)/2.0)*std::abs(psiStep2);
      if((Eg2+incr)>supplyE[i]) {
        psiStep2 = psiStep2*0.5;
      } else {
        psiPlant = psiPlant + psiStep2;
        Eg2 = Eg2+incr;
        wPrev = w;
      }
    }
    supplyPsiPlant[i] = psiPlant;
    if(supplyPsiPlant[i]<psiMax) supplyPsiPlant[i] = psiMax;
    supplyFittedE[i] = std::max(supplyFittedE[i-1], Eg2); //Ensure non-decreasing function
    // supplyFittedE[i] = Eg2;
    if(i==1) {
      supplydEdp[0] = (supplyFittedE[1]-supplyFittedE[0])/(supplyPsiPlant[0]-supplyPsiPlant[1]); 
      if((supplyPsiPlant[0]-supplyPsiPlant[1])==0.0) supplydEdp[0] = 0.0;
    }
    else if(i>1) {
      supplydEdp[i-1] = 0.5*(supplyFittedE[i-1]-supplyFittedE[i-2])/(supplyPsiPlant[i-2]-supplyPsiPlant[i-1])+0.5*(supplyFittedE[i]-supplyFittedE[i-1])/(supplyPsiPlant[i-1]-supplyPsiPlant[i]); 
      if((supplyPsiPlant[i-2]-supplyPsiPlant[i-1])==0.0) supplydEdp[i-1] = 0.0;
      else if((supplyPsiPlant[i-1]-supplyPsiPlant[i])==0.0) supplydEdp[i-1] = 0.0;
    }
    if(i==(maxNsteps-1)) {
      supplydEdp[i] = (supplyFittedE[i]-supplyFittedE[i-1])/(supplyPsiPlant[i-1]-supplyPsiPlant[i]); 
      if((supplyPsiPlant[i-1]-supplyPsiPlant[i])==0.0) supplydEdp[i] = 0.0;
    }
  }
  return(List::create(Named("E") = supplyE,
                      Named("FittedE") = supplyFittedE,
                      Named("psiRootCrown")=supplyPsiRoot, 
                      Named("psiPlant")=supplyPsiPlant,
                      Named("dEdP")=supplydEdp));
}

//' @rdname hydraulics_supplyfunctions
// [[Rcpp::export("hydraulics_supplyFunctionThreeElements")]]
List supplyFunctionThreeElements(double Emax, double psiSoil, double krhizomax, double kxylemmax, double kleafmax, 
                                 double n, double alpha, 
                                 double stemc, double stemd, 
                                 double leafc, double leafd,
                                 double psiCav = 0.0, 
                                 double dE = 0.1, double psiMax = -10.0) {
  dE = std::min(dE,Emax/5.0);
  int maxNsteps = round(Emax/dE)+1;
  NumericVector supplyE(maxNsteps);
  NumericVector supplydEdp(maxNsteps);
  NumericVector supplyFittedE(maxNsteps);
  NumericVector supplyPsiRoot(maxNsteps);
  NumericVector supplyPsiStem(maxNsteps);
  NumericVector supplyPsiLeaf(maxNsteps);
  double Eg1 = 0.0;
  double Eg2 = 0.0;
  double Eg3 = 0.0;
  double psiStep1 = -0.1;
  double psiStep2 = -0.1;
  double psiStep3 = -0.1;
  double psiRoot = psiSoil;
  double psiStem = psiSoil;
  double psiLeaf = psiSoil;
  double vgPrev = vanGenuchtenConductance(psiRoot, krhizomax, n, alpha);
  double vg = 0.0;
  //conductance can decrease if psiCav < psiStem/psiLeaf
  double wPrevStem = xylemConductance(std::min(psiCav,psiStem), kxylemmax, stemc, stemd); 
  double wPrevLeaf = xylemConductance(std::min(psiCav,psiLeaf), kleafmax, leafc, leafd); 
  double w = 0.0;
  double incr = 0.0;
  supplyPsiRoot[0] = psiSoil;
  supplyPsiStem[0] = psiSoil;
  supplyPsiLeaf[0] = psiSoil;
  supplyE[0] = 0.0;
  double psiPrec = -0.000001;
  for(int i=1;i<maxNsteps;i++) {
    supplyE[i] = supplyE[i-1]+dE;
    psiStep1 = -0.01;
    
    // Root
    psiRoot = supplyPsiRoot[i-1];
    vgPrev = vanGenuchtenConductance(psiRoot, krhizomax, n, alpha);
    while((psiStep1<psiPrec) && (psiRoot>psiMax))  {
      vg = vanGenuchtenConductance(psiRoot+psiStep1, krhizomax, n, alpha);
      incr = ((vg+vgPrev)/2.0)*std::abs(psiStep1);
      if((Eg1+incr)>supplyE[i]) {
        psiStep1 = psiStep1*0.5;
      } else {
        psiRoot = psiRoot + psiStep1;
        Eg1 = Eg1+incr;
        vgPrev = vg;
      }
    }
    supplyPsiRoot[i] = psiRoot;
    if(supplyPsiRoot[i]<psiMax) supplyPsiRoot[i] = psiMax;
    
    //Stem
    psiStep2 = -0.01;
    Eg2 = 0.0;
    psiStem = psiRoot;
    wPrevStem = xylemConductance(std::min(psiCav,psiStem), kxylemmax, stemc, stemd);
    while((psiStep2<psiPrec) && (psiStem>psiMax))  {
      w = xylemConductance(std::min(psiCav,psiStem+psiStep2), kxylemmax, stemc, stemd);
      incr = ((w+wPrevStem)/2.0)*std::abs(psiStep2);
      if((Eg2+incr)>supplyE[i]) {
        psiStep2 = psiStep2*0.5;
      } else {
        psiStem = psiStem + psiStep2;
        Eg2 = Eg2+incr;
        wPrevStem = w;
      }
    }
    supplyPsiStem[i] = psiStem;
    if(supplyPsiStem[i]<psiMax) supplyPsiStem[i] = psiMax;
    
    //Leaf
    psiStep3 = -0.01;
    Eg3 = 0.0;
    psiLeaf = psiStem;
    wPrevLeaf = xylemConductance(psiLeaf, kleafmax,  leafc, leafd);
    while((psiStep3<psiPrec) && (psiLeaf>psiMax))  {
      w = xylemConductance(psiLeaf+psiStep3, kleafmax,  leafc, leafd);
      incr = ((w+wPrevLeaf)/2.0)*std::abs(psiStep3);
      if((Eg3+incr)>supplyE[i]) {
        psiStep3 = psiStep3*0.5;
      } else {
        psiLeaf = psiLeaf + psiStep3;
        Eg3 = Eg3+incr;
        wPrevLeaf = w;
      }
    }
    supplyPsiLeaf[i] = psiLeaf;
    if(supplyPsiLeaf[i]<psiMax) supplyPsiLeaf[i] = psiMax;
    
    
    //Ensure non-decreasing function
    supplyFittedE[i] = std::max(supplyFittedE[i-1], Eg3); 
    
    
    // supplyFittedE[i] = Eg2;
    if(i==1) {
      supplydEdp[0] = (supplyFittedE[1]-supplyFittedE[0])/(supplyPsiLeaf[0]-supplyPsiLeaf[1]); 
      if((supplyPsiLeaf[0]-supplyPsiLeaf[1])==0.0) supplydEdp[0] = 0.0;
    }
    else if(i>1) {
      supplydEdp[i-1] = 0.5*(supplyFittedE[i-1]-supplyFittedE[i-2])/(supplyPsiLeaf[i-2]-supplyPsiLeaf[i-1])+0.5*(supplyFittedE[i]-supplyFittedE[i-1])/(supplyPsiLeaf[i-1]-supplyPsiLeaf[i]); 
      if((supplyPsiLeaf[i-2]-supplyPsiLeaf[i-1])==0.0) supplydEdp[i-1] = 0.0;
      else if((supplyPsiLeaf[i-1]-supplyPsiLeaf[i])==0.0) supplydEdp[i-1] = 0.0;
    }
    if(i==(maxNsteps-1)) {
      supplydEdp[i] = (supplyFittedE[i]-supplyFittedE[i-1])/(supplyPsiLeaf[i-1]-supplyPsiLeaf[i]); 
      if((supplyPsiLeaf[i-1]-supplyPsiLeaf[i])==0.0) supplydEdp[i] = 0.0;
    }
  }
  return(List::create(Named("E") = supplyE,
                      Named("FittedE") = supplyFittedE,
                      Named("psiRootCrown")=supplyPsiRoot, 
                      Named("psiStem")=supplyPsiStem,
                      Named("psiLeaf")=supplyPsiLeaf,
                      Named("dEdP")=supplydEdp));
}

//' @rdname hydraulics_supplyfunctions
// [[Rcpp::export("hydraulics_supplyFunctionBelowground")]]
List supplyFunctionBelowground(List hydraulicNetwork, 
                           double minFlow = 0.0, int maxNsteps=400, 
                           int ntrial = 10, double psiTol = 0.0001, double ETol = 0.0001,
                           double pCrit = 0.001) {
  NumericVector psiSoil = hydraulicNetwork["psisoil"]; 
  int nlayers = psiSoil.size();
  NumericVector supplyE(maxNsteps);
  NumericVector supplydEdp(maxNsteps);
  NumericMatrix supplyERhizo(maxNsteps,nlayers);
  NumericMatrix supplyPsiRhizo(maxNsteps,nlayers);
  NumericVector supplyPsiRoot(maxNsteps);
  List sol = E2psiBelowground(minFlow, hydraulicNetwork,
                          NumericVector::create(0),
                          ntrial,psiTol, ETol);
  NumericVector solERhizo =sol["ERhizo"];
  NumericVector solPsiRhizo = sol["psiRhizo"];
  supplyERhizo(0,_) = solERhizo;
  supplyPsiRhizo(0,_) = solPsiRhizo;
  supplyE[0] = sol["E"];
  supplyPsiRoot[0] = sol["psiRootCrown"];

  //Calculate initial slope
  List solI = E2psiBelowground(minFlow+ETol*2.0, hydraulicNetwork,
                           sol["x"],
                           ntrial,psiTol, ETol);
  double psiRootI = solI["psiRootCrown"];
  double maxdEdp = (ETol*2.0)/std::abs(psiRootI - supplyPsiRoot[0]);

  int nsteps = 1;
  double dE = std::min(0.05,maxdEdp*0.05);
  for(int i=1;i<maxNsteps;i++) {
    // if(i==3) stop("kk");
    supplyE[i] = supplyE[i-1]+dE;
    sol = E2psiBelowground(supplyE[i], hydraulicNetwork,
                       sol["x"],
                       ntrial,psiTol, ETol);
    solERhizo =sol["ERhizo"];
    solPsiRhizo = sol["psiRhizo"];
    supplyERhizo(i,_) =  solERhizo;
    supplyPsiRhizo(i,_) = solPsiRhizo;
    supplyPsiRoot[i] = sol["psiRootCrown"];

    if(!NumericVector::is_na(supplyPsiRoot[i])) {
      if(i==1) {
        supplydEdp[0] = (supplyE[1]-supplyE[0])/std::abs(supplyPsiRoot[1] - supplyPsiRoot[0]);
      } else {
        double d1 = (supplyE[i-1]-supplyE[i-2])/std::abs(supplyPsiRoot[i-1] - supplyPsiRoot[i-2]);
        double d2 = (supplyE[i]-supplyE[i-1])/std::abs(supplyPsiRoot[i] - supplyPsiRoot[i-1]);
        supplydEdp[i-1] = (d1+d2)/2.0;
      }
      if(supplyE[i]>0.1) dE = std::min(0.1,supplydEdp[i-1]*0.05);
      nsteps++;
      if((supplydEdp[i-1]<(pCrit*maxdEdp)) && (i>5)) break;
    } else {
      break;
    }
  }
  //Calculate last dEdP
  if(nsteps>1) supplydEdp[nsteps-1] = (supplyE[nsteps-1]-supplyE[nsteps-2])/std::abs(supplyPsiRoot[nsteps-1]-supplyPsiRoot[nsteps-2]);
  //Copy values tp nsteps
  NumericVector supplyEDef(nsteps);
  NumericVector supplydEdpDef(nsteps);
  NumericMatrix supplyERhizoDef(nsteps,nlayers);
  NumericMatrix supplyPsiRhizoDef(nsteps,nlayers);
  NumericVector supplyPsiRootDef(nsteps);
  for(int i=0;i<nsteps;i++) {
    supplyEDef[i] = supplyE[i];
    supplydEdpDef[i] = supplydEdp[i];
    supplyPsiRootDef[i] = supplyPsiRoot[i];
    supplyERhizoDef(i,_) = supplyERhizo(i,_);
    supplyPsiRhizoDef(i,_) = supplyPsiRhizo(i,_);
  }
  return(List::create(Named("E") = supplyEDef,
                      Named("ERhizo") = supplyERhizoDef,
                      Named("psiRhizo")=supplyPsiRhizoDef,
                      Named("psiRootCrown")=supplyPsiRootDef,
                      Named("dEdP")=supplydEdpDef));
  
}

//' @rdname hydraulics_supplyfunctions
// [[Rcpp::export("hydraulics_supplyFunctionAboveground")]]
List supplyFunctionAboveground(NumericVector Erootcrown, NumericVector psiRootCrown, 
                               List hydraulicNetwork) {

  

  NumericVector PLCstem = hydraulicNetwork["PLCstem"];
  
  int nStemSegments = PLCstem.size();
  int maxNsteps = Erootcrown.size();
  NumericVector supplyE(maxNsteps);
  NumericVector supplydEdp(maxNsteps);
  NumericVector supplyPsiLeaf(maxNsteps);
  NumericVector supplyKterm(maxNsteps);
  NumericMatrix supplyPsiStem(maxNsteps,nStemSegments);

  
  
  int Nsteps = 0;
  for(int i=0;i<maxNsteps;i++) {
    List sol = E2psiAboveground(Erootcrown[i], psiRootCrown[i],                        
                                hydraulicNetwork);
    NumericVector solPsiStem = sol["psiStem"];
    supplyPsiStem(i,_) = solPsiStem; 
    supplyPsiLeaf[i] = sol["psiLeaf"];
    if(NumericVector::is_na(supplyPsiLeaf[i])) {
      break; 
    } else {
      Nsteps = Nsteps + 1;
    }
    supplyKterm[i] = sol["kterm"];
    supplyE[i] = sol["E"];
    
    if(i==1) {
      supplydEdp[0] = (supplyE[1]-supplyE[0])/std::abs(supplyPsiLeaf[1]-supplyPsiLeaf[0]);
    } else {
      double d1 = (supplyE[i-1]-supplyE[i-2])/std::abs(supplyPsiLeaf[i-1]-supplyPsiLeaf[i-2]);
      double d2 = (supplyE[i]-supplyE[i-1])/std::abs(supplyPsiLeaf[i]-supplyPsiLeaf[i-1]);
      supplydEdp[i-1] = (d1+d2)/2.0;
    }
    
  }
  //Calculate last dEdP
  if(Nsteps>1) supplydEdp[Nsteps-1] = (supplyE[Nsteps-1]-supplyE[Nsteps-2])/std::abs(supplyPsiLeaf[Nsteps-1]-supplyPsiLeaf[Nsteps-2]);
  
  // Rcout<<Nsteps;
  //Copy values tp nsteps
  NumericVector supplyKtermDef(Nsteps);
  NumericVector supplyEDef(Nsteps);
  NumericVector supplydEdpDef(Nsteps);
  NumericVector supplyPsiLeafDef(Nsteps);
  NumericMatrix supplyPsiStemDef(Nsteps,nStemSegments);

  for(int i=0;i<Nsteps;i++) {
    supplyEDef[i] = supplyE[i];
    supplydEdpDef[i] = supplydEdp[i];
    supplyKtermDef[i] = supplyKterm[i];
    supplyPsiLeafDef[i] = supplyPsiLeaf[i];
    supplyPsiStemDef(i,_) = supplyPsiStem(i,_); 
  }
  
  return(List::create(Named("E") = supplyEDef,
                      Named("psiStem")=supplyPsiStemDef,
                      Named("psiLeaf")=supplyPsiLeafDef,
                      Named("dEdP")=supplydEdpDef,
                      Named("kterm") = supplyKtermDef));
  
}


List supplyFunctionAbovegroundCapacitance(NumericVector Erootcrown, NumericVector psiRootCrown,
                                  NumericVector psiStemPrev, NumericVector PLCstemPrev,
                                  double psiLeafPrev, 
                                  double kstemmax, double stemc, double stemd,
                                  double kleafmax, double leafc, double leafd,
                                  double Vsapwood, double stemfapo, double stempi0, double stemeps,
                                  double Vleaf, double leaffapo, double leafpi0, double leafeps,
                                  double tstep = 3600.0) {
  int nnodes = psiStemPrev.size(); // stem nodes + leaf
  int maxNsteps = Erootcrown.size();
  NumericVector supplydEdp(maxNsteps);
  NumericVector supplyE(maxNsteps);
  NumericVector supplyPsiLeaf(maxNsteps);
  NumericVector supplyKterm(maxNsteps);
  NumericMatrix supplyEStem(maxNsteps,nnodes);
  NumericMatrix supplyPsiStem(maxNsteps,nnodes);

  
  int Nsteps = 0;
  for(int i=0;i<maxNsteps;i++) {
    List sol = E2psiAbovegroundCapacitance(Erootcrown[i], psiRootCrown[i],
                                   psiStemPrev, PLCstemPrev,
                                   psiLeafPrev, 
                                   kstemmax, stemc, stemd,
                                   kleafmax, leafc, leafd,
                                   Vsapwood, stemfapo, stempi0, stemeps,
                                   Vleaf, leaffapo, leafpi0, leafeps,
                                   tstep);
    NumericVector solEStem = sol["EStem"];
    NumericVector solNewPsiStem = sol["psiStem"];
    supplyPsiStem(i,_) = solNewPsiStem;
    supplyEStem(i,_) = solEStem;
    supplyE[i] = sol["ELeaf"];
    supplyPsiLeaf[i] = sol["psiLeaf"];
    if(NumericVector::is_na(supplyPsiLeaf[i])) {
      break;
    } else {
      Nsteps = Nsteps + 1;
    }
    supplyKterm[i] = sol["kterm"];

    if(i==1) {
      supplydEdp[0] = (supplyE[1]-supplyE[0])/std::abs(supplyPsiLeaf[1]-supplyPsiLeaf[0]);
    } else {
      double d1 = (supplyE[i-1]-supplyE[i-2])/std::abs(supplyPsiLeaf[i-1]-supplyPsiLeaf[i-2]);
      double d2 = (supplyE[i]-supplyE[i-1])/std::abs(supplyPsiLeaf[i]-supplyPsiLeaf[i-1]);
      supplydEdp[i-1] = (d1+d2)/2.0;
    }
    
  }
  //Calculate last dEdP
  if(Nsteps>1) supplydEdp[Nsteps-1] = (supplyE[Nsteps-1]-supplyE[Nsteps-2])/std::abs(supplyPsiLeaf[Nsteps-1]-supplyPsiLeaf[Nsteps-2]);
  
  // Rcout<<Nsteps;
  //Copy values tp nsteps
  NumericVector supplyKtermDef(Nsteps);
  NumericVector supplyEDef(Nsteps);
  NumericVector supplydEdpDef(Nsteps);
  NumericVector supplyPsiLeafDef(Nsteps);
  NumericMatrix supplyEStemDef(Nsteps,nnodes);
  NumericMatrix supplyPsiStemDef(Nsteps,nnodes);

  for(int i=0;i<Nsteps;i++) {
    supplyEDef[i] = supplyE[i];
    supplydEdpDef[i] = supplydEdp[i];
    supplyKtermDef[i] = supplyKterm[i];
    supplyPsiLeafDef[i] = supplyPsiLeaf[i];
    supplyPsiStemDef(i,_) = supplyPsiStem(i,_);
    supplyEStemDef(i,_) = supplyEStem(i,_);
  }
  
  return(List::create(Named("E") = supplyEDef,
                      Named("EStem")=supplyEStemDef,
                      Named("psiStem")=supplyPsiStemDef,
                      Named("psiLeaf")=supplyPsiLeafDef,
                      Named("dEdP")=supplydEdpDef,
                      Named("kterm") = supplyKtermDef));
  
}


// List supplyFunctionAbovegroundCapacitance(NumericVector Erootcrown, NumericVector psiRootCrown,
//                                double EPrev, double psiRootCrownPrev,
//                                NumericVector psiStemPrev, NumericVector PLCstemPrev, NumericVector RWCsympstemPrev,
//                                double psiLeafPrev, double RWCsympleafPrev,
//                                double kstemmax, double stemc, double stemd,
//                                double kleafmax, double leafc, double leafd,
//                                double Vsapwood, double stemfapo, double stempi0, double stemeps,
//                                double Vleaf, double leaffapo, double leafpi0, double leafeps,
//                                double klat,
//                                double tstep = 3600.0, int nSubSteps = 1000) {
//   int nnodes = psiStemPrev.size(); // stem nodes + leaf
//   int maxNsteps = Erootcrown.size();
//   NumericVector supplyE(maxNsteps);
//   NumericVector supplyEdif(maxNsteps);
//   NumericVector supplydEdp(maxNsteps);
//   NumericVector supplyPsiLeaf(maxNsteps);
//   NumericVector supplyRWCsympleaf(maxNsteps);
//   NumericVector supplyKterm(maxNsteps);
//   NumericMatrix supplyPsiStem(maxNsteps,nnodes);
//   NumericMatrix supplyPLCstem(maxNsteps,nnodes);
//   NumericMatrix supplyRWCsympstem(maxNsteps,nnodes);
// 
// 
// 
//   int Nsteps = 0;
//   for(int i=0;i<maxNsteps;i++) {
//     List sol = E2psiAbovegroundCapacitance(Erootcrown[i], psiRootCrown[i],
//                                 EPrev, psiRootCrownPrev,
//                                 psiStemPrev, PLCstemPrev, RWCsympstemPrev,
//                                 psiLeafPrev, RWCsympleafPrev,
//                                 kstemmax, stemc, stemd,
//                                 kleafmax, leafc, leafd,
//                                 Vsapwood, stemfapo, stempi0, stemeps,
//                                 Vleaf, leaffapo, leafpi0, leafeps,
//                                 klat, 
//                                 tstep, nSubSteps);
//     NumericVector solNewPsiStem = sol["psiStem"];
//     NumericVector solNewPLC = sol["PLCstem"];
//     NumericVector solNewRWCsympstem = sol["RWCsympstem"];
//     supplyPsiStem(i,_) = solNewPsiStem;
//     supplyPLCstem(i,_) = solNewPLC;
//     supplyRWCsympstem(i,_) = solNewRWCsympstem;
//     supplyRWCsympleaf[i] = sol["RWCsympleaf"];
//     supplyPsiLeaf[i] = sol["psiLeaf"];
//     supplyEdif[i] = sol["Edif"];
//     if(NumericVector::is_na(supplyPsiLeaf[i])) {
//       break;
//     } else {
//       Nsteps = Nsteps + 1;
//     }
//     supplyKterm[i] = sol["kterm"];
//     supplyE[i] = sol["E"];
// 
//     if(i==1) {
//       supplydEdp[0] = (supplyE[1]-supplyE[0])/std::abs(supplyPsiLeaf[1]-supplyPsiLeaf[0]);
//     } else {
//       double d1 = (supplyE[i-1]-supplyE[i-2])/std::abs(supplyPsiLeaf[i-1]-supplyPsiLeaf[i-2]);
//       double d2 = (supplyE[i]-supplyE[i-1])/std::abs(supplyPsiLeaf[i]-supplyPsiLeaf[i-1]);
//       supplydEdp[i-1] = (d1+d2)/2.0;
//     }
// 
//   }
//   //Calculate last dEdP
//   if(Nsteps>1) supplydEdp[Nsteps-1] = (supplyE[Nsteps-1]-supplyE[Nsteps-2])/std::abs(supplyPsiLeaf[Nsteps-1]-supplyPsiLeaf[Nsteps-2]);
// 
//   // Rcout<<Nsteps;
//   //Copy values tp nsteps
//   NumericVector supplyKtermDef(Nsteps);
//   NumericVector supplyEDef(Nsteps);
//   NumericVector supplyEdifDef(Nsteps);
//   NumericVector supplydEdpDef(Nsteps);
//   NumericVector supplyPsiLeafDef(Nsteps);
//   NumericVector supplyRWCsympleafDef(Nsteps);
//   NumericMatrix supplyPsiStemDef(Nsteps,nnodes);
//   NumericMatrix supplyPLCstemDef(Nsteps,nnodes);
//   NumericMatrix supplyRWCsympstemDef(Nsteps,nnodes);
// 
//   for(int i=0;i<Nsteps;i++) {
//     supplyEDef[i] = supplyE[i];
//     supplyEdifDef[i] = supplyEdif[i];
//     supplydEdpDef[i] = supplydEdp[i];
//     supplyKtermDef[i] = supplyKterm[i];
//     supplyPsiLeafDef[i] = supplyPsiLeaf[i];
//     supplyRWCsympleafDef[i]= supplyRWCsympleaf[i];
//     supplyPsiStemDef(i,_) = supplyPsiStem(i,_);
//     supplyPLCstemDef(i,_) = supplyPLCstem(i,_);
//     supplyRWCsympstemDef(i,_) = supplyRWCsympstem(i,_);
//   }
// 
//   return(List::create(Named("E") = supplyEDef,
//                       Named("Edif") = supplyEdifDef,
//                       Named("PLCstem")=supplyPLCstemDef,
//                       Named("RWCsympstem")=supplyRWCsympstemDef,
//                       Named("RWCsympleaf")=supplyRWCsympleafDef,
//                       Named("psiStem")=supplyPsiStemDef,
//                       Named("psiLeaf")=supplyPsiLeafDef,
//                       Named("dEdP")=supplydEdpDef,
//                       Named("kterm") = supplyKtermDef));
// 
// }

//' @rdname hydraulics_supplyfunctions
// [[Rcpp::export("hydraulics_supplyFunctionFineRootLeaf")]]
List supplyFunctionFineRootLeaf(double psiFineRoot,
                                List hydraulicNetwork,
                                double minFlow = 0.0, int maxNsteps=400, 
                                double ETol = 0.0001, double pCrit = 0.001) {
  NumericVector supplyE(maxNsteps);
  NumericVector supplydEdp(maxNsteps);
  NumericVector supplyPsiRootCrown(maxNsteps);
  NumericVector supplyPsiStem1(maxNsteps);
  NumericVector supplyPsiStem2(maxNsteps);
  NumericVector supplyPsiLeaf(maxNsteps);
  NumericVector supplyKterm(maxNsteps);
  
  List sol = E2psiFineRootLeaf(minFlow, psiFineRoot, 
                               hydraulicNetwork);
  supplyPsiRootCrown[0] = sol["psiRootCrown"];
  supplyPsiStem1[0] = sol["psiStem1"];
  supplyPsiStem2[0] = sol["psiStem2"];
  supplyPsiLeaf[0] = sol["psiLeaf"];
  supplyKterm[0] = sol["kterm"];
  supplyE[0] = minFlow;
  
  //Calculate initial slope
  List solI = E2psiFineRootLeaf(minFlow+ETol*2.0, psiFineRoot, 
                                hydraulicNetwork);
  double psiLeafI = solI["psiLeaf"];
  double maxdEdp = (ETol*2.0)/std::abs(psiLeafI - supplyPsiLeaf[0]);
  
  int nsteps = 1;
  double dE = std::min(0.0005,maxdEdp*0.05);
  for(int i=1;i<maxNsteps;i++) {
    // if(i==3) stop("kk");
    supplyE[i] = supplyE[i-1]+dE;
    sol = E2psiFineRootLeaf(supplyE[i], psiFineRoot, 
                            hydraulicNetwork);
    supplyPsiRootCrown[i] = sol["psiRootCrown"];
    supplyPsiStem1[i] = sol["psiStem1"];
    supplyPsiStem2[i] = sol["psiStem2"];
    supplyPsiLeaf[i] = sol["psiLeaf"];
    supplyKterm[i] = sol["kterm"];
    
    if(!NumericVector::is_na(supplyPsiLeaf[i])) {
      if(i==1) {
        supplydEdp[0] = (supplyE[1]-supplyE[0])/std::abs(supplyPsiLeaf[1] - supplyPsiLeaf[0]);
      } else {
        double d1 = (supplyE[i-1]-supplyE[i-2])/std::abs(supplyPsiLeaf[i-1] - supplyPsiLeaf[i-2]);
        double d2 = (supplyE[i]-supplyE[i-1])/std::abs(supplyPsiLeaf[i] - supplyPsiLeaf[i-1]);
        supplydEdp[i-1] = (d1+d2)/2.0;
      }
      if(supplyE[i]>0.1) dE = std::min(0.05,supplydEdp[i-1]*0.05);
      else if(supplyE[i]>0.05) dE = std::min(0.01,supplydEdp[i-1]*0.05);
      else if(supplyE[i]>0.01) dE = std::min(0.005,supplydEdp[i-1]*0.05);
      nsteps++;
      if(supplydEdp[i-1]<(pCrit*maxdEdp)) break;
    } else {
      break;
    }
  }
  //Calculate last dEdP
  if(nsteps>1) supplydEdp[nsteps-1] = (supplyE[nsteps-1]-supplyE[nsteps-2])/std::abs(supplyPsiLeaf[nsteps-1] - supplyPsiLeaf[nsteps-2]);
  //Copy values tp nsteps
  NumericVector supplyKtermDef(nsteps);
  NumericVector supplyEDef(nsteps);
  NumericVector supplydEdpDef(nsteps);
  NumericVector supplyPsiRootCrownDef(nsteps);
  NumericVector supplyPsiStem1Def(nsteps);
  NumericVector supplyPsiStem2Def(nsteps);
  NumericVector supplyPsiLeafDef(nsteps);
  NumericVector supplyPsiRootDef(nsteps);
  for(int i=0;i<nsteps;i++) {
    if(NumericVector::is_na(supplyE[i])) stop("NA E in supplyFunctionNetwork");
    supplyEDef[i] = supplyE[i];
    supplydEdpDef[i] = supplydEdp[i];
    supplyKtermDef[i] = supplyKterm[i];
    supplyPsiRootCrownDef[i] = supplyPsiRootCrown[i];
    supplyPsiStem1Def[i] = supplyPsiStem1[i];
    supplyPsiStem2Def[i] = supplyPsiStem2[i];
    supplyPsiLeafDef[i] = supplyPsiLeaf[i];
  }
  return(List::create(Named("E") = supplyEDef,
                      Named("psiRootCrown")=supplyPsiRootCrownDef,
                      Named("psiStem1")=supplyPsiStem2Def,
                      Named("psiStem2")=supplyPsiStem2Def,
                      Named("psiLeaf")=supplyPsiLeafDef,
                      Named("dEdP")=supplydEdpDef,
                      Named("kterm") = supplyKtermDef));
  
}

//' @rdname hydraulics_supplyfunctions
// [[Rcpp::export("hydraulics_supplyFunctionNetworkStem1")]]
List supplyFunctionNetworkStem1(List hydraulicNetwork,
                           double minFlow = 0.0, int maxNsteps=400, 
                           int ntrial = 200, double psiTol = 0.0001, double ETol = 0.0001,
                           double pCrit = 0.001) {
  NumericVector psiSoil = hydraulicNetwork["psisoil"]; 
  int nlayers = psiSoil.size();
  NumericVector supplyE(maxNsteps);
  NumericVector supplydEdp(maxNsteps);
  NumericMatrix supplyERhizo(maxNsteps,nlayers);
  NumericMatrix supplyPsiRhizo(maxNsteps,nlayers);
  NumericVector supplyPsiRoot(maxNsteps);
  NumericVector supplyPsiStem1(maxNsteps);

  List sol = E2psiNetworkStem1(minFlow, hydraulicNetwork,
                               NumericVector::create(0),
                               ntrial,psiTol, ETol);
  NumericVector solERhizo = sol["ERhizo"];
  NumericVector solPsiRhizo = sol["psiRhizo"];
  supplyERhizo(0,_) = solERhizo;
  supplyPsiRhizo(0,_) = solPsiRhizo;
  supplyPsiStem1[0] = sol["psiStem1"];
  supplyPsiRoot[0] = sol["psiRootCrown"];
  supplyE[0] = minFlow;
  
  //Calculate initial slope
  List solI = E2psiNetworkStem1(minFlow+ETol*2.0, hydraulicNetwork,
                                sol["x"],
                                ntrial,psiTol, ETol);
  double psiStem1I = solI["psiStem1"];
  double maxdEdp = (ETol*2.0)/std::abs(psiStem1I - supplyPsiStem1[0]);
  
  int nsteps = 1;
  double dE = std::min(0.0005,maxdEdp*0.05);
  for(int i=1;i<maxNsteps;i++) {
    // if(i==3) stop("kk");
    supplyE[i] = supplyE[i-1]+dE;
    sol = E2psiNetworkStem1(supplyE[i], hydraulicNetwork,
                            sol["x"],
                            ntrial,psiTol, ETol);
    solERhizo = sol["ERhizo"];
    solPsiRhizo = sol["psiRhizo"];
    supplyERhizo(i,_) = solERhizo;
    supplyPsiRhizo(i,_) = solPsiRhizo;
    supplyPsiStem1[i] = sol["psiStem1"];
    supplyPsiRoot[i] = sol["psiRootCrown"];

    if(!NumericVector::is_na(supplyPsiStem1[i])) {
      if(i==1) {
        supplydEdp[0] = (supplyE[1]-supplyE[0])/std::abs(supplyPsiStem1[1] - supplyPsiStem1[0]);
      } else {
        double d1 = (supplyE[i-1]-supplyE[i-2])/std::abs(supplyPsiStem1[i-1] - supplyPsiStem1[i-2]);
        double d2 = (supplyE[i]-supplyE[i-1])/std::abs(supplyPsiStem1[i] - supplyPsiStem1[i-1]);
        supplydEdp[i-1] = (d1+d2)/2.0;
      }
      if(supplyE[i]>0.1) dE = std::min(0.05,supplydEdp[i-1]*0.05);
      else if(supplyE[i]>0.05) dE = std::min(0.01,supplydEdp[i-1]*0.05);
      else if(supplyE[i]>0.01) dE = std::min(0.005,supplydEdp[i-1]*0.05);
      nsteps++;
      if(supplydEdp[i-1]<(pCrit*maxdEdp)) break;
    } else {
      break;
    }
  }
  //Calculate last dEdP
  if(nsteps>1) supplydEdp[nsteps-1] = (supplyE[nsteps-1]-supplyE[nsteps-2])/std::abs(supplyPsiStem1[nsteps-1] - supplyPsiStem1[nsteps-2]);
  //Copy values tp nsteps
  NumericVector supplyEDef(nsteps);
  NumericVector supplydEdpDef(nsteps);
  NumericMatrix supplyERhizoDef(nsteps,nlayers);
  NumericMatrix supplyPsiRhizoDef(nsteps,nlayers);
  NumericVector supplyPsiStem1Def(nsteps);
  NumericVector supplyPsiRootDef(nsteps);
  for(int i=0;i<nsteps;i++) {
    if(NumericVector::is_na(supplyE[i])) stop("NA E in supplyFunctionNetwork");
    supplyEDef[i] = supplyE[i];
    supplydEdpDef[i] = supplydEdp[i];
    supplyPsiRootDef[i] = supplyPsiRoot[i];
    supplyERhizoDef(i,_) = supplyERhizo(i,_);
    supplyPsiRhizoDef(i,_) = supplyPsiRhizo(i,_);
    supplyPsiStem1Def[i] = supplyPsiStem1[i];
  }
  return(List::create(Named("E") = supplyEDef,
                      Named("ERhizo") = supplyERhizoDef,
                      Named("psiRhizo")=supplyPsiRhizoDef,
                      Named("psiRootCrown")=supplyPsiRootDef,
                      Named("psiStem1")=supplyPsiStem1Def,
                      Named("dEdP")=supplydEdpDef));
  
}

//' @rdname hydraulics_supplyfunctions
// [[Rcpp::export("hydraulics_supplyFunctionNetwork")]]
List supplyFunctionNetwork(List hydraulicNetwork,
                           double minFlow = 0.0, int maxNsteps=400, 
                           int ntrial = 200, double psiTol = 0.0001, double ETol = 0.0001,
                           double pCrit = 0.001) {
  NumericVector psiSoil = hydraulicNetwork["psisoil"]; 
  NumericVector PLCstem = hydraulicNetwork["PLCstem"]; 
  int nlayers = psiSoil.size();
  int nStemSegments = PLCstem.size();
  NumericVector supplyE(maxNsteps);
  NumericVector supplydEdp(maxNsteps);
  NumericMatrix supplyERhizo(maxNsteps,nlayers);
  NumericMatrix supplyPsiRhizo(maxNsteps,nlayers);
  NumericVector supplyPsiRoot(maxNsteps);
  NumericMatrix supplyPsiStem(maxNsteps,nStemSegments);
  NumericVector supplyPsiLeaf(maxNsteps);
  NumericVector supplyKterm(maxNsteps);
  
  List sol = E2psiNetwork(minFlow, hydraulicNetwork,
                          NumericVector::create(0),
                          ntrial,psiTol, ETol);
  NumericVector solERhizo = sol["ERhizo"];
  NumericVector solPsiRhizo = sol["psiRhizo"];
  NumericVector solPsiStem = sol["psiStem"];
  supplyERhizo(0,_) = solERhizo;
  supplyPsiRhizo(0,_) = solPsiRhizo;
  supplyPsiStem(0,_) = solPsiStem;
  supplyPsiLeaf[0] = sol["psiLeaf"];
  supplyPsiRoot[0] = sol["psiRootCrown"];
  supplyKterm[0] = sol["kterm"];
  supplyE[0] = minFlow;
  
  //Calculate initial slope
  List solI = E2psiNetwork(minFlow+ETol*2.0, hydraulicNetwork,
                           sol["x"],
                              ntrial,psiTol, ETol);
  double psiLeafI = solI["psiLeaf"];
  double maxdEdp = (ETol*2.0)/std::abs(psiLeafI - supplyPsiLeaf[0]);
  
  int nsteps = 1;
  double dE = std::min(0.0005,maxdEdp*0.05);
  for(int i=1;i<maxNsteps;i++) {
    // if(i==3) stop("kk");
    supplyE[i] = supplyE[i-1]+dE;
    sol = E2psiNetwork(supplyE[i], hydraulicNetwork,
                       sol["x"],
                          ntrial,psiTol, ETol);
    solERhizo = sol["ERhizo"];
    solPsiRhizo = sol["psiRhizo"];
    solPsiStem = sol["psiStem"];
    supplyERhizo(i,_) = solERhizo;
    supplyPsiRhizo(i,_) = solPsiRhizo;
    supplyPsiStem(i,_) = solPsiStem;
    supplyPsiLeaf[i] = sol["psiLeaf"];
    supplyPsiRoot[i] = sol["psiRootCrown"];
    supplyKterm[i] = sol["kterm"];
    
    if(!NumericVector::is_na(supplyPsiLeaf[i])) {
      if(i==1) {
        supplydEdp[0] = (supplyE[1]-supplyE[0])/std::abs(supplyPsiLeaf[1] - supplyPsiLeaf[0]);
      } else {
        double d1 = (supplyE[i-1]-supplyE[i-2])/std::abs(supplyPsiLeaf[i-1] - supplyPsiLeaf[i-2]);
        double d2 = (supplyE[i]-supplyE[i-1])/std::abs(supplyPsiLeaf[i] - supplyPsiLeaf[i-1]);
        supplydEdp[i-1] = (d1+d2)/2.0;
      }
      if(supplyE[i]>0.1) dE = std::min(0.05,supplydEdp[i-1]*0.05);
      else if(supplyE[i]>0.05) dE = std::min(0.01,supplydEdp[i-1]*0.05);
      else if(supplyE[i]>0.01) dE = std::min(0.005,supplydEdp[i-1]*0.05);
      nsteps++;
      if(supplydEdp[i-1]<(pCrit*maxdEdp)) break;
    } else {
      break;
    }
  }
  //Calculate last dEdP
  if(nsteps>1) supplydEdp[nsteps-1] = (supplyE[nsteps-1]-supplyE[nsteps-2])/std::abs(supplyPsiLeaf[nsteps-1] - supplyPsiLeaf[nsteps-2]);
  //Copy values tp nsteps
  NumericVector supplyKtermDef(nsteps);
  NumericVector supplyEDef(nsteps);
  NumericVector supplydEdpDef(nsteps);
  NumericMatrix supplyERhizoDef(nsteps,nlayers);
  NumericMatrix supplyPsiRhizoDef(nsteps,nlayers);
  NumericMatrix supplyPsiStemDef(nsteps,nStemSegments);
  NumericVector supplyPsiLeafDef(nsteps);
  NumericVector supplyPsiRootDef(nsteps);
  for(int i=0;i<nsteps;i++) {
    if(NumericVector::is_na(supplyE[i])) stop("NA E in supplyFunctionNetwork");
    supplyEDef[i] = supplyE[i];
    supplydEdpDef[i] = supplydEdp[i];
    supplyKtermDef[i] = supplyKterm[i];
    supplyPsiRootDef[i] = supplyPsiRoot[i];
    supplyERhizoDef(i,_) = supplyERhizo(i,_);
    supplyPsiRhizoDef(i,_) = supplyPsiRhizo(i,_);
    supplyPsiStemDef(i,_) = supplyPsiStem(i,_);
    supplyPsiLeafDef[i] = supplyPsiLeaf[i];
  }
  return(List::create(Named("E") = supplyEDef,
                      Named("ERhizo") = supplyERhizoDef,
                      Named("psiRhizo")=supplyPsiRhizoDef,
                      Named("psiRootCrown")=supplyPsiRootDef,
                      Named("psiStem")=supplyPsiStemDef,
                      Named("psiLeaf")=supplyPsiLeafDef,
                      Named("dEdP")=supplydEdpDef,
                      Named("kterm") = supplyKtermDef));
  
}

List supplyFunctionNetworkCapacitance(NumericVector psiSoil, 
                                      NumericVector psiStemPrev, NumericVector PLCstemPrev,
                                      double psiLeafPrev, 
                                      NumericVector krhizomax, NumericVector nsoil, NumericVector alphasoil,
                                      NumericVector krootmax, double rootc, double rootd, 
                                      double kstemmax, double stemc, double stemd,
                                      double kleafmax, double leafc, double leafd,
                                      double Vsapwood, double stemfapo, double stempi0, double stemeps,
                                      double Vleaf, double leaffapo, double leafpi0, double leafeps,
                                      double tstep = 3600.0,
                                      double minFlow = 0.0, int maxNsteps=400, 
                                      int ntrial = 200, double psiTol = 0.0001, double ETol = 0.0001,
                                      double pCrit = 0.001) {
  int nlayers = psiSoil.size();
  int nStemSegments = PLCstemPrev.size();
  NumericVector supplyE(maxNsteps);
  NumericVector supplyELeaf(maxNsteps);
  NumericVector supplydEdp(maxNsteps);
  NumericMatrix supplyERhizo(maxNsteps,nlayers);
  NumericMatrix supplyPsiRhizo(maxNsteps,nlayers);
  NumericVector supplyPsiRoot(maxNsteps);
  NumericMatrix supplyPsiStem(maxNsteps,nStemSegments);
  NumericMatrix supplyEStem(maxNsteps,nStemSegments);
  NumericVector supplyPsiLeaf(maxNsteps);
  NumericVector supplyKterm(maxNsteps);
  
  supplyE[0] = minFlow;
  List sol = E2psiNetworkCapacitance(minFlow, psiSoil,
                                     psiStemPrev, PLCstemPrev,
                                     psiLeafPrev,  
                                     krhizomax,  nsoil,  alphasoil,
                                     krootmax,  rootc,  rootd,
                                     kstemmax,  stemc,  stemd,
                                     kleafmax,  leafc,  leafd,
                                     Vsapwood, stemfapo, stempi0, stemeps,
                                     Vleaf, leaffapo, leafpi0, leafeps,
                                     tstep,
                                     NumericVector::create(0),
                                     ntrial,psiTol, ETol);
  NumericVector solERhizo = sol["ERhizo"];
  NumericVector solPsiRhizo = sol["psiRhizo"];
  NumericVector solPsiStem = sol["psiStem"];
  NumericVector solEStem = sol["psiStem"];
  supplyERhizo(0,_) = solERhizo;
  supplyPsiRhizo(0,_) = solPsiRhizo;
  supplyPsiStem(0,_) = solPsiStem;
  supplyEStem(0,_) = solEStem;
  supplyPsiLeaf[0] = sol["psiLeaf"];
  supplyPsiRoot[0] = sol["psiRootCrown"];
  supplyKterm[0] = sol["kterm"];
  supplyELeaf[0] =sol["ELeaf"];
  
  //Calculate initial slope
  List solI = E2psiNetworkCapacitance(minFlow+ETol*2.0, psiSoil,
                                      psiStemPrev, PLCstemPrev,
                                      psiLeafPrev,  
                                      krhizomax,  nsoil,  alphasoil,
                                      krootmax,  rootc,  rootd,
                                      kstemmax,  stemc,  stemd,
                                      kleafmax,  leafc,  leafd,
                                      Vsapwood, stemfapo, stempi0, stemeps,
                                      Vleaf, leaffapo, leafpi0, leafeps,
                                      tstep,
                                      sol["x"],
                                      ntrial,psiTol, ETol);
  double psiLeafI = solI["psiLeaf"];
  double EI = solI["ELeaf"];
  double maxdEdp = (EI - supplyELeaf[0])/std::abs(psiLeafI - supplyPsiLeaf[0]);
  
  int nsteps = 1;
  double dE = std::min(0.05,maxdEdp*0.05);
  for(int i=1;i<maxNsteps;i++) {
    // if(i==3) stop("kk");
    supplyE[i] = supplyE[i-1]+dE;
    sol = E2psiNetworkCapacitance(supplyE[i], psiSoil,
                                  psiStemPrev, PLCstemPrev,
                                  psiLeafPrev,  
                                  krhizomax,  nsoil,  alphasoil,
                                  krootmax,  rootc,  rootd,
                                  kstemmax,  stemc,  stemd,
                                  kleafmax,  leafc,  leafd,
                                  Vsapwood, stemfapo, stempi0, stemeps,
                                  Vleaf, leaffapo, leafpi0, leafeps,
                                  tstep,
                                  sol["x"],
                                  ntrial,psiTol, ETol);
    supplyELeaf[i] = sol["ELeaf"];
    solERhizo = sol["ERhizo"];
    solPsiRhizo = sol["psiRhizo"];
    solPsiStem = sol["psiStem"];
    solEStem = sol["EStem"];
    supplyERhizo(i,_) = solERhizo;
    supplyPsiRhizo(i,_) = solPsiRhizo;
    supplyPsiStem(i,_) = solPsiStem;
    supplyEStem(i,_) = solEStem;
    supplyPsiLeaf[i] = sol["psiLeaf"];
    supplyPsiRoot[i] = sol["psiRootCrown"];
    supplyKterm[i] = sol["kterm"];
    
    // Rcout<<supplyPsiLeaf[i]<<"\n";
    if(!NumericVector::is_na(supplyPsiLeaf[i])) {
      if(i==1) {
        supplydEdp[0] = (supplyELeaf[1]-supplyELeaf[0])/std::abs(supplyPsiLeaf[1] - supplyPsiLeaf[0]);
      } else {
        double d1 = (supplyELeaf[i-1]-supplyELeaf[i-2])/std::abs(supplyPsiLeaf[i-1] - supplyPsiLeaf[i-2]);
        double d2 = (supplyELeaf[i]-supplyELeaf[i-1])/std::abs(supplyPsiLeaf[i] - supplyPsiLeaf[i-1]);
        supplydEdp[i-1] = (d1+d2)/2.0;
      }
      if(supplyELeaf[i]>0.1) dE = std::min(0.05,supplydEdp[i-1]*0.05);
      nsteps++;
      // Rcout<<supplydEdp[i-1]<<" "<<pCrit*maxdEdp<<" "<<maxNsteps<<"\n";
      if(supplydEdp[i-1]<(pCrit*maxdEdp)) break;
    } else {
      break;
    }
  }
  //Calculate last dEdP
  if(nsteps>1) supplydEdp[nsteps-1] = (supplyELeaf[nsteps-1]-supplyELeaf[nsteps-2])/std::abs(supplyPsiLeaf[nsteps-1] - supplyPsiLeaf[nsteps-2]);
  //Copy values tp nsteps
  NumericVector supplyKtermDef(nsteps);
  NumericVector supplyEDef(nsteps);
  NumericVector supplyELeafDef(nsteps);
  NumericVector supplydEdpDef(nsteps);
  NumericMatrix supplyERhizoDef(nsteps,nlayers);
  NumericMatrix supplyPsiRhizoDef(nsteps,nlayers);
  NumericMatrix supplyEStemDef(nsteps,nStemSegments);
  NumericMatrix supplyPsiStemDef(nsteps,nStemSegments);
  NumericVector supplyPsiLeafDef(nsteps);
  NumericVector supplyPsiRootDef(nsteps);
  for(int i=0;i<nsteps;i++) {
    if(NumericVector::is_na(supplyE[i])) stop("NA E in supplyFunctionNetwork");
    supplyEDef[i] = supplyE[i];
    supplyELeafDef[i] = supplyELeaf[i];
    supplydEdpDef[i] = supplydEdp[i];
    supplyKtermDef[i] = supplyKterm[i];
    supplyPsiRootDef[i] = supplyPsiRoot[i];
    supplyERhizoDef(i,_) = supplyERhizo(i,_);
    supplyPsiRhizoDef(i,_) = supplyPsiRhizo(i,_);
    supplyPsiStemDef(i,_) = supplyPsiStem(i,_);
    supplyEStemDef(i,_) = supplyEStem(i,_);
    supplyPsiLeafDef[i] = supplyPsiLeaf[i];
  }
  return(List::create(Named("E") = supplyEDef,
                      Named("ELeaf") = supplyELeafDef,
                      Named("EStem") = supplyEStemDef,
                      Named("ERhizo") = supplyERhizoDef,
                      Named("psiRhizo")=supplyPsiRhizoDef,
                      Named("psiRootCrown")=supplyPsiRootDef,
                      Named("psiStem")=supplyPsiStemDef,
                      Named("psiLeaf")=supplyPsiLeafDef,
                      Named("dEdP")=supplydEdpDef,
                      Named("kterm") = supplyKtermDef));
  
}


//' @rdname hydraulics_supplyfunctions
// [[Rcpp::export("hydraulics_regulatedPsiXylem")]]
NumericVector regulatedPsiXylem(double E, double psiUpstream, double kxylemmax, double c, double d, double psiStep = -0.01) {
  //If Ein > Ecrit then set Ein to Ecrit
  double psiUnregulated = E2psiXylem(E, psiUpstream, kxylemmax, c, d, 0.0);
  double Ec = ECrit(psiUpstream, kxylemmax,c,d);
  double Ein = E;
  if(Ein > Ec) {
    Ein = Ec;
    psiUnregulated = psiCrit(c,d);
  }
  double deltaPsiUnregulated = psiUnregulated - psiUpstream;
  double kp = xylemConductance(psiUpstream, kxylemmax, c, d);
  double deltaPsiRegulated = deltaPsiUnregulated*(xylemConductance(psiUnregulated, kxylemmax, c, d)/kp);
  //replace by maximum if found for lower psi values
  // Rcout <<"Initial "<<psiUnregulated << " "<< deltaPsiRegulated <<"\n";
  for(double psi = psiUpstream; psi > psiUnregulated; psi +=psiStep) {
    double deltaPsi = (psi-psiUpstream)*(xylemConductance(psi, kxylemmax, c,d)/kp);
    // Rcout <<psi << " "<< deltaPsi<< " "<< deltaPsiRegulated <<"\n";
    if(NumericVector::is_na(deltaPsiRegulated)) deltaPsiRegulated = deltaPsi;
    else if(deltaPsi < deltaPsiRegulated) deltaPsiRegulated = deltaPsi;
  }
  double psiRegulated = psiUpstream + deltaPsiRegulated;
  double Efin = EXylem(psiRegulated, psiUpstream, kxylemmax, c, d);
  double relativeConductance1 = Efin/Ein;
  double relativeConductance2 = Efin/E;
  return(NumericVector::create(psiUnregulated, psiRegulated, Ein, Efin, relativeConductance1, relativeConductance2));
}

//' @rdname hydraulics_supplyfunctions
// [[Rcpp::export("hydraulics_regulatedPsiTwoElements")]]
NumericVector regulatedPsiTwoElements(double Emax, double psiSoil, double krhizomax, double kxylemmax, double n, double alpha, double c, double d, double dE = 0.1, double psiMax = -10.0) {
  List s = supplyFunctionTwoElements(Emax, psiSoil, krhizomax, kxylemmax, n, alpha, c, d, 0.0, dE,psiMax);
  NumericVector supplyPsi = s["PsiPlant"];
  NumericVector Efitted = s["FittedE"];
  NumericVector dEdP = s["dEdP"];
  int maxNsteps = Efitted.size();
  double deltaPsiRegulated = 0.0;
  double deltaPsiRegulatedi=0.0;
  double dEdP0 = dEdP[0];
  for(int i=1;i<maxNsteps;i++) {
    if(supplyPsi[i]>psiMax) {
      deltaPsiRegulatedi = (supplyPsi[i] - psiSoil)*std::min(1.0, dEdP[i]/dEdP0);
      // Rcout<<supplydEdp <<" "<<deltaPsiRegulatedi<<"\n";
      if(deltaPsiRegulatedi < deltaPsiRegulated) {
        deltaPsiRegulated = deltaPsiRegulatedi;
      }
    }
  }
  //Regulated potential
  double psiRegulated = psiSoil + deltaPsiRegulated;
  //Find transpiration corresponding to regulated potential
  double ERegulated = 0.0, dEdPRegulated = 0.0;
  for(int i=1;i<maxNsteps;i++) {
    if((supplyPsi[i-1] >= psiRegulated) && (supplyPsi[i]<psiRegulated)) {
      ERegulated = Efitted[i]*std::abs((supplyPsi[i-1]-psiRegulated)/(supplyPsi[i-1]-supplyPsi[i])) + Efitted[i-1]*std::abs((supplyPsi[i]-psiRegulated)/(supplyPsi[i-1]-supplyPsi[i]));
      ERegulated = std::min(ERegulated, Emax);
      psiRegulated = supplyPsi[i]*std::abs((supplyPsi[i-1]-psiRegulated)/(supplyPsi[i-1]-supplyPsi[i])) + supplyPsi[i-1]*std::abs((supplyPsi[i]-psiRegulated)/(supplyPsi[i-1]-supplyPsi[i]));
      dEdPRegulated = dEdP[i]*std::abs((supplyPsi[i-1]-psiRegulated)/(supplyPsi[i-1]-supplyPsi[i])) + dEdP[i-1]*std::abs((supplyPsi[i]-psiRegulated)/(supplyPsi[i-1]-supplyPsi[i]));
      if((supplyPsi[i-1]-supplyPsi[i])==0.0) dEdPRegulated = 0.0;
      // Rcout<<dEdP[i]<< " "<<dEdP[i-1]<< " "<<dEdPRegulated<<"\n";
      break;
    }
  }
  return(NumericVector::create(supplyPsi[maxNsteps-1], psiRegulated, Efitted[maxNsteps-1], ERegulated, dEdPRegulated));
}



//' Scaling from conductivity to conductance
//' 
//' Functions used to scale from tissue conductivity to conductance of different elements of the continuum.
//' 
//' @param psiSoil Soil water potential (in MPa). A scalar or a vector depending on the function.
//' @param psiRhizo Water potential (in MPa) in the rhizosphere (root surface).
//' @param psiStem Water potential (in MPa) in the stem.
//' @param psiLeaf Water potential (in MPa) in the leaf.
//' @param PLCstem Percent loss of conductance (in \%) in the stem.
//' @param L Vector with the length of coarse roots (mm) for each soil layer.
//' @param V Vector with the proportion [0-1] of fine roots within each soil layer.
//' @param krhizomax Maximum rhizosphere hydraulic conductance (defined as flow per leaf surface unit and per pressure drop).
//' @param kleafmax Maximum leaf hydraulic conductance (defined as flow per leaf surface unit and per pressure drop).
//' @param kstemmax Maximum stem xylem hydraulic conductance (defined as flow per leaf surface unit and per pressure drop).
//' @param krootmax Maximum root xylem hydraulic conductance (defined as flow per leaf surface unit and per pressure drop).
//' @param psiStep Water potential precision (in MPa).
//' @param rootc,rootd Parameters of the Weibull function for roots (root xylem vulnerability curve).
//' @param stemc,stemd Parameters of the Weibull function for stems (stem xylem vulnerability curve).
//' @param leafc,leafd Parameters of the Weibull function for leaves (leaf vulnerability curve).
//' @param n,alpha Parameters of the Van Genuchten function (rhizosphere vulnerability curve).
//' @param averageResistancePercent Average (across water potential values) resistance percent of the rhizosphere, with respect to total resistance (rhizosphere + root xylem + stem xylem).
//' @param initialValue Initial value of rhizosphere conductance.
//' @param xylemConductivity Xylem conductivity as flow per length of conduit and pressure drop (in kg·m-1·s-1·MPa-1).
//' @param Al2As Leaf area to sapwood area (in m2·m-2).
//' @param height Plant height (in cm).
//' @param refheight Reference plant height of measurement of xylem conductivity (in cm).
//' @param taper A boolean flag to indicate correction by taper of xylem conduits (Christoffersen et al. 2017).
//' 
//' @details Details of the hydraulic model are given in the medfate book
//' 
//' @return
//' Values returned for each function are:
//' \itemize{
//'   \item{\code{hydraulics_maximumSoilPlantConductance}: The maximum soil-plant conductance, in the same units as the input segment conductances.}
//'   \item{\code{hydraulics_averageRhizosphereResistancePercent}: The average percentage of resistance due to the rhizosphere, calculated across water potential values.}
//'   \item{\code{hydraulics_findRhizosphereMaximumConductance}: The maximum rhizosphere conductance value given an average rhizosphere resistance and the vulnerability curves of rhizosphere, root and stem elements.}
//'   \item{\code{hydraulics_taperFactorSavage}: Taper factor according to Savage et al. (2010).}
//' }
//' 
//' @references
//' Christoffersen, B. O., M. Gloor, S. Fauset, N. M. Fyllas, D. R. Galbraith, T. R. Baker, L. Rowland, R. A. Fisher, O. J. Binks, S. A. Sevanto, C. Xu, S. Jansen, B. Choat, M. Mencuccini, N. G. McDowell, and P. Meir. 2016. Linking hydraulic traits to tropical forest function in a size-structured and trait-driven model (TFS v.1-Hydro). Geoscientific Model Development Discussions 9: 4227–4255.
//' 
//' Savage, V. M., L. P. Bentley, B. J. Enquist, J. S. Sperry, D. D. Smith, P. B. Reich, and E. I. von Allmen. 2010. Hydraulic trade-offs and space filling enable better predictions of vascular structure and function in plants. Proceedings of the National Academy of Sciences of the United States of America 107:22722–7.
//' 
//' Olson, M.E., Anfodillo, T., Rosell, J.A., Petit, G., Crivellaro, A., Isnard, S., \enc{León-Gómez}{Leon-Gomez}, C., \enc{Alvarado-Cárdenas}{Alvarado-Cardenas}, L.O., and Castorena, M. 2014. Universal hydraulics of the flowering plants: Vessel diameter scales with stem length across angiosperm lineages, habits and climates. Ecology Letters 17: 988–997.
//' 
//' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
//' 
//' @seealso
//' \code{\link{hydraulics_psi2K}}, \code{\link{hydraulics_supplyFunctionPlot}}, \code{\link{spwb}}, \code{\link{soil}}
//' 
//' @name hydraulics_scalingconductance
// [[Rcpp::export("hydraulics_maximumSoilPlantConductance")]]
double maximumSoilPlantConductance(NumericVector krhizomax, NumericVector krootmax, 
                                   double kstemmax, double kleafmax) {
  int nlayers = krhizomax.length();
  double krhizo = 0.0;
  double kroot = 0.0;
  for(int i=0;i<nlayers;i++) {
    krhizo = krhizo + krhizomax[i];
    kroot = kroot + krootmax[i];
  }
  double rrhizo = 1.0/krhizo;
  double rroot = 1.0/kroot;
  double rstem = 1.0/kstemmax;
  double rleaf = 1.0/kleafmax;
  return(1.0/(rrhizo+rroot+rstem+rleaf));
}

//' @rdname hydraulics_scalingconductance
// [[Rcpp::export("hydraulics_soilPlantResistances")]]
NumericVector soilPlantResistances(NumericVector psiSoil, NumericVector psiRhizo, 
                                   NumericVector psiStem, NumericVector PLCstem,
                                   double psiLeaf, 
                                   NumericVector krhizomax, NumericVector n, NumericVector alpha,
                                   NumericVector krootmax, double rootc, double rootd, 
                                   double kstemmax, double stemc, double stemd,
                                   double kleafmax, double leafc, double leafd) {
  int nlayers = psiSoil.length();
  double krhizo = 0.0;
  double kroot = 0.0;
  for(int i=0;i<nlayers;i++) {
    krhizo = krhizo + vanGenuchtenConductance(psiSoil[i], krhizomax[i], n[i], alpha[i]);
    kroot = kroot + xylemConductance(psiRhizo[i], krootmax[i], rootc, rootd);
  }
  double rrhizo = 1.0/krhizo;
  double rroot = 1.0/kroot;
  int nStemSegments = psiStem.length();
  double kxsegmax = kstemmax*((double) nStemSegments);
  double rstem = 0.0;
  double plcCond = NA_REAL;
  for(int i=0;i<nStemSegments;i++) {
    plcCond = (1.0-PLCstem[i]);
    rstem = rstem + 1.0/(kxsegmax*std::min(plcCond, xylemConductance(psiStem[i], 1.0, stemc, stemd)));
  }
  double rleaf = 1.0/xylemConductance(psiLeaf, kleafmax, leafc, leafd);
  NumericVector resistances = NumericVector::create(rrhizo, rroot, rstem, rleaf);
  return(resistances);
}

/*
 * Parametrization of rhizosphere conductance
 */
double rhizosphereResistancePercent(double psiSoil, 
                                    double krhizomax, double n, double alpha,
                                    double krootmax, double rootc, double rootd,
                                    double kstemmax, double stemc, double stemd,
                                    double kleafmax, double leafc, double leafd) {
  double krhizo = vanGenuchtenConductance(psiSoil, krhizomax, n, alpha);
  double kroot = xylemConductance(psiSoil, krootmax, rootc, rootd);
  double kstem = xylemConductance(psiSoil, kstemmax, stemc, stemd);
  double kleaf = xylemConductance(psiSoil, kleafmax, leafc, leafd);
  return(100.0*(1.0/krhizo)/((1.0/kroot)+(1.0/kstem)+(1.0/kleaf)+(1.0/krhizo)));
}

//' @rdname hydraulics_scalingconductance
// [[Rcpp::export("hydraulics_averageRhizosphereResistancePercent")]]
double averageRhizosphereResistancePercent(double krhizomax, double n, double alpha,
                                           double krootmax, double rootc, double rootd,
                                           double kstemmax, double stemc, double stemd, 
                                           double kleafmax, double leafc, double leafd,
                                           double psiStep = -0.01){
  double psiC = psiCrit(stemc, stemd);
  double cnt = 0.0;
  double sum = 0.0;
  for(double psi=0.0; psi>psiC;psi += psiStep) {
    sum +=rhizosphereResistancePercent(psi, krhizomax, n,alpha,krootmax, rootc, rootd,
                                       kstemmax, stemc,stemd,
                                       kleafmax, leafc,leafd);
    cnt+=1.0;
  }
  return(sum/cnt);
}




//' @rdname hydraulics_scalingconductance
// [[Rcpp::export("hydraulics_findRhizosphereMaximumConductance")]]
double findRhizosphereMaximumConductance(double averageResistancePercent, double n, double alpha,
                                         double krootmax, double rootc, double rootd,
                                         double kstemmax, double stemc, double stemd,
                                         double kleafmax, double leafc, double leafd,
                                         double initialValue = 0.0) {
  double step = 1.0;
  double fTol = 0.1;
  // Rcout<<exp(initialValue)<<"\n";
  double krhizomaxlog = initialValue;
  int nsteps = 0;
  int max_nsteps = 100;
  double f = averageRhizosphereResistancePercent(exp(krhizomaxlog), n,alpha,krootmax, rootc, rootd,
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
    f = averageRhizosphereResistancePercent(exp(krhizomaxlog), n,alpha,krootmax, rootc, rootd,
                                            kstemmax, stemc,stemd,
                                            kleafmax, leafc,leafd);
    nsteps++;
  }
  // Rcout<< nsteps<< " "<<exp(krhizomaxlog) << " "<< f << " "<< averageResistancePercent<< " "<< step<<"\n";
  return(exp(krhizomaxlog));
}



/**
 * BIOMECHANICS
 * 
 * Savage, V. M., L. P. Bentley, B. J. Enquist, J. S. Sperry, D. D. Smith, P. B. Reich, and E. I. von Allmen. 2010. Hydraulic trade-offs and space filling enable better predictions of vascular structure and function in plants. Proceedings of the National Academy of Sciences of the United States of America 107:22722–7.
 * Christoffersen, B. O., M. Gloor, S. Fauset, N. M. Fyllas, D. R. Galbraith, T. R. Baker, L. Rowland, R. A. Fisher, O. J. Binks, S. A. Sevanto, C. Xu, S. Jansen, B. Choat, M. Mencuccini, N. G. McDowell, and P. Meir. 2016. Linking hydraulic traits to tropical forest function in a size-structured and trait-driven model (TFS v.1-Hydro). Geoscientific Model Development Discussions 0:1–60.
 * 
 * height - Tree height in cm
 */
//' @rdname hydraulics_scalingconductance
// [[Rcpp::export("hydraulics_taperFactorSavage")]]
double taperFactorSavage(double height) {
  double b_p0 = 1.32, b_p13 = 1.85; //normalizing constants (p = 1/3)
  double a_p0 = 7.20E-13, a_p13 = 6.67E-13;
  double n_ext = 2.0; //Number of daughter branches per parent
  double N = ((3.0*log(1.0-(height/4.0)*(1.0-pow(n_ext, 1.0/3.0))))/log(n_ext))-1.0;
  double K_0 = a_p0*pow(pow(n_ext, N/2.0),b_p0);
  double K_13 = a_p13*pow(pow(n_ext, N/2.0),b_p13);
  return(K_13/K_0);
}

/**
 *  Returns the terminal conduit radius (in micras)
 *  
 *  height - plant height in cm
 */
//' @rdname hydraulics_scalingconductance
// [[Rcpp::export("hydraulics_terminalConduitRadius")]]
double terminalConduitRadius(double height) {
  double dh  = pow(10,1.257 +  0.24*log10(height/100.0));//Olson, M.E., Anfodillo, T., Rosell, J.A., Petit, G., Crivellaro, A., Isnard, S., León-Gómez, C., Alvarado-Cárdenas, L.O., & Castorena, M. 2014. Universal hydraulics of the flowering plants: Vessel diameter scales with stem length across angiosperm lineages, habits and climates. Ecology Letters 17: 988–997.
  return(dh/2.0);
}


//' @rdname hydraulics_scalingconductance
// [[Rcpp::export("hydraulics_referenceConductivityHeightFactor")]]
double referenceConductivityHeightFactor(double refheight, double height) {
  double rhref  = terminalConduitRadius(refheight);
  double rh  = terminalConduitRadius(height);
  double df = pow(rh/rhref,2.0);
  return(df);
}


/**
 * Calculate maximum leaf-specific stem hydraulic conductance (in mmol·m-2·s-1·MPa-1)
 * 
 * xylemConductivity - Sapwood-specific conductivity of stem xylem (in kg·m-1·s-1·MPa-1), 
 *                     assumed to be measured at distal twigs
 * refheight - Reference plant height (on which xylem conductivity was measured)
 * Al2As - Leaf area to sapwood area ratio (in m2·m-2)
 * height - plant height (in cm)
 * taper - boolean to apply taper
 */
//' @rdname hydraulics_scalingconductance
// [[Rcpp::export("hydraulics_maximumStemHydraulicConductance")]]
double maximumStemHydraulicConductance(double xylemConductivity, double refheight, double Al2As, double height, 
                                       bool taper = false) {
  
  
  // Christoffersen, B. O., M. Gloor, S. Fauset, N. M. Fyllas, D. R. Galbraith, T. R. Baker, L. Rowland, R. A. Fisher, O. J. Binks, S. A. Sevanto, C. Xu, S. Jansen, B. Choat, M. Mencuccini, N. G. McDowell, and P. Meir. 2016. Linking hydraulic traits to tropical forest function in a size-structured and trait-driven model (TFS v.1-Hydro). Geoscientific Model Development Discussions 0:1–60.
  double kmax = 0.0;
  if(!taper) {
    double xylemConductivityCorrected = xylemConductivity*referenceConductivityHeightFactor(refheight, height);
    kmax =   (1000.0/0.018)*(xylemConductivityCorrected/Al2As)*(100.0/height);
  } else {
    double petioleConductivity = xylemConductivity*referenceConductivityHeightFactor(refheight, 100.0);
    // Correct reference conductivity in relation to the reference plant height in which it was measured
    kmax =   (1000.0/0.018)*(petioleConductivity/Al2As)*(100.0/height)*(taperFactorSavage(height)/(taperFactorSavage(100.0)));
  } 
  return(kmax); 
}

/**
 * Proportions of root xylem conductance
 * 
 * Calculates the proportion of total xylem conductance that corresponds to each layer in a network of 
 * parallel xylem resistances.
 * 
 * Sperry, J. S., Y. Wang, B. T. Wolfe, D. S. Mackay, W. R. L. Anderegg, N. G. Mcdowell, and W. T. Pockman. 2016. 
 * Pragmatic hydraulic theory predicts stomatal responses to climatic water deficits. 
 * New Phytologist 212:577–589.
 * 
 */
//' @rdname hydraulics_scalingconductance
// [[Rcpp::export("hydraulics_rootxylemConductanceProportions")]]
NumericVector rootxylemConductanceProportions(NumericVector L, NumericVector V) {
  int nlayers = L.size();
  //Weights
  NumericVector w(nlayers, 0.0);
  double wsum=0.0;
  for(int i=0;i<nlayers;i++) {
    if(L[i]>0.0) {
      w[i]= V[i]*(1.0/L[i]);
      wsum +=w[i];
    }
  }
  for(int i=0;i<nlayers;i++) w[i] = w[i]/wsum;
  return(w);
}

/**
 * Calculate maximum leaf-specific root hydraulic conductance (in mmol·m-2·s-1·MPa-1)
 * 
 * xylemConductivity - Sapwood-specific conductivity of root xylem (in kg·m-1·s-1·MPa-1)
 * Al2As - Leaf area to sapwood area ratio (in m2·m-2)
 * v - proportion of fine roots in each soil layer
 * widths - soil layer depths (in mm)
 */
// double maximumRootHydraulicConductance(double xylemConductivity, double Al2As, NumericVector v, 
//                                        NumericVector widths, double depthWidthRatio = 1.0){
//   NumericVector rl = coarseRootLengths(v,widths, depthWidthRatio);
//   NumericVector w = xylemConductanceProportions(v,widths, depthWidthRatio);
//   int nlayers = v.length();
//   double kmax = 0.0;
//   for(int i=0;i<nlayers;i++) {
//     if(rl[i]>0.0) kmax = kmax + w[i]*(1000.0/0.018)*(xylemConductivity/((rl[i]/1000.0)*Al2As)); 
//   }
//   return(kmax); 
// }



