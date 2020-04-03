#include <Rcpp.h>
#include <numeric>
#include <math.h>
using namespace Rcpp;

double Rn = 0.008314472; // The perfect gas constant MPa·l/K·mol = kJ/K·mol
double waterViscosity;
double sucroseMolarWeight = 342.3; //g*mol-1

/**
 * On the pressure dependence of the viscosity of aqueous sugar solutions
 * Rheol Acta (2002) 41: 369–374 DOI 10.1007/s00397-002-0238-y
 * 
 *  x - sugar concentration (mol/l)
 *  temp - temperature (degrees C)
 */
// [[Rcpp::export("moisture_relativeSapViscosity")]]
double relativeSapViscosity(double conc, double temp) {
  double x = conc*sucroseMolarWeight/1e3; //from mol/l to g*cm-3
  double Tkelvin = temp + 273.15;
  double q0a = 1.12; //g*cm-3
  double q1 = -0.248;
  double Ea = 2.61; //kJ*mol-1 energy of activation
  double va = x/(q0a*exp(-1.0*Ea/(Rn*Tkelvin)));
  double relVisc = exp(va/(1.0 + q1*va)); // relative viscosity
  double relWat = exp(-3.7188+(578.919/(-137.546+ Tkelvin))); // Vogel equation for liquid dynamic viscosity (= 1 for 25ºC)
  return(relWat*relVisc);
}
/**
 * Van 't Hoff equation
 *  conc - mol/l 
 *  temp - deg C
 *  wp - MPa
 */
// [[Rcpp::export("moisture_osmoticWaterPotential")]]
double osmoticWaterPotential(double conc, double temp) {
  return(- conc*Rn*(temp + 273.15));
}
// [[Rcpp::export("moisture_sugarConcentration")]]
double sugarConcentration(double osmoticWP, double temp) {
  return(- osmoticWP/(Rn*(temp + 273.15)));
}

/**
 *  Turgor (MPa)
 *  conc - mol/l 
 *  temp - deg C
 *  psi - water potential (MPa)
 */
// [[Rcpp::export("moisture_turgor")]]
double turgor(double psi, double conc, double temp) {
  return(std::max(0.0, psi-osmoticWaterPotential(conc,temp)));
}

/**
 * floem flow (Holtta et al. 2017)
 *  psiUpstream, psiDownstream - water potential upstream (leaves)  and downstream
 *  concUpstream, concDownstream - sugar concentration upstream (leaves) and downstream (stem)
 *  k_f - floem conductance per leaf area basis (l*m-2*MPa-1*s-1)
 *  
 *  out mol*s-1
 */
// [[Rcpp::export("moisture_floemFlow")]]
double floemFlow(double psiUpstream, double psiDownstream,
                 double concUpstream, double concDownstream,
                 double temp, double k_f = 3.0e-5) {
  double turgor_up = turgor(psiUpstream, concUpstream, temp);
  double turgor_down = turgor(psiDownstream, concDownstream, temp);
  double concMean = (concUpstream+concDownstream)/2.0;
  double relVisc = relativeSapViscosity(concMean, temp);
  if(temp < 0.0) k_f = 0.0; // No floem flow if temperature below zero
  return(k_f*concMean*(turgor_up - turgor_down)/relVisc);
}
/**
* Calculates symplastic relative water content from tissue water potential
* 
*  Bartlett, M. K., C. Scoffoni, and L. Sack. 2012. 
*  The determinants of leaf turgor loss point and prediction of drought tolerance of species and biomes: a global meta-analysis. 
*  Ecology letters 15:393–405.
*  
*  psi - Water potential (MPa)
*  pi0 - Full turgor osmotic potential (MPa)
*  epsilon - Bulk modulus elasticity (MPa)
*  
*  Returns symplastic RWC as proportion of maximum hydration 
*/
// [[Rcpp::export("moisture_symplasticRWC")]]
double symplasticRelativeWaterContent(double psiSym, double pi0, double epsilon) {
  double psi_tl = (pi0*epsilon)/(pi0+epsilon);
  double rwc = 0;
  if(psiSym< psi_tl) {
    rwc = (-std::abs(pi0))/psiSym;
  } else {
    double c = std::abs(pi0);
    double b = psiSym+epsilon - c;
    double a = -epsilon;
    rwc = ((-b)-sqrt(pow(b,2.0)-4.0*a*c))/(2.0*a);
  }
  return(rwc);
}

// [[Rcpp::export("moisture_symplasticPsi")]]
double symplasticWaterPotential(double RWC, double pi0, double epsilon) {
  double pt = std::max(0.0,-pi0 +(-1.0)*epsilon*(1.0 - RWC));
  double ps = pi0/RWC;
  return(pt+ps);
}
/**
 * Calculates apoplastic relative water content from water potential
 * 
 *  psi - Leaf water potential (MPa)
 *  c, d - Parameters of the vulnerability curve
 *  
 *  Returns Apoplastic RWC as proportion of maximum hydration 
 */
// [[Rcpp::export("moisture_apoplasticRWC")]]
double apoplasticRelativeWaterContent(double psiApo, double c, double d) {
  if(psiApo>=0.0) return(1.0);
  return(exp(-pow(psiApo/d,c)));
}

// [[Rcpp::export("moisture_apoplasticPsi")]]
double apoplasticWaterPotential(double RWC, double c, double d) {
  return(d*pow(-1.0*log(RWC),1.0/c));
}
/**
 * Calculates leaf relative water content from leaf water potential
 * 
 *  Bartlett, M. K., C. Scoffoni, and L. Sack. 2012. 
 *  The determinants of leaf turgor loss point and prediction of drought tolerance of species and biomes: a global meta-analysis. 
 *  Ecology letters 15:393–405.
 *  
 *  psiApo - Apoplastic water potential (MPa)
 *  pi0 - Full turgor osmotic potential (MPa)
 *  epsilon - bulk modulus elasticity (MPa)
 *  af - Apoplastic fraction (proportion)
 *  femb - Fraction of embolized conduits
 *  
 *  Returns tissue RWC as proportion of maximum hydration (= g H2O / g H2O sat = m3 H2O / m3 H2O sat)
 */
// [[Rcpp::export("moisture_tissueRWC")]]
double tissueRelativeWaterContent(double psiSym, double pi0, double epsilon, 
                                  double psiApo, double c, double d, 
                                  double af, double femb = 0.0) {
  double sym_rwc = symplasticRelativeWaterContent(psiSym, pi0, epsilon);
  double apo_rwc = std::min(1.0-femb, apoplasticRelativeWaterContent(psiApo, c, d)); //Water content cannot be higher than the fraction of non-embolized conduits
  return(sym_rwc*(1.0-af)+apo_rwc*af);
}

// [[Rcpp::export("moisture_tissueFMC")]]
double tissueFMC(double RWC, double density, double d0 = 1.54) {
  return(100*RWC*((1.0/density) - (1.0/d0)));
}
  