#include <Rcpp.h>
#include <numeric>
#include <math.h>
using namespace Rcpp;


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
// [[Rcpp::export("moisture.symplasticRWC")]]
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

// [[Rcpp::export("moisture.symplasticPsi")]]
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
// [[Rcpp::export("moisture.apoplasticRWC")]]
double apoplasticRelativeWaterContent(double psiApo, double c, double d) {
  if(psiApo>=0.0) return(1.0);
  return(exp(-pow(psiApo/d,c)));
}

// [[Rcpp::export("moisture.apoplasticPsi")]]
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
 *  
 *  Returns tissue RWC as proportion of maximum hydration (= g H2O / g H2O sat = m3 H2O / m3 H2O sat)
 */
// [[Rcpp::export("moisture.tissueRWC")]]
double tissueRelativeWaterContent(double psiSym, double pi0, double epsilon, 
                                  double psiApo, double c, double d, 
                                  double af) {
  double sym_rwc = symplasticRelativeWaterContent(psiSym, pi0, epsilon);
  double apo_rwc = apoplasticRelativeWaterContent(psiApo, c, d);
  return(sym_rwc*(1.0-af)+apo_rwc*af);
}

// [[Rcpp::export("moisture.tissueFMC")]]
double tissueFMC(double RWC, double density, double d0 = 1.54) {
  return(100*RWC*((1.0/density) - (1.0/d0)));
}
  