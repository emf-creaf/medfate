#include <Rcpp.h>
#include <numeric>
#include <math.h>
#include "carbon.h"
#include "forestutils.h"
#include "paramutils.h"
using namespace Rcpp;


/**
 * Calculate water capacity of stem per leaf area (in mm = l·m-2)
 * 
 * Al2As - Leaf area to sapwood area ratio (in m2·m-2)
 * height - plant height (in cm)
 * wd - wood density (in g·cm-3)
 * http://www.fao.org/forestry/17109/en/
 */
// [[Rcpp::export("moisture_sapwoodWaterCapacity")]]
double sapwoodWaterCapacity(double Al2As, double height, NumericVector V, NumericVector L, double wd) {
  int nlayers = V.size();
  double woodPorosity = (1.0- (wd/1.54));
  double vAbove = 1000*(height/(Al2As*100.0))*woodPorosity;
  double vBelow = 0.0;
  for(int i=0;i<nlayers;i++) {
    vBelow += 1000*(V[i]*(L[i]/10.0)/(Al2As*100.0))*woodPorosity;
  }
  return(vAbove+vBelow);
  // 
  // return(1000*(height/(Al2As*100.0))*(1.0- (wd/1.54)));
}

/**
 * Calculate  capacity of leaves per leaf area (in mm = l·m-2)
 * 
 * SLA - Specific leaf area (in m2/kg)
 * ld - leaf density (in g/cm3 = 1000 kg/m3)
 */
// [[Rcpp::export("moisture_leafWaterCapacity")]]
double leafWaterCapacity(double SLA, double ld) {
  return(1000.0/(1000.0*ld*SLA))*(1.0- (ld/1.54));
}


// [[Rcpp::export("moisture_turgorLossPoint")]]
double turgorLossPoint(double pi0, double epsilon) {
  return((pi0*epsilon)/(pi0+epsilon));
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
  double psi_tl = turgorLossPoint(pi0,epsilon);
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
