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
double symplasticRelativeWaterContent(double psi, double pi0, double epsilon) {
  double psi_tl = (pi0*epsilon)/(pi0+epsilon);
  double rwc = 0;
  if(psi< psi_tl) {
    rwc = (-std::abs(pi0))/psi;
  } else {
    double c = std::abs(pi0);
    double b = psi+epsilon - c;
    double a = -epsilon;
    rwc = ((-b)-sqrt(pow(b,2.0)-4.0*a*c))/(2.0*a);
  }
  return(rwc);
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
double apoplasticRelativeWaterContent(double psi, double c, double d, double cellWallFraction = 0.07) {
  if(psi>=0.0) return(1.0);
  return(cellWallFraction+(exp(-pow(psi/d,c))*(1.0 - cellWallFraction)));
}

/**
 * Calculates leaf relative water content from leaf water potential
 * 
 *  Bartlett, M. K., C. Scoffoni, and L. Sack. 2012. 
 *  The determinants of leaf turgor loss point and prediction of drought tolerance of species and biomes: a global meta-analysis. 
 *  Ecology letters 15:393–405.
 *  
 *  psi - Leaf water potential (MPa)
 *  pi0 - Full turgor osmotic potential (MPa)
 *  epsilon - bulk modulus elasticity (MPa)
 *  af - Apoplastic fraction (proportion)
 *  
 *  Returns Leaf RWC as proportion of maximum hydration (= g H2O / g H2O sat = m3 H2O / m3 H2O sat)
 */
// [[Rcpp::export("moisture.leafRWC")]]
double leafRelativeWaterContent(double psi, double pi0, double epsilon, double af) {
  return((1.0-af)*symplasticRelativeWaterContent(psi,pi0,epsilon)+af);
}

/**
 * Calculates branch relative water content from leaf water potential
 * Assumes there is no heartwood (all tissue is either sapwood or living tissue)
 * 
 * Christoffersen, B.O., Gloor, M., Fauset, S., Fyllas, N.M., Galbraith, D.R., Baker, T.R., Rowland, L., Fisher, R.A., Binks, O.J., Sevanto, S.A., Xu, C., Jansen, S., Choat, B., Mencuccini, M., McDowell, N.G., & Meir, P. 2016. Linking hydraulic traits to tropical forest function in a size-structured and trait-driven model (TFS v.1-Hydro). Geoscientific Model Development Discussions 0: 1–60.
 */
// [[Rcpp::export("moisture.branchRWC")]]
double branchRelativeWaterContent(double psi, double wd, double c, double d, double af = 0.80) {
  double pi0_sap = 0.52 - 4.16*wd;
  double eps_sap = sqrt(1.02*exp(8.5*wd)-2.89);
  double sym_rwc = symplasticRelativeWaterContent(psi, pi0_sap, eps_sap);
  double apo_rwc = apoplasticRelativeWaterContent(psi, c, d);
  return(sym_rwc*(1.0-af)+apo_rwc*af);
}
  
/**
 * Calculates fine fuel (i.e. leaves + small branches) moisture content content from water potential
 * 
 *  psi - water potential (MPa)
 *  leaf_pi0 - Leaf full turgor osmotic potential (MPa)
 *  leaf_eps - Leaf bulk modulus elasticity (MPa)
 *  leaf_af - Leaf apoplastic fraction (proportion)
 *  wd - wood density (g/cm3)
 *  c,d - parameters of the vulnerability curve
 *  r635 - Ratio of biomass of foliar + small branches (<6.35) to foliar biomass
 *  
 *  Returns Fine fuel moisture content as percentage of dry weight (= g H2O / g dry weight)
 */
// [[Rcpp::export("moisture.fineFuelRWC")]]
double fineFuelRelativeWaterContent(double psi, double leaf_pi0, double leaf_eps, double leaf_af, 
                               double wd, double c, double d, double r635) {
  double leafRWC = leafRelativeWaterContent(psi, leaf_pi0, leaf_eps, leaf_af);
  double branchRWC = branchRelativeWaterContent(psi, wd, c, d, 0.8);
  double leafProp = (1.0/r635);
  return(leafRWC*leafProp + branchRWC*(1.0-leafProp));
}



