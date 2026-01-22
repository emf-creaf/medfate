#include <numeric>
#include <math.h>


/**
 * Calculate  capacity of leaves per leaf area (in mm = l·m-2)
 * 
 * SLA - Specific leaf area (in m2/kg)
 * ld - leaf density (in g/cm3 = 1000 kg/m3)
 */
//' @rdname moisture
//' @keywords internal
// [[Rcpp::export("moisture_leafWaterCapacity")]]
double leafWaterCapacity_c(double SLA, double ld) {
  return(1000.0/(1000.0*ld*SLA))*(1.0- (ld/1.54));
}


//' @rdname moisture
//' @keywords internal
// [[Rcpp::export("moisture_turgorLossPoint")]]
double turgorLossPoint_c(double pi0, double epsilon) {
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
//' @rdname moisture
//' @keywords internal
// [[Rcpp::export("moisture_symplasticRWC")]]
double symplasticRelativeWaterContent_c(double psiSym, double pi0, double epsilon) {
  double psi_tl = turgorLossPoint_c(pi0,epsilon);
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

//' @rdname moisture
//' @keywords internal
// [[Rcpp::export("moisture_symplasticPsi")]]
double symplasticWaterPotential_c(double RWC, double pi0, double epsilon) {
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
//' @rdname moisture
//' @keywords internal
// [[Rcpp::export("moisture_apoplasticRWC")]]
double apoplasticRelativeWaterContent_c(double psiApo, double c, double d) {
  if(psiApo>=0.0) return(1.0);
  return(exp(-pow(psiApo/d,c)));
}

//' @rdname moisture
//' @keywords internal
// [[Rcpp::export("moisture_apoplasticPsi")]]
double apoplasticWaterPotential_c(double RWC, double c, double d) {
  double psi = d*pow(-1.0*log(RWC),1.0/c);
  if( psi< -40.0) psi = -40.0; //Minimum value
  return(psi);
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
//' @rdname moisture
//' @keywords internal
// [[Rcpp::export("moisture_tissueRWC")]]
double tissueRelativeWaterContent_c(double psiSym, double pi0, double epsilon, 
                                  double psiApo, double c, double d, 
                                  double af) {
  double sym_rwc = symplasticRelativeWaterContent_c(psiSym, pi0, epsilon);
  double apo_rwc = apoplasticRelativeWaterContent_c(psiApo, c, d); //Water content cannot be higher than the fraction of non-embolized conduits
  return(sym_rwc*(1.0-af)+apo_rwc*af);
}
