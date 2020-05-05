#include <Rcpp.h>
#include <numeric>
#include <math.h>
#include "carbon.h"
using namespace Rcpp;


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

// [[Rcpp::export("moisture_tissueFMC")]]
double tissueFMC(double RWC, double density, double d0 = 1.54) {
  return(100*RWC*((1.0/density) - (1.0/d0)));
}



/**
 * Translates soil water balance results to fuel moisture content of plant cohorts.
 * 
 *   spwb - The output of spwb() or growth()
 */
// [[Rcpp::export("moisture_cohortFMC")]]
List cohortFMC(List spwb) {
  List x;
  if(spwb.containsElementNamed("spwbInput")) x = spwb["spwbInput"];
  else if(spwb.containsElementNamed("growthInput")) x = spwb["growthInput"];
  else stop("wrong input");
  
  //Draw cohort-based variables
  DataFrame cohorts = Rcpp::as<Rcpp::DataFrame>(x["cohorts"]);
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  DataFrame paramsAnatomy = Rcpp::as<Rcpp::DataFrame>(x["paramsAnatomy"]);
  DataFrame paramsTranspiration = Rcpp::as<Rcpp::DataFrame>(x["paramsTranspiration"]);
  
  //Anatomy parameters
  NumericVector r635 = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["r635"]);
  NumericVector WoodDens = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["WoodDensity"]);
  NumericVector LeafDens = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["LeafDensity"]);
  
  //Transpiration parameters
  NumericVector VCstem_c = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["VCstem_c"]);
  NumericVector VCstem_d = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["VCstem_d"]);
  NumericVector VCleaf_c = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["VCleaf_c"]);
  NumericVector VCleaf_d = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["VCleaf_d"]);
  
  //Water storage parameters
  DataFrame paramsWaterStorage = Rcpp::as<Rcpp::DataFrame>(x["paramsWaterStorage"]);
  NumericVector StemAF = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["StemAF"]);
  NumericVector LeafAF = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["LeafAF"]);
  
  List plants = spwb["Plants"];
  NumericMatrix psiapoleaf = Rcpp::as<Rcpp::NumericMatrix>(plants["LeafPsiMin"]);
  NumericMatrix StemPLC = Rcpp::as<Rcpp::NumericMatrix>(plants["StemPLC"]);
  NumericMatrix RWCsymleaf = Rcpp::as<Rcpp::NumericMatrix>(plants["LeafRWC"]);
  NumericMatrix RWCsymstem = Rcpp::as<Rcpp::NumericMatrix>(plants["StemRWC"]);
  List l = psiapoleaf.attr("dimnames");
  CharacterVector days = l[0];
  CharacterVector cohNames = l[1];
  int numDays = psiapoleaf.nrow();
  int numCohorts = psiapoleaf.ncol();
  
  NumericMatrix leafFMC(numDays, numCohorts);
  NumericMatrix twigFMC(numDays, numCohorts);
  NumericMatrix fineFMC(numDays, numCohorts);
  leafFMC.attr("dimnames") = l;
  twigFMC.attr("dimnames") = l;
  fineFMC.attr("dimnames") = l;
  for(int c=0;c<numCohorts;c++) {
    double f_apo_leaf = LeafAF[c];
    double f_apo_stem = StemAF[c];
    double density_leaf = LeafDens[c];
    double density_stem = WoodDens[c];
    double leafc = VCleaf_c[c];
    double leafd = VCleaf_d[c];
    double p_leaves = 1.0/r635[c];
    for(int d=0;d<numDays;d++) {
      double rwc_apo_leaf = apoplasticRelativeWaterContent(psiapoleaf(d,c), leafc, leafd);
      double rwc_sym_leaf = RWCsymleaf(d,c);
      double rwc_leaf = rwc_apo_leaf*f_apo_leaf + rwc_sym_leaf*(1.0 - f_apo_leaf);
      leafFMC(d,c) = tissueFMC(rwc_leaf, density_leaf);
      double rwc_apo_stem = 1.0-StemPLC(d,c);
      double rwc_sym_stem = RWCsymstem(d,c);
      double rwc_stem = rwc_apo_stem*f_apo_stem + rwc_sym_stem*(1.0 - f_apo_stem);
      twigFMC(d,c) = tissueFMC(rwc_stem, density_stem);
      
      fineFMC(d,c)  =leafFMC(d,c)*p_leaves + twigFMC(d,c)*(1.0 - p_leaves);
    }
  }
  return(List::create(Named("LeafFMC") = leafFMC, Named("TwigFMC") = twigFMC, Named("FineFMC") = fineFMC));
}

  