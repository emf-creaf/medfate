#include <Rcpp.h>
#include <numeric>
#include <math.h>
#include "carbon.h"
#include "forestutils.h"
#include "paramutils.h"
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



/**
 * Translates soil water balance results to fuel moisture content of plant cohorts.
 * 
 *   spwb - The output of spwb() or growth()
 *   SpParams - Data frame of (additional) species parameters (specially FMCmax)
 */
// [[Rcpp::export("moisture_cohortFMC")]]
NumericMatrix cohortFMC(List spwb, DataFrame SpParams) {
  List x;
  if(spwb.containsElementNamed("spwbInput")) x = spwb["spwbInput"];
  else if(spwb.containsElementNamed("growthInput")) x = spwb["growthInput"];
  else stop("wrong input");

  List plants = spwb["Plants"];
  NumericMatrix StemPLC = Rcpp::as<Rcpp::NumericMatrix>(plants["StemPLC"]);
  List l = StemPLC.attr("dimnames");
  StringVector days = l[0];
  StringVector cohNames = l[1];
  int numDays = StemPLC.nrow();
  int numCohorts = StemPLC.ncol();
  
  //Draw cohort-based variables
  DataFrame cohorts = Rcpp::as<Rcpp::DataFrame>(x["cohorts"]);
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  List control = x["control"];
  String transpirationMode = control["transpirationMode"];
  
  //FMCmax
  IntegerVector SP = cohorts["SP"];
  NumericVector FMCmax = speciesNumericParameter(SP, SpParams, "maxFMC");
  
  if(transpirationMode == "Granier") {
    NumericVector LeafPI0 = speciesNumericParameter(SP, SpParams, "LeafPI0");
    NumericVector LeafEPS = speciesNumericParameter(SP, SpParams, "LeafEPS");
    NumericVector LeafAF = speciesNumericParameter(SP, SpParams, "LeafAF");
    NumericMatrix PlantPsi = Rcpp::as<Rcpp::NumericMatrix>(plants["PlantPsi"]);

    NumericMatrix leafFMC(numDays, numCohorts);
    leafFMC.attr("dimnames") = l;
    //Average values for Mediterranean climate species
    for(int c=0;c<numCohorts;c++) {
      if(NumericVector::is_na(LeafPI0[c])) {
        String s = cohNames[c];
        LeafPI0[c] = -2.0;  
        warning("Default LeafPI0 = %f for cohort %s", LeafPI0[c],s.get_cstring());
      }
      if(NumericVector::is_na(LeafEPS[c])) {
        String s = cohNames[c];
        LeafEPS[c] = 17.0; 
        warning("Default LeafEPS = %f for cohort %s",LeafEPS[c],s.get_cstring());
      }
      if(NumericVector::is_na(LeafAF[c])) {
        String s = cohNames[c];
        LeafAF[c] = 0.29; 
        warning("Default LeafAF = %f for cohort %s",LeafAF[c],s.get_cstring());
      }
      for(int d=0;d<numDays;d++) {
        double rwc = symplasticRelativeWaterContent(PlantPsi(d,c), LeafPI0[c], LeafEPS[c]);
        leafFMC(d,c) = (rwc*(1.0-LeafAF[c])+ (1.0 - StemPLC(d,c))*LeafAF[c])*FMCmax[c];
      }
    }
    return(leafFMC);
  } else {
    DataFrame paramsAnatomy = Rcpp::as<Rcpp::DataFrame>(x["paramsAnatomy"]);
    DataFrame paramsTranspiration = Rcpp::as<Rcpp::DataFrame>(x["paramsTranspiration"]);
    //Anatomy parameters
    // NumericVector r635 = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["r635"]);
    // NumericVector WoodDens = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["WoodDensity"]);
    // NumericVector LeafDens = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["LeafDensity"]);
    //Transpiration parameters
    // NumericVector VCstem_c = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["VCstem_c"]);
    // NumericVector VCstem_d = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["VCstem_d"]);
    NumericVector VCleaf_c = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["VCleaf_c"]);
    NumericVector VCleaf_d = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["VCleaf_d"]);
    
    //Water storage parameters
    DataFrame paramsWaterStorage = Rcpp::as<Rcpp::DataFrame>(x["paramsWaterStorage"]);
    // NumericVector StemAF = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["StemAF"]);
    NumericVector LeafAF = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["LeafAF"]);
    
    NumericMatrix psiapoleaf = Rcpp::as<Rcpp::NumericMatrix>(plants["LeafPsiMin"]);
    NumericMatrix RWCsymleaf = Rcpp::as<Rcpp::NumericMatrix>(plants["LeafSympRWC"]);
    // NumericMatrix RWCsymstem = Rcpp::as<Rcpp::NumericMatrix>(plants["StemSympRWC"]);
    
    NumericMatrix leafFMC(numDays, numCohorts);
    // NumericMatrix twigFMC(numDays, numCohorts);
    // NumericMatrix fineFMC(numDays, numCohorts);
    leafFMC.attr("dimnames") = l;
    // twigFMC.attr("dimnames") = l;
    // fineFMC.attr("dimnames") = l;
    for(int c=0;c<numCohorts;c++) {
      double f_apo_leaf = LeafAF[c];
      // double f_apo_stem = StemAF[c];
      // double density_leaf = LeafDens[c];
      // double density_stem = WoodDens[c];
      double leafc = VCleaf_c[c];
      double leafd = VCleaf_d[c];
      // double p_leaves = 1.0/r635[c];
      for(int d=0;d<numDays;d++) {
        double rwc_apo_leaf = apoplasticRelativeWaterContent(psiapoleaf(d,c), leafc, leafd);
        double rwc_sym_leaf = RWCsymleaf(d,c);
        double rwc_leaf = rwc_apo_leaf*f_apo_leaf + rwc_sym_leaf*(1.0 - f_apo_leaf);
        leafFMC(d,c) = rwc_leaf*FMCmax[c];
        // double rwc_apo_stem = 1.0-StemPLC(d,c);
        // double rwc_sym_stem = RWCsymstem(d,c);
        // double rwc_stem = rwc_apo_stem*f_apo_stem + rwc_sym_stem*(1.0 - f_apo_stem);
        // twigFMC(d,c) = tissueFMC(rwc_stem, density_stem);
        // 
        // fineFMC(d,c)  =leafFMC(d,c)*p_leaves + twigFMC(d,c)*(1.0 - p_leaves);
      }
    }
    // return(List::create(Named("LeafFMC") = leafFMC, Named("TwigFMC") = twigFMC, Named("FineFMC") = fineFMC));
    return(leafFMC);
  }
}


// [[Rcpp::export("moisture_cohortFMCDay")]]
NumericVector cohortFMCDay(List spwb_day, List x, DataFrame SpParams) {
  DataFrame plants = Rcpp::as<Rcpp::DataFrame>(spwb_day["Plants"]);
  
  //Draw cohort-based variables
  DataFrame cohorts = Rcpp::as<Rcpp::DataFrame>(x["cohorts"]);
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  StringVector cohNames = above.attr("row.names");
  int numCohorts = cohNames.size();
  List control = x["control"];
  String transpirationMode = control["transpirationMode"];
  
  
  NumericVector leafFMC(numCohorts,0.0);
  leafFMC.attr("names") = cohNames;
  
  //FMCmax
  IntegerVector SP = cohorts["SP"];
  NumericVector FMCmax = speciesNumericParameter(SP, SpParams, "maxFMC");
  
  if(transpirationMode == "Granier") {
    NumericVector LeafPI0 = speciesNumericParameter(SP, SpParams, "LeafPI0");
    NumericVector LeafEPS = speciesNumericParameter(SP, SpParams, "LeafEPS");
    NumericVector LeafAF = speciesNumericParameter(SP, SpParams, "LeafAF");
    NumericVector PlantPsi = Rcpp::as<Rcpp::NumericVector>(plants["PlantPsi"]);
    NumericVector StemPLC = Rcpp::as<Rcpp::NumericVector>(plants["StemPLC"]);
    
    //Average values for Mediterranean climate species
    for(int c=0;c<numCohorts;c++) {
      if(NumericVector::is_na(LeafPI0[c])) {
        String s = cohNames[c];
        LeafPI0[c] = -2.0;  
        warning("Default LeafPI0 = %f for cohort %s", LeafPI0[c],s.get_cstring());
      }
      if(NumericVector::is_na(LeafEPS[c])) {
        String s = cohNames[c];
        LeafEPS[c] = 17.0; 
        warning("Default LeafEPS = %f for cohort %s",LeafEPS[c],s.get_cstring());
      }
      if(NumericVector::is_na(LeafAF[c])) {
        String s = cohNames[c];
        LeafAF[c] = 0.29; 
        warning("Default LeafAF = %f for cohort %s",LeafAF[c],s.get_cstring());
      }
      double rwc = symplasticRelativeWaterContent(PlantPsi(c), LeafPI0[c], LeafEPS[c]);
      leafFMC[c] = (rwc*(1.0-LeafAF[c])+ (1.0 - StemPLC[c])*LeafAF[c])*FMCmax[c];
    }
    return(leafFMC);
  } else {
    DataFrame paramsAnatomy = Rcpp::as<Rcpp::DataFrame>(x["paramsAnatomy"]);
    DataFrame paramsTranspiration = Rcpp::as<Rcpp::DataFrame>(x["paramsTranspiration"]);
    NumericVector VCleaf_c = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["VCleaf_c"]);
    NumericVector VCleaf_d = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["VCleaf_d"]);
    
    //Water storage parameters
    DataFrame paramsWaterStorage = Rcpp::as<Rcpp::DataFrame>(x["paramsWaterStorage"]);
    NumericVector LeafAF = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["LeafAF"]);
    
    NumericVector StemPLC = Rcpp::as<Rcpp::NumericVector>(plants["StemPLC"]);
    NumericVector psiapoleaf = Rcpp::as<Rcpp::NumericVector>(plants["LeafPsiMin"]);
    NumericVector RWCsymleaf = Rcpp::as<Rcpp::NumericVector>(plants["LeafSympRWC"]);
    for(int c=0;c<numCohorts;c++) {
      double f_apo_leaf = LeafAF[c];
      double leafc = VCleaf_c[c];
      double leafd = VCleaf_d[c];
      double rwc_apo_leaf = apoplasticRelativeWaterContent(psiapoleaf[c], leafc, leafd);
      double rwc_sym_leaf = RWCsymleaf[c];
      double rwc_leaf = rwc_apo_leaf*f_apo_leaf + rwc_sym_leaf*(1.0 - f_apo_leaf);
      leafFMC[c] = rwc_leaf*FMCmax[c];
    }
    return(leafFMC);
  }
}
