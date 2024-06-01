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
//' Tissue moisture functions
//' 
//' Set of functions used to calculate tissue moisture from water potential and viceversa.
//' 
//' @param psiSym,psiApo Symplastic or apoplastic water potential (MPa).
//' @param RWC Relative water content [0-1].
//' @param pi0 Full turgor osmotic potential (MPa).
//' @param epsilon Bulk modulus of elasticity (MPa).
//' @param c,d Parameters of the xylem vulnerability curve.
//' @param af Apoplastic fraction (proportion) in the segment (e.g. leaf or stem).
//' @param L Vector with the length of coarse roots (mm) for each soil layer.
//' @param V Vector with the proportion [0-1] of fine roots within each soil layer.
//' @param Al2As Leaf area to sapwood area (in m2·m-2).
//' @param height Plant height (in cm).
//' @param SLA Specific leaf area (mm2·mg-1).
//' @param wd Wood density (g·cm-3).
//' @param ld Leaf tissue density (g·cm-3).
//' 
//' @return
//' Values returned for each function are:
//' \itemize{
//'   \item{\code{moisture_symplasticRWC}: Relative water content [0-1] of the symplastic fraction.}
//'   \item{\code{moisture_apoplasticRWC}: Relative water content [0-1] of the apoplastic fraction.}
//'   \item{\code{moisture_symplasticWaterPotential}: Water potential (in MPa) of the symplastic fraction.}
//'   \item{\code{moisture_apoplasticWaterPotential}: Water potential (in MPa) of the apoplastic fraction.}
//'   \item{\code{moisture_turgorLossPoint}: Water potential (in MPa) corresponding to turgor loss point.}
//'   \item{\code{moisture_segmentRWC}: Segment relative water content [0-1].}
//'   \item{\code{water_plant}: A vector of water content (mm) per plant cohort.}
//' }
//' 
//' @references
//' Bartlett, M.K., Scoffoni, C., Sack, L. 2012. The determinants of leaf turgor loss point and prediction of drought tolerance of species and biomes: a global meta-analysis. Ecology Letters 15: 393–405.
//' 
//' \enc{Hölttä}{Holtta}, T., Cochard, H., Nikinmaa, E., Mencuccini, M. 2009. Capacitive effect of cavitation in xylem conduits: Results from a dynamic model. Plant, Cell and Environment 32: 10–21.
//' 
//' Martin-StPaul, N., Delzon, S., Cochard, H. 2017. Plant resistance to drought depends on timely stomatal closure. Ecology Letters 20: 1437–1447.
//' 
//' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
//' 
//' @seealso
//' \code{\link{hydraulics_psi2K}}, \code{\link{hydraulics_supplyFunctionPlot}}, \code{\link{spwb}}, \code{\link{soil}}
//' 
//' @examples
//' psi = seq(-10,0, by=0.1)
//' rwc_s = rep(NA, length(psi))
//' for(i in 1:length(psi)) rwc_s[i] = moisture_symplasticRWC(psi[i],-3,12)
//' plot(psi, rwc_s, type="l", xlab="Water potential (MPa)", ylab = "Symplasmic RWC")
//' 
//' @name moisture
//' @keywords internal
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
//' @rdname moisture
//' @keywords internal
// [[Rcpp::export("moisture_leafWaterCapacity")]]
double leafWaterCapacity(double SLA, double ld) {
  return(1000.0/(1000.0*ld*SLA))*(1.0- (ld/1.54));
}


//' @rdname moisture
//' @keywords internal
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
//' @rdname moisture
//' @keywords internal
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

//' @rdname moisture
//' @keywords internal
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
//' @rdname moisture
//' @keywords internal
// [[Rcpp::export("moisture_apoplasticRWC")]]
double apoplasticRelativeWaterContent(double psiApo, double c, double d) {
  if(psiApo>=0.0) return(1.0);
  return(exp(-pow(psiApo/d,c)));
}

//' @rdname moisture
//' @keywords internal
// [[Rcpp::export("moisture_apoplasticPsi")]]
double apoplasticWaterPotential(double RWC, double c, double d) {
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
double tissueRelativeWaterContent(double psiSym, double pi0, double epsilon, 
                                  double psiApo, double c, double d, 
                                  double af) {
  double sym_rwc = symplasticRelativeWaterContent(psiSym, pi0, epsilon);
  double apo_rwc = apoplasticRelativeWaterContent(psiApo, c, d); //Water content cannot be higher than the fraction of non-embolized conduits
  return(sym_rwc*(1.0-af)+apo_rwc*af);
}

//' @rdname moisture
//' @keywords internal
// [[Rcpp::export("plant_water")]]
NumericVector plantWaterContent(List x) {
  List control = x["control"];
  String transpirationMode = control["transpirationMode"];
  DataFrame cohorts = Rcpp::as<Rcpp::DataFrame>(x["cohorts"]);
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  NumericVector LAIlive = Rcpp::as<Rcpp::NumericVector>(above["LAI_live"]);
  NumericVector vol(LAIlive.size(), NA_REAL);
  
  DataFrame paramsTransp = Rcpp::as<Rcpp::DataFrame>(x["paramsTranspiration"]);
  DataFrame paramsWaterStorage = Rcpp::as<Rcpp::DataFrame>(x["paramsWaterStorage"]);
  NumericVector StemPI0 = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["StemPI0"]);
  NumericVector StemEPS = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["StemEPS"]);
  NumericVector StemAF = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["StemAF"]);
  NumericVector Vsapwood = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["Vsapwood"]); //l·m-2 = mm
  NumericVector LeafPI0 = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["LeafPI0"]);
  NumericVector LeafEPS = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["LeafEPS"]);
  NumericVector LeafAF = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["LeafAF"]);
  NumericVector Vleaf = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["Vleaf"]); //l·m-2 = mm
  
  DataFrame internalWater = Rcpp::as<Rcpp::DataFrame>(x["internalWater"]);
  if(transpirationMode == "Granier") {
    NumericVector VCstem_c = Rcpp::as<Rcpp::NumericVector>(paramsTransp["VCstem_c"]);
    NumericVector VCstem_d = Rcpp::as<Rcpp::NumericVector>(paramsTransp["VCstem_d"]);
    NumericVector PlantPsi = clone(Rcpp::as<Rcpp::NumericVector>(internalWater["PlantPsi"]));
    for(int c = 0;c< LAIlive.size();c++) {
      double leafrwc = tissueRelativeWaterContent(PlantPsi[c], LeafPI0[c], LeafEPS[c], 
                                                  PlantPsi[c], VCstem_c[c], VCstem_d[c], 
                                                  LeafAF[c]);
      
      double stemrwc = tissueRelativeWaterContent(PlantPsi[c], StemPI0[c], StemEPS[c], 
                                                  PlantPsi[c], VCstem_c[c], VCstem_d[c], 
                                                  StemAF[c]);
      vol[c] = ((Vleaf[c] * leafrwc) + (Vsapwood[c] * stemrwc))*LAIlive[c];
    }
  } else {
    NumericVector LeafPLCVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["LeafPLC"]);
    NumericVector StemPLCVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["StemPLC"]);
    NumericVector LeafSympPsiVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["LeafSympPsi"]);
    NumericVector StemSympPsiVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["StemSympPsi"]);
    for(int c = 0;c< LAIlive.size();c++) {
      double leaf_sym_rwc = symplasticRelativeWaterContent(LeafSympPsiVEC[c], LeafPI0[c], LeafEPS[c]);
      double leaf_apo_rwc = (1.0 - LeafPLCVEC[c]); 
      double leafrwc =(leaf_sym_rwc*(1.0-LeafAF[c])+leaf_apo_rwc*LeafAF[c]);
      
      double stem_sym_rwc = symplasticRelativeWaterContent(StemSympPsiVEC[c], StemPI0[c], StemEPS[c]);
      double stem_apo_rwc = (1.0 - StemPLCVEC[c]); 
      double stemrwc =(stem_sym_rwc*(1.0-StemAF[c])+stem_apo_rwc*StemAF[c]);
      
      vol[c] = ((Vleaf[c] * leafrwc) + (Vsapwood[c] * stemrwc))*LAIlive[c];
    }
  }
    
  vol.attr("names") = above.attr("row.names");
  return(vol);
}