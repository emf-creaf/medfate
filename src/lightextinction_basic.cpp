#include <Rcpp.h>
#include "forestutils.h"
#include "paramutils.h"
#include <math.h>
#include <meteoland.h>
using namespace Rcpp;

/***
 * FUNCTIONS FOR LIGHT EXTINCTION (BASIC)
 */
double availableLight(double h, NumericVector H, NumericVector LAI_expanded, NumericVector LAI_dead, NumericVector k, NumericVector CR) {
  double s= 0.0, p=0.0;
  for(int j=0; j< H.size(); j++) {
    double cbh = H[j]*(1.0-CR[j]);
    // p = (H[j]-h)/(H[j]*CR[j]);
    p = leafAreaProportion(h, H[j], cbh, H[j]);
    if(p<0.0) p = 0.0;
    else if(p>1.0) p=1.0;
    s = s + k[j]*p*(LAI_expanded[j]+LAI_dead[j]);
  }
  return(100*exp((-1)*s));
}

NumericVector parcohortC(NumericVector H, NumericVector LAI_expanded,  NumericVector LAI_dead, NumericVector k, NumericVector CR){
  int n = H.size();
  NumericVector ci(n);
  for(int i=0; i<n;i++) ci[i] = availableLight( H[i]*(1.0-(1.0-CR[i])/2.0), H, LAI_expanded, LAI_dead, k, CR);
  ci.attr("names") = H.attr("names");
  return(ci);
}

// [[Rcpp::export(".parcohort")]]
NumericVector parcohort(IntegerVector SP, NumericVector H, NumericVector CR, NumericVector LAI, DataFrame SpParams){
  int n = SP.size();
  NumericVector LAI_dead(n, 0.0);
  NumericVector kPAR = speciesNumericParameterWithImputation(SP, SpParams, "kPAR", true, true);  
  return(parcohortC(H,LAI,LAI_dead,kPAR,CR));
}

//' Radiation extinction functions used in basic transpiration sub-model
//' 
//' @param x An object of class \code{\link{forest}}
//' @param SpParams A data frame with species parameters (see \code{\link{SpParamsMED}}).
//' @param z A numeric vector with height values.
//' @param gdd Growth degree days.
//' 
//' 
//' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
//' 
//' @seealso  \code{\link{spwb}}, \code{\link{light_advanced}}
//' 
//' @name light_basic
//' @keywords internal
// [[Rcpp::export("light_PARcohort")]]
NumericVector PARcohort(List x, DataFrame SpParams, double gdd = NA_REAL) {
  DataFrame above = forest2aboveground(x, SpParams, gdd, false);
  IntegerVector SP = above["SP"];
  NumericVector H = above["H"];
  NumericVector LAI = above["LAI_expanded"];
  NumericVector CR = above["CR"];
  NumericVector pc = parcohort(SP, H, CR, LAI, SpParams);
  pc.attr("names") = cohortIDs(x, SpParams);
  return(pc);
}
NumericVector parheight(NumericVector heights, IntegerVector SP, NumericVector H, NumericVector CR, NumericVector LAI, DataFrame SpParams){
  int n = SP.size();
  NumericVector kPAR = speciesNumericParameterWithImputation(SP, SpParams, "kPAR", true, true);  
  NumericVector LAI_dead(n, 0.0);
  NumericVector AL(heights.size());
  for(int i=0; i<heights.size();i++) AL[i] = availableLight(heights[i], H,LAI,LAI_dead, kPAR,CR);
  return(AL);
}

NumericVector swrheight(NumericVector heights, IntegerVector SP, NumericVector H, NumericVector CR, NumericVector LAI, DataFrame SpParams){
  int n = SP.size();
  NumericVector kPAR = speciesNumericParameterWithImputation(SP, SpParams, "kPAR", true, true);  
  NumericVector kSWR(n), LAI_dead(n);
  for(int i=0; i<n;i++) {
    kSWR[i] = kPAR[i]/1.35;
    LAI_dead[i]=0.0;
  }
  NumericVector AL(heights.size());
  for(int i=0; i<heights.size();i++) AL[i] = availableLight(heights[i], H,LAI, LAI_dead, kSWR,CR);
  return(AL);
}

// [[Rcpp::export(".parheight")]]
NumericVector parheight(NumericVector z, List x, DataFrame SpParams, double gdd = NA_REAL) {
  DataFrame above = forest2aboveground(x, SpParams, gdd, false);
  IntegerVector SP = above["SP"];
  NumericVector H = above["H"];
  NumericVector LAI = above["LAI_expanded"];
  NumericVector CR = above["CR"];
  return(parheight(z, SP, H, CR, LAI, SpParams));
}

//' @rdname light_basic
//' @keywords internal
// [[Rcpp::export("light_PARground")]]
double PARground(List x, DataFrame SpParams, double gdd = NA_REAL) {
  DataFrame above = forest2aboveground(x, SpParams, gdd, false);
  NumericVector LAIlive = above["LAI_live"];
  double woodyLAI = sum(LAIlive); //For herb LAI correction
  NumericVector LAIphe = above["LAI_expanded"];
  NumericVector LAIdead = above["LAI_dead"];
  IntegerVector SP = above["SP"];
  NumericVector kPAR = speciesNumericParameterWithImputation(SP, SpParams, "kPAR", true, true);
  int numCohorts = LAIphe.size();
  double s = 0.0;
  for(int c=0;c<numCohorts;c++) {
    s += (kPAR[c]*(LAIphe[c]+LAIdead[c]));
  }
  //Herb layer effects on light extinction and interception
  s += 0.5*herbLAIAllometric(x["herbCover"], x["herbHeight"], woodyLAI);
  //Percentage of irradiance reaching the ground
  double LgroundPAR = 100.0*exp((-1.0)*s);
  return(LgroundPAR);
}

// [[Rcpp::export(".swrheight")]]
NumericVector swrheight(NumericVector z, List x, DataFrame SpParams, double gdd = NA_REAL) {
  DataFrame above = forest2aboveground(x, SpParams, gdd, false);
  IntegerVector SP = above["SP"];
  NumericVector H = above["H"];
  NumericVector LAI = above["LAI_expanded"];
  NumericVector CR = above["CR"];
  return(swrheight(z, SP, H, CR, LAI, SpParams));
}

//' @rdname light_basic
//' @keywords internal
// [[Rcpp::export("light_SWRground")]]
double SWRground(List x, DataFrame SpParams, double gdd = NA_REAL) {
  DataFrame above = forest2aboveground(x, SpParams, gdd, false);
  NumericVector LAIphe = above["LAI_expanded"];
  NumericVector LAIlive = above["LAI_live"];
  double woodyLAI = sum(LAIlive);
  NumericVector LAIdead = above["LAI_dead"];
  IntegerVector SP = above["SP"];
  NumericVector kPAR = speciesNumericParameterWithImputation(SP, SpParams, "kPAR", true, true);
  int numCohorts = LAIphe.size();
  double s = 0.0;
  for(int c=0;c<numCohorts;c++) {
    s += (kPAR[c]*(LAIphe[c]+LAIdead[c]));
  }
  //Herb layer effects on light extinction and interception
  s += 0.5*herbLAIAllometric(x["herbCover"], x["herbHeight"], woodyLAI);
  //Percentage of irradiance reaching the ground
  double LgroundSWR = 100.0*exp((-1.0)*s/1.35);
  return(LgroundSWR);
}

// [[Rcpp::export(".parExtinctionProfile")]]
NumericVector parExtinctionProfile(NumericVector z, List x, DataFrame SpParams, double gdd = NA_REAL, bool includeHerbs = false) {
  DataFrame above = forest2aboveground(x, SpParams, gdd, false);
  IntegerVector SP = above["SP"];
  NumericVector H = above["H"];
  NumericVector LAI = above["LAI_expanded"];
  NumericVector LAIlive = above["LAI_live"];
  double woodyLAI = sum(LAIlive);
  NumericVector CR = above["CR"];
  if(includeHerbs) {
    SP.push_back(0);
    H.push_back(x["herbHeight"]);
    LAI.push_back(herbLAIAllometric(x["herbCover"], x["herbHeight"], woodyLAI));
    CR.push_back(1.0);
  }
  return(parheight(z, SP, H, CR, LAI, SpParams));
}

// [[Rcpp::export(".swrExtinctionProfile")]]
NumericVector swrExtinctionProfile(NumericVector z, List x, DataFrame SpParams, double gdd = NA_REAL, bool includeHerbs = false) {
  DataFrame above = forest2aboveground(x, SpParams,  gdd, false);
  IntegerVector SP = above["SP"];
  NumericVector H = above["H"];
  NumericVector LAI = above["LAI_expanded"];
  NumericVector LAIlive = above["LAI_live"];
  double woodyLAI = sum(LAIlive);
  NumericVector CR = above["CR"];
  if(includeHerbs) {
    SP.push_back(0);
    H.push_back(x["herbHeight"]);
    LAI.push_back(herbLAIAllometric(x["herbCover"], x["herbHeight"], woodyLAI));
    CR.push_back(1.0);
  }
  return(swrheight(z, SP, H, CR, LAI, SpParams));
}



/**
 * Fraction of the incident SWR radiation in a layer that is absorbed
 */
NumericVector layerAbsorbedSWRFractionIncident(NumericMatrix LAIme, NumericMatrix LAImd, NumericVector kSWR) {
  int nlayer = LAIme.nrow();
  int ncoh = LAIme.ncol();
  double s;
  NumericVector f(nlayer);
  for(int l = 0;l<nlayer;l++) {
    s = 0.0;
    for(int c=0;c<ncoh; c++) s+=kSWR[c]*(LAIme(l,c)+LAImd(l,c));
    f[l] = 1.0 - exp(-1.0*s);
  }
  return(f);
}

/**
 * Fraction of the incident SWR radiation in a layer that is absorbed by live leaves of each cohort
 */
NumericMatrix cohortLayerAbsorbedSWRFractionIncident(NumericVector fi, NumericMatrix LAIme, NumericMatrix LAImd, NumericVector kSWR) {
  int nlayer = LAIme.nrow();
  int ncoh = LAIme.ncol();
  NumericMatrix fij(nlayer, ncoh); 
  double s = 0.0;
  for(int l = 0;l<nlayer;l++) {
    s = 0.0;
    for(int c=0;c<ncoh; c++) s+=kSWR[c]*(LAIme(l,c)+LAImd(l,c));
    if(s>0.0) {
      for(int c=0;c<ncoh; c++) fij(l,c) = fi[l]*kSWR[c]*LAIme(l,c)/s; 
    }
  }
  return(fij);
}

/**
 * Fraction of the SWR radiation that is absorbed by each cohort
 */
NumericVector cohortAbsorbedSWRFraction(NumericMatrix LAIme, NumericMatrix LAImd, NumericVector kSWR) {
  NumericVector fi = layerAbsorbedSWRFractionIncident(LAIme, LAImd, kSWR);
  NumericVector fij = cohortLayerAbsorbedSWRFractionIncident(fi, LAIme, LAImd, kSWR);
  int ncoh = LAIme.ncol();
  int nlayer = LAIme.nrow();
  NumericVector fj(ncoh);
  NumericVector rem(nlayer);
  for(int i = 0;i<nlayer;i++) {
    rem[i] = 1.0;
    for(int h = (nlayer-1);h>i;h--) {
      rem[i] = rem[i]*(1.0-fi[h]);
    }
  }
  double s;
  for(int j=0;j<ncoh;j++) {
    s = 0.0;
    for(int i = 0;i<nlayer;i++) {
      s = s + fij(i,j)*rem[i];
    }
    fj[j] = s;
  }
  return(fj);
}

NumericVector cohortAbsorbedSWRFraction(NumericVector z, NumericVector LAI_expanded, NumericVector LAI_dead, NumericVector H, NumericVector CR, NumericVector kPAR) {
  NumericMatrix LAIme =  LAIdistributionVectors(z, LAI_expanded, H, CR);
  NumericMatrix LAImd =  LAIdistributionVectors(z, LAI_dead, H, CR);
  NumericVector kSWR(kPAR.size());
  for(int i=0;i<kPAR.size();i++) kSWR[i] = kPAR[i]/1.35;
  return(cohortAbsorbedSWRFraction(LAIme, LAImd, kSWR));
}

//' @rdname light_basic
//' @keywords internal
// [[Rcpp::export("light_cohortAbsorbedSWRFraction")]]
NumericVector cohortAbsorbedSWRFraction(NumericVector z, List x, DataFrame SpParams, double gdd = NA_REAL) {
   NumericMatrix LAIme =  LAIdistribution(z, x, SpParams, gdd);
   NumericMatrix LAImd(LAIme.nrow(), LAIme.ncol());
   int nlayer = LAIme.nrow();
   int ncoh = LAIme.ncol();
   for(int i=0;i<nlayer;i++) for(int j=0;j<ncoh;j++) LAImd(i,j)=0.0; 
   NumericVector kPAR = cohortNumericParameterWithImputation(x, SpParams, "kPAR");
   NumericVector kSWR(kPAR.size());
   for(int i=0;i<kPAR.size();i++) kSWR[i] = kPAR[i]/1.35;
   NumericVector caswrf = cohortAbsorbedSWRFraction(LAIme, LAImd, kSWR);
   caswrf.attr("names") = cohortIDs(x, SpParams);
   return(caswrf);
 }
