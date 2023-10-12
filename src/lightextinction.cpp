#include <Rcpp.h>
#include "forestutils.h"
#include "paramutils.h"
#include <meteoland.h>
using namespace Rcpp;

const double SIGMA_Wm2 = 5.67*1e-8;
/***
 * FUNCTIONS FOR LIGHT EXTINCTION
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
  NumericVector kPAR = speciesNumericParameterWithImputation(SP, SpParams, "kPAR", true);  
  return(parcohortC(H,LAI,LAI_dead,kPAR,CR));
}

//' Radiation transfer functions
//' 
//' Functions \code{light_layerIrradianceFraction} and \code{light_layerIrradianceFractionBottomUp} calculate 
//' the fraction of above-canopy irradiance (and the soil irradiance, respectively) reaching each vegetation layer. 
//' Function \code{light_layerSunlitFraction} calculates the proportion of sunlit leaves in each vegetation layer. 
//' Function \code{light_cohortSunlitShadeAbsorbedRadiation} calculates the amount of radiation absorved 
//' by cohort and vegetation layers, while differentiating between sunlit and shade leaves.
//' 
//' @param LAIme A numeric matrix of live expanded LAI values per vegetation layer (row) and cohort (column).
//' @param LAImd A numeric matrix of dead LAI values per vegetation layer (row) and cohort (column).
//' @param LAImx A numeric matrix of maximum LAI values per vegetation layer (row) and cohort (column).
//' @param k A vector of light extinction coefficients.
//' @param kb A vector of direct light extinction coefficients.
//' @param kd A vector of diffuse light extinction coefficients.
//' @param Ib0 Above-canopy direct incident radiation.
//' @param Id0 Above-canopy diffuse incident radiation.
//' @param Ibf Fraction of above-canopy direct radiation reaching each vegetation layer.
//' @param Idf Fraction of above-canopy diffuse radiation reaching each vegetation layer.
//' @param alpha A vecfor of leaf absorbance by species.
//' @param beta Solar elevation (in radians).
//' @param gamma Vector of canopy reflectance values.
//' @param kPAR A vector of visible light extinction coefficients for each cohort.
//' @param alphaSWR A vecfor of hort-wave absorbance coefficients for each cohort.
//' @param gammaSWR A vector of short-wave reflectance coefficients (albedo) for each cohort.
//' @param ddd A dataframe with direct and diffuse radiation for different subdaily time steps (see function \code{radiation_directDiffuseDay} in package meteoland).
//' @param ntimesteps Number of subdaily time steps.
//' @param trunkExtinctionFraction Fraction of extinction due to trunks (for winter deciduous forests).
//' @param LWRatm Atmospheric downward long-wave radiation (W/m2).
//' @param Tsoil Soil temperature (Celsius).
//' @param Tair Canopy layer air temperature vector (Celsius).
//' @param x An object of class \code{\link{forest}}
//' @param SpParams A data frame with species parameters (see \code{\link{SpParamsMED}}).
//' @param z A numeric vector with height values.
//' @param gdd Growth degree days.
//' 
//' @details
//' Functions for short-wave radiation are adapted from Anten & Bastiaans (2016), 
//' whereas long-wave radiation balance follows Flerchinger et al. (2009). 
//' Vegetation layers are assumed to be ordered from bottom to top.
//' 
//' @return
//' Functions \code{light_layerIrradianceFraction}, \code{light_layerIrradianceFractionBottomUp}  and \code{light_layerSunlitFraction} 
//' return a numeric vector of length equal to the number of vegetation layers. 
//' 
//' Function \code{light_cohortSunlitShadeAbsorbedRadiation} returns a list with 
//' two elements (matrices): \code{I_sunlit} and \code{I_shade}.
//' 
//' @references
//' Anten, N.P.R., Bastiaans, L., 2016. The use of canopy models to analyze light competition among plants, in: Hikosaka, K., Niinemets, U., Anten, N.P.R. (Eds.), Canopy Photosynthesis: From Basics to Application. Springer, pp. 379–398.
//' 
//' Flerchinger, G. N., Xiao, W., Sauer, T. J., Yu, Q. 2009. Simulation of within-canopy radiation exchange. NJAS - Wageningen Journal of Life Sciences 57 (1): 5–15. https://doi.org/10.1016/j.njas.2009.07.004.
//' 
//' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
//' 
//' @seealso  \code{\link{spwb}}
//' 
//' @examples
//' LAI <- 2
//' nlayer <- 10
//' LAIlayerlive <- matrix(rep(LAI/nlayer,nlayer),nlayer,1)
//' LAIlayerdead <- matrix(0,nlayer,1)
//' kb <- 0.8
//' kd_PAR <- 0.5
//' kd_SWR <- kd_PAR/1.35
//' alpha_PAR <- 0.9
//' gamma_PAR <- 0.04
//' gamma_SWR <- 0.05
//' alpha_SWR <- 0.7
//' 
//' Ibfpar <- light_layerIrradianceFraction(LAIlayerlive,LAIlayerdead,LAIlayerlive,kb, alpha_PAR)
//' Idfpar <- light_layerIrradianceFraction(LAIlayerlive,LAIlayerdead,LAIlayerlive,kd_PAR, alpha_PAR)
//' Ibfswr <- light_layerIrradianceFraction(LAIlayerlive,LAIlayerdead,LAIlayerlive,kb, alpha_SWR)
//' Idfswr <- light_layerIrradianceFraction(LAIlayerlive,LAIlayerdead,LAIlayerlive,kd_SWR, alpha_SWR)
//' fsunlit <- light_layerSunlitFraction(LAIlayerlive, LAIlayerdead, kb)
//' SHarea <- (1-fsunlit)*LAIlayerlive[,1] 
//' SLarea <- fsunlit*LAIlayerlive[,1] 
//' 
//' oldpar <- par(mar=c(4,4,1,1), mfrow=c(1,2))
//' plot(Ibfpar*100, 1:nlayer,type="l", ylab="Layer", 
//'      xlab="Percentage of irradiance", xlim=c(0,100), ylim=c(1,nlayer), col="dark green")
//' lines(Idfpar*100, 1:nlayer, col="dark green", lty=2)
//' lines(Ibfswr*100, 1:nlayer, col="red")
//' lines(Idfswr*100, 1:nlayer, col="red", lty=2)
//'   
//' plot(fsunlit*100, 1:nlayer,type="l", ylab="Layer", 
//'      xlab="Percentage of leaves", xlim=c(0,100), ylim=c(1,nlayer))
//' lines((1-fsunlit)*100, 1:nlayer, lty=2)
//' par(oldpar)  
//'   
//' solarElevation <- 0.67
//' SWR_direct <- 1100
//' SWR_diffuse <- 300
//' PAR_direct <- 550
//' PAR_diffuse <- 150
//' 
//' abs_PAR <- light_cohortSunlitShadeAbsorbedRadiation(PAR_direct, PAR_diffuse,
//'                         Ibfpar, Idfpar, beta = solarElevation, 
//'                         LAIlayerlive, LAIlayerdead, kb, kd_PAR, alpha_PAR, gamma_PAR)
//' abs_SWR <- light_cohortSunlitShadeAbsorbedRadiation(SWR_direct, SWR_diffuse,
//'                          Ibfswr, Idfswr, beta = solarElevation, 
//'                          LAIlayerlive, LAIlayerdead, kb, kd_SWR, alpha_SWR, gamma_SWR)
//' oldpar <- par(mar=c(4,4,1,1), mfrow=c(1,2))
//' absRadSL <- abs_SWR$I_sunlit[,1]
//' absRadSH <- abs_SWR$I_shade[,1]
//' lambda <- 546.6507
//' QSL <- abs_PAR$I_sunlit[,1]*lambda*0.836*0.01
//' QSH <- abs_PAR$I_shade[,1]*lambda*0.836*0.01
//' plot(QSL, 1:nlayer,type="l", ylab="Layer", 
//'    xlab="Absorbed PAR quantum flux per leaf area", ylim=c(1,nlayer), col="dark green", 
//'    xlim=c(0,max(QSL)))
//' lines(QSH, 1:nlayer, col="dark green", lty=2)
//' plot(absRadSL, 1:nlayer,type="l", ylab="Layer", 
//'    xlab="Absorbed SWR per leaf area (W/m2)", ylim=c(1,nlayer), col="red", 
//'    xlim=c(0,max(absRadSL)))
//' lines(absRadSH, 1:nlayer, col="red", lty=2)
//' par(oldpar)
//'   
//' @name light
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
  NumericVector kPAR = speciesNumericParameterWithImputation(SP, SpParams, "kPAR", true);  
  NumericVector LAI_dead(n, 0.0);
  NumericVector AL(heights.size());
  for(int i=0; i<heights.size();i++) AL[i] = availableLight(heights[i], H,LAI,LAI_dead, kPAR,CR);
  return(AL);
}

NumericVector swrheight(NumericVector heights, IntegerVector SP, NumericVector H, NumericVector CR, NumericVector LAI, DataFrame SpParams){
  int n = SP.size();
  NumericVector kPAR = speciesNumericParameterWithImputation(SP, SpParams, "kPAR", true);  
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

//' @rdname light
// [[Rcpp::export("light_PARground")]]
double PARground(List x, DataFrame SpParams, double gdd = NA_REAL) {
  DataFrame above = forest2aboveground(x, SpParams, gdd, false);
  NumericVector LAIlive = above["LAI_live"];
  double woodyLAI = sum(LAIlive); //For herb LAI correction
  NumericVector LAIphe = above["LAI_expanded"];
  NumericVector LAIdead = above["LAI_dead"];
  IntegerVector SP = above["SP"];
  NumericVector kPAR = speciesNumericParameterWithImputation(SP, SpParams, "kPAR", true);
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

//' @rdname light
// [[Rcpp::export("light_SWRground")]]
double SWRground(List x, DataFrame SpParams, double gdd = NA_REAL) {
  DataFrame above = forest2aboveground(x, SpParams, gdd, false);
  NumericVector LAIphe = above["LAI_expanded"];
  NumericVector LAIlive = above["LAI_live"];
  double woodyLAI = sum(LAIlive);
  NumericVector LAIdead = above["LAI_dead"];
  IntegerVector SP = above["SP"];
  NumericVector kPAR = speciesNumericParameterWithImputation(SP, SpParams, "kPAR", true);
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

//' @rdname light
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

//' @rdname light
// [[Rcpp::export("light_layerIrradianceFraction")]]
NumericVector layerIrradianceFraction(NumericMatrix LAIme, NumericMatrix LAImd,NumericMatrix LAImx, NumericVector k, NumericVector alpha, double trunkExtinctionFraction = 0.1) {
  int nlayer = LAIme.nrow();
  int ncoh = LAIme.ncol();
  NumericVector Ifraction(nlayer);
  double s = 0.0;
  for(int i=nlayer-1;i>=0;i--) { //Start from top layer
    Ifraction[i] = exp(-1.0*s);
    //for subsequent layers increase s
    //Extinction is the maximum between the sum of dead(standing) and expanded leaves and a fraction of maximum live leaves corresponding to trunks (for winter-deciduous stands)
    for(int j =0;j<ncoh;j++) s = s + (k[j]*pow(alpha[j],0.5)*(std::max(LAIme(i,j)+LAImd(i,j), trunkExtinctionFraction*LAImx(i,j))));
  }
  return(Ifraction);
}

//' @rdname light
// [[Rcpp::export("light_layerIrradianceFractionBottomUp")]]
NumericVector layerIrradianceFractionBottomUp(NumericMatrix LAIme, NumericMatrix LAImd,NumericMatrix LAImx, NumericVector k, NumericVector alpha, double trunkExtinctionFraction = 0.1) {
  int nlayer = LAIme.nrow();
  int ncoh = LAIme.ncol();
  NumericVector Ifraction(nlayer);
  double s = 0.0;
  for(int i=0;i<nlayer;i++) { //Start from bottom layer
    Ifraction[i] = exp(-1.0*s);
    //for subsequent layers increase s
    //Extinction is the maximum between the sum of dead(standing) and expanded leaves and a fraction of maximum live leaves corresponding to trunks (for winter-deciduous stands)
    for(int j =0;j<ncoh;j++) s = s + (k[j]*pow(alpha[j],0.5)*(std::max(LAIme(i,j)+LAImd(i,j), trunkExtinctionFraction*LAImx(i,j))));
  }
  return(Ifraction);
}


double groundDiffuseIrradianceFraction(NumericMatrix LAIme, NumericMatrix LAImd, NumericMatrix LAImx, NumericVector k, double trunkExtinctionFraction = 0.1){
  int nlayer = LAIme.nrow();
  int ncoh = LAIme.ncol();
  double s = 0.0;
  for(int i=nlayer-1;i>=0;i--) { //Start from top layer
    //for subsequent layers increase s
    for(int j =0;j<ncoh;j++) s = s + (k[j]*(std::max(LAIme(i,j)+LAImd(i,j), trunkExtinctionFraction*LAImx(i,j))));
  }
  return(exp(-1.0*s));
}
double groundDirectIrradianceFraction(NumericMatrix LAIme, NumericMatrix LAImd, NumericMatrix LAImx, NumericVector k,  NumericVector alpha, double trunkExtinctionFraction = 0.1){
  int nlayer = LAIme.nrow();
  int ncoh = LAIme.ncol();
  double s = 0.0;
  for(int i=nlayer-1;i>=0;i--) { //Start from top layer
    //for subsequent layers increase s
    for(int j =0;j<ncoh;j++) s = s + (k[j]*pow(alpha[j],0.5)*(std::max(LAIme(i,j)+LAImd(i,j), trunkExtinctionFraction*LAImx(i,j))));
  }
  return(exp(-1.0*s));
}

/**
 *  Diffuse light absorbed radiation per unit of leaf area, given diffuse light level at the top of the layer
 *  I_{da, ij}
 */
NumericMatrix cohortDiffuseAbsorbedRadiation(double Id0, NumericVector Idf, NumericMatrix LAIme, NumericMatrix LAImd, NumericVector kd, NumericVector alpha, NumericVector gamma) {
  int ncoh = alpha.size();
  int nlayer = Idf.size();
  NumericMatrix Ida(nlayer, ncoh);
  for(int i = 0;i<nlayer;i++) {
    double s = 0.0;
    for(int j = 0;j<ncoh;j++) s += kd[j]*pow(alpha[j],0.5)*(LAIme(i,j)/2.0+LAImd(i,j)/2.0);
    for(int j = 0;j<ncoh;j++) {
      Ida(i,j) = Id0*(1.0-gamma[j])*Idf[i]*pow(alpha[j],0.5)*kd[j]*exp(-1.0*s);
    }
  }
  return(Ida);
}


/**
 *  Scattered light absorbed radiation per unit of leaf area, given direct light level at the top of the layer
 *  I_{bsa, ij}
 */
NumericMatrix cohortScatteredAbsorbedRadiation(double Ib0, NumericVector Ibf, NumericMatrix LAIme, NumericMatrix LAImd, NumericVector kb, NumericVector alpha, NumericVector gamma) {
  int ncoh = alpha.size();
  int nlayer = Ibf.size();
  NumericMatrix Ibsa(nlayer, ncoh);
  for(int i = 0;i<nlayer;i++) {
    double s1 = 0.0, s2=0.0;
    for(int j = 0;j<ncoh;j++) {
      s1 += kb[j]*alpha[j]*(LAIme(i,j)/2.0+LAImd(i,j)/2.0);
      s2 += kb[j]*(LAIme(i,j)/2.0+LAImd(i,j)/2.0);
    }
    for(int j = 0;j<ncoh;j++) {
      Ibsa(i,j) = Ib0*(1.0-gamma[j])*Ibf[i]*kb[j]*(pow(alpha[j], 0.5)*exp(-1.0*s1) - (alpha[j]/(1.0-gamma[j]))*exp(-1.0*s2));
    }
  }
  return(Ibsa);
}


/**
 * I_{SU,ij}
 * I_{SH,ij}
 */
//' @rdname light
// [[Rcpp::export("light_cohortSunlitShadeAbsorbedRadiation")]]
List cohortSunlitShadeAbsorbedRadiation(double Ib0, double Id0, NumericVector Ibf, NumericVector Idf, double beta,
                             NumericMatrix LAIme, NumericMatrix LAImd, 
                             NumericVector kb,  NumericVector kd, NumericVector alpha, NumericVector gamma) {
  NumericMatrix Ida = cohortDiffuseAbsorbedRadiation(Id0, Idf, LAIme, LAImd, kd, alpha, gamma);
  NumericMatrix Ibsa = cohortScatteredAbsorbedRadiation(Ib0, Ibf, LAIme, LAImd, kb, alpha, gamma);
  int ncoh = alpha.size();
  int nlayer = Ibf.size();
  // double sinb = sin(beta);
  NumericMatrix Ish(nlayer,ncoh); 
  NumericMatrix Isu(nlayer, ncoh);
  // Rcout<<Ib0<<" "<<beta<<" "<<sinb <<" "<<Ib0*alpha[0]*(0.5/sinb)<<"\n";
  for(int j = 0;j<ncoh;j++) {
    for(int i = 0;i<nlayer;i++) {
      Ish(i,j) = Ida(i,j)+Ibsa(i,j); //Absorved radiation in shade leaves (i.e. diffuse+scatter)
      Isu(i,j) = Ish(i,j)+Ib0*alpha[j]; //Absorved radiation in sunlit leaves (i.e. diffuse+scatter+direct)
    }
  }
  List s = List::create(Named("I_sunlit")=Isu, Named("I_shade") = Ish);
  return(s);
}


/**
 *  Sunlit leaf fraction per layer
 *  f_{SL, ij}
 */
//' @rdname light
// [[Rcpp::export("light_layerSunlitFraction")]]
NumericVector layerSunlitFraction(NumericMatrix LAIme, NumericMatrix LAImd, NumericVector kb) {
  int ncoh = kb.size();
  int nlayer = LAIme.nrow();
  NumericVector fSL(nlayer);
  double s1=0.0;
  for(int i = nlayer-1;i>=0;i--) {//Start from top layer
    double s2=0.0;
    for(int j = 0;j<ncoh;j++) {
      s1 += kb[j]*(LAIme(i,j)+LAImd(i,j));
      s2 += kb[j]*(LAIme(i,j)/2.0+LAImd(i,j)/2.0);
    }
    fSL[i] = exp(-1.0*s1)*exp(-1.0*s2);
  }
  return(fSL);
}

/*
 * Calculates the amount of radiation absorved by each cohort
 */
//' @rdname light
// [[Rcpp::export("light_instantaneousLightExtinctionAbsortion")]]
List instantaneousLightExtinctionAbsortion(NumericMatrix LAIme, NumericMatrix LAImd, NumericMatrix LAImx, 
                                           NumericVector kDIR, NumericVector kPAR, NumericVector alphaSWR, NumericVector gammaSWR,
                                           DataFrame ddd, int ntimesteps = 24, double trunkExtinctionFraction = 0.1) {

  int numCohorts = LAIme.ncol();
  int nz = LAIme.nrow();
  
  NumericVector solarElevation = ddd["SolarElevation"]; //in radians
  NumericVector SWR_direct = ddd["SWR_direct"]; //in kW·m-2
  NumericVector SWR_diffuse = ddd["SWR_diffuse"]; //in kW·m-2
  NumericVector PAR_direct = ddd["PAR_direct"]; //in kW·m-2
  NumericVector PAR_diffuse = ddd["PAR_diffuse"]; //in kW·m-2
  
  //Light PAR/SWR coefficients
  NumericVector gammaPAR(numCohorts); //PAR albedo 
  NumericVector gammaLWR(numCohorts, 0.03); //3% albedo of LWR
  NumericVector alphaPAR(numCohorts), alphaLWR(numCohorts);
  NumericVector kSWR(numCohorts), kLWR(numCohorts);
  for(int c=0;c<numCohorts;c++) {
    kSWR[c] = kPAR[c]/1.35;
    alphaPAR[c] = alphaSWR[c]*1.35;
    gammaPAR[c] = gammaSWR[c]*0.8; // (PAR albedo 80% of SWR albedo)
    alphaLWR[c] = 0.97; //Longwave coefficients
    // kLWR[c] = 0.8;
  }
  
  //Average radiation extinction fractions for direct and diffuse PAR/SWR radiation and LWR radiation
  //Include extinction from trunks in winter
  NumericVector Ibfpar = layerIrradianceFraction(LAIme,LAImd,LAImx, kDIR, alphaPAR, trunkExtinctionFraction);
  NumericVector Idfpar = layerIrradianceFraction(LAIme,LAImd,LAImx, kPAR, alphaPAR, trunkExtinctionFraction);
  NumericVector Ibfswr = layerIrradianceFraction(LAIme,LAImd,LAImx, kDIR, alphaSWR,trunkExtinctionFraction);
  NumericVector Idfswr = layerIrradianceFraction(LAIme,LAImd,LAImx, kSWR, alphaSWR, trunkExtinctionFraction);
  // NumericVector Idflwr_td = layerIrradianceFraction(LAIme,LAImd,LAImx, kLWR, alphaLWR, trunkExtinctionFraction);
  // NumericVector Idflwr_bu = layerIrradianceFractionBottomUp(LAIme,LAImd,LAImx, kLWR, alphaLWR, trunkExtinctionFraction);
  
  //Fraction of incoming diffuse/direct SWR radiation and LWR radiation reaching the ground
  double gbf = groundDirectIrradianceFraction(LAIme,LAImd,LAImx, kDIR, alphaSWR, trunkExtinctionFraction);
  double gdf = groundDiffuseIrradianceFraction(LAIme,LAImd,LAImx, kSWR, trunkExtinctionFraction);
  // double glwr = groundDiffuseIrradianceFraction(LAIme,LAImd,LAImx, kLWR, trunkExtinctionFraction);
  // 
  //Average sunlit fraction
  NumericVector fsunlit = layerSunlitFraction(LAIme, LAImd, kDIR);

  // Rcout<<rad<<" "<< solarElevation[0]<<" "<<SWR_direct[0]<<"\n";

  List abs_PAR_SL_COH_list(ntimesteps);
  List abs_SWR_SL_COH_list(ntimesteps);
  List abs_PAR_SH_COH_list(ntimesteps);
  List abs_SWR_SH_COH_list(ntimesteps);
  // List abs_LWR_SL_COH_list(ntimesteps);
  // List abs_LWR_SH_COH_list(ntimesteps);
  List abs_PAR_SL_ML_list(ntimesteps);
  List abs_SWR_SL_ML_list(ntimesteps);
  List abs_PAR_SH_ML_list(ntimesteps);
  List abs_SWR_SH_ML_list(ntimesteps);
  // List abs_LWR_SL_ML_list(ntimesteps);
  // List abs_LWR_SH_ML_list(ntimesteps);
  NumericVector abs_SWR_can(ntimesteps,0.0), abs_LWR_can(ntimesteps,0.0);
  NumericVector abs_SWR_soil(ntimesteps,0.0), abs_LWR_soil(ntimesteps,0.0);
  for(int n=0;n<ntimesteps;n++) {
    
    //Calculate PAR absorved radiation for sunlit and shade leaves
    List abs_PAR = cohortSunlitShadeAbsorbedRadiation(PAR_direct[n]*1000.0, PAR_diffuse[n]*1000.0, 
                                                      Ibfpar, Idfpar, solarElevation[n],
                                                      LAIme, LAImd, 
                                                      kDIR,  kPAR, alphaPAR, gammaPAR);
    //Calculate sWR absorved radiation for sunlit and shade leaves
    List abs_SWR = cohortSunlitShadeAbsorbedRadiation(SWR_direct[n]*1000.0, SWR_diffuse[n]*1000.0, 
                                                      Ibfswr, Idfswr, solarElevation[n],
                                                      LAIme, LAImd, 
                                                      kDIR,  kSWR, alphaSWR, gammaSWR);
    NumericMatrix mswrsl = abs_SWR["I_sunlit"];
    NumericMatrix mswrsh = abs_SWR["I_shade"];
    NumericMatrix mparsl = abs_PAR["I_sunlit"];
    NumericMatrix mparsh = abs_PAR["I_shade"];

    // Rcout << SWR_direct[n]*1000.0 << " "<< SWR_diffuse[n]*1000.0 << " " << LWR_diffuse[n] << "\n";
    //Multiple layer
    abs_PAR_SL_ML_list[n] = mparsl;
    abs_PAR_SH_ML_list[n] = mparsh;
    abs_SWR_SL_ML_list[n] = mswrsl;
    abs_SWR_SH_ML_list[n] = mswrsh;
    //Aggregate light (PAR, SWR, LWR) for sunlit leaves and shade leaves
    NumericVector vparsl(numCohorts,0.0), vparsh(numCohorts,0.0);
    NumericVector vswrsl(numCohorts,0.0), vswrsh(numCohorts,0.0);
    for(int c=0;c<numCohorts;c++){
      for(int i=0;i<nz;i++){
        mparsl(i,c)=mparsl(i,c)*LAIme(i,c)*fsunlit[i];
        mparsh(i,c)=mparsh(i,c)*LAIme(i,c)*(1.0-fsunlit[i]);
        mswrsl(i,c)=mswrsl(i,c)*LAIme(i,c)*fsunlit[i];
        mswrsh(i,c)=mswrsh(i,c)*LAIme(i,c)*(1.0-fsunlit[i]);
        vparsl[c]+=mparsl(i,c);
        vparsh[c]+=mparsh(i,c);
        vswrsl[c]+=mswrsl(i,c);
        vswrsh[c]+=mswrsh(i,c);
      }
      // Rcout<<"Hola "<<vswrsl[c]<<" "<<vswrsh[c]<<" "<<vparsl[c]<<" "<<vparsh[c]<<"\n";
    }
    
    abs_PAR_SL_COH_list[n] = vparsl;
    abs_PAR_SH_COH_list[n] = vparsh;
    abs_SWR_SL_COH_list[n] = vswrsl;
    abs_SWR_SH_COH_list[n] = vswrsh;
    //Calculate canopy absorbed radiation (includes absortion by trunks in winter)
    double abs_dir_swr = SWR_direct[n]*1000.0*(1.0 - gbf); //W/m2
    double abs_dif_swr = SWR_diffuse[n]*1000.0*(1.0 - gdf); //W/m2
    abs_SWR_can[n] = abs_dir_swr+abs_dif_swr;
    // Rcout<<n<<" "<< abs_SWR_can[n]<< " "<<abs_LWR_can[n]<<"\n";
    //Calculate soil absorved radiation
    abs_SWR_soil[n] = 0.90*((gbf*SWR_direct[n]*1000.0)+(gdf*SWR_diffuse[n]*1000.0)); //10% reflectance for SWR (Geiger, The climate near the ground)

    // Rcout<<n<<" PAR : "<<(PAR_direct[n]*1000.0)+(PAR_diffuse[n]*1000.0)<<" SWR: "<<(SWR_direct[n]*1000.0)+(SWR_diffuse[n]*1000.0) <<" can: "<< abs_SWR_can[n]<< " soil: "<< abs_SWR_soil[n]<<" LWR: "<< (LWR_diffuse[n]) <<" can: "<< abs_LWR_can[n]<< " soil: "<< abs_LWR_soil[n]<<"\n";
  }
  List multilayer = List::create(_["PAR_SL"] = abs_PAR_SL_ML_list,
                                 _["PAR_SH"] = abs_PAR_SH_ML_list,
                                 _["SWR_SL"] = abs_SWR_SL_ML_list,
                                 _["SWR_SH"] = abs_SWR_SH_ML_list);
  List sunshade = List::create(_["PAR_SL"] = abs_PAR_SL_COH_list,
                               _["PAR_SH"] = abs_PAR_SH_COH_list,
                               _["SWR_SL"] = abs_SWR_SL_COH_list,
                               _["SWR_SH"] = abs_SWR_SH_COH_list);
  List res = List::create(_["fsunlit"] = fsunlit,
                          _["multilayer"] = multilayer,
                          _["sunshade"] = sunshade,
                          _["SWR_can"] = abs_SWR_can,
                          _["SWR_soil"] = abs_SWR_soil,
                          _["gbf"] = gbf, //ground direct SWR fraction
                          _["gdf"] = gdf); //ground diffuse LWR fraction
  return(res);
}

/**
 *  LWR model of Ma and Liu (2019), based on Flerchinger et al (2009)
 *  
 *  Ma Y, Liu H (2019) An Advanced Multiple-Layer Canopy Model in the WRF Model With Large-Eddy Simulations to Simulate Canopy Flows and Scalar Transport Under Different Stability Conditions. J Adv Model Earth Syst 11:2330–2351. https://doi.org/10.1029/2018MS001347
 *  Flerchinger GN, Xiao W, Sauer TJ, Yu Q (2009) Simulation of within-canopy radiation exchange. NJAS - Wageningen J Life Sci 57:5–15. https://doi.org/10.1016/j.njas.2009.07.004
 */
//' @rdname light
// [[Rcpp::export("light_longwaveRadiationSHAW")]]
List longwaveRadiationSHAW(NumericMatrix LAIme, NumericMatrix LAImd, NumericMatrix LAImx, 
                           double LWRatm, double Tsoil, NumericVector Tair, double trunkExtinctionFraction = 0.1) {
  int ncoh = LAIme.ncol();
  int ncanlayers = Tair.size();
  NumericVector Lup(ncanlayers), Ldown(ncanlayers), Lnet(ncanlayers);
  NumericVector tau(ncanlayers), sumTauComp(ncanlayers);
  NumericMatrix lai_ij(ncanlayers, ncoh);
  NumericMatrix tauM(ncanlayers, ncoh);
  NumericMatrix LnetM(ncanlayers, ncoh);
  if(ncoh>0) LnetM.attr("dimnames") = List::create(seq(1,ncanlayers), seq(1,ncoh));
  
  double Kdlw = 0.7815; //Extinction coefficient fo LWR
  double eps_c = 0.97;
  double eps_g = 0.97;
  //Transmissivity
  for(int i=0;i<ncanlayers;i++) {
    double lai_layer = 0.0;
    sumTauComp[i] = 0.0;
    for(int j=0;j<ncoh;j++) {
      lai_ij(i,j) = std::max(LAIme(i,j)+LAImd(i,j), trunkExtinctionFraction*LAImx(i,j));
      tauM(i,j) = exp(-Kdlw*lai_ij(i,j));
      sumTauComp[i] += (1.0-tauM(i,j)); 
      lai_layer +=lai_ij(i,j);
    }
    tau[i] = exp(-Kdlw*lai_layer);
  }
  //Downwards
  for(int i=(ncanlayers-1);i>=0;i--) {
    double Ldown_upper = 0.0;
    if(i==(ncanlayers-1)) Ldown_upper = LWRatm;
    else Ldown_upper = Ldown[i+1];
    Ldown[i] = tau[i]*Ldown_upper + (1.0 - tau[i])*eps_c*SIGMA_Wm2*pow(Tair[i]+273.16,4.0);
  }
  //Upwards
  double Lup_g = (1.0 - eps_g)*Ldown[0] + eps_g*SIGMA_Wm2*pow(Tsoil+273.16,4.0);
  for(int i=0;i<ncanlayers;i++) {
    double Lup_lower = 0.0;
    if(i==0) Lup_lower = Lup_g;
    else Lup_lower = Lup[i-1];
    Lup[i] = tau[i]*Lup_lower + (1.0 - tau[i])*eps_c*SIGMA_Wm2*pow(Tair[i]+273.16,4.0);
  }
  //Net
  for(int i=0;i<ncanlayers;i++) {
    double Lup_lower = 0.0;
    if(i==0) Lup_lower = Lup_g;
    else Lup_lower = Lup[i-1];
    Lnet[i] = eps_c*(1.0 - tau[i])*(Ldown[i]+Lup_lower - 2.0*SIGMA_Wm2*pow(Tair[i]+273.16,4.0));
    for(int j=0;j<ncoh;j++) {
      LnetM(i,j) =0.0;
      if(LAIme(i,j)>0.0) {
        LnetM(i,j) =  Lnet[i]*((1.0-tauM(i,j))/sumTauComp[i]);
        //Correct for the fact that extinction included all leaves and energy balance is on expanded leaves
        LnetM(i,j) = LnetM(i,j)*(LAIme(i,j)/lai_ij(i,j)); 
        
      } 
    }
  }
  double Lnet_g = eps_g*(Ldown[0] - SIGMA_Wm2*pow(Tsoil+273.16,4.0));
  double Lnet_c = sum(Lnet);
  DataFrame LWR = DataFrame::create(_["Ldown"] = Ldown, 
                                    _["Lup"] = Lup,
                                    _["Lnet"] = Lnet);
  return(List::create(_["LWR_layer"] = LWR,
                      _["Ldown_ground"] = Ldown[0],
                      _["Lup_ground"] = Lup_g,
                      _["Lnet_ground"] = Lnet_g,
                      _["Ldown_canopy"] = LWRatm,
                      _["Lup_canopy"] = Lup[(ncanlayers-1)],
                      _["Lnet_canopy"] = Lnet_c,
                      _["Lnet_cohort_layer"] = LnetM));
}