#include <Rcpp.h>
#include "forestutils.h"
#include "paramutils.h"
#include "incbeta.h"
#include <math.h>
#include <meteoland.h>
using namespace Rcpp;

const double SIGMA_Wm2 = 5.67*1e-8;

//' Advanced radiation transfer functions
//' 
//' Functions \code{light_layerDirectIrradianceFraction} and \code{light_layerDiffuseIrradianceFraction} calculate 
//' the fraction of above-canopy direct and diffuse radiation reaching each vegetation layer. 
//' Function \code{light_layerSunlitFraction} calculates the proportion of sunlit leaves in each vegetation layer. 
//' Function \code{light_cohortSunlitShadeAbsorbedRadiation} calculates the amount of radiation absorbed 
//' by cohort and vegetation layers, while differentiating between sunlit and shade leaves.
//' 
//' @param LAIme A numeric matrix of live expanded LAI values per vegetation layer (row) and cohort (column).
//' @param LAImd A numeric matrix of dead LAI values per vegetation layer (row) and cohort (column).
//' @param LAImx A numeric matrix of maximum LAI values per vegetation layer (row) and cohort (column).
//' @param k A vector of light extinction coefficients.
//' @param kb,kDIR A vector of direct light extinction coefficients.
//' @param kd A vector of diffuse light extinction coefficients.
//' @param Ib0 Above-canopy direct incident radiation.
//' @param Id0 Above-canopy diffuse incident radiation.
//' @param leafAngle Average leaf inclination angle (in radians)
//' @param p,q Parameters of the beta distribution for leaf angles
//' @param ClumpingIndex The extent to which foliage has a nonrandom spatial distribution.
//' @param alpha A vector of leaf absorbance by species.
//' @param gamma A vector of leaf reflectance values.
//' @param beta Solar elevation (in radians).
//' @param alphaSWR A vecfor of hort-wave absorbance coefficients for each cohort.
//' @param gammaSWR A vector of short-wave reflectance coefficients (albedo) for each cohort.
//' @param ddd A dataframe with direct and diffuse radiation for different subdaily time steps (see function \code{radiation_directDiffuseDay} in package meteoland).
//' @param ntimesteps Number of subdaily time steps.
//' @param trunkExtinctionFraction Fraction of extinction due to trunks (for winter deciduous forests).
//' @param LWRatm Atmospheric downward long-wave radiation (W/m2).
//' @param Tsoil Soil temperature (Celsius).
//' @param Tair Canopy layer air temperature vector (Celsius).
//' 
//' @details
//' Functions for short-wave radiation are adapted from Anten & Bastiaans (2016), 
//' whereas long-wave radiation balance follows Flerchinger et al. (2009). 
//' Vegetation layers are assumed to be ordered from bottom to top.
//' 
//' @return
//' Functions \code{light_layerDirectIrradianceFraction}, \code{light_layerDiffuseIrradianceFraction}
//' and \code{light_layerSunlitFraction} return a numeric vector of length equal to the number of vegetation layers. 
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
//' solarElevation <- 0.67 # in radians
//' SWR_direct <- 1100
//' SWR_diffuse <- 300
//' PAR_direct <- 550
//' PAR_diffuse <- 150
//' 
//' LAI <- 2
//' nlayer <- 10
//' LAIlayerlive <- matrix(rep(LAI/nlayer,nlayer),nlayer,1)
//' LAIlayerdead <- matrix(0,nlayer,1)
//' leafAngle <- 60 # in degrees
//' 
//' ## Extinction coefficients
//' kb <- light_directExtinctionCoefficient(leafAngle*(pi/180), solarElevation)
//' kd_PAR <- 0.5
//' kd_SWR <- kd_PAR/1.35
//' 
//' ## Absorption/Reflection coefficients
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
//' 
//' abs_PAR <- light_cohortSunlitShadeAbsorbedRadiation(PAR_direct, PAR_diffuse,
//'                         Ibfpar, Idfpar, 
//'                         LAIlayerlive, LAIlayerdead, kb, kd_PAR, alpha_PAR, gamma_PAR)
//' abs_SWR <- light_cohortSunlitShadeAbsorbedRadiation(SWR_direct, SWR_diffuse,
//'                          Ibfswr, Idfswr, 
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
//' @name light_advanced
// [[Rcpp::export("light_leafAngleCDF")]]
double leafAngleCDF(double leafAngle, double p, double q) {
  double theta = leafAngle*(2.0/M_PI);
  return(incbeta(p,q,theta));
}

double G_function1(double leafAngle, double solarElevation) {
  double G = NA_REAL;
  if(solarElevation > leafAngle) {
    G = sin(solarElevation)*cos(leafAngle);
  } else {
    double a = sin(solarElevation)*cos(leafAngle)*asin(tan(solarElevation)/tan(leafAngle));
    double b = sqrt(pow(sin(leafAngle), 2.0) - pow(sin(solarElevation), 2.0));
    G = (2.0/M_PI)*(a + b);
  }
  return(G);
}

//' @rdname light_advanced
// [[Rcpp::export("light_directionalExtinctionCoefficient")]]
double directionalExtinctionCoefficient(double p, double q, double solarElevation) {
  double mu = sin(solarElevation);
  double G = 0.0;
  for(int i=0;i<9;i++){
    double a1 = i*10.0*(M_PI/180.0);
    double a2 = (i + 1)*10.0*(M_PI/180.0);
    double Fi = leafAngleCDF(a2, p, q) - leafAngleCDF(a1, p, q);
    G = G + G_function1((a1+a2)/2.0, solarElevation)*Fi;
  }
  return(G/mu);
}


//' @rdname light_advanced
// [[Rcpp::export("light_layerDirectIrradianceFraction")]]
NumericVector layerDirectIrradianceFraction(NumericMatrix LAIme, NumericMatrix LAImd,NumericMatrix LAImx, 
                                            NumericVector kb, NumericVector ClumpingIndex, 
                                            NumericVector alpha, NumericVector gamma, double trunkExtinctionFraction = 0.1) {
  int nlayer = LAIme.nrow();
  int ncoh = LAIme.ncol();
  NumericVector Ifraction(nlayer);
  double s = 0.0;
  double gamma_i = 0.0;
  for(int i=nlayer-1;i>=0;i--) { //Start from top layer
    // I fraction for the current layer (no extinction for top layer)
    Ifraction[i] = (1.0 - gamma_i)*exp(-1.0*s);
    //for subsequent layers increase s
    //Extinction is the maximum between the sum of dead(standing) and expanded leaves and a fraction of maximum live leaves corresponding to trunks (for winter-deciduous stands)
    double gsum = 0.0;
    double lsum = 0.0;
    for(int j =0;j<ncoh;j++) {
      gsum = gamma[j]*(LAIme(i,j)+LAImd(i,j));
      lsum = LAIme(i,j)+LAImd(i,j);
      s = s + (kb[j]*pow(alpha[j],0.5)*ClumpingIndex[j]*(std::max(LAIme(i,j)+LAImd(i,j), trunkExtinctionFraction*LAImx(i,j))));
    }
    //Average gamma according to LAI
    gamma_i = gsum/lsum;
    if(lsum == 0.0) gamma_i = 0.0;
  }
  return(Ifraction);
}

double groundDirectIrradianceFraction(NumericMatrix LAIme, NumericMatrix LAImd,NumericMatrix LAImx, 
                                            NumericVector kb, NumericVector ClumpingIndex, 
                                            NumericVector alpha, double trunkExtinctionFraction = 0.1) {
  int nlayer = LAIme.nrow();
  int ncoh = LAIme.ncol();
  double s = 0.0;
  for(int i=nlayer-1;i>=0;i--) { //Start from top layer
    //for subsequent layers increase s
    //Extinction is the maximum between the sum of dead(standing) and expanded leaves and a fraction of maximum live leaves corresponding to trunks (for winter-deciduous stands)
    for(int j =0;j<ncoh;j++) {
      s = s + (kb[j]*pow(alpha[j],0.5)*ClumpingIndex[j]*(std::max(LAIme(i,j)+LAImd(i,j), trunkExtinctionFraction*LAImx(i,j))));
    }
  }
  return(exp(-1.0*s));
}



//' @rdname light_advanced
// [[Rcpp::export("light_layerDiffuseIrradianceFraction")]]
NumericMatrix layerDiffuseIrradianceFraction(NumericMatrix LAIme, NumericMatrix LAImd,NumericMatrix LAImx, 
                                             NumericMatrix K, NumericVector ClumpingIndex, NumericVector ZF,
                                             NumericVector alpha, NumericVector gamma, double trunkExtinctionFraction = 0.1) {
   int nlayer = LAIme.nrow();
   int ncoh = LAIme.ncol();
   int nZ = ZF.length();
  
  NumericMatrix Ifraction(nZ,nlayer); //Fraction of irradiance from direction k in each layer
   for(int k = 0;k<nZ;k++) { //Sky fractions
     double s = 0.0;
     double gamma_i = 0.0;
     for(int i=nlayer-1;i>=0;i--) { //Start from top layer
       // I fraction for the current layer (no extinction for top layer) and sky fraction
       Ifraction(k,i) = ZF[k]*(1.0 - gamma_i)*exp(-1.0*s);
       //for subsequent layers increase s
       //Extinction is the maximum between the sum of dead(standing) and expanded leaves and a fraction of maximum live leaves corresponding to trunks (for winter-deciduous stands)
       double gsum = 0.0;
       double lsum = 0.0;
       for(int j =0;j<ncoh;j++) {
         gsum = gamma[j]*(LAIme(i,j)+LAImd(i,j));
         lsum = LAIme(i,j)+LAImd(i,j);
         s = s + (K(k,j)*pow(alpha[j],0.5)*ClumpingIndex[j]*(std::max(LAIme(i,j)+LAImd(i,j), trunkExtinctionFraction*LAImx(i,j))));
       }
       //Average gamma according to LAI
       gamma_i = gsum/lsum;
       if(lsum == 0.0) gamma_i = 0.0;
     }
   }
   return(Ifraction);
 }

double groundDiffuseIrradianceFraction(NumericMatrix LAIme, NumericMatrix LAImd,NumericMatrix LAImx, 
                                       NumericMatrix K, NumericVector ClumpingIndex, NumericVector ZF,
                                       NumericVector alpha, double trunkExtinctionFraction = 0.1) {
  int nlayer = LAIme.nrow();
  int ncoh = LAIme.ncol();
  int nZ = ZF.length();
  
  double gif = 0.0; //Overall fraction of ground irradiance
  for(int k = 0;k<nZ;k++) { //Sky fractions
    double s = 0.0;
    for(int i=nlayer-1;i>=0;i--) { //Start from top layer
      //for subsequent layers increase s
      //Extinction is the maximum between the sum of dead(standing) and expanded leaves and a fraction of maximum live leaves corresponding to trunks (for winter-deciduous stands)
      for(int j =0;j<ncoh;j++) {
        s = s + (K(k,j)*pow(alpha[j],0.5)*ClumpingIndex[j]*(std::max(LAIme(i,j)+LAImd(i,j), trunkExtinctionFraction*LAImx(i,j))));
      }
    }
    // I fraction for the current layer (no extinction for top layer)
    gif = gif + ZF[k]*exp(-1.0*s);
  }
  return(gif);
}

/**
 *  Diffuse light absorbed radiation per unit of leaf area, given diffuse light level at the top of the layer
 *  I_{da, ij}
 */
NumericMatrix cohortDiffuseAbsorbedRadiation(double Id0, NumericMatrix Idf, 
                                             NumericMatrix LAIme, NumericMatrix LAImd, 
                                             NumericMatrix K, NumericVector ClumpingIndex,
                                             NumericVector alpha, NumericVector gamma) {
  int ncoh = alpha.length();
  int nlayer = LAIme.nrow();
  int nZ = K.nrow();
  NumericMatrix Ida(nlayer, ncoh);
  for(int i = 0;i<nlayer;i++) {
    //Initialize to zero for all cohorts
    for(int j = 0;j<ncoh;j++) Ida(i,j) = 0.0;
    for(int k = 0;k<nZ;k++) { //Over all sky zones
      double s = 0.0;
      for(int j = 0;j<ncoh;j++) s += K(k,i)*pow(alpha[j],0.5)*ClumpingIndex[j]*(LAIme(i,j)+LAImd(i,j));
      for(int j = 0;j<ncoh;j++) {
        Ida(i,j) = Ida(i,j) + Id0*(1.0-gamma[j])*Idf(k,i)*pow(alpha[j],0.5)*K(k,i)*exp(-1.0*s);
      }
    }
  }
  return(Ida);
}


/**
 *  Scattered light absorbed radiation per unit of leaf area, given direct light level at the top of the layer
 *  I_{bsa, ij}
 */
NumericMatrix cohortScatteredAbsorbedRadiation(double Ib0, NumericVector Ibf, 
                                               NumericMatrix LAIme, NumericMatrix LAImd, 
                                               NumericVector kb, NumericVector ClumpingIndex,
                                               NumericVector alpha, NumericVector gamma) {
  int ncoh = alpha.size();
  int nlayer = Ibf.size();
  NumericMatrix Ibsa(nlayer, ncoh);
  for(int i = 0;i<nlayer;i++) {
    double s1 = 0.0, s2=0.0;
    for(int j = 0;j<ncoh;j++) {
      s1 += kb[j]*alpha[j]*ClumpingIndex[j]*(LAIme(i,j)+LAImd(i,j));
      s2 += kb[j]*ClumpingIndex[j]*(LAIme(i,j)/2.0+LAImd(i,j)/2.0);
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
//' @rdname light_advanced
// [[Rcpp::export("light_cohortSunlitShadeAbsorbedRadiation")]]
List cohortSunlitShadeAbsorbedRadiation(double Ib0, double Id0,
                                        NumericMatrix LAIme, NumericMatrix LAImd, NumericMatrix LAImx,
                                        NumericVector kb, NumericMatrix K, NumericVector ZF, NumericVector ClumpingIndex, 
                                        NumericVector alpha, NumericVector gamma, double trunkExtinctionFraction = 0.1) {
  int ncoh = alpha.size();
  int nlayer = LAIme.nrow();
  
  NumericVector Ibf = layerDirectIrradianceFraction(LAIme,LAImd,LAImx, 
                                                    kb, ClumpingIndex, 
                                                    alpha, gamma, trunkExtinctionFraction);
  NumericMatrix Idf = layerDiffuseIrradianceFraction(LAIme,LAImd,LAImx, 
                                                     K, ClumpingIndex, ZF, 
                                                     alpha, gamma, trunkExtinctionFraction);
  
  NumericMatrix Ida = cohortDiffuseAbsorbedRadiation(Id0, Idf, 
                                                     LAIme, LAImd,
                                                     K, ClumpingIndex, 
                                                     alpha, gamma);
  NumericMatrix Ibsa = cohortScatteredAbsorbedRadiation(Ib0, Ibf, 
                                                        LAIme, LAImd, 
                                                        kb, ClumpingIndex,
                                                        alpha, gamma);
  
  NumericMatrix Ish(nlayer,ncoh); 
  NumericMatrix Isu(nlayer, ncoh);
  // Rcout<<Ib0<<" "<<beta<<" "<<sinb <<" "<<Ib0*alpha[0]*(0.5/sinb)<<"\n";
  for(int j = 0;j<ncoh;j++) {
    for(int i = 0;i<nlayer;i++) {
      Ish(i,j) = Ida(i,j)+Ibsa(i,j); //Absorved radiation in shade leaves (i.e. diffuse+scatter)
      Isu(i,j) = Ish(i,j)+Ib0*kb[j]*alpha[j]; //Absorved radiation in sunlit leaves (i.e. diffuse+scatter+direct)
    }
  }
  List s = List::create(Named("I_sunlit")=Isu, Named("I_shade") = Ish);
  return(s);
}


/**
 *  Sunlit leaf fraction per layer
 *  f_{SL, ij}
 */
//' @rdname light_advanced
// [[Rcpp::export("light_layerSunlitFraction")]]
NumericVector layerSunlitFraction(NumericMatrix LAIme, NumericMatrix LAImd, 
                                  NumericVector kb, NumericVector ClumpingIndex) {
  int ncoh = kb.size();
  int nlayer = LAIme.nrow();
  NumericVector fSL(nlayer);
  double s1=0.0;
  for(int i = nlayer-1;i>=0;i--) {//Start from top layer
    double s2=0.0;
    for(int j = 0;j<ncoh;j++) {
      s1 += kb[j]*ClumpingIndex[j]*(LAIme(i,j)+LAImd(i,j));
      s2 += kb[j]*ClumpingIndex[j]*(LAIme(i,j)+LAImd(i,j));
    }
    fSL[i] = exp(-1.0*s1)*exp(-1.0*s2);
  }
  return(fSL);
}

/*
 * Calculates the amount of radiation absorbed by each cohort
 */
//' @rdname light_advanced
// [[Rcpp::export("light_instantaneousLightExtinctionAbsortion")]]
List instantaneousLightExtinctionAbsortion(NumericMatrix LAIme, NumericMatrix LAImd, NumericMatrix LAImx, 
                                           NumericVector p, NumericVector q, NumericVector ClumpingIndex, 
                                           NumericVector alphaSWR, NumericVector gammaSWR,
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
  NumericVector kDIR(numCohorts, NA_REAL);
  NumericMatrix K9DIR(9, numCohorts);
  NumericVector ZF(numCohorts, 0.1111111); //Uniform overcast sky
  for(int k=0;k<9;k++){
    double zk = ((double)(k*10 + 5))*(M_PI/180.0); // 5, 15, 25,...
    for(int c=0;c<numCohorts;c++) {
      K9DIR(k,c) = directionalExtinctionCoefficient(p[c], q[c], zk);
    }
  }
  for(int c=0;c<numCohorts;c++) {
    alphaPAR[c] = alphaSWR[c]*1.35;
    gammaPAR[c] = gammaSWR[c]*0.8; // (PAR albedo 80% of SWR albedo)
    alphaLWR[c] = 0.97; //Longwave coefficients
  }

  // Rcout<<rad<<" "<< solarElevation[0]<<" "<<SWR_direct[0]<<"\n";
  List fsunlit_list(ntimesteps);
  
  List abs_PAR_SL_COH_list(ntimesteps);
  List abs_SWR_SL_COH_list(ntimesteps);
  List abs_PAR_SH_COH_list(ntimesteps);
  List abs_SWR_SH_COH_list(ntimesteps);
  List abs_PAR_SL_ML_list(ntimesteps);
  List abs_SWR_SL_ML_list(ntimesteps);
  List abs_PAR_SH_ML_list(ntimesteps);
  List abs_SWR_SH_ML_list(ntimesteps);
  
  NumericVector abs_SWR_can(ntimesteps,0.0), abs_LWR_can(ntimesteps,0.0);
  NumericVector abs_SWR_soil(ntimesteps,0.0), abs_LWR_soil(ntimesteps,0.0);
  NumericVector gbf(ntimesteps,0.0), gdf(ntimesteps,0.0);
  for(int n=0;n<ntimesteps;n++) {
    //Calculate direct beam extinction coefficients
    for(int c=0;c<numCohorts;c++) {
      kDIR[c] = directionalExtinctionCoefficient(p[c], q[c], solarElevation[n]);
    }

    //Average sunlit fraction
    NumericVector fsunlit = layerSunlitFraction(LAIme, LAImd, kDIR, ClumpingIndex);
    fsunlit_list[n] = fsunlit;
    
    //Fraction of incoming diffuse/direct SWR radiation and LWR radiation reaching the ground
    gbf[n] = groundDirectIrradianceFraction(LAIme,LAImd,LAImx, 
                                            kDIR, ClumpingIndex, 
                                            alphaSWR, trunkExtinctionFraction);
    gdf[n] = groundDiffuseIrradianceFraction(LAIme,LAImd,LAImx, 
                                             K9DIR, ClumpingIndex, ZF, 
                                             alphaSWR, trunkExtinctionFraction);
    
    //Calculate PAR absorbed radiation for sunlit and shade leaves
    List abs_PAR = cohortSunlitShadeAbsorbedRadiation(PAR_direct[n]*1000.0, PAR_diffuse[n]*1000.0, 
                                                      LAIme, LAImd, LAImx,
                                                      kDIR, K9DIR, ZF, ClumpingIndex,
                                                      alphaPAR, gammaPAR);
    //Calculate SWR absorbed radiation for sunlit and shade leaves
    List abs_SWR = cohortSunlitShadeAbsorbedRadiation(SWR_direct[n]*1000.0, SWR_diffuse[n]*1000.0, 
                                                      LAIme, LAImd, LAImx,
                                                      kDIR, K9DIR, ZF, ClumpingIndex,
                                                      alphaSWR, gammaSWR);
    
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
    double abs_dir_swr = SWR_direct[n]*1000.0*(1.0 - gbf[n]); //W/m2
    double abs_dif_swr = SWR_diffuse[n]*1000.0*(1.0 - gdf[n]); //W/m2
    abs_SWR_can[n] = abs_dir_swr+abs_dif_swr;
    // Rcout<<n<<" "<< abs_SWR_can[n]<< " "<<abs_LWR_can[n]<<"\n";
    //Calculate soil absorved radiation
    abs_SWR_soil[n] = 0.90*((gbf[n]*SWR_direct[n]*1000.0)+(gdf[n]*SWR_diffuse[n]*1000.0)); //10% reflectance for SWR (Geiger, The climate near the ground)

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
  List res = List::create(_["fsunlit"] = fsunlit_list,
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
//' @rdname light_advanced
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