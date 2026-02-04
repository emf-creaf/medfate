#include <RcppArmadillo.h>
#include "medfate.h"
#include "incbeta_c.h"
#include "modelInput_c.h"
#include "biophysicsutils_c.h"
#include "lightextinction_advanced_c.h"
#include "radiation_c.h"


//' @name light_advanced
//' @keywords internal
// [[Rcpp::export("light_leafAngleCDF")]]
double leafAngleCDF_c(double leafAngle, double p, double q) {
  double theta = leafAngle*(2.0/M_PI);
  return(incbeta_c(p,q,theta));
}

double G_function1_c(double leafAngle, double solarElevation) {
  double G = NA_REAL;
  if(solarElevation > leafAngle) {
    G = sin(solarElevation)*cos(leafAngle);
  } else {
    double a = sin(solarElevation)*cos(leafAngle)*asin(tan(solarElevation)/tan(leafAngle));
    double b = sqrt(pow(sin(leafAngle), 2.0) - pow(sin(solarElevation), 2.0));
    if (std::isnan(b)) b = 0;
    G = (2.0/M_PI)*(a + b);
  }
  return(G);
}


//' @rdname light_advanced
//' @keywords internal
// [[Rcpp::export("light_directionalExtinctionCoefficient")]]
double directionalExtinctionCoefficient_c(double p, double q, double solarElevation) {
  double mu = sin(solarElevation);
  double G = 0.0;
  for(int i=0;i<9;i++){
    double a1 = i*10.0*(M_PI/180.0);
    double a2 = (i + 1)*10.0*(M_PI/180.0);
    double Fi = leafAngleCDF_c(a2, p, q) - leafAngleCDF_c(a1, p, q);
    G = G + G_function1_c((a1+a2)/2.0, solarElevation)*Fi;
  }
  return(G/mu);
}


void layerDirectIrradianceFraction_c(std::vector<double>& Ifraction,
                                     const arma::mat& LAIme, const arma::mat& LAImd, const arma::mat& LAImx, 
                                     const std::vector<double>& kb, const std::vector<double>& ClumpingIndex, 
                                     const std::vector<double>& alpha, const std::vector<double>& gamma, 
                                     double trunkExtinctionFraction) {
  int nlayer = LAIme.n_rows;
  int ncoh = LAIme.n_cols;
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
      s = s + (kb[j]*std::sqrt(alpha[j])*ClumpingIndex[j]*(std::max(LAIme(i,j)+LAImd(i,j), trunkExtinctionFraction*LAImx(i,j))));
    }
    //Average gamma according to LAI
    gamma_i = gsum/lsum;
    if(lsum == 0.0) gamma_i = 0.0;
  }
}

void layerDiffuseIrradianceFraction_c(arma::mat& Ifraction,
                                      const arma::mat& LAIme, const arma::mat& LAImd, const arma::mat& LAImx, 
                                      const arma::mat& K, const std::vector<double>& ClumpingIndex, const std::vector<double>& ZF,
                                      const std::vector<double>& alpha, const std::vector<double>& gamma, 
                                      double trunkExtinctionFraction) {
  int nlayer = LAIme.n_rows;
  int ncoh = LAIme.n_cols;
  int nZ = ZF.size();
  
  // Rcout << nlayer << " " << ncoh << " " << nZ<<"\n";
  for(int k = 0;k<nZ;k++) { //Sky fractions
    double s = 0.0;
    double gamma_i = 0.0;
    for(int i=nlayer-1;i>=0;i--) { //Start from top layer
      // I fraction for the current layer (no extinction for top layer) and sky fraction
      Ifraction(k,i) = ZF[k]*(1.0 - gamma_i)*exp(-1.0*s);
      // Rcout<< k<< " "<< i << " "<< s<< " "<< Ifraction(k,i)<<"\n";
      //for subsequent layers increase s
      //Extinction is the maximum between the sum of dead(standing) and expanded leaves and a fraction of maximum live leaves corresponding to trunks (for winter-deciduous stands)
      double gsum = 0.0;
      double lsum = 0.0;
      for(int j =0;j<ncoh;j++) {
        gsum = gamma[j]*(LAIme(i,j)+LAImd(i,j));
        lsum = LAIme(i,j)+LAImd(i,j);
        s = s + (K(k,j)*sqrt(alpha[j])*ClumpingIndex[j]*(std::max(LAIme(i,j)+LAImd(i,j), trunkExtinctionFraction*LAImx(i,j))));
      }
      //Average gamma according to LAI
      gamma_i = gsum/lsum;
      if(lsum == 0.0) gamma_i = 0.0;
    }
  }
}

double groundDirectIrradianceFraction_c(const arma::mat& LAIme, const arma::mat& LAImd, const arma::mat& LAImx, 
                                        const std::vector<double>& kb, const std::vector<double>& ClumpingIndex, 
                                        const std::vector<double>& alpha,
                                        double trunkExtinctionFraction = 0.1) {
  int nlayer = LAIme.n_rows;
  int ncoh = LAIme.n_cols;
  double s = 0.0;
  for(int i=nlayer-1;i>=0;i--) { //Start from top layer
    //for subsequent layers increase s
    //Extinction is the maximum between the sum of dead(standing) and expanded leaves and a fraction of maximum live leaves corresponding to trunks (for winter-deciduous stands)
    for(int j =0;j<ncoh;j++) {
      s = s + (kb[j]*std::sqrt(alpha[j])*ClumpingIndex[j]*(std::max(LAIme(i,j)+LAImd(i,j), trunkExtinctionFraction*LAImx(i,j))));
    }
  }
  return(exp(-1.0*s));
}

double groundDiffuseIrradianceFraction_c(const arma::mat& LAIme, const arma::mat& LAImd, const arma::mat& LAImx, 
                                         const arma::mat& K, const std::vector<double>& ClumpingIndex, const std::vector<double>& ZF,
                                         const std::vector<double>& alpha, double trunkExtinctionFraction = 0.1) {
  int nlayer = LAIme.n_rows;
  int ncoh = LAIme.n_cols;
  int nZ = ZF.size();
  
  double gif = 0.0; //Overall fraction of ground irradiance
  for(int k = 0;k<nZ;k++) { //Sky fractions
    double s = 0.0;
    for(int i=nlayer-1;i>=0;i--) { //Start from top layer
      //for subsequent layers increase s
      //Extinction is the maximum between the sum of dead(standing) and expanded leaves and a fraction of maximum live leaves corresponding to trunks (for winter-deciduous stands)
      for(int j =0;j<ncoh;j++) {
        s = s + (K(k,j)*std::sqrt(alpha[j])*ClumpingIndex[j]*(std::max(LAIme(i,j)+LAImd(i,j), trunkExtinctionFraction*LAImx(i,j))));
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
void cohortDiffuseAbsorbedRadiation_c(arma::mat& Ida, double Id0, const arma::mat& Idf, 
                                      const arma::mat& LAIme, const arma::mat& LAImd, 
                                      const arma::mat& K, const std::vector<double>& ClumpingIndex,
                                      const std::vector<double>& alpha, const std::vector<double>& gamma) {
  int ncoh = alpha.size();
  int nlayer = LAIme.n_rows;
  int nZ = K.n_rows;
  // Rcout << ncoh << " "<<nlayer<< " "<< nZ<<"\n";
  for(int i = 0;i<nlayer;i++) {
    //Initialize to zero for all cohorts
    for(int j = 0;j<ncoh;j++) Ida(i,j) = 0.0;
    for(int k = 0;k<nZ;k++) { //Over all sky zones
      if(std::isnan(Idf(k,i))) throw medfate::MedfateInternalError("NA Idf");
      double s = 0.0;
      for(int j = 0;j<ncoh;j++) s += K(k,j)*std::sqrt(alpha[j])*ClumpingIndex[j]*(LAIme(i,j)+LAImd(i,j));
      for(int j = 0;j<ncoh;j++) {
        Ida(i,j) = Ida(i,j) + Id0*(1.0-gamma[j])*Idf(k,i)*std::sqrt(alpha[j])*K(k,j)*exp(-1.0*s);
      }
    }
  }
}

/**
 *  Scattered light absorbed radiation per unit of leaf area, given direct light level at the top of the layer
 *  I_{bsa, ij}
 */
void cohortScatteredAbsorbedRadiation_c(arma::mat& Ibsa, double Ib0, const std::vector<double>& Ibf, 
                                        const arma::mat& LAIme, const arma::mat& LAImd, 
                                        const std::vector<double>& kb, const std::vector<double>& ClumpingIndex,
                                        const std::vector<double>& alpha, const std::vector<double>& gamma) {
  int ncoh = alpha.size();
  int nlayer = Ibf.size();
  for(int i = 0;i<nlayer;i++) {
    double s1 = 0.0, s2=0.0;
    for(int j = 0;j<ncoh;j++) {
      s1 += kb[j]*std::sqrt(alpha[j])*ClumpingIndex[j]*(LAIme(i,j)+LAImd(i,j));
      s2 += kb[j]*ClumpingIndex[j]*(LAIme(i,j)+LAImd(i,j));
    }
    for(int j = 0;j<ncoh;j++) {
      double diff = std::sqrt(alpha[j])*exp(-1.0*s1) - alpha[j]*exp(-1.0*s2);
      // double diff = exp(-1.0*s1) - exp(-1.0*s2);
      // Rcout<< i << " "<< j << " "<< s1 << " " << s2 << " " << diff<<"\n";
      Ibsa(i,j) = Ib0*Ibf[i]*std::sqrt(alpha[j])*kb[j]*diff;
    }
  }
}

void cohortSunlitShadeAbsorbedRadiation_c(arma::mat& I_sunlit, arma::mat& I_shade,
                                          double Ib0, double Id0,
                                          const arma::mat& LAIme, const arma::mat& LAImd, const arma::mat& LAImx,
                                          const std::vector<double>& kb, const arma::mat& K, const std::vector<double>& ClumpingIndex, const std::vector<double>& ZF, 
                                          const std::vector<double>& alpha, const std::vector<double>& gamma, double trunkExtinctionFraction) {
  int ncoh = alpha.size();
  int nlayer = LAIme.n_rows;
  int nZ = ZF.size();
  
  std::vector<double> Ibf(nlayer, 0.0); 
  arma::mat Idf(nZ,nlayer);
  arma::mat Ida(nlayer, ncoh);
  arma::mat Ibsa(nlayer, ncoh);
  layerDirectIrradianceFraction_c(Ibf, LAIme,LAImd,LAImx, 
                                  kb, ClumpingIndex, 
                                  alpha, gamma, trunkExtinctionFraction);
  layerDiffuseIrradianceFraction_c(Idf, LAIme,LAImd,LAImx, 
                                   K, ClumpingIndex, ZF, 
                                   alpha, gamma, trunkExtinctionFraction);
  cohortDiffuseAbsorbedRadiation_c(Ida, Id0, Idf, 
                                   LAIme, LAImd,
                                   K, ClumpingIndex, 
                                   alpha, gamma);
  cohortScatteredAbsorbedRadiation_c(Ibsa, Ib0, Ibf, 
                                     LAIme, LAImd, 
                                     kb, ClumpingIndex,
                                     alpha, gamma);
  
  // Rcout << Id0 << " " << Ib0 <<"\n";
  // Rcout<<Ib0<<" "<<beta<<" "<<sinb <<" "<<Ib0*alpha[0]*(0.5/sinb)<<"\n";
  for(int i = 0;i<nlayer;i++) {
    for(int j = 0;j<ncoh;j++) {
      if(std::isnan(Ida(i,j))) throw medfate::MedfateInternalError("NA Ida");
      if(std::isnan(Ibsa(i,j))) throw medfate::MedfateInternalError("NA Ibsa");
      // Ish(i,j) = Ida(i,j)+Ibsa(i,j); //Absorbed radiation in shade leaves (i.e. diffuse+scatter)
      I_shade(i,j) = Ida(i,j)+Ibsa(i,j);
      I_sunlit(i,j) = I_shade(i,j)+Ib0*kb[j]*alpha[j]; //Absorbed radiation in sunlit leaves (i.e. diffuse+scatter+direct)
      // Rcout<<i<< " "<< j<<" "<< Ida(i,j)<<" "<< Ibsa(i,j)<<"\n";
    }
  }
}


void layerSunlitFraction_c(std::vector<double>& fSL,
                           const arma::mat& LAIme, const arma::mat& LAImd, 
                           const std::vector<double>& kb, const std::vector<double>& ClumpingIndex) {
  int ncoh = kb.size();
  int nlayer = LAIme.n_rows;
  double s1=0.0;
  for(int i = nlayer-1;i>=0;i--) {//Start from top layer
    double s2=0.0;
    for(int j = 0;j<ncoh;j++) {
      s1 += kb[j]*ClumpingIndex[j]*(LAIme(i,j)+LAImd(i,j));
      s2 += kb[j]*ClumpingIndex[j]*0.5*(LAIme(i,j)+LAImd(i,j));
    }
    fSL[i] = exp(-1.0*s1)*exp(-1.0*s2);
    if(fSL[i]<0.00001) fSL[i] = 0.0; //Avoids overestimation of absorbed radiation per leaf area when sunlit fraction is close to zero
  }
}

void instantaneousLightExtinctionAbsortion_c(InstantaneousLightExtinctionAbsortion_RESULT& res,
                                             const arma::mat& LAIme, const arma::mat& LAImd, const arma::mat& LAImx, 
                                             const std::vector<double>& p, const std::vector<double>& q, const std::vector<double>& ClumpingIndex, 
                                             const std::vector<double>& alphaSWR, const std::vector<double>& gammaSWR,
                                             const DirectDiffuseDay_RESULT& ddd, int ntimesteps = 24, double trunkExtinctionFraction = 0.1) {
  
  int numCohorts = LAIme.n_cols;
  int ncanlayers = LAIme.n_rows;
  
  //Light PAR/SWR coefficients
  std::vector<double> gammaPAR(numCohorts); //PAR albedo 
  std::vector<double> gammaLWR(numCohorts, 0.03); //3% albedo of LWR
  std::vector<double> alphaPAR(numCohorts), alphaLWR(numCohorts);
  std::vector<double> kDIR(numCohorts, NA_REAL);
  std::vector<double> ZF = {0.178, 0.514, 0.308}; //Standard overcast sky
  std::vector<double> Zangles = {15.0, 45.0, 75.0}; 
  int nZ = ZF.size(); 
  arma::mat K9DIR(nZ, numCohorts); //Goudriaan 1988
  for(int k=0;k<nZ;k++){
    for(int c=0;c<numCohorts;c++) {
      K9DIR(k,c) = directionalExtinctionCoefficient_c(p[c], q[c], Zangles[k]*(M_PI/180.0));
    }
  }
  for(int c=0;c<numCohorts;c++) {
    alphaPAR[c] = alphaSWR[c]*1.35;
    gammaPAR[c] = gammaSWR[c]*0.8; // (PAR albedo 80% of SWR albedo)
    alphaLWR[c] = 0.97; //Longwave coefficients
  }
  
  for(int n=0;n<ntimesteps;n++) {
    arma::mat& abs_PAR_sunlit = res.multilayer.PAR_SL[n];
    arma::mat& abs_PAR_shade = res.multilayer.PAR_SH[n];
    arma::mat& abs_SWR_sunlit = res.multilayer.SWR_SL[n];
    arma::mat& abs_SWR_shade = res.multilayer.SWR_SH[n];
    std::vector<double>& vparsl = res.sunshade.PAR_SL[n];
    std::vector<double>& vparsh = res.sunshade.PAR_SH[n];
    std::vector<double>& vswrsl = res.sunshade.SWR_SL[n];
    std::vector<double>& vswrsh = res.sunshade.SWR_SH[n];
    std::vector<double>& fsunlit = res.fsunlit[n];
    
    // n = 8;
    //Calculate direct beam extinction coefficients
    for(int c=0;c<numCohorts;c++) {
      kDIR[c] = directionalExtinctionCoefficient_c(p[c], q[c], std::max(0.0001, ddd.SolarElevation[n]));
    }
    
    //Average sunlit fraction
    layerSunlitFraction_c(fsunlit, LAIme, LAImd, kDIR, ClumpingIndex);

    //Fraction of incoming diffuse/direct SWR radiation and LWR radiation reaching the ground
    res.gbf[n] = groundDirectIrradianceFraction_c(LAIme,LAImd,LAImx, 
                                                  kDIR, ClumpingIndex, 
                                                  alphaSWR, trunkExtinctionFraction);
    res.gdf[n] = groundDiffuseIrradianceFraction_c(LAIme,LAImd,LAImx, 
                                                   K9DIR, ClumpingIndex, ZF, 
                                                   alphaSWR, trunkExtinctionFraction);
    
    //Calculate PAR absorbed radiation for sunlit and shade leaves
    cohortSunlitShadeAbsorbedRadiation_c(abs_PAR_sunlit, abs_PAR_shade, 
                                         ddd.PAR_direct[n]*1000.0, ddd.PAR_diffuse[n]*1000.0, 
                                         LAIme, LAImd, LAImx,
                                         kDIR, K9DIR, ClumpingIndex, ZF,
                                         alphaPAR, gammaPAR);
    
    //Calculate SWR absorbed radiation for sunlit and shade leaves
    cohortSunlitShadeAbsorbedRadiation_c(abs_SWR_sunlit, abs_SWR_shade, 
                                         ddd.SWR_direct[n]*1000.0, ddd.SWR_diffuse[n]*1000.0, 
                                         LAIme, LAImd, LAImx,
                                         kDIR, K9DIR, ClumpingIndex, ZF, 
                                         alphaSWR, gammaSWR);

    //Aggregate light (PAR, SWR, LWR) for sunlit leaves and shade leaves
    for(int c=0;c<numCohorts;c++){
      vparsl[c]  = 0.0;
      vparsh[c]  = 0.0;
      vswrsl[c]  = 0.0;
      vswrsh[c]  = 0.0;
      for(int i=0;i<ncanlayers;i++){
        abs_PAR_sunlit(i,c)=abs_PAR_sunlit(i,c)*LAIme(i,c)*fsunlit[i];
        abs_PAR_shade(i,c)=abs_PAR_shade(i,c)*LAIme(i,c)*(1.0-fsunlit[i]);
        abs_SWR_sunlit(i,c)=abs_SWR_sunlit(i,c)*LAIme(i,c)*fsunlit[i];
        abs_SWR_shade(i,c)=abs_SWR_shade(i,c)*LAIme(i,c)*(1.0-fsunlit[i]);
        vparsl[c]+=abs_PAR_sunlit(i,c);
        vparsh[c]+=abs_PAR_shade(i,c);
        vswrsl[c]+=abs_SWR_sunlit(i,c);
        vswrsh[c]+=abs_SWR_shade(i,c);
      }
    }
    //Calculate canopy absorbed radiation (includes absortion by trunks in winter)
    double abs_dir_swr = ddd.SWR_direct[n]*1000.0*(1.0 - res.gbf[n]); //W/m2
    double abs_dif_swr = ddd.SWR_diffuse[n]*1000.0*(1.0 - res.gdf[n]); //W/m2
    res.SWR_can[n] = abs_dir_swr+abs_dif_swr;
    //Calculate soil absorved radiation
    res.SWR_soil[n] = 0.90*((res.gbf[n]*ddd.SWR_direct[n]*1000.0)+(res.gdf[n]*ddd.SWR_diffuse[n]*1000.0)); //10% reflectance for SWR (Geiger, The climate near the ground)
  }
}

Rcpp::List copyInstantaneousLightExtinctionAbsortionResult_c(const InstantaneousLightExtinctionAbsortion_RESULT& res) {
  int ntimesteps = res.fsunlit.size();
  Rcpp::List fsunlit_list(ntimesteps);
  Rcpp::List abs_PAR_SL_COH_list(ntimesteps);
  Rcpp::List abs_SWR_SL_COH_list(ntimesteps);
  Rcpp::List abs_PAR_SH_COH_list(ntimesteps);
  Rcpp::List abs_SWR_SH_COH_list(ntimesteps);
  Rcpp::List abs_PAR_SL_ML_list(ntimesteps);
  Rcpp::List abs_SWR_SL_ML_list(ntimesteps);
  Rcpp::List abs_PAR_SH_ML_list(ntimesteps);
  Rcpp::List abs_SWR_SH_ML_list(ntimesteps);
  for(int n=0;n<ntimesteps;n++) {
    fsunlit_list[n] = Rcpp::wrap(res.fsunlit[n]);
    abs_PAR_SL_ML_list[n] = Rcpp::wrap(res.multilayer.PAR_SL[n]);
    abs_PAR_SH_ML_list[n] = Rcpp::wrap(res.multilayer.PAR_SH[n]);
    abs_SWR_SL_ML_list[n] = Rcpp::wrap(res.multilayer.SWR_SL[n]);
    abs_SWR_SH_ML_list[n] = Rcpp::wrap(res.multilayer.SWR_SH[n]);
    abs_PAR_SL_COH_list[n] = Rcpp::wrap(res.sunshade.PAR_SL[n]);
    abs_PAR_SH_COH_list[n] = Rcpp::wrap(res.sunshade.PAR_SH[n]);
    abs_SWR_SL_COH_list[n] = Rcpp::wrap(res.sunshade.SWR_SL[n]);
    abs_SWR_SH_COH_list[n] = Rcpp::wrap(res.sunshade.SWR_SH[n]);
  }
  Rcpp::List multilayer = Rcpp::List::create(Rcpp::Named("PAR_SL") = abs_PAR_SL_ML_list,
                                             Rcpp::Named("PAR_SH") = abs_PAR_SH_ML_list,
                                             Rcpp::Named("SWR_SL") = abs_SWR_SL_ML_list,
                                             Rcpp::Named("SWR_SH") = abs_SWR_SH_ML_list);
  Rcpp::List sunshade = Rcpp::List::create(Rcpp::Named("PAR_SL") = abs_PAR_SL_COH_list,
                                           Rcpp::Named("PAR_SH") = abs_PAR_SH_COH_list,
                                           Rcpp::Named("SWR_SL") = abs_SWR_SL_COH_list,
                                           Rcpp::Named("SWR_SH") = abs_SWR_SH_COH_list);
  Rcpp::List r_l = Rcpp::List::create(Rcpp::Named("fsunlit") = fsunlit_list,
                                      Rcpp::Named("multilayer") = multilayer,
                                      Rcpp::Named("sunshade") = sunshade,
                                      Rcpp::Named("SWR_can") = Rcpp::wrap(res.SWR_can),
                                      Rcpp::Named("SWR_soil") = Rcpp::wrap(res.SWR_soil),
                                      Rcpp::Named("gbf") = Rcpp::wrap(res.gbf), //ground direct SWR fraction
                                      Rcpp::Named("gdf") = Rcpp::wrap(res.gdf)); //ground diffuse LWR fraction
  return(r_l);
}

void longwaveRadiationSHAW_inner_c(LongWaveRadiation_RESULT& res, 
                                   const arma::mat& LAIme, const arma::mat& LAImd, const arma::mat& LAImx, 
                                   double LWRatm, double Tsoil, const std::vector<double>& Tair, double trunkExtinctionFraction) {
  int ncoh = LAIme.n_cols;
  int ncanlayers = Tair.size();

  std::vector<double>& Lup = res.LWR_layer.Lup;
  std::vector<double>& Ldown = res.LWR_layer.Ldown;
  std::vector<double>& Lnet = res.LWR_layer.Lnet;
  std::vector<double>& tau = res.LWR_layer.tau;
  std::vector<double>& sumTauComp = res.LWR_layer.sumTauComp;
  
  arma::mat& LnetM = res.Lnet_cohort_layer;
  
  arma::mat lai_ij(ncanlayers, ncoh);
  arma::mat tauM(ncanlayers, ncoh);
  
  double Kdlw = 0.7815; //Extinction coefficient fo LWR
  double eps_c = 0.97;
  double eps_g = 0.97;
  
  double sigma_pow_Tsoil = SIGMA_Wm2*pow(Tsoil+273.16,4.0);
  
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
  double Lup_g = (1.0 - eps_g)*Ldown[0] + eps_g*sigma_pow_Tsoil;
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
  double Lnet_g = eps_g*(Ldown[0] - sigma_pow_Tsoil);
  double Lnet_c = 0.0;
  for(int i=0;i<ncanlayers;i++) Lnet_c += Lnet[i];
  
  res.Ldown_ground = Ldown[0];
  res.Lup_ground = Lup_g;
  res.Lnet_ground = Lnet_g;
  res.Ldown_canopy = LWRatm;
  res.Lup_canopy = Lup[(ncanlayers-1)];
  res.Lnet_canopy = Lnet_c;
}


Rcpp::List copyLongWaveRadiationResult_c(const LongWaveRadiation_RESULT& res) {
  Rcpp::NumericMatrix LnetM = Rcpp::wrap(res.Lnet_cohort_layer);
  int ncanlayers = res.LWR_layer.Ldown.size();
  int ncoh = LnetM.ncol();
  if(ncoh>0) LnetM.attr("dimnames") =  Rcpp::List::create(Rcpp::seq(1,ncanlayers), Rcpp::seq(1,ncoh));
  
  Rcpp::DataFrame LWR_layer = Rcpp::DataFrame::create(Rcpp::Named("tau") = Rcpp::wrap(res.LWR_layer.tau),
                                                      Rcpp::Named("sumTauComp") = Rcpp::wrap(res.LWR_layer.sumTauComp),
                                                      Rcpp::Named("Ldown") = Rcpp::wrap(res.LWR_layer.Ldown), 
                                                      Rcpp::Named("Lup") = Rcpp::wrap(res.LWR_layer.Lup),
                                                      Rcpp::Named("Lnet") = Rcpp::wrap(res.LWR_layer.Lnet));
  Rcpp::List lwr_struct = Rcpp::List::create(Rcpp::Named("LWR_layer") = LWR_layer,
                                             Rcpp::Named("Ldown_ground") = res.Ldown_ground,
                                             Rcpp::Named("Lup_ground") = res.Lup_ground,
                                             Rcpp::Named("Lnet_ground") = res.Lnet_ground,
                                             Rcpp::Named("Ldown_canopy") = res.Ldown_canopy,
                                             Rcpp::Named("Lup_canopy") = res.Lup_canopy,
                                             Rcpp::Named("Lnet_canopy") = res.Lnet_canopy,
                                             Rcpp::Named("Lnet_cohort_layer") = LnetM);
  return(lwr_struct);
}