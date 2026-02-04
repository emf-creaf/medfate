#include <RcppArmadillo.h>
#include "medfate.h"
#include "incbeta_c.h"
#include "lightextinction_advanced_c.h"


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

double groundDiffuseIrradianceFraction(const arma::mat& LAIme, const arma::mat& LAImd, const arma::mat& LAImx, 
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

void cohortSunlitShadeAbsorbedRadiation_c(CohortSunlitShadeAbsorbedRadiation_RESULT& res,
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
      res.I_shade(i,j) = Ida(i,j)+Ibsa(i,j);
      res.I_sunlit(i,j) = res.I_shade(i,j)+Ib0*kb[j]*alpha[j]; //Absorbed radiation in sunlit leaves (i.e. diffuse+scatter+direct)
      // Rcout<<i<< " "<< j<<" "<< Ida(i,j)<<" "<< Ibsa(i,j)<<"\n";
    }
  }
}
Rcpp::List copyCohortSunlitShadeAbsorvedRadiationResult_c(CohortSunlitShadeAbsorbedRadiation_RESULT& res) {
  Rcpp::List s = Rcpp::List::create(
    Rcpp::Named("I_sunlit") = Rcpp::wrap(res.I_sunlit), 
    Rcpp::Named("I_shade") = Rcpp::wrap(res.I_shade)
  );
  return(s);
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