#include <RcppArmadillo.h>
#include "forestutils.h"
#include "paramutils.h"
#include "incbeta_c.h"
#include "biophysicsutils_c.h"
#include "lightextinction_advanced_c.h"
#include "radiation_c.h"
#include <math.h>
#include <meteoland.h>
using namespace Rcpp;

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
//' @param K A vector of light extinction coefficients.
//' @param kb A vector of direct light extinction coefficients.
//' @param ZF Fraction of sky angles.
//' @param Ib0 Above-canopy direct incident radiation.
//' @param Id0 Above-canopy diffuse incident radiation.
//' @param leafAngle Average leaf inclination angle (in radians).
//' @param leafAngleSD Standard deviation of leaf inclination angle (in radians).
//' @param p,q Parameters of the beta distribution for leaf angles
//' @param ClumpingIndex The extent to which foliage has a nonrandom spatial distribution.
//' @param alpha A vector of leaf absorbance by species.
//' @param gamma A vector of leaf reflectance values.
//' @param solarElevation Solar elevation (in radians).
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
//' @seealso  \code{\link{spwb}}, \code{\link{light_basic}}
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
//' meanLeafAngle <- 60 # in degrees
//' sdLeafAngle <- 20
//' 
//' beta <- light_leafAngleBetaParameters(meanLeafAngle*(pi/180), sdLeafAngle*(pi/180))
//' 
//' ## Extinction coefficients
//' kb <- light_directionalExtinctionCoefficient(beta["p"], beta["q"], solarElevation)
//' kd_PAR <- 0.5
//' kd_SWR <- kd_PAR/1.35
//' @name light_advanced
//' @keywords internal
// [[Rcpp::export("light_leafAngleBetaParameters")]]
NumericVector leafAngleBetaParameters(double leafAngle, double leafAngleSD) {
  double pow_sum = (leafAngleSD*leafAngleSD) + (leafAngle*leafAngle);
  double p_num = 1.0 - pow_sum/(leafAngle*M_PI/2.0);
  double p_den = pow_sum/(leafAngle*leafAngle) - 1.0;
  double p = p_num/p_den;
  double q = (M_PI/(2.0*leafAngle) - 1.0)*p;
  return(NumericVector::create(_["p"] = p,
                               _["q"] = q));
}




//' @rdname light_advanced
//' @keywords internal
// [[Rcpp::export("light_layerDirectIrradianceFraction")]]
NumericVector layerDirectIrradianceFraction(NumericMatrix LAIme, NumericMatrix LAImd,NumericMatrix LAImx, 
                                            NumericVector kb, NumericVector ClumpingIndex, 
                                            NumericVector alpha, NumericVector gamma, double trunkExtinctionFraction = 0.1) {
  int nlayer = LAIme.nrow();
  std::vector<double> Ifraction_vec(nlayer, 0.0);
  layerDirectIrradianceFraction_c(Ifraction_vec,
                                  as<arma::mat>(LAIme), as<arma::mat>(LAImd), as<arma::mat>(LAImx),
                                  as<std::vector<double>>(kb),as<std::vector<double>>(ClumpingIndex),
                                  as<std::vector<double>>(alpha),as<std::vector<double>>(gamma),
                                  trunkExtinctionFraction);
  NumericVector Ifraction = Rcpp::wrap(Ifraction_vec);
  return(Ifraction);
}


//' @rdname light_advanced
//' @keywords internal
// [[Rcpp::export("light_layerDiffuseIrradianceFraction")]]
NumericMatrix layerDiffuseIrradianceFraction(NumericMatrix LAIme, NumericMatrix LAImd,NumericMatrix LAImx, 
                                             NumericMatrix K, NumericVector ClumpingIndex, NumericVector ZF,
                                             NumericVector alpha, NumericVector gamma, double trunkExtinctionFraction = 0.1) {
   int nlayer = LAIme.nrow();
   int nZ = ZF.size();
  
  arma::mat Ifraction_mat(nZ,nlayer); //Fraction of irradiance from direction k in each layer
  layerDiffuseIrradianceFraction_c(Ifraction_mat,
                                  as<arma::mat>(LAIme), as<arma::mat>(LAImd), as<arma::mat>(LAImx),
                                  as<arma::mat>(K), as<std::vector<double>>(ClumpingIndex), as<std::vector<double>>(ZF),
                                  as<std::vector<double>>(alpha), as<std::vector<double>>(gamma),
                                  trunkExtinctionFraction);
  
  // Rcout << nlayer << " " << ncoh << " " << nZ<<"\n";
   NumericMatrix Ifraction = Rcpp::wrap(Ifraction_mat);
   return(Ifraction);
 }




/**
 * I_{SU,ij}
 * I_{SH,ij}
 */
//' @rdname light_advanced
//' @keywords internal
// [[Rcpp::export("light_cohortSunlitShadeAbsorbedRadiation")]]
List cohortSunlitShadeAbsorbedRadiation(double Ib0, double Id0,
                                        NumericMatrix LAIme, NumericMatrix LAImd, NumericMatrix LAImx,
                                        NumericVector kb, NumericMatrix K, NumericVector ClumpingIndex, NumericVector ZF, 
                                        NumericVector alpha, NumericVector gamma, double trunkExtinctionFraction = 0.1) {
  int ncoh = alpha.size();
  int nlayer = LAIme.nrow();
  arma::mat I_sunlit(nlayer, ncoh);
  arma::mat I_shade(nlayer, ncoh);
  cohortSunlitShadeAbsorbedRadiation_c(I_sunlit, I_shade, 
                                       Ib0, Id0,
                                       as<arma::mat>(LAIme), as<arma::mat>(LAImd), as<arma::mat>(LAImx),
                                       as<std::vector<double>>(kb), as<arma::mat>(K), as<std::vector<double>>(ClumpingIndex), as<std::vector<double>>(ZF), 
                                       as<std::vector<double>>(alpha), as<std::vector<double>>(gamma), 
                                       trunkExtinctionFraction);
  return(Rcpp::List::create(
      Rcpp::Named("I_sunlit") = Rcpp::wrap(I_sunlit), 
      Rcpp::Named("I_shade") = Rcpp::wrap(I_shade)
  ));
}


/**
 *  Sunlit leaf fraction per layer
 *  f_{SL, ij}
 */
//' @rdname light_advanced
//' @keywords internal
// [[Rcpp::export("light_layerSunlitFraction")]]
NumericVector layerSunlitFraction(NumericMatrix LAIme, NumericMatrix LAImd, 
                                  NumericVector kb, NumericVector ClumpingIndex) {
  int nlayer = LAIme.nrow();
  std::vector<double> fSL_vec(nlayer, 0.0);
  layerSunlitFraction_c(fSL_vec,
                        as<arma::mat>(LAIme), as<arma::mat>(LAImd),
                        as<std::vector<double>>(kb), as<std::vector<double>>(ClumpingIndex));
  NumericVector fSL = Rcpp::wrap(fSL_vec);
  return(fSL);
}

/*
 * Calculates the amount of radiation absorbed by each cohort
 */
//' @rdname light_advanced
//' @keywords internal
// [[Rcpp::export("light_instantaneousLightExtinctionAbsortion")]]
List instantaneousLightExtinctionAbsortion(NumericMatrix LAIme, NumericMatrix LAImd, NumericMatrix LAImx, 
                                           NumericVector p, NumericVector q, NumericVector ClumpingIndex, 
                                           NumericVector alphaSWR, NumericVector gammaSWR,
                                           DataFrame ddd, int ntimesteps = 24, double trunkExtinctionFraction = 0.1) {

  int numCohorts = LAIme.ncol();
  int ncanlayers = LAIme.nrow();
  
  InstantaneousLightExtinctionAbsortion_RESULT lightExtinctionAbsortion(numCohorts, ncanlayers, ntimesteps);
  DirectDiffuseDay_RESULT directDiffuseDay(ntimesteps);
  
  directDiffuseDay.SolarElevation = as<std::vector<double>> (ddd["SolarElevation"]); //in radians
  directDiffuseDay.SWR_direct = as<std::vector<double>> (ddd["SWR_direct"]); //in kW·m-2
  directDiffuseDay.SWR_diffuse = as<std::vector<double>> (ddd["SWR_diffuse"]); //in kW·m-2
  directDiffuseDay.PAR_direct = as<std::vector<double>> (ddd["PAR_direct"]); //in kW·m-2
  directDiffuseDay.PAR_diffuse = as<std::vector<double>> (ddd["PAR_diffuse"]); //in kW·m-2
  
  instantaneousLightExtinctionAbsortion_c(lightExtinctionAbsortion,
                                          as<arma::mat>(LAIme), as<arma::mat>(LAImd), as<arma::mat>(LAImx),
                                          as<std::vector<double>>(p), as<std::vector<double>>(q), as<std::vector<double>>(ClumpingIndex),
                                          as<std::vector<double>>(alphaSWR), as<std::vector<double>>(gammaSWR),
                                          directDiffuseDay, ntimesteps, trunkExtinctionFraction);
  
  return(copyInstantaneousLightExtinctionAbsortionResult_c(lightExtinctionAbsortion));
}

void longwaveRadiationSHAW_inner(List internalLWR, NumericMatrix LAIme, NumericMatrix LAImd, NumericMatrix LAImx, 
                                 double LWRatm, double Tsoil, NumericVector Tair, double trunkExtinctionFraction = 0.1) {

  int ncanlayers = Tair.size();
  int numCohorts = LAIme.ncol();
  LongWaveRadiation_RESULT res(ncanlayers, numCohorts);
  
  longwaveRadiationSHAW_inner_c(res, 
                                as<arma::mat>(LAIme), as<arma::mat>(LAImd), as<arma::mat>(LAImx),
                                LWRatm, Tsoil, as<std::vector<double>>(Tair), trunkExtinctionFraction);
  Rcpp::NumericMatrix LnetM = Rcpp::wrap(res.Lnet_cohort_layer);
  if(numCohorts>0) LnetM.attr("dimnames") =  Rcpp::List::create(Rcpp::seq(1,ncanlayers), Rcpp::seq(1,numCohorts));
  Rcpp::DataFrame LWR_layer = Rcpp::DataFrame::create(Rcpp::Named("tau") = Rcpp::wrap(res.LWR_layer.tau),
                                                      Rcpp::Named("sumTauComp") = Rcpp::wrap(res.LWR_layer.sumTauComp),
                                                      Rcpp::Named("Ldown") = Rcpp::wrap(res.LWR_layer.Ldown), 
                                                      Rcpp::Named("Lup") = Rcpp::wrap(res.LWR_layer.Lup),
                                                      Rcpp::Named("Lnet") = Rcpp::wrap(res.LWR_layer.Lnet));
  internalLWR["Lnet_cohort_layer"] = LnetM;
  internalLWR["LWR_layer"] = LWR_layer;
  internalLWR["Ldown_ground"] = res.Ldown_ground;
  internalLWR["Lup_ground"] = res.Lup_ground;
  internalLWR["Lnet_ground"] = res.Lnet_ground;
  internalLWR["Ldown_canopy"] = res.Ldown_canopy;
  internalLWR["Lup_canopy"] = res.Lup_canopy;
  internalLWR["Lnet_canopy"] = res.Lnet_canopy;
}


/**
 *  LWR model of Ma and Liu (2019), based on Flerchinger et al (2009)
 *  
 *  Ma Y, Liu H (2019) An Advanced Multiple-Layer Canopy Model in the WRF Model With Large-Eddy Simulations to Simulate Canopy Flows and Scalar Transport Under Different Stability Conditions. J Adv Model Earth Syst 11:2330–2351. https://doi.org/10.1029/2018MS001347
 *  Flerchinger GN, Xiao W, Sauer TJ, Yu Q (2009) Simulation of within-canopy radiation exchange. NJAS - Wageningen J Life Sci 57:5–15. https://doi.org/10.1016/j.njas.2009.07.004
 */
//' @rdname light_advanced
//' @keywords internal
// [[Rcpp::export("light_longwaveRadiationSHAW")]]
List longwaveRadiationSHAW(NumericMatrix LAIme, NumericMatrix LAImd, NumericMatrix LAImx, 
                            double LWRatm, double Tsoil, NumericVector Tair, double trunkExtinctionFraction = 0.1) {
  int ncanlayers = Tair.size();
  int numCohorts = LAIme.ncol();
  LongWaveRadiation_RESULT LWRres(ncanlayers, numCohorts);

  longwaveRadiationSHAW_inner_c(LWRres, 
                                as<arma::mat>(LAIme), as<arma::mat>(LAImd), as<arma::mat>(LAImx),
                                LWRatm, Tsoil, as<std::vector<double>>(Tair), trunkExtinctionFraction);
  
  return(copyLongWaveRadiationResult_c(LWRres));
}

