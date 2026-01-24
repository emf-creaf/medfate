#include <cmath>
#include <vector>
#include "soil_c.h"
#include "Rcpp.h"


/*=============================================================================
 * Implementation of Soil class
 *=============================================================================**/
Soil::Soil(int nlayersIn, 
           std::string& modelIn,
           std::vector<double>& widthsIn,
           std::vector<double>& clayIn,
           std::vector<double>& sandIn,
           std::vector<double>& omIn,
           std::vector<double>& nitrogenIn,
           std::vector<double>& phIn,
           std::vector<double>& bdIn,
           std::vector<double>& rfcIn,
           std::vector<double>& macroIn,
           std::vector<double>& KsatIn,
           std::vector<double>& VG_alphaIn,
           std::vector<double>& VG_nIn,
           std::vector<double>& VG_theta_resIn,
           std::vector<double>& VG_theta_satIn,
           std::vector<std::string>& usda_typeIn,
           std::vector<double>& theta_SATIn,
           std::vector<double>& theta_FCIn,
           std::vector<double>& WIn,
           std::vector<double>& psiIn,
           std::vector<double>& thetaIn,
           std::vector<double>& TempIn,
           ClappHornberger& clapp_hornbergerIn)  { 
  nlayers = nlayersIn;
  model = modelIn;
  widths = widthsIn;
  clay = clayIn;
  sand = sandIn;
  om = omIn;
  nitrogen = nitrogenIn;
  ph = phIn;
  bd = bdIn;
  rfc = rfcIn;
  macro = macroIn;
  Ksat = KsatIn;
  VG_alpha = VG_alphaIn;
  VG_n = VG_nIn;
  VG_theta_res = VG_theta_resIn;
  VG_theta_sat = VG_theta_satIn;
  usda_type = usda_typeIn;
  theta_SAT = theta_SATIn;
  theta_FC = theta_FCIn;
  W = WIn;
  psi = psiIn;
  theta = thetaIn;
  Temp = TempIn;
  clapp_hornberger = clapp_hornbergerIn;
}

Soil::Soil(Rcpp::DataFrame x, Rcpp::String model = "VG") {
  nlayers = x.nrow();
  model = model.get_cstring();
  
  widths = Rcpp::as< std::vector<double> >(x["widths"]);
  clay = Rcpp::as< std::vector<double> >(x["clay"]);
  sand = Rcpp::as< std::vector<double> >(x["sand"]);
  om = Rcpp::as< std::vector<double> >(x["om"]);
  ph = std::vector<double>(nlayers);
  if(x.containsElementNamed("ph")) ph = Rcpp::as< std::vector<double> >(x["ph"]);
  nitrogen = std::vector<double>(nlayers);
  if(x.containsElementNamed("nitrogen")) nitrogen = Rcpp::as< std::vector<double> >(x["nitrogen"]);
  bd = Rcpp::as< std::vector<double> >(x["bd"]);
  rfc = Rcpp::as< std::vector<double> >(x["rfc"]);
  macro = Rcpp::as< std::vector<double> >(x["macro"]);
  Ksat = Rcpp::as< std::vector<double> >(x["Ksat"]);
  VG_alpha = Rcpp::as< std::vector<double> >(x["VG_alpha"]);
  VG_n = Rcpp::as< std::vector<double> >(x["VG_n"]);
  VG_theta_res = Rcpp::as< std::vector<double> >(x["VG_theta_res"]);
  VG_theta_sat = Rcpp::as< std::vector<double> >(x["VG_theta_sat"]);
  usda_type = std::vector<std::string>(nlayers);
  theta_FC = std::vector<double>(nlayers);
  theta_SAT = std::vector<double>(nlayers);
  W  = Rcpp::as< std::vector<double> >(x["W"]);
  psi = std::vector<double>(nlayers);
  theta = std::vector<double>(nlayers);
  Temp= Rcpp::as< std::vector<double> >(x["Temp"]);
  for(int l=0;l<nlayers;l++) {
    setW(l, W[l]); // this fills psi and theta based on W
    usda_type[l] = USDAType_c(clay[l], sand[l]);
    if(model=="SX") {
      theta_SAT[l] = thetaSATSaxton_c(clay[l], sand[l], om[l]); 
      theta_FC[l] = psi2thetaSaxton_c(clay[l], sand[l], fieldCapacityPsi, om[l]); 
    } else if(model=="VG") {
      theta_SAT[l] = VG_theta_sat[l]; 
      theta_FC[l] = psi2thetaVanGenuchten_c(VG_n[l], VG_alpha[l], VG_theta_res[l], VG_theta_sat[l], fieldCapacityPsi); 
    }
  }
  clapp_hornberger = ClappHornberger(usda_type[0]);
}
double Soil::getW(int layer) {return W[layer];}
int Soil::getNlayers() {return nlayers; }
std::string Soil::getModel() {return model; }
ClappHornberger Soil::getClappHornberger() {return clapp_hornberger; }
double Soil::getWidth(int layer) {return widths[layer]; }  
double Soil::getClay(int layer) {return clay[layer]; }
double Soil::getSand(int layer) {return sand[layer]; }
double Soil::getOM(int layer) {return om[layer]; }
double Soil::getNitrogen(int layer) {return nitrogen[layer]; }
double Soil::getPH(int layer) {return ph[layer]; }
double Soil::getBD(int layer) {return bd[layer]; }
double Soil::getRFC(int layer) {return rfc[layer]; }
double Soil::getMacro(int layer) {return macro[layer]; }
double Soil::getKsat(int layer) {return Ksat[layer]; }
double Soil::getVG_alpha(int layer) {return VG_alpha[layer]; }
double Soil::getVG_n(int layer) {return VG_n[layer]; }
double Soil::getVG_theta_res(int layer) {return VG_theta_res[layer]; }
double Soil::getVG_theta_sat(int layer) {return VG_theta_sat[layer]; }
std::string Soil::getUSDAType(int layer) {return usda_type[layer]; }
double Soil::getThetaSAT(int layer) {return theta_SAT[layer]; }
double Soil::getThetaFC(int layer) {return theta_FC[layer]; }
double Soil::getPsi(int layer) {return psi[layer]; }
double Soil::getTheta(int layer) {return theta[layer]; }
double Soil::getWaterSAT(int layer) {
  double water_SAT = widths[layer]*theta_SAT[layer]*(1.0-(rfc[layer]/100.0));
  return(water_SAT);
}
double Soil::getWaterFC(int layer) {
  double water_FC = widths[layer]*theta_FC[layer]*(1.0-(rfc[layer]/100.0));
  return(water_FC);
}
double Soil::getTemp(int layer) {return Temp[layer]; }
void Soil::setPsi(int layer, double value) {
  psi[layer] = value; 
  if(model=="VG") {
    theta[layer] = psi2thetaVanGenuchten_c(VG_n[layer], VG_alpha[layer], VG_theta_res[layer], VG_theta_sat[layer], psi[layer]);
  } else {
    theta[layer] = psi2thetaSaxton_c(clay[layer], sand[layer], psi[layer], om[layer]);
  }
  W[layer] = theta[layer]/theta_FC[layer];
}
void Soil::setTheta(int layer, double value) {
  theta[layer] = value; 
  W[layer] = theta[layer]/theta_FC[layer];
  if(model=="VG") {
    psi[layer] = theta2psiVanGenuchten_c(VG_n[layer], VG_alpha[layer], VG_theta_res[layer], VG_theta_sat[layer], theta[layer]);
  } else {
    psi[layer] = theta2psiSaxton_c(clay[layer], sand[layer], theta[layer], om[layer]);
  }
}
void Soil::setW(int layer, double value) {
  W[layer] = value; 
  theta[layer] = W[layer]*theta_FC[layer];
  if(model=="VG") {
    psi[layer] = theta2psiVanGenuchten_c(VG_n[layer], VG_alpha[layer], VG_theta_res[layer], VG_theta_sat[layer], theta[layer]);
  } else {
    psi[layer] = theta2psiSaxton_c(clay[layer], sand[layer], theta[layer], om[layer]);
  }
}
void Soil::setTemp(int layer, double value) {
  Temp[layer] = value; 
}


/*=============================================================================
 * Implementation of soil routines using C++ code
 *=============================================================================*/


std::string USDAType_c(double clay, double sand) {
  double silt = 100 - clay - sand;
  if((silt+1.5*clay)<15) return("Sand");
  else if(((silt+1.5*clay)>=15) && ((silt + 2.0*clay)<30)) return("Loamy sand");
  else if(((clay>=7) && (clay<20) && (sand>52) && ((silt + 2.0*clay)>=30)) || ((clay < 7) && (silt < 50) && ((silt + 2.0*clay)>=30))) return("Sandy loam");
  else if(((clay>=7) && (clay<27)) && ((silt>=28) && (silt<50)) && (sand<=52)) return("Loam");
  else if(((silt>=50) && ((clay>=12) && (clay<27))) || ((silt>=50) && (silt<80) && (clay <12))) return("Silt loam");
  else if((silt>=80) && (clay<12)) return("Silt");
  else if(((clay>=20) && (clay<35)) && (silt<28) && (sand>45)) return("Sandy clay loam");
  else if(((clay>=27) && (clay<40)) && ((sand>20) && (sand<=45))) return("Clay loam");
  else if(((clay>=27) && (clay<40)) && (sand<=20)) return("Silty clay loam");
  else if((clay>=35) && (sand>45)) return("Sandy clay");
  else if((clay>=40) && (silt>=40)) return("Silty clay");
  else if((clay>=40) && (sand<=45) && (silt<40)) return("Clay");
  return("Unknown");
}
/**
 * Saturated conductivity (mmolH20·m-1·s-1·MPa-1)
 */
//' Soil texture and hydraulics
//' 
//' Low-level functions relating soil texture with soil hydraulics and soil water content.
//'
//' @param clay Percentage of clay (in percent weight).
//' @param sand Percentage of sand (in percent weight).
//' @param n,alpha,theta_res,theta_sat Parameters of the Van Genuchten-Mualem model (m = 1 - 1/n).
//' @param psi Water potential (in MPa).
//' @param theta Relative water content (in percent volume).
//' @param om Percentage of organic matter (optional, in percent weight).
//' @param mmol Boolean flag to indicate that saturated conductivity units should be returned in mmol/m/s/MPa. If \code{mmol = FALSE} then units are cm/day.
//' @param bd Bulk density (in g/cm3).
//' @param topsoil A boolean flag to indicate topsoil layer.
//' @param soilType A string indicating the soil type.
//' @param soil Initialized soil object (returned by function \code{\link{soil}}).
//' @param model Either 'SX' or 'VG' for Saxton's or Van Genuchten's water retention models.
//' @param minPsi Minimum water potential (in MPa) to calculate the amount of extractable water.
//' @param pWeight Percentage of corresponding to rocks, in weight.
//' @param bulkDensity Bulk density of the soil fraction (g/cm3).
//' @param rockDensity Rock density (g/cm3).
//' 
//' @details
//' \itemize{
//' \item{\code{soil_psi2thetaSX()} and \code{soil_theta2psiSX()} calculate water potentials (MPa) and water contents (theta) using texture data the formulae of Saxton et al. (1986) or Saxton & Rawls (2006) depending on whether organic matter is available.}
//' \item{\code{soil_psi2thetaVG()} and \code{soil_theta2psiVG()} to the same calculations as before, but using the Van Genuchten - Mualem equations (\enc{Wösten}{Wosten} & van Genuchten 1988). }
//' \item{\code{soil_saturatedConductivitySX()} returns the saturated conductivity of the soil (in mmol/m/s/MPa or cm/day), estimated from formulae of Saxton et al. (1986) or Saxton & Rawls (2006) depending on whether organic matter is available.}
//' \item{\code{soil_unsaturatedConductivitySX()} returns the unsaturated conductivity of the soil (in mmol/m/s/MPa or cm/day), estimated from formulae of Saxton et al. (1986) or Saxton & Rawls (2006) depending on whether organic matter is available.}
//' \item{\code{soil_USDAType()} returns the USDA type (a string) for a given texture.}
//' \item{\code{soil_vanGenuchtenParamsCarsel()} gives parameters for van Genuchten-Mualem equations (alpha, n, theta_res and theta_sat, where alpha is in MPa-1) for a given texture type (Leij et al. 1996) }
//' \item{\code{soil_vanGenuchtenParamsToth()} gives parameters for van Genuchten-Mualem equations (alpha, n, theta_res and theta_sat, where alpha is in MPa-1) for a given texture, organic matter and bulk density (Toth et al. 2015).}
//' \item{\code{soil_psi()} returns the water potential (MPa) of each soil layer, according to its water retention model.}
//' \item{\code{soil_theta()} returns the moisture content (as percent of soil volume) of each soil layer, according to its water retention model.}
//' \item{\code{soil_water()} returns the water volume (mm) of each soil layer, according to its water retention model.}
//' \item{\code{soil_conductivity()} returns the conductivity of each soil layer (in mmol/m/s/MPa or cm/day).}
//' \item{\code{soil_waterExtractable()} returns the water volume (mm) extractable from the soil according to its water retention curves and up to a given soil water potential.}
//' \item{\code{soil_waterFC()} and \code{soil_thetaFC()} calculate the water volume (in mm) and moisture content (as percent of soil volume) of each soil layer at field capacity, respectively.}
//' \item{\code{soil_waterWP()} and \code{soil_thetaWP()} calculate the water volume (in mm) and moisture content (as percent of soil volume) of each soil layer at wilting point (-1.5 MPa), respectively. }
//' \item{\code{soil_waterSAT()}, \code{soil_thetaSATSX()} and \code{soil_thetaSAT()} calculate the saturated water volume (in mm) and moisture content (as percent of soil volume) of each soil layer.}
//' \item{\code{soil_saturatedWaterDepth()} returns the depth to saturation in mm from surface.}
//' \item{\code{soil_rockWeight2Volume()} transforms rock percentage from weight to volume basis.}
//' }
//' 
//' @return Depends on the function (see details).
//' 
//' @references
//' Leij, F.J., Alves, W.J., Genuchten, M.T. Van, Williams, J.R., 1996. The UNSODA Unsaturated Soil Hydraulic Database User’s Manual Version 1.0.
//' 
//' Saxton, K.E., Rawls, W.J., Romberger, J.S., Papendick, R.I., 1986. Estimating generalized soil-water characteristics from texture. Soil Sci. Soc. Am. J. 50, 1031–1036.
//' 
//' Saxton, K.E., Rawls, W.J., 2006. Soil water characteristic estimates by texture and organic matter for hydrologic solutions. Soil Sci. Soc. Am. J. 70, 1569. doi:10.2136/sssaj2005.0117
//' 
//' \enc{Wösten}{Wosten}, J.H.M., & van Genuchten, M.T. 1988. Using texture and other soil properties to predict the unsaturated soil hydraulic functions. Soil Science Society of America Journal 52: 1762–1770.
//' 
//' \enc{Tóth}{Toth}, B., Weynants, M., Nemes, A., \enc{Makó}{Mako}, A., Bilas, G., and \enc{Tóth}{Toth}, G. 2015. New generation of hydraulic pedotransfer functions for Europe. European Journal of Soil Science 66: 226–238.
//' 
//' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
//' 
//' @seealso  \code{\link{soil}}
//' 
//' @examples
//' #Determine USDA soil texture type
//' type = soil_USDAType(clay=40, sand=10)
//' type
//' 
//' #Van Genuchten's params (bulk density = 1.3 g/cm)
//' vg = soil_vanGenuchtenParamsToth(40,10,1,1.3,TRUE)
//' vg
//' 
//' # Define soil with default params
//' soil_df <- defaultSoilParams(4)
//' soil_df
//' 
//' # Initialize soil parameters and state variables
//' s = soil(soil_df)
//' 
//' # Plot Saxton's and Van Genuchten's water retention curves
//' plot(s, model="both")
//' 
//' @name soil_texture
//' @keywords internal
// [[Rcpp::export("soil_saturatedConductivitySX")]]
double saturatedConductivitySaxton_c(double clay, double sand, double bd, double om, bool mmol = true) {
  double Ksat;
  //If organic matter is missing use Saxton et al (1986)
  //Otherwise use Saxton & Rawls (2006)
  if(std::isnan(om)) {
    double theta_sat = 0.332 - 7.251E-4*sand + 0.1276*log10(clay);
    Ksat = 2.778e-6*exp(12.012+-7.55e-2*sand+(-3.8950 + 3.671e-2*sand - 0.1103*clay + 8.7546e-4*pow(clay,2.0))/theta_sat);
    //m/s to cm/day
    Ksat = Ksat*100.0*86400.0;
  } else {
    sand = sand/100.0;
    clay = clay/100.0;
    //om = om/100.0; //OM should be in percentage in Saxton's 2006
    double theta33t = (-0.251*sand) + (0.195*clay) + (0.011*om) + (0.006*(sand*om)) - (0.027*(clay*om)) + (0.452*(sand*clay)) + 0.299;
    double theta33 = theta33t + (1.283*pow(theta33t,2.0) - 0.374 * theta33t - 0.015);
    double theta_S33t = (0.278*sand) + (0.034*clay)+ (0.022*om) - (0.018*(sand*om)) - (0.027*(clay*om)) - (0.584*(sand*clay)) + 0.078;
    double theta_S33 = theta_S33t + (0.636*theta_S33t-0.107);
    double theta_sat = theta33+theta_S33 - (0.097*sand) + 0.043;
    double theta1500t = -0.024*sand + 0.487*clay+0.006*om + 0.005*(sand*om) - 0.013*(clay*om) + 0.068*(sand*clay) + 0.031;
    double theta1500 = theta1500t + (0.14*theta1500t - 0.02);
    double B = 3.816712/(log(theta33)-log(theta1500)); //3.816712 = log(1500) - log(33)
    double lambda = 1.0/B;
    Ksat = 1930.0*pow(theta_sat - theta33, 3.0 - lambda);
    //mm/h to cm/day
    Ksat = Ksat*0.1*24.0;
  }
  //Correct for bulk density
  double bdsoil = 2.73; //Density of soil particles
  double bdref = 1.2; //Reference bulk density for Ksat
  Ksat = Ksat*std::pow((bdsoil - bd)/(bdsoil - bdref),3.0);
  //cm/day to mmolH20·m-1·s-1·MPa-1
  if(mmol) Ksat = Ksat*cmdTOmmolm2sMPa;
  return(Ksat);
}


//' @rdname soil_texture
//' @keywords internal
// [[Rcpp::export("soil_unsaturatedConductivitySX")]]
double unsaturatedConductivitySaxton_c(double theta, double clay, double sand, double bd, double om, bool mmol = true) {
  double Kunsat;
  //If organic matter is missing use Saxton et al (1986)
  //Otherwise use Saxton & Rawls (2006)
  if(std::isnan(om)) {
    Kunsat = 2.778e-6*exp(12.012+-7.55e-2*sand+(-3.8950 + 3.671e-2*sand - 0.1103*clay + 8.7546e-4*pow(clay,2.0))/theta);
    //m/s to cm/day
    Kunsat = Kunsat*100.0*86400.0;
  } else {
    sand = sand/100.0;
    clay = clay/100.0;
    //om = om/100.0; //OM should be in percentage in Saxton's 2006
    double theta33t = (-0.251*sand) + (0.195*clay) + (0.011*om) + (0.006*(sand*om)) - (0.027*(clay*om)) + (0.452*(sand*clay)) + 0.299;
    double theta33 = theta33t + (1.283*pow(theta33t,2.0) - 0.374 * theta33t - 0.015);
    double theta_S33t = (0.278*sand) + (0.034*clay)+ (0.022*om) - (0.018*(sand*om)) - (0.027*(clay*om)) - (0.584*(sand*clay)) + 0.078;
    double theta_S33 = theta_S33t + (0.636*theta_S33t-0.107);
    double theta_sat = theta33+theta_S33 - (0.097*sand) + 0.043;
    double theta1500t = -0.024*sand + 0.487*clay+0.006*om + 0.005*(sand*om) - 0.013*(clay*om) + 0.068*(sand*clay) + 0.031;
    double theta1500 = theta1500t + (0.14*theta1500t - 0.02);
    double B = 3.816712/(log(theta33)-log(theta1500)); //3.816712 = log(1500) - log(33)
    double lambda = 1.0/B;
    double Ksat = 1930.0*pow(theta_sat - theta33, 3.0 - lambda);
    Kunsat = Ksat*pow(theta/theta_sat, 3.0 + (2.0/lambda));
    //mm/h to cm/day
    Kunsat = Kunsat*0.1*24.0;
  }
  //Correct for bulk density
  double bdsoil = 2.73; //Density of soil particles
  double bdref = 1.2; //Reference bulk density for Ksat
  Kunsat = Kunsat*std::pow((bdsoil - bd)/(bdsoil - bdref),3.0);
  //cm/day to mmolH20·m-1·s-1·MPa-1
  if(mmol) Kunsat = Kunsat*cmdTOmmolm2sMPa;
  return(Kunsat);
}

/**
 *  Returns water content (% volume) at saturation according to Saxton's pedotransfer model
 */
//' @rdname soil_texture
//' @keywords internal
// [[Rcpp::export("soil_thetaSATSX")]]
double thetaSATSaxton_c(double clay, double sand, double om) {
  double theta_sat;
  //If organic matter is missing use Saxton et al (1986)
  //Otherwise use Saxton & Rawls (2006)
  if(std::isnan(om)) {
    theta_sat = 0.332 - 7.251E-4*sand + 0.1276*log10(clay);
  } else {
    sand = sand/100.0;
    clay = clay/100.0;
    //om = om/100.0; // om as percentage in Saxton's 2006
    double theta33t = (-0.251*sand) + (0.195*clay) + (0.011*om) + (0.006*(sand*om)) - (0.027*(clay*om)) + (0.452*(sand*clay)) + 0.299;
    double theta33 = theta33t + (1.283*pow(theta33t,2.0) - 0.374 * theta33t - 0.015);
    double theta_S33t = (0.278*sand) + (0.034*clay)+ (0.022*om) - (0.018*(sand*om)) - (0.027*(clay*om)) - (0.584*(sand*clay)) + 0.078;
    double theta_S33 = theta_S33t + (0.636*theta_S33t-0.107);
    theta_sat = theta33+theta_S33 - (0.097*sand) + 0.043;
  }
  return(theta_sat);
}

/**
 * Returns soil water potential (in MPa) according to Saxton's pedotransfer model
 * theta - soil water content (in % volume)
 */
//' @rdname soil_texture
//' @keywords internal
// [[Rcpp::export("soil_theta2psiSX")]]
double theta2psiSaxton_c(double clay, double sand, double theta, double om) {
  double A;
  double B;
  double psi;
  //If organic matter is missing use Saxton et al (1986)
  //Otherwise use Saxton & Rawls (2006)
  if(std::isnan(om)) {
    A = -0.1 * exp(-4.396 - (0.0715*clay)-(0.0004880*pow(sand,2.0)) - (0.00004285*pow(sand,2.0)*clay));
    B = -3.140 - (0.00222*pow(clay,2.0)) - (0.00003484*pow(sand,2.0)*(clay));
    psi = A*pow(theta,B);
    if(psi > -0.01) { // If calculated psi > -10 KPa use linear part
      double theta_sat = thetaSATSaxton_c(clay, sand, om);
      double psi_e = -0.1*(-0.108+(0.341*theta_sat));//air-entry tension in MPa
      double theta_10 = pow(-0.01/A, 1.0/B);//exp((2.302-log(A))/B);
      psi = -0.01 - ((theta-theta_10)*(-0.01 - psi_e)/(theta_sat - theta_10));
      psi = std::min(psi,psi_e); //Truncate to air entry tension
    }
  } else {
    sand = sand/100.0;
    clay = clay/100.0;
    //om = om/100.0;//OM should be in percentage in Saxton's 2006
    double theta1500t = -0.024*sand + 0.487*clay+0.006*om + 0.005*(sand*om) - 0.013*(clay*om) + 0.068*(sand*clay) + 0.031;
    double theta1500 = theta1500t + (0.14*theta1500t - 0.02);
    if(theta1500<0.00001) theta1500 = 0.00001;//Truncate theta1500 to avoid NaN when taking logarithms
    double theta33t = -0.251*sand + 0.195*clay + 0.011*om + 0.006*(sand*om) - 0.027*(clay*om) + 0.452*(sand*clay) + 0.299;
    double theta33 = theta33t + (1.283*pow(theta33t,2.0) - 0.374 * theta33t - 0.015);
    if(theta33<0.00001) theta33 = 0.00001;//Truncate theta33 to avoid NaN when taking logarithms
    B = 3.816712/(log(theta33)-log(theta1500)); //3.816712 = log(1500) - log(33)
    A = exp(3.496508 + B*log(theta33)); // 3.496508 = log(33)
    psi = -0.001*(A*pow(theta,-1.0*B));
    // Rcout<<" "<<theta33<<" "<< theta1500<<" "<<A<<" "<< B <<" "<<psi<<"\n";
    if(psi > fieldCapacityPsi) { // If calculated psi > -33 KPa use linear part
      double theta_S33t = (0.278*sand) + (0.034*clay)+ (0.022*om) - (0.018*(sand*om)) - (0.027*(clay*om)) - (0.584*(sand*clay)) + 0.078;
      double theta_S33 = theta_S33t + (0.636*theta_S33t-0.107);
      double theta_sat = theta33+theta_S33 - (0.097*sand) + 0.043;
      double psi_et = -(21.67*sand) - (27.93*clay) - (81.97*theta_S33)+(71.12*(sand*theta_S33))+(8.29*(clay*theta_S33))+(14.05*(sand*clay))+27.16;
      double psi_e = -0.001*(psi_et + ((0.02*pow(psi_et, 2.0)) - (0.113*psi_et) - 0.70));//air-entry tension in MPa
      // Rcout<<psi_et<<" "<< psi_e<<"\n";
      if(psi_e>0.0) psi_e = 0.0;
      psi = fieldCapacityPsi - ((theta-theta33)*(fieldCapacityPsi - psi_e)/(theta_sat - theta33));
      psi = std::min(psi,psi_e); //Truncate to air entry tension
      // Rcout<<psi<<"\n";
    }
  }
  if(psi < minimumPsi) psi = minimumPsi;
  if(theta==0.0) psi = minimumPsi;
  if(psi>0.0) psi = 0.0;
  return(psi);
}


/**
 *  Returns water content (% volume) according to Saxton's pedotransfer model
 *  psi - Soil water potential (in MPa)
 */
//' @rdname soil_texture
//' @keywords internal
// [[Rcpp::export("soil_psi2thetaSX")]]
double psi2thetaSaxton_c(double clay, double sand, double psi, double om) {
  double A;
  double B;
  double theta;
  //If organic matter is missing use Saxton et al (1986)
  //Otherwise use Saxton & Rawls (2006)
  if(std::isnan(om)) { // less than -10 kPa = -0.01 MPa
    A = -0.1 * exp(-4.396 - (0.0715*clay)-(0.0004880*pow(sand,2.0)) - (0.00004285*pow(sand,2.0)*clay));
    B = -3.140 - (0.00222*pow(clay,2.0)) - (0.00003484*pow(sand,2.0)*(clay));
    if(psi< -0.01) {
      theta = pow(psi/A, 1.0/B);
    } else { //Linear part of the relationship (from -10 kPa to air entry tension)
      double theta_sat = thetaSATSaxton_c(clay, sand, om);
      double psi_e = -0.1*(-0.108+(0.341*theta_sat));//air-entry tension in MPa
      double theta_10 = pow(-0.01/A, 1.0/B);//exp((2.302-log(A))/B);
      psi = std::min(psi,psi_e); //Truncate to air entry tension
      theta = theta_10+(((-0.01-psi)*(theta_sat - theta_10))/(-0.01-psi_e));
    }
  } else {
    sand = sand/100.0;
    clay = clay/100.0;
    //om = om/100.0; // OM should be in percentage in Saxton's 2006
    double theta1500t = (-0.024*sand) + (0.487*clay) + (0.006*om) + (0.005*(sand*om)) - (0.013*(clay*om)) + (0.068*(sand*clay)) + 0.031;
    double theta1500 = theta1500t + ((0.14*theta1500t) - 0.02);
    if(theta1500<0.00001) theta1500 = 0.00001;//Truncate theta1500 to avoid NaN when taking logarithms
    double theta33t = (-0.251*sand) + (0.195*clay) + (0.011*om) + (0.006*(sand*om)) - (0.027*(clay*om)) + (0.452*(sand*clay)) + 0.299;
    double theta33 = theta33t + (1.283*pow(theta33t,2.0) - 0.374 * theta33t - 0.015);
    if(theta33<0.00001) theta33 = 0.00001;//Truncate theta33 to avoid NaN when taking logarithms
    B = 3.816712/(log(theta33)-log(theta1500)); //3.816712 = log(1500) - log(33)
    A = exp(3.496508 + B*log(theta33)); // 3.496508 = log(33)
    // Rcout<<theta1500t<<" "<<theta1500<<" "<<theta33t<<" "<<theta33<<" "<< A<<" "<<B<<" "<< psi<<"\n";
    if(psi< fieldCapacityPsi) {
      psi = psi*(-1000.0);
      theta = pow(psi/A, -1.0/B);
    } else {//Linear part of the relationship (from -10 kPa to air entry tension)
      double theta_S33t = (0.278*sand) + (0.034*clay)+ (0.022*om) - (0.018*(sand*om)) - (0.027*(clay*om)) - (0.584*(sand*clay)) + 0.078;
      double theta_S33 = theta_S33t + (0.636*theta_S33t-0.107);
      double theta_sat = theta33+theta_S33 - (0.097*sand) + 0.043;
      double psi_et = -(21.67*sand) - (27.93*clay) - (81.97*theta_S33)+(71.12*(sand*theta_S33))+(8.29*(clay*theta_S33))+(14.05*(sand*clay))+27.16;
      double psi_e = -0.001*(psi_et + ((0.02*pow(psi_et, 2.0)) - (0.113*psi_et) - 0.70));//air-entry tension in MPa
      if(psi_e>0.0) psi_e = 0.0;
      psi = std::min(psi,psi_e); //Truncate to air entry tension
      theta = theta33+(((fieldCapacityPsi-psi)*(theta_sat - theta33))/(fieldCapacityPsi-psi_e));
    }
  }
  return(theta);
}

/**
 *  Returns soil hydraulic conductivity according to Van Genuchten's pedotransfer model (m = 1 - 1/n)
 */
//' @rdname soil_texture
//' @param ksat saturated hydraulic conductance
//' @keywords internal
// [[Rcpp::export("soil_psi2kVG")]]
double psi2kVanGenuchten_c(double ksat, double n, double alpha, double theta_res, double theta_sat, double psi){
  double m = 1.0 - (1.0/n);
  double Se = pow(1.0 + pow(alpha*std::abs(psi),n),-m);
  double k = ksat*pow(Se,0.5)*pow(1.0 - pow(1.0 - pow(Se, 1.0/m), m), 2.0);
  return(k);
}

double psi2kVanGenuchtenMicropores_c(double k_b, double n, double alpha, double theta_res, double theta_sat, 
                                   double psi, double psi_b){
  double m = 1.0 - (1.0/n);
  double Se = pow(1.0 + pow(alpha*std::abs(psi),n),-m);
  double Se_b = pow(1.0 + pow(alpha*std::abs(psi_b),n),-m);
  //For pressure heads above psi_b, micropore conductivity set to k_b (MACRO 5.0 Larsbo & Jarvis)
  Se = std::min(Se, Se_b);
  double k = k_b*pow(Se/Se_b,0.5)*std::pow(1.0 - std::pow(1.0 - std::pow(Se, 1.0/m), m), 2.0)/std::pow(1.0 - std::pow(1.0 - std::pow(Se_b, 1.0/m), m), 2.0);
  return(k);
}

//From Larsbo et al. (2005) eq. 8
double psi2DVanGenuchten_c(double k_sat, double n, double alpha, double theta_res, double theta_sat, 
                           double psi){
  double m = 1.0 - (1.0/n);
  double l = 0.5;
  double Se = pow(1.0 + pow(alpha*std::abs(psi),n),-m);
  double f_1 = ((1.0 - m)*k_sat)/(alpha*m*(theta_sat - theta_res));
  double f_2 = pow(Se, l - 1.0/m);
  double f_3 = pow(1.0 - pow(Se, 1.0/m),-1.0*m) + pow(1.0 - pow(Se, 1.0/m), m) - 2.0;
  return(f_1*f_2*f_3);
}



//' @rdname soil_texture
//' @keywords internal
// [[Rcpp::export("soil_psi2cVG")]]
double psi2cVanGenuchten_c(double n, double alpha, double theta_res, double theta_sat, double psi){
  double m = 1.0 - (1.0/n);
  double num = alpha*m*n*(theta_sat - theta_res)*pow(alpha*std::abs(psi), n - 1.0);
  double den = pow(1.0 + pow(alpha*std::abs(psi),n),m + 1.0);
  double c = num/den;
  return(c);
}

/**
 *  Returns water content (% volume) according to Van Genuchten's pedotransfer model (m = 1 - 1/n)
 *  psi - Soil water potential (in MPa)
 */
//' @rdname soil_texture
//' @keywords internal
// [[Rcpp::export("soil_psi2thetaVG")]]
double psi2thetaVanGenuchten_c(double n, double alpha, double theta_res, double theta_sat, double psi) {
  double m = 1.0 - (1.0/n);
  double T = pow(1.0  + pow(alpha*std::abs(psi),n),-m);
  return(theta_res+T*(theta_sat-theta_res));
}

/**
 *  Returns  soil water potential (in MPa) according to Van Genuchten's pedotransfer model (m = 1 - 1/n)
 *  theta - soil water content (in % volume)
 */
//' @rdname soil_texture
//' @keywords internal
// [[Rcpp::export("soil_theta2psiVG")]]
double theta2psiVanGenuchten_c(double n, double alpha, double theta_res, double theta_sat, double theta) {
  //if theta > theta_sat then psi = 0
  theta = std::min(theta, theta_sat);
  theta = std::max(theta, theta_res);
  double T = (theta-theta_res)/(theta_sat-theta_res); //content relative
  double m = 1.0 - (1.0/n);
  // double T = pow(pow(alpha*std::abs(psi),n)+1.0,-m);
  double psi = -(1.0/alpha)*pow(pow(T,-1.0/m)-1.0,1.0/n);
  if(psi < minimumPsi) psi = minimumPsi;
  return(psi);
}

