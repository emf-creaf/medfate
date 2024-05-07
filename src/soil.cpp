// [[Rcpp::interfaces(r,cpp)]]

#include <Rcpp.h>
using namespace Rcpp;

/**
 * Conversion factor from conductivity in cm·day-1 to molH20·m-2·MPa-1·s-1
 *  1 day = 86400 sec
 *  1 mol H20 = 18.01528 g
 */
const double cmdTOmmolm2sMPa = 655.2934; //100.0/(18.01528*86400.0*0.00009804139432); 
/**
 * Conversion factor from cm to MPa
 */
const double cmTOMPa = 0.00009804139432; 

/**
 * Conversion factor from m to MPa
 */
const double mTOMPa = 0.009804139432; //1/9.804139*0.000001; 

CharacterVector layerNames(int nlayers) {
  CharacterVector ln(nlayers);
  for(int l=0;l<nlayers;l++){
    String s("");
    s += (l+1);
    ln[l] = s;
  }
  return(ln);
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
//' @param soil Soil object (returned by function \code{\link{soil}}).
//' @param model Either 'SX' or 'VG' for Saxton's or Van Genuchten's water retention models; or 'both' to plot both retention models.
//' @param minPsi Minimum water potential (in MPa) to calculate the amount of extractable water.
//' @param pWeight Percentage of corresponding to rocks, in weight.
//' @param bulkDensity Bulk density of the soil fraction (g/cm3).
//' @param rockDensity Rock density (g/cm3).
//' 
//' @details
//' \itemize{
//' \item{\code{soil_psi2thetaSX()} and \code{soil_theta2psiSX()} calculate water potentials (MPa) and water contents (theta) using texture data the formulae of Saxton et al. (1986) or Saxton & Rawls (2006) depending on whether organic matter is available.}
//' \item{\code{soil_psi2thetaVG()} and \code{soil_theta2psiVG()} to the same calculations as before, but using the Van Genuchten - Mualem equations (\enc{Wösten}{Wosten} & van Genuchten 1988). }
//' \item{\code{soil_saturatedConductivitySX()} returns the saturated conductivity of the soil (in cm/day or mmol/m/s/MPa), estimated from formulae of Saxton et al. (1986) or Saxton & Rawls (2006) depending on whether organic matter is available.}
//' \item{\code{soil_unsaturatedConductivitySX()} returns the unsaturated conductivity of the soil (in cm/day or mmol/m/s/MPa), estimated from formulae of Saxton et al. (1986) or Saxton & Rawls (2006) depending on whether organic matter is available.}
//' \item{\code{soil_USDAType()} returns the USDA type (a string) for a given texture.}
//' \item{\code{soil_vanGenuchtenParamsCarsel()} gives parameters for van Genuchten-Mualem equations (alpha, n, theta_res and theta_sat, where alpha is in MPa-1) for a given texture type (Leij et al. 1996) }
//' \item{\code{soil_vanGenuchtenParamsToth()} gives parameters for van Genuchten-Mualem equations (alpha, n, theta_res and theta_sat, where alpha is in MPa-1) for a given texture, organic matter and bulk density (Toth et al. 2015).}
//' \item{\code{soil_psi()} returns the water potential (MPa) of each soil layer, according to its water retention model.}
//' \item{\code{soil_theta()} returns the moisture content (as percent of soil volume) of each soil layer, according to its water retention model.}
//' \item{\code{soil_water()} returns the water volume (mm) of each soil layer, according to its water retention model.}
//' \item{\code{soil_conductivity()} returns the conductivity of each soil layer (mmol/m/s/MPa), according the Saxton model.}
//' \item{\code{soil_waterExtractable()} returns the water volume (mm) extractable from the soil according to its water retention curves and up to a given soil water potential.}
//' \item{\code{soil_waterFC()} and \code{soil_thetaFC()} calculate the water volume (in mm) and moisture content (as percent of soil volume) of each soil layer at field capacity, respectively.}
//' \item{\code{soil_waterWP()} and \code{soil_thetaWP()} calculate the water volume (in mm) and moisture content (as percent of soil volume) of each soil layer at wilting point (-1.5 MPa), respectively. }
//' \item{\code{soil_waterSAT()}, \code{soil_thetaSATSX()} and \code{soil_thetaSAT()} calculate the saturated water volume (in mm) and moisture content (as percent of soil volume) of each soil layer.}
//' \item{\code{soil_saturatedWaterDepth()} returns the depth to saturation in mm from surface.}
//' \item{\code{soil_rockWeight2Volume()} transforms rock percentage from weight to volume basis.}
//' \item{\code{soil_retentionCurvePlot()} allows ploting the water retention curve of a given soil layer.}
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
//' # Initialize soil object with default params
//' s = soil(defaultSoilParams())
//' 
//' # Plot Saxton's and Van Genuchten's water retention curves
//' soil_retentionCurvePlot(s, model="both")
//' 
//' @name soil_texture
// [[Rcpp::export("soil_saturatedConductivitySX")]]
double saturatedConductivitySaxton(double clay, double sand, double bd, double om = NA_REAL, bool mmol = true) {
  double Ksat = NA_REAL;
  //If organic matter is missing use Saxton et al (1986)
  //Otherwise use Saxton & Rawls (2006)
  if(NumericVector::is_na(om)) {
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
// [[Rcpp::export("soil_unsaturatedConductivitySX")]]
double unsaturatedConductivitySaxton(double theta, double clay, double sand, double bd, double om = NA_REAL, bool mmol = true) {
  double Kunsat = NA_REAL;
  //If organic matter is missing use Saxton et al (1986)
  //Otherwise use Saxton & Rawls (2006)
  if(NumericVector::is_na(om)) {
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
// [[Rcpp::export("soil_thetaSATSX")]]
double thetaSATSaxton(double clay, double sand, double om = NA_REAL) {
  double theta_sat = NA_REAL;
  //If organic matter is missing use Saxton et al (1986)
  //Otherwise use Saxton & Rawls (2006)
  if(NumericVector::is_na(om)) {
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
// [[Rcpp::export("soil_theta2psiSX")]]
double theta2psiSaxton(double clay, double sand, double theta, double om = NA_REAL) {
  double A = NA_REAL;
  double B = NA_REAL;
  double psi = NA_REAL;
  //If organic matter is missing use Saxton et al (1986)
  //Otherwise use Saxton & Rawls (2006)
  if(NumericVector::is_na(om)) {
    A = -0.1 * exp(-4.396 - (0.0715*clay)-(0.0004880*pow(sand,2.0)) - (0.00004285*pow(sand,2.0)*clay));
    B = -3.140 - (0.00222*pow(clay,2.0)) - (0.00003484*pow(sand,2.0)*(clay));
    psi = A*pow(theta,B);
    if(psi > -0.01) { // If calculated psi > -10 KPa use linear part
      double theta_sat = thetaSATSaxton(clay, sand, om);
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
    if(psi > -0.033) { // If calculated psi > -33 KPa use linear part
      double theta_S33t = (0.278*sand) + (0.034*clay)+ (0.022*om) - (0.018*(sand*om)) - (0.027*(clay*om)) - (0.584*(sand*clay)) + 0.078;
      double theta_S33 = theta_S33t + (0.636*theta_S33t-0.107);
      double theta_sat = theta33+theta_S33 - (0.097*sand) + 0.043;
      double psi_et = -(21.67*sand) - (27.93*clay) - (81.97*theta_S33)+(71.12*(sand*theta_S33))+(8.29*(clay*theta_S33))+(14.05*(sand*clay))+27.16;
      double psi_e = -0.001*(psi_et + ((0.02*pow(psi_et, 2.0)) - (0.113*psi_et) - 0.70));//air-entry tension in MPa
      // Rcout<<psi_et<<" "<< psi_e<<"\n";
      if(psi_e>0.0) psi_e = 0.0;
      psi = -0.033 - ((theta-theta33)*(-0.033 - psi_e)/(theta_sat - theta33));
      psi = std::min(psi,psi_e); //Truncate to air entry tension
      // Rcout<<psi<<"\n";
    }
  }
  if(psi < -40.0) psi = -40.0;
  if(theta==0.0) psi = -40.0;
  if(psi>0.0) psi = 0.0;
  return(psi);
}

/**
 *  Returns water content (% volume) according to Saxton's pedotransfer model
 *  psi - Soil water potential (in MPa)
 */
//' @rdname soil_texture
// [[Rcpp::export("soil_psi2thetaSX")]]
double psi2thetaSaxton(double clay, double sand, double psi, double om = NA_REAL) {
  double A = NA_REAL;
  double B = NA_REAL;
  double theta = NA_REAL;
  //If organic matter is missing use Saxton et al (1986)
  //Otherwise use Saxton & Rawls (2006)
  if(NumericVector::is_na(om)) { // less than -10 kPa = -0.01 MPa
    A = -0.1 * exp(-4.396 - (0.0715*clay)-(0.0004880*pow(sand,2.0)) - (0.00004285*pow(sand,2.0)*clay));
    B = -3.140 - (0.00222*pow(clay,2.0)) - (0.00003484*pow(sand,2.0)*(clay));
    if(psi< -0.01) {
      theta = pow(psi/A, 1.0/B);
    } else { //Linear part of the relationship (from -10 kPa to air entry tension)
      double theta_sat = thetaSATSaxton(clay, sand, om);
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
    if(psi< -0.033) {
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
      theta = theta33+(((-0.033-psi)*(theta_sat - theta33))/(-0.033-psi_e));
    }
  }
  return(theta);
}

/**
 *  Returns soil hydraulic conductivity according to Van Genuchten's pedotransfer model (m = 1 - 1/n)
 */
//' @rdname soil_texture
//' @param ksat saturated hydraulic conductance
// [[Rcpp::export("soil_psi2kVG")]]
double psi2kVanGenuchten(double ksat, double n, double alpha, double theta_res, double theta_sat, double psi){
  double m = 1.0 - (1.0/n);
  double Se = pow(1.0 + pow(alpha*std::abs(psi),n),-m);
  double k = ksat*pow(Se,0.5)*pow(1.0 - pow(1.0 - pow(Se, 1.0/m), m), 2.0);
  return(k);
}

double psi2kVanGenuchtenMicropores(double k_b, double n, double alpha, double theta_res, double theta_sat, 
                                   double psi, double psi_b){
  double m = 1.0 - (1.0/n);
  double Se = pow(1.0 + pow(alpha*std::abs(psi),n),-m);
  double Se_b = pow(1.0 + pow(alpha*std::abs(psi_b),n),-m);
  //For pressure heads above psi_b, micropore conductivity set to k_b (MACRO 5.0 Larsbo & Jarvis)
  Se = std::min(Se, Se_b);
  double k = k_b*pow(Se/Se_b,0.5)*pow(1.0 - pow(1.0 - pow(Se, 1.0/m), m), 2.0)/pow(1.0 - pow(1.0 - pow(Se_b, 1.0/m), m), 2.0);
  return(k);
}

//From Larsbo et al. (2005) eq. 8
double psi2DVanGenuchten(double k_sat, double n, double alpha, double theta_res, double theta_sat, 
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
// [[Rcpp::export("soil_psi2cVG")]]
double psi2cVanGenuchten(double n, double alpha, double theta_res, double theta_sat, double psi){
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
// [[Rcpp::export("soil_psi2thetaVG")]]
double psi2thetaVanGenuchten(double n, double alpha, double theta_res, double theta_sat, double psi) {
  double m = 1.0 - (1.0/n);
  double T = pow(1.0  + pow(alpha*std::abs(psi),n),-m);
  return(theta_res+T*(theta_sat-theta_res));
}

/**
 *  Returns  soil water potential (in MPa) according to Van Genuchten's pedotransfer model (m = 1 - 1/n)
 *  theta - soil water content (in % volume)
 */
//' @rdname soil_texture
// [[Rcpp::export("soil_theta2psiVG")]]
double theta2psiVanGenuchten(double n, double alpha, double theta_res, double theta_sat, double theta) {
  //if theta > theta_sat then psi = 0
  theta = std::min(theta, theta_sat);
  theta = std::max(theta, theta_res);
  double T = (theta-theta_res)/(theta_sat-theta_res); //content relative
  double m = 1.0 - (1.0/n);
  // double T = pow(pow(alpha*std::abs(psi),n)+1.0,-m);
  double psi = -(1.0/alpha)*pow(pow(T,-1.0/m)-1.0,1.0/n);
  if(psi < -40.0) psi = -40.0;
  return(psi);
}


//' @rdname soil_texture
// [[Rcpp::export("soil_USDAType")]]
String USDAType(double clay, double sand) {
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

NumericVector psi2thetasoil(List soil, NumericVector psi, String model="SX") {
  NumericVector SD = soil["dVec"];
  int nlayers = SD.size();
  NumericVector theta(nlayers);
  if(model=="SX") {
    NumericVector clay =soil["clay"];
    NumericVector sand = soil["sand"];
    NumericVector om = soil["om"];
    for(int l=0;l<nlayers;l++) {
      theta[l] = psi2thetaSaxton(clay[l], sand[l], psi[l], om[l]); 
    }
  } else if(model=="VG") {
    NumericVector n =soil["VG_n"];
    NumericVector alpha = soil["VG_alpha"];
    NumericVector theta_res = soil["VG_theta_res"];
    NumericVector theta_sat = soil["VG_theta_sat"];
    for(int l=0;l<nlayers;l++) {
      theta[l] = psi2thetaVanGenuchten(n[l],alpha[l],theta_res[l], theta_sat[l], psi[l]); 
    }
  }
  return(theta);
}

NumericVector psi2thetasoil(List soil, double psi, String model="SX") {
  NumericVector SD = soil["dVec"];
  int nlayers = SD.size();
  NumericVector theta(nlayers);
  if(model=="SX") {
    NumericVector clay =soil["clay"];
    NumericVector sand = soil["sand"];
    NumericVector om = soil["om"];
    for(int l=0;l<nlayers;l++) {
      theta[l] = psi2thetaSaxton(clay[l], sand[l], psi, om[l]); 
    }
  } else if(model=="VG") {
    NumericVector n =soil["VG_n"];
    NumericVector alpha = soil["VG_alpha"];
    NumericVector theta_res = soil["VG_theta_res"];
    NumericVector theta_sat = soil["VG_theta_sat"];
    for(int l=0;l<nlayers;l++) {
      theta[l] = psi2thetaVanGenuchten(n[l],alpha[l],theta_res[l], theta_sat[l], psi); 
    }
  }
  return(theta);
}


/**
 * Returns water content in volume per soil volume at field capacity, according to the given pedotransfer model
 */
//' @rdname soil_texture
// [[Rcpp::export("soil_thetaFC")]]
NumericVector thetaFC(List soil, String model="SX") {
  return(psi2thetasoil(soil, -0.033, model));
}


/**
 * Returns water content in volume per soil volume at field capacity, according to the given pedotransfer model
 */
//' @rdname soil_texture
// [[Rcpp::export("soil_thetaWP")]]
NumericVector thetaWP(List soil, String model="SX") {
  return(psi2thetasoil(soil, -1.5, model));
}

/**
 * Returns water content in volume per soil volume at saturation, according to the given pedotransfer model
 */
//' @rdname soil_texture
// [[Rcpp::export("soil_thetaSAT")]]
NumericVector thetaSAT(List soil, String model="SX") {
  NumericVector SD = soil["dVec"];
  int nlayers = SD.size();
  NumericVector Theta_Sat(nlayers);
  if(model=="SX") {
    NumericVector clay =soil["clay"];
    NumericVector sand = soil["sand"];
    NumericVector om = soil["om"];
    for(int l=0;l<nlayers;l++) {
      Theta_Sat[l] = thetaSATSaxton(clay[l], sand[l], om[l]); 
    }
  } else if(model=="VG") {
    NumericVector theta_sat = soil["VG_theta_sat"];
    for(int l=0;l<nlayers;l++) {
      Theta_Sat[l] = theta_sat[l]; 
    }
  }
  return(Theta_Sat);
}

/**
 * Returns water content in mm at field capacity, according to the given pedotransfer model
 */
//' @rdname soil_texture
// [[Rcpp::export("soil_waterFC")]]
NumericVector waterFC(List soil, String model="SX") {
  NumericVector dVec = soil["dVec"];
  NumericVector Theta_FC = thetaFC(soil, model);
  NumericVector rfc = soil["rfc"];
  int nlayers = dVec.size();
  NumericVector Water_FC(nlayers);
  for(int i=0;i<nlayers;i++) Water_FC[i] = dVec[i]*Theta_FC[i]*(1.0-(rfc[i]/100.0));
  return(Water_FC);
}

/**
 * Returns water content in mm at saturation, according to the given pedotransfer model
 */
//' @rdname soil_texture
// [[Rcpp::export("soil_waterSAT")]]
NumericVector waterSAT(List soil, String model="SX") {
  NumericVector dVec = soil["dVec"];
  NumericVector Theta_SAT = thetaSAT(soil, model);
  NumericVector rfc = soil["rfc"];
  int nlayers = dVec.size();
  NumericVector Water_SAT(nlayers);
  for(int i=0;i<nlayers;i++) Water_SAT[i] = dVec[i]*Theta_SAT[i]*(1.0-(rfc[i]/100.0));
  return(Water_SAT);
}

//' @rdname soil_texture
// [[Rcpp::export("soil_waterWP")]]
NumericVector waterWP(List soil, String model="SX") {
  NumericVector dVec = soil["dVec"];
  NumericVector theta_WP = thetaWP(soil, model);
  NumericVector rfc = soil["rfc"];
  int nlayers = dVec.size();
  NumericVector Water_WP(nlayers);
  for(int i=0;i<nlayers;i++) Water_WP[i] = dVec[i]*theta_WP[i]*(1.0-(rfc[i]/100.0));
  return(Water_WP);
}

//' @rdname soil_texture
 // [[Rcpp::export("soil_waterPsi")]]
 NumericVector waterPsi(List soil, double psi, String model="SX") {
   NumericVector dVec = soil["dVec"];
   NumericVector theta_psi = psi2thetasoil(soil, psi, model);
   NumericVector rfc = soil["rfc"];
   int nlayers = dVec.size();
   NumericVector Water_psi(nlayers);
   for(int i=0;i<nlayers;i++) Water_psi[i] = dVec[i]*theta_psi[i]*(1.0-(rfc[i]/100.0));
   return(Water_psi);
 }

//' @rdname soil_texture
// [[Rcpp::export("soil_waterExtractable")]]
NumericVector waterExtractable(List soil, String model="SX", double minPsi = -5.0) {
  NumericVector dVec = soil["dVec"];
  NumericVector theta_FC = thetaFC(soil, model);
  NumericVector theta_Min = psi2thetasoil(soil, minPsi, model);
  NumericVector rfc = soil["rfc"];
  int nlayers = dVec.size();
  NumericVector Water_Extr(nlayers);
  for(int i=0;i<nlayers;i++) Water_Extr[i] = dVec[i]*((theta_FC[i]- theta_Min[i])*(1.0-(rfc[i]/100.0)));
  return(Water_Extr);
}

/**
 * Returns current water content (in prop. volume), according to the given pedotransfer model
 */
//' @rdname soil_texture
// [[Rcpp::export("soil_theta")]]
NumericVector theta(List soil, String model="SX") {
  NumericVector Theta_FC = thetaFC(soil, model);
  NumericVector W = soil["W"];
  NumericVector Theta = Theta_FC * W;
  return(Theta);
}

//' @rdname soil_texture
// [[Rcpp::export("soil_water")]]
NumericVector water(List soil, String model="SX") {
  NumericVector dVec = soil["dVec"];
  NumericVector Theta = theta(soil, model);
  NumericVector rfc = soil["rfc"];
  int nlayers = dVec.size();
  NumericVector Water(nlayers);
  for(int i=0;i<nlayers;i++) Water[i] = dVec[i]*Theta[i]*(1.0-(rfc[i]/100.0));
  return(Water);
}

//' @rdname soil_texture
// [[Rcpp::export("soil_rockWeight2Volume")]]
double rockWeight2Volume(double pWeight, double bulkDensity, double rockDensity = 2.3) {
  double rVolume = pWeight/rockDensity;
  double sVolume = (100.0-pWeight)/bulkDensity;
  return(100.0*rVolume/(rVolume+sVolume));
}
/**
 * Returns current water potential, according to the given pedotransfer model
 */
//' @rdname soil_texture
// [[Rcpp::export("soil_psi")]]
NumericVector psi(List soil, String model="SX") {
  NumericVector Theta = theta(soil, model);
  int nlayers = Theta.size();
  NumericVector psi(nlayers);
  if(model=="SX") {
    NumericVector clay =soil["clay"];
    NumericVector sand = soil["sand"];
    NumericVector om = soil["om"];
    for(int l=0;l<nlayers;l++) {
      psi[l] = theta2psiSaxton(clay[l], sand[l], Theta[l], om[l]);
    }
  } else if(model=="VG") {
    NumericVector n =soil["VG_n"];
    NumericVector alpha = soil["VG_alpha"];
    NumericVector theta_res = soil["VG_theta_res"];
    NumericVector theta_sat = soil["VG_theta_sat"];
    for(int l=0;l<nlayers;l++) {
      psi[l] = theta2psiVanGenuchten(n[l],alpha[l],theta_res[l], theta_sat[l], Theta[l]); 
    }
  }
  return(psi);
}

//' @rdname soil_texture
// [[Rcpp::export("soil_conductivity")]]
NumericVector conductivity(List soil, String model="SX") {
  NumericVector W = soil["W"];
  int nlayers = W.size();
  NumericVector K(nlayers);
  if(model=="SX") {
    NumericVector Theta = theta(soil, model);
    NumericVector clay =soil["clay"];
    NumericVector sand = soil["sand"];
    NumericVector bd = soil["bd"];
    NumericVector om = soil["om"];
    for(int l=0;l<nlayers;l++) {
      K[l] = unsaturatedConductivitySaxton(Theta[l], clay[l], sand[l], bd[l], om[l]);
    }
  } else {
    NumericVector psiSoil = psi(soil, model);
    NumericVector Ksat = soil["Ksat"];
    for(int l=0;l<nlayers;l++) {
      NumericVector n =soil["VG_n"];
      NumericVector alpha = soil["VG_alpha"];
      NumericVector theta_res = soil["VG_theta_res"];
      NumericVector theta_sat = soil["VG_theta_sat"];
      K[l] = psi2kVanGenuchten(Ksat[l], n[l], alpha[l], theta_res[l], theta_sat[l], psiSoil[l]);
    }
  } 
  return(K);
}

//' @rdname soil_texture
 // [[Rcpp::export("soil_capacitance")]]
 NumericVector capacitance(List soil, String model="SX") {
   NumericVector W = soil["W"];
   int nlayers = W.size();
   NumericVector C(nlayers);
   if(model=="SX") {
     stop("Capacitance not available for model 'SX'");
   } else {
     NumericVector psiSoil = psi(soil, model);
     NumericVector Ksat = soil["Ksat"];
     for(int l=0;l<nlayers;l++) {
       NumericVector n =soil["VG_n"];
       NumericVector alpha = soil["VG_alpha"];
       NumericVector theta_res = soil["VG_theta_res"];
       NumericVector theta_sat = soil["VG_theta_sat"];
       C[l] = psi2cVanGenuchten(n[l], alpha[l], theta_res[l], theta_sat[l], psiSoil[l]);
     }
   } 
   return(C);
 }

//' @rdname soil_texture
// [[Rcpp::export("soil_saturatedWaterDepth")]]
double saturatedWaterDepth(List soil, String model = "SX") {
  NumericVector dVec = soil["dVec"];
  NumericVector W = soil["W"];
  NumericVector Theta_FC = thetaFC(soil, model);
  NumericVector Theta_SAT = thetaSAT(soil, model);
  int nlayers = W.length();
  double z = 0.0;
  int nunsaturated = 0;
  for(int l=0;l<nlayers;l++) {
    if(W[l]>1.0) {
      z = z + dVec[l]*(Theta_SAT[l]-Theta_FC[l]*W[l])/(Theta_SAT[l]-Theta_FC[l]);
    } else {
      z = z + dVec[l];
      nunsaturated++;
    }
  }
  if(nunsaturated==nlayers) z = NA_REAL;
  return(z);
}


/* 
 * Parameters for the Van Genuchten-Mualem equations, taken from:
 * Leij, F.J., Alves, W.J., Genuchten, M.T. Van, Williams, J.R., 1996. The UNSODA Unsaturated Soil Hydraulic Database User’s Manual Version 1.0.
 * after Carsel, R.F., & Parrish, R.S. 1988. Developing joint probability distributions of soil water retention characteristics. Water Resources Research 24: 755–769.
 * 
 * Parameter 'alpha' was transformed from pressure in cm to pressure in MPa
 * Textural parameters (1 cm = 0.00009804139432 MPa)
 * 
 *  0 - alpha
 *  1 - n
 *  2 - residual volumetric water content
 *  3 - saturated water content 
 *  4 - saturated soil conductivity (mmol·m-2·s-1·MPa-1)
 */
//' @rdname soil_texture
// [[Rcpp::export("soil_vanGenuchtenParamsCarsel")]]
NumericVector vanGenuchtenParamsCarsel(String soilType) {
  NumericVector vg(5,NA_REAL);
  if(soilType=="Sand") {vg[0]=1478.967; vg[1]=2.68; vg[2] = 0.045; vg[3]=0.43; vg[4] = 712.80;}
  else if(soilType=="Loamy sand") {vg[0]=1264.772; vg[1]=2.28;vg[2] = 0.057; vg[3]=0.41; vg[4] = 350.16;}
  else if(soilType=="Sandy loam") {vg[0]=764.983; vg[1]=1.89; vg[2] = 0.065; vg[3]=0.41; vg[4] = 106.08;}
  else if(soilType=="Loam") {vg[0]=367.1918; vg[1]=1.56; vg[2] = 0.078; vg[3]=0.43; vg[4] = 24.96;}
  else if(soilType=="Silt") {vg[0]=163.1964; vg[1]=1.37; vg[2] = 0.034; vg[3]=0.46; vg[4] = 6.00;}
  else if(soilType=="Silt loam") {vg[0]=203.9955; vg[1]=1.41; vg[2] = 0.067; vg[3]=0.45; vg[4]=10.80;}
  else if(soilType=="Sandy clay loam") {vg[0]=601.7866; vg[1]=1.48; vg[2] = 0.100; vg[3]=0.39;vg[4]=31.44;}
  else if(soilType=="Clay loam") {vg[0]=193.7957; vg[1]=1.31; vg[2] = 0.095; vg[3]=0.41;vg[4]=6.24;}
  else if(soilType=="Silty clay loam") {vg[0]=101.9977; vg[1]=1.23; vg[2] = 0.089; vg[3]=0.43;vg[4]=1.68;}
  else if(soilType=="Sandy clay") {vg[0]=275.3939; vg[1]=1.23; vg[2] = 0.100; vg[3]=0.38;vg[4]=2.88;}
  else if(soilType=="Silty clay") {vg[0]=50.99887; vg[1]=1.09; vg[2] = 0.070; vg[3]=0.36;vg[4]=0.48;}
  else if(soilType=="Clay") {vg[0]=81.59819; vg[1]=1.09; vg[2] = 0.068; vg[3]=0.38;vg[4]=4.80;}
  vg[4] = vg[4]*cmdTOmmolm2sMPa;
  vg.attr("names") = CharacterVector::create("alpha", "n", "theta_res", "theta_sat", "Ks");
  return(vg);
}

//' @rdname soil_texture
// [[Rcpp::export("soil_campbellParamsClappHornberger")]]
NumericVector campbellParamsClappHornberger(String soilType) {
   NumericVector cp(4,NA_REAL);
   if(soilType=="Sand") {cp[0]=0.395; cp[1]=-12.1; cp[2] = 4.05; cp[3]=63.36;}
   else if(soilType=="Loamy sand") {cp[0]=0.410; cp[1]=-9.1;cp[2] = 4.38; cp[3]=56.28;}
   else if(soilType=="Sandy loam") {cp[0]=0.435; cp[1]=-21.8; cp[2] = 4.90; cp[3]=12.48;}
   else if(soilType=="Silt loam") {cp[0]=0.485; cp[1]=-78.6; cp[2] = 5.30; cp[3]=2.59;}
   else if(soilType=="Loam") {cp[0]=0.451; cp[1]=-47.8; cp[2] = 5.39; cp[3]=2.50;}
   else if(soilType=="Silt") {cp[0]=0.485; cp[1]=-78.6; cp[2] = 5.30; cp[3]=2.59;} // EQUAL TO SILT LOAM
   else if(soilType=="Sandy clay loam") {cp[0]=0.420; cp[1]=-29.9; cp[2] = 7.12; cp[3]=2.27;}
   else if(soilType=="Silty clay loam") {cp[0]=0.477; cp[1]=-35.6; cp[2] = 7.75; cp[3]=0.61;}
   else if(soilType=="Clay loam") {cp[0]=0.476; cp[1]=-63.0; cp[2] = 8.52; cp[3]=0.88;}
   else if(soilType=="Sandy clay") {cp[0]=0.426; cp[1]=-15.3; cp[2] = 10.4; cp[3]=0.38;}
   else if(soilType=="Silty clay") {cp[0]=0.492; cp[1]=-49.0; cp[2] = 10.4; cp[3]=0.37;}
   else if(soilType=="Clay") {cp[0]=0.482; cp[1]=-40.5; cp[2] = 11.4; cp[3]=0.46;}
   cp.attr("names") = CharacterVector::create("theta_sat", "psi_sat_cm", "b", "K_sat_cm_h");
   return(cp);
 }
/* 
 * Parameters for the Van Genuchten-Mualem equations, taken from:
 * Tóth, B., Weynants, M., Nemes, A., Makó, A., Bilas, G., & Tóth, G. 2015. New generation of hydraulic pedotransfer functions for Europe. European Journal of Soil Science 66: 226–238.
 * Parameter 'alpha' was transformed from pressure in cm to pressure in MPa
 * Textural parameters (1 MPa = 0.00009804139432 cm)
 * Model #21
 * 
 *  0 - alpha
 *  1 - n
 *  2 - residual volumetric water content
 *  3 - saturated water content 
 */
//' @rdname soil_texture
// [[Rcpp::export("soil_vanGenuchtenParamsToth")]]
NumericVector vanGenuchtenParamsToth(double clay, double sand, double om, double bd, bool topsoil) {
  double silt = 100.0 - clay - sand;
  double ts = 1.0;
  if(!topsoil) ts = 0.0;
  if(NumericVector::is_na(om)) om = 0.0;
  NumericVector vg(4,NA_REAL);
  //Theta_res
  if(sand>=2.0) vg[2] = 0.041;
  else vg[2] = 0.179;
  //Theta_sat
  vg[3] = 0.83080 - 0.28217*bd+0.0002728*clay + 0.000187*silt; 
  //Alpha
  vg[0] = (1.0/0.00009804139432)*pow(10.0,(-0.43348 - 0.41729*bd - 0.04762*om+0.21810*ts - 0.01582*clay - 0.01207*silt));
  //N
  vg[1] = 1.0 + pow(10.0, 0.22236 - 0.30189*bd - 0.05558*ts - 0.005306*clay - 0.003084*silt - 0.01072*om);
  vg.attr("names") = CharacterVector::create("alpha", "n", "theta_res", "theta_sat");
  return(vg);
}


//' Soil initialization
//'
//' Initializes soil parameters and state variables for its use in simulations.
//' 
//' @param SoilParams A data frame of soil parameters (see an example in \code{\link{defaultSoilParams}}).
//' @param VG_PTF Pedotransfer functions to obtain parameters for the van Genuchten-Mualem equations. Either \code{"Carsel"} (Carsel and Parrish 1988) or \code{"Toth"} (Toth et al. 2015).
//' @param W A numerical vector with the initial relative water content of each soil layer.
//' @param SWE Initial snow water equivalent of the snow pack on the soil surface (mm).
//' 
//' @return
//' Function \code{soil} returns a list of class \code{soil} with the following elements:
//' \itemize{
//'   \item{\code{W}: State variable with relative water content of each layer (in as proportion relative to FC).}
//'   \item{\code{SWE}: Initial snow water equivalent of the snow pack on the soil surface (mm).}
//'   \item{\code{Temp}: State variable with temperature (in ºC) of each layer.}
//'   \item{\code{Gsoil}: Gamma parameter for bare soil evaporation (see \code{\link{hydrology_soilEvaporationAmount}}).}
//'   \item{\code{dVec}: Width of soil layers (in mm).}
//'   \item{\code{sand}: Sand percentage for each layer (in percent volume).}
//'   \item{\code{clay}: Clay percentage for each layer (in percent volume).}
//'   \item{\code{om}: Organic matter percentage for each layer (in percent volume).}
//'   \item{\code{nitrogen}: Sum of total nitrogen (ammonia, organic and reduced nitrogen) for each layer (in g/kg).}
//'   \item{\code{VG_alpha}, \code{VG_n}, \code{VG_theta_res}, \code{VG_theta_sat}: Parameters for van Genuchten's pedotransfer functions, for each layer, corresponding to the USDA texture type.}
//'   \item{\code{Ksat}: Saturated soil conductivity for each layer (estimated using function \code{\link{soil_saturatedConductivitySX}}.}
//'   \item{\code{macro}: Macroporosity for each layer (estimated using Stolf et al. 2011).}
//'   \item{\code{rfc}: Percentage of rock fragment content for each layer.}
//' }
//' 
//' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
//' 
//' @details 
//' Function \code{print} prompts a description of soil characteristics and state variables (water content and temperature) 
//' according to a water retention curve (either Saxton's or Van Genuchten's). 
//' Volume at field capacity is calculated assuming a soil water potential equal to -0.033 MPa. 
//' Parameter \code{Temp} is initialized as missing for all soil layers. 
//' 
//' If available, the user can specify columns \code{VG_alpha}, \code{VG_n}, \code{VG_theta_res}, \code{VG_theta_sat} and \code{K_sat},
//' to override Van Genuchten parameters an saturated conductivity estimated from pedotransfer functions when calling function \code{soil}. 
//' 
//' @references
//' Carsel, R.F., and Parrish, R.S. 1988. Developing joint probability distributions of soil water retention characteristics. Water Resources Research 24: 755–769.
//' 
//' \enc{Tóth}{Toth}, B., Weynants, M., Nemes, A., \enc{Makó}{Mako}, A., Bilas, G., and \enc{Tóth}{Toth}, G. 2015. New generation of hydraulic pedotransfer functions for Europe. European Journal of Soil Science 66: 226–238.
//' 
//' Stolf, R., Thurler, A., Oliveira, O., Bacchi, S., Reichardt, K., 2011. Method to estimate soil macroporosity and microporosity based on sand content and bulk density. Rev. Bras. Ciencias do Solo 35, 447–459.
//' 
//' @seealso   \code{\link{soil_redefineLayers}}, \code{\link{soil_psi2thetaSX}}, \code{\link{soil_psi2thetaVG}}, \code{\link{spwb}}, \code{\link{defaultSoilParams}}
//' 
//' @examples
//' # Default parameters
//' df_soil <- defaultSoilParams()
//' 
//' # Initializes soil
//' s = soil(df_soil)
//' 
//' # Prints soil characteristics according to Saxton's water retention curve
//' print(s, model="SX")
//' 
//' # Prints soil characteristics according to Van Genuchten's water retention curve
//' print(s, model="VG")
//' 
//' # Add columns 'VG_theta_sat' and 'VG_theta_res' with custom values
//' df_soil$VG_theta_sat <- 0.400 
//' df_soil$VG_theta_res <- 0.040 
//' 
//' # Reinitialize soil (should override estimations)
//' s2 = soil(df_soil)
//' print(s2, model="VG")
//' @name soil
// [[Rcpp::export("soil")]]
List soil(DataFrame SoilParams, String VG_PTF = "Toth", 
          NumericVector W = NumericVector::create(1.0), 
          double SWE = 0.0) {
  NumericVector dVec = clone(as<NumericVector>(SoilParams["widths"]));
  int nlayers = dVec.size();

  if(W.size()==1) {
    double w0 = W[0];
    W = NumericVector(nlayers);
    for(int l=0;l<nlayers;l++) W[l] = w0; 
  } else {
    W = clone(W);
  }
  
  //Soil parameters related to physical structure
  NumericVector clay = clone(as<NumericVector>(SoilParams["clay"]));
  NumericVector sand = clone(as<NumericVector>(SoilParams["sand"]));
  NumericVector bd = clone(as<NumericVector>(SoilParams["bd"]));
  NumericVector rfc = clone(as<NumericVector>(SoilParams["rfc"]));
  
  if(any(is_na(clay))) stop("Missing values in soil 'clay'");
  if(any(is_na(sand))) stop("Missing values in soil 'sand'");
  if(any(is_na(bd))) stop("Missing values in soil 'bd'");
  if(any(is_na(rfc))) stop("Missing values in soil 'rfc'");

  //Optional
  NumericVector om(nlayers, NA_REAL);
  NumericVector nitrogen(nlayers, NA_REAL);
  if(SoilParams.containsElementNamed("om")) om = clone(as<NumericVector>(SoilParams["om"]));
  if(SoilParams.containsElementNamed("nitrogen")) nitrogen = clone(as<NumericVector>(SoilParams["nitrogen"]));
  
  //Parameters to be calculated and state variables
  NumericVector macro(nlayers, NA_REAL);
  NumericVector temperature(nlayers, NA_REAL);
  CharacterVector usda_Type(nlayers);
  NumericVector VG_alpha(nlayers, NA_REAL);
  NumericVector VG_n(nlayers, NA_REAL);
  NumericVector VG_theta_res(nlayers, NA_REAL);
  NumericVector VG_theta_sat(nlayers, NA_REAL);
  NumericVector Ksat(nlayers, NA_REAL);
  
  //Get parameters from input if specified
  if(SoilParams.containsElementNamed("VG_alpha")) {
    VG_alpha = clone(as<NumericVector>(SoilParams["VG_alpha"]));
  }
  if(SoilParams.containsElementNamed("VG_n")) {
    VG_n = clone(as<NumericVector>(SoilParams["VG_n"]));
  }
  if(SoilParams.containsElementNamed("VG_theta_res")) {
    VG_theta_res = clone(as<NumericVector>(SoilParams["VG_theta_res"]));
  }
  if(SoilParams.containsElementNamed("VG_theta_sat")) {
    VG_theta_sat = clone(as<NumericVector>(SoilParams["VG_theta_sat"]));
  }
  if(SoilParams.containsElementNamed("Ksat")) {
    Ksat = clone(as<NumericVector>(SoilParams["Ksat"]));
  }
  for(int l=0;l<nlayers;l++) {
    usda_Type[l] = USDAType(clay[l],sand[l]);
    NumericVector vgl;
    if(VG_PTF=="Carsel") {
      if(NumericVector::is_na(Ksat[l])) {
        vgl = vanGenuchtenParamsCarsel(usda_Type[l]); 
        Ksat[l] = vgl[4]; //Use Carsel estimate for Ksat
      }
    } else if(VG_PTF=="Toth") {
      if(NumericVector::is_na(Ksat[l])) {
        if(!SoilParams.containsElementNamed("bd")) stop("bd missing in SoilParams");
        NumericVector bd = as<NumericVector>(SoilParams["bd"]);
        // if(l==0) vgl = vanGenuchtenParamsToth(clay[l], sand[l], om[l], bd[l], TRUE);
        //Use non-top soil equation for all layers
        vgl = vanGenuchtenParamsToth(clay[l], sand[l], om[l], bd[l], FALSE);
        // Stolf, R., Thurler, A., Oliveira, O., Bacchi, S., Reichardt, K., 2011. Method to estimate soil macroporosity and microporosity based on sand content and bulk density. Rev. Bras. Ciencias do Solo 35, 447–459.
        macro[l] = std::max(0.0,0.693 - 0.465*bd[l] + 0.212*(sand[l]/100.0));
        Ksat[l] = saturatedConductivitySaxton(clay[l], sand[l], bd[l], om[l]); 
      }
    } else {
      stop("Wrong value for 'VG_PTF'");
    }
    if(NumericVector::is_na(VG_alpha[l])) VG_alpha[l] = vgl[0];
    if(NumericVector::is_na(VG_n[l])) VG_n[l] = vgl[1];
    if(NumericVector::is_na(VG_theta_res[l])) VG_theta_res[l] = vgl[2];
    if(NumericVector::is_na(VG_theta_sat[l])) VG_theta_sat[l] = vgl[3];
  }
  double Gsoil = 0.5; //TO DO, implement pedotransfer functions for Gsoil
  List l = List::create(_["W"] = W, 
                       _["SWE"] = SWE,
                       _["Temp"] = temperature,
                       _["Gsoil"] = Gsoil,
                       _["dVec"] = dVec,
                       _["sand"] = sand, _["clay"] = clay, _["om"] = om, _["nitrogen"] = nitrogen,
                       _["VG_alpha"] = VG_alpha,_["VG_n"] = VG_n, 
                       _["VG_theta_res"] = VG_theta_res,_["VG_theta_sat"] = VG_theta_sat,
                       _["Ksat"] = Ksat,
                       _["macro"] = macro,
                       _["bd"] = bd,
                       _["rfc"] = rfc);
  l.attr("class") = CharacterVector::create("soil","list");
  return(l);
}


// [[Rcpp::export(".modifySoilLayerParam")]]
void modifySoilLayerParam(List soil, String paramName, int layer, double newValue, 
                          String VG_PTF = "Toth") {
  
  //Perform modification
  NumericVector paramVec = as<NumericVector>(soil[paramName]);
  paramVec[layer] = newValue;
  
  //Recalculate necessary soil parameters
  NumericVector clay = as<NumericVector>(soil["clay"]);
  NumericVector sand = as<NumericVector>(soil["sand"]);
  NumericVector om = as<NumericVector>(soil["om"]);
  NumericVector bd = as<NumericVector>(soil["bd"]);
  NumericVector rfc =as<NumericVector>(soil["rfc"]);
  NumericVector dVec = as<NumericVector>(soil["dVec"]);
  NumericVector macro = as<NumericVector>(soil["macro"]);
  NumericVector VG_alpha = as<NumericVector>(soil["VG_alpha"]);
  NumericVector VG_n = as<NumericVector>(soil["VG_n"]);
  NumericVector VG_theta_res = as<NumericVector>(soil["VG_theta_res"]);
  NumericVector VG_theta_sat = as<NumericVector>(soil["VG_theta_sat"]);
  NumericVector Ksat = as<NumericVector>(soil["Ksat"]);
  
  int nlayers = dVec.size();

  //Parameters to be re-calculated
  CharacterVector usda_Type(nlayers);
  for(int l=0;l<nlayers;l++) {
    usda_Type[l] = USDAType(clay[l],sand[l]);
    NumericVector vgl;
    if(VG_PTF=="Carsel") {
      vgl = vanGenuchtenParamsCarsel(usda_Type[l]); 
    } else if(VG_PTF=="Toth") {
      if(l==0) vgl = vanGenuchtenParamsToth(clay[l], sand[l], om[l], bd[l], TRUE);
      else vgl = vanGenuchtenParamsToth(clay[l], sand[l], om[l], bd[l], FALSE);
    } else {
      stop("Wrong value for 'VG_PTF'");
    }
    VG_alpha[l] = vgl[0];
    VG_n[l] = vgl[1];
    VG_theta_res[l] = vgl[2];
    VG_theta_sat[l] = vgl[3];
    // Stolf, R., Thurler, A., Oliveira, O., Bacchi, S., Reichardt, K., 2011. Method to estimate soil macroporosity and microporosity based on sand content and bulk density. Rev. Bras. Ciencias do Solo 35, 447–459.
    macro[l] = std::max(0.0,0.693 - 0.465*bd[l] + 0.212*(sand[l]/100.0));
    Ksat[l] = saturatedConductivitySaxton(clay[l], sand[l], bd[l], om[l]);
  }
}