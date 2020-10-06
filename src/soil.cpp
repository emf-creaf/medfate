// [[Rcpp::interfaces(r,cpp)]]

#include <Rcpp.h>
using namespace Rcpp;

// sand 1.7-2.9 W·m-1·K-1, clay 0.8-6.3 W·m-1·K-1 [Geiger et al. The Climate near the Ground]
const double thermalConductivitySand = 1.57025; //W·m-1·K-1 From Dharssi et al. 2009
const double thermalConductivitySilt = 1.57025; //W·m-1·K-1 From Dharssi et al. 2009
const double thermalConductivityClay = 1.16025; //W·m-1·K-1 From Dharssi et al. 2009
const double thermalConductivityAir = 0.025; //W·m-1·K-1 From Dharssi et al. 2009
const double capacitySand = 1.25*1e6; //kg·m-3 
const double capacitySilt = 1.19*1e6; //kg·m-3 
const double capacityClay = 1.23*1e6; //kg·m-3 

/**
 * Conversion factor from conductivity in cm·day-1 to molH20·m-1·MPa-1·s-1
 */
const double cmdTOmmolm2sMPa = 655.2934; //100.0/(18.01528*86400.0*0.00009804139432); 


CharacterVector layerNames(int nlayers) {
  CharacterVector ln(nlayers);
  for(int l=0;l<nlayers;l++){
    char Result[16]; 
    sprintf(Result, "%d", l+1);
    ln[l] = Result;
  }
  return(ln);
}
/**
 * Returs saturated conductivity (mmolH20·m-1·s-1·MPa-1)
 */
// [[Rcpp::export("soil_saturatedConductivitySX")]]
double saturatedConductivitySaxton(double clay, double sand, double om = NA_REAL, bool mmol = true) {
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
    Ksat = 1930.0*pow(theta_sat - theta33, 3.0 - 2.0/B);
    //mm/h to cm/day
    Ksat = Ksat*0.1*24.0;
  }
  //cm/day to mmolH20·m-1·s-1·MPa-1
  if(mmol) Ksat = Ksat*cmdTOmmolm2sMPa;
  return(Ksat);
}
/**
 *  Returns water content (% volume) at saturation according to Saxton's pedotransfer model
 */
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
 *  Returns water content (% volume) according to Van Genuchten's pedotransfer model (m = 1 - 1/n)
 *  psi - Soil water potential (in MPa)
 */
// [[Rcpp::export("soil_psi2thetaVG")]]
double psi2thetaVanGenuchten(double n, double alpha, double theta_res, double theta_sat, double psi) {
  double m = 1.0 - (1.0/n);
  double T = pow(pow(alpha*std::abs(psi),n)+1.0,-m);
  return(theta_res+T*(theta_sat-theta_res));
}
/**
 *  Returns  soil water potential (in MPa) according to Van Genuchten's pedotransfer model (m = 1 - 1/n)
 *  theta - soil water content (in % volume)
 */
// [[Rcpp::export("soil_theta2psiVG")]]
double theta2psiVanGenuchten(double n, double alpha, double theta_res, double theta_sat, double theta) {
  double T = (theta-theta_res)/(theta_sat-theta_res); //content relative
  double m = 1.0 - (1.0/n);
  // double T = pow(pow(alpha*std::abs(psi),n)+1.0,-m);
  double psi = -(1.0/alpha)*pow(pow(T,-1.0/m)-1.0,1.0/n);
  if(psi < -40.0) psi = -40.0;
  return(psi);
}



// [[Rcpp::export("soil_USDAType")]]
String USDAType(double clay, double sand) {
  double silt = 100 - clay - sand;
  if((silt+1.5*clay)<15) return("Sand");
  else if(((silt+1.5*clay)>=15) & ((silt + 2.0*clay)<30)) return("Loamy sand");
  else if(((clay>=7) & (clay<20) & (sand>52) & ((silt + 2.0*clay)>=30)) | ((clay < 7) & (silt < 50) & ((silt + 2.0*clay)>=30))) return("Sandy loam");
  else if(((clay>=7) & (clay<27)) & ((silt>=28) & (silt<50)) & (sand<=52)) return("Loam");
  else if(((silt>=50) & ((clay>=12) & (clay<27))) | ((silt>=50) & (silt<80) & (clay <12))) return("Silt loam");
  else if((silt>=80) & (clay<12)) return("Silt");
  else if(((clay>=20) & (clay<35)) & (silt<28) & (sand>45)) return("Sandy clay loam");
  else if(((clay>=27) & (clay<40)) & ((sand>20) & (sand<=45))) return("Clay loam");
  else if(((clay>=27) & (clay<40)) & (sand<=20)) return("Silty clay loam");
  else if((clay>=35) & (sand>45)) return("Sandy clay");
  else if((clay>=40) & (silt>=40)) return("Silty clay");
  else if((clay>=40) & (sand<=45) &(silt<40)) return("Clay");
  return("Unknown");
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
  
/**
 * Soil thermal conductivity 
 *
 * Dharssi, I., Vidale, P.L., Verhoef, A., MacPherson, B., Jones, C., & Best, M. 2009. New soil physical properties implemented in the Unified Model at PS18. 9–12.
 * Best et al. 2011
 */
NumericVector layerThermalConductivity(NumericVector sand, NumericVector clay, NumericVector W, NumericVector Theta_FC) {
  int nlayers = sand.length();
  NumericVector thermalCond(nlayers,0.0);
  for(int l=0;l<nlayers;l++) {
    double silt = 100 - sand[l] - clay[l];
    double lambda_m = ((thermalConductivitySand*sand[l])+(thermalConductivitySilt*silt)+(thermalConductivityClay*clay[l]))/(silt+sand[l]+clay[l]);
    double lambda_dry = pow(thermalConductivityAir, Theta_FC[l])*pow(lambda_m, (1.0-Theta_FC[l]));
    double Ke = 0.0;
    if(W[l]>=0.1) Ke = log10(W[l]) + 1.0;
    double lambda_s = std::max(1.58,std::min(2.2,1.58 + 12.4*(lambda_dry-0.25)));
    thermalCond[l] = (lambda_s-lambda_dry)*Ke + lambda_dry;
  }
  return(thermalCond);
}


/**
 * Soil thermal capacity. Simplified from:
 * 
 *  returns - J·m-3·K-1
 * Cox, P.M., Betts, R.A., Bunton, C.B., Essery, R.L.H., Rowntree, P.R., & Smith, J. 1999. The impact of new land surface physics on the GCM simulation of climate and climate sensitivity. Climate Dynamics 15: 183–203.
 */
NumericVector layerThermalCapacity(NumericVector sand, NumericVector clay, NumericVector W, NumericVector Theta_FC) {
  int nlayers = sand.length();
  NumericVector thermalCap(nlayers,0.0);
  for(int l=0;l<nlayers;l++) {
    thermalCap[l] = ((sand[l]*capacitySand)+(clay[l]*capacityClay) + ((100.0-clay[l]-sand[l])*capacitySilt))/100.0;
    thermalCap[l] = thermalCap[l] + 4.19*1e3*1000.0*Theta_FC[l]*W[l];//Add water
  }
  return(thermalCap);
}



/**
 * Calculates midpoints of soil layers
 */
NumericVector midpoints(NumericVector dVec) {
  int nlayers = dVec.length();
  double sumz = 0.0;
  NumericVector midZ(nlayers);
  for(int l = 0;l<nlayers; l++) {
    midZ[l] = sumz + dVec[l]/2.0;
    sumz = sumz + dVec[l];
  }
  return(midZ);
}

/**
 * Soil temperature gradient (in ºC/m)
 */
// [[Rcpp::export("soil_temperatureGradient")]]
NumericVector temperatureGradient(NumericVector dVec, NumericVector Temp) {
  NumericVector midZ = midpoints(dVec);
  int nlayers = Temp.length();
  NumericVector gradTemp(nlayers,0.0);
  if(nlayers>1) {
    for(int l = 0;l<nlayers-1; l++) {
      gradTemp[l] = (Temp[l+1]-Temp[l])/(0.001*(midZ[l+1]-midZ[l]));
    }
  }
  gradTemp[nlayers-1] = (15.5-Temp[nlayers-1])/(0.001*(10000.0-midZ[nlayers-1])); //15.5º at 10 m
  return(gradTemp);
}


// [[Rcpp::export("soil_temperatureChange")]]
NumericVector temperatureChange(NumericVector dVec, NumericVector Temp,
                                    NumericVector sand, NumericVector clay, 
                                    NumericVector W, NumericVector Theta_FC,
                                    double Gdown) {
  NumericVector lambda = layerThermalConductivity(sand, clay, W, Theta_FC);
  NumericVector Ca = layerThermalCapacity(sand, clay, W, Theta_FC);
  int nlayers = Temp.length();
  NumericVector gradTemp = temperatureGradient(dVec, Temp);
  NumericVector midZ = midpoints(dVec);
  double Gup = -Gdown; //Gdown > 0 when net flux is in the direction of soil
  double Gi;
  NumericVector tempch(nlayers);
  for(int l = 0;l<nlayers; l++) {
    Gi = lambda[l]*gradTemp[l]; //Gi < 0 when net flux is downward
    tempch[l] = (Gi-Gup)/(Ca[l]*0.001*dVec[l]);
    Gup = Gi;
  }
  return(tempch);
}

// [[Rcpp::export("soil")]]
List soil(DataFrame SoilParams, String VG_PTF = "Toth", 
          NumericVector W = NumericVector::create(1.0), 
          double SWE = 0.0) {
  double SoilDepth = 0.0;
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
  NumericVector om = clone(as<NumericVector>(SoilParams["om"]));
  NumericVector bd = clone(as<NumericVector>(SoilParams["bd"]));
  NumericVector rfc = clone(as<NumericVector>(SoilParams["rfc"]));


  //Parameters to be calculated and state variables
  NumericVector macro(nlayers, NA_REAL);
  NumericVector temperature(nlayers, NA_REAL);
  NumericVector Theta_FC(nlayers);
  CharacterVector usda_Type(nlayers);
  NumericVector VG_alpha(nlayers);
  NumericVector VG_n(nlayers);
  NumericVector VG_theta_res(nlayers);
  NumericVector VG_theta_sat(nlayers);
  NumericVector Ksat(nlayers);
  for(int l=0;l<nlayers;l++) {
    usda_Type[l] = USDAType(clay[l],sand[l]);
    NumericVector vgl;
    if(VG_PTF=="Carsel") {
      vgl = vanGenuchtenParamsCarsel(usda_Type[l]); 
    } else if(VG_PTF=="Toth") {
      if(!SoilParams.containsElementNamed("bd")) stop("bd missing in SoilParams");
      NumericVector bd = as<NumericVector>(SoilParams["bd"]);
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
    Ksat[l] = saturatedConductivitySaxton(clay[l], sand[l], om[l]);
    SoilDepth +=dVec[l];
  }
  // Saturated vertical hydraulic conductivity (mm/day) 
  double Kdrain = 0.05*saturatedConductivitySaxton(clay[nlayers-1], sand[nlayers-1], om[nlayers-1], false);
  double Ksoil = 0.05;
  double Gsoil = 0.5; //TO DO, implement pedotransfer functions for Gsoil
  List l = List::create(_["SoilDepth"] = SoilDepth,
                      _["W"] = W, 
                      _["SWE"] = SWE,
                      _["Temp"] = temperature,
                      _["Ksoil"] = Ksoil, _["Gsoil"] = Gsoil,
                      _["dVec"] = dVec,
                      _["sand"] = sand, _["clay"] = clay, _["om"] = om,
                      _["VG_alpha"] = VG_alpha,_["VG_n"] = VG_n, 
                      _["VG_theta_res"] = VG_theta_res,_["VG_theta_sat"] = VG_theta_sat,
                      _["Ksat"] = Ksat,
                      _["Kdrain"] = Kdrain,
                      _["macro"] = macro,
                      _["bd"] = bd,
                      _["rfc"] = rfc);
  l.attr("class") = CharacterVector::create("soil","list");
  return(l);
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
// [[Rcpp::export("soil_thetaFC")]]
NumericVector thetaFC(List soil, String model="SX") {
  return(psi2thetasoil(soil, -0.033, model));
}


/**
 * Returns water content in volume per soil volume at field capacity, according to the given pedotransfer model
 */
// [[Rcpp::export("soil_thetaWP")]]
NumericVector thetaWP(List soil, String model="SX") {
  return(psi2thetasoil(soil, -1.5, model));
}

/**
 * Returns water content in volume per soil volume at saturation, according to the given pedotransfer model
 */
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
// [[Rcpp::export("soil_theta")]]
NumericVector theta(List soil, String model="SX") {
  NumericVector Theta_FC = thetaFC(soil, model);
  NumericVector W = soil["W"];
  NumericVector Theta = Theta_FC * W;
  return(Theta);
}

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

// [[Rcpp::export("soil_rockWeight2Volume")]]
double rockWeight2Volume(double pWeight, double bulkDensity, double rockDensity = 2.3) {
  double rVolume = pWeight/rockDensity;
  double sVolume = (100.0-pWeight)/bulkDensity;
  return(100.0*rVolume/(rVolume+sVolume));
}
/**
 * Returns current water potential, according to the given pedotransfer model
 */
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


// [[Rcpp::export("soil_waterTableDepth")]]
double waterTableDepth(List soil, String model = "SX") {
  NumericVector dVec = soil["dVec"];
  NumericVector W = soil["W"];
  NumericVector Theta_FC = thetaFC(soil, model);
  NumericVector Theta_SAT = thetaSAT(soil, model);
  int nlayers = W.length();
  double z = 0.0;
  for(int l=0;l<nlayers;l++) {
    if(W[l]>1.0) {
      z = z + dVec[l]*(Theta_SAT[l]-Theta_FC[l]*W[l])/(Theta_SAT[l]-Theta_FC[l]);
    } else {
      z = z + dVec[l];
    }
  }
  return(z);
}

// [[Rcpp::export("soil_thermalCapacity")]]
NumericVector thermalCapacity(List soil, String model = "SX") {
  NumericVector sand = soil["sand"];
  NumericVector clay = soil["clay"];
  NumericVector W = soil["W"];
  NumericVector Theta_FC = thetaFC(soil, model);
  return(layerThermalCapacity(sand, clay, W, Theta_FC));
}

// [[Rcpp::export("soil_thermalConductivity")]]
NumericVector thermalConductivity(List soil, String model = "SX") {
  NumericVector sand = soil["sand"];
  NumericVector clay = soil["clay"];
  NumericVector W = soil["W"];
  NumericVector Theta_FC = thetaFC(soil, model);
  return(layerThermalConductivity(sand, clay, W, Theta_FC));
}
