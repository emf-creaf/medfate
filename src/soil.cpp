#include <Rcpp.h>
using namespace Rcpp;

// sand 1.7-2.9 W·m-1·K-1, clay 0.8-6.3 W·m-1·K-1 [Geiger et al. The Climate near the Ground]
const double thermalConductivitySand = 1.57025; //W·m-1·K-1 From Dharssi et al. 2009
const double thermalConductivitySilt = 1.57025; //W·m-1·K-1 From Dharssi et al. 2009
const double thermalConductivityClay = 1.16025; //W·m-1·K-1 From Dharssi et al. 2009
const double thermalConductivityAir = 0.025; //W·m-1·K-1 From Dharssi et al. 2009
const double capacitySand = 1.25*pow(10.0,6.0); //kg·m-3 
const double capacitySilt = 1.19*pow(10.0,6.0); //kg·m-3 
const double capacityClay = 1.23*pow(10.0,6.0); //kg·m-3 

/**
 * Returns Soil water potential (in MPa)
 * theta - soil water content (in % volume)
 */
// [[Rcpp::export("soil.theta2psi")]]
double theta2psi(double clay, double sand, double theta, double om = NA_REAL) {
  double A = NA_REAL;
  double B = NA_REAL;
  double psi = NA_REAL;
  //If organic matter is missing use Saxton et al (1986)
  //Otherwise use Saxton & Rawls (2006)
  if(NumericVector::is_na(om)) {
    A = -0.1 * exp(-4.396 - (0.0715*clay)-(0.0004880*pow(sand,2.0)) - (0.00004285*pow(sand,2.0)*clay));
    B = -3.140 - (0.00222*pow(clay,2.0)) - (0.00003484*pow(sand,2.0)*(clay));
    psi = A*pow(theta,B);
  } else {
    sand = sand/100.0;
    clay = clay/100.0;
    om = om/100.0;
    double theta1500t = -0.024*sand + 0.487*clay+0.006*om + 0.005*(sand*om) - 0.013*(clay*om) + 0.068*(sand*clay) + 0.031;
    double theta1500 = theta1500t + (0.14*theta1500t - 0.02);
    double theta33t = -0.251*sand + 0.195*clay + 0.011*om + 0.006*(sand*om) - 0.027*(clay*om) + 0.452*(sand*clay) + 0.299;
    double theta33 = theta33t + (1.283*pow(theta33t,2.0) - 0.374 * theta33t - 0.015);
    B = 3.816712/(log(theta33)-log(theta1500)); //3.816712 = log(1500) - log(33)
    A = exp(3.496508 + B*log(theta33)); // 3.496508 = log(33)
    psi = -0.001*(A*pow(theta,-1.0*B));
  }
  if(psi < -40.0) psi = -40.0;
  if(theta==0.0) psi = -40.0;
  return(psi);
}
/**
 *  Returns water content (% volume)
 *  psi - Soil water potential (in MPa)
 */
// [[Rcpp::export("soil.psi2theta")]]
double psi2theta(double clay, double sand, double psi, double om = NA_REAL) {
  double A = NA_REAL;
  double B = NA_REAL;
  double theta = NA_REAL;
  //If organic matter is missing use Saxton et al (1986)
  //Otherwise use Saxton & Rawls (2006)
  if(NumericVector::is_na(om)) {
    A = -0.1 * exp(-4.396 - (0.0715*clay)-(0.0004880*pow(sand,2.0)) - (0.00004285*pow(sand,2.0)*clay));
    B = -3.140 - (0.00222*pow(clay,2.0)) - (0.00003484*pow(sand,2.0)*(clay));
    theta = pow(psi/A, 1.0/B);
  } else {
    sand = sand/100.0;
    clay = clay/100.0;
    om = om/100.0;
    double theta1500t = (-0.024*sand) + (0.487*clay) + (0.006*om) + (0.005*(sand*om)) - (0.013*(clay*om)) + (0.068*(sand*clay)) + 0.031;
    double theta1500 = theta1500t + ((0.14*theta1500t) - 0.02);
    double theta33t = (-0.251*sand) + (0.195*clay) + (0.011*om) + (0.006*(sand*om)) - (0.027*(clay*om)) + (0.452*(sand*clay)) + 0.299;
    double theta33 = theta33t + (1.283*pow(theta33t,2.0) - 0.374 * theta33t - 0.015);
    B = 3.816712/(log(theta33)-log(theta1500)); //3.816712 = log(1500) - log(33)
    A = exp(3.496508 + B*log(theta33)); // 3.496508 = log(33)
    psi = psi*(-1000.0);
    // Rcout<<theta1500t<<" "<<theta1500<<" "<<theta33t<<" "<<theta33<<" "<< A<<" "<<B<<" "<< psi<<"\n";
    theta = pow(psi/A, -1.0/B);
  }
  return(theta);
}

// [[Rcpp::export("soil.USDAType")]]
String soilUSDAType(double clay, double sand) {
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
 * Leij, F.J., Alves, W.J., Genuchten, M.T. Van, Williams, J.R., 1996. The UNSODA Unsaturated Soil Hydraulic Database User’s Manual Version 1.0.
 * Textural parameters (1 MPa = 0.00009804139432 cm)
 */
// [[Rcpp::export("soil.vanGenuchtenParams")]]
NumericVector vanGenuchtenParams(String soilType) {
  NumericVector vg(2,NA_REAL);
  if(soilType=="Sand") {vg[0]=1478.967; vg[1]=2.68;}
  else if(soilType=="Loamy sand") {vg[0]=1264.772; vg[1]=2.28;}
  else if(soilType=="Sandy loam") {vg[0]=764.983; vg[1]=1.89;}
  else if(soilType=="Loam") {vg[0]=367.1918; vg[1]=1.56;}
  else if(soilType=="Silt") {vg[0]=163.1964; vg[1]=1.37;}
  else if(soilType=="Silt loam") {vg[0]=203.9955; vg[1]=1.41;}
  else if(soilType=="Sandy clay loam") {vg[0]=601.7866; vg[1]=1.48;}
  else if(soilType=="Clay loam") {vg[0]=193.7957; vg[1]=1.31;}
  else if(soilType=="Silty clay loam") {vg[0]=101.9977; vg[1]=1.23;}
  else if(soilType=="Sandy clay") {vg[0]=275.3939; vg[1]=1.23;}
  else if(soilType=="Silty clay") {vg[0]=50.99887; vg[1]=1.09;}
  else if(soilType=="Clay") {vg[0]=81.59819; vg[1]=1.09;}
  return(vg);
}

/**
 * Soil thermal conductivity 
 *
 * Dharssi, I., Vidale, P.L., Verhoef, A., MacPherson, B., Jones, C., & Best, M. 2009. New soil physical properties implemented in the Unified Model at PS18. 9–12.
 * Best et al. 2011
 */
NumericVector layerthermalconductivity(NumericVector sand, NumericVector clay, NumericVector W, NumericVector Theta_FC) {
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

// [[Rcpp::export("soil.thermalconductivity")]]
NumericVector soilthermalconductivity(List soil) {
  NumericVector sand = soil["sand"];
  NumericVector clay = soil["clay"];
  NumericVector W = soil["W"];
  NumericVector Theta_FC = soil["Theta_FC"];
  return(layerthermalconductivity(sand, clay, W, Theta_FC));
}

/**
 * Soil thermal capacity. Simplified from:
 * 
 *  returns - J·m-3·K-1
 * Cox, P.M., Betts, R.A., Bunton, C.B., Essery, R.L.H., Rowntree, P.R., & Smith, J. 1999. The impact of new land surface physics on the GCM simulation of climate and climate sensitivity. Climate Dynamics 15: 183–203.
 */
NumericVector layerthermalcapacity(NumericVector sand, NumericVector clay, NumericVector W, NumericVector Theta_FC) {
  int nlayers = sand.length();
  NumericVector thermalCap(nlayers,0.0);
  for(int l=0;l<nlayers;l++) {
    thermalCap[l] = ((sand[l]*capacitySand)+(clay[l]*capacityClay) + ((100.0-clay[l]-sand[l])*capacitySilt))/100.0;
    thermalCap[l] = thermalCap[l] + 4.19*pow(10.0,3.0)*1000.0*Theta_FC[l]*W[l];//Add water
  }
  return(thermalCap);
}

// [[Rcpp::export("soil.thermalcapacity")]]
NumericVector soilthermalcapacity(List soil) {
  NumericVector sand = soil["sand"];
  NumericVector clay = soil["clay"];
  NumericVector W = soil["W"];
  NumericVector Theta_FC = soil["Theta_FC"];
  return(layerthermalcapacity(sand, clay, W, Theta_FC));
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
// [[Rcpp::export("soil.temperaturegradient")]]
NumericVector soilTemperatureGradient(NumericVector dVec, NumericVector Temp) {
  NumericVector midZ = midpoints(dVec);
  int nlayers = Temp.length();
  NumericVector gradTemp(nlayers-1,0.0);
  if(nlayers>1) {
    for(int l = 0;l<nlayers-1; l++) {
      gradTemp[l] = (Temp[l+1]-Temp[l])/(0.001*(midZ[l+1]-midZ[l]));
    }
  }
  return(gradTemp);
}


// [[Rcpp::export("soil.temperaturechange")]]
NumericVector soilTemperatureChange(NumericVector dVec, NumericVector Temp,
                                    NumericVector sand, NumericVector clay, 
                                    NumericVector W, NumericVector Theta_FC,
                                    double Gdown) {
  NumericVector lambda = layerthermalconductivity(sand, clay, W, Theta_FC);
  NumericVector Ca = layerthermalcapacity(sand, clay, W, Theta_FC);
  int nlayers = Temp.length();
  NumericVector gradTemp = soilTemperatureGradient(dVec, Temp);
  NumericVector midZ = midpoints(dVec);
  double Gup = -Gdown; //Gdown > 0 when net flux is in the direction of soil
  double Gi;
  NumericVector tempch(nlayers);
  for(int l = 0;l<nlayers; l++) {
    if(l<(nlayers-1)) {
      Gi = lambda[l]*gradTemp[l]; //Gi < 0 when net flux is downward
      tempch[l] = (Gi-Gup)/(Ca[l]*0.001*dVec[l]);
      Gup = Gi;
    } else {
      tempch[l] = (-Gup)/(Ca[l]*0.001*dVec[l]);
    }
  }
  return(tempch);
}

// [[Rcpp::export("soil")]]
List soil(List SoilParams, NumericVector W = NumericVector::create(1.0)) {
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
  //Soil parameters related to texture
  NumericVector clay = clone(as<NumericVector>(SoilParams["clay"]));
  NumericVector sand = clone(as<NumericVector>(SoilParams["sand"]));
  NumericVector om = clone(as<NumericVector>(SoilParams["om"]));
  NumericVector macro = clone(as<NumericVector>(SoilParams["macro"]));
  NumericVector rfc = clone(as<NumericVector>(SoilParams["rfc"]));


  NumericVector temperature(nlayers, NA_REAL);
  NumericVector Theta_FC(nlayers);
  NumericVector psi(nlayers);
  CharacterVector usda_Type(nlayers);
  NumericVector VG_alpha(nlayers);
  NumericVector VG_n(nlayers);
  for(int l=0;l<nlayers;l++) {
    Theta_FC[l] = psi2theta(clay[l], sand[l], -0.033, om[l]); //FC to -33 kPa = -0.033 MPa
    usda_Type[l] = soilUSDAType(clay[l],sand[l]);
    psi[l] = theta2psi(clay[l], sand[l], W[l]*Theta_FC[l], om[l]);
    NumericVector vgl = vanGenuchtenParams(usda_Type[l]);
    VG_alpha[l] = vgl[0];
    VG_n[l] = vgl[1];
    SoilDepth +=dVec[l];
  }
  List l = List::create(_["SoilDepth"] = SoilDepth,
                      _["W"] = W, _["psi"] = psi, _["Temp"] = temperature,
                      _["Ksoil"] = SoilParams["Ksoil"], _["Gsoil"] = SoilParams["Gsoil"],
                      _["dVec"] = dVec,
                      _["sand"] = sand, _["clay"] = clay, _["om"] = om,
                      _["usda_Type"] = usda_Type,
                      _["VG_alpha"] = VG_alpha,_["VG_n"] = VG_n, 
                      _["macro"] = macro, _["rfc"] = rfc,
                      _["Theta_FC"] = Theta_FC);
  l.attr("class") = CharacterVector::create("soil","list");
  return(l);
}


// [[Rcpp::export("soil.waterFC")]]
NumericVector waterFC(List soil) {
  NumericVector dVec = soil["dVec"];
  NumericVector Theta_FC = soil["Theta_FC"];
  NumericVector rfc = soil["rfc"];
  int nlayers = dVec.size();
  NumericVector Water_FC(nlayers);
  for(int i=0;i<nlayers;i++) Water_FC[i] = dVec[i]*Theta_FC[i]*(1.0-(rfc[i]/100.0));
  return(Water_FC);
}