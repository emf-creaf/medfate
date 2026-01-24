// [[Rcpp::interfaces(r,cpp)]]

#include <Rcpp.h>
#include <string>
#include <vector>
#include "soil.h"
#include "soil_c.h"
using namespace Rcpp;

CharacterVector layerNames(int nlayers) {
  CharacterVector ln(nlayers);
  for(int l=0;l<nlayers;l++){
    String s("");
    s += (l+1);
    ln[l] = s;
  }
  return(ln);
}


//' @rdname soil_texture
//' @keywords internal
// [[Rcpp::export("soil_USDAType")]]
String USDAType(double clay, double sand) {
  return(String(USDAType_c(clay, sand)));
}

NumericVector psi2thetasoil(DataFrame soil, NumericVector psi, String model="SX") {
  NumericVector SD = soil["widths"];
  int nlayers = SD.size();
  NumericVector theta(nlayers);
  if(model=="SX") {
    NumericVector clay =soil["clay"];
    NumericVector sand = soil["sand"];
    NumericVector om = soil["om"];
    for(int l=0;l<nlayers;l++) {
      theta[l] = psi2thetaSaxton_c(clay[l], sand[l], psi[l], om[l]); 
    }
  } else if(model=="VG") {
    NumericVector n =soil["VG_n"];
    NumericVector alpha = soil["VG_alpha"];
    NumericVector theta_res = soil["VG_theta_res"];
    NumericVector theta_sat = soil["VG_theta_sat"];
    for(int l=0;l<nlayers;l++) {
      theta[l] = psi2thetaVanGenuchten_c(n[l],alpha[l],theta_res[l], theta_sat[l], psi[l]); 
    }
  }
  return(theta);
}

NumericVector psi2thetasoil(DataFrame soil, double psi, String model="SX") {
  NumericVector SD = soil["widths"];
  int nlayers = SD.size();
  NumericVector theta(nlayers);
  if(model=="SX") {
    NumericVector clay =soil["clay"];
    NumericVector sand = soil["sand"];
    NumericVector om = soil["om"];
    for(int l=0;l<nlayers;l++) {
      theta[l] = psi2thetaSaxton_c(clay[l], sand[l], psi, om[l]); 
    }
  } else if(model=="VG") {
    NumericVector n =soil["VG_n"];
    NumericVector alpha = soil["VG_alpha"];
    NumericVector theta_res = soil["VG_theta_res"];
    NumericVector theta_sat = soil["VG_theta_sat"];
    for(int l=0;l<nlayers;l++) {
      theta[l] = psi2thetaVanGenuchten_c(n[l],alpha[l],theta_res[l], theta_sat[l], psi); 
    }
  }
  return(theta);
}


/**
 * Returns water content in volume per soil volume at field capacity, according to the given pedotransfer model
 */
//' @rdname soil_texture
//' @keywords internal
// [[Rcpp::export("soil_thetaFC")]]
NumericVector thetaFC(DataFrame soil, String model="SX") {
  if(!soil.inherits("soil")) {
    if(soil.inherits("data.frame")) stop("Please, initialize soil parameters using function `soil()`");
    else stop("Wrong class for `soil`.");
  }
  return(psi2thetasoil(soil, fieldCapacityPsi, model));
}


/**
 * Returns water content in volume per soil volume at field capacity, according to the given pedotransfer model
 */
//' @rdname soil_texture
//' @keywords internal
// [[Rcpp::export("soil_thetaWP")]]
NumericVector thetaWP(DataFrame soil, String model="SX") {
  if(!soil.inherits("soil")) {
    if(soil.inherits("data.frame")) stop("Please, initialize soil parameters using function `soil()`");
    else stop("Wrong class for `soil`.");
  }
  return(psi2thetasoil(soil, -1.5, model));
}

/**
 * Returns water content in volume per soil volume at saturation, according to the given pedotransfer model
 */
//' @rdname soil_texture
// [[Rcpp::export("soil_thetaSAT")]]
NumericVector thetaSAT(DataFrame soil, String model="SX") {
  if(!soil.inherits("soil")) {
    if(soil.inherits("data.frame")) stop("Please, initialize soil parameters using function `soil()`");
    else stop("Wrong class for `soil`.");
  }
  NumericVector SD = soil["widths"];
  int nlayers = SD.size();
  NumericVector Theta_Sat(nlayers);
  if(model=="SX") {
    NumericVector clay =soil["clay"];
    NumericVector sand = soil["sand"];
    NumericVector om = soil["om"];
    for(int l=0;l<nlayers;l++) {
      Theta_Sat[l] = thetaSATSaxton_c(clay[l], sand[l], om[l]); 
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
NumericVector waterFC(DataFrame soil, String model="SX") {
  if(!soil.inherits("soil")) {
    if(soil.inherits("data.frame")) stop("Please, initialize soil parameters using function `soil()`");
    else stop("Wrong class for `soil`.");
  }
  NumericVector widths = soil["widths"];
  NumericVector Theta_FC = thetaFC(soil, model);
  NumericVector rfc = soil["rfc"];
  int nlayers = widths.size();
  NumericVector Water_FC(nlayers);
  for(int i=0;i<nlayers;i++) Water_FC[i] = widths[i]*Theta_FC[i]*(1.0-(rfc[i]/100.0));
  return(Water_FC);
}

/**
 * Returns water content in mm at saturation, according to the given pedotransfer model
 */
//' @rdname soil_texture
//' @keywords internal
// [[Rcpp::export("soil_waterSAT")]]
NumericVector waterSAT(DataFrame soil, String model="SX") {
  if(!soil.inherits("soil")) {
    if(soil.inherits("data.frame")) stop("Please, initialize soil parameters using function `soil()`");
    else stop("Wrong class for `soil`.");
  }
  NumericVector widths = soil["widths"];
  NumericVector Theta_SAT = thetaSAT(soil, model);
  NumericVector rfc = soil["rfc"];
  int nlayers = widths.size();
  NumericVector Water_SAT(nlayers);
  for(int i=0;i<nlayers;i++) Water_SAT[i] = widths[i]*Theta_SAT[i]*(1.0-(rfc[i]/100.0));
  return(Water_SAT);
}

//' @rdname soil_texture
//' @keywords internal
// [[Rcpp::export("soil_waterWP")]]
NumericVector waterWP(DataFrame soil, String model="SX") {
  if(!soil.inherits("soil")) {
    if(soil.inherits("data.frame")) stop("Please, initialize soil parameters using function `soil()`");
    else stop("Wrong class for `soil`.");
  }
  NumericVector widths = soil["widths"];
  NumericVector theta_WP = thetaWP(soil, model);
  NumericVector rfc = soil["rfc"];
  int nlayers = widths.size();
  NumericVector Water_WP(nlayers);
  for(int i=0;i<nlayers;i++) Water_WP[i] = widths[i]*theta_WP[i]*(1.0-(rfc[i]/100.0));
  return(Water_WP);
}

//' @rdname soil_texture
//' @keywords internal
// [[Rcpp::export("soil_waterPsi")]]
 NumericVector waterPsi(DataFrame soil, double psi, String model="SX") {
   if(!soil.inherits("soil")) {
     if(soil.inherits("data.frame")) stop("Please, initialize soil parameters using function `soil()`");
     else stop("Wrong class for `soil`.");
   }
   NumericVector widths = soil["widths"];
   NumericVector theta_psi = psi2thetasoil(soil, psi, model);
   NumericVector rfc = soil["rfc"];
   int nlayers = widths.size();
   NumericVector Water_psi(nlayers);
   for(int i=0;i<nlayers;i++) Water_psi[i] = widths[i]*theta_psi[i]*(1.0-(rfc[i]/100.0));
   return(Water_psi);
 }

//' @rdname soil_texture
//' @keywords internal
// [[Rcpp::export("soil_waterExtractable")]]
NumericVector waterExtractable(DataFrame soil, String model="SX", double minPsi = -5.0) {
  if(!soil.inherits("soil")) {
    if(soil.inherits("data.frame")) stop("Please, initialize soil parameters using function `soil()`");
    else stop("Wrong class for `soil`.");
  }
  NumericVector widths = soil["widths"];
  NumericVector theta_FC = thetaFC(soil, model);
  NumericVector theta_Min = psi2thetasoil(soil, minPsi, model);
  NumericVector rfc = soil["rfc"];
  int nlayers = widths.size();
  NumericVector Water_Extr(nlayers);
  for(int i=0;i<nlayers;i++) Water_Extr[i] = widths[i]*((theta_FC[i]- theta_Min[i])*(1.0-(rfc[i]/100.0)));
  return(Water_Extr);
}

/**
 * Returns current water content (in prop. volume), according to the given pedotransfer model
 */
//' @rdname soil_texture
//' @keywords internal
// [[Rcpp::export("soil_theta")]]
NumericVector theta(DataFrame soil, String model="SX") {
  if(!soil.inherits("soil")) {
    if(soil.inherits("data.frame")) stop("Please, initialize soil parameters using function `soil()`");
    else stop("Wrong class for `soil`.");
  }
  NumericVector Theta_FC = thetaFC(soil, model);
  NumericVector W = soil["W"];
  NumericVector Theta = Theta_FC * W;
  return(Theta);
}

//' @rdname soil_texture
//' @keywords internal
// [[Rcpp::export("soil_water")]]
NumericVector water(DataFrame soil, String model="SX") {
  if(!soil.inherits("soil")) {
    if(soil.inherits("data.frame")) stop("Please, initialize soil parameters using function `soil()`");
    else stop("Wrong class for `soil`.");
  }
  NumericVector widths = soil["widths"];
  NumericVector Theta = theta(soil, model);
  NumericVector rfc = soil["rfc"];
  int nlayers = widths.size();
  NumericVector Water(nlayers);
  for(int i=0;i<nlayers;i++) Water[i] = widths[i]*Theta[i]*(1.0-(rfc[i]/100.0));
  return(Water);
}

//' @rdname soil_texture
//' @keywords internal
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
//' @keywords internal
// [[Rcpp::export("soil_psi")]]
NumericVector psi(DataFrame soil, String model="SX") {
  NumericVector Theta = theta(soil, model);
  int nlayers = Theta.size();
  NumericVector psi(nlayers);
  if(model=="SX") {
    NumericVector clay =soil["clay"];
    NumericVector sand = soil["sand"];
    NumericVector om = soil["om"];
    for(int l=0;l<nlayers;l++) {
      psi[l] = theta2psiSaxton_c(clay[l], sand[l], Theta[l], om[l]);
    }
  } else if(model=="VG") {
    NumericVector n =soil["VG_n"];
    NumericVector alpha = soil["VG_alpha"];
    NumericVector theta_res = soil["VG_theta_res"];
    NumericVector theta_sat = soil["VG_theta_sat"];
    for(int l=0;l<nlayers;l++) {
      psi[l] = theta2psiVanGenuchten_c(n[l],alpha[l],theta_res[l], theta_sat[l], Theta[l]); 
    }
  }
  return(psi);
}

//' @rdname soil_texture
//' @keywords internal
// [[Rcpp::export("soil_conductivity")]]
NumericVector conductivity(DataFrame soil, String model="SX", bool mmol = true) {
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
      K[l] = unsaturatedConductivitySaxton_c(Theta[l], clay[l], sand[l], bd[l], om[l], mmol);
    }
  } else {
    NumericVector psiSoil = psi(soil, model);
    NumericVector Ksat = soil["Ksat"];
    for(int l=0;l<nlayers;l++) {
      NumericVector n =soil["VG_n"];
      NumericVector alpha = soil["VG_alpha"];
      NumericVector theta_res = soil["VG_theta_res"];
      NumericVector theta_sat = soil["VG_theta_sat"];
      K[l] = psi2kVanGenuchten_c(Ksat[l], n[l], alpha[l], theta_res[l], theta_sat[l], psiSoil[l]);
      // mmolH20·m-1·s-1·MPa-1 to cm/day
      if(!mmol) K[l] = K[l]/cmdTOmmolm2sMPa;
    }
  } 
  return(K);
}

//' @rdname soil_texture
//' @keywords internal
// [[Rcpp::export("soil_capacitance")]]
 NumericVector capacitance(DataFrame soil, String model="SX") {
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
       C[l] = psi2cVanGenuchten_c(n[l], alpha[l], theta_res[l], theta_sat[l], psiSoil[l]);
     }
   } 
   return(C);
 }

//' @rdname soil_texture
//' @keywords internal
// [[Rcpp::export("soil_saturatedWaterDepth")]]
double saturatedWaterDepth(DataFrame soil, String model = "SX") {
  NumericVector widths = soil["widths"];
  NumericVector W = soil["W"];
  NumericVector Theta_FC = thetaFC(soil, model);
  NumericVector Theta_SAT = thetaSAT(soil, model);
  int nlayers = W.length();
  double z = 0.0;
  int nunsaturated = 0;
  for(int l=0;l<nlayers;l++) {
    if(W[l]>1.0) {
      z = z + widths[l]*(Theta_SAT[l]-Theta_FC[l]*W[l])/(Theta_SAT[l]-Theta_FC[l]);
    } else {
      z = z + widths[l];
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
//' @keywords internal
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
//' @keywords internal
// [[Rcpp::export("soil_campbellParamsClappHornberger")]]
NumericVector campbellParamsClappHornberger(String soilType) {
   NumericVector cp(4,NA_REAL);
   ClappHornberger clapp = ClappHornberger(soilType.get_cstring());
   cp[0] = clapp.theta_sat;
   cp[1] = clapp.psi_sat_cm;
   cp[2] = clapp.b;
   cp[3] = clapp.K_sat_cm_h;
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
//' @keywords internal
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
//' @param x A data frame of soil parameters (see an example in \code{\link{defaultSoilParams}}).
//' @param VG_PTF Pedotransfer functions to obtain parameters for the van Genuchten-Mualem equations. Either \code{"Carsel"} (Carsel and Parrish 1988) or \code{"Toth"} (Toth et al. 2015).
//' 
//' @return
//' Function \code{soil} returns a data frame of class \code{soil} with the following columns:
//' \itemize{
//'   \item{\code{widths}: Width of soil layers (in mm).}
//'   \item{\code{sand}: Sand percentage for each layer (in percent volume).}
//'   \item{\code{clay}: Clay percentage for each layer (in percent volume).}
//'   \item{\code{om}: Organic matter percentage for each layer (in percent volume).}
//'   \item{\code{nitrogen}: Sum of total nitrogen (ammonia, organic and reduced nitrogen) for each layer (in g/kg).}
//'   \item{\code{ph}: pH in water of each layer (0-14).}
//'   \item{\code{rfc}: Percentage of rock fragment content for each layer.}
//'   \item{\code{macro}: Macroporosity for each layer (estimated using Stolf et al. 2011).}
//'   \item{\code{Ksat}: Saturated soil conductivity for each layer (in mmol·m-1·s-1·MPa-1, estimated using function \code{\link{soil_saturatedConductivitySX}}.}
//'   \item{\code{VG_alpha}, \code{VG_n}, \code{VG_theta_res}, \code{VG_theta_sat}: Parameters for van Genuchten's pedotransfer functions, for each layer, corresponding to the USDA texture type.}
//'   \item{\code{W}: State variable with relative water content of each layer (in as proportion relative to FC).}
//'   \item{\code{Temp}: State variable with temperature (in Celsius) of each layer.}
//' }
//' 
//' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
//' 
//' @details 
//' Function \code{summary} prompts a description of soil characteristics and state variables (water content and temperature) 
//' according to a water retention curve (either Saxton's or Van Genuchten's). 
//' Volume at field capacity is calculated assuming a soil water potential equal to fieldCapacityPsi MPa. 
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
//' @seealso \code{\link{plot.soil}}, \code{\link{soil_redefineLayers}}, \code{\link{soil_psi2thetaSX}}, \code{\link{soil_psi2thetaVG}}, \code{\link{spwb}}, \code{\link{defaultSoilParams}}
//' 
//' @examples
//' # Default parameters
//' df_soil <- defaultSoilParams()
//' 
//' # Initializes soil
//' s <- soil(df_soil)
//' s
//' 
//' # Prints soil characteristics according to Saxton's water retention curve
//' summary(s, model="SX")
//' 
//' # Prints soil characteristics according to Van Genuchten's water retention curve
//' summary(s, model="VG")
//' 
//' # Add columns 'VG_theta_sat' and 'VG_theta_res' with custom values
//' df_soil$VG_theta_sat <- 0.400 
//' df_soil$VG_theta_res <- 0.040 
//' 
//' # Reinitialize soil (should override estimations)
//' s2 <- soil(df_soil)
//' s2
//' summary(s2, model="VG")
//' @name soil
// [[Rcpp::export("soil")]]
DataFrame soilInit(DataFrame x, String VG_PTF = "Toth") {
  int nlayers = x.nrow();
  if((VG_PTF!="Toth") && (VG_PTF!="Carsel")) stop("Wrong VG_PTF (should be either 'Toth' or 'Carsel')");
  //Soil parameters related to physical structure
  NumericVector widths = clone(as<NumericVector>(x["widths"]));
  NumericVector clay = clone(as<NumericVector>(x["clay"]));
  NumericVector sand = clone(as<NumericVector>(x["sand"]));
  NumericVector bd = clone(as<NumericVector>(x["bd"]));
  NumericVector rfc = clone(as<NumericVector>(x["rfc"]));
  
  if(any(is_na(widths))) stop("Missing values in soil 'widths'");
  if(any(is_na(clay))) stop("Missing values in soil 'clay'");
  if(any(is_na(sand))) stop("Missing values in soil 'sand'");
  if(any(is_na(bd))) stop("Missing values in soil 'bd'");
  if(any(is_na(rfc))) stop("Missing values in soil 'rfc'");

  NumericVector W(nlayers, 1.0);
  if(x.containsElementNamed("W")) W = clone(as<NumericVector>(x["W"]));
  if(any(is_na(W))) stop("Missing values in soil 'W'");
  
  //Optional
  NumericVector om(nlayers, NA_REAL);
  NumericVector nitrogen(nlayers, NA_REAL);
  NumericVector ph(nlayers, NA_REAL);
  if(x.containsElementNamed("om")) om = clone(as<NumericVector>(x["om"]));
  if(x.containsElementNamed("nitrogen")) nitrogen = clone(as<NumericVector>(x["nitrogen"]));
  if(x.containsElementNamed("ph")) ph = clone(as<NumericVector>(x["ph"]));
  
  //Parameters to be calculated and state variables
  NumericVector macro(nlayers, NA_REAL);
  NumericVector temperature(nlayers, NA_REAL);
  CharacterVector usda_Type(nlayers, NA_STRING);
  NumericVector VG_alpha(nlayers, NA_REAL);
  NumericVector VG_n(nlayers, NA_REAL);
  NumericVector VG_theta_res(nlayers, NA_REAL);
  NumericVector VG_theta_sat(nlayers, NA_REAL);
  NumericVector Ksat(nlayers, NA_REAL);
  
  //Get parameters from input if specified
  if(x.containsElementNamed("VG_alpha")) VG_alpha = clone(as<NumericVector>(x["VG_alpha"]));
  if(x.containsElementNamed("VG_n")) VG_n = clone(as<NumericVector>(x["VG_n"]));
  if(x.containsElementNamed("VG_theta_res")) VG_theta_res = clone(as<NumericVector>(x["VG_theta_res"]));
  if(x.containsElementNamed("VG_theta_sat")) VG_theta_sat = clone(as<NumericVector>(x["VG_theta_sat"]));
  if(x.containsElementNamed("Ksat")) Ksat = clone(as<NumericVector>(x["Ksat"]));
  if(x.containsElementNamed("macro")) macro = clone(as<NumericVector>(x["macro"]));
  if(x.containsElementNamed("usda")) usda_Type = clone(as<CharacterVector>(x["usda"]));
  
  for(int l=0;l<nlayers;l++) {
    if(CharacterVector::is_na(usda_Type[l]))  usda_Type[l] = USDAType(clay[l],sand[l]);
    // Stolf, R., Thurler, A., Oliveira, O., Bacchi, S., Reichardt, K., 2011. Method to estimate soil macroporosity and microporosity based on sand content and bulk density. Rev. Bras. Ciencias do Solo 35, 447–459.
    if(NumericVector::is_na(macro[l])) macro[l] = std::max(0.0,0.693 - 0.465*bd[l] + 0.212*(sand[l]/100.0));
    
    NumericVector vgl;
    if(VG_PTF=="Carsel") {
      vgl = vanGenuchtenParamsCarsel(usda_Type[l]); 
      if(NumericVector::is_na(Ksat[l])) {
        Ksat[l] = vgl[4]; //Use Carsel estimate for Ksat
      }
    } else if(VG_PTF=="Toth") {
      // if(l==0) vgl = vanGenuchtenParamsToth(clay[l], sand[l], om[l], bd[l], TRUE);
      //Use non-top soil equation for all layers
      vgl = vanGenuchtenParamsToth(clay[l], sand[l], om[l], bd[l], false);
      if(NumericVector::is_na(Ksat[l])) {
        Ksat[l] = saturatedConductivitySaxton_c(clay[l], sand[l], bd[l], om[l], true); 
      }
    } else {
      stop("Wrong value for 'VG_PTF'");
    }
    if(NumericVector::is_na(VG_alpha[l])) VG_alpha[l] = vgl[0];
    if(NumericVector::is_na(VG_n[l])) VG_n[l] = vgl[1];
    if(NumericVector::is_na(VG_theta_res[l])) VG_theta_res[l] = vgl[2];
    if(NumericVector::is_na(VG_theta_sat[l])) VG_theta_sat[l] = vgl[3];
  }
  DataFrame l = DataFrame::create(_["widths"] = widths,
                        _["sand"] = sand, 
                        _["clay"] = clay, 
                        _["usda"] = usda_Type,
                        _["om"] = om, 
                        _["nitrogen"] = nitrogen,
                        _["ph"] = ph,
                        _["bd"] = bd,
                        _["rfc"] = rfc,
                        _["macro"] = macro,
                        _["Ksat"] = Ksat,
                        _["VG_alpha"] = VG_alpha,
                        _["VG_n"] = VG_n, 
                        _["VG_theta_res"] = VG_theta_res,
                        _["VG_theta_sat"] = VG_theta_sat,
                        _["W"] = W, 
                        _["Temp"] = temperature);
  l.attr("class") = CharacterVector::create("soil","data.frame");
  return(l);
}


// [[Rcpp::export(".modifySoilLayerParam")]]
void modifySoilLayerParam(DataFrame soil, String paramName, int layer, double newValue, 
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
  NumericVector widths = as<NumericVector>(soil["widths"]);
  NumericVector macro = as<NumericVector>(soil["macro"]);
  NumericVector VG_alpha = as<NumericVector>(soil["VG_alpha"]);
  NumericVector VG_n = as<NumericVector>(soil["VG_n"]);
  NumericVector VG_theta_res = as<NumericVector>(soil["VG_theta_res"]);
  NumericVector VG_theta_sat = as<NumericVector>(soil["VG_theta_sat"]);
  NumericVector Ksat = as<NumericVector>(soil["Ksat"]);
  
  int nlayers = widths.size();

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
    Ksat[l] = saturatedConductivitySaxton_c(clay[l], sand[l], bd[l], om[l], true);
  }
}

// [[Rcpp::export(.testSoilDataFrameToStructure)]]
NumericVector testSoilDataFrameToStructure(DataFrame x, String model = "VG") {
  Soil soil(x, model);
  NumericVector sizes = {sizeof(x),sizeof(soil)};
  return(sizes);
}
