// [[Rcpp::interfaces(r,cpp)]]

#include <RcppArmadillo.h>
#include "soil.h"
#include "numerical_solving.h"
#include "communication_structures.h"
#include "soil_thermodynamics.h"
#include "soil_thermodynamics_c.h"
using namespace Rcpp;


/**
 * Soil thermal conductivity 
 *
 * Dharssi, I., Vidale, P.L., Verhoef, A., MacPherson, B., Jones, C., & Best, M. 2009. New soil physical properties implemented in the Unified Model at PS18. 9–12.
 * Best et al. 2011
 */
NumericVector layerThermalConductivity(NumericVector sand, NumericVector clay, 
                                       NumericVector W, NumericVector Theta_SAT, NumericVector Theta_FC,
                                       NumericVector Temp) {
  int nlayers = sand.length();
  std::vector<double> thermalCond_vec(nlayers);
  layerThermalConductivity_c(thermalCond_vec,
                             as<std::vector<double>>(sand), as<std::vector<double>>(clay),
                             as<std::vector<double>>(W), as<std::vector<double>>(Theta_SAT), as<std::vector<double>>(Theta_FC),
                             as<std::vector<double>>(Temp));
  NumericVector thermalCond = Rcpp::wrap(thermalCond_vec);
  return(thermalCond);
}


/**
 * Soil thermal capacity. Simplified from:
 * 
 *  returns - J·m-3·K-1
 * Cox, P.M., Betts, R.A., Bunton, C.B., Essery, R.L.H., Rowntree, P.R., & Smith, J. 1999. The impact of new land surface physics on the GCM simulation of climate and climate sensitivity. Climate Dynamics 15: 183–203.
 */
NumericVector layerThermalCapacity(NumericVector sand, NumericVector clay, 
                                   NumericVector W, NumericVector Theta_SAT, NumericVector Theta_FC,
                                   NumericVector Temp) {
  int nlayers = sand.length();
  std::vector<double> thermalCap_vec(nlayers);
  layerThermalCapacity_c(thermalCap_vec,
                             as<std::vector<double>>(sand), as<std::vector<double>>(clay),
                             as<std::vector<double>>(W), as<std::vector<double>>(Theta_SAT), as<std::vector<double>>(Theta_FC),
                             as<std::vector<double>>(Temp));
  NumericVector thermalCap = Rcpp::wrap(thermalCap_vec);
  return(thermalCap);
}


//' Soil thermodynamic functions
//' 
//' Functions \code{soil_thermalConductivity} and \code{soil_thermalCapacity} calculate thermal conductivity and thermal capacity 
//' for each soil layer, given its texture and water content. Functions \code{soil_temperatureGradient} and \code{soil_temperatureChange} 
//' are used to calculate soil temperature gradients (in ºC/m) and temporal temperature change (in ºC/s) 
//' given soil layer texture and water content (and possibly including heat flux from above).
//' 
//' @param soil Soil object (returned by function \code{\link{soil}}).
//' @param model Either 'SX' or 'VG' for Saxton's or Van Genuchten's pedotransfer models.
//' @param widths Width of soil layers (in mm).
//' @param Temp Temperature (in ºC) for each soil layer.
//' @param clay Percentage of clay (in percent weight) for each layer.
//' @param sand Percentage of sand (in percent weight) for each layer.
//' @param W Soil moisture (in percent of field capacity) for each layer.
//' @param Theta_SAT Relative water content (in percent volume) at saturation for each layer.
//' @param Theta_FC Relative water content (in percent volume) at field capacity for each layer.
//' @param Gdown Downward heat flux from canopy to soil (in W·m-2).
//' @param tstep Time step (interval) in seconds.
//' 
//' @return 
//' Function \code{soil_thermalConductivity} returns a vector with values of thermal conductivity (W/m/ºK) for each soil layer. 
//' 
//' Function \code{soil_thermalCapacity} returns a vector with values of heat storage capacity (J/m3/ºK) for each soil layer. 
//' 
//' Function \code{soil_temperatureGradient} returns a vector with values of temperature gradient between consecutive soil layers. 
//' 
//' Function \code{soil_temperatureChange} returns a vector with values of instantaneous temperature change (ºC/s) for each soil layer.
//' 
//' @references
//' Cox, P.M., Betts, R.A., Bunton, C.B., Essery, R.L.H., Rowntree, P.R., and Smith, J. 1999. The impact of new land surface physics on the GCM simulation of climate and climate sensitivity. Climate Dynamics 15: 183–203.
//' 
//' Dharssi, I., Vidale, P.L., Verhoef, A., MacPherson, B., Jones, C., and Best, M. 2009. New soil physical properties implemented in the Unified Model at PS18. 9–12.
//' 
//' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
//' 
//' @seealso \code{\link{soil}}
//' 
//' @examples
//' #Define soil and complete parameters
//' examplesoil = soil(defaultSoilParams(4))
//' 
//' soil_thermalConductivity(examplesoil)
//' soil_thermalCapacity(examplesoil)
//' 
//' #Values change when altering water content (drier layers have lower conductivity and capacity)
//' examplesoil$W = c(0.1, 0.4, 0.7, 1.0)
//' soil_thermalConductivity(examplesoil)
//' soil_thermalCapacity(examplesoil)
//' 
//' @name soil_thermodynamics
//' @keywords internal
// [[Rcpp::export("soil_thermalCapacity")]]
NumericVector thermalCapacity(DataFrame soil, String model = "SX") {
  NumericVector sand = soil["sand"];
  NumericVector clay = soil["clay"];
  NumericVector W = soil["W"];
  NumericVector Temp = soil["Temp"];
  NumericVector Theta_FC = thetaFC(soil, model);
  NumericVector Theta_SAT = thetaSAT(soil, model);
  return(layerThermalCapacity(sand, clay, 
                              W, Theta_SAT, Theta_FC,
                              Temp));
}

//' @rdname soil_thermodynamics
//' @keywords internal
// [[Rcpp::export("soil_thermalConductivity")]]
NumericVector thermalConductivity(DataFrame soil, String model = "SX") {
  NumericVector sand = soil["sand"];
  NumericVector clay = soil["clay"];
  NumericVector W = soil["W"];
  NumericVector Theta_FC = thetaFC(soil, model);
  NumericVector Theta_SAT = thetaSAT(soil, model);
  NumericVector Temp = soil["Temp"];
  return(layerThermalConductivity(sand, clay, 
                                  W, Theta_SAT, Theta_FC,
                                  Temp));
}


/**
 * Soil temperature gradient (in ºC/m)
 */
//' @name soil_thermodynamics
//' @keywords internal
// [[Rcpp::export("soil_temperatureGradient")]]
NumericVector temperatureGradient(NumericVector widths, NumericVector Temp) {
  int nlayers = Temp.length();
  std::vector<double> gradTemp_vec(nlayers);
  temperatureGradient_c(gradTemp_vec,
                        as<std::vector<double>>(widths), 
                        as<std::vector<double>>(Temp));
  NumericVector gradTemp = Rcpp::wrap(gradTemp_vec);
  return(gradTemp);
}

NumericVector temperatureChange_inner(List SEBcommunication, NumericVector widths, NumericVector Temp,
                                      NumericVector sand, NumericVector clay,
                                      NumericVector W, NumericVector Theta_SAT, NumericVector Theta_FC,
                                      double Gdown, double tstep) {
  NumericVector k = layerThermalConductivity(sand, clay, 
                                             W, Theta_SAT, Theta_FC, 
                                             Temp);
  NumericVector Cv = layerThermalCapacity(sand, clay, 
                                          W, Theta_SAT, Theta_FC, 
                                          Temp);
  int nlayers = Temp.length();
  
  NumericVector a = SEBcommunication[SOILEBCOM_a];
  NumericVector b = SEBcommunication[SOILEBCOM_b];
  NumericVector c = SEBcommunication[SOILEBCOM_c];
  NumericVector d = SEBcommunication[SOILEBCOM_d];
  NumericVector e = SEBcommunication[SOILEBCOM_e];
  NumericVector f = SEBcommunication[SOILEBCOM_f];
  NumericVector dZ_m = SEBcommunication[SOILEBCOM_dZ_m];
  NumericVector dZUp = SEBcommunication[SOILEBCOM_dZUp];
  NumericVector dZDown = SEBcommunication[SOILEBCOM_dZDown];
  NumericVector Zcent = SEBcommunication[SOILEBCOM_Zcent];
  NumericVector Zup = SEBcommunication[SOILEBCOM_Zup];
  NumericVector Zdown = SEBcommunication[SOILEBCOM_Zdown];
  NumericVector k_up = SEBcommunication[SOILEBCOM_k_up];
  NumericVector k_down = SEBcommunication[SOILEBCOM_k_down];
  NumericVector tempch = SEBcommunication[SOILEBCOM_tempch];
  
  //Estimate layer interfaces
  for(int l=0;l<nlayers;l++) dZ_m[l] = widths[l]*0.001; //mm to m
  
  for(int l=0;l<nlayers;l++) {
    if(l==0) { //first layer
      dZUp[l] = dZ_m[0]/2.0; //Distance from ground to mid-layer
      Zcent[l] = -1.0*dZ_m[0]/2.0; //Center of the layer (negative downwards)
      Zup[l] = 0.0; // Upper limit
      Zdown[l] = -1.0*dZ_m[0]; //Lower limit (negative downwards)
    } else {
      dZUp[l] = (dZ_m[l - 1]/2.0) + (dZ_m[l]/2.0); // Distance between mid-layers
      Zcent[l] = Zcent[l-1] - (dZ_m[l]/2.0); //Center of the layer (negative downwards)
      Zup[l] = Zdown[l-1]; //Upper limit
      Zdown[l] = Zdown[l-1] - dZ_m[l]; //Lower limit (negative downwards)
    }
    if(l<(nlayers - 1)) {
      dZDown[l] = (dZ_m[l]/2.0) + (dZ_m[l + 1]/2.0);
    } else { //last layer
      dZDown[l] = dZ_m[l]/2.0;
    }
  }
  for(int l=0;l<nlayers;l++) {
    k_up[l] = 0.0;
    k_down[l] = 0.0;
    if(l==0) { //first layer
      k_down[l] = (k[l]*k[l+1]*(Zcent[l]-Zcent[l+1]))/(k[l]*(Zdown[l] - Zcent[l+1]) + k[l+1]*(Zcent[l] - Zdown[l]));
      a[l] = 0.0;
      c[l] = -1.0*(k_down[l]/dZDown[l]);
      b[l] = (Cv[l]*dZ_m[l]/tstep) - c[l];
      d[l] = Gdown - (k_down[l]/dZDown[l])*(Temp[l] - Temp[l+1]);
    } else if(l<(nlayers - 1)) {
      k_up[l] = k_down[l-1];
      k_down[l] = (k[l]*k[l+1]*(Zcent[l]-Zcent[l+1]))/(k[l]*(Zdown[l] - Zcent[l+1]) + k[l+1]*(Zcent[l] - Zdown[l]));
      a[l] = -1.0*(k_up[l]/dZUp[l]);
      c[l] = -1.0*(k_down[l]/dZDown[l]);
      b[l] = (Cv[l]*dZ_m[l]/tstep) - a[l] - c[l];
      d[l] = (k_up[l]/dZUp[l])*(Temp[l-1] - Temp[l]) - (k_down[l]/dZDown[l])*(Temp[l] - Temp[l+1]);
    } else { // last layer
      a[l] = -1.0*(k_up[l]/dZUp[l]);
      c[l] = 0.0;
      b[l] = (Cv[l]*dZ_m[l]/tstep) - a[l];
      d[l] = (k_up[l]/dZUp[l])*(Temp[l-1] - Temp[l]);
    }
  }
  tridiagonalSolving(a,b,c,d, e,f, tempch);
  return(tempch);
}
//' @name soil_thermodynamics
//' @keywords internal
// [[Rcpp::export("soil_temperatureChange")]]
NumericVector temperatureChange(NumericVector widths, NumericVector Temp,
                                 NumericVector sand, NumericVector clay,
                                 NumericVector W, NumericVector Theta_SAT, NumericVector Theta_FC,
                                 double Gdown, double tstep) {
   List SEBcommunication = communicationSoilEnergyBalance(widths.size());
   return(temperatureChange_inner(SEBcommunication, widths, Temp,
                                  sand, clay,
                                  W, Theta_SAT, Theta_FC,
                                  Gdown, tstep));
}


