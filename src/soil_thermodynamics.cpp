// [[Rcpp::interfaces(r,cpp)]]

#include <Rcpp.h>
#include "soil.h"
#include "numerical_solving.h"
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
 * Calculates midpoints of soil layers
 */
NumericVector midpoints(NumericVector widths) {
  int nlayers = widths.length();
  double sumz = 0.0;
  NumericVector midZ(nlayers);
  for(int l = 0;l<nlayers; l++) {
    midZ[l] = sumz + widths[l]/2.0;
    sumz = sumz + widths[l];
  }
  return(midZ);
}

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
NumericVector layerThermalCapacity(NumericVector sand, NumericVector clay, 
                                   NumericVector W, NumericVector Theta_SAT, NumericVector Theta_FC,
                                   NumericVector Temp) {
  int nlayers = sand.length();
  NumericVector thermalCap(nlayers,0.0);
  for(int l=0;l<nlayers;l++) {
    thermalCap[l] = ((sand[l]*capacitySand)+(clay[l]*capacityClay) + ((100.0-clay[l]-sand[l])*capacitySilt))/100.0;
    thermalCap[l] = thermalCap[l] + 4.19*1e3*1000.0*Theta_FC[l]*W[l];//Add water
  }
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
  NumericVector midZ = midpoints(widths);
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

//' @name soil_thermodynamics
//' @keywords internal
// [[Rcpp::export("soil_temperatureChange")]]
NumericVector temperatureChange(NumericVector widths, NumericVector Temp,
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

  NumericVector dZ_m = widths*0.001; //mm to m

  //Estimate layer interfaces
  NumericVector dZUp(nlayers), dZDown(nlayers), Zcent(nlayers), Zup(nlayers), Zdown(nlayers);
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
  NumericVector k_up(nlayers), k_down(nlayers);
  NumericVector a(nlayers),  b(nlayers),  c(nlayers),  d(nlayers);
  for(int l=0;l<nlayers;l++) {
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
  NumericVector tempch = tridiagonalSolving(a,b,c,d);
  return(tempch);
}
// NumericVector temperatureChange(NumericVector widths, NumericVector Temp,
//                                 NumericVector sand, NumericVector clay,
//                                 NumericVector W, NumericVector Theta_FC,
//                                 double Gdown, double tstep) {
//   NumericVector lambda = layerThermalConductivity(sand, clay, W, Theta_FC);
//   NumericVector Ca = layerThermalCapacity(sand, clay, W, Theta_FC);
//   int nlayers = Temp.length();
//   NumericVector gradTemp = temperatureGradient(widths, Temp);
//   NumericVector midZ = midpoints(widths);
//   double Gup = -Gdown; //Gdown > 0 when net flux is in the direction of soil
//   double Gi;
//   NumericVector tempch(nlayers);
//   for(int l = 0;l<nlayers; l++) {
//     Gi = lambda[l]*gradTemp[l]; //Gi < 0 when net flux is downward
//     tempch[l] = (Gi-Gup)/(Ca[l]*0.001*widths[l])*tstep;
//     Gup = Gi;
//   }
//   return(tempch);
// }

