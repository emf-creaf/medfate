// [[Rcpp::interfaces(r,cpp)]]

#include <RcppArmadillo.h>
#include "soil.h"
#include "numerical_solving_c.h"
#include "soil_thermodynamics_c.h"
using namespace Rcpp;


/**
 * Soil thermal conductivity 
 *
 * Dharssi, I., Vidale, P.L., Verhoef, A., MacPherson, B., Jones, C., & Best, M. 2009. New soil physical properties implemented in the Unified Model at PS18. 9–12.
 * Best et al. 2011
 */
void layerThermalConductivity_c(std::vector<double>& thermalCond,
                                const std::vector<double>& sand, const std::vector<double>& clay, 
                                const std::vector<double>& W, const std::vector<double>& Theta_SAT, const std::vector<double>& Theta_FC,
                                const std::vector<double>& Temp) {
  int nlayers = sand.size();
  for(int l=0;l<nlayers;l++) {
    double silt = 100 - sand[l] - clay[l];
    double lambda_m = ((thermalConductivitySand*sand[l])+(thermalConductivitySilt*silt)+(thermalConductivityClay*clay[l]))/(silt+sand[l]+clay[l]);
    double lambda_dry = pow(thermalConductivityAir, Theta_FC[l])*pow(lambda_m, (1.0-Theta_FC[l]));
    double Ke = 0.0;
    if(W[l]>=0.1) Ke = log10(W[l]) + 1.0;
    double lambda_s = std::max(1.58,std::min(2.2,1.58 + 12.4*(lambda_dry-0.25)));
    thermalCond[l] = (lambda_s-lambda_dry)*Ke + lambda_dry;
  }
}


/**
 * Soil thermal capacity. Simplified from:
 * 
 *  returns - J·m-3·K-1
 * Cox, P.M., Betts, R.A., Bunton, C.B., Essery, R.L.H., Rowntree, P.R., & Smith, J. 1999. The impact of new land surface physics on the GCM simulation of climate and climate sensitivity. Climate Dynamics 15: 183–203.
 */
void layerThermalCapacity_c(std::vector<double>& thermalCap,
                            const std::vector<double>& sand, const std::vector<double>& clay, 
                            const std::vector<double>& W, const std::vector<double>& Theta_SAT, const std::vector<double>& Theta_FC,
                            const std::vector<double>& Temp) {
  int nlayers = sand.size();
  for(int l=0;l<nlayers;l++) {
    thermalCap[l] = ((sand[l]*capacitySand)+(clay[l]*capacityClay) + ((100.0-clay[l]-sand[l])*capacitySilt))/100.0;
    thermalCap[l] = thermalCap[l] + 4.19*1e3*1000.0*Theta_FC[l]*W[l];//Add water
  }
}


/**
 * Calculates midpoints of soil layers
 */
void midpoints_c(std::vector<double> midZ, const std::vector<double>& widths) {
  int nlayers = widths.size();
  double sumz = 0.0;
  for(int l = 0;l<nlayers; l++) {
    midZ[l] = sumz + widths[l]/2.0;
    sumz = sumz + widths[l];
  }
}


/**
 * Soil temperature gradient (in ºC/m)
 */
//' @name soil_thermodynamics
//' @keywords internal
// [[Rcpp::export("soil_temperatureGradient")]]
void temperatureGradient_c(std::vector<double>& gradTemp, const std::vector<double>& widths, const std::vector<double>& Temp) {
  int nlayers = Temp.size();
  std::vector<double> midZ(nlayers);
  midpoints_c(midZ, widths);
  if(nlayers>1) {
    for(int l = 0;l<nlayers-1; l++) {
      gradTemp[l] = (Temp[l+1]-Temp[l])/(0.001*(midZ[l+1]-midZ[l]));
    }
  }
  gradTemp[nlayers-1] = (15.5-Temp[nlayers-1])/(0.001*(10000.0-midZ[nlayers-1])); //15.5º at 10 m
}

void temperatureChange_inner_c(SoilEnergyBalance_COMM& SEBcomm, 
                               const std::vector<double>& widths, const std::vector<double>& Temp,
                               const std::vector<double>& sand, const std::vector<double>& clay,
                               const std::vector<double>& W, const std::vector<double>& Theta_SAT, const std::vector<double>& Theta_FC,
                               double Gdown, double tstep) {
  int nlayers = Temp.size();
  
  std::vector<double> k(nlayers);
  std::vector<double> Cv(nlayers);
  layerThermalConductivity_c(k , sand, clay, 
                             W, Theta_SAT, Theta_FC, 
                             Temp);
  layerThermalCapacity_c(Cv, sand, clay, 
                         W, Theta_SAT, Theta_FC, 
                         Temp);
  
  std::vector<double>& a = SEBcomm.a;
  std::vector<double>& b = SEBcomm.b;
  std::vector<double>& c = SEBcomm.c;
  std::vector<double>& d = SEBcomm.d;
  std::vector<double>& e = SEBcomm.e;
  std::vector<double>& f = SEBcomm.f;
  
  std::vector<double>& dZ_m = SEBcomm.dZ_m;
  std::vector<double>& dZUp = SEBcomm.dZUp;
  std::vector<double>& dZDown = SEBcomm.dZDown;

  std::vector<double>& Zcent = SEBcomm.Zcent;
  std::vector<double>& Zup = SEBcomm.Zup;
  std::vector<double>& Zdown = SEBcomm.Zdown;
  
  std::vector<double>& k_up = SEBcomm.k_up;
  std::vector<double>& k_down = SEBcomm.k_down;
  std::vector<double>& tempch = SEBcomm.tempch;

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
  tridiagonalSolving_c(a,b,c,d, e,f, tempch);
}


