#include <RcppArmadillo.h>
#include "medfate.h"

#ifndef SOIL_THERMODYNAMICS_C_H
#define SOIL_THERMODYNAMICS_C_H


// sand 1.7-2.9 W·m-1·K-1, clay 0.8-6.3 W·m-1·K-1 [Geiger et al. The Climate near the Ground]
const double thermalConductivitySand = 1.57025; //W·m-1·K-1 From Dharssi et al. 2009
const double thermalConductivitySilt = 1.57025; //W·m-1·K-1 From Dharssi et al. 2009
const double thermalConductivityClay = 1.16025; //W·m-1·K-1 From Dharssi et al. 2009
const double thermalConductivityAir = 0.025; //W·m-1·K-1 From Dharssi et al. 2009
const double capacitySand = 1.25*1e6; //kg·m-3 
const double capacitySilt = 1.19*1e6; //kg·m-3 
const double capacityClay = 1.23*1e6; //kg·m-3 

// ----------------------------------------------------------------------------
// Soil Energy Balance Communication Structure (15 fields)
// ----------------------------------------------------------------------------

struct SoilEnergyBalance_COMM {
  
  // Geometry
  std::vector<double> dZ_m;
  std::vector<double> dZUp;
  std::vector<double> dZDown;
  std::vector<double> Zup;
  std::vector<double> Zdown;
  std::vector<double> Zcent;
  
  // Thermal conductivity polynomial coefficients
  std::vector<double> a;
  std::vector<double> b;
  std::vector<double> c;
  std::vector<double> d;
  std::vector<double> e;
  std::vector<double> f;
  
  // Thermal conductivity
  std::vector<double> k_up;
  std::vector<double> k_down;
  
  // Temperature change
  std::vector<double> tempch;
  
  // Constructor
  SoilEnergyBalance_COMM(size_t nlayers = 0) {
    dZ_m = std::vector<double>(nlayers, medfate::NA_DOUBLE);
    dZUp = std::vector<double>(nlayers, medfate::NA_DOUBLE);
    dZDown = std::vector<double>(nlayers, medfate::NA_DOUBLE);
    Zup = std::vector<double>(nlayers, medfate::NA_DOUBLE);
    Zdown = std::vector<double>(nlayers, medfate::NA_DOUBLE);
    Zcent = std::vector<double>(nlayers, medfate::NA_DOUBLE);
    
    k_up = std::vector<double>(nlayers, medfate::NA_DOUBLE);
    k_down = std::vector<double>(nlayers, medfate::NA_DOUBLE);
    tempch = std::vector<double>(nlayers, medfate::NA_DOUBLE);
    
    a = std::vector<double>(nlayers, medfate::NA_DOUBLE);
    b = std::vector<double>(nlayers, medfate::NA_DOUBLE);
    c = std::vector<double>(nlayers, medfate::NA_DOUBLE);
    d = std::vector<double>(nlayers, medfate::NA_DOUBLE);
    e = std::vector<double>(nlayers, medfate::NA_DOUBLE);
    f = std::vector<double>(nlayers, medfate::NA_DOUBLE);
  }
};


void layerThermalConductivity_c(std::vector<double>& thermalCond,
                                const std::vector<double>& sand, const std::vector<double>& clay, 
                                const std::vector<double>& W, const std::vector<double>& Theta_SAT, const std::vector<double>& Theta_FC,
                                const std::vector<double>& Temp);

void layerThermalCapacity_c(std::vector<double>& thermalCap,
                            const std::vector<double>& sand, const std::vector<double>& clay, 
                            const std::vector<double>& W, const std::vector<double>& Theta_SAT, const std::vector<double>& Theta_FC,
                            const std::vector<double>& Temp);

void temperatureGradient_c(std::vector<double>& gradTemp, const std::vector<double>& widths, const std::vector<double>& Temp);

void temperatureChange_inner_c(SoilEnergyBalance_COMM& SEBcomm, 
                               const std::vector<double>& widths, const std::vector<double>& Temp,
                               const std::vector<double>& sand, const std::vector<double>& clay,
                               const std::vector<double>& W, const std::vector<double>& Theta_SAT, const std::vector<double>& Theta_FC,
                               double Gdown, double tstep);
#endif
