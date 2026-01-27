#include "medfate.h"
#include <RcppArmadillo.h>
#include <vector>
#include <string>
#include <limits>
#include <cmath>

#ifndef COMMUNICATION_STRUCTURES_C_H
#define COMMUNICATION_STRUCTURES_C_H


// ----------------------------------------------------------------------------
// Soil Water Balance Communication Structure (50 fields)
// ----------------------------------------------------------------------------

struct SoilWaterBalance_COMM {
  size_t nlayers;
  
  // Layer dimensions
  std::vector<double> dZ_m;
  std::vector<double> dZUp;
  std::vector<double> dZDown;
  
  // Soil hydraulic properties
  std::vector<double> lambda;
  
  // Water content states
  std::vector<double> theta_micro;
  std::vector<double> theta_b;
  std::vector<double> theta_macro;
  std::vector<double> theta_sat_fict;
  
  // Saturated conductivity
  std::vector<double> Ksat_b;
  std::vector<double> Ksat_b_ms;
  std::vector<double> Ksat;
  std::vector<double> Ksat_ms;
  
  // Water potential
  std::vector<double> Psi;
  std::vector<double> Psi_m;
  std::vector<double> Psi_step;
  std::vector<double> Psi_step_m;
  
  // Conductivity
  std::vector<double> K;
  std::vector<double> K_ms;
  std::vector<double> Kbc;
  std::vector<double> Kbc_ms;
  std::vector<double> K_step;
  std::vector<double> K_step_ms;
  
  // Capacitance
  std::vector<double> C;
  std::vector<double> C_m;
  std::vector<double> C_step;
  std::vector<double> C_step_m05;
  
  // Macro-porosity
  std::vector<double> S_macro;
  std::vector<double> e_macro;
  std::vector<double> S_macro_step;
  std::vector<double> theta_macro_step;
  
  // Polytelnyic coefficients
  std::vector<double> a;
  std::vector<double> b;
  std::vector<double> c;
  std::vector<double> d;
  std::vector<double> e;
  std::vector<double> f;
  
  // Additional properties
  std::vector<double> Kmacro_ms;
  std::vector<double> K_step_ms05;
  std::vector<double> waterFluidity;
  std::vector<double> finalSourceSinks_m3s;
  std::vector<double> capill_below;
  std::vector<double> drain_above;
  std::vector<double> drain_below;
  std::vector<double> theta_micro_step;
  std::vector<double> lateral_flows_step_mm;
  
  // Constructor
  SoilWaterBalance_COMM(size_t n = 0) : nlayers(n) {
    dZ_m = std::vector<double>(n, medfate::NA_DOUBLE);
    dZUp = std::vector<double>(n, medfate::NA_DOUBLE);
    dZDown = std::vector<double>(n, medfate::NA_DOUBLE);
    lambda = std::vector<double>(n, medfate::NA_DOUBLE);
    theta_micro = std::vector<double>(n, medfate::NA_DOUBLE);
    theta_b = std::vector<double>(n, medfate::NA_DOUBLE);
    theta_macro = std::vector<double>(n, medfate::NA_DOUBLE);
    theta_sat_fict = std::vector<double>(n, medfate::NA_DOUBLE);
    Ksat_b = std::vector<double>(n, medfate::NA_DOUBLE);
    Ksat_b_ms = std::vector<double>(n, medfate::NA_DOUBLE);
    Ksat = std::vector<double>(n, medfate::NA_DOUBLE);
    Ksat_ms = std::vector<double>(n, medfate::NA_DOUBLE);
    Psi = std::vector<double>(n, medfate::NA_DOUBLE);
    Psi_m = std::vector<double>(n, medfate::NA_DOUBLE);
    Psi_step = std::vector<double>(n, medfate::NA_DOUBLE);
    Psi_step_m = std::vector<double>(n, medfate::NA_DOUBLE);
    K = std::vector<double>(n, medfate::NA_DOUBLE);
    K_ms = std::vector<double>(n, medfate::NA_DOUBLE);
    Kbc = std::vector<double>(n, medfate::NA_DOUBLE);
    Kbc_ms = std::vector<double>(n, medfate::NA_DOUBLE);
    K_step = std::vector<double>(n, medfate::NA_DOUBLE);
    K_step_ms = std::vector<double>(n, medfate::NA_DOUBLE);
    C = std::vector<double>(n, medfate::NA_DOUBLE);
    C_m = std::vector<double>(n, medfate::NA_DOUBLE);
    C_step = std::vector<double>(n, medfate::NA_DOUBLE);
    C_step_m05 = std::vector<double>(n, medfate::NA_DOUBLE);
    S_macro = std::vector<double>(n, medfate::NA_DOUBLE);
    e_macro = std::vector<double>(n, medfate::NA_DOUBLE);
    S_macro_step = std::vector<double>(n, medfate::NA_DOUBLE);
    theta_macro_step = std::vector<double>(n, medfate::NA_DOUBLE);
    a = std::vector<double>(n, medfate::NA_DOUBLE);
    b = std::vector<double>(n, medfate::NA_DOUBLE);
    c = std::vector<double>(n, medfate::NA_DOUBLE);
    d = std::vector<double>(n, medfate::NA_DOUBLE);
    e = std::vector<double>(n, medfate::NA_DOUBLE);
    f = std::vector<double>(n, medfate::NA_DOUBLE);
    Kmacro_ms = std::vector<double>(n, medfate::NA_DOUBLE);
    K_step_ms05 = std::vector<double>(n, medfate::NA_DOUBLE);
    waterFluidity = std::vector<double>(n, medfate::NA_DOUBLE);
    finalSourceSinks_m3s = std::vector<double>(n, medfate::NA_DOUBLE);
    capill_below = std::vector<double>(n, medfate::NA_DOUBLE);
    drain_above = std::vector<double>(n, medfate::NA_DOUBLE);
    drain_below = std::vector<double>(n, medfate::NA_DOUBLE);
    theta_micro_step = std::vector<double>(n, medfate::NA_DOUBLE);
    lateral_flows_step_mm = std::vector<double>(n, medfate::NA_DOUBLE);
                                        
  }
};

// ----------------------------------------------------------------------------
// Soil Energy Balance Communication Structure (15 fields)
// ----------------------------------------------------------------------------

struct SoilEnergyBalance {
  size_t nlayers;
  
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
  SoilEnergyBalance(size_t n = 0) : nlayers(n) {
    auto na = std::vector<double>(n, medfate::NA_DOUBLE);
    dZ_m = dZUp = dZDown = na;
    Zup = Zdown = Zcent = na;
    a = b = c = d = e = f = na;
    k_up = k_down = tempch = na;
  }
};

#endif