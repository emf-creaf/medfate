#include "medfate.h"
#include "modelInput_c.h"
#include <RcppArmadillo.h>
#include <vector>
#include <string>
#include <limits>
#include <cmath>

#ifndef COMMUNICATION_STRUCTURES_C_H
#define COMMUNICATION_STRUCTURES_C_H

Rcpp::NumericMatrix copyNumericMatrix_c(arma::mat comm, int rows, int cols);

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