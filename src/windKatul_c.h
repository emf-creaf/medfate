#include <RcppArmadillo.h>
#include "medfate.h"

#ifndef WINDKATUL_C_H
#define WINDKATUL_C_H


/* K-Epsilon Models of Katul et al. (2004)
 * Code adapted from Matlab to Rcpp
 * www: https://nicholas.duke.edu/people/faculty/katul/k_epsilon_model.htm
 */

/* (A)
 *   
 * Second order turbulence model closure constants for canopy-layer simulation.
 * Revised by csz: Fri Mar 28 2003
 */

//Rate of MKE lost by drag converted into TKE (in [0, 1])
// (roughly-) Fitted value for canopy-layer velocity profile simulation
const double Bp=1.0;
// Mixing lenght constant (Seginer, 1974; Massman and Weil, 1999)
// alpha=0.09 (Katul and Chang, 1999)
double alphaCNT=0.05;

// ASL Values for sigma_u/u*, sigma_v/u*, sigma_w/u*
// (set upper boundary condition on TKE)
const double AAu=2.3;
const double AAv=2.1;
const double AAw=1.25;
const double Aq=0.5*((AAu*AAu) + (AAv*AAv) + (AAw*AAw));

// Determine the Kolmogorov constant from Au, Av, Aw
// (constant for the turbulent viscocity)
const double Cu=1.0/(Aq*Aq);

// Von karman constant
const double kv=0.4;

// Constants for the dissipation budget (Launder and Spalding, 1974)
const double Ce1=1.44;
const double Ce2=1.92;

// Prandtl numbers
// Prandtl number (PrTKE = 1) for TKE hardwired in solver routines.
const double PrTKE=1.0;
// dissipation budget (Detering and Etling, 1974)
const double Pr=(kv*kv)/(sqrt(Cu)*(Ce2-Ce1));

// Wake TKE budget coefficients (Sanz, 2003)
const double Cg=pow(2.0/alphaCNT,2.0/3.0);

// Shortcircuit in Dissipation (linear in Bp)
const double Bd=sqrt(Cu)*Cg*Bp + 3.0/PrTKE;

// Dissipation budget (constant with alpha)
const double Ce4=PrTKE*(2.0/Pr-sqrt(Cu)/6.0*Cg*(Ce2-Ce1));                

// Ce5[=1.5]=Ce4 (Green, 1992), not [=0.6]!=Ce4 (Liu et al., 1996).
// *** Proof: see ``A note on k-epsilon modelling...'' (Sanz, 2003)
// *** The '1.5' value above for engineering (e.g. wind tunnel) flows (Cu=0.09)
const double Ce5=Ce4;

struct CanopyTurbulenceModel_RESULT {
  std::vector<double> z1;
  std::vector<double> U1;
  std::vector<double> dU1;
  std::vector<double> epsilon1;
  std::vector<double> k1;
  std::vector<double> uw1;
  std::vector<double> Lmix1;
  CanopyTurbulenceModel_RESULT(size_t N) {
    z1 = std::vector<double>(N, medfate::NA_DOUBLE);
    U1 = std::vector<double>(N, medfate::NA_DOUBLE);
    dU1 = std::vector<double>(N, medfate::NA_DOUBLE);
    epsilon1 = std::vector<double>(N, medfate::NA_DOUBLE);
    uw1 = std::vector<double>(N, medfate::NA_DOUBLE);
    Lmix1 = std::vector<double>(N, medfate::NA_DOUBLE);
  }
};
Rcpp::DataFrame copyCanopyTurbulenceModelResult_c(const CanopyTurbulenceModel_RESULT& canopyTurbulenceModel);

struct CanopyTurbulence_RESULT {
  std::vector<double> zmid;
  std::vector<double> u;
  std::vector<double> du;
  std::vector<double> epsilon;
  std::vector<double> k;
  std::vector<double> uw;
  CanopyTurbulence_RESULT(size_t N) {
    zmid = std::vector<double>(N, medfate::NA_DOUBLE);
    u = std::vector<double>(N, medfate::NA_DOUBLE);
    du = std::vector<double>(N, medfate::NA_DOUBLE);
    epsilon = std::vector<double>(N, medfate::NA_DOUBLE);
    k = std::vector<double>(N, medfate::NA_DOUBLE);
    uw = std::vector<double>(N, medfate::NA_DOUBLE);
  }
};
Rcpp::DataFrame copyCanopyTurbulenceResult_c(const CanopyTurbulence_RESULT& canopyTurbulence);

void windCanopyTurbulence_inner_c(CanopyTurbulence_RESULT& canopyTurbulence, CanopyTurbulenceModel_RESULT& canopyTurbulenceModel, 
                                  const std::vector<double>& zmid, 
                                  const std::vector<double>& LAD, 
                                  double canopyHeight,
                                  double u, double windMeasurementHeight, 
                                  std::string model);

#endif
