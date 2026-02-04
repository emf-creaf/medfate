#include <RcppArmadillo.h>

#ifndef PHOTOSYNTHESIS_C_H
#define PHOTOSYNTHESIS_C_H

const double R_gas = 8.314; //(J/mol/ºK) Universal gas constant
const double O2_conc = 209.0; //mmol*mol-1 (Collatz et al. 2001)
const double quantumYield = 0.3; //mol photon * mol-1 electron
const double lightResponseCurvature = 0.9;

struct Photo{
  double Ci;
  double A;
};

struct BaldocchiPhoto{
  double Gsw, Cs, Ci, An, Ag;
};

double gammaTemp_c(double Tleaf);
double KmTemp_c(double Tleaf, double Oi = 209.0);
double VmaxTemp_c(double Vmax298, double Tleaf);
double JmaxTemp_c(double Jmax298, double Tleaf);

double gLeafBoundary_c(double u, double leafWidth, double gBound0 = 0.397);
double gCrown_c(double u, double gCrown0 = 0.150);

double electronLimitedPhotosynthesis_c(double Q, double Ci, double GT, double Jmax);
double electronLimitedPhotosynthesisDerivative_c(double Q, double Ci, double GT, double Jmax);
double rubiscoLimitedPhotosynthesis_c(double Ci, double GT, double Km, double Vmax);
double rubiscoLimitedPhotosynthesisDerivative_c(double Ci, double GT, double Km, double Vmax);
double photosynthesis_Ci_c(double Q, double Ci, double GT, double Km, double Vmax, double Jmax);

double f_c(double x, double Q, double Ca, double Gc, double GT, double Km, double Vmax, double Jmax);
double fder_c(double x, double Q, double Ca, double Gc, double GT, double Km, double Vmax, double Jmax);

void leafphotosynthesis_inner_c(Photo& photo, 
                                double Q, double Catm, double Gc, double Tleaf, double Vmax298, double Jmax298);

void photosynthesisBaldocchi_inner_c(BaldocchiPhoto &photoOut,
                                     double Q, 
                                     double Catm, 
                                     double Tleaf, 
                                     double u,
                                     double Vmax298, 
                                     double Jmax298, 
                                     double leafWidth,
                                     double Gsw_AC_slope,
                                     double Gsw_AC_intercept);
#endif

