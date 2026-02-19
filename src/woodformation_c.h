#ifndef WOODFORMATION_C_H
#define WOODFORMATION_C_H
/// Constants
const double Rn = 8.314; // The perfect gas constant
const double T0 = -273.15; // Absolute 0 temperature in degC
const double Tref = 15.0; // Reference temperature in degC

double _metR_c(double Tc, double DHa, double DSd, double DHd);
double _pi2n_c(double pi, double V, double Tc);
double _n2pi_c(double n, double V, double Tc);
double _microT_c(double Tc, double inflection, double scale);
double _divide_c(double psi, double Tc,
                 double Nc, double phi0, double pi0,
                 double Y_P, double Y_T);
double relative_expansion_rate_c(double psi, double Tc, double pi, double phi, double Y_P, double Y_T);
#endif