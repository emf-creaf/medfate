#include <vector>

#ifndef INCGAMMA_C_H
#define INCGAMMA_C_H

const double dwarf=0.0000001;
const double giant=999999999.9;
const double explow = -300.0;
const double machtol = 1.0E-15;
const double sqrttwopi=2.5066282746310005024;
const double lnsqrttwopi=0.9189385332046727418;
const double twopi=6.2831853071795864769;
const double oneoversqrtpi=0.5641895835477562869;
const double epss=1.0E-15;

std::vector<double> incgam_c(double a, double x);
double invincgam_c(double a, double p, double q);
double errorfunction_c(double x, bool erfcc, bool expo);
#endif
