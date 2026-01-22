#ifndef INCBETA_C_H
#define INCBETA_C_H

const int MAXIT = 100;
const double EPS = 3.0e-7;
const double FPMIN = 1.0e-30;

double betacf_c(double a, double b, double x);
double incbeta_c(double a, double b, double x);
#endif
