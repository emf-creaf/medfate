#include <Rcpp.h>

#ifndef WOODFORMATION_H
#define WOODFORMATION_H
#endif
using namespace Rcpp;

List initialize_ring(double phi0=0.13, double pi0=-0.8, double CRD0=8.3);
void grow_ring(List ring, double psi, double Tc, 
               double Nc=8.85, double phi0=0.13, double pi0=-0.8, double CRD0=8.3,
               double Y_P=0.05, double Y_T=8, double h=0.043*1.8, double s=1.8);