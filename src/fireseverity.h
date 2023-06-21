#include <Rcpp.h>

#ifndef FIRESEVERITY_H
#define FIRESEVERITY_H
#endif
using namespace Rcpp;

double plumeTemperature(double Ib_surf, double z, double T_air = 25.0, double rho_air = 1.169);
double barkThermalDiffusivity(double fmc_bark, double rho_bark = 500.0, double T_air = 25.0);
double radialBoleNecrosis(double Ib_surf, double t_res, double bark_diffusivity,
                          double T_air = 25.0, double rho_air = 1.169, double T_necrosis = 60.0);
double leafThermalFactor(double SLA, double h = 130.0, double c = 2500.0);
double necrosisCriticalTemperature(double t_res, double tissue_factor, double T_air = 25.0, double T_necrosis = 60.0);
double necrosisHeight(double Ib_surf, double t_res, double tissue_factor, 
                      double T_air = 25.0, double rho_air = 1.169, double T_necrosis = 60.0);
