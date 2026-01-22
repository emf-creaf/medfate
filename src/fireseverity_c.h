#ifndef FIRESEVERITY_C_H
#define FIRESEVERITY_C_H

const double temperatureNecrosis = 60.0;

double plumeTemperature_c(double Ib_surf, double z, double T_air, double rho_air);
double barkThermalDiffusivity_c(double fmc_bark, double rho_bark, double T_air);
double radialBoleNecrosis_c(double Ib_surf, double t_res, double bark_diffusivity,
                          double T_air, double rho_air);
double leafThermalFactor_c(double SLA, double h, double c);
double necrosisCriticalTemperature_c(double t_res, double tissue_factor, double T_air);
double necrosisHeight_c(double Ib_surf, double t_res, double tissue_factor, 
                      double T_air, double rho_air);

#endif
