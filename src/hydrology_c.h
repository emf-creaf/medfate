#include "medfate.h"
#include "soil_c.h"

#ifndef HYDROLOGY_C_H
#define HYDROLOGY_C_H

double soilEvaporationAmount_c(double DEF,double PETs, double Gsoil);
double soilEvaporation_c(Soil& soil,  
                         double snowpack, 
                         double pet, double LgroundSWR,
                         bool modifySoil);
void herbaceousTranspiration_c(std::vector<double>& EherbVec, 
                               Soil& soil, 
                               double pet, double LherbSWR, 
                               double herbLAI,
                               const std::vector<double> V,
                               bool modifySoil = true);
double interceptionGashDay_c(double Rainfall, double Cm, double p, double ER=0.05);
double interceptionLiuDay_c(double Rainfall, double Cm, double p, double ER=0.05);
double snowMelt_c(double tday, double rad, double LgroundSWR, double elevation);
double rainfallIntensity_c(int month, double prec, const std::vector<double>& rainfallIntensityPerMonth);
double infiltrationBoughton_c(double input, double Ssoil);
double infitrationGreenAmpt_c(double t, double Psi_w, double Ksat, double theta_sat, double theta_dry);
void infiltrationRepartition_c(double I, 
                               std::vector<double> &Ivec, 
                               const std::vector<double> &widths, 
                               const std::vector<double> &macro, 
                               double a, double b);
double infiltrationAmount_c(double rainfallInput, double rainfallIntensity, Soil& soil, 
                            std::string model, double K_correction);
double microporeImbibitionRate_c(double theta_b, double theta_micro, 
                                 double D_theta_b, double D_theta_micro,
                                 double S_macro);
double rootFindingMacropores_c(double S_t, double K_up, double Ksat_ms, double Ksat_b_ms, double kin_exp,
                               double e_macro, double lambda, double dZ_m, double sourceSink_macro_m3s, double tstep, 
                               int Nmax = 100);

#endif
