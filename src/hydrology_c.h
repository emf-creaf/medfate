#include "medfate.h"
#include "modelInput_c.h"
#include "soil_c.h"
#include "communication_structures_c.h"

#ifndef HYDROLOGY_C_H
#define HYDROLOGY_C_H

struct SoilWaterBalance_RESULT {
  double localSourceSinks_mm;
  double lateralSourceSinks_mm;
  double infiltration_mm;
  double infiltrationExcess_mm;
  double saturationExcess_mm;
  double runoff_mm;
  double deepDrainage_mm;
  double capillarityRise_mm;
};

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
                               bool modifySoil);
double interceptionGashDay_c(double Rainfall, double Cm, double p, double ER);
double interceptionLiuDay_c(double Rainfall, double Cm, double p, double ER);
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
void waterInputs_c(std::vector<double>& waterInputs,
                   ModelInput& x,
                   double prec, double rainfallIntensity,
                   double pet, double tday, double rad, double elevation,
                   double Cm, double LgroundPAR, double LgroundSWR, 
                   bool modifyInput);
double microporeImbibitionRate_c(double theta_b, double theta_micro, 
                                 double D_theta_b, double D_theta_micro,
                                 double S_macro);
double rootFindingMacropores_c(double S_t, double K_up, double Ksat_ms, double Ksat_b_ms, double kin_exp,
                               double e_macro, double lambda, double dZ_m, double sourceSink_macro_m3s, double tstep, 
                               int Nmax);

void soilWaterBalance_inner_c(SoilWaterBalance_RESULT &SWBres, SoilWaterBalance_COMM &SWBcomm, Soil &soil, 
                              double rainfallInput, double rainfallIntensity, double snowmelt, const std::vector<double> &sourceSink, 
                              double runon, const std::vector<double> &lateralFlows, double waterTableDepth,
                              std::string infiltrationMode, double infiltrationCorrection, 
                              std::string soilDomains, 
                              int nsteps, int max_nsubsteps);

#endif
