#include "medfate.h"
#include "modelInput_c.h"
#include "soil_c.h"

#ifndef HYDROLOGY_C_H
#define HYDROLOGY_C_H



struct WaterInputs_COMM {
  double rain;
  double snow;
  double interception;
  double netrain;
  double melt;
};

// ----------------------------------------------------------------------------
// Soil Water Balance Communication Structure
// ----------------------------------------------------------------------------

struct SoilWaterBalance_COMM {

  // Layer dimensions
  std::vector<double> dZ_m;
  std::vector<double> dZUp;
  std::vector<double> dZDown;
  
  // Soil hydraulic properties
  std::vector<double> lambda;
  
  // Water content states
  std::vector<double> Theta;
  std::vector<double> theta_micro;
  std::vector<double> theta_b;
  std::vector<double> theta_macro;
  std::vector<double> theta_sat_fict;
  std::vector<double> prop_saturated;
  
  // Saturated conductivity
  std::vector<double> Ksat_b;
  std::vector<double> Ksat_b_ms;
  std::vector<double> Ksat;
  std::vector<double> Ksat_ms;
  
  // Water potential
  std::vector<double> Psi;
  std::vector<double> Psi_m;
  std::vector<double> Psi_step;
  std::vector<double> Psi_step_m;
  std::vector<double> Psi_step_t1;
  std::vector<double> Psi_step_t05;
  
  // Conductivity
  std::vector<double> K;
  std::vector<double> K_ms;
  std::vector<double> Kbc;
  std::vector<double> Kbc_ms;
  std::vector<double> K_step;
  std::vector<double> K_step_ms;
  
  // Capacitance
  std::vector<double> C;
  std::vector<double> C_m;
  std::vector<double> C_step;
  std::vector<double> C_step_m;
  std::vector<double> C_step_m05;
  
  // Macro-porosity
  std::vector<double> S_macro;
  std::vector<double> e_macro;
  std::vector<double> S_macro_step;
  std::vector<double> theta_macro_step;
  std::vector<double> Kmacro_step_ms;
  
  // Polytelnyic coefficients
  std::vector<double> a;
  std::vector<double> b;
  std::vector<double> c;
  std::vector<double> d;
  std::vector<double> e;
  std::vector<double> f;
  
  // Inputs/Outputs
  std::vector<double> source_sink_def_mm;
  std::vector<double> lateral_flows_step_mm;
  std::vector<double> IVec_mm;
  std::vector<double> finalSourceSinks_m3s;
  std::vector<double> source_sink_def_m3s;
  std::vector<double> matrixImbibition_m3s;
  std::vector<double> matrixExcess_m3s;
  std::vector<double> saturated_matrix_correction_m3s;
  std::vector<double> saturated_macropore_correction_m3s;
  std::vector<double> matrix_macropore_flows_mm;
  
  // Additional properties
  std::vector<double> Kmacro_ms;
  std::vector<double> K_step_ms05;
  std::vector<double> waterFluidity;
  std::vector<double> capill_below;
  std::vector<double> drain_above;
  std::vector<double> drain_below;
  std::vector<double> theta_micro_step;
  
  // Constructor
  SoilWaterBalance_COMM(size_t n = 0) {
    dZ_m = std::vector<double>(n, medfate::NA_DOUBLE);
    dZUp = std::vector<double>(n, medfate::NA_DOUBLE);
    dZDown = std::vector<double>(n, medfate::NA_DOUBLE);
    lambda = std::vector<double>(n, medfate::NA_DOUBLE);
    
    // Water content states
    Theta = std::vector<double>(n, medfate::NA_DOUBLE);
    theta_micro = std::vector<double>(n, medfate::NA_DOUBLE);
    theta_b = std::vector<double>(n, medfate::NA_DOUBLE);
    theta_macro = std::vector<double>(n, medfate::NA_DOUBLE);
    theta_sat_fict = std::vector<double>(n, medfate::NA_DOUBLE);
    prop_saturated = std::vector<double>(n, medfate::NA_DOUBLE);
    
    // Saturated conductivity
    Ksat_b = std::vector<double>(n, medfate::NA_DOUBLE);
    Ksat_b_ms = std::vector<double>(n, medfate::NA_DOUBLE);
    Ksat = std::vector<double>(n, medfate::NA_DOUBLE);
    Ksat_ms = std::vector<double>(n, medfate::NA_DOUBLE);
    
    // Water potential
    Psi = std::vector<double>(n, medfate::NA_DOUBLE);
    Psi_m = std::vector<double>(n, medfate::NA_DOUBLE);
    Psi_step = std::vector<double>(n, medfate::NA_DOUBLE);
    Psi_step_m = std::vector<double>(n, medfate::NA_DOUBLE);
    Psi_step_t1 = std::vector<double>(n, medfate::NA_DOUBLE);
    Psi_step_t05 = std::vector<double>(n, medfate::NA_DOUBLE);
    
    // Conductivity
    K = std::vector<double>(n, medfate::NA_DOUBLE);
    K_ms = std::vector<double>(n, medfate::NA_DOUBLE);
    Kbc = std::vector<double>(n, medfate::NA_DOUBLE);
    Kbc_ms = std::vector<double>(n, medfate::NA_DOUBLE);
    K_step = std::vector<double>(n, medfate::NA_DOUBLE);
    K_step_ms = std::vector<double>(n, medfate::NA_DOUBLE);
    
    //Capacitance
    C = std::vector<double>(n, medfate::NA_DOUBLE);
    C_m = std::vector<double>(n, medfate::NA_DOUBLE);
    C_step = std::vector<double>(n, medfate::NA_DOUBLE);
    C_step_m = std::vector<double>(n, medfate::NA_DOUBLE);
    C_step_m05 = std::vector<double>(n, medfate::NA_DOUBLE);
    
    // Macro-porosity
    S_macro = std::vector<double>(n, medfate::NA_DOUBLE);
    e_macro = std::vector<double>(n, medfate::NA_DOUBLE);
    S_macro_step = std::vector<double>(n, medfate::NA_DOUBLE);
    Kmacro_step_ms = std::vector<double>(n, medfate::NA_DOUBLE);
    theta_macro_step = std::vector<double>(n, medfate::NA_DOUBLE);
    
    // Polytelnyic coefficients
    a = std::vector<double>(n, medfate::NA_DOUBLE);
    b = std::vector<double>(n, medfate::NA_DOUBLE);
    c = std::vector<double>(n, medfate::NA_DOUBLE);
    d = std::vector<double>(n, medfate::NA_DOUBLE);
    e = std::vector<double>(n, medfate::NA_DOUBLE);
    f = std::vector<double>(n, medfate::NA_DOUBLE);
    
    // Additional properties
    Kmacro_ms = std::vector<double>(n, medfate::NA_DOUBLE);
    K_step_ms05 = std::vector<double>(n, medfate::NA_DOUBLE);
    waterFluidity = std::vector<double>(n, medfate::NA_DOUBLE);
    capill_below = std::vector<double>(n, medfate::NA_DOUBLE);
    drain_above = std::vector<double>(n, medfate::NA_DOUBLE);
    drain_below = std::vector<double>(n, medfate::NA_DOUBLE);
    theta_micro_step = std::vector<double>(n, medfate::NA_DOUBLE);
    
    // Inputs/Outputs
    source_sink_def_mm = std::vector<double>(n, medfate::NA_DOUBLE);
    lateral_flows_step_mm = std::vector<double>(n, medfate::NA_DOUBLE);
    IVec_mm = std::vector<double>(n, medfate::NA_DOUBLE);
    finalSourceSinks_m3s = std::vector<double>(n, medfate::NA_DOUBLE);
    source_sink_def_m3s = std::vector<double>(n, medfate::NA_DOUBLE);
    matrixImbibition_m3s = std::vector<double>(n, medfate::NA_DOUBLE);
    matrixExcess_m3s = std::vector<double>(n, medfate::NA_DOUBLE);
    saturated_matrix_correction_m3s = std::vector<double>(n, medfate::NA_DOUBLE);
    saturated_macropore_correction_m3s = std::vector<double>(n, medfate::NA_DOUBLE);
    matrix_macropore_flows_mm = std::vector<double>(n, medfate::NA_DOUBLE);
    
  }
};


// ----------------------------------------------------------------------------
// Soil Water Balance Output Structure
// ----------------------------------------------------------------------------

struct SoilWaterBalance_RESULT {
  double localSourceSinks_mm;
  double lateralSourceSinks_mm;
  double infiltration_mm;
  double infiltrationMatrix_mm;
  double infiltrationMacropores_mm;
  double infiltrationExcess_mm;
  double infiltrationExcessMatrix_mm;
  double infiltrationExcessMacropores_mm;
  double saturationExcess_mm;
  double saturationExcessMatrix_mm;
  double saturationExcessMacropores_mm;
  double runoff_mm;
  double deepDrainage_mm;
  double drainageMatrix_mm;
  double drainageMacropores_mm;
  double capillarityRise_mm;
  double capillarityMatrix_mm;
  double capillarityMacropores_mm;
  double correction_mm;
  double correctionMatrix_mm;
  double correctionMacropores_mm;
  double volumeChange_mm;
  double matrixVolumeChange_mm;
  double macroporeVolumeChange_mm;
  double matrixMacroporeFlow_mm;
  int substeps;
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
void waterInputs_c(WaterInputs_COMM& waterInputs,
                   ModelInput& x,
                   double prec, double rainfallIntensity,
                   double pet, double tday, double rad, double elevation,
                   double Cm, double LgroundPAR, double LgroundSWR, 
                   bool modifyInput);
void agricultureWaterInputs_c(WaterInputs_COMM& waterInputs,
                              AgricultureModelInput& x,
                              double prec, double tday, double rad, double elevation,
                              double LgroundSWR, 
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
