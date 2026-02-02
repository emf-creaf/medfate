#include <RcppArmadillo.h>
#include "modelInput_c.h"

#ifndef INNER_SUREAU_C_H
#define INNER_SUREAU_C_H

struct SureauOpt {
  double Lsym;
  double Ssym;
  double CLapo;
  double CTapo;
  double Lcav;
  double Scav;
  double Eord;
};
struct SureauParams { 
  int npools;
  double TPhase_gmin;
  double Q10_1_gmin;
  double Q10_2_gmin;
  double Gsw_AC_slope;
  double fTRBToLeaf;
  double C_SApoInit;
  double C_LApoInit;
  double k_SLApoInit;
  double k_CSApoInit; 
  double* k_RCApoInit;
  double* VGrhizo_kmax;
  double slope_gs;
  double P50_gs;
  double Tgs_optim;
  double Tgs_sens;
  double JarvisPAR;
  double gmin20;
  double gsMax;
  double gmin_S;
  double gsNight;
  double VCleaf_P50;
  double VCleaf_slope;
  double VCstem_P50;
  double VCstem_slope;
  double VCroot_P50;
  double VCroot_slope;
  double PiFullTurgor_Leaf;
  double epsilonSym_Leaf;
  double PiFullTurgor_Stem;
  double epsilonSym_Stem;
};

struct SureauNetwork {
  SureauParams params;
  double LAI;
  double Psi_LApo;
  double Psi_LSym;
  double Psi_RCApo;
  double Psi_SApo;
  double Psi_SSym;
  double Psi_SApo_cav;
  double Psi_LApo_cav;
  double PLC_Stem;
  double PLC_Leaf;                   
  double C_SApo;
  double C_LApo;
  double C_SSym;
  double C_LSym;
  double k_SLApo;                     
  double k_CSApo;
  double k_SSym;
  double k_LSym;
  double* PsiSoil;
  double* k_RSApo;
  double* k_SoilToStem;
  double* k_Soil;
  double k_Plant;
  double Q_SApo_sat_mmol_perLeafArea;
  double Q_LApo_sat_mmol_perLeafArea;
  double Q_SSym_sat_mmol_perLeafArea;
  double Q_LSym_sat_mmol_perLeafArea;
  double Einst;
  double Einst_SL;
  double Einst_SH;
  double Elim;
  double Elim_SL;
  double Elim_SH;
  double Emin_L;
  double Emin_L_SL;
  double Emin_L_SH;
  double Emin_S;
  int Diag_nwhile_cavit;
  double Diag_deltaRegulMax;
  double Diag_deltaPLCMax;
  double Diag_timeStepInSeconds; 
};

double PLC_derivative_c(double plc, double slope);
double PLC_c(double Pmin, double slope, double P50);
double invPLC_c(double plc, double slope, double P50);
double RWC_c(double PiFT, double Esymp, double Pmin);
double Emin_c(double gmin, double gBL, double gCrown, 
              double VPD, double airPressure =101.3);
double Turgor_c(double PiFT, double Esymp, double Rstemp);
void update_conductances_c(SureauNetwork &network);
void update_capacitances_c(SureauNetwork &network);
double gsJarvis_c(SureauParams &params, double PAR, double Temp, int option = 1);

void semi_implicit_integration_inner_c(SureauNetwork& network,
                                       double dt, 
                                       const SureauOpt& opt, 
                                       const std::string& stemCavitationRecovery, 
                                       const std::string& leafCavitationRecovery);

void copyParams_c(SureauParams& params, SureauParams& sinkParams);
void copyNetwork_c(SureauNetwork& network, SureauNetwork& sinkNetwork);
void deleteSureauNetworkPointers_c(SureauNetwork &network);

void initSureauParams_inner_c(SureauParams& params, int c,
                              InternalWater& internalWater, 
                              TranspirationParams& paramsTranspiration, 
                              WaterStorageParams& paramsWaterStorage,
                              std::vector<double>& VCroot_kmax, 
                              std::vector<double>& VGrhizo_kmax,
                              ControlParameters& control, 
                              double sapFluidityDay);

void initSureauNetwork_inner_c(SureauNetwork& network, int c, 
                               std::vector<double>& LAIphe,
                               InternalWater& internalWater, 
                               AnatomyParams& paramsAnatomy, 
                               TranspirationParams& paramsTranspiration, 
                               WaterStorageParams& paramsWaterStorage,
                               std::vector<double>& VCroot_kmax, std::vector<double>& VGrhizo_kmax,
                               std::vector<double>& PsiSoil, std::vector<double>& VG_n, std::vector<double>& VG_alpha,
                               ControlParameters& control, double sapFluidityDay);

#endif
