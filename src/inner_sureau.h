#include <Rcpp.h>

#ifndef INNER_COCHARD_H
#define INNER_COCHARD_H
#endif
using namespace Rcpp;

struct SureauParams { 
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
  double Elim;
  double Elim_SL;
  double Elim_SH;
  double Emin_L;
  double Emin_L_SL;
  double Emin_L_SH;
  double Emin_S;
  int Diag_nwhile_cavit;
  int Diag_deltaRegulMax;
  int Diag_deltaPLCMax;
  int Diag_timeStepInSeconds; 
};
void initSureauNetwork_inner(SureauNetwork &network, int c, NumericVector LAIphe,
                             DataFrame internalWater, 
                             DataFrame paramsAnatomy, DataFrame paramsTranspiration, DataFrame paramsWaterStorage,
                             NumericVector VCroot_kmax, NumericVector VGrhizo_kmax,
                             NumericVector PsiSoil, NumericVector VG_n, NumericVector VG_alpha,
                             List control, double sapFluidityDay = 1.0);
void deleteSureauNetworkPointers(SureauNetwork &network);


void innerSureau(List x, SureauNetwork* networks, List input, List output, int n, double tstep, 
                 bool verbose = false);
