#define STRICT_R_HEADERS
#include <Rcpp.h>
#include "struct_sureau.h"
#include "photosynthesis.h"
#include "biophysicsutils.h"
#include "hydraulics.h"
#include "soil.h"
#include "tissuemoisture.h"
#include <meteoland.h>
using namespace Rcpp;

//Function to clone SureauParams object
void copyParams(SureauParams params, SureauParams &sinkParams) {
  sinkParams.npools = params.npools;
  sinkParams.TPhase_gmin = params.TPhase_gmin;
  sinkParams.Q10_1_gmin = params.Q10_1_gmin;
  sinkParams.Q10_2_gmin = params.Q10_2_gmin;
  sinkParams.Gsw_AC_slope = params.Gsw_AC_slope;
  sinkParams.fTRBToLeaf = params.fTRBToLeaf;
  sinkParams.C_SApoInit = params.C_SApoInit;
  sinkParams.C_LApoInit = params.C_LApoInit;
  sinkParams.k_SLApoInit = params.k_SLApoInit;
  sinkParams.k_CSApoInit = params.k_CSApoInit;

  for(int i=0;i<params.npools;i++) {
    sinkParams.k_RCApoInit[i] = params.k_RCApoInit[i];
  }
  
  sinkParams.slope_gs = params.slope_gs;
  sinkParams.P50_gs = params.P50_gs;
  sinkParams.Tgs_optim = params.Tgs_optim;
  sinkParams.Tgs_sens = params.Tgs_sens;
  sinkParams.JarvisPAR = params.JarvisPAR;
  
  sinkParams.gmin20 = params.gmin20;
  sinkParams.gsMax = params.gsMax;
  sinkParams.gmin_S = params.gmin_S;
  sinkParams.gsNight = params.gsNight;

  sinkParams.VCleaf_P50 = params.VCleaf_P50; 
  sinkParams.VCleaf_slope = params.VCleaf_slope; 
  sinkParams.VCstem_P50 = params.VCstem_P50; 
  sinkParams.VCstem_slope = params.VCstem_slope; 
  sinkParams.VCroot_P50 = params.VCroot_P50; 
  sinkParams.VCroot_slope = params.VCroot_slope; 
  sinkParams.PiFullTurgor_Leaf = params.PiFullTurgor_Leaf; 
  sinkParams.epsilonSym_Leaf = params.epsilonSym_Leaf;
  sinkParams.PiFullTurgor_Stem = params.PiFullTurgor_Stem; 
  sinkParams.epsilonSym_Stem = params.epsilonSym_Stem;
}

//Function to copy SureauNetwork object
void copyNetwork(SureauNetwork network, SureauNetwork &sinkNetwork) {
  copyParams(network.params, sinkNetwork.params);
  sinkNetwork.LAI = network.LAI;
  sinkNetwork.Psi_LApo = network.Psi_LApo;
  sinkNetwork.Psi_LSym = network.Psi_LSym;
  sinkNetwork.Psi_RCApo = network.Psi_RCApo;
  sinkNetwork.Psi_SApo = network.Psi_SApo;
  sinkNetwork.Psi_SSym = network.Psi_SSym;
  sinkNetwork.Psi_SApo_cav = network.Psi_SApo_cav;
  sinkNetwork.Psi_LApo_cav = network.Psi_LApo_cav;
  sinkNetwork.PLC_Stem = network.PLC_Stem;
  sinkNetwork.PLC_Leaf = network.PLC_Leaf;                   
  sinkNetwork.C_SApo = network.C_SApo;
  sinkNetwork.C_LApo = network.C_LApo;
  sinkNetwork.C_SSym = network.C_SSym;
  sinkNetwork.C_LSym = network.C_LSym;
  sinkNetwork.k_SLApo = network.k_SLApo;                     
  sinkNetwork.k_CSApo = network.k_CSApo;
  sinkNetwork.k_SSym = network.k_SSym;
  sinkNetwork.k_LSym = network.k_LSym;
  for(int i=0;i<network.params.npools;i++) {
    sinkNetwork.k_RSApo[i] = network.k_RSApo[i];
    sinkNetwork.k_SoilToStem[i] = network.k_SoilToStem[i];
    sinkNetwork.k_Soil[i] = network.k_Soil[i];
    sinkNetwork.PsiSoil[i] = network.PsiSoil[i];
  }
  
  sinkNetwork.k_Plant = network.k_Plant;
  
  //Water content (mmol m-2)
  sinkNetwork.Q_SApo_sat_mmol_perLeafArea = network.Q_SApo_sat_mmol_perLeafArea;
  sinkNetwork.Q_LApo_sat_mmol_perLeafArea = network.Q_LApo_sat_mmol_perLeafArea;
  sinkNetwork.Q_SSym_sat_mmol_perLeafArea = network.Q_SSym_sat_mmol_perLeafArea;
  sinkNetwork.Q_LSym_sat_mmol_perLeafArea = network.Q_LSym_sat_mmol_perLeafArea;

  sinkNetwork.Einst = network.Einst;
  sinkNetwork.Einst_SL = network.Einst_SL;
  sinkNetwork.Einst_SH = network.Einst_SH;
  sinkNetwork.Elim = network.Elim;
  sinkNetwork.Elim_SL = network.Elim_SL;
  sinkNetwork.Elim_SH = network.Elim_SH;
  sinkNetwork.Emin_L = network.Emin_L;
  sinkNetwork.Emin_L_SL = network.Emin_L_SL;
  sinkNetwork.Emin_L_SH = network.Emin_L_SH;
  sinkNetwork.Emin_S = network.Emin_S;
  sinkNetwork.Diag_nwhile_cavit = network.Diag_nwhile_cavit;
  sinkNetwork.Diag_deltaRegulMax = network.Diag_deltaRegulMax;
  sinkNetwork.Diag_deltaPLCMax = network.Diag_deltaPLCMax;
  sinkNetwork.Diag_timeStepInSeconds = network.Diag_timeStepInSeconds;
}

//Function to delete pointers
void deleteSureauNetworkPointers(SureauNetwork &network) {
  delete[] network.params.k_RCApoInit;
  delete[] network.k_Soil;
  delete[] network.PsiSoil;
  delete[] network.k_RSApo;
  delete[] network.k_SoilToStem;
}

List structToList(SureauNetwork snetwork) {
  
  List network = List::create();
  //Params
  List params = List::create();
  // CONSTANTS (Control variables)
  params.push_back(snetwork.params.npools, "npools"); 
  params.push_back(snetwork.params.TPhase_gmin, "TPhase_gmin"); 
  params.push_back(snetwork.params.Q10_1_gmin, "Q10_1_gmin"); 
  params.push_back(snetwork.params.Q10_2_gmin, "Q10_2_gmin"); 
  params.push_back(snetwork.params.Tgs_optim, "Tgs_optim"); 
  params.push_back(snetwork.params.Tgs_sens, "Tgs_sens"); 
  params.push_back(snetwork.params.JarvisPAR, "JarvisPAR"); 
  params.push_back(snetwork.params.Gsw_AC_slope, "Gsw_AC_slope");
  params.push_back(snetwork.params.fTRBToLeaf, "fTRBToLeaf");
  params.push_back(snetwork.params.C_SApoInit, "C_SApoInit");
  params.push_back(snetwork.params.C_LApoInit, "C_LApoInit");
  params.push_back(snetwork.params.k_SLApoInit, "k_SLApoInit");
  params.push_back(snetwork.params.k_CSApoInit, "k_CSApoInit");
  NumericVector k_RCApoInit(snetwork.params.npools, NA_REAL);
  for(int l=0;l < snetwork.params.npools; l++) {
    k_RCApoInit[l] = snetwork.params.k_RCApoInit[l];
  }
  params.push_back(k_RCApoInit, "k_RCApoInit"); 
  params.push_back(snetwork.params.slope_gs, "slope_gs");
  params.push_back(snetwork.params.P50_gs, "P50_gs");
  params.push_back(snetwork.params.gmin20, "gmin20"); 
  params.push_back(snetwork.params.gsMax, "gsMax"); 
  params.push_back(snetwork.params.gmin_S, "gmin_S");
  params.push_back(snetwork.params.gsNight, "gsNight"); 
  params.push_back(snetwork.params.VCleaf_P50, "VCleaf_P50"); 
  params.push_back(snetwork.params.VCleaf_slope, "VCleaf_slope"); 
  params.push_back(snetwork.params.VCstem_P50, "VCstem_P50"); 
  params.push_back(snetwork.params.VCstem_slope, "VCstem_slope"); 
  params.push_back(snetwork.params.VCroot_P50, "VCroot_P50"); 
  params.push_back(snetwork.params.VCroot_slope, "VCroot_slope"); 
  params.push_back(snetwork.params.PiFullTurgor_Leaf, "PiFullTurgor_Leaf"); 
  params.push_back(snetwork.params.epsilonSym_Leaf, "epsilonSym_Leaf"); 
  params.push_back(snetwork.params.PiFullTurgor_Stem, "PiFullTurgor_Stem"); 
  params.push_back(snetwork.params.epsilonSym_Stem, "epsilonSym_Stem"); 
  
  network.push_back(params, "params");
  
  network.push_back(snetwork.LAI, "LAI");
  network.push_back(snetwork.Psi_LApo, "Psi_LApo"); 
  network.push_back(snetwork.Psi_LSym, "Psi_LSym"); 
  network.push_back(snetwork.Psi_RCApo, "Psi_RCApo");
  network.push_back(snetwork.Psi_SApo, "Psi_SApo"); 
  network.push_back(snetwork.Psi_SSym, "Psi_SSym");
  network.push_back(snetwork.Psi_SApo_cav, "Psi_SApo_cav"); 
  network.push_back(snetwork.Psi_LApo_cav, "Psi_LApo_cav"); 
  network.push_back(snetwork.PLC_Stem, "PLC_Stem"); 
  network.push_back(snetwork.PLC_Leaf, "PLC_Leaf"); 
  network.push_back(snetwork.C_SApo, "C_SApo"); 
  network.push_back(snetwork.C_LApo, "C_LApo"); 
  network.push_back(snetwork.C_SSym, "C_SSym"); 
  network.push_back(snetwork.C_LSym, "C_LSym");
  network.push_back(snetwork.k_SLApo, "k_SLApo");
  network.push_back(snetwork.k_CSApo, "k_CSApo");
  network.push_back(snetwork.k_SSym, "k_SSym"); 
  network.push_back(snetwork.k_LSym, "k_LSym"); 
  
  NumericVector k_RSApo(snetwork.params.npools, NA_REAL);
  NumericVector k_SoilToStem(snetwork.params.npools, NA_REAL);
  NumericVector k_Soil(snetwork.params.npools, NA_REAL);
  NumericVector PsiSoil(snetwork.params.npools, NA_REAL);
  for(int l=0;l < snetwork.params.npools; l++) {
    k_RSApo[l] = snetwork.k_RSApo[l];
    k_SoilToStem[l] = snetwork.k_SoilToStem[l];
    k_Soil[l] = snetwork.k_Soil[l];
    PsiSoil[l] = snetwork.PsiSoil[l];
  }
  network.push_back(k_RSApo, "k_RSApo"); 
  network.push_back(k_SoilToStem, "k_SoilToStem"); 
  network.push_back(PsiSoil, "PsiSoil");
  network.push_back(k_Soil, "k_Soil"); 
  network.push_back(snetwork.k_Plant, "k_Plant"); 
  network.push_back(snetwork.Q_SApo_sat_mmol_perLeafArea, "Q_SApo_sat_mmol_perLeafArea"); 
  network.push_back(snetwork.Q_LApo_sat_mmol_perLeafArea, "Q_LApo_sat_mmol_perLeafArea"); 
  network.push_back(snetwork.Q_SSym_sat_mmol_perLeafArea, "Q_SSym_sat_mmol_perLeafArea"); 
  network.push_back(snetwork.Q_LSym_sat_mmol_perLeafArea, "Q_LSym_sat_mmol_perLeafArea"); 
  network.push_back(snetwork.Einst, "Einst"); 
  network.push_back(snetwork.Einst_SL, "Einst_SL");
  network.push_back(snetwork.Einst_SH, "Einst_SH"); 
  network.push_back(snetwork.Elim, "Elim"); 
  network.push_back(snetwork.Elim_SL, "Elim_SL");
  network.push_back(snetwork.Elim_SH, "Elim_SH"); 
  network.push_back(snetwork.Emin_L, "Emin_L"); 
  network.push_back(snetwork.Emin_L_SL, "Emin_L_SL"); 
  network.push_back(snetwork.Emin_L_SH, "Emin_L_SH"); 
  network.push_back(snetwork.Emin_S, "Emin_S"); 
  
  //Diagnostics
  network.push_back(snetwork.Diag_nwhile_cavit, "Diag_nwhile_cavit");
  network.push_back(snetwork.Diag_deltaRegulMax, "Diag_deltaRegulMax");
  network.push_back(snetwork.Diag_deltaPLCMax, "Diag_deltaPLCMax");
  network.push_back(snetwork.Diag_timeStepInSeconds, "Diag_timeStepInSeconds");
  return(network);
}
SureauNetwork listToStruct(List network) {
  List params = network["params"];
  //Copy from List to SureauNetwork
  SureauNetwork snetwork;
  snetwork.params.npools = params["npools"];
  snetwork.params.TPhase_gmin = params["TPhase_gmin"];
  snetwork.params.Q10_1_gmin = params["Q10_1_gmin"]; 
  snetwork.params.Q10_2_gmin = params["Q10_2_gmin"]; 
  snetwork.params.Tgs_optim = params["Tgs_optim"];
  snetwork.params.Tgs_sens = params["Tgs_sens"];
  snetwork.params.JarvisPAR = params["JarvisPAR"];
  snetwork.params.Gsw_AC_slope = params["Gsw_AC_slope"];
  snetwork.params.fTRBToLeaf = params["fTRBToLeaf"];
  snetwork.params.C_SApoInit = params["C_SApoInit"];
  snetwork.params.C_LApoInit = params["C_LApoInit"];
  snetwork.params.k_SLApoInit = params["k_SLApoInit"];
  snetwork.params.k_CSApoInit = params["k_CSApoInit"];
  snetwork.params.slope_gs = params["slope_gs"];
  snetwork.params.P50_gs = params["P50_gs"];
  snetwork.params.gmin20 = params["gmin20"];
  snetwork.params.gsMax = params["gsMax"];
  snetwork.params.gmin_S = params["gmin_S"];
  snetwork.params.gsNight = params["gsNight"];
  snetwork.params.VCleaf_P50 = params["VCleaf_P50"];
  snetwork.params.VCleaf_slope = params["VCleaf_slope"];
  snetwork.params.VCstem_P50 = params["VCstem_P50"];
  snetwork.params.VCstem_slope = params["VCstem_slope"];
  snetwork.params.VCroot_P50 = params["VCroot_P50"];
  snetwork.params.VCroot_slope = params["VCroot_slope"];
  snetwork.params.PiFullTurgor_Leaf = params["PiFullTurgor_Leaf"];
  snetwork.params.epsilonSym_Leaf = params["epsilonSym_Leaf"];
  snetwork.params.PiFullTurgor_Stem = params["PiFullTurgor_Stem"];
  snetwork.params.epsilonSym_Stem = params["epsilonSym_Stem"];
  NumericVector k_RCApoInit = Rcpp::as<Rcpp::NumericVector>(params["k_RCApoInit"]);
  snetwork.params.k_RCApoInit = new double[k_RCApoInit.size()];
  for(int l=0;l < k_RCApoInit.size(); l++) {
    snetwork.params.k_RCApoInit[l] = k_RCApoInit[l];
  }
  
  snetwork.LAI = network["LAI"];
  snetwork.Psi_LApo = network["Psi_LApo"];
  snetwork.Psi_LSym = network["Psi_LSym"];
  snetwork.Psi_RCApo = network["Psi_RCApo"];
  snetwork.Psi_SApo = network["Psi_SApo"];
  snetwork.Psi_SSym = network["Psi_SSym"];
  snetwork.Psi_SApo_cav = network["Psi_SApo_cav"];
  snetwork.Psi_LApo_cav = network["Psi_LApo_cav"];
  snetwork.PLC_Stem = network["PLC_Stem"];
  snetwork.PLC_Leaf = network["PLC_Leaf"];
  snetwork.C_SApo = network["C_SApo"];
  snetwork.C_LApo = network["C_LApo"];
  snetwork.C_SSym = network["C_SSym"];
  snetwork.C_LSym = network["C_LSym"];
  snetwork.k_SLApo = network["k_SLApo"];
  snetwork.k_CSApo = network["k_CSApo"];
  snetwork.k_SSym = network["k_SSym"];
  snetwork.k_LSym = network["k_LSym"];

  NumericVector k_RSApo = Rcpp::as<Rcpp::NumericVector>(network["k_RSApo"]);
  NumericVector k_SoilToStem = Rcpp::as<Rcpp::NumericVector>(network["k_SoilToStem"]);
  NumericVector k_Soil = Rcpp::as<Rcpp::NumericVector>(network["k_Soil"]);
  NumericVector PsiSoil = Rcpp::as<Rcpp::NumericVector>(network["PsiSoil"]);
  snetwork.k_RSApo = new double[k_RSApo.size()];
  snetwork.k_SoilToStem = new double[k_RSApo.size()];
  snetwork.k_Soil = new double[k_RSApo.size()];
  snetwork.PsiSoil = new double[k_RSApo.size()];
  for(int l=0;l < snetwork.params.npools; l++) {
    snetwork.k_RSApo[l] = k_RSApo[l];
    snetwork.k_SoilToStem[l] = k_SoilToStem[l]; 
    snetwork.k_Soil[l] = k_Soil[l];
    snetwork.PsiSoil[l] = PsiSoil[l];
  }
  
  snetwork.k_Plant = network["k_Plant"];
  snetwork.Q_SApo_sat_mmol_perLeafArea = network["Q_SApo_sat_mmol_perLeafArea"];
  snetwork.Q_LApo_sat_mmol_perLeafArea = network["Q_LApo_sat_mmol_perLeafArea"];
  snetwork.Q_SSym_sat_mmol_perLeafArea = network["Q_SSym_sat_mmol_perLeafArea"];
  snetwork.Q_LSym_sat_mmol_perLeafArea = network["Q_LSym_sat_mmol_perLeafArea"];
  
  snetwork.Einst = network["Einst"];
  snetwork.Einst_SL = network["Einst_SL"];
  snetwork.Einst_SH = network["Einst_SH"];
  snetwork.Elim = network["Elim"];
  snetwork.Elim_SL = network["Elim_SL"];
  snetwork.Elim_SH = network["Elim_SH"];
  snetwork.Emin_L = network["Emin_L"];
  snetwork.Emin_L_SL = network["Emin_L_SL"];
  snetwork.Emin_L_SH = network["Emin_L_SH"];
  snetwork.Emin_S = network["Emin_S"];

  //Diagnostics
  snetwork.Diag_nwhile_cavit = network["Diag_nwhile_cavit"];
  snetwork.Diag_deltaRegulMax = network["Diag_deltaRegulMax"];
  snetwork.Diag_deltaPLCMax = network["Diag_deltaPLCMax"];
  snetwork.Diag_timeStepInSeconds = network["Diag_timeStepInSeconds"];

  return(snetwork);
}

double PLC_derivative(double plc, double slope) {
  return(-1.0*slope/25.0 * plc/100 * (1.0 - plc/100));
}
double PLC(double Pmin, double slope, double P50){
  return (100.0 / (1.0 + exp(slope / 25.0 * (Pmin - P50))));
}
double invPLC(double plc, double slope, double P50){
  return(P50 + log((100.0/plc)-1.0)*(25.0/slope));
}
double averagePsiSigmoid(NumericVector psi, NumericVector v, double slope, double P50) {
  int nlayers = psi.size();
  NumericVector K(nlayers);
  for(int l=0;l<nlayers;l++) K[l]= PLC(psi[l], slope, P50);
  double psires =  invPLC(sum(K*v), slope, P50);
  psires = std::max(psires, -40.0); //Limits plant water potential to -40 MPa
  return(psires);
}
double RWC(double PiFT, double Esymp, double Pmin) {
  double A = std::max((-1.0 * (Pmin + PiFT - Esymp) - sqrt(pow(Pmin + PiFT - Esymp, 2.0) + 4.0 * (Pmin * Esymp))) / (2.0 * Esymp), 1.0 - PiFT / Pmin);
  return(A);
}

double Emin(double gmin, double gBL, double gCrown, 
            double VPD, double airPressure =101.3) {
  double gmintot = 1.0/(1.0/gmin+ 1.0/gBL + 1.0/gCrown);
  return(gmintot * VPD /airPressure); 
}

// # Update plant conductances
void update_conductances(SureauNetwork &network) {
  
  double* k_RCApoInit = network.params.k_RCApoInit;
  double k_CSApoInit = network.params.k_CSApoInit;
  double k_SLApoInit = network.params.k_SLApoInit;
  double k_LSym = network.k_LSym;
  double plc_leaf = network.PLC_Leaf;
  double plc_stem = network.PLC_Stem;
  
  double*  k_RSApo = network.k_RSApo;
  double*  k_SoilToStem = network.k_SoilToStem;
  double*  k_Soil = network.k_Soil;
  double*  k_RCApo = new double[network.params.npools];
  
  double k_SLApo = k_SLApoInit * (1.0 - (plc_leaf/100.0)); //Conductance from stem apo to leaf apo
  network.k_SLApo  = k_SLApo;
  double k_CSApo = k_CSApoInit * (1.0 - (plc_stem/100.0)); //conductance from root crown to stem apo
  network.k_CSApo = k_CSApo;
  double sum_k_RCApo = 0.0;
  //# calculate k_SoilToStem and k_RSApo with cavitation
  for(int i = 0;i<network.params.npools;i++) {
    k_RCApo[i] = k_RCApoInit[i] * (1.0 - (plc_stem/100.0)); // conductance from root surface to root crown
    k_RSApo[i] = 1.0/((1.0/k_RCApo[i]) + (1.0/k_CSApo)); // conductance from root surface to stem
    //# Root from root length
    k_SoilToStem[i] = 1.0/((1.0/k_Soil[i]) + (1.0/k_RSApo[i])); // # conductance from soil to stem
    sum_k_RCApo += k_RCApo[i];
  }
  delete[] k_RCApo;

  // Compute k_plant (from root to leaf) for diagnostic only
  // Rcout << " PLCstem "<< ((double) network["PLC_Stem"]) << " PLCleaf "<< ((double) network["PLC_Leaf"]) << " "<< sum(k_RSApo) << " Leaf " << k_SLApoInit << "/"<<k_SLApo << " " << k_LSym<<"\n";
  network.k_Plant =  1.0/ (1.0 /sum_k_RCApo + 1.0/k_CSApo + 1.0/k_SLApo + 1.0/k_LSym);
}

// # update symplasmic plant capacitances for Trunk and leaves
void update_capacitances(SureauNetwork &network) {
  SureauParams params = network.params;
  double dbxmin = 1.0e-100; //# NM minimal double to avoid-INF
  
  double LAI = network.LAI;
  double Psi_SSym = network.Psi_SSym;
  double Psi_LSym = network.Psi_LSym;
  double Q_LSym_sat_mmol_perLeafArea = network.Q_LSym_sat_mmol_perLeafArea;
  // double Q_SSym_sat_mmol_perLeafArea = network["Q_SSym_sat_mmol_perLeafArea"]; //not used
  
  double epsilonSym_Leaf = params.epsilonSym_Leaf;
  double PiFullTurgor_Leaf = params.PiFullTurgor_Leaf;
  double epsilonSym_Stem = params.epsilonSym_Stem;
  double PiFullTurgor_Stem = params.PiFullTurgor_Stem;
  double PsiTLP_Leaf = turgorLossPoint(PiFullTurgor_Leaf, epsilonSym_Leaf);
  double PsiTLP_Stem = turgorLossPoint(PiFullTurgor_Stem, epsilonSym_Stem);

  //#----Compute the relative water content of the symplasm----
  double RWC_LSym = 1.0 - RWC(PiFullTurgor_Leaf, epsilonSym_Leaf, Psi_LSym - dbxmin);
  //#----Compute the derivative of the relative water content of the symplasm----
  double RWC_LSym_prime;
  if(Psi_LSym > PsiTLP_Leaf) { //# FP derivative of -Pi0- Eps(1-RWC)+Pi0/RWC
    RWC_LSym_prime = RWC_LSym / (-1.0*PiFullTurgor_Leaf - Psi_LSym - epsilonSym_Leaf + 2.0 * epsilonSym_Leaf * RWC_LSym);
  } else {
    RWC_LSym_prime = -1.0*PiFullTurgor_Leaf / pow(Psi_LSym, 2.0);// # FP derivative of Pi0/Psi
  }
  //# Compute the leaf capacitance (mmol/MPa/m2_sol)
  if (LAI==0){
    network.C_LSym = 0.0;
  } else {
    network.C_LSym = Q_LSym_sat_mmol_perLeafArea * RWC_LSym_prime;
  } //# changed 25/10/2021 by NM
  
  
  //#----Stem symplasmic canopy water content----
  double RWC_SSym = 1.0 - RWC(PiFullTurgor_Stem, epsilonSym_Stem, Psi_SSym - dbxmin);
  
  //#----Compute the derivative of the relative water content of the symplasm----
  double RWC_SSym_prime;
  if (Psi_SSym > PsiTLP_Stem) {
    RWC_SSym_prime = RWC_SSym / (-1.0* PiFullTurgor_Stem - Psi_SSym - epsilonSym_Stem + 2.0 * epsilonSym_Stem * RWC_SSym);
  } else {
    RWC_SSym_prime = -1.0* PiFullTurgor_Stem / pow(Psi_SSym, 2.0);
  }
  //# Compute the capacitance (mmol/MPa/m2_leaf)
  network.C_SSym = Q_LSym_sat_mmol_perLeafArea * RWC_SSym_prime; // #  changed 25/10/2021 by NM. --> Stem capacitance per leaf area can only decrease with LAI (cannot increase when LAI<1 )
  // Rcout << "update "<< Psi_SSym <<  " " << Q_LSym_sat_mmol_perLeafArea << " " << RWC_SSym_prime << " C_SSym: " << network.C_SSym<<"\n";
  // MIQUEL - we could use instead: network["C_SSym"] = Q_SSym_sat_mmol_perLeafArea * RWC_SSym_prime; 
  network.C_SApo = params.C_SApoInit; //MIQUEL: Why are these not scaled?
  network.C_LApo = params.C_LApoInit;
}

double Turgor(double PiFT, double Esymp, double Rstemp) {
  return(-1.0*PiFT - Esymp * Rstemp);
}

//Currently not called for optimization purposes
// double regulFact(double psi, List params, String regulationType = "Sigmoid") {
//   // double stomatalClosure = NA_REAL;
//   double regulFact = NA_REAL;
//   // double regulFactPrime = NA_REAL;
//   
//   if(regulationType == "PiecewiseLinear") {
//     double PsiStartClosing = params["PsiStartClosing"];
//     double PsiClose = params["PsiClose"];
//     if (psi > PsiStartClosing) {
//       // stomatalClosure = 0.0;
//       regulFact = 1.0;
//       // regulFactPrime = 0.0;
//     } else if (psi > PsiClose) {
//       // stomatalClosure = 1.0;
//       regulFact = (psi - PsiClose) / (PsiStartClosing - PsiClose);
//       // regulFactPrime = 1.0 / (PsiStartClosing - PsiClose);
//     } else {
//       // stomatalClosure = 2.0;
//       regulFact = 0.0;
//       // regulFactPrime = 0.0;
//     }
//   } else if (regulationType == "Sigmoid") {
//     double slope_gs = params["slope_gs"];
//     double P50_gs = params["P50_gs"];
//     double PL_gs = 1.0 / (1.0 + exp(slope_gs / 25.0 * (psi - P50_gs)));
//     regulFact = 1.0 - PL_gs;
//     // double al = slope_gs / 25.0;
//     // regulFactPrime = al * PL_gs * regulFact;
//   } else if (regulationType == "Turgor") {
//     double turgorPressureAtGsMax = params["turgorPressureAtGsMax"];
//     double epsilonSym_Leaf = params["epsilonSym_Leaf"];
//     double PiFullTurgor_Leaf = params["PiFullTurgor_Leaf"];
//     double rs = RWC(PiFullTurgor_Leaf, epsilonSym_Leaf, psi);
//     double turgor = Turgor(PiFullTurgor_Leaf, epsilonSym_Leaf, rs);
//     turgorPressureAtGsMax = Turgor(PiFullTurgor_Leaf, epsilonSym_Leaf, 0.0);
//     regulFact = std::max(0.0, std::min(1.0, turgor / turgorPressureAtGsMax));
//     // Rcout<< psi << " "<< turgor<< " "<< regulFact << "\n";
//   }
//   return(regulFact);
// }

//# stomatal conductance calculation with Jarvis type formulations
double gsJarvis(SureauParams &params, double PAR, double Temp, int option = 1){
  double JarvisPAR = params.JarvisPAR;
  double gsMax = params.gsMax;
  double gsNight = params.gsNight;
  double Tgs_optim = params.Tgs_optim;
  double Tgs_sens = params.Tgs_sens;
  double gsMax2, gsNight2;
  if (option == 1) { //# temperature effect on gs
    double tempEff = 1.0/(1.0 + pow((Temp - Tgs_optim)/Tgs_sens, 2.0));
    gsMax2    = std::max(0.0, gsMax * tempEff);
    gsNight2  = std::max(0.0, gsNight * tempEff);
  } else {
    gsMax2    = std::max(0.0, gsMax);
    gsNight2  = std::max(0.0, gsNight);
  }
  double Q = irradianceToPhotonFlux(PAR); //From W m-2 to micromol s-1 m-2
  double gs_bound = gsNight2 + (gsMax2 - gsNight2) * (1.0 - exp(-1.0*JarvisPAR*Q));
  return(gs_bound);
}

void initSureauParams_inner(SureauParams &params, int c,
                            DataFrame internalWater, 
                            DataFrame paramsTranspiration, DataFrame paramsWaterStorage,
                            NumericVector VCroot_kmax, 
                            List control, double sapFluidityDay = 1.0) {
  String stomatalSubmodel = control["stomatalSubmodel"];
  bool soilDisconnection = control["soilDisconnection"];
  
  //This will change depending on the layers connected
  params.npools = VCroot_kmax.size();
  params.TPhase_gmin = control["TPhase_gmin"];
  params.Q10_1_gmin = control["Q10_1_gmin"];
  params.Q10_2_gmin = control["Q10_2_gmin"];
  if(stomatalSubmodel=="Jarvis") {
    NumericVector Gs_Toptim = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["Gs_Toptim"]);
    NumericVector Gs_Tsens = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["Gs_Tsens"]);
    params.Tgs_optim = Gs_Toptim[c];
    params.Tgs_sens = Gs_Tsens[c];
    params.JarvisPAR = control["JarvisPAR"];
  } else {
    NumericVector Gsw_AC_slope = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["Gsw_AC_slope"]);
    params.Gsw_AC_slope = Gsw_AC_slope[c];
  }
  params.fTRBToLeaf = control["fTRBToLeaf"];//ratio of bark area to leaf area
  params.C_SApoInit = control["C_SApoInit"]; //Maximum capacitance of the stem apoplasm
  params.C_LApoInit = control["C_LApoInit"]; //Maximum capacitance of the leaf apoplasm
  NumericVector VCleafapo_kmax = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["VCleafapo_kmax"]);
  NumericVector VCstem_kmax = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["VCstem_kmax"]);
  params.k_SLApoInit = sapFluidityDay*VCleafapo_kmax[c]; //Maximum conductance from trunk apoplasm to leaf apoplasm
  params.k_CSApoInit = sapFluidityDay*VCstem_kmax[c]; //Maximum conductance from root crown to stem apoplasm
  // Rcout << "par init "<< c<< " VCstem_kmax[c]: "<<VCstem_kmax[c] << " sapfluidity: "<< sapFluidityDay<<" k_CSApoInit: " << params.k_CSApoInit<< "\n";
  params.k_RCApoInit = new double [params.npools];
  for(int l = 0;l<params.npools;l++) params.k_RCApoInit[l] = sapFluidityDay*VCroot_kmax[l];
  NumericVector Gs_P50 = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["Gs_P50"]);
  NumericVector Gs_slope = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["Gs_slope"]);
  params.slope_gs = Gs_slope[c];
  params.P50_gs = Gs_P50[c];
  
  //PLANT RELATED PARAMETERS
  NumericVector Gswmax = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["Gswmax"]);
  NumericVector Gswmin = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["Gswmin"]);
  NumericVector LeafPI0 = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["LeafPI0"]);
  NumericVector LeafEPS = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["LeafEPS"]);
  NumericVector StemPI0 = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["StemPI0"]);
  NumericVector StemEPS = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["StemEPS"]);
  NumericVector StemSympPsiVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["StemSympPsi"]);
  NumericVector LeafSympPsiVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["LeafSympPsi"]);
  NumericVector VCleaf_P50 = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["VCleaf_P50"]);
  NumericVector VCleaf_slope = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["VCleaf_slope"]);
  NumericVector VCstem_P50 = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["VCstem_P50"]);
  NumericVector VCstem_slope = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["VCstem_slope"]);
  NumericVector VCroot_P50 = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["VCroot_P50"]);
  NumericVector VCroot_slope = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["VCroot_slope"]);
  
  double gmin20 = Gswmin[c]*1000.0; //Leaf cuticular transpiration
  bool leafCuticularTranspiration = control["leafCuticularTranspiration"];
  if(soilDisconnection) gmin20 = gmin20*RWC(LeafPI0[c], LeafEPS[c], LeafSympPsiVEC[c]);
  if(!leafCuticularTranspiration) gmin20 = 0.0;
  params.gmin20 = gmin20; //mol to mmol 
  params.gsMax = Gswmax[c]*1000.0; //mol to mmol 
  bool stemCuticularTranspiration = control["stemCuticularTranspiration"];
  double gmin_S = Gswmin[c]*1000.0;  // gmin for stem equal to leaf gmin, in mmol
  if(soilDisconnection) gmin_S = gmin_S*RWC(StemPI0[c], StemEPS[c], StemSympPsiVEC[c]);
  if(!stemCuticularTranspiration) gmin_S = 0.0;
  params.gmin_S = gmin_S;
  double gs_NightFrac = control["gs_NightFrac"];
  params.gsNight = gs_NightFrac*Gswmax[c]*1000.0; 
  
  params.VCleaf_P50 = VCleaf_P50[c]; 
  params.VCleaf_slope = VCleaf_slope[c]; 
  params.VCstem_P50 = VCstem_P50[c]; 
  params.VCstem_slope = VCstem_slope[c]; 
  params.VCroot_P50 = VCroot_P50[c]; 
  params.VCroot_slope = VCroot_slope[c]; 
  params.PiFullTurgor_Leaf = LeafPI0[c]; 
  params.epsilonSym_Leaf = LeafEPS[c]; 
  params.PiFullTurgor_Stem = StemPI0[c]; 
  params.epsilonSym_Stem = StemEPS[c]; 
  
}
void initSureauNetwork_inner(SureauNetwork &network, int c, NumericVector LAIphe,
                             DataFrame internalWater, 
                             DataFrame paramsAnatomy, DataFrame paramsTranspiration, DataFrame paramsWaterStorage,
                             NumericVector VCroot_kmax, NumericVector VGrhizo_kmax,
                             NumericVector PsiSoil, NumericVector VG_n, NumericVector VG_alpha,
                             List control, double sapFluidityDay = 1.0) {
  
  String stomatalSubmodel = control["stomatalSubmodel"];
  
  //Root distribution input
  NumericVector Einst = Rcpp::as<Rcpp::NumericVector>(internalWater["Einst"]);
  NumericVector Elim = Rcpp::as<Rcpp::NumericVector>(internalWater["Elim"]);
  NumericVector Emin_L = Rcpp::as<Rcpp::NumericVector>(internalWater["Emin_L"]);
  NumericVector Emin_S = Rcpp::as<Rcpp::NumericVector>(internalWater["Emin_S"]);
  NumericVector StemPLCVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["StemPLC"]);
  NumericVector LeafPLCVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["LeafPLC"]);
  NumericVector StemPsiVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["StemPsi"]);
  NumericVector LeafPsiVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["LeafPsi"]);
  NumericVector StemSympPsiVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["StemSympPsi"]);
  NumericVector LeafSympPsiVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["LeafSympPsi"]);
  NumericVector RootCrownPsiVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["RootCrownPsi"]);
  
  NumericVector Vmax298 = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["Vmax298"]);
  NumericVector Jmax298 = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["Jmax298"]);
  NumericVector VCleafapo_kmax = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["VCleafapo_kmax"]);
  NumericVector VCleaf_P50 = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["VCleaf_P50"]);
  NumericVector VCleaf_slope = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["VCleaf_slope"]);
  NumericVector VCstem_kmax = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["VCstem_kmax"]);
  NumericVector VCstem_P50 = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["VCstem_P50"]);
  NumericVector VCstem_slope = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["VCstem_slope"]);
  NumericVector VCroot_P50 = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["VCroot_P50"]);
  NumericVector VCroot_slope = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["VCroot_slope"]);
  NumericVector kstem_symp = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["kstem_symp"]);
  NumericVector kleaf_symp = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["kleaf_symp"]);
  NumericVector Gswmax = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["Gswmax"]);
  NumericVector Gswmin = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["Gswmin"]);
  
  NumericVector Gs_P50 = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["Gs_P50"]);
  NumericVector Gs_slope = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["Gs_slope"]);
  
  
  NumericVector StemPI0 = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["StemPI0"]);
  NumericVector StemEPS = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["StemEPS"]);
  NumericVector StemAF = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["StemAF"]);
  NumericVector Vsapwood = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["Vsapwood"]); //l·m-2 = mm
  NumericVector LeafPI0 = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["LeafPI0"]);
  NumericVector LeafEPS = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["LeafEPS"]);
  NumericVector LeafAF = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["LeafAF"]);
  NumericVector Vleaf = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["Vleaf"]); //l·m-2 = mm
  
  bool soilDisconnection = control["soilDisconnection"];
  
  //Params
  initSureauParams_inner(network.params, c, internalWater, 
                         paramsTranspiration, paramsWaterStorage,
                         VCroot_kmax, 
                         control, sapFluidityDay);
  
  //LAI
  network.LAI = LAIphe[c];
  
  //Water potentials
  network.Psi_LApo = LeafPsiVEC[c]; 
  network.Psi_LSym = LeafSympPsiVEC[c];
  network.Psi_RCApo = RootCrownPsiVEC[c];

  network.Psi_SApo = StemPsiVEC[c]; 
  network.Psi_SSym = StemSympPsiVEC[c];
  network.Psi_SApo_cav = std::min(0.0, invPLC(StemPLCVEC[c]*100.0, VCstem_slope[c], VCstem_P50[c])); //Sureau operates with %
  network.Psi_LApo_cav = std::min(0.0, invPLC(LeafPLCVEC[c]*100.0, VCleaf_slope[c], VCleaf_P50[c])); //Sureau operates with %
  
  //PLC levels
  network.PLC_Stem = StemPLCVEC[c]*100.0; //Sureau operates with %
  network.PLC_Leaf = LeafPLCVEC[c]*100.0; //Sureau operates with %
  //Capacitances (mmol m-2 MPa-1)
  network.C_SApo = NA_REAL; //Capacitance of the stem apoplasm
  network.C_LApo = NA_REAL; //Capacitance of the leaf apoplasm
  network.C_SSym = NA_REAL; //Capacitance of the stem symplasm (HOW TO ESTIMATE THEM?)
  network.C_LSym = NA_REAL; //Capacitance of the leaf symplasm (HOW TO ESTIMATE THEM?)
  //Conductances (mmol m-2 MPa-1 s-1)
  network.k_SLApo = NA_REAL; //Conductance from trunk apoplasm to leaf apoplasm
  network.k_CSApo = NA_REAL; //Conductance from root crown to trunk apoplasm
  network.k_SSym = sapFluidityDay*kstem_symp[c]; //Conductance from trunk apoplasm to trunk symplasm (CONTROL PARAMETER?)
  network.k_LSym = sapFluidityDay*kleaf_symp[c]; //Conductance from leaf apoplasm to leaf symplasm
  
  network.k_RSApo = new double[network.params.npools];
  network.k_SoilToStem = new double[network.params.npools];
  network.k_Soil = new double[network.params.npools];
  network.PsiSoil = new double[network.params.npools];
  for(int l=0;l < network.params.npools; l++) {
    network.PsiSoil[l] = PsiSoil[l];
    network.k_SoilToStem[l] = NA_REAL; //Conductance from soil to trunk apoplasm
    network.k_RSApo[l] = NA_REAL; //Conductance from rhizosphere surface to trunk apoplasm
    double k_soil_cl = VGrhizo_kmax[l];
    if(soilDisconnection) k_soil_cl = k_soil_cl*(1.0 - PLC(PsiSoil[l], Gs_slope[c], Gs_P50[c])/100.0);
    network.k_Soil[l] = vanGenuchtenConductance(PsiSoil[l],
                                        k_soil_cl, 
                                        VG_n[l], VG_alpha[l]); 
  }
  network.k_Plant = NA_REAL; //Whole-plant conductance  = Plant_kmax
  
  //Water content (mmol m-2)
  double l2mmol = 1.0e6/18.0;
  network.Q_SApo_sat_mmol_perLeafArea = Vsapwood[c]*StemAF[c]*l2mmol; //Water content in stem apoplasm
  network.Q_LApo_sat_mmol_perLeafArea = Vleaf[c]*LeafAF[c]*l2mmol; //Water content in leaf apoplasm
  network.Q_SSym_sat_mmol_perLeafArea = Vsapwood[c]*(1.0 - StemAF[c])*l2mmol; //Water content in stem symplasm
  network.Q_LSym_sat_mmol_perLeafArea = Vleaf[c]*(1.0 - LeafAF[c])*l2mmol; //Water content in leaf symplasm
  
  //Flows (mmol m-2 s-1)
  network.Einst = Einst[c]; //Total transpiration
  network.Einst_SL = NA_REAL; //Total transpiration (sunlit leaves)
  network.Einst_SH = NA_REAL; //Total transpiration (shade leaves)
  network.Elim = Elim[c]; //Stomatal transpiration
  network.Elim_SL = NA_REAL; //Stomatal transpiration (sunlit leaves)
  network.Elim_SH = NA_REAL; //Stomatal transpiration (shade leaves)
  network.Emin_L = Emin_L[c]; //Leaf cuticular transpiration
  network.Emin_L_SL = NA_REAL; //Leaf cuticular transpiration (sunlit leaves)
  network.Emin_L_SH = NA_REAL; //Leaf cuticular transpiration (shade leaves)
  network.Emin_S = Emin_S[c]; //Stem cuticular transpiration
  
  //Diagnostics
  network.Diag_nwhile_cavit = NA_INTEGER;
  network.Diag_deltaRegulMax = NA_REAL;
  network.Diag_deltaPLCMax = NA_REAL;
  network.Diag_timeStepInSeconds = NA_REAL;
  
  // Update plant conductances and capacitances according to network status
  update_conductances(network);
  update_capacitances(network);
}


List initSureauNetwork(int c, NumericVector LAIphe,
                       DataFrame internalWater, 
                       DataFrame paramsAnatomy, DataFrame paramsTranspiration, DataFrame paramsWaterStorage,
                       NumericVector VCroot_kmax, NumericVector VGrhizo_kmax,
                       NumericVector PsiSoil, NumericVector VG_n, NumericVector VG_alpha,
                       List control, double sapFluidityDay = 1.0) {
  
  SureauNetwork snetwork;
  initSureauNetwork_inner(snetwork, c, LAIphe,
                          internalWater, 
                          paramsAnatomy, paramsTranspiration, paramsWaterStorage,
                          VCroot_kmax, VGrhizo_kmax,
                          PsiSoil, VG_n, VG_alpha,
                          control, sapFluidityDay);
  List network = structToList(snetwork);
  deleteSureauNetworkPointers(snetwork);
  return(network);
}

//' Sureau-ECOS inner functions for testing only
//' 
//' Function \code{initSureauNetworks} initializes hydraulic networks for all plant cohorts in x
//' Function \code{semi_implicit_integration} updates water potentials and cavitation across the hydraulic network
//' 
//' @param x An object of class \code{\link{spwbInput}} or \code{\link{growthInput}} created using \code{transpirationMode = "Sureau"}.
//'  
//' @return Function \code{initSureauNetworks} returns a vector of length equal to the number of cohorts. Each element is a list with Sureau-ECOS parameters.
//' Function \code{semi_implicit_integration} does not return anything, but modifies input parameter \code{network}.
//' 
//' @author
//' \itemize{
//'   \item{Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF}
//'   \item{Nicolas Martin-StPaul, URFM-INRAE}
//' }
//' 
//' @references
//' Ruffault J, Pimont F, Cochard H, Dupuy JL, Martin-StPaul N (2022) 
//' SurEau-Ecos v2.0: a trait-based plant hydraulics model for simulations of plant water status and drought-induced mortality at the ecosystem level.
//' Geoscientific Model Development 15, 5593-5626 (doi:10.5194/gmd-15-5593-2022).
//' 
//' 
//' @seealso  \code{\link{spwb}}
//' 
//' @name sureau_ecos
//' @keywords internal
// [[Rcpp::export("initSureauNetworks")]]
List initSureauNetworks(List x) {
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  NumericVector LAIphe = Rcpp::as<Rcpp::NumericVector>(above["LAI_expanded"]);
  
  List belowLayers = Rcpp::as<Rcpp::List>(x["belowLayers"]);
  NumericMatrix VCroot_kmax= Rcpp::as<Rcpp::NumericMatrix>(belowLayers["VCroot_kmax"]);
  NumericMatrix VGrhizo_kmax= Rcpp::as<Rcpp::NumericMatrix>(belowLayers["VGrhizo_kmax"]);
  DataFrame internalWater = Rcpp::as<Rcpp::DataFrame>(x["internalWater"]);
  DataFrame paramsAnatomy = Rcpp::as<Rcpp::DataFrame>(x["paramsAnatomy"]);
  DataFrame paramsTranspiration = Rcpp::as<Rcpp::DataFrame>(x["paramsTranspiration"]);
  DataFrame paramsWaterStorage = Rcpp::as<Rcpp::DataFrame>(x["paramsWaterStorage"]);
  
  List control = x["control"];
  DataFrame soil = Rcpp::as<Rcpp::DataFrame>(x["soil"]);
  NumericVector psiSoil = psi(soil, "VG");
  NumericVector VG_n = Rcpp::as<Rcpp::NumericVector>(soil["VG_n"]);
  NumericVector VG_alpha = Rcpp::as<Rcpp::NumericVector>(soil["VG_alpha"]);
  
  int numCohorts = internalWater.nrow();
  List networks(numCohorts);
  for(int c = 0;c<numCohorts;c++) {
    networks[c] = initSureauNetwork(c, LAIphe,
                                     internalWater, 
                                     paramsAnatomy, paramsTranspiration, paramsWaterStorage,
                                     VCroot_kmax(c,_), VGrhizo_kmax(c,_),
                                     psiSoil, VG_n, VG_alpha, 
                                     control, 1.0);
  }
  networks.attr("names") = above.attr("row.names");
  return(networks);
}


void calculateRhizoPsi_inner(int c, int nlayers,
                             SureauNetwork network, NumericMatrix RhizoPsiMAT,
                             LogicalMatrix layerConnected, 
                             List RHOP, List layerConnectedPools,
                             bool plantWaterPools) {
  int numCohorts = layerConnected.nrow();
  double* k_SoilToStem = network.k_SoilToStem;
  double* k_Soil = network.k_Soil;
  double* PsiSoil = network.PsiSoil;
  SureauParams params = network.params;
  double VCroot_slope = params.VCroot_slope;
  double VCroot_P50 = params.VCroot_P50;
  
  double Psi_SApo = network.Psi_SApo;
  if(!plantWaterPools) {
    int cl = 0;
    for(int l=0;l<nlayers;l++) {
      if(layerConnected(c,l)) {
        double fluxSoilToStem_mmolm2s = k_SoilToStem[cl]*(PsiSoil[cl] - Psi_SApo);
        RhizoPsiMAT(c,l) = PsiSoil[cl] - fluxSoilToStem_mmolm2s/k_Soil[cl];
        // Rcout<< c << l <<" "<< fluxSoilToStem_mmolm2s<<" " << PsiSoil[cl]<< " "<< RhizoPsiMAT(c,l)<<"\n";
        cl++;
      } 
    }
  } else {
    NumericMatrix RHOPcoh = Rcpp::as<Rcpp::NumericMatrix>(RHOP[c]);
    LogicalMatrix layerConnectedCoh = Rcpp::as<Rcpp::LogicalMatrix>(layerConnectedPools[c]);
    NumericVector rplv(numCohorts,NA_REAL);
    NumericVector vplv(numCohorts,NA_REAL);
    int cl = 0;
    for(int l=0;l<nlayers;l++) {
      int clj = 0;
      for(int j=0;j<numCohorts;j++) {
        if(layerConnectedCoh(j,l)) {
          double fluxSoilToStem_mmolm2s = k_SoilToStem[cl]*(PsiSoil[cl] - Psi_SApo);
          rplv[clj] = PsiSoil[cl] - fluxSoilToStem_mmolm2s/k_Soil[cl];
          vplv[clj] = RHOPcoh(j,l);
          cl++;
          clj++;
        }
      }
      NumericVector pv(clj,NA_REAL);
      NumericVector vv(clj,NA_REAL);
      for(int j=0;j<clj;j++) {
        pv[j] = rplv[j];
        vv[j] = vplv[j];
      }
      RhizoPsiMAT(c,l) = averagePsiSigmoid(pv, vv, VCroot_slope, VCroot_P50);
    }
  }
}

void semi_implicit_integration_inner(SureauNetwork &network,
                                     double dt, NumericVector opt, 
                                     String stemCavitationRecovery = "annual", 
                                     String leafCavitationRecovery = "total") {
  
  SureauParams params = network.params;
  double* PsiSoil = network.PsiSoil;
  
  // Step 1. Initializing current time step according to computation options (FP)
  double dbxmin = 1.0e-100; // FP minimal double to avoid 0/0
  double Psi_LApo_n = network.Psi_LApo;
  double Psi_SApo_n = network.Psi_SApo;
  double Psi_LSym_n = network.Psi_LSym;
  double Psi_SSym_n = network.Psi_SSym;
  double Psi_LApo_cav = network.Psi_LApo_cav;
  double Psi_SApo_cav = network.Psi_SApo_cav;
  
  //Conductances
  double K_SL = network.k_SLApo;
  double k_SSym = network.k_SSym;
  double k_LSym = network.k_LSym;
  double c_LSym = network.C_LSym;
  double c_SSym = network.C_SSym;
  double c_LApo = network.C_LApo;
  double c_SApo = network.C_SApo;
  double PLC_Leaf = network.PLC_Leaf;
  double PLC_Stem = network.PLC_Stem;
  double Q_LApo_sat_mmol_perLeafArea = network.Q_LApo_sat_mmol_perLeafArea;
  double Q_SApo_sat_mmol_perLeafArea = network.Q_SApo_sat_mmol_perLeafArea;
  
  //Modifiers
  double Lsym = opt["Lsym"];
  double Ssym = opt["Ssym"];
  double CLapo = opt["CLapo"];
  double CTapo = opt["CTapo"];
  double Lcav = opt["Lcav"];
  double Scav = opt["Scav"];
  
  
  //Apply modifiers
  double K_LSym = Lsym * k_LSym;   
  double K_SSym = Ssym * k_SSym;   
  double C_LSym = Lsym * c_LSym;   
  double C_SSym = Ssym * c_SSym;   
  double C_LApo = CLapo * c_LApo; 
  double C_SApo = CTapo * c_SApo; 
  
  double E_nph = network.Elim; // Leaf stomatal transpiration
  double Emin_L_nph = network.Emin_L; //Leaf cuticular transpiration
  double Emin_S_nph = network.Emin_S; //Stem cuticular transpiration
  
  
  
  
  double VCleaf_slope = params.VCleaf_slope;
  double VCstem_slope = params.VCstem_slope;
  double VCleaf_P50 = params.VCleaf_P50;
  double VCstem_P50 = params.VCstem_P50;
  
  //Compute K_L_Cav et K_S_Cav
  double PLC_prime_L = PLC_derivative(PLC_Leaf, VCleaf_slope);
  double K_L_Cav = -1.0 * Lcav * Q_LApo_sat_mmol_perLeafArea * PLC_prime_L / dt;  // avec WBveg$Q_LSym_sat en l/m2 sol # changed by NM (25/10/2021)
  double PLC_prime_S = PLC_derivative(PLC_Stem, VCstem_slope);
  double K_S_Cav = -1.0 * Scav * Q_SApo_sat_mmol_perLeafArea * PLC_prime_S / dt;  // opt$Scav * WBveg$K_S_Cav #FP corrected a bug sign herehanged by NM (25/10/2021)
  // Rcout<< "0 "<< PLC_prime_L << " "<<K_L_Cav<<" "<<PLC_prime_S<< " "<< K_S_Cav<<"\n";
  
  // Step 2. While loop in order to decide if cavitation or not :
  //  In order to account for the cavitation that occurs only when potentials go below their lowest value "cav" (formerly called "mem" in an earlier version)
  // the following computations are done trying sequentially the resolutions of LApo and TApo eventually activating
  // the appropriate cavitation events when needed (starting assuming no cavit at all...)
  // in case of computational problem, the last case assume no cavitation flux
  bool LcavitWellComputed = false; //initialized to false
  bool ScavitWellComputed = false;
  
  NumericVector delta_L_cavs, delta_S_cavs;
  if ((Lcav==0.0) && (Scav==0.0)) { // no cavitation flux computed
    delta_L_cavs = NumericVector::create(0.0);
    delta_S_cavs = NumericVector::create(0.0);
  } else if ((Lcav==0.0) && (Scav==1.0)) {// Scav only
    delta_L_cavs=NumericVector::create(0.0,0.0,0.0);
    delta_S_cavs=NumericVector::create(0.0,1.0,0.0);
  } else if ((Lcav==1.0) && (Scav==0.0)) {// Lcav only
    delta_L_cavs=NumericVector::create(0.0,1.0,0.0);
    delta_S_cavs=NumericVector::create(0.0,0.0,0.0);
  } else { //#Lcav=1 and Scav=1
    delta_L_cavs=NumericVector::create(0.0,1.0,0.0,1.0,0.0); // the fifth case is here in case no solution with others...
    delta_S_cavs=NumericVector::create(0.0,0.0,1.0,1.0,0.0);
  }
  
  double k_CSApo = network.k_CSApo;
  double* k_SoilToStem = network.k_SoilToStem;
  
  double alpha, Psi_td, Psi_RCApo_np1, Psi_LApo_np1, Psi_SApo_np1, Psi_LSym_np1, Psi_SSym_np1;
  double psirefL, psirefS;
  
  int nwhilecomp = 0; // # count the number of step in while loop (if more than 4 no solution and warning)
  while (((!LcavitWellComputed)||(!ScavitWellComputed)) && (nwhilecomp<delta_L_cavs.size())) {
    double delta_L_cav = delta_L_cavs[nwhilecomp];
    double delta_S_cav = delta_S_cavs[nwhilecomp];
    
    //# 2.1 LApo
    alpha = exp(-1.0*(K_SL+K_LSym+delta_L_cav*K_L_Cav)/C_LApo*dt);
    Psi_td = (K_SL*Psi_SApo_n + K_LSym*Psi_LSym_n + delta_L_cav*K_L_Cav*Psi_LApo_cav)/(K_SL + K_LSym+delta_L_cav*K_L_Cav + dbxmin);// # dbxmin to avoid 0/0
    Psi_LApo_np1 = alpha * Psi_LApo_n + (1.0 - alpha) * Psi_td;
    
    //# 2.2. SApo
    double k_sum = 0.0, kpsi_sum = 0.0;
    for(int l=0;l<params.npools;l++) {
      k_sum +=k_SoilToStem[l];
      kpsi_sum += (k_SoilToStem[l] * PsiSoil[l]);
    }
    alpha = exp(-1.0*(K_SL+K_SSym + k_sum +delta_S_cav*K_S_Cav)/C_SApo*dt);
    Psi_td = (K_SL*Psi_LApo_n + K_SSym*Psi_SSym_n + kpsi_sum + delta_S_cav*K_S_Cav*Psi_SApo_cav)/(K_SL + K_SSym + k_sum + delta_S_cav*K_S_Cav + dbxmin);// # dbxmin to avoid 0/0
    Psi_SApo_np1 = alpha * Psi_SApo_n + (1.0 - alpha) * Psi_td;
    
    //Update Psi_RCApo from flow and conductance
    double fluxSoilToStem_mmolm2s = 0.0;
    for(int cl=0;cl<params.npools;cl++) fluxSoilToStem_mmolm2s += k_SoilToStem[cl]*(PsiSoil[cl] - Psi_SApo_np1);
    Psi_RCApo_np1 = Psi_SApo_np1 + fluxSoilToStem_mmolm2s/k_CSApo;
    
    // Rcout<< kpsi_sum<<" "<< k_sum << " " << alpha << " " << Psi_td << "\n";
    //# 2.3 Compute Psi_SApo_np1 (Explicit approach only)
    //# 2.4 check if cavitation is well computed according to delta_cav, np1 and "cav"
    LcavitWellComputed = (delta_L_cav==(Psi_LApo_np1 < Psi_LApo_cav)) || (Lcav==0.0);
    ScavitWellComputed = (delta_S_cav==(Psi_SApo_np1 < Psi_SApo_cav)) || (Scav==0.0);
    
    nwhilecomp = nwhilecomp + 1;
    
    // Rcout<< nwhilecomp << " "<< Psi_LApo_np1 << " "<<Psi_LApo_cav<<" "<<Psi_SApo_np1<< " "<< Psi_SApo_cav<<"\n";
    
    if ((delta_L_cavs.size() > 1) && (nwhilecomp==delta_L_cavs.size())) { //# we tried the normal cases and the computation is still not ok so we have done a last one desactivating cavitation water source (delta_cav=0)
      Rcerr << "water flux due to Cavitation ignored with time step, no solution from the implicit solver="<<dt<<"\n";
    }
  } //# end of the while loop with check on cavitation options
  
  network.Diag_nwhile_cavit = nwhilecomp;  // # Diagnostic step to track cavit event and eventual errors (corresponding to nwhilecomp==5)
  
  //# Step 3. Compute Psi_Symp_np1 (L and S)
  alpha = exp(-1.0*K_LSym/C_LSym*dt);
  Psi_td = (K_LSym*Psi_LApo_n - (E_nph + Emin_L_nph))/(K_LSym + dbxmin);// # dbxmin to avoid 0/0
  Psi_LSym_np1 = alpha * Psi_LSym_n +(1.0 - alpha) * Psi_td;
  alpha = exp(-1.0*K_SSym/C_SSym*dt);
  Psi_td = (K_SSym*Psi_SApo_n - Emin_S_nph)/(K_SSym + dbxmin); // # dbxmin to avoid 0/0
  Psi_SSym_np1 = alpha * Psi_SSym_n +(1.0 - alpha) * Psi_td;
  // Rcout<< "integr. "<< Psi_td << " C_SSym:"<< C_SSym<< "  K_SSym: "<< K_SSym<<" alpha: "<<alpha << " "<< Psi_SSym_np1<<"\n";
  
  //#Step 4 : set computed values in network and update Psi_cav, PLC and Psi_AllSoil
  network.Psi_LApo = std::min(-0.00001, Psi_LApo_np1);
  network.Psi_SApo = std::min(-0.00001,Psi_SApo_np1);
  network.Psi_RCApo = std::min(-0.00001,Psi_RCApo_np1);
  network.Psi_LSym = std::min(-0.00001,Psi_LSym_np1);
  network.Psi_SSym = std::min(-0.00001,Psi_SSym_np1);
  
  //# Cavitation
  psirefL = network.Psi_LApo;  //# the reference is at current time step for other modes  (implicit, explicit)
  psirefS = network.Psi_SApo;  //# The reference is at current time step for other modes (implicit, explicit)
  if(stemCavitationRecovery!="total") {
    if (psirefS < Psi_SApo_cav) {
      network.Psi_SApo_cav = psirefS;
      network.PLC_Stem = PLC(psirefS, VCstem_slope, VCstem_P50);
    }
  } else { //Immediate refilling
    network.Psi_SApo_cav = psirefS;
    network.PLC_Stem = PLC(psirefS, VCstem_slope, VCstem_P50);
  }
  if(leafCavitationRecovery!="total") {
    if(psirefL < Psi_LApo_cav) {
      network.Psi_LApo_cav = psirefL;
      network.PLC_Leaf = PLC(psirefL, VCleaf_slope, VCleaf_P50);
    }
  } else { //Immediate refilling
    network.Psi_LApo_cav = psirefL;
    network.PLC_Leaf = PLC(psirefL, VCleaf_slope, VCleaf_P50);
  }
}

void innerSureau(List x, SureauNetwork* networks, List input, List output, int n, double tstep, 
                 bool verbose = false) {
  
  // Communication structures
  BaldocchiPhoto PB_SL, PB_SH;
  
  IntegerVector nsmalltimesteps = IntegerVector::create(6,12, 24, 60);
  // IntegerVector nsmalltimesteps = IntegerVector::create(2,4, 8, 16);
  NumericVector opt = NumericVector::create(_["Lsym"] = 1.0,
                                            _["Ssym"] = 1.0,
                                            _["Eord"] = 1.0,
                                            _["Lcav"] = 1.0,
                                            _["Scav"] = 1.0,
                                            _["CLapo"] = 1.0,
                                            _["CTapo"] = 1.0);
  
  // Extract control variables
  List control = x["control"];
  String soilFunctions = control["soilFunctions"];
  String stemCavitationRecovery = control["stemCavitationRecovery"];
  String leafCavitationRecovery = control["leafCavitationRecovery"];
  String rhizosphereOverlap = control["rhizosphereOverlap"];
  bool plantWaterPools = (rhizosphereOverlap!="total");
  bool plantCapacitance = control["plantCapacitance"];
  bool cavitationFlux = control["cavitationFlux"];
  String stomatalSubmodel = control["stomatalSubmodel"];
  bool sunlitShade = control["sunlitShade"];
  if(!cavitationFlux) {
    opt["Lcav"] = 0.0;
    opt["Scav"] = 0.0;
  }
  if(!plantCapacitance) {
    opt["CLapo"] = 0.0;
    opt["CTapo"] = 0.0;
    opt["Lsym"] = 0.0;
    opt["Ssym"] = 0.0;
  }
  
  DataFrame cohorts = Rcpp::as<Rcpp::DataFrame>(x["cohorts"]);
  int numCohorts = cohorts.nrow();
  DataFrame soil = Rcpp::as<Rcpp::DataFrame>(x["soil"]);
  int nlayers = soil.nrow();
  
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  NumericVector LAIphe = Rcpp::as<Rcpp::NumericVector>(above["LAI_expanded"]);
  NumericVector LAIlive = Rcpp::as<Rcpp::NumericVector>(above["LAI_live"]);
  
  // Extract parameters
  // Rcout<<"params\n";
  DataFrame canopyParams = Rcpp::as<Rcpp::DataFrame>(x["canopy"]);
  NumericVector zlow = canopyParams["zlow"];
  NumericVector zmid = canopyParams["zmid"];
  NumericVector zup = canopyParams["zup"];
  NumericVector Tair = canopyParams["Tair"];
  NumericVector VPair = canopyParams["VPair"];
  NumericVector Cair = canopyParams["Cair"];
  
  DataFrame paramsAnatomy = Rcpp::as<Rcpp::DataFrame>(x["paramsAnatomy"]);
  NumericVector LeafWidth = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["LeafWidth"]);
  
  DataFrame paramsTranspiration = Rcpp::as<Rcpp::DataFrame>(x["paramsTranspiration"]);
  NumericVector VCleaf_P50 = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["VCleaf_P50"]);
  NumericVector VCleaf_slope = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["VCleaf_slope"]);
  
  DataFrame paramsWaterStorage = Rcpp::as<Rcpp::DataFrame>(x["paramsWaterStorage"]);
  NumericVector StemPI0 = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["StemPI0"]);
  NumericVector StemEPS = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["StemEPS"]);
  NumericVector StemAF = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["StemAF"]);
  NumericVector Vsapwood = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["Vsapwood"]); //l·m-2 = mm
  NumericVector LeafPI0 = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["LeafPI0"]);
  NumericVector LeafEPS = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["LeafEPS"]);
  NumericVector LeafAF = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["LeafAF"]);
  NumericVector Vleaf = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["Vleaf"]); //l·m-2 = mm
  
  DataFrame internalWater = Rcpp::as<Rcpp::DataFrame>(x["internalWater"]);
  NumericVector EinstVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["Einst"]);
  NumericVector ElimVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["Elim"]);
  NumericVector Emin_LVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["Emin_L"]);
  NumericVector Emin_SVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["Emin_S"]);
  NumericVector StemPLCVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["StemPLC"]);
  NumericVector LeafPLCVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["LeafPLC"]);
  NumericVector LeafPsiVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["LeafPsi"]);
  NumericVector LeafSympPsiVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["LeafSympPsi"]);
  NumericVector StemPsiVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["StemPsi"]);
  NumericVector StemSympPsiVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["StemSympPsi"]);
  NumericVector RootCrownPsiVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["RootCrownPsi"]);
  
  // //Water pools
  DataFrame belowdf = Rcpp::as<Rcpp::DataFrame>(x["below"]);
  List belowLayers = Rcpp::as<Rcpp::List>(x["belowLayers"]);
  NumericMatrix Wpool = Rcpp::as<Rcpp::NumericMatrix>(belowLayers["Wpool"]);
  List RHOP;
  NumericVector poolProportions(numCohorts);
  if(plantWaterPools) {
    RHOP = belowLayers["RHOP"];
    poolProportions = belowdf["poolProportions"];
  }
  NumericMatrix RhizoPsiMAT = Rcpp::as<Rcpp::NumericMatrix>(belowLayers["RhizoPsi"]);
  
  //Extract output to be filled
  
  NumericMatrix SoilWaterExtract = Rcpp::as<Rcpp::NumericMatrix>(output["Extraction"]);
  List ExtractionPools = Rcpp::as<Rcpp::List>(output["ExtractionPools"]);
  NumericMatrix soilLayerExtractInst = Rcpp::as<Rcpp::NumericMatrix>(output["ExtractionInst"]);
  NumericMatrix minPsiRhizo = Rcpp::as<Rcpp::NumericMatrix>(output["RhizoPsi"]);
  List Plants = output["Plants"];
  NumericVector PWB = Plants["WaterBalance"];
  NumericVector Eplant = Plants["Transpiration"];
  NumericVector Agplant = Plants["GrossPhotosynthesis"];
  NumericVector Anplant = Plants["NetPhotosynthesis"];
  NumericVector minLeafPsi = Plants["LeafPsiMin"];
  NumericVector maxLeafPsi = Plants["LeafPsiMax"];
  NumericVector minStemPsi = Plants["StemPsi"];
  NumericVector minRootPsi = Plants["RootPsi"];
  
  List PlantsInst = output["PlantsInst"];
  NumericMatrix Einst = Rcpp::as<Rcpp::NumericMatrix>(PlantsInst["E"]);
  NumericMatrix Aginst = Rcpp::as<Rcpp::NumericMatrix>(PlantsInst["Ag"]);
  NumericMatrix Aninst = Rcpp::as<Rcpp::NumericMatrix>(PlantsInst["An"]);
  NumericMatrix dEdPInst = Rcpp::as<Rcpp::NumericMatrix>(PlantsInst["dEdP"]);
  NumericMatrix PWBinst = Rcpp::as<Rcpp::NumericMatrix>(PlantsInst["PWB"]);
  NumericMatrix StemSympRWCInst = Rcpp::as<Rcpp::NumericMatrix>(PlantsInst["StemSympRWC"]);
  NumericMatrix LeafSympRWCInst = Rcpp::as<Rcpp::NumericMatrix>(PlantsInst["LeafSympRWC"]);
  NumericMatrix StemRWCInst = Rcpp::as<Rcpp::NumericMatrix>(PlantsInst["StemRWC"]);
  NumericMatrix LeafRWCInst = Rcpp::as<Rcpp::NumericMatrix>(PlantsInst["LeafRWC"]);
  NumericMatrix StemPsiInst = Rcpp::as<Rcpp::NumericMatrix>(PlantsInst["StemPsi"]);
  NumericMatrix LeafPsiInst = Rcpp::as<Rcpp::NumericMatrix>(PlantsInst["LeafPsi"]);
  NumericMatrix RootPsiInst = Rcpp::as<Rcpp::NumericMatrix>(PlantsInst["RootPsi"]);
  NumericMatrix StemPLC = Rcpp::as<Rcpp::NumericMatrix>(PlantsInst["StemPLC"]);
  NumericMatrix LeafPLC = Rcpp::as<Rcpp::NumericMatrix>(PlantsInst["LeafPLC"]);
  NumericMatrix StemSympPsiInst = Rcpp::as<Rcpp::NumericMatrix>(PlantsInst["StemSympPsi"]);
  NumericMatrix LeafSympPsiInst = Rcpp::as<Rcpp::NumericMatrix>(PlantsInst["LeafSympPsi"]);
  List ShadeInst = output["ShadeLeavesInst"];
  NumericMatrix LAI_SH = Rcpp::as<Rcpp::NumericMatrix>(ShadeInst["LAI"]);
  NumericMatrix Vmax298_SH = Rcpp::as<Rcpp::NumericMatrix>(ShadeInst["Vmax298"]);
  NumericMatrix Jmax298_SH = Rcpp::as<Rcpp::NumericMatrix>(ShadeInst["Jmax298"]);
  NumericMatrix SWR_SH = Rcpp::as<Rcpp::NumericMatrix>(ShadeInst["Abs_SWR"]);
  NumericMatrix PAR_SH = Rcpp::as<Rcpp::NumericMatrix>(ShadeInst["Abs_PAR"]);
  NumericMatrix LWR_SH = Rcpp::as<Rcpp::NumericMatrix>(ShadeInst["Net_LWR"]);
  NumericMatrix Ag_SH = Rcpp::as<Rcpp::NumericMatrix>(ShadeInst["Ag"]);
  NumericMatrix An_SH = Rcpp::as<Rcpp::NumericMatrix>(ShadeInst["An"]);
  NumericMatrix E_SH = Rcpp::as<Rcpp::NumericMatrix>(ShadeInst["E"]);
  NumericMatrix VPD_SH = Rcpp::as<Rcpp::NumericMatrix>(ShadeInst["VPD"]);
  NumericMatrix Psi_SH = Rcpp::as<Rcpp::NumericMatrix>(ShadeInst["Psi"]);
  NumericMatrix Temp_SH = Rcpp::as<Rcpp::NumericMatrix>(ShadeInst["Temp"]);
  NumericMatrix GSW_SH = Rcpp::as<Rcpp::NumericMatrix>(ShadeInst["Gsw"]);
  NumericMatrix Ci_SH = Rcpp::as<Rcpp::NumericMatrix>(ShadeInst["Ci"]);
  List SunlitInst = output["SunlitLeavesInst"];
  NumericMatrix LAI_SL = Rcpp::as<Rcpp::NumericMatrix>(SunlitInst["LAI"]);
  NumericMatrix Vmax298_SL = Rcpp::as<Rcpp::NumericMatrix>(SunlitInst["Vmax298"]);
  NumericMatrix Jmax298_SL = Rcpp::as<Rcpp::NumericMatrix>(SunlitInst["Jmax298"]);
  NumericMatrix SWR_SL = Rcpp::as<Rcpp::NumericMatrix>(SunlitInst["Abs_SWR"]);
  NumericMatrix PAR_SL = Rcpp::as<Rcpp::NumericMatrix>(SunlitInst["Abs_PAR"]);
  NumericMatrix LWR_SL = Rcpp::as<Rcpp::NumericMatrix>(SunlitInst["Net_LWR"]);
  NumericMatrix Ag_SL = Rcpp::as<Rcpp::NumericMatrix>(SunlitInst["Ag"]);
  NumericMatrix An_SL = Rcpp::as<Rcpp::NumericMatrix>(SunlitInst["An"]);
  NumericMatrix E_SL = Rcpp::as<Rcpp::NumericMatrix>(SunlitInst["E"]);
  NumericMatrix VPD_SL = Rcpp::as<Rcpp::NumericMatrix>(SunlitInst["VPD"]);
  NumericMatrix Psi_SL = Rcpp::as<Rcpp::NumericMatrix>(SunlitInst["Psi"]);
  NumericMatrix Temp_SL = Rcpp::as<Rcpp::NumericMatrix>(SunlitInst["Temp"]);
  NumericMatrix GSW_SL = Rcpp::as<Rcpp::NumericMatrix>(SunlitInst["Gsw"]);
  NumericMatrix Ci_SL = Rcpp::as<Rcpp::NumericMatrix>(SunlitInst["Ci"]);
  
  //Extract input  
  // Rcout<<"input\n";
  NumericVector zWind = input["zWind"];
  double Patm = input["Patm"];
  double f_dry = input["f_dry"];
  
  IntegerVector iLayerCohort = input["iLayerCohort"];
  IntegerVector iLayerSunlit = input["iLayerSunlit"];
  IntegerVector iLayerShade = input["iLayerShade"];
  IntegerVector nlayerscon = input["nlayerscon"];
  LogicalMatrix layerConnected = input["layerConnected"];
  List layerConnectedPools = input["layerConnectedPools"];
  

  for(int c=0;c<numCohorts;c++) { //Plant cohort loop
    
    if(LAIphe[c]>0.0 && LeafPLCVEC[c] < 0.999) {
      // Rcout << "\n*** HOUR STEP " << n << " cohort " << c << "***\n";
      //Init temporary variables
      SureauNetwork network_n;
      network_n.params.k_RCApoInit = new double[networks[c].params.npools];
      network_n.k_RSApo = new double[networks[c].params.npools];
      network_n.k_SoilToStem = new double[networks[c].params.npools];
      network_n.k_Soil = new double[networks[c].params.npools];
      network_n.PsiSoil = new double[networks[c].params.npools];
      NumericVector ElayersVEC(networks[c].params.npools,0.0); //Instantaneous flow rate
      NumericVector fluxSoilToStem_mm(networks[c].params.npools, 0.0); //Cummulative flow
      
      // # A. LOOP ON THE IMPLICIT SOLVER IN PSI, trying different time steps until results are OK
      bool regulationWellComputed = false;
      bool cavitationWellComputed = false;
      
      double gmin_S = networks[c].params.gmin_S;
      double gmin20 = networks[c].params.gmin20;
      double TPhase_gmin = networks[c].params.TPhase_gmin;
      double Q10_1_gmin = networks[c].params.Q10_1_gmin;
      double Q10_2_gmin = networks[c].params.Q10_2_gmin;
      double fTRBToLeaf = networks[c].params.fTRBToLeaf;
      double Gsw_AC_slope = networks[c].params.Gsw_AC_slope;
      double gsNight = networks[c].params.gsNight;
      double slope_gs = networks[c].params.slope_gs;
      double P50_gs = networks[c].params.P50_gs;
      double LAI = networks[c].LAI;
      
      int nwhilecomp = 0;
      
      double Agsum = 0.0, Ansum = 0.0;
      
      while ((!regulationWellComputed || !cavitationWellComputed) && (nwhilecomp<nsmalltimesteps.size())) { //# LOOP TO TRY DIFFERENT TIME STEPS
        //Copy values to temporary network
        copyNetwork(networks[c], network_n); 
        
        regulationWellComputed = false;
        cavitationWellComputed = false;
        double deltaRegulMax = 1.0e-100;
        double deltaPLCMax = 1.0e-100;
        
        //Reset output fluxes to zero
        Agsum = 0.0;
        Ansum = 0.0;
        EinstVEC[c] = 0.0;
        ElimVEC[c] = 0.0;
        Emin_LVEC[c] = 0.0;
        Emin_SVEC[c] = 0.0;
        for(int i=0;i < networks[c].params.npools;i++) {
          ElayersVEC[i] = 0.0;
          fluxSoilToStem_mm[i] = 0.0; 
        }
        
        int nts = nsmalltimesteps[nwhilecomp];// # number of small time steps
        double dt = tstep / ((double) nts); //Determine number of seconds of small time steps
        // Rcout<< " Attempt #" << nwhilecomp<<" nts "<< nts << " dt " << dt << "\n";
        for(int its = 1; its <= nts; its++) { //#INTERNAL LOOP ON SMALL TIME STEPS
          
          //Current leaf water potential (same for sunlit and shade leaves)
          double Psi_LSym = network_n.Psi_LSym;
          // Current stomatal regulation ("Sigmoid")
          double regul_ini = 1.0 - (1.0 / (1.0 + exp(slope_gs / 25.0 * (Psi_LSym - P50_gs))));
          
          //Leaf temperature for sunlit and shade leaves
          double Elim_SL = network_n.Elim_SL;
          double Elim_SH = network_n.Elim_SH;
          double Elim = network_n.Elim;
          if(NumericVector::is_na(Elim_SL)) Elim_SL = Elim * (LAI_SL(c,n)/LAI);
          if(NumericVector::is_na(Elim_SH)) Elim_SH = Elim * (LAI_SH(c,n)/LAI);
          if(!sunlitShade) Elim_SH = Elim_SL;
          
          Temp_SL(c,n) = leafTemperature2(SWR_SL(c,n)/LAI_SL(c,n), LWR_SL(c,n)/LAI_SL(c,n), 
                  Tair[iLayerSunlit[c]], zWind[iLayerSunlit[c]], 
                                              Elim_SL,  LeafWidth[c]);
          Temp_SH(c,n) = leafTemperature2(SWR_SH(c,n)/LAI_SH(c,n), LWR_SH(c,n)/LAI_SH(c,n), 
                  Tair[iLayerShade[c]], zWind[iLayerShade[c]], 
                                             Elim_SH,  LeafWidth[c]);
          if(!sunlitShade) Temp_SH(c,n) = Temp_SL(c,n);
          
          //VPD
          double VPD_air = meteoland::utils_saturationVP(Tair[iLayerCohort[c]]) - VPair[iLayerCohort[c]];
          VPD_SL(c,n) = std::max(0.0,leafVapourPressure(Temp_SL(c,n), Psi_LSym) - VPair[iLayerSunlit[c]]);
          VPD_SH(c,n) = std::max(0.0,leafVapourPressure(Temp_SH(c,n), Psi_LSym) - VPair[iLayerShade[c]]);
          if(!sunlitShade) VPD_SH(c,n) = VPD_SL(c,n);
          // Rcout<< "  AirT "<< Tair[iLayerCohort[c]] << " LT_SL "<< Temp_SL(c,n)<< " LT_SH "<< Temp_SH(c,n)<<"\n";
          // Rcout<< "  VPD_air "<< VPD_air << " VPD_SL "<< VPD_SL(c,n)<< " VPD_SH "<< VPD_SH(c,n)<<"\n";
          
          //gCR = g Crown
          double gCR = 1000.0*gCrown(zWind[iLayerCohort[c]]); 
          //Assumes well coupled canopy (for compatibility with Sperry and leaf temperature balance)
          //gBL = g Boundary Layer
          double gBL = 1000.0*gLeafBoundary(zWind[iLayerCohort[c]], LeafWidth[c]); // mmol boundary layer conductance
          
          //# Leaf cuticular conductances and cuticular transpiration
          double gmin_SL = gmin(Temp_SL(c,n), gmin20, TPhase_gmin, Q10_1_gmin, Q10_2_gmin);
          double gmin_SH = gmin(Temp_SH(c,n), gmin20, TPhase_gmin, Q10_1_gmin, Q10_2_gmin);
          double Emin_L_SL = Emin(gmin_SL, gBL, gCR, VPD_SL(c,n), Patm)*f_dry; //Add f_dry to decrease transpiration in rainy days
          double Emin_L_SH = Emin(gmin_SH, gBL, gCR, VPD_SH(c,n), Patm)*f_dry;
          double Emin_L = ((Emin_L_SL*LAI_SL(c,n)) + (Emin_L_SH*LAI_SH(c,n)))/LAI; 
          network_n.Emin_L = Emin_L;
          
          //Compute stem cuticular transpiration
          double Emin_S = fTRBToLeaf * Emin(gmin_S, gBL, gCR, VPD_air, Patm);
          network_n.Emin_S =  Emin_S*f_dry; //Add f_dry to decrease transpiration in rainy days
          // Rcout<< "  Emin_S "<< Emin_S<<" Emin_L_SL "<< Emin_L_SL<<" Emin_L_SH "<< Emin_L_SH<<" Emin_L "<< Emin_L<<"\n";
          
          // Current stomatal regulation ("Sigmoid")
          double regul = 1.0 - (1.0 / (1.0 + exp(slope_gs / 25.0 * (Psi_LSym - P50_gs))));
          
          double gs_SL, gs_SH;
          if(stomatalSubmodel=="Jarvis") {
            gs_SL = gsJarvis(networks[c].params, PAR_SL(c,n), Temp_SL(c,n));
            gs_SH = gsJarvis(networks[c].params, PAR_SH(c,n), Temp_SH(c,n));
            //Rcout<< "  PAR_SL "<< PAR_SL(c,n)<<"  gs_SL "<< gs_SL<<"  PAR_SH "<< PAR_SH(c,n)<<" gs_SH "<< gs_SH<<"\n";
            gs_SL = gs_SL * regul;
            gs_SH = gs_SH * regul;
          } else {
            photosynthesisBaldocchi_inner(PB_SL, 
                                          irradianceToPhotonFlux(PAR_SL(c,n))/LAI_SL(c,n), 
                                          Cair[iLayerSunlit[c]], 
                                          std::max(0.0,Temp_SL(c,n)), 
                                          zWind[iLayerCohort[c]],
                                          Vmax298_SL(c,n)/LAI_SL(c,n), 
                                          Jmax298_SL(c,n)/LAI_SL(c,n), 
                                          LeafWidth[c],
                                          Gsw_AC_slope,
                                          gsNight/1000.0);
            gs_SL = PB_SL.Gsw*1000.0; //From mmol to mol 
            gs_SL = std::max(gsNight, gs_SL)*regul;
            // Rcout<<c << " "<<n << " Bald gs: "<< PB_SL.Gsw << " regul: "<< regul << " gs_SL: "<< gs_SL<<"\n";
            photosynthesisBaldocchi_inner(PB_SH, 
                                          irradianceToPhotonFlux(PAR_SH(c,n))/LAI_SH(c,n), 
                                          Cair[iLayerSunlit[c]], 
                                          std::max(0.0,Temp_SH(c,n)), 
                                          zWind[iLayerCohort[c]],
                                          Vmax298_SH(c,n)/LAI_SH(c,n), 
                                          Jmax298_SH(c,n)/LAI_SH(c,n), 
                                          LeafWidth[c],
                                          Gsw_AC_slope,
                                          gsNight/1000.0);
            gs_SH = PB_SH.Gsw*1000.0; //From mmol to mol
            gs_SH = std::max(gsNight, gs_SH)*regul;
          }
          if(!sunlitShade) gs_SH = gs_SL;
          
          // Store stomatal conductance          
          GSW_SL(c,n) = gs_SL/1000.0; // From mmol to mol
          GSW_SH(c,n) = gs_SH/1000.0; // From mmol to mol
          // Stomatal transpiration
          double Gwdiff_SL = 1.0/(1.0/gCR + 1.0/gs_SL + 1.0/gBL); 
          double Gwdiff_SH = 1.0/(1.0/gCR + 1.0/gs_SH + 1.0/gBL); 
          Elim_SL = Gwdiff_SL * (VPD_SL(c,n)/Patm)*f_dry; //Add f_dry to decrease transpiration in rainy days
          Elim_SH = Gwdiff_SH * (VPD_SH(c,n)/Patm)*f_dry;
          
          //Photosynthesis
          double Gwdiff_all_SL = 1.0/(1.0/gCR + 1.0/(gs_SL + gmin_SL) + 1.0/gBL); 
          double Gwdiff_all_SH = 1.0/(1.0/gCR + 1.0/(gs_SH + gmin_SH) + 1.0/gBL); 
          NumericVector LP_SL = leafphotosynthesis(irradianceToPhotonFlux(PAR_SL(c,n))/LAI_SL(c,n), 
                                                   Cair[iLayerSunlit[c]], Gwdiff_all_SL/(1000.0*1.6), //From mmol to mol 
                                                   std::max(0.0,Temp_SL(c,n)), 
                                                   Vmax298_SL(c,n)/LAI_SL(c,n), Jmax298_SL(c,n)/LAI_SL(c,n));
          NumericVector LP_SH = leafphotosynthesis(irradianceToPhotonFlux(PAR_SH(c,n))/LAI_SH(c,n), 
                                                   Cair[iLayerShade[c]], Gwdiff_all_SH/(1000.0*1.6), //From mmol to mol
                                                   std::max(0.0,Temp_SH(c,n)), 
                                                   Vmax298_SH(c,n)/LAI_SH(c,n), Jmax298_SH(c,n)/LAI_SH(c,n));
          if(!sunlitShade) LP_SH = LP_SL;
          Ci_SL(c,n) = LP_SL[0];
          Ci_SH(c,n) = LP_SH[0];
          Ag_SL(c,n) = LP_SL[1];
          Ag_SH(c,n) = LP_SH[1];
          An_SL(c,n) = Ag_SL(c,n) - 0.015*VmaxTemp(Vmax298_SL(c,n)/LAI_SL(c,n), Temp_SL(c,n));
          An_SH(c,n) = Ag_SH(c,n) - 0.015*VmaxTemp(Vmax298_SH(c,n)/LAI_SH(c,n), Temp_SH(c,n));
          
          Agsum += Ag_SL(c,n)*LAI_SL(c,n) + Ag_SH(c,n)*LAI_SH(c,n);
          Ansum += An_SL(c,n)*LAI_SL(c,n) + An_SH(c,n)*LAI_SH(c,n);
          
          network_n.Elim_SL = Elim_SL;
          network_n.Elim_SH = Elim_SH;
          Elim = ((Elim_SL*LAI_SL(c,n)) + (Elim_SH*LAI_SH(c,n)))/LAI; 
          network_n.Elim = Elim;
          // Rcout<< "  Elim_SL "<< Elim_SL<<"  Elim_SH "<< Elim_SH<<"  Elim "<< Elim<<"\n";
          
          //Add transpiration sources
          network_n.Einst = Elim + Emin_S + Emin_L;
          network_n.Einst_SL = Elim_SL + Emin_L_SL; //For sunlit photosynthesis/transpiration
          network_n.Einst_SH = Elim_SH + Emin_L_SH; //For shade photosynthesis/transpiration
          
          //Effects on water potentials and flows
          semi_implicit_integration_inner(network_n, 
                                          dt, opt, stemCavitationRecovery, leafCavitationRecovery);
          update_conductances(network_n);
          update_capacitances(network_n);
          
          // # QUANTITIES TO CHECK IF THE RESOLUTION IS OK
          // # 1. delta regulation between n and np1 (MIQUEL: Only Psi_LSym changes between the two calculations, params should be the same)
          deltaRegulMax = std::max(deltaRegulMax,std::abs(regul - regul_ini));
          
          // # 2. PLC at n and np1
          deltaPLCMax = std::max(deltaPLCMax, network_n.PLC_Leaf - network_n.PLC_Leaf);
          deltaPLCMax = std::max(deltaPLCMax, network_n.PLC_Stem - network_n.PLC_Stem);
          
          // # 3. update of soil on small time step (done by FP in version 16)
          double Psi_SApo = network_n.Psi_SApo;
          double* k_SoilToStem = network_n.k_SoilToStem; 
          double* PsiSoil = network_n.PsiSoil;
          for(int l=0;l < networks[0].params.npools;l++) {
            double fluxSoilToStem_mmolm2s = k_SoilToStem[l]*(PsiSoil[l] - Psi_SApo);
            ElayersVEC[l] += fluxSoilToStem_mmolm2s;
            fluxSoilToStem_mm[l] += (fluxSoilToStem_mmolm2s*0.001*0.01802*LAIphe[c]*dt);
          }
          //MIQUEL (27/04/2024): Changed network to network_n
          EinstVEC[c] += network_n.Einst;
          ElimVEC[c] += network_n.Elim;
          Emin_LVEC[c] += network_n.Emin_L;
          Emin_SVEC[c] += network_n.Emin_S;
          
        } //# end loop small time step
        
        //Divide average fluxes by time steps
        for(int l=0;l < networks[c].params.npools;l++) ElayersVEC[l] = ElayersVEC[l]/((double) nts);
        EinstVEC[c] = EinstVEC[c]/((double) nts);
        ElimVEC[c] = ElimVEC[c]/((double) nts);
        Emin_LVEC[c] = Emin_LVEC[c]/((double) nts);
        Emin_SVEC[c] = Emin_SVEC[c]/((double) nts);
        Agsum = Agsum/((double) nts);
        Ansum = Ansum/((double) nts);
        
        // # TESTS ON RESOLUTION
        network_n.Diag_deltaRegulMax = deltaRegulMax;
        regulationWellComputed = (deltaRegulMax<0.05);
        network_n.Diag_deltaPLCMax = deltaPLCMax;
        cavitationWellComputed = (deltaPLCMax<1.0);// # 1%
        network_n.Diag_timeStepInSeconds = dt;
        nwhilecomp = nwhilecomp + 1;
      } //# end while

      // # B. SAVING SOLUTION AT NEXT TIME STEP IN ORIGINAL NETWORK
      copyNetwork(network_n, networks[c]);
      
      //Store leaf values (final substep)
      E_SL(c,n) = networks[c].Einst_SL;
      E_SH(c,n) = networks[c].Einst_SH;
      Psi_SH(c,n) = networks[c].Psi_LSym;
      Psi_SL(c,n) = networks[c].Psi_LSym;
      dEdPInst(c,n) = networks[c].k_Plant;
      
      //Store state
      LeafPsiVEC[c] = networks[c].Psi_LApo;
      LeafSympPsiVEC[c] = networks[c].Psi_LSym;
      StemPsiVEC[c] = networks[c].Psi_SApo;
      StemSympPsiVEC[c] = networks[c].Psi_SSym;
      RootCrownPsiVEC[c] = networks[c].Psi_RCApo;
      StemPLCVEC[c] = networks[c].PLC_Stem/100.0;
      LeafPLCVEC[c] = networks[c].PLC_Leaf/100.0;
      
      // Rcout<<iPMSunlit[c]<<" "<<iPMShade[c] <<" "<<GwSunlit[iPMSunlit[c]]<<" "<<GwShade[iPMShade[c]]<<" "<<fittedE[iPMSunlit[c]]<<" "<<fittedE[iPMShade[c]]<<"\n";
      //Get leaf status
      
      //Scale photosynthesis
      Aginst(c,n) = (1e-6)*12.01017*Agsum*tstep;
      Aninst(c,n) = (1e-6)*12.01017*Ansum*tstep;
      
      //Scale from instantaneous flow to water volume in the time step
      Einst(c,n) = EinstVEC[c]*0.001*0.01802*LAIphe[c]*tstep;
      
      
      //Calculate and copy RhizoPsi from connected layers to RhizoPsi from soil layers
      calculateRhizoPsi_inner(c, nlayers,
                              networks[c], RhizoPsiMAT,
                              layerConnected,
                              RHOP, layerConnectedPools,
                              plantWaterPools);
      
      //Balance between extraction and transpiration
      PWBinst(c,n) = sum(fluxSoilToStem_mm) - Einst(c,n);
      
      //Add step transpiration to daily plant cohort transpiration
      Eplant[c] += Einst(c,n);
      Anplant[c] += Aninst(c,n);
      Agplant[c] += Aginst(c,n);
      //Add PWB
      PWB[c] += PWBinst(c,n);
      
      
      //Copy transpiration and from connected layers to transpiration from soil layers
      //And update soil water content (soil water potential will not be updated until next day!)
      if(!plantWaterPools) {
        int cl = 0;
        for(int l=0;l<nlayers;l++) {
          if(layerConnected(c,l)) {
            SoilWaterExtract(c,l) += fluxSoilToStem_mm[cl]; //Add to cummulative transpiration from layers
            soilLayerExtractInst(l,n) += fluxSoilToStem_mm[cl];
            cl++;
          }
        }
      } else {
        NumericMatrix RHOPcoh = Rcpp::as<Rcpp::NumericMatrix>(RHOP[c]);
        LogicalMatrix layerConnectedCoh = Rcpp::as<Rcpp::LogicalMatrix>(layerConnectedPools[c]);
        int cl = 0;
        for(int j = 0;j<numCohorts;j++) {
          NumericMatrix ExtractionPoolsCoh = Rcpp::as<Rcpp::NumericMatrix>(ExtractionPools[j]);
          for(int l=0;l<nlayers;l++) {
            if(layerConnectedCoh(j,l)) {
              SoilWaterExtract(c,l) += fluxSoilToStem_mm[cl]; //Add to cummulative transpiration from layers
              soilLayerExtractInst(l,n) += fluxSoilToStem_mm[cl];
              ExtractionPoolsCoh(c,l) += fluxSoilToStem_mm[cl];
              cl++;
            }
          }
        }
      }
      //Delete pointers
      deleteSureauNetworkPointers(network_n);
    } else if(LAIlive[c]>0.0) { //Cohorts with living individuals but no LAI (or completely embolized)
      E_SL(c,n) = 0.0;
      E_SH(c,n) = 0.0;
      Psi_SH(c,n) = networks[c].Psi_LSym;
      Psi_SL(c,n) = networks[c].Psi_LSym;
      dEdPInst(c,n) = networks[c].k_Plant;
      LeafPsiVEC[c] = networks[c].Psi_LApo;
      LeafSympPsiVEC[c] = networks[c].Psi_LSym;
      StemPsiVEC[c] = networks[c].Psi_SApo;
      StemSympPsiVEC[c] = networks[c].Psi_SSym;
      RootCrownPsiVEC[c] = networks[c].Psi_RCApo;
      StemPLCVEC[c] = networks[c].PLC_Stem/100.0;
      LeafPLCVEC[c] = networks[c].PLC_Leaf/100.0;
      Aginst(c,n) = 0.0;
      Aninst(c,n) = 0.0;
      Einst(c,n) = 0.0;
      PWBinst(c,n) = 0.0;
    }
  }
  
}

//' @rdname sureau_ecos
//' @param network A hydraulic network element of the list returned by \code{initSureauNetworks}
//' @param dt Smallest time step (seconds)
//' @param opt Option flag vector
//' @param stemCavitationRecovery,leafCavitationRecovery A string indicating how refilling of embolized conduits is done:
//'           \itemize{
//'             \item{"none" - no refilling.}
//'             \item{"annual" - every first day of the year.}
//'             \item{"rate" - following a rate of new sapwood formation.}
//'             \item{"total" - instantaneous complete refilling.}
//'           }
//' @keywords internal
// [[Rcpp::export("semi_implicit_integration")]]
List semi_implicit_integration(List network, double dt, NumericVector opt, 
                               String stemCavitationRecovery = "annual", String leafCavitationRecovery = "total") {

  //Copy values from List to SureauNetwork
  SureauNetwork snetwork = listToStruct(network);
  //Integration
  semi_implicit_integration_inner(snetwork,
                                  dt, opt, 
                                  stemCavitationRecovery, 
                                  leafCavitationRecovery);
  //Copy back values from SureauNetwork to List
  network = structToList(snetwork);
  //Free memory
  deleteSureauNetworkPointers(snetwork);
  return(network);
}

