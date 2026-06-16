#define STRICT_R_HEADERS
#include <RcppArmadillo.h>
#include "inner_sureau_c.h"
#include "photosynthesis.h"
#include "photosynthesis_c.h"
#include "biophysicsutils_c.h"
#include "meteoland/utils_c.hpp"
#include "hydraulics.h"
#include "hydraulics_c.h"
#include "soil.h"
#include "soil_c.h"
#include "tissuemoisture_c.h"
#include <meteoland.h>
using namespace Rcpp;


List structToList(SureauNetwork& snetwork) {
  
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
  NumericVector VGrhizo_kmax(snetwork.params.npools, NA_REAL);
  NumericVector k_RCApoInit(snetwork.params.npools, NA_REAL);
  for(int l=0;l < snetwork.params.npools; l++) {
    k_RCApoInit[l] = snetwork.params.k_RCApoInit[l];
    VGrhizo_kmax[l] = snetwork.params.VGrhizo_kmax[l];
  }
  params.push_back(k_RCApoInit, "k_RCApoInit"); 
  params.push_back(VGrhizo_kmax, "VGrhizo_kmax"); 
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
  network.push_back(snetwork.LAImistletoe, "LAImistletoe");
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
  network.push_back(snetwork.Emist, "Emist"); 
  network.push_back(snetwork.Emist_SL, "Emist_SL");
  network.push_back(snetwork.Emist_SH, "Emist_SH"); 
  
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
  NumericVector VGrhizo_kmax = Rcpp::as<Rcpp::NumericVector>(params["VGrhizo_kmax"]);
  snetwork.params.k_RCApoInit = new double[k_RCApoInit.size()];
  snetwork.params.VGrhizo_kmax = new double[VGrhizo_kmax.size()];
  for(int l=0;l < k_RCApoInit.size(); l++) {
    snetwork.params.k_RCApoInit[l] = k_RCApoInit[l];
    snetwork.params.VGrhizo_kmax[l] = VGrhizo_kmax[l];
  }
  
  snetwork.LAI = network["LAI"];
  snetwork.LAImistletoe = network["LAImistletoe"];
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
  snetwork.Emist = network["Emist"];
  snetwork.Emist_SL = network["Emist_SL"];
  snetwork.Emist_SH = network["Emist_SH"];
  
  //Diagnostics
  snetwork.Diag_nwhile_cavit = network["Diag_nwhile_cavit"];
  snetwork.Diag_deltaRegulMax = network["Diag_deltaRegulMax"];
  snetwork.Diag_deltaPLCMax = network["Diag_deltaPLCMax"];
  snetwork.Diag_timeStepInSeconds = network["Diag_timeStepInSeconds"];

  return(snetwork);
}


void initSureauParams_inner(SureauParams &params, int c,
                            DataFrame internalWater, 
                            DataFrame paramsTranspiration, DataFrame paramsWaterStorage,
                            NumericVector VCroot_kmax, NumericVector VGrhizo_kmax,
                            List control, double sapFluidityDay = 1.0) {
  String stomatalSubmodel = control["stomatalSubmodel"];

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
  params.VGrhizo_kmax = new double [params.npools];
  for(int l = 0;l<params.npools;l++)  {
    params.k_RCApoInit[l] = sapFluidityDay*VCroot_kmax[l];
    params.VGrhizo_kmax[l] = VGrhizo_kmax[l];
  }
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
  if(!leafCuticularTranspiration) gmin20 = 0.0;
  params.gmin20 = gmin20; //mol to mmol 
  params.gsMax = Gswmax[c]*1000.0; //mol to mmol 
  bool stemCuticularTranspiration = control["stemCuticularTranspiration"];
  double gmin_S = Gswmin[c]*1000.0;  // gmin for stem equal to leaf gmin, in mmol
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
void initSureauNetwork_inner(SureauNetwork &network, int c, NumericVector LAIphe, NumericVector LAImistletoe,
                             DataFrame internalWater, 
                             DataFrame paramsAnatomy, DataFrame paramsTranspiration, DataFrame paramsWaterStorage,
                             NumericVector VCroot_kmax, NumericVector VGrhizo_kmax,
                             NumericVector PsiSoil, NumericVector VG_n, NumericVector VG_alpha,
                             List control, double sapFluidityDay = 1.0) {
  
  String stomatalSubmodel = control["stomatalSubmodel"];
  
  //Root distribution input
  NumericVector Einst = Rcpp::as<Rcpp::NumericVector>(internalWater["Einst"]);
  NumericVector Emist = Rcpp::as<Rcpp::NumericVector>(internalWater["Emist"]);
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
  

  //Params
  initSureauParams_inner(network.params, c, internalWater, 
                         paramsTranspiration, paramsWaterStorage,
                         VCroot_kmax, VGrhizo_kmax,
                         control, sapFluidityDay);
  
  //LAI
  network.LAI = LAIphe[c];
  network.LAImistletoe = LAImistletoe[c];
  
  //Water potentials
  network.Psi_LApo = LeafPsiVEC[c]; 
  network.Psi_LSym = LeafSympPsiVEC[c];
  network.Psi_RCApo = RootCrownPsiVEC[c];

  network.Psi_SApo = StemPsiVEC[c]; 
  network.Psi_SSym = StemSympPsiVEC[c];
  network.Psi_SApo_cav = std::min(0.0, invPLC_c(StemPLCVEC[c]*100.0, VCstem_slope[c], VCstem_P50[c])); //Sureau operates with %
  network.Psi_LApo_cav = std::min(0.0, invPLC_c(LeafPLCVEC[c]*100.0, VCleaf_slope[c], VCleaf_P50[c])); //Sureau operates with %
  
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
    network.k_Soil[l] = vanGenuchtenConductance_c(PsiSoil[l],
                                                VGrhizo_kmax[l], 
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
  network.Emist = Emist[c]; //Mistletoe transpiration
  network.Emist_SL = NA_REAL; //Mistletoe transpiration (sunlit leaves)
  network.Emist_SH = NA_REAL; //Mistletoe transpiration (shade leaves)
  
  //Diagnostics
  network.Diag_nwhile_cavit = NA_INTEGER;
  network.Diag_deltaRegulMax = NA_REAL;
  network.Diag_deltaPLCMax = NA_REAL;
  network.Diag_timeStepInSeconds = NA_REAL;
  
  // Update plant conductances and capacitances according to network status
  update_conductances_c(network);
  update_capacitances_c(network);
}


List initSureauNetwork(int c, NumericVector LAIphe, NumericVector LAImistletoe,
                       DataFrame internalWater, 
                       DataFrame paramsAnatomy, DataFrame paramsTranspiration, DataFrame paramsWaterStorage,
                       NumericVector VCroot_kmax, NumericVector VGrhizo_kmax,
                       NumericVector PsiSoil, NumericVector VG_n, NumericVector VG_alpha,
                       List control, double sapFluidityDay = 1.0) {
  
  SureauNetwork snetwork;
  initSureauNetwork_inner(snetwork, c, LAIphe, LAImistletoe,
                          internalWater, 
                          paramsAnatomy, paramsTranspiration, paramsWaterStorage,
                          VCroot_kmax, VGrhizo_kmax,
                          PsiSoil, VG_n, VG_alpha,
                          control, sapFluidityDay);
  
  List network = structToList(snetwork);

  deleteSureauNetworkPointers_c(snetwork);
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
  NumericVector LAImistletoe = Rcpp::as<Rcpp::NumericVector>(above["LAI_mistletoe"]);
  
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
    networks[c] = initSureauNetwork(c, LAIphe, LAImistletoe,
                                     internalWater, 
                                     paramsAnatomy, paramsTranspiration, paramsWaterStorage,
                                     VCroot_kmax(c,_), VGrhizo_kmax(c,_),
                                     psiSoil, VG_n, VG_alpha,
                                     control, 1.0);
  }
  networks.attr("names") = above.attr("row.names");
  return(networks);
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
  
  SureauOpt opt_c;
  opt_c.Lsym = opt["Lsym"];
  opt_c.Ssym = opt["Ssym"];
  opt_c.Lcav = opt["Lcav"];
  opt_c.Scav = opt["Scav"];
  opt_c.CLapo = opt["CLapo"];
  opt_c.CTapo = opt["CTapo"];
  
  //Integration
  std::string stemCavitationRecovery_str = stemCavitationRecovery.get_cstring();
  std::string leafCavitationRecovery_str = leafCavitationRecovery.get_cstring();
  
  semi_implicit_integration_inner_c(snetwork,
                                    dt, opt_c, 
                                    stemCavitationRecovery_str, 
                                    leafCavitationRecovery_str);
  
  //Copy back values from SureauNetwork to List
  network = structToList(snetwork);
  //Free memory
  deleteSureauNetworkPointers_c(snetwork);
  return(network);
}

