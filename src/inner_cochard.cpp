#define STRICT_R_HEADERS
#include <Rcpp.h>
#include "photosynthesis.h"
#include "biophysicsutils.h"
#include "hydraulics.h"
#include "soil.h"
#include "tissuemoisture.h"
using namespace Rcpp;

double PLC_derivative(double plc, double slope) {
  return(-1.0*slope/25.0 * plc/100 * (1.0 - plc/100));
}
double PLC(double Pmin, double slope, double P50){
  return (100.0 / (1.0 + exp(slope / 25.0 * (Pmin - P50))));
}
double invPLC(double plc, double slope, double P50){
  return(P50 + log(100.0/plc)*(25.0/slope));
}
// [[Rcpp::export]]
List initHydraulicArchitecture(List x) {
  //Root distribution input
  DataFrame belowdf = Rcpp::as<Rcpp::DataFrame>(x["below"]);
  List belowLayers = Rcpp::as<Rcpp::List>(x["belowLayers"]);
  NumericMatrix VCroot_kmax= Rcpp::as<Rcpp::NumericMatrix>(belowLayers["VCroot_kmax"]);
  NumericMatrix VGrhizo_kmax= Rcpp::as<Rcpp::NumericMatrix>(belowLayers["VGrhizo_kmax"]);
  
  DataFrame internalWater = Rcpp::as<Rcpp::DataFrame>(x["internalWater"]);
  NumericVector Einst = Rcpp::as<Rcpp::NumericVector>(internalWater["Einst"]);
  NumericVector Elim = Rcpp::as<Rcpp::NumericVector>(internalWater["Elim"]);
  NumericVector Emin_L = Rcpp::as<Rcpp::NumericVector>(internalWater["Emin_L"]);
  NumericVector Emin_S = Rcpp::as<Rcpp::NumericVector>(internalWater["Emin_S"]);
  NumericVector StemPLCVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["StemPLC"]);
  NumericVector LeafPLCVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["LeafPLC"]);
  NumericVector StemPsiVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["StemPsi"]);
  NumericVector LeafPsiVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["LeafPsi"]);
  NumericVector RootCrownPsiVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["RootCrownPsi"]);
  NumericVector StemSympPsiVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["StemSympPsi"]);

  DataFrame paramsTranspiration = Rcpp::as<Rcpp::DataFrame>(x["paramsTranspiration"]);
  NumericVector Plant_kmax = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["Plant_kmax"]);
  NumericVector VCleaf_kmax = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["VCleaf_kmax"]);
  NumericVector VCleaf_P50 = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["VCleaf_P50"]);
  NumericVector VCleaf_slope = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["VCleaf_slope"]);
  NumericVector VCstem_kmax = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["VCstem_kmax"]);
  NumericVector VCstem_P50 = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["VCstem_P50"]);
  NumericVector VCstem_slope = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["VCstem_slope"]);
  NumericVector VCroot_P50 = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["VCroot_P50"]);
  NumericVector VCroot_slope = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["VCroot_slope"]);
  
  DataFrame paramsWaterStorage = Rcpp::as<Rcpp::DataFrame>(x["paramsWaterStorage"]);
  NumericVector Vleaf = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["Vleaf"]);
  NumericVector LeafAF = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["LeafAF"]);
  NumericVector Vsapwood = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["Vsapwood"]);
  NumericVector StemAF = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["StemAF"]);
  
  int numCohorts = StemPLCVEC.size();
  List HAs(numCohorts);
  for(int c = 0;c<numCohorts;c++) {
    List HA = List::create();
    //Params
    List params = List::create();
    params.push_back(VCleaf_P50[c], "VCleaf_P50"); 
    params.push_back(VCleaf_slope[c], "VCleaf_slope"); 
    params.push_back(VCstem_P50[c], "VCstem_P50"); 
    params.push_back(VCstem_slope[c], "VCstem_slope"); 
    params.push_back(VCroot_P50[c], "VCrootP50"); 
    params.push_back(VCroot_slope[c], "VCroot_slope"); 
    
    HA.push_back(params, "params");
    
    //Water potentials
    HA.push_back(StemPsiVEC[c], "Psi_LApo"); // Leaf apo psi  = Stem psi
    HA.push_back(LeafPsiVEC[c], "Psi_LSym"); // Leaf sym psi = Leaf psi
    HA.push_back(RootCrownPsiVEC[c], "Psi_SApo"); // Stem apo psi = root crown psi
    HA.push_back(StemSympPsiVEC[c], "Psi_SSym");
    HA.push_back(invPLC(StemPLCVEC[c]*100.0, VCstem_slope[c], VCstem_P50[c]), "Psi_SApo_cav"); //Sureau operates with %
    HA.push_back(invPLC(LeafPLCVEC[c]*100.0, VCleaf_slope[c], VCleaf_P50[c]), "Psi_LApo_cav"); //Sureau operates with %
    //PLC levels
    HA.push_back(StemPLCVEC[c]*100.0, "PLC_Stem"); //Sureau operates with %
    HA.push_back(LeafPLCVEC[c]*100.0, "PLC_Leaf"); //Sureau operates with %
    //Capacitances (mmol m-2 MPa-1)
    HA.push_back(2e-05, "C_SApo"); //Capacitance of the stem apoplasm
    HA.push_back(1e-05, "C_LApo"); //Capacitance of the leaf apoplasm
    HA.push_back(4761.905, "C_SSym"); //Capacitance of the stem symplasm (HOW TO ESTIMATE THEM?)
    HA.push_back(361.9048, "C_LSym"); //Capacitance of the leaf symplasm (HOW TO ESTIMATE THEM?)
    //Conductances (mmol m-2 MPa-1 s-1)
    HA.push_back(Plant_kmax[c], "k_Plant"); //Whole-plant conductance  = Plant_kmax
    HA.push_back(VCstem_kmax[c], "k_SLApo"); //Conductance from trunk apoplasm to leaf apoplasm = VCstem_kmax
    HA.push_back(0.26, "k_SSym"); //Conductance from stem apoplasm to stem symplasm (CONTROL PARAMETER?)
    HA.push_back(VCleaf_kmax[c], "k_LSym"); //Conductance from leaf apoplasm to leaf symplasm = VCleaf_kmax
    HA.push_back(VCroot_kmax(c,_), "k_SoilToStem"); //Conductance from rhizosphere to root crown?
    //Water content (mmol m-2)
    double l2mmol = 1.0e6/18.0;
    HA.push_back(Vsapwood[c]*StemAF[c]*l2mmol, "Q_SApo_sat_mmol_perLeafArea"); //Water content in stem apoplasm
    HA.push_back(Vleaf[c]*LeafAF[c]*l2mmol, "Q_LApo_sat_mmol_perLeafArea"); //Water content in leaf apoplasm
    HA.push_back(Vsapwood[c]*(1.0 - StemAF[c])*l2mmol, "Q_SSym_sat_mmol_perLeafArea"); //Water content in stem symplasm
    HA.push_back(Vleaf[c]*(1.0 - LeafAF[c])*l2mmol, "Q_LSym_sat_mmol_perLeafArea"); //Water content in leaf symplasm
    //Flows (mmol m-2 s-1)
    HA.push_back(Einst[c], "Einst"); //Total transpiration
    HA.push_back(Elim[c], "Elim"); //Stomatal transpiration
    HA.push_back(Emin_L[c], "Emin_L"); //Leaf cuticular transpiration
    HA.push_back(Emin_S[c], "Emin_S"); //Stem cuticular transpiration
    
    //Diagnostics
    HA.push_back(NA_INTEGER, "Diag_nwhile_cavit");
    

    HAs[c] = HA;
  }
  
  return(HAs);
}

// [[Rcpp::export]]
void semi_implicit_integration(List HA, List soil, 
                               double dt, int nsmalltimesteps, NumericVector opt) {
  
  List params = as<Rcpp::List>(HA["params"]);
  
  // Step 1. Initializing current time step according to computation options (FP)
  double dbxmin = 1.0e-100; // FP minimal double to avoid 0/0
  double Psi_LApo_n = HA["Psi_LApo"];
  double Psi_SApo_n = HA["Psi_SApo"];
  double Psi_LSym_n = HA["Psi_LSym"];
  double Psi_SSym_n = HA["Psi_SSym"];
  double Psi_LApo_cav = HA["Psi_LApo_cav"];
  double Psi_SApo_cav = HA["Psi_SApo_cav"];
  
  //Conductances
  double K_SL = HA["k_SLApo"];
  double k_SSym = HA["k_SSym"];
  double k_LSym = HA["k_LSym"];
  double c_LSym = HA["C_LSym"];
  double c_SSym = HA["C_SSym"];
  double c_LApo = HA["C_LApo"];
  double c_SApo = HA["C_SApo"];
  double PLC_Leaf = HA["PLC_Leaf"];
  double PLC_Stem = HA["PLC_Stem"];
  double Q_LApo_sat_mmol_perLeafArea = HA["Q_LApo_sat_mmol_perLeafArea"];
  double Q_SApo_sat_mmol_perLeafArea = HA["Q_SApo_sat_mmol_perLeafArea"];
  
  //Modifiers
  double Lsym = opt["Lsym"];
  double Ssym = opt["Ssym"];
  double CLapo = opt["CLapo"];
  double CTapo = opt["CTapo"];
  // double Eord = opt["Eord"]; MIQUEL: NOT USED
  double Lcav = opt["Lcav"];
  double Scav = opt["Scav"];

  
  //Apply modifiers
  double K_LSym = Lsym * k_LSym;   
  double K_SSym = Ssym * k_SSym;   
  double C_LSym = Lsym * c_LSym;   
  double C_SSym = Ssym * c_SSym;   
  double C_LApo = CLapo * c_LApo; 
  double C_SApo = CTapo * c_SApo; 
  
  double E_nph = HA["Elim"]; // Leaf stomatal transpiration
  double Emin_L_nph = HA["Emin_L"]; //Leaf cuticular transpiration
  double Emin_S_nph = HA["Emin_S"]; //Stem cuticular transpiration
  
  
  double VCleaf_slope = params["VCleaf_slope"];
  double VCstem_slope = params["VCstem_slope"];
  double VCleaf_P50 = params["VCleaf_P50"];
  double VCstem_P50 = params["VCstem_P50"];

  //Compute K_L_Cav et K_S_Cav
  double PLC_prime_L = PLC_derivative(PLC_Leaf, VCleaf_slope);
  double K_L_Cav = -1.0 * Lcav * Q_LApo_sat_mmol_perLeafArea * PLC_prime_L / dt;  // avec WBveg$Q_LSym_sat en l/m2 sol # changed by NM (25/10/2021)
  double PLC_prime_S = PLC_derivative(PLC_Stem, VCstem_slope);
  double K_S_Cav = -1.0 * Scav * Q_SApo_sat_mmol_perLeafArea * PLC_prime_S / dt;  // opt$Scav * WBveg$K_S_Cav #FP corrected a bug sign herehanged by NM (25/10/2021)


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

  NumericVector k_SoilToStem = HA["k_SoilToStem"];
  NumericVector PsiSoil = psi(soil, "VG");

  double alpha, Psi_td, Psi_LApo_np1, Psi_SApo_np1, Psi_LSym_np1, Psi_SSym_np1;
  double psiref;

  int nwhilecomp = 0; // # count the number of step in while loop (if more than 4 no solution and warning)
  while (((!LcavitWellComputed)||(!ScavitWellComputed)) && (nwhilecomp<delta_L_cavs.size())) {
    nwhilecomp = nwhilecomp + 1;
    double delta_L_cav = delta_L_cavs[nwhilecomp];
    double delta_S_cav = delta_S_cavs[nwhilecomp];

    //# 2.1 LApo
    alpha = exp(-1.0*(K_SL+K_LSym+delta_L_cav*K_L_Cav)/C_LApo*dt);
    Psi_td = (K_SL*Psi_SApo_n + K_LSym*Psi_LSym_n + delta_L_cav*K_L_Cav*Psi_LApo_cav)/(K_SL + K_LSym+delta_L_cav*K_L_Cav + dbxmin);// # dbxmin to avoid 0/0
    Psi_LApo_np1 = alpha * Psi_LApo_n + (1.0 - alpha) * Psi_td;

    //# 2.2. SApo
    alpha = exp(-1.0*(K_SL+K_SSym + sum(k_SoilToStem)+delta_S_cav*K_S_Cav)/C_SApo*dt);
    Psi_td = (K_SL*Psi_LApo_n + K_SSym*Psi_SSym_n + sum(k_SoilToStem * PsiSoil)+ delta_S_cav*K_S_Cav*Psi_SApo_cav)/(K_SL + K_SSym+sum(k_SoilToStem)+delta_S_cav*K_S_Cav + dbxmin);// # dbxmin to avoid 0/0
    Psi_SApo_np1 = alpha * Psi_SApo_n + (1.0 - alpha) * Psi_td;

    //# 2.3 Compute Psi_SApo_np1 (Explicit approach only)
    //# 2.4 check if cavitation is well computed according to delta_cav, np1 and "cav"
    LcavitWellComputed = (delta_L_cav==(Psi_LApo_np1 < Psi_LApo_cav)) || (Lcav==0.0);
    ScavitWellComputed = (delta_S_cav==(Psi_SApo_np1 < Psi_SApo_cav)) || (Scav==0.0);
    if ((delta_L_cavs.size() > 1) && (nwhilecomp==delta_L_cavs.size())) { //# we tried the normal cases and the computation is still not ok so we have done a last one desactivating cavitation water source (delta_cav=0)
      Rcerr << "water flux due to Cavitation ignored with time step, no solution from the implicit solver="<<dt<<"\n";
    }
  } //# end of the while loop with check on cavitation options
   
  HA["Diag_nwhile_cavit"] = nwhilecomp;  // # Diagnostic step to track cavit event and eventual errors (corresponding to nwhilecomp==5)

  //# Step 3. Compute Psi_Symp_np1 (L and S)
  alpha = exp(-1.0*K_LSym/C_LSym*dt);
  Psi_td = (K_LSym*Psi_LApo_n - (E_nph + Emin_L_nph))/(K_LSym + dbxmin);// # dbxmin to avoid 0/0
  Psi_LSym_np1 = alpha * Psi_LSym_n +(1.0 - alpha) * Psi_td;
  alpha = exp(-1.0*K_SSym/C_SSym*dt);
  Psi_td = (K_SSym*Psi_SApo_n - Emin_S_nph)/(K_SSym + dbxmin); // # dbxmin to avoid 0/0
  Psi_SSym_np1 = alpha * Psi_SSym_n +(1.0 - alpha) * Psi_td;

  //#Step 4 : set computed values in HA and update Psi_cav, PLC and Psi_AllSoil
  HA["Psi_LApo"] = std::min(-0.00001, Psi_LApo_np1);
  HA["Psi_SApo"] = std::min(-0.00001,Psi_SApo_np1);
  HA["Psi_LSym"] = std::min(-0.00001,Psi_LSym_np1);
  HA["Psi_SSym"] = std::min(-0.00001,Psi_SSym_np1);

  //# Cavitation
  psiref = HA["Psi_LApo"];  //# the reference is at current time step for other modes  (implicit, explicit)
  if(psiref < Psi_LApo_cav) {
    HA["Psi_LApo_cav"] = psiref;
    HA["PLC_Leaf"] = PLC(psiref, VCleaf_slope, VCleaf_P50);
  }

  psiref = HA["Psi_SApo"];  //# The reference is at current time step for other modes (implicit, explicit)
  if (psiref < Psi_SApo_cav) {
    HA["Psi_SApo_cav"] = psiref;
    HA["PLC_Stem"] = PLC(psiref, VCstem_slope, VCstem_P50);
  }


  // WBveg["Psi_AllSoil"] = sum(k_SoilToStem * PsiSoil)/sum(k_SoilToStem);
  // 
}

void innerCochard(List x, List input, List output, int n, double tstep, 
                 bool verbose = false, bool modifyInput = true) {
  
  // Extract control variables
  List control = x["control"];
  String soilFunctions = control["soilFunctions"];
  String cavitationRefill = control["cavitationRefill"];
  String rhizosphereOverlap = control["rhizosphereOverlap"];
  bool plantWaterPools = (rhizosphereOverlap!="total");

  DataFrame cohorts = Rcpp::as<Rcpp::DataFrame>(x["cohorts"]);
  int numCohorts = cohorts.nrow();
  
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  NumericVector LAIphe = Rcpp::as<Rcpp::NumericVector>(above["LAI_expanded"]);
  NumericVector N = Rcpp::as<Rcpp::NumericVector>(above["N"]);
  
  List soil = x["soil"];
  NumericVector dVec = soil["dVec"];
  NumericVector Water_FC = waterFC(soil, soilFunctions);
  // NumericVector Theta_FC = thetaFC(soil, soilFunctions);
  // NumericVector VG_n = Rcpp::as<Rcpp::NumericVector>(soil["VG_n"]);
  // NumericVector VG_alpha = Rcpp::as<Rcpp::NumericVector>(soil["VG_alpha"]);
  NumericVector Tsoil = soil["Temp"]; 
  int nlayers = Tsoil.length();
  
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
  NumericVector leafWidth = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["LeafWidth"]);
  
  DataFrame paramsTransp = Rcpp::as<Rcpp::DataFrame>(x["paramsTranspiration"]);
  NumericVector Gswmin = Rcpp::as<Rcpp::NumericVector>(paramsTransp["Gswmin"]);
  NumericVector Gswmax = Rcpp::as<Rcpp::NumericVector>(paramsTransp["Gswmax"]);
  NumericVector VCstem_kmax = Rcpp::as<Rcpp::NumericVector>(paramsTransp["VCstem_kmax"]);
  NumericVector VCstem_c = Rcpp::as<Rcpp::NumericVector>(paramsTransp["VCstem_c"]);
  NumericVector VCstem_d = Rcpp::as<Rcpp::NumericVector>(paramsTransp["VCstem_d"]);
  NumericVector VCroot_kmax_sum = Rcpp::as<Rcpp::NumericVector>(paramsTransp["VCroot_kmax"]);
  NumericVector VCroot_c = paramsTransp["VCroot_c"];
  NumericVector VCroot_d = paramsTransp["VCroot_d"];
  NumericVector VCleaf_kmax = Rcpp::as<Rcpp::NumericVector>(paramsTransp["VCleaf_kmax"]);
  NumericVector VCleaf_c = Rcpp::as<Rcpp::NumericVector>(paramsTransp["VCleaf_c"]);
  NumericVector VCleaf_d = Rcpp::as<Rcpp::NumericVector>(paramsTransp["VCleaf_d"]);
  
  DataFrame paramsWaterStorage = Rcpp::as<Rcpp::DataFrame>(x["paramsWaterStorage"]);
  NumericVector StemPI0 = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["StemPI0"]);
  NumericVector StemEPS = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["StemEPS"]);
  NumericVector StemAF = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["StemAF"]);
  NumericVector Vsapwood = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["Vsapwood"]); //l·m-2 = mm
  NumericVector LeafPI0 = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["LeafPI0"]);
  NumericVector LeafEPS = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["LeafEPS"]);
  NumericVector LeafAF = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["LeafAF"]);
  NumericVector Vleaf = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["Vleaf"]); //l·m-2 = mm
  
  //Extract internal variables
  // Rcout<<"internal\n";
  DataFrame internalWater = Rcpp::as<Rcpp::DataFrame>(x["internalWater"]);
  NumericVector StemPLCVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["StemPLC"]);
  NumericVector Stem1PsiVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["Stem1Psi"]);
  NumericVector Stem2PsiVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["Stem2Psi"]);
  NumericVector EinstVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["Einst"]);
  NumericVector LeafPsiVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["LeafPsi"]);
  NumericVector RootCrownPsiVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["RootCrownPsi"]);
  NumericVector StemSympPsiVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["StemSympPsi"]);
  NumericVector LeafSympPsiVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["LeafSympPsi"]);
  
  //Water pools
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
  
  // Rcout<<"EB\n";
  List EB = output["EnergyBalance"];
  DataFrame Tinst = Rcpp::as<Rcpp::DataFrame>(EB["Temperature"]);
  DataFrame CEBinst = Rcpp::as<Rcpp::DataFrame>(EB["CanopyEnergyBalance"]);
  DataFrame SEBinst = Rcpp::as<Rcpp::DataFrame>(EB["SoilEnergyBalance"]);
  
  NumericMatrix SoilWaterExtract = Rcpp::as<Rcpp::NumericMatrix>(output["Extraction"]);
  NumericMatrix soilLayerExtractInst = Rcpp::as<Rcpp::NumericMatrix>(output["ExtractionInst"]);
  
  NumericMatrix minPsiRhizo = Rcpp::as<Rcpp::NumericMatrix>(output["RhizoPsi"]);
  
  // Rcout<<"Plants\n";
  List Plants = output["Plants"];
  NumericVector PWB = Plants["WaterBalance"];
  NumericVector Eplant = Plants["Transpiration"];
  NumericVector Agplant = Plants["GrossPhotosynthesis"];
  NumericVector Anplant = Plants["NetPhotosynthesis"];
  NumericVector minLeafPsi = Plants["LeafPsiMin"];
  NumericVector maxLeafPsi = Plants["LeafPsiMax"];
  NumericVector minStemPsi = Plants["StemPsi"];
  NumericVector minRootPsi = Plants["RootPsi"];
  
  // Rcout<<"Leaves\n";
  List Sunlit = output["SunlitLeaves"];
  List Shade = output["ShadeLeaves"];
  NumericVector LAI_SL = Sunlit["LAI"];
  NumericVector Vmax298SL = Sunlit["Vmax298"];
  NumericVector Jmax298SL = Sunlit["Jmax298"];
  NumericVector maxGSW_SL = Sunlit["GSWMax"];
  NumericVector minGSW_SL = Sunlit["GSWMin"];
  NumericVector minTemp_SL = Sunlit["TempMin"];
  NumericVector maxTemp_SL = Sunlit["TempMax"];
  NumericVector minLeafPsi_SL = Sunlit["LeafPsiMin"];
  NumericVector maxLeafPsi_SL = Sunlit["LeafPsiMax"];
  
  NumericVector LAI_SH = Shade["LAI"];
  NumericVector Vmax298SH = Shade["Vmax298"];
  NumericVector Jmax298SH = Shade["Jmax298"];
  NumericVector maxGSW_SH = Shade["GSWMax"];
  NumericVector minGSW_SH = Shade["GSWMin"];
  NumericVector minTemp_SH = Shade["TempMin"];
  NumericVector maxTemp_SH = Shade["TempMax"];
  NumericVector minLeafPsi_SH = Shade["LeafPsiMin"];
  NumericVector maxLeafPsi_SH = Shade["LeafPsiMax"];
  
  // Rcout<<"PlantsInst\n";
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
  NumericMatrix PLC = Rcpp::as<Rcpp::NumericMatrix>(PlantsInst["StemPLC"]);
  NumericMatrix StemSympPsiInst = Rcpp::as<Rcpp::NumericMatrix>(PlantsInst["StemSympPsi"]);
  NumericMatrix LeafSympPsiInst = Rcpp::as<Rcpp::NumericMatrix>(PlantsInst["LeafSympPsi"]);
  
  List ShadeInst = output["ShadeLeavesInst"];
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
  IntegerVector iLayerCohort = input["iLayerCohort"];
  IntegerVector iLayerSunlit = input["iLayerSunlit"];
  IntegerVector iLayerShade = input["iLayerShade"];
  IntegerVector nlayerscon = input["nlayerscon"];
  LogicalMatrix layerConnected = input["layerConnected"];
  List layerConnectedPools = input["layerConnectedPools"];

  for(int c=0;c<numCohorts;c++) { //Plant cohort loop
    // 
    // 
    // // # A. LOOP ON THE IMPLICIT SOLVER IN PSI, trying different time steps until results are OK
    // bool regulationWellComputed = false;
    // bool cavitationWellComputed = false;
    // List params  = WBveg["params"];
    // NumericVector kSoil = WBsoil["kSoil"];
    // 
    // List WBveg_np1;
    // 
    // double LAI = (double) WBveg["LAI"];
    // 
    // int nwhilecomp = 0;
    // NumericVector fluxSoilToStemLargeTimeStep(kSoil.size(), 0.0);
    // 
    // while ((!regulationWellComputed || !cavitationWellComputed) && (nwhilecomp<nsmalltimesteps.size())) { //# LOOP TO TRY DIFFERENT TIME STEPS
    //   List WBveg_n = clone(WBveg); // # initial value of WBveg
    //   List WBsoil_n = clone(WBsoil); // # initial value of WBsoil
    //   
    //   regulationWellComputed = false;
    //   cavitationWellComputed = false;
    //   int nwhilecomp = nwhilecomp + 1;
    //   double deltaRegulMax = 1.0e-100;
    //   double deltaPLCMax = 1.0e-100;
    //   
    //   double fluxEvaporationSoilLargeTimeStep = 0.0;
    //   for(int i=0;i < kSoil.size();i++) fluxSoilToStemLargeTimeStep[i] = 0.0;
    //   
    //   int nts = nsmalltimesteps[nwhilecomp];// # number of small time steps
    //   for(int its = 1; its <= nts; its++) { //#INTERNAL LOOP ON SMALL TIME STEPS
    //     double p = (((double) its ) - 0.5)/((double) nts);
    //     List WBclim = interp_WBclim(WBclim_current, WBclim_next, p); // # climate at nph
    //     double ETPr = WBveg_n["ETPr"];
    //     double Tair_mean = WBclim["Tair_mean"];
    //     double RHair = WBclim["RHair"];
    //     compute_evaporationG(WBsoil_n, RHair, Tair_mean,
    //                          ((double) Nhours)/((double) nts), LAI,
    //                          ETPr, (double) params["K"]);
    //     fluxEvaporationSoilLargeTimeStep = fluxEvaporationSoilLargeTimeStep + ((double) WBsoil_n["E_Soil3"])/((double) nts);
    //     
    //     WBveg_np1 = clone(WBveg_n); // Clone WBveg object
    //     compute_transpiration(WBveg_np1, WBclim, Nhours, opt, stomatalRegFormulation);// # transpi with climate at nph
    //     semi_implicit_temporal_integration(WBveg_np1,  WBsoil_n, ((double) Nhours) * 3600.0 / ((double) nts), nts, opt);
    //     update_kplant(WBveg_np1, WBsoil_n);
    //     update_capacitancesApoAndSym(WBveg_np1);
    //     
    //     // # QUANTITIES TO CHECK IF THE RESOLUTION IS OK
    //     // # 1. delta regulation between n and np1 (MIQUEL: Only Psi_LSym changes between the two calculations, params should be the same)
    //     NumericVector regul_np1 = regulFact_comp(WBveg_np1["Psi_LSym"], params, stomatalRegFormulation);
    //     NumericVector regul_n = regulFact_comp(WBveg_n["Psi_LSym"], params, stomatalRegFormulation); //# TODO check why recomputed? should be in WBveg_tmp
    //     
    //     deltaRegulMax = std::max(deltaRegulMax,std::abs((double) regul_np1["regulFact"] - (double) regul_n["regulFact"]));
    //     
    //     // # 2. PLC at n and np1
    //     deltaPLCMax = std::max(deltaPLCMax, (double) WBveg_np1["PLC_Leaf"] - (double) WBveg_n["PLC_Leaf"]);
    //     deltaPLCMax = std::max(deltaPLCMax, (double) WBveg_np1["PLC_Stem"] - (double) WBveg_n["PLC_Stem"]);
    //     WBveg_n = WBveg_np1; //# Update WBveg_n
    //     
    //     // # 3. update of soil on small time step (done by FP in version 16)
    //     NumericVector fluxSoilToStem_mm = WBveg_np1["fluxSoilToStem"];
    //     double Psi_SApo = WBveg_np1["Psi_SApo"];
    //     NumericVector k_SoilToStem = WBveg["k_SoilToStem"]; //MIQUEL: Why WBveg here?
    //     NumericVector PsiSoil = WBsoil_n["PsiSoil"];
    //     for(int i=0;i < kSoil.size();i++) {
    //       double fluxSoilToStem_mmolm2s = k_SoilToStem[i]*(PsiSoil[i] - Psi_SApo);
    //       fluxSoilToStemLargeTimeStep[i] = fluxSoilToStemLargeTimeStep[i] + fluxSoilToStem_mmolm2s/((double) nts);// # mean flux over one large time step
    //       fluxSoilToStem_mm[i] = convertFluxFrom_mmolm2s_To_mm(fluxSoilToStem_mmolm2s, ((double) Nhours)/((double) nts), LAI); // # Quantity from each soil layer to the below part
    //     }
    //     // # NB the time step for fluxSoilToStem_mm is Nhours/nts!
    //     update_soilWater(WBsoil_n, fluxSoilToStem_mm);
    //   } //# end loop small time step
    //   
    //   // # TESTS ON RESOLUTION
    //   WBveg_np1["Diag_deltaRegulMax"] = deltaRegulMax;
    //   regulationWellComputed = (deltaRegulMax<0.05);
    //   WBveg_np1["Diag_deltaPLCMax"] = deltaPLCMax;
    //   cavitationWellComputed = (deltaPLCMax<1.0);// # 1%
    //   WBveg_np1["Diag_timeStepInHours"] = ((double) Nhours)/((double) nts);
    // } //# end while
    // 
    // // # B. SAVING SOLUTION AT NEXT TIME STEP IN WBveg
    // WBveg = WBveg_np1;
    // // # final update of transpiration at clim_next (useful for consistency in outputs, but not required for the computations)
    // compute_transpiration(WBveg, WBclim_next, Nhours, opt, stomatalRegFormulation);
    // // # C. UPDATING FLUX FROM SOIL (WBveg$fluxSoilToStem is used as input in UpdateSoilWater.WBsoil)
    // // #TODO FP suggests moving the computation of  fluxSoilToStem in the main loop, as it is the coupling between the two models...
    // // # mean soil quantities on large time steps
    // WBveg["Emin_mm"]  = convertFluxFrom_mmolm2s_To_mm((double) WBveg["Emin"], (double) Nhours, LAI); // # Flux from each soil layer to the below part
    // WBveg["Emin_S_mm"] = convertFluxFrom_mmolm2s_To_mm((double) WBveg["Emin_S"], (double) Nhours, LAI); // # Flux from each soil layer to the below part
    // 
    // double SumFluxSoilToStem = sum(fluxSoilToStemLargeTimeStep); // # flux total en mmol/m2/s / used for Tleaf
    // WBveg["SumFluxSoilToStem"] = SumFluxSoilToStem;
    // NumericVector fluxSoilToStem_mm = WBveg["fluxSoilToStem_mm"];
    // for(int i=0;i< fluxSoilToStem_mm.size();i++) {
    //   fluxSoilToStem_mm[i] = convertFluxFrom_mmolm2s_To_mm(fluxSoilToStemLargeTimeStep[i], (double) Nhours, LAI); // # Flux from each soil layer to the below part  in mm
    // }
    // WBveg["transpiration_mm"] = convertFluxFrom_mmolm2s_To_mm(((double) WBveg["Emin"]) + ((double) WBveg["Emin_S"]) + ((double) WBveg["Elim"]),
    //                                            (double) Nhours, LAI); //# total flux in mm
    
  }

}
