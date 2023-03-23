#define STRICT_R_HEADERS
#include <Rcpp.h>
#include "photosynthesis.h"
#include "biophysicsutils.h"
#include "hydraulics.h"
#include "soil.h"
#include "tissuemoisture.h"
#include <meteoland.h>
using namespace Rcpp;

double PLC_derivative(double plc, double slope) {
  return(-1.0*slope/25.0 * plc/100 * (1.0 - plc/100));
}
double PLC(double Pmin, double slope, double P50){
  return (100.0 / (1.0 + exp(slope / 25.0 * (Pmin - P50))));
}
double invPLC(double plc, double slope, double P50){
  return(P50 + log((100.0/plc)-1.0)*(25.0/slope));
}
double RWC(double PiFT, double Esymp, double Pmin) {
  double A = std::max((-1.0 * (Pmin + PiFT - Esymp) - sqrt(pow(Pmin + PiFT - Esymp, 2.0) + 4.0 * (Pmin * Esymp))) / (2.0 * Esymp), 1.0 - PiFT / Pmin);
  return(A);
}
double gCrown(double gCrown0, double windSpeed){
  windSpeed=  std::max(0.1, windSpeed); //# to avoid very high conductance values 
  return(gCrown0*pow(windSpeed,0.6));
}
double Emin(double gmin, double gBL, double gCrown, 
            double VPD, double airPressure =101.3) {
  double gmintot = 1.0/(1.0/gmin+ 1.0/gBL + 1.0/gCrown);
  return(gmintot * VPD /airPressure); 
}
double gmin(double leafTemperature, double gmin_20, 
            double TPhase, double Q10_1, double Q10_2) {
  double gmin = NA_REAL;
  if (leafTemperature<= TPhase) {
    gmin = gmin_20 * pow(Q10_1,(leafTemperature - 20.0) / 10.0);
  } else if (leafTemperature > TPhase) {
    gmin = gmin_20 * pow(Q10_1, (TPhase - 20.0) / 10.0) * pow(Q10_2, (leafTemperature- TPhase) / 10.0);
  }
  return(gmin);
}

// # Update plant conductances
void update_conductances(List network) {
  List params = as<Rcpp::List>(network["params"]);
  
  NumericVector k_RSApoInit = params["k_RSApoInit"];
  
  NumericVector k_RSApo = network["k_RSApo"];
  NumericVector k_SoilToStem = network["k_SoilToStem"];
  NumericVector k_Soil = network["k_Soil"];
  
  network["k_SLApo"] = ((double) params["k_SLApoInit"]) * (1.0 - ((double) network["PLC_Leaf"])/100.0);

  for(int i = 0;i<k_RSApo.size();i++) {
    //# calculate k_RSApo and k_SLApo with cavitation
    k_RSApo[i] = k_RSApoInit[i] * (1.0 - ((double) network["PLC_Stem"])/100.0);
    //# Root from root length
    k_SoilToStem[i] = 1.0/((1.0/k_Soil[i]) + (1.0/k_RSApo[i])); // # conductance from soil to collar (two resistances in series Rsoil and Rroot)
  }

  // Compute k_plant (from root to leaf) for diagnostic only
  network["k_Plant"] =  1.0/ (1.0 /sum(k_RSApo) + 1.0/((double) network["k_SLApo"]) + 1.0/((double) network["k_LSym"]));
}

// # update symplasmic plant capacitances for Trunk and leaves
void update_capacitances(List network) {
  List params = as<Rcpp::List>(network["params"]);
  double dbxmin = 1.0e-100; //# NM minimal double to avoid-INF
  
  double LAI = network["LAI"];
  double Psi_SSym = network["Psi_SSym"];
  double Psi_LSym = network["Psi_LSym"];
  double Q_LSym_sat_mmol_perLeafArea = network["Q_LSym_sat_mmol_perLeafArea"];
  
  double epsilonSym_Leaf = params["epsilonSym_Leaf"];
  double PiFullTurgor_Leaf = params["PiFullTurgor_Leaf"];
  double epsilonSym_Stem = params["epsilonSym_Stem"];
  double PiFullTurgor_Stem = params["PiFullTurgor_Stem"];
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
    network["C_LSym"] = 0.0;
  } else {
    network["C_LSym"] = Q_LSym_sat_mmol_perLeafArea * RWC_LSym_prime;
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
  network["C_SSym"] = Q_LSym_sat_mmol_perLeafArea * RWC_SSym_prime; // #  changed 25/10/2021 by NM. --> Stem capacitance per leaf area can only decrease with LAI (cannot increase when LAI<1 )
  network["C_SApo"] = params["C_SApoInit"];
  network["C_LApo"] = params["C_LApoInit"];
}

double Turgor(double PiFT, double Esymp, double Rstemp) {
  return(-1.0*PiFT - Esymp * Rstemp);
}

NumericVector regulFact(double psi, List params, String regulationType = "Turgor") {
  double stomatalClosure = NA_REAL;
  double regulFact = NA_REAL;
  double regulFactPrime = NA_REAL;
  
  if(regulationType == "PiecewiseLinear") {
    double PsiStartClosing = params["PsiStartClosing"];
    double PsiClose = params["PsiClose"];
    if (psi > PsiStartClosing) {
      stomatalClosure = 0.0;
      regulFact = 1.0;
      regulFactPrime = 0.0;
    } else if (psi > PsiClose) {
      stomatalClosure = 1.0;
      regulFact = (psi - PsiClose) / (PsiStartClosing - PsiClose);
      regulFactPrime = 1.0 / (PsiStartClosing - PsiClose);
    } else {
      stomatalClosure = 2.0;
      regulFact = 0.0;
      regulFactPrime = 0.0;
    }
  } else if (regulationType == "Sigmoid") {
    double slope_gs = params["slope_gs"];
    double P50_gs = params["P50_gs"];
    double PL_gs = 1.0 / (1.0 + exp(slope_gs / 25.0 * (psi - P50_gs)));
    regulFact = 1.0 - PL_gs;
    double al = slope_gs / 25.0;
    regulFactPrime = al * PL_gs * regulFact;
  } else if (regulationType == "Turgor") {
    double turgorPressureAtGsMax = params["turgorPressureAtGsMax"];
    double epsilonSym_Leaf = params["epsilonSym_Leaf"];
    double PiFullTurgor_Leaf = params["PiFullTurgor_Leaf"];
    double rs = RWC(PiFullTurgor_Leaf, epsilonSym_Leaf, psi);
    double turgor = Turgor(PiFullTurgor_Leaf, epsilonSym_Leaf, rs);
    regulFact = std::max(0.0, std::min(1.0, turgor / turgorPressureAtGsMax));
    if((regulFact == 1.0) || (regulFact == 0.0)) {
      regulFactPrime = 0.0;
    } else {
      regulFactPrime = 0.0;// # TODO insert the derivative to compute regulFactPrime
    }
  }
  NumericVector res = NumericVector::create(_["regulFact"] = regulFact, 
                                            _["stomatalClosure"] = stomatalClosure, 
                                            _["regulFactPrime"] = regulFactPrime);
  return(res);
}

//# stomatal conductance calculation with Jarvis type formulations
double gsJarvis(List params, double PAR, double Temp, int option = 1){
  double JarvisPAR = params["JarvisPAR"];
  double gsMax = params["gsMax"];
  double gsNight = params["gsNight"];
  double Tgs_optim = params["Tgs_optim"];
  double Tgs_sens = params["Tgs_sens"];
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

List initCochardNetwork(int c, NumericVector LAIphe,
                       DataFrame internalWater, 
                       DataFrame paramsAnatomy, DataFrame paramsTranspiration, DataFrame paramsWaterStorage,
                       NumericVector PsiSoil, NumericVector VCroot_kmax, NumericVector VGrhizo_kmax,
                       double sapFluidityDay = 1.0) {
  //Root distribution input
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

  NumericVector Plant_kmax = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["Plant_kmax"]);
  NumericVector VCleaf_kmax = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["VCleaf_kmax"]);
  NumericVector VCleaf_P50 = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["VCleaf_P50"]);
  NumericVector VCleaf_slope = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["VCleaf_slope"]);
  NumericVector VCstem_kmax = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["VCstem_kmax"]);
  NumericVector VCstem_P50 = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["VCstem_P50"]);
  NumericVector VCstem_slope = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["VCstem_slope"]);
  NumericVector VCroot_P50 = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["VCroot_P50"]);
  NumericVector VCroot_slope = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["VCroot_slope"]);
  NumericVector Gswmax = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["Gswmax"]);
  NumericVector Gswmin = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["Gswmin"]);
  
  
  NumericVector StemPI0 = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["StemPI0"]);
  NumericVector StemEPS = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["StemEPS"]);
  NumericVector StemAF = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["StemAF"]);
  NumericVector Vsapwood = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["Vsapwood"]); //l·m-2 = mm
  NumericVector LeafPI0 = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["LeafPI0"]);
  NumericVector LeafEPS = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["LeafEPS"]);
  NumericVector LeafAF = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["LeafAF"]);
  NumericVector Vleaf = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["Vleaf"]); //l·m-2 = mm
  
  
  List network = List::create();
  //Params
  List params = List::create();
  // CONSTANTS TO BE REVISED
  params.push_back(37.5, "TPhase_gmin"); 
  params.push_back(1.2, "Q10_1_gmin"); 
  params.push_back(4.8, "Q10_2_gmin"); 
  
  params.push_back(1.15, "turgorPressureAtGsMax"); 
  
  params.push_back(0.003, "JarvisPAR"); 
  params.push_back(20.0, "gsNight"); //mmol
  params.push_back(25.0, "Tgs_optim"); 
  params.push_back(17.0, "Tgs_sens"); 
  
  params.push_back(150.0, "gCrown0"); 
  params.push_back(2.0, "gmin_S"); 
  params.push_back(0.8, "fTRBToLeaf");
  
  //PLANT RELATED PARAMETERS
  params.push_back(Gswmin[c]*1000.0, "gmin20"); //mol to mmol 
  params.push_back(Gswmax[c]*1000.0, "gsMax"); //mol to mmol 
  
  params.push_back(VCleaf_P50[c], "VCleaf_P50"); 
  params.push_back(VCleaf_slope[c], "VCleaf_slope"); 
  params.push_back(VCstem_P50[c], "VCstem_P50"); 
  params.push_back(VCstem_slope[c], "VCstem_slope"); 
  params.push_back(VCroot_P50[c], "VCrootP50"); 
  params.push_back(VCroot_slope[c], "VCroot_slope"); 
  params.push_back(sapFluidityDay*VCstem_kmax[c], "k_SLApoInit"); //Maximum conductance from trunk apoplasm to leaf apoplasm = VCstem_kmax
  params.push_back(sapFluidityDay*VCroot_kmax, "k_RSApoInit"); //Maximum conductance from rhizosphere surface to root crown
  params.push_back(LeafPI0[c], "PiFullTurgor_Leaf"); 
  params.push_back(LeafEPS[c], "epsilonSym_Leaf"); 
  params.push_back(StemPI0[c], "PiFullTurgor_Stem"); 
  params.push_back(StemEPS[c], "epsilonSym_Stem"); 
  params.push_back(2e-05, "C_SApoInit"); //Maximum capacitance of the stem apoplasm
  params.push_back(1e-05, "C_LApoInit"); //Maximum capacitance of the leaf apoplasm
  
  network.push_back(params, "params");
  
  //LAI
  network.push_back(LAIphe[c], "LAI");
  
  //Water potentials
  network.push_back(StemPsiVEC[c], "Psi_LApo"); // Leaf apo psi  = Stem psi
  network.push_back(LeafPsiVEC[c], "Psi_LSym"); // Leaf sym psi = Leaf psi
  network.push_back(RootCrownPsiVEC[c], "Psi_SApo"); // Stem apo psi = root crown psi
  network.push_back(StemSympPsiVEC[c], "Psi_SSym");
  network.push_back(std::min(0.0, invPLC(StemPLCVEC[c]*100.0, VCstem_slope[c], VCstem_P50[c])), "Psi_SApo_cav"); //Sureau operates with %
  network.push_back(std::min(0.0, invPLC(LeafPLCVEC[c]*100.0, VCleaf_slope[c], VCleaf_P50[c])), "Psi_LApo_cav"); //Sureau operates with %
  network.push_back(PsiSoil, "PsiSoil"); // Leaf apo psi  = Stem psi
  //PLC levels
  network.push_back(StemPLCVEC[c]*100.0, "PLC_Stem"); //Sureau operates with %
  network.push_back(LeafPLCVEC[c]*100.0, "PLC_Leaf"); //Sureau operates with %
  //Capacitances (mmol m-2 MPa-1)
  network.push_back(NA_REAL, "C_SApo"); //Capacitance of the stem apoplasm
  network.push_back(NA_REAL, "C_LApo"); //Capacitance of the leaf apoplasm
  network.push_back(NA_REAL, "C_SSym"); //Capacitance of the stem symplasm (HOW TO ESTIMATE THEM?)
  network.push_back(NA_REAL, "C_LSym"); //Capacitance of the leaf symplasm (HOW TO ESTIMATE THEM?)
  //Conductances (mmol m-2 MPa-1 s-1)
  network.push_back(NA_REAL, "k_Plant"); //Whole-plant conductance  = Plant_kmax
  network.push_back(NA_REAL, "k_SLApo"); //Conductance from trunk apoplasm to leaf apoplasm = VCstem_kmax
  network.push_back(sapFluidityDay*0.26, "k_SSym"); //Conductance from stem apoplasm to stem symplasm (CONTROL PARAMETER?)
  network.push_back(sapFluidityDay*VCleaf_kmax[c], "k_LSym"); //Conductance from leaf apoplasm to leaf symplasm = VCleaf_kmax
  NumericVector k_RSApo(VGrhizo_kmax.size(), NA_REAL);
  network.push_back(k_RSApo, "k_RSApo"); //Conductance from rhizosphere surface to root crown?
  NumericVector k_SoilToStem(VGrhizo_kmax.size(), NA_REAL);
  network.push_back(k_SoilToStem, "k_SoilToStem"); //Conductance from rhizosphere surface to root crown?
  network.push_back(VGrhizo_kmax, "k_Soil"); //Conductance in the rhizosphere
  //Water content (mmol m-2)
  double l2mmol = 1.0e6/18.0;
  network.push_back(Vsapwood[c]*StemAF[c]*l2mmol, "Q_SApo_sat_mmol_perLeafArea"); //Water content in stem apoplasm
  network.push_back(Vleaf[c]*LeafAF[c]*l2mmol, "Q_LApo_sat_mmol_perLeafArea"); //Water content in leaf apoplasm
  network.push_back(Vsapwood[c]*(1.0 - StemAF[c])*l2mmol, "Q_SSym_sat_mmol_perLeafArea"); //Water content in stem symplasm
  network.push_back(Vleaf[c]*(1.0 - LeafAF[c])*l2mmol, "Q_LSym_sat_mmol_perLeafArea"); //Water content in leaf symplasm
  //Flows (mmol m-2 s-1)
  network.push_back(Einst[c], "Einst"); //Total transpiration
  network.push_back(Elim[c], "Elim"); //Stomatal transpiration
  network.push_back(NA_REAL, "Elim_SL"); //Stomatal transpiration (sunlit leaves)
  network.push_back(NA_REAL, "Elim_SH"); //Stomatal transpiration (shade leaves)
  network.push_back(Emin_L[c], "Emin_L"); //Leaf cuticular transpiration
  network.push_back(NA_REAL, "Emin_L_SL"); //Leaf cuticular transpiration (sunlit leaves)
  network.push_back(NA_REAL, "Emin_L_SH"); //Leaf cuticular transpiration (shade leaves)
  network.push_back(Emin_S[c], "Emin_S"); //Stem cuticular transpiration
  
  //Diagnostics
  network.push_back(NA_INTEGER, "Diag_nwhile_cavit");
  network.push_back(NA_INTEGER, "Diag_deltaRegulMax");
  network.push_back(NA_INTEGER, "Diag_deltaPLCMax");
  network.push_back(NA_INTEGER, "Diag_timeStepInSeconds");
  
  // Update plant conductances and capacitances according to network status
  update_conductances(network);
  update_capacitances(network);
  return(network);
}

// Initializes network for all plant cohorts in x
// [[Rcpp::export]]
List initCochardNetworks(List x) {
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  NumericVector LAIphe = Rcpp::as<Rcpp::NumericVector>(above["LAI_expanded"]);
  
  List belowLayers = Rcpp::as<Rcpp::List>(x["belowLayers"]);
  NumericMatrix VCroot_kmax= Rcpp::as<Rcpp::NumericMatrix>(belowLayers["VCroot_kmax"]);
  NumericMatrix VGrhizo_kmax= Rcpp::as<Rcpp::NumericMatrix>(belowLayers["VGrhizo_kmax"]);
  DataFrame internalWater = Rcpp::as<Rcpp::DataFrame>(x["internalWater"]);
  DataFrame paramsAnatomy = Rcpp::as<Rcpp::DataFrame>(x["paramsAnatomy"]);
  DataFrame paramsTranspiration = Rcpp::as<Rcpp::DataFrame>(x["paramsTranspiration"]);
  DataFrame paramsWaterStorage = Rcpp::as<Rcpp::DataFrame>(x["paramsWaterStorage"]);
  
  List soil = x["soil"];
  NumericVector psiSoil = psi(soil, "VG");
  int numCohorts = internalWater.nrow();
  List networks(numCohorts);
  for(int c = 0;c<numCohorts;c++) {
    networks[c] = initCochardNetwork(c, LAIphe,
                                     internalWater, 
                                     paramsAnatomy, paramsTranspiration, paramsWaterStorage,
                                     psiSoil, VCroot_kmax(c,_), VGrhizo_kmax(c,_));
  }
  return(networks);
}


// dt - Smallest time step (seconds)
// opt - Option flag vector
// [[Rcpp::export]]
void semi_implicit_integration(List network, double dt, NumericVector opt) {
  
  List params = as<Rcpp::List>(network["params"]);
  NumericVector PsiSoil = network["PsiSoil"];
  
  // Step 1. Initializing current time step according to computation options (FP)
  double dbxmin = 1.0e-100; // FP minimal double to avoid 0/0
  double Psi_LApo_n = network["Psi_LApo"];
  double Psi_SApo_n = network["Psi_SApo"];
  double Psi_LSym_n = network["Psi_LSym"];
  double Psi_SSym_n = network["Psi_SSym"];
  double Psi_LApo_cav = network["Psi_LApo_cav"];
  double Psi_SApo_cav = network["Psi_SApo_cav"];
  
  //Conductances
  double K_SL = network["k_SLApo"];
  double k_SSym = network["k_SSym"];
  double k_LSym = network["k_LSym"];
  double c_LSym = network["C_LSym"];
  double c_SSym = network["C_SSym"];
  double c_LApo = network["C_LApo"];
  double c_SApo = network["C_SApo"];
  double PLC_Leaf = network["PLC_Leaf"];
  double PLC_Stem = network["PLC_Stem"];
  double Q_LApo_sat_mmol_perLeafArea = network["Q_LApo_sat_mmol_perLeafArea"];
  double Q_SApo_sat_mmol_perLeafArea = network["Q_SApo_sat_mmol_perLeafArea"];
  
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
  
  double E_nph = network["Elim"]; // Leaf stomatal transpiration
  double Emin_L_nph = network["Emin_L"]; //Leaf cuticular transpiration
  double Emin_S_nph = network["Emin_S"]; //Stem cuticular transpiration
  
  
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

  NumericVector k_SoilToStem = network["k_SoilToStem"];

  double alpha, Psi_td, Psi_LApo_np1, Psi_SApo_np1, Psi_LSym_np1, Psi_SSym_np1;
  double psiref;

  int nwhilecomp = 0; // # count the number of step in while loop (if more than 4 no solution and warning)
  while (((!LcavitWellComputed)||(!ScavitWellComputed)) && (nwhilecomp<delta_L_cavs.size())) {
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
    
    nwhilecomp = nwhilecomp + 1;
    if ((delta_L_cavs.size() > 1) && (nwhilecomp==delta_L_cavs.size())) { //# we tried the normal cases and the computation is still not ok so we have done a last one desactivating cavitation water source (delta_cav=0)
      Rcerr << "water flux due to Cavitation ignored with time step, no solution from the implicit solver="<<dt<<"\n";
    }
  } //# end of the while loop with check on cavitation options
   
  network["Diag_nwhile_cavit"] = nwhilecomp;  // # Diagnostic step to track cavit event and eventual errors (corresponding to nwhilecomp==5)

  //# Step 3. Compute Psi_Symp_np1 (L and S)
  alpha = exp(-1.0*K_LSym/C_LSym*dt);
  Psi_td = (K_LSym*Psi_LApo_n - (E_nph + Emin_L_nph))/(K_LSym + dbxmin);// # dbxmin to avoid 0/0
  Psi_LSym_np1 = alpha * Psi_LSym_n +(1.0 - alpha) * Psi_td;
  alpha = exp(-1.0*K_SSym/C_SSym*dt);
  Psi_td = (K_SSym*Psi_SApo_n - Emin_S_nph)/(K_SSym + dbxmin); // # dbxmin to avoid 0/0
  Psi_SSym_np1 = alpha * Psi_SSym_n +(1.0 - alpha) * Psi_td;

  //#Step 4 : set computed values in network and update Psi_cav, PLC and Psi_AllSoil
  network["Psi_LApo"] = std::min(-0.00001, Psi_LApo_np1);
  network["Psi_SApo"] = std::min(-0.00001,Psi_SApo_np1);
  network["Psi_LSym"] = std::min(-0.00001,Psi_LSym_np1);
  network["Psi_SSym"] = std::min(-0.00001,Psi_SSym_np1);

  //# Cavitation
  psiref = network["Psi_LApo"];  //# the reference is at current time step for other modes  (implicit, explicit)
  if(psiref < Psi_LApo_cav) {
    network["Psi_LApo_cav"] = psiref;
    network["PLC_Leaf"] = PLC(psiref, VCleaf_slope, VCleaf_P50);
  }

  psiref = network["Psi_SApo"];  //# The reference is at current time step for other modes (implicit, explicit)
  if (psiref < Psi_SApo_cav) {
    network["Psi_SApo_cav"] = psiref;
    network["PLC_Stem"] = PLC(psiref, VCstem_slope, VCstem_P50);
  }
}

void innerCochard(List x, List input, List output, int n, double tstep, 
                 bool verbose = false, bool modifyInput = true) {
  
  // Extract hydraulic networks
  List networks = input["networks"];
  
  IntegerVector nsmalltimesteps = IntegerVector::create(2,4, 8, 16);
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
  int nlayers = dVec.length();
  
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
  NumericVector StemPsiVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["StemPsi"]);
  NumericVector LeafPsiVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["LeafPsi"]);
  NumericVector RootCrownPsiVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["RootCrownPsi"]);
  NumericVector StemSympPsiVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["StemSympPsi"]);
  
  // //Water pools
  // DataFrame belowdf = Rcpp::as<Rcpp::DataFrame>(x["below"]);
  // List belowLayers = Rcpp::as<Rcpp::List>(x["belowLayers"]);
  // NumericMatrix Wpool = Rcpp::as<Rcpp::NumericMatrix>(belowLayers["Wpool"]);
  // List RHOP;
  // NumericVector poolProportions(numCohorts);
  // if(plantWaterPools) {
  //   RHOP = belowLayers["RHOP"];
  //   poolProportions = belowdf["poolProportions"];
  // }
  // NumericMatrix RhizoPsiMAT = Rcpp::as<Rcpp::NumericMatrix>(belowLayers["RhizoPsi"]);
  
  //Extract output to be filled
  
  NumericMatrix SoilWaterExtract = Rcpp::as<Rcpp::NumericMatrix>(output["Extraction"]);
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
  List PlantsInst = output["PlantsInst"];
  NumericMatrix Einst = Rcpp::as<Rcpp::NumericMatrix>(PlantsInst["E"]);
  NumericMatrix Aginst = Rcpp::as<Rcpp::NumericMatrix>(PlantsInst["Ag"]);
  NumericMatrix Aninst = Rcpp::as<Rcpp::NumericMatrix>(PlantsInst["An"]);
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
    // # A. LOOP ON THE IMPLICIT SOLVER IN PSI, trying different time steps until results are OK
    bool regulationWellComputed = false;
    bool cavitationWellComputed = false;
    List network = networks[c];
    
    List params = network["params"];
    double gmin_S = params["gmin_S"];
    double gCrown0 = params["gCrown0"];
    double gmin20 = params["gmin20"];
    double TPhase_gmin = params["TPhase_gmin"];
    double Q10_1_gmin = params["Q10_1_gmin"];
    double Q10_2_gmin = params["Q10_2_gmin"];
    double fTRBToLeaf = params["fTRBToLeaf"];
    
    NumericVector kSoil = network["k_Soil"];

    double LAI = network["LAI"];
    Rcout << "\n*** HOUR STEP " << n << " cohort " << c << "***\n";

    int nwhilecomp = 0;
    NumericVector fluxSoilToStemLargeTimeStep(kSoil.size(), 0.0);

    List network_n;
    
    double Gwdiff_SL, Gwdiff_SH;
    while ((!regulationWellComputed || !cavitationWellComputed) && (nwhilecomp<nsmalltimesteps.size())) { //# LOOP TO TRY DIFFERENT TIME STEPS
      network_n = clone(network); // # initial value of WBveg
    //   List WBsoil_n = clone(WBsoil); // # initial value of WBsoil
    //   
      regulationWellComputed = false;
      cavitationWellComputed = false;
      double deltaRegulMax = 1.0e-100;
      double deltaPLCMax = 1.0e-100;

      //Reset soil flux to zero
      for(int i=0;i < kSoil.size();i++) fluxSoilToStemLargeTimeStep[i] = 0.0;

      int nts = nsmalltimesteps[nwhilecomp];// # number of small time steps
      double dt = tstep / ((double) nts); //Determine number of seconds of small time steps
      Rcout<< " Attempt #" << nwhilecomp<<" nts "<< nts << " dt " << dt << "\n";
      for(int its = 1; its <= nts; its++) { //#INTERNAL LOOP ON SMALL TIME STEPS
        // double p = (((double) its ) - 0.5)/((double) nts);
        // List WBclim = interp_WBclim(WBclim_current, WBclim_next, p); // # climate at nph
        // double ETPr = WBveg_n["ETPr"];
        // double Tair_mean = WBclim["Tair_mean"];
        // double RHair = WBclim["RHair"];
        // compute_evaporationG(WBsoil_n, RHair, Tair_mean,
        //                      ((double) Nhours)/((double) nts), LAI,
        //                      ETPr, (double) params["K"]);
        // fluxEvaporationSoilLargeTimeStep = fluxEvaporationSoilLargeTimeStep + ((double) WBsoil_n["E_Soil3"])/((double) nts);

        // network_np1 = clone(network_n); // Clone WBveg object
    //     compute_transpiration(WBveg_np1, WBclim, Nhours, opt, stomatalRegFormulation);// # transpi with climate at nph
    
        //Current leaf water potential (same for sunlit and shade leaves)
        double Psi_LSym = network_n["Psi_LSym"];
        NumericVector regul_ini = regulFact(Psi_LSym, params);
        
        //Leaf temperature for sunlit and shade leaves
        double Elim_SL = network_n["Elim_SL"];
        double Elim_SH = network_n["Elim_SH"];
        double Elim = network_n["Elim"];
        if(NumericVector::is_na(Elim_SL)) Elim_SL = Elim * (LAI_SL[c]/LAI);
        if(NumericVector::is_na(Elim_SH)) Elim_SH = Elim * (LAI_SH[c]/LAI);
        Temp_SL(c,n) = leafTemperature2(SWR_SL(c,n)/LAI_SL[c], LWR_SL(c,n)/LAI_SL[c], 
                                   Tair[iLayerSunlit[c]], zWind[iLayerSunlit[c]], 
                                   Elim_SL,  LeafWidth[c]);
        Temp_SH(c,n) = leafTemperature2(SWR_SH(c,n)/LAI_SH[c], LWR_SH(c,n)/LAI_SH[c], 
                                   Tair[iLayerShade[c]], zWind[iLayerShade[c]], 
                                   Elim_SH,  LeafWidth[c]);
        
        //VPD
        double VPD_air = meteoland::utils_saturationVP(Tair[iLayerCohort[c]]) - VPair[iLayerCohort[c]];
        VPD_SL(c,n) = std::max(0.0,leafVapourPressure(Temp_SL(c,n), Psi_LSym) - VPair[iLayerSunlit[c]]);
        VPD_SH(c,n) = std::max(0.0,leafVapourPressure(Temp_SH(c,n), Psi_LSym) - VPair[iLayerShade[c]]);
        Rcout<< "  AirT "<< Tair[iLayerCohort[c]] << " LT_SL "<< Temp_SL(c,n)<< " LT_SH "<< Temp_SH(c,n)<<"\n";
        Rcout<< "  VPD_air "<< VPD_air << " VPD_SL "<< VPD_SL(c,n)<< " VPD_SH "<< VPD_SH(c,n)<<"\n";
        
        //gCR = g Crown
        double gCR = gCrown(gCrown0, zWind[iLayerCohort[c]]);
        //gBL = g Boundary Layer
        double gBL = 397.0*pow(zWind[iLayerCohort[c]]/(LeafWidth[c]*0.0072), 0.5); // mmol boundary layer conductance
        
        //# Leaf cuticular conductances and cuticular transpiration
        double gmin_SL = gmin(Temp_SL(c,n), gmin20, TPhase_gmin, Q10_1_gmin, Q10_2_gmin);
        double gmin_SH = gmin(Temp_SH(c,n), gmin20, TPhase_gmin, Q10_1_gmin, Q10_2_gmin);
        double Emin_L_SL = Emin(gmin_SL, gBL, gCR, VPD_SL(c,n), Patm);
        double Emin_L_SH = Emin(gmin_SH, gBL, gCR, VPD_SH(c,n), Patm);
        double Emin_L = ((Emin_L_SL*LAI_SL[c]) + (Emin_L_SH*LAI_SH[c]))/LAI; 
        network_n["Emin_L"] = Emin_L;
        
        //Compute stem cuticular transpiration
        double Emin_S = fTRBToLeaf * Emin(gmin_S, gBL, gCR, VPD_air, Patm);
        network_n["Emin_S"] =  Emin_S;
        Rcout<< "  Emin_S "<< Emin_S<<" Emin_L_SL "<< Emin_L_SL<<" Emin_L_SH "<< Emin_L_SH<<" Emin_L "<< Emin_L<<"\n";
        
        // Current stomatal regulation
        NumericVector regul = regulFact(Psi_LSym, params);
        double gs_SL = gsJarvis(params, PAR_SL(c,n), Temp_SL(c,n));
        double gs_SH = gsJarvis(params, PAR_SH(c,n), Temp_SH(c,n));
        Rcout<< "  PAR_SL "<< PAR_SL(c,n)<<"  gs_SL "<< gs_SL<<"  PAR_SH "<< PAR_SH(c,n)<<" gs_SH "<< gs_SH<<"\n";
        GSW_SL(c,n) = gs_SL * regul["regulFact"];
        GSW_SH(c,n) = gs_SH * regul["regulFact"];
        
        // Stomatal transpiration
        Gwdiff_SL = 1.0/(1.0/gCR + 1.0/GSW_SL(c,n) + 1.0/gBL); 
        Gwdiff_SH = 1.0/(1.0/gCR + 1.0/GSW_SH(c,n) + 1.0/gBL); 
        Elim_SL = Gwdiff_SL * VPD_SL(c,n)/Patm;
        Elim_SH = Gwdiff_SH * VPD_SL(c,n)/Patm;
        network_n["Elim_SL"] = Elim_SL;
        network_n["Elim_SH"] = Elim_SH;
        Elim = ((Elim_SL*LAI_SL[c]) + (Elim_SH*LAI_SH[c]))/LAI; 
        network_n["Elim"] = Elim;
        Rcout<< "  Elim_SL "<< Elim_SL<<"  Elim_SH "<< Elim_SH<<"  Elim "<< Elim<<"\n";
        
        //Add transpiration sources
        network_n["Einst"] = Elim + Emin_S + Emin_L;
        
        //Effects on water potentials and flows
        semi_implicit_integration(network_n, dt, opt);
        update_conductances(network_n);
        update_capacitances(network_n);

        // # QUANTITIES TO CHECK IF THE RESOLUTION IS OK
        // # 1. delta regulation between n and np1 (MIQUEL: Only Psi_LSym changes between the two calculations, params should be the same)
        deltaRegulMax = std::max(deltaRegulMax,std::abs(((double) regul["regulFact"]) - ((double) regul_ini["regulFact"])));
     
        // # 2. PLC at n and np1
        deltaPLCMax = std::max(deltaPLCMax, (double) network_n["PLC_Leaf"] - (double) network_n["PLC_Leaf"]);
        deltaPLCMax = std::max(deltaPLCMax, (double) network_n["PLC_Stem"] - (double) network_n["PLC_Stem"]);
        // network_n = network_np1; //# Update network_n

        // # 3. update of soil on small time step (done by FP in version 16)
        // NumericVector fluxSoilToStem_mm = WBveg_np1["fluxSoilToStem"];
        // double Psi_SApo = WBveg_np1["Psi_SApo"];
    //     NumericVector k_SoilToStem = WBveg["k_SoilToStem"]; //MIQUEL: Why WBveg here?
    //     NumericVector PsiSoil = WBsoil_n["PsiSoil"];
    //     for(int i=0;i < kSoil.size();i++) {
    //       double fluxSoilToStem_mmolm2s = k_SoilToStem[i]*(PsiSoil[i] - Psi_SApo);
    //       fluxSoilToStemLargeTimeStep[i] = fluxSoilToStemLargeTimeStep[i] + fluxSoilToStem_mmolm2s/((double) nts);// # mean flux over one large time step
    //       fluxSoilToStem_mm[i] = convertFluxFrom_mmolm2s_To_mm(fluxSoilToStem_mmolm2s, ((double) Nhours)/((double) nts), LAI); // # Quantity from each soil layer to the below part
    //     }
    //     // # NB the time step for fluxSoilToStem_mm is Nhours/nts!
    //     update_soilWater(WBsoil_n, fluxSoilToStem_mm);
      } //# end loop small time step

      // # TESTS ON RESOLUTION
      network_n["Diag_deltaRegulMax"] = deltaRegulMax;
      regulationWellComputed = (deltaRegulMax<0.05);
      network_n["Diag_deltaPLCMax"] = deltaPLCMax;
      cavitationWellComputed = (deltaPLCMax<1.0);// # 1%
      network_n["Diag_timeStepInSeconds"] = dt;
      nwhilecomp = nwhilecomp + 1;
    } //# end while

    // # B. SAVING SOLUTION AT NEXT TIME STEP IN WBveg
    networks[c] = network_n;
    network = network_n;
    
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
    
    //Store leaf values (final substep)
    E_SL(c,n) = network["Elim_SL"];
    E_SH(c,n) = network["Elim_SH"];
    Psi_SH(c,n) = network["Psi_LSym"];
    Psi_SL(c,n) = network["Psi_LSym"];
    
    //Sunlit/shade photosynthesis
    NumericVector LP_SL = leafphotosynthesis(irradianceToPhotonFlux(PAR_SL(c,n))/LAI_SL[c], 
                                             Cair[iLayerSunlit[c]], Gwdiff_SL/1.6, 
                                             std::max(0.0,Temp_SL(c,n)), 
                                             Vmax298SL[c]/LAI_SL[c], Jmax298SL[c]/LAI_SL[c]);
    NumericVector LP_SH = leafphotosynthesis(irradianceToPhotonFlux(PAR_SH(c,n))/LAI_SH[c], 
                                             Cair[iLayerShade[c]], Gwdiff_SH/1.6, 
                                             std::max(0.0,Temp_SH(c,n)), 
                                             Vmax298SH[c]/LAI_SH[c], Jmax298SH[c]/LAI_SH[c]);
    Ci_SL(c,n) = LP_SL[0];
    Ci_SH(c,n) = LP_SH[0];
    Ag_SL(c,n) = LP_SL[1];
    Ag_SH(c,n) = LP_SH[1];
    An_SL(c,n) = Ag_SL(c,n) - 0.015*VmaxTemp(Vmax298SL[c]/LAI_SL[c], Temp_SL(c,n));
    An_SH(c,n) = Ag_SH(c,n) - 0.015*VmaxTemp(Vmax298SH[c]/LAI_SH[c], Temp_SH(c,n));
    
    //Store state
    StemPsiVEC[c] = network["Psi_LApo"];
    LeafPsiVEC[c] = network["Psi_LSym"];
    RootCrownPsiVEC[c] = network["Psi_SApo"];
    StemPLCVEC[c] = ((double) network["PLC_Stem"])/100.0;
    LeafPLCVEC[c] = ((double) network["PLC_Leaf"])/100.0;
    EinstVEC[c] = network["Einst"];
    ElimVEC[c] = network["Elim"];
    Emin_LVEC[c] = network["Emin_L"];
    Emin_SVEC[c] = network["Emin_S"];
    
    // Rcout<<iPMSunlit[c]<<" "<<iPMShade[c] <<" "<<GwSunlit[iPMSunlit[c]]<<" "<<GwShade[iPMShade[c]]<<" "<<fittedE[iPMSunlit[c]]<<" "<<fittedE[iPMShade[c]]<<"\n";
    //Get leaf status


    
    //Scale photosynthesis
    double Agsum = Ag_SL(c,n)*LAI_SL[c] + Ag_SH(c,n)*LAI_SH[c];
    double Ansum = An_SL(c,n)*LAI_SL[c] + An_SH(c,n)*LAI_SH[c];
    Aginst(c,n) = (1e-6)*12.01017*Agsum*tstep;
    Aninst(c,n) = (1e-6)*12.01017*Ansum*tstep;

    //Scale from instantaneous flow to water volume in the time step
    Einst(c,n) = EinstVEC[c]*0.001*0.01802*LAIphe[c]*tstep;
    

    // 
    // NumericVector Esoilcn(nlayerscon[c],0.0);
    // NumericVector ElayersVEC(nlayerscon[c],0.0);
    
    
    // //Get info from sFunctionBelow (this will be different depending on wether capacitance is considered)
    // NumericMatrix ERhizo = Rcpp::as<Rcpp::NumericMatrix>(sFunctionBelow["ERhizo"]);
    // NumericMatrix RhizoPsi = Rcpp::as<Rcpp::NumericMatrix>(sFunctionBelow["psiRhizo"]);
    
    //Store steady state stem and rootcrown and root surface water potential values
    // NumericMatrix newStemPsi = Rcpp::as<Rcpp::NumericMatrix>(sFunctionAbove["psiStem"]);
    // Stem1PsiVEC[c] = newStemPsi(iPM,0); 
    // Stem2PsiVEC[c] = newStemPsi(iPM,1);
    // for(int lc=0;lc<nlayerscon[c];lc++) {
    //   ElayersVEC[lc] = ERhizo(iPM,lc)*tstep; //Scale according to the time step
    // }
    
    //Copy RhizoPsi and from connected layers to RhizoPsi from soil layers
    // copyRhizoPsi(c,iPM, 
    //              RhizoPsi, RhizoPsiMAT,
    //              layerConnected, 
    //              RHOP, layerConnectedPools,
    //              VCroot_c, VCroot_d,  
    //              plantWaterPools);
    
    // Store the PLC corresponding to stem1 water potential
    // if(cavitationRefill!="total") {
    //   StemPLCVEC[c] = std::max(StemPLCVEC[c], 1.0 - xylemConductance(Stem1PsiVEC[c], 1.0, VCstem_c[c], VCstem_d[c])); 
    // } else { //Immediate refilling
    //   StemPLCVEC[c] = 1.0 - xylemConductance(Stem1PsiVEC[c], 1.0, VCstem_c[c], VCstem_d[c]); 
    // }
    
    //Scale soil water extracted from leaf to cohort level
    // for(int lc=0;lc<nlayerscon[c];lc++) {
    //   Esoilcn[lc] = ElayersVEC[lc]*0.001*0.01802*LAIphe[c]; //Scale from flow to water volume in the time step
    // }
    
    //Balance between extraction and transpiration
    // PWBinst(c,n) = sum(Esoilcn) - Einst(c,n);
    
    //Add step transpiration to daily plant cohort transpiration
    Eplant[c] += Einst(c,n);
    Anplant[c] += Aninst(c,n);
    Agplant[c] += Aginst(c,n);
    //Add PWB
    // PWB[c] += PWBinst(c,n); 
    
    
    
    //Copy transpiration and from connected layers to transpiration from soil layers
    //And update soil water content (soil water potential will not be updated until next day!)
    // if(!plantWaterPools) {
    //   int cl = 0;
    //   for(int l=0;l<nlayers;l++) {
    //     if(layerConnected(c,l)) {
    //       SoilWaterExtract(c,l) += Esoilcn[cl]; //Add to cummulative transpiration from layers
    //       soilLayerExtractInst(l,n) += Esoilcn[cl];
    //       //Apply extraction to soil layer
    //       if(modifyInput) Ws[l] = std::max(Ws[l] - (Esoilcn[cl]/Water_FC[l]),0.0);
    //       cl++;
    //     } 
    //   }
    // } else {
    //   NumericMatrix RHOPcoh = Rcpp::as<Rcpp::NumericMatrix>(RHOP[c]);
    //   LogicalMatrix layerConnectedCoh = Rcpp::as<Rcpp::LogicalMatrix>(layerConnectedPools[c]);
    //   int cl = 0;
    //   for(int j = 0;j<numCohorts;j++) {
    //     for(int l=0;l<nlayers;l++) {
    //       if(layerConnectedCoh(j,l)) {
    //         SoilWaterExtract(c,l) += Esoilcn[cl]; //Add to cummulative transpiration from layers
    //         soilLayerExtractInst(l,n) += Esoilcn[cl];
    //         //Apply extraction to soil layer
    //         if(modifyInput) Wpool(j,l) = Wpool(j,l) - (Esoilcn[cl]/(Water_FC[l]*poolProportions[j])); //Apply extraction from pools
    //         cl++;
    //       }
    //     }
    //   }
    // }
    if(N[c]>0.0) {
      //Store (for output) instantaneous leaf, stem and root potential, plc and rwc values
      PLC(c,n) = StemPLCVEC[c];
      StemSympRWCInst(c,n) = symplasticRelativeWaterContent(StemSympPsiVEC[c], StemPI0[c], StemEPS[c]);
      LeafSympRWCInst(c,n) = symplasticRelativeWaterContent(LeafPsiVEC[c], LeafPI0[c], LeafEPS[c]);
      StemRWCInst(c,n) = StemSympRWCInst(c,n)*(1.0 - StemAF[c]) + (1.0 - StemPLCVEC[c])*StemAF[c];
      LeafRWCInst(c,n) = LeafSympRWCInst(c,n)*(1.0 - LeafAF[c]) + (1.0 - LeafPLCVEC[c])*LeafAF[c];
      StemPsiInst(c,n) = StemPsiVEC[c]; 
      LeafPsiInst(c,n) = LeafPsiVEC[c]; //Store instantaneous (average) leaf potential
      RootPsiInst(c,n) = RootCrownPsiVEC[c]; //Store instantaneous root crown potential
      StemSympPsiInst(c,n) = StemSympPsiVEC[c];
      LeafSympPsiInst(c,n) = LeafPsiVEC[c];
      
      //Store the minimum water potential of the day (i.e. mid-day)
      minGSW_SL[c] = std::min(minGSW_SL[c], GSW_SL(c,n));
      minGSW_SH[c] = std::min(minGSW_SH[c], GSW_SH(c,n));
      maxGSW_SL[c] = std::max(maxGSW_SL[c], GSW_SL(c,n));
      maxGSW_SH[c] = std::max(maxGSW_SH[c], GSW_SH(c,n));
      minTemp_SL[c] = std::min(minTemp_SL[c], Temp_SL(c,n));
      minTemp_SH[c] = std::min(minTemp_SH[c], Temp_SH(c,n));
      maxTemp_SL[c] = std::max(maxTemp_SL[c], Temp_SL(c,n));
      maxTemp_SH[c] = std::max(maxTemp_SH[c], Temp_SH(c,n));
      minLeafPsi_SL[c] = std::min(minLeafPsi_SL[c], Psi_SL(c,n));
      minLeafPsi_SH[c] = std::min(minLeafPsi_SH[c], Psi_SH(c,n));
      maxLeafPsi_SL[c] = std::max(maxLeafPsi_SL[c], Psi_SL(c,n));
      maxLeafPsi_SH[c] = std::max(maxLeafPsi_SH[c], Psi_SH(c,n));
      minLeafPsi[c] = std::min(minLeafPsi[c], LeafPsiInst(c,n));
      maxLeafPsi[c] = std::max(maxLeafPsi[c], LeafPsiInst(c,n));
      minStemPsi[c] = std::min(minStemPsi[c], StemPsiInst(c,n));
      minRootPsi[c] = std::min(minRootPsi[c], RootPsiInst(c,n));
      // for(int l=0;l<nlayers;l++) {
      //   minPsiRhizo(c,l) = std::min(minPsiRhizo(c,l), RhizoPsiMAT(c,l));
      // }
    }
  }
}
