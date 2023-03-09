#include "tissuemoisture.h"
#include <Rcpp.h>
using namespace Rcpp;

double convertFluxFrom_mmolm2s_To_mm(double x, double Nhours, double LAI=1.0) {
  return(x * (LAI * Nhours * 3600.0 * 18.0) / 1.0e6);
}
double temp2SVP(double TD) {
  if(!NumericVector::is_na(TD)) return(0.61078*exp((17.269*TD)/(237.3+TD)));
  return(NA_REAL);
}
double vapourPressureFromRH(double T, double RH) {
  return(temp2SVP(T)*(RH/100.0));
}

//#interpolate climate values between two WBclim objects
List interp_WBclim(List clim1, List clim2, double p=0.5) {
  List res =  clone(clim1);
  res["Tair_mean"] = (1.0-p)*((double) clim1["Tair_mean"]) + p*((double) clim2["Tair_mean"]);
  res["RG"] = (1.0-p)*((double) clim1["RG"]) + p*((double) clim2["RG"]);
  res["WS"] = (1.0-p)*((double) clim1["WS"]) + p*((double) clim2["WS"]);
  res["VPD"] = (1.0-p)*((double) clim1["VPD"]) + p*((double) clim2["VPD"]);
  res["RHair_mean"] = (1.0-p)*((double) clim1["RHair_mean"]) + p*((double) clim2["RHair_mean"]);
  res["ETP"] = (1.0-p)*((double) clim1["ETP"]) + p*((double) clim2["ETP"]);
  return(res);
}

void update_soilConductanceAndPsi(List WBsoil) {
  List params = WBsoil["params"];
  NumericVector soilWaterStock = WBsoil["soilWaterStock"];
  NumericVector REW = WBsoil["REW"];
  NumericVector PsiSoil = WBsoil["PsiSoil"];
  NumericVector kSoil = WBsoil["kSoil"];
  
  NumericVector V_field_capacity = params["V_field_capacity"];
  NumericVector V_saturation_capacity_vg = params["V_saturation_capacity_vg"];
  NumericVector V_residual_capacity_vg = params["V_residual_capacity_vg"];
  NumericVector m = params["m"];
  NumericVector n = params["n"];
  NumericVector alpha_vg = params["alpha_vg"];
  NumericVector Ksat_vg = params["Ksat_vg"];
  NumericVector B_GC = params["B_GC"];
  NumericVector I_vg = params["I_vg"];
  double offSetPsoil = params["offSetPsoil"];
  
  //# Compute soil hydraulic conductivity with Van Genuchten
  //# Soil water holding capacity  (volumetric)
  NumericVector totalavailwater = V_saturation_capacity_vg - V_residual_capacity_vg;
  NumericVector totalavailwater_FC = (V_field_capacity - V_residual_capacity_vg);
  // Compute the relative water content (m3 water/m3 soil) based on Water Reserve and soil volume
  NumericVector actualavailwater = (soilWaterStock - V_residual_capacity_vg);
  for(int i = 0;i<REW.size();i++) {
    REW[i] = actualavailwater[i] / totalavailwater[i]; //  #/ numeric [ncouches
    if(REW[i] <=0.0001) REW[i] = 0.0001;
    else if(REW[i]>1.0) REW[i] = 1.0;
    // KSoil_temp <- REW^(WBsoil$params$I_vg) * (1 - (1 - REW^(1 / WBsoil$params$m))^WBsoil$params$m)^2
    double KSoil_temp = pow(REW[i], I_vg[i]) * pow(1.0 - pow(1.0 - pow(REW[i],(1.0 / m[i])),m[i]),2.0);
    // PsiSoil <- (-1 * ((((1 / REW)^(1 / WBsoil$params$m)) - 1)^(1 / WBsoil$params$n)) / WBsoil$params$alpha_vg / 10000) - WBsoil$params$offSetPsoil  # diviser par 10000 pour passer de cm à MPa
    PsiSoil[i] = (-1.0 * (pow((pow(1.0 / REW[i], 1.0 / m[i])) - 1.0, 1.0 / n[i])) / (alpha_vg[i] * 10000.0)) - offSetPsoil;//  # diviser par 10000 pour passer de cm à MPa
    // # Compute Soil conductance
    kSoil[i] = 1000.0 * Ksat_vg[i] * B_GC[i] * KSoil_temp;
  }

  double REW_tot_fc = sum(actualavailwater) / sum(totalavailwater_FC); // #/ numeric sum of all layers
  WBsoil["REW_tot_fc"] = std::min(REW_tot_fc,1.0);
  WBsoil["REW_tot_sat"] = sum(actualavailwater) / sum(totalavailwater); //  #/ numeric sum of all layers
}

// # Update soil water reservoirs, potential and psi according to wate fluxes from the plant
void update_soilWater(List WBsoil, NumericVector fluxEvap) {
  NumericVector soilWaterStock = WBsoil["soilWaterStock"];
  for(int i=0;i<soilWaterStock.size();i++) soilWaterStock[i] = soilWaterStock[i] - fluxEvap[i];
  update_soilConductanceAndPsi(WBsoil);
}


void update_kplant(List WBveg, List WBsoil) {
  List params = as<Rcpp::List>(WBveg["params"]);
  
  NumericVector k_RSApoInit = params["k_RSApoInit"];
  
  NumericVector k_RSApo = WBveg["k_RSApo"];
  NumericVector k_SoilToStem = WBveg["k_SoilToStem"];
  NumericVector kSoil = WBsoil["kSoil"];
  
  WBveg["k_SLApo"] = ((double) params["k_SLApoInit"]) * (1.0 - ((double) WBveg["PLC_Leaf"])/100.0);
  
  for(int i = 0;i<k_RSApo.size();i++) {
    //# calculate k_RSApo and k_SLApo with cavitation
    k_RSApo[i] = k_RSApoInit[i] * (1.0 - ((double) WBveg["PLC_Stem"])/100.0);
    //# Root from root length
    k_SoilToStem[i] = 1.0/((1.0/kSoil[i]) + (1.0/k_RSApo[i])); // # conductance from soil to collar (two resistances in series Rsoil and Rroot)
  }
    
  // Compute k_plant (from root to leaf) for diagnostic only
  WBveg["k_Plant"] =  1.0/ (1.0 /sum(k_RSApo) + 1.0/((double) WBveg["k_SLApo"]) + 1.0/((double) WBveg["k_LSym"]));
        
}


//# compute Evaporation from ETP and Gsoil and update SWS, Psi and in each soil layer
void compute_evaporationG(List WBsoil, double RHair, double Tair, 
                          double Nhours, double LAI, 
                          double ETP, double K) {
  //# created 03/01/2021 by JR / based on SurEau.C with gsoil0
  //# such as Esoil = gSoil0 * REW1 * VPDsoil/Patm

  List params = WBsoil["params"];
  NumericVector REW = WBsoil["REW"];
  NumericVector Evaporation = WBsoil["Evaporation"];
  NumericVector soilWaterStock = WBsoil["soilWaterStock"];
  
  //#TODO: improve this relation....
  double Tsoil = 0.6009*Tair+3.59;// # from relation fitted on O3HP
  
  double VPDsoil = vapourPressureFromRH(Tsoil, RHair);
  
  double evaporation = 0.0;
  if (Tsoil > 0.0) { //# no evaporation from frozen soil
    double g_Soil = ((double) params["gSoil0"]) * ((double) REW[0]);
    double E_Soil1 = g_Soil * VPDsoil / 101.3; // #  VPD effect
    double ETP_mmol_s = 1.0e6 * ETP / (3600.0 * Nhours * 18.0);
    double E_Soil2 = (g_Soil / ((double) params["gSoil0"])) * ETP_mmol_s * exp(-1.0*K * LAI);// # limitation by ETP depending on radiation reaching the soil
    double E_Soil3 = std::min(E_Soil1, E_Soil2);
    evaporation = convertFluxFrom_mmolm2s_To_mm(E_Soil3, Nhours); // # Conversion from mmol/m2/s to mm
  }
    
  WBsoil["EvaporationSum"] = evaporation;
  WBsoil["Evaporation"] = evaporation;
  
  soilWaterStock[0] = soilWaterStock[0] - evaporation;
  update_soilConductanceAndPsi(WBsoil);
}

double turgor_comp(double PiFT, double Esymp, double Rstemp) {
  return(-1.0*PiFT - Esymp * Rstemp);
}

List compute_regulFact(double psi, List params, String regulationType) {
  
  
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
   //# regulFactPrime = al*exp(-al*(psi-P50_gs)) / ((1+exp(-al*(psi-P50_gs)))^2) # formulation in psi (same but slower to compute)
  } else if (regulationType == "Turgor") {
    double turgorPressureAtGsMax = params["turgorPressureAtGsMax"];
    double epsilonSym_Leaf = params["epsilonSym_Leaf"];
    double PiFullTurgor_Leaf = params["PiFullTurgor_Leaf"];
    double rs = symplasticRelativeWaterContent(psi, PiFullTurgor_Leaf, epsilonSym_Leaf);
    double turgor = turgor_comp(PiFullTurgor_Leaf, epsilonSym_Leaf, rs);
    regulFact = std::max(0.0, std::min(1.0, turgor / turgorPressureAtGsMax));
    if((regulFact == 1.0) || (regulFact == 0.0)) {
      regulFactPrime = 0.0;
    } else {
      regulFactPrime = 0.0;// # TODO insert the derivative to compute regulFactPrime
    }
  }
  List res = List::create(_["regulFact"] = regulFact, 
                          _["stomatalClosure"] = stomatalClosure, 
                          _["regulFactPrime"] = regulFactPrime);
  return(res);
}


double PLCPrime_comp(double plc, double slope) {
  return(-1.0*slope/25.0 * plc/100 * (1.0 - plc/100));
}
double PLC_comp(double Pmin, double slope, double P50){
  return (100.0 / (1.0 + exp(slope / 25.0 * (Pmin - P50))));
}

void semi_implicit_temporal_integration(List WBveg, List WBsoil, double dt, int nsmalltimesteps, NumericVector opt) {
  
  List params = as<Rcpp::List>(WBveg["params"]);
  
  // Step 1. Initializing current time step according to computation options (FP)
  double dbxmin = 1e-100; // FP minimal double to avoid 0/0
  double Psi_LApo_n = WBveg["Psi_LApo"];
  double Psi_SApo_n = WBveg["Psi_SApo"];
  double Psi_LSym_n = WBveg["Psi_LSym"];
  double Psi_SSym_n = WBveg["Psi_SSym"];
  double Psi_LApo_cav = WBveg["Psi_LApo_cav"];
  double Psi_SApo_cav = WBveg["Psi_SApo_cav"];
  
  //Modifiers
  double Lsym = opt["Lsym"];
  double Ssym = opt["Ssym"];
  double CLapo = opt["CLapo"];
  double CTapo = opt["CTapo"];
  double Eord = opt["Eord"];
  double Lcav = opt["Lcav"];
  double Scav = opt["Scav"];
  
  double k_SSym = WBveg["k_SSym"];
  double k_LSym = WBveg["k_LSym"];
  double c_LSym = WBveg["C_LSym"];
  double c_SSym = WBveg["C_SSym"];
  double c_LApo = WBveg["C_LApo"];
  double c_SApo = WBveg["C_SApo"];
  
  double K_LSym = Lsym * k_LSym;   //
  double K_SSym = Ssym * k_SSym;   //
  double C_LSym = Lsym * c_LSym;   //
  double C_SSym = Ssym * c_SSym;   //
  double C_LApo = CLapo * c_LApo; //
  double C_SApo = CTapo * c_SApo; //
  
  
  double K_SL = WBveg["k_SLApo"];
  
  double E_nph = WBveg["Elim"]; // COMPUTED WITH WBveg$Psi_LSym AND INTERPOLATE CLIMATE
  double Eprime = WBveg["Eprime"];
  double Eprime_nph = Eord * Eprime;
  
  double Emin_L_nph = WBveg["Emin"];
  double Emin_S_nph = WBveg["Emin_S"];
  
  double PLC_Leaf = WBveg["PLC_Leaf"];
  double PLC_Stem = WBveg["PLC_Stem"];
  double Q_LApo_sat_mmol_perLeafArea = WBveg["Q_LApo_sat_mmol_perLeafArea"];
  double Q_SApo_sat_mmol_perLeafArea = WBveg["Q_SApo_sat_mmol_perLeafArea"];
  
  double slope_VC_Leaf = params["slope_VC_Leaf"];
  double slope_VC_Stem = params["slope_VC_Stem"];
  double P50_VC_Leaf = params["P50_VC_Leaf"];
  double P50_VC_Stem = params["P50_VC_Stem"];
  
  //Compute K_L_Cav et K_S_Cav
  double PLC_prime_L = PLCPrime_comp(PLC_Leaf, slope_VC_Leaf);
  double K_L_Cav = -1.0 * Lcav * Q_LApo_sat_mmol_perLeafArea * PLC_prime_L / dt;  // avec WBveg$Q_LSym_sat en l/m2 sol # changed by NM (25/10/2021)
  double PLC_prime_S = PLCPrime_comp(PLC_Stem, slope_VC_Stem);
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
    
  NumericVector k_SoilToStem = WBveg["k_SoilToStem"];
  NumericVector PsiSoil = WBsoil["PsiSoil"];

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
  
  WBveg["Diag_nwhile_cavit"] = nwhilecomp;  // # Diagnostic step to track cavit event and eventual errors (corresponding to nwhilecomp==5)

  //# Step 3. Compute Psi_Symp_np1 (L and S)
  alpha = exp(-1.0*K_LSym/C_LSym*dt);
  Psi_td = (K_LSym*Psi_LApo_n - (E_nph + Emin_L_nph))/(K_LSym + dbxmin);// # dbxmin to avoid 0/0
  Psi_LSym_np1 = alpha * Psi_LSym_n +(1.0 - alpha) * Psi_td;
  alpha = exp(-1.0*K_SSym/C_SSym*dt);
  Psi_td = (K_SSym*Psi_SApo_n - Emin_S_nph)/(K_SSym + dbxmin); // # dbxmin to avoid 0/0
  Psi_SSym_np1 = alpha * Psi_SSym_n +(1.0 - alpha) * Psi_td;

  //#Step 4 : set computed values in WBveg and update Psi_cav, PLC and Psi_AllSoil
  WBveg["Psi_LApo"] = std::min(-0.00001, Psi_LApo_np1);
  WBveg["Psi_SApo"] = std::min(-0.00001,Psi_SApo_np1);
  WBveg["Psi_LSym"] = std::min(-0.00001,Psi_LSym_np1);
  WBveg["Psi_SSym"] = std::min(-0.00001,Psi_SSym_np1);

  //# Cavitation
  psiref = WBveg["Psi_LApo"];  //# the reference is at current time step for other modes  (implicit, explicit)
  if(psiref < Psi_LApo_cav) {
    WBveg["Psi_LApo_cav"] = psiref;
    WBveg["PLC_Leaf"] = PLC_comp(psiref, slope_VC_Leaf, P50_VC_Leaf);
  }

  psiref = WBveg["Psi_SApo"];  //# The reference is at current time step for other modes (implicit, explicit)
  if (psiref < Psi_SApo_cav) {
    WBveg["Psi_SApo_cav"] = psiref;
    WBveg["PLC_Stem"] = PLC_comp(psiref, slope_VC_Stem, P50_VC_Stem);
  }


  WBveg["Psi_AllSoil"] = sum(k_SoilToStem * PsiSoil)/sum(k_SoilToStem);

}

// [[Rcpp::export]]
List compute_plantNextTimeStep(List WBveg, List WBsoil, List WBclim_current, List WBclim_next, 
                               int Nhours, 
                               IntegerVector nsmalltimesteps, NumericVector opt) {
  
  // # A. LOOP ON THE IMPLICIT SOLVER IN PSI, trying different time steps until results are OK
  bool regulationWellComputed = false;
  bool cavitationWellComputed = false;
  
  int nwhilecomp = 0;
  NumericVector fluxSoilToStemLargeTimeStep;
  while ((!regulationWellComputed || !cavitationWellComputed) && (nwhilecomp<nsmalltimesteps.size())) { //# LOOP TO TRY DIFFERENT TIME STEPS
    List WBveg_n = clone(WBveg); // # initial value of WBveg
    List WBsoil_n = clone(WBsoil); // # initial value of WBsoil
     
    double LAI = (double) WBveg["LAI"];
     
    List params  = WBveg_n["params"];
    
    bool regulationWellComputed = false;
    bool cavitationWellComputed = false;
    int nwhilecomp = nwhilecomp + 1;
    double deltaRegulMax = 1.0e-100;
    double deltaPLCMax = 1.0e-100;
    
    int nts = nsmalltimesteps[nwhilecomp];// # number of small time steps
    double fluxEvaporationSoilLargeTimeStep = 0.0;
    for(int its = 1; its <= nts; its++) { //#INTERNAL LOOP ON SMALL TIME STEPS
      double p = (((double) its ) - 0.5)/((double) nts);
      List WBclim = interp_WBclim(WBclim_current, WBclim_next, p); // # climate at nph
      double ETPr = WBveg_n["ETPr"];
      double Tair_mean = WBclim["Tair_mean"];
      double RHair = WBclim["RHair"];
      compute_evaporationG(WBsoil_n, RHair, Tair_mean, 
                           ((double) Nhours)/((double) nts), (double) WBveg_n["LAI"],
                           ETPr, (double) params["K"]);
      double fluxEvaporationSoilLargeTimeStep = fluxEvaporationSoilLargeTimeStep + ((double) WBsoil_n["E_Soil3"])/((double) nts);
      
      List WBveg_np1 = clone(WBveg_n);
      compute_transpiration(WBveg_np1, WBclim, Nhours, modeling_options);// # transpi with climate at nph
      semi_implicit_temporal_integration(Wveg_np1,  WBsoil_n, dt = Nhours * 3600 / nts, opt = opt);
      update_kplant(WBveg_np1,WBsoil_n);
      update_capacitancesApoAndSym(WBveg_np1);

      // # QUANTITIES TO CHECK IF THE RESOLUTION IS OK
      // # 1. delta regulation between n and np1
      List regul_np1 = compute_regulFact(psi = WBveg_np1$Psi_LSym, params = WBveg_np1$params,regulationType=modeling_options$stomatalRegFormulation);
      List regul_n   = compute_regulFact(psi = WBveg_n$Psi_LSym  , params = WBveg_n$params  ,regulationType=modeling_options$stomatalRegFormulation); //# TODO check why recomputed? should be in WBveg_tmp
      
      deltaRegulMax = std::max(deltaRegulMax,std::abs(((double) regul_np1["regulFact"]) - ((double) regul_n["regulFact"])));
      
      // # 2. PLC at n and np1
      deltaPLCMax = max(deltaPLCMax,WBveg_np1$PLC_Leaf-WBveg_n$PLC_Leaf,WBveg_np1$PLC_Stem-WBveg_n$PLC_Stem)
      WBveg_n = WBveg_np1; //# Update WBveg_n

      // # 3. update of soil on small time step (done by FP in version 16)
      NumericVector fluxSoilToStem_mm = WBveg_np1["fluxSoilToStem"];
      double Psi_SApo = WBveg_np1["Psi_SApo"];
      NumericVector k_SoilToStem = WBveg["k_SoilToStem"]; //MIQUEL: Why WBveg here?
      NumericVector PsiSoil = WBsoil_n["PsiSoil"]; 
      for(int i=0;i < fluxSoilToStem_mm.size();i++) {
        double fluxSoilToStem_mmolm2s = k_SoilToStem[i]*(PsiSoil[i] - Psi_SApo);
        fluxSoilToStemLargeTimeStep[i] = fluxSoilToStemLargeTimeStep[i] + fluxSoilToStem_mmolm2s/((double) nts);// # mean flux over one large time step
        fluxSoilToStem_mm[i] = convertFluxFrom_mmolm2s_To_mm(fluxSoilToStem_mmolm2s, ((double) Nhours)/((double) nts), LAI); // # Quantity from each soil layer to the below part
      }
      // # NB the time step for fluxSoilToStem_mm is Nhours/nts!
      update_soilWater(WBsoil_n, fluxSoilToStem_mm);
    } //# end loop small time step
    
    // # TESTS ON RESOLUTION
    WBveg_np1["Diag_deltaRegulMax"] = deltaRegulMax;
    regulationWellComputed = (deltaRegulMax<0.05);
    WBveg_np1["Diag_deltaPLCMax"] = deltaPLCMax;
    cavitationWellComputed = (deltaPLCMax<1.0);// # 1%
    WBveg_np1["Diag_timeStepInHours"] = ((double) Nhours)/((double) nts);
  } //# end while
  
  // # B. SAVING SOLUTION AT NEXT TIME STEP IN WBveg
  WBveg = WBveg_np1;
  // # final update of transpiration at clim_next (useful for consistency in outputs, but not required for the computations)
  compute_transpiration(WBveg, WBclim_next, Nhours,modeling_options=modeling_options);
  // # C. UPDATING FLUX FROM SOIL (WBveg$fluxSoilToStem is used as input in UpdateSoilWater.WBsoil)
  // #TODO FP suggests moving the computation of  fluxSoilToStem in the main loop, as it is the coupling between the two models...
  // # mean soil quantities on large time steps
  WBveg["Emin_mm"]  = convertFluxFrom_mmolm2s_To_mm((double) WBveg["Emin"], (double) Nhours, LAI); // # Flux from each soil layer to the below part
  WBveg["Emin_S_mm"] = convertFluxFrom_mmolm2s_To_mm((double) WBveg["Emin_S"], (double) Nhours, LAI); // # Flux from each soil layer to the below part
     
  WBveg["SumFluxSoilToStem"] = sum(fluxSoilToStemLargeTimeStep);// # flux total en mmol/m2/s / used for Tleaf
  NumericVector fluxSoilToStem_mm = WBveg["fluxSoilToStem_mm"];
  for(int i=0;i< fluxSoilToStem_mm.size();i++) {
    fluxSoilToStem_mm[i] = convertFluxFrom_mmolm2s_To_mm(fluxSoilToStemLargeTimeStep[i], (double) Nhours, LAI); // # Flux from each soil layer to the below part  in mm
  }
  WBveg["transpiration_mm"] = convertFluxFrom_mmolm2s_To_mm((double)(WBveg["Emin"] + WBveg["Emin_S"] + WBveg["Elim"]), (double) Nhours, LAI); //# total flux in mm
     
  return(WBveg)
}