// [[Rcpp::interfaces(r,cpp)]]
#define STRICT_R_HEADERS
#include <meteoland.h>
#include <Rcpp.h>
using namespace Rcpp;

double convertFluxFrom_mmolm2s_To_mm(double x, double timeStep, double LAI=1) {
  return(x * (LAI * timeStep * 3600.0 * 18.0) / 1.0e6);
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

//# compute Evaporation from ETP and Gsoil and update SWS, Psi and in each soil layer
void compute_evaporationG(List WBsoil, double RHair, double Tair, int Nhours, double LAI, double ETP, double K) {
  //# created 03/01/2021 by JR / based on SurEau.C with gsoil0
  //# such as Esoil = gSoil0 * REW1 * VPDsoil/Patm

  List params = WBsoil["params"];
  NumericVector REW = WBsoil["REW"];
  NumericVector Evaporation = WBsoil["Evaporation"];
  NumericVector soilWaterStock = WBsoil["soilWaterStock"];
  
  //#TODO: improve this relation....
  double Tsoil = 0.6009*Tair+3.59;// # from relation fitted on O3HP
  
  // VPDsoil <- compute.VPDfromRHandT(RHair, Tsoil)
  double VPDsoil = meteoland::vapourPressureFromRH(Tsoil, RHair);
  
  if (Tsoil < 0.0) { //# no evaporation from frozen soil
    for(int i=0;i< Evaporation.size(); i++) Evaporation[i] = 0.0;
  } else {
    double g_Soil = ((double) params["gSoil0"]) * ((double) REW[0]);
    double E_Soil1 = g_Soil * VPDsoil / 101.3; // #  VPD effect
    double ETP_mmol_s = 1.0e6 * ETP / ((double) (3600 * Nhours * 18));
    double E_Soil2 = (g_Soil / ((double) params["gSoil0"])) * ETP_mmol_s * exp(-1.0*K * LAI);// # limitation by ETP depending on radiation reaching the soil
    double E_Soil3 = std::min(E_Soil1, E_Soil2);
    WBsoil["Evaporation"] = convertFluxFrom_mmolm2s_To_mm(E_Soil3, timeStep=Nhours); // # Conversion from mmol/m2/s to mm
    Evaporation = WBsoil["Evaporation"];
  }
    
  WBsoil["EvaporationSum"] = sum(Evaporation);
      
  soilWaterStock[0] = soilWaterStock[0] - WBsoil$Evaporation
    WBsoil <- compute.soilConductanceAndPsi.WBsoil(WBsoil)
      
}



double PLCPrime_comp(double plc, double slope) {
  return(-1.0*slope/25.0 * plc/100 * (1.0 - plc/100));
}
double PLC_comp(double Pmin, double slope, double P50){
  return (100.0 / (1.0 + exp(slope / 25.0 * (Pmin - P50))));
}

void SemiImplicitTemporalIntegration(List WBveg, List WBsoil, double dt, int nsmalltimesteps, NumericVector opt) {
  
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
void compute_plantNextTimeStep(List WBveg, List WBsoil, List WBclim_current, List WBclim_next, 
                               int Nhours, 
                               IntegerVector nsmalltimesteps, NumericVector opt) {
  
  // # A. LOOP ON THE IMPLICIT SOLVER IN PSI, trying different time steps until results are OK
  bool regulationWellComputed = false;
  bool cavitationWellComputed = false;
  
  int nwhilecomp = 0;
  while ((!regulationWellComputed || !cavitationWellComputed) && (nwhilecomp<nsmalltimesteps.size())) { //# LOOP TO TRY DIFFERENT TIME STEPS
    List WBveg_n = clone(WBveg); // # initial value of WBveg
    List WBsoil_n = clone(WBsoil); // # initial value of WBsoil
     
    bool regulationWellComputed = false;
    bool cavitationWellComputed = false;
    int nwhilecomp = nwhilecomp + 1;
    double deltaRegulMax = 1.0e-100;
    double deltaPLCMax = 1.0e-100;
    
    int nts = nsmalltimesteps[nwhilecomp];// # number of small time steps
    double fluxSoilToStemLargeTimeStep = 0.0;
    double fluxEvaporationSoilLargeTimeStep = 0.0;
    for(int its = 1; i <= nts; i++) { //#INTERNAL LOOP ON SMALL TIME STEPS
      double p = (((double) its ) - 0.5)/((double) nts);
      List WBclim = interp_WBclim(WBclim_current, WBclim_next, p); // # climate at nph
      List WBsoil_n = compute_evaporationG(WBsoil = WBsoil_n,ETP = WBveg_n$ETPr,Tair = WBclim$Tair_mean, RHair = WBclim$RHair,K = WBveg_n$params$K,LAI = WBveg_n$LAI,Nhours = Nhours/nts)
//       fluxEvaporationSoilLargeTimeStep = fluxEvaporationSoilLargeTimeStep + WBsoil_n$E_Soil3/nts
// #WBclim = lapply(seq_along(WBclim_current),function(i) unlist(0.5*(WBclim_current[i])+unlist(WBclim_next[i])))
//       WBveg_tmp <- compute.transpiration.WBveg(WBveg_n, WBclim, Nhours, modeling_options) # transpi with climate at nph
//       WBveg_np1 <- implicit.temporal.integration.atnp1(WBveg_tmp,  WBsoil_n, dt = Nhours * 3600 / nts, opt = opt)
//       WBveg_np1 <- update.kplant.WBveg(WBveg_np1,WBsoil_n)
//       WBveg_np1 <- update.capacitancesApoAndSym.WBveg(WBveg_np1)
//       
// #browser()
// # QUANTITIES TO CHECK IF THE RESOLUTION IS OK
// # 1. delta regulation between n and np1
//       regul_np1 = compute.regulFact(psi = WBveg_np1$Psi_LSym, params = WBveg_np1$params,regulationType=modeling_options$stomatalRegFormulation)
//         regul_n   = compute.regulFact(psi = WBveg_n$Psi_LSym  , params = WBveg_n$params  ,regulationType=modeling_options$stomatalRegFormulation)# TODO check why recomputed? should be in WBveg_tmp
//       deltaRegulMax = max(deltaRegulMax,abs(regul_np1$regulFact-regul_n$regulFact))
// # 2. PLC at n and np1
//         deltaPLCMax = max(deltaPLCMax,WBveg_np1$PLC_Leaf-WBveg_n$PLC_Leaf,WBveg_np1$PLC_Stem-WBveg_n$PLC_Stem)
//         WBveg_n = WBveg_np1 # Update WBveg_n
//       
// # 3. update of soil on small time step (done by FP in version 16)
//       fluxSoilToStem = WBveg$k_SoilToStem*(WBsoil_n$PsiSoil-WBveg_np1$Psi_SApo)
// # NB the time step for fluxSoilToStem is Nhours/nts!
//         WBveg_np1$fluxSoilToStem = convertFluxFrom_mmolm2s_To_mm(fluxSoilToStem, LAI = WBveg$LAI, timeStep = Nhours/nts) # Quantity from each soil layer to the below part
//       WBsoil_n <- update.soilWater.WBsoil(WBsoil = WBsoil_n, fluxEvap = WBveg_np1$fluxSoilToStem)
//         fluxSoilToStemLargeTimeStep = fluxSoilToStemLargeTimeStep + fluxSoilToStem/nts # mean flux over one large time step
//       
// # if (opt$numericalScheme == "Explicit" ) {
// #   write.WBoutput(Date = NA, WBoutput = WBoutput, WBsoil = WBsoil, WBveg = WBveg_n, WBclim = WBclim_next)
// # }
//       
//       
//     } # end loop small time step
// # TESTS ON RESOLUTION
//     WBveg_np1$Diag_deltaRegulMax = deltaRegulMax
//     regulationWellComputed = (deltaRegulMax<0.05)
//       WBveg_np1$Diag_deltaPLCMax = deltaPLCMax
//     cavitationWellComputed = (deltaPLCMax<1) # 1%
//       WBveg_np1$Diag_timeStepInHours = Nhours/nts
//     
// # if (nwhilecomp==length(opt$nsmalltimesteps)&deltaRegulMax>0.05) {
// #   warning(paste0('regulation inacurate(deltaRegulMax=',signif(deltaRegulMax, digits = 3),'; please reduce the time step, currently=',Nhours/nts))
// # }
// # if (nwhilecomp==length(opt$nsmalltimesteps)&deltaPLCMax>1) {# 1%
// #   warning(paste0('water release from cavitation inacurate(deltaPLCMax(%)=',signif(deltaPLCMax, digits = 3),'; please reduce the time step, currently=',Nhours/nts))
// # }
  } # end while
// # B. SAVING SOLUTION AT NEXT TIME STEP IN WBveg
//     WBveg = WBveg_np1
//   
//   WBveg <- compute.transpiration.WBveg(WBveg, WBclim_next, Nhours,modeling_options=modeling_options) # final update of transpiration at clim_next (useful for consistency in outputs, but not required for the computations)
//     
//     
// # C. UPDATING FLUX FROM SOIL (WBveg$fluxSoilToStem is used as input in UpdateSoilWater.WBsoil)
// #TODO FP suggests moving the computation of  fluxSoilToStem in the main loop, as it is the coupling between the two models...
//     
// # mean soil quantities on large time steps
//     WBveg$Emin_mm  = convertFluxFrom_mmolm2s_To_mm(WBveg$Emin, LAI = WBveg$LAI, timeStep = Nhours) # Flux from each soil layer to the below part
//     WBveg$Emin_S_mm = convertFluxFrom_mmolm2s_To_mm(WBveg$Emin_S, LAI = WBveg$LAI, timeStep = Nhours) # Flux from each soil layer to the below part
//     
//     WBveg$SumFluxSoilToStem <- sum(fluxSoilToStemLargeTimeStep) # flux total en mmol/m2/s / used for Tleaf
//     WBveg$fluxSoilToStem_mm  <- convertFluxFrom_mmolm2s_To_mm(fluxSoilToStemLargeTimeStep, LAI = WBveg$LAI, timeStep = Nhours) # Flux from each soil layer to the below part  in mm
//     WBveg$transpiration_mm     <- convertFluxFrom_mmolm2s_To_mm((WBveg$Emin + WBveg$Emin_S + WBveg$Elim),LAI = WBveg$LAI, timeStep = Nhours) # total flux in mm
//     
//     return(WBveg)
}