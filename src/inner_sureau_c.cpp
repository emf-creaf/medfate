#define STRICT_R_HEADERS
#include <RcppArmadillo.h>
#include "medfate.h"
#include "inner_sureau_c.h"
#include "biophysicsutils_c.h"
#include "hydraulics_c.h"
#include "modelInput_c.h"
#include "soil_c.h"
#include "tissuemoisture_c.h"
#include <meteoland.h>
using namespace Rcpp;


double PLC_derivative_c(double plc, double slope) {
  return(-1.0*slope/25.0 * plc/100 * (1.0 - plc/100));
}
double PLC_c(double Pmin, double slope, double P50){
  return (100.0 / (1.0 + exp(slope / 25.0 * (Pmin - P50))));
}
double invPLC_c(double plc, double slope, double P50){
  return(P50 + log((100.0/plc)-1.0)*(25.0/slope));
}

double RWC_c(double PiFT, double Esymp, double Pmin) {
  double A = std::max((-1.0 * (Pmin + PiFT - Esymp) - sqrt(pow(Pmin + PiFT - Esymp, 2.0) + 4.0 * (Pmin * Esymp))) / (2.0 * Esymp), 1.0 - PiFT / Pmin);
  return(A);
}

double Emin_c(double gmin, double gBL, double gCrown, 
            double VPD, double airPressure) {
  double gmintot = 1.0/(1.0/gmin+ 1.0/gBL + 1.0/gCrown);
  return(gmintot * VPD /airPressure); 
}
double Turgor_c(double PiFT, double Esymp, double Rstemp) {
  return(-1.0*PiFT - Esymp * Rstemp);
}


// # Update plant conductances
void update_conductances_c(SureauNetwork &network) {
  
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
void update_capacitances_c(SureauNetwork &network) {
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
  double PsiTLP_Leaf = turgorLossPoint_c(PiFullTurgor_Leaf, epsilonSym_Leaf);
  double PsiTLP_Stem = turgorLossPoint_c(PiFullTurgor_Stem, epsilonSym_Stem);
  
  //#----Compute the relative water content of the symplasm----
  double RWC_LSym = 1.0 - RWC_c(PiFullTurgor_Leaf, epsilonSym_Leaf, Psi_LSym - dbxmin);
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
  double RWC_SSym = 1.0 - RWC_c(PiFullTurgor_Stem, epsilonSym_Stem, Psi_SSym - dbxmin);
  
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


//# stomatal conductance calculation with Jarvis type formulations
double gsJarvis_c(SureauParams &params, double PAR, double Temp, int option){
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
  double Q = irradianceToPhotonFlux_c(PAR, defaultLambda); //From W m-2 to micromol s-1 m-2
  double gs_bound = gsNight2 + (gsMax2 - gsNight2) * (1.0 - exp(-1.0*JarvisPAR*Q));
  return(gs_bound);
}


void semi_implicit_integration_inner_c(SureauNetwork& network,
                                       double dt, 
                                       const SureauOpt& opt, 
                                       const std::string& stemCavitationRecovery, 
                                       const std::string& leafCavitationRecovery) {
  
  const SureauParams& params = network.params;
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
  

  //Apply modifiers
  double K_LSym = opt.Lsym * k_LSym;   
  double K_SSym = opt.Ssym * k_SSym;   
  double C_LSym = opt.Lsym * c_LSym;   
  double C_SSym = opt.Ssym * c_SSym;   
  double C_LApo = opt.CLapo * c_LApo; 
  double C_SApo = opt.CTapo * c_SApo; 
  
  double E_nph = network.Elim; // Leaf stomatal transpiration
  double Emin_L_nph = network.Emin_L; //Leaf cuticular transpiration
  double Emin_S_nph = network.Emin_S; //Stem cuticular transpiration
  
  
  
  
  double VCleaf_slope = params.VCleaf_slope;
  double VCstem_slope = params.VCstem_slope;
  double VCleaf_P50 = params.VCleaf_P50;
  double VCstem_P50 = params.VCstem_P50;
  
  //Compute K_L_Cav et K_S_Cav
  double PLC_prime_L = PLC_derivative_c(PLC_Leaf, VCleaf_slope);
  double K_L_Cav = -1.0 * opt.Lcav * Q_LApo_sat_mmol_perLeafArea * PLC_prime_L / dt;  // avec WBveg$Q_LSym_sat en l/m2 sol # changed by NM (25/10/2021)
  double PLC_prime_S = PLC_derivative_c(PLC_Stem, VCstem_slope);
  double K_S_Cav = -1.0 * opt.Scav * Q_SApo_sat_mmol_perLeafArea * PLC_prime_S / dt;  // opt$Scav * WBveg$K_S_Cav #FP corrected a bug sign herehanged by NM (25/10/2021)
  // Rcout<< "0 "<< PLC_prime_L << " "<<K_L_Cav<<" "<<PLC_prime_S<< " "<< K_S_Cav<<"\n";
  
  // Step 2. While loop in order to decide if cavitation or not :
  //  In order to account for the cavitation that occurs only when potentials go below their lowest value "cav" (formerly called "mem" in an earlier version)
  // the following computations are done trying sequentially the resolutions of LApo and TApo eventually activating
  // the appropriate cavitation events when needed (starting assuming no cavit at all...)
  // in case of computational problem, the last case assume no cavitation flux
  bool LcavitWellComputed = false; //initialized to false
  bool ScavitWellComputed = false;
  
  std::vector<double> delta_L_cavs, delta_S_cavs;
  if ((opt.Lcav==0.0) && (opt.Scav==0.0)) { // no cavitation flux computed
    delta_L_cavs = {0.0};
    delta_S_cavs = {0.0};
  } else if ((opt.Lcav==0.0) && (opt.Scav==1.0)) {// Scav only
    delta_L_cavs= {0.0,0.0,0.0};
    delta_S_cavs= {0.0,1.0,0.0};
  } else if ((opt.Lcav==1.0) && (opt.Scav==0.0)) {// Lcav only
    delta_L_cavs= {0.0,1.0,0.0};
    delta_S_cavs={0.0,0.0,0.0};
  } else { //#Lcav=1 and Scav=1
    delta_L_cavs= {0.0,1.0,0.0,1.0,0.0}; // the fifth case is here in case no solution with others...
    delta_S_cavs= {0.0,0.0,1.0,1.0,0.0};
  }
  
  double k_CSApo = network.k_CSApo;
  double* k_SoilToStem = network.k_SoilToStem;
  
  double alpha, Psi_td, Psi_LApo_np1, Psi_SApo_np1, Psi_LSym_np1, Psi_SSym_np1;
  double psirefL, psirefS;
  
  double Psi_RCApo_np1 = 0.0;
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
    LcavitWellComputed = (delta_L_cav==(Psi_LApo_np1 < Psi_LApo_cav)) || (opt.Lcav==0.0);
    ScavitWellComputed = (delta_S_cav==(Psi_SApo_np1 < Psi_SApo_cav)) || (opt.Scav==0.0);
    
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
      network.PLC_Stem = PLC_c(psirefS, VCstem_slope, VCstem_P50);
    }
  } else { //Immediate refilling
    network.Psi_SApo_cav = psirefS;
    network.PLC_Stem = PLC_c(psirefS, VCstem_slope, VCstem_P50);
  }
  if(leafCavitationRecovery!="total") {
    if(psirefL < Psi_LApo_cav) {
      network.Psi_LApo_cav = psirefL;
      network.PLC_Leaf = PLC_c(psirefL, VCleaf_slope, VCleaf_P50);
    }
  } else { //Immediate refilling
    network.Psi_LApo_cav = psirefL;
    network.PLC_Leaf = PLC_c(psirefL, VCleaf_slope, VCleaf_P50);
  }
}


//Function to clone SureauParams object
void copyParams_c(SureauParams& params, SureauParams& sinkParams) {
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
    sinkParams.VGrhizo_kmax[i] = params.VGrhizo_kmax[i];
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
void copyNetwork_c(SureauNetwork& network, SureauNetwork& sinkNetwork) {
  copyParams_c(network.params, sinkNetwork.params);
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
void deleteSureauNetworkPointers_c(SureauNetwork &network) {
  delete[] network.params.k_RCApoInit;
  delete[] network.params.VGrhizo_kmax;
  delete[] network.k_Soil;
  delete[] network.PsiSoil;
  delete[] network.k_RSApo;
  delete[] network.k_SoilToStem;
}



void initSureauParams_inner_c(SureauParams& params, int c,
                              InternalWater& internalWater, 
                              TranspirationParams& paramsTranspiration, 
                              WaterStorageParams& paramsWaterStorage,
                              std::vector<double>& VCroot_kmax, 
                              std::vector<double>& VGrhizo_kmax,
                              ControlParameters& control, 
                              double sapFluidityDay) {

  //This will change depending on the layers connected
  params.npools = VCroot_kmax.size();
  params.TPhase_gmin = control.advancedWB.TPhase_gmin;
  params.Q10_1_gmin = control.advancedWB.Q10_1_gmin;
  params.Q10_2_gmin = control.advancedWB.Q10_2_gmin;
  if(control.sureau.stomatalSubmodel=="Jarvis") {
    params.Tgs_optim = paramsTranspiration.Gs_Toptim[c];
    params.Tgs_sens = paramsTranspiration.Gs_Tsens[c];
    params.JarvisPAR = control.sureau.JarvisPAR;
  } else {
    params.Gsw_AC_slope = paramsTranspiration.Gsw_AC_slope[c];
  }
  params.fTRBToLeaf = control.sureau.fTRBToLeaf;//ratio of bark area to leaf area
  params.C_SApoInit = control.sureau.C_SApoInit; //Maximum capacitance of the stem apoplasm
  params.C_LApoInit = control.sureau.C_LApoInit; //Maximum capacitance of the leaf apoplasm
  params.k_SLApoInit = sapFluidityDay*paramsTranspiration.VCleafapo_kmax[c]; //Maximum conductance from trunk apoplasm to leaf apoplasm
  params.k_CSApoInit = sapFluidityDay*paramsTranspiration.VCstem_kmax[c]; //Maximum conductance from root crown to stem apoplasm
  // Rcout << "par init "<< c<< " VCstem_kmax[c]: "<<VCstem_kmax[c] << " sapfluidity: "<< sapFluidityDay<<" k_CSApoInit: " << params.k_CSApoInit<< "\n";
  params.k_RCApoInit = new double [params.npools];
  params.VGrhizo_kmax = new double [params.npools];
  for(int l = 0;l<params.npools;l++)  {
    params.k_RCApoInit[l] = sapFluidityDay*VCroot_kmax[l];
    params.VGrhizo_kmax[l] = VGrhizo_kmax[l];
  }
  params.slope_gs = paramsTranspiration.Gs_slope[c];
  params.P50_gs = paramsTranspiration.Gs_P50[c];
  
  //PLANT RELATED PARAMETERS
  double gmin20 = paramsTranspiration.Gswmin[c]*1000.0; //Leaf cuticular transpiration
  bool leafCuticularTranspiration = control.sureau.leafCuticularTranspiration;
  if(!leafCuticularTranspiration) gmin20 = 0.0;
  params.gmin20 = gmin20; //mol to mmol 
  params.gsMax = paramsTranspiration.Gswmax[c]*1000.0; //mol to mmol 
  bool stemCuticularTranspiration = control.sureau.stemCuticularTranspiration;
  double gmin_S = paramsTranspiration.Gswmin[c]*1000.0;  // gmin for stem equal to leaf gmin, in mmol
  if(!stemCuticularTranspiration) gmin_S = 0.0;
  params.gmin_S = gmin_S;
  double gs_NightFrac = control.sureau.gs_NightFrac;
  params.gsNight = gs_NightFrac*paramsTranspiration.Gswmax[c]*1000.0; 
  
  params.VCleaf_P50 = paramsTranspiration.VCleaf_P50[c]; 
  params.VCleaf_slope = paramsTranspiration.VCleaf_slope[c]; 
  params.VCstem_P50 = paramsTranspiration.VCstem_P50[c]; 
  params.VCstem_slope = paramsTranspiration.VCstem_slope[c]; 
  params.VCroot_P50 = paramsTranspiration.VCroot_P50[c]; 
  params.VCroot_slope = paramsTranspiration.VCroot_slope[c]; 
  params.PiFullTurgor_Leaf = paramsWaterStorage.LeafPI0[c]; 
  params.epsilonSym_Leaf = paramsWaterStorage.LeafEPS[c]; 
  params.PiFullTurgor_Stem = paramsWaterStorage.StemPI0[c]; 
  params.epsilonSym_Stem = paramsWaterStorage.StemEPS[c]; 
  
}

void initSureauNetwork_inner_c(SureauNetwork& network, int c, 
                               std::vector<double>& LAIphe,
                               InternalWater& internalWater, 
                               AnatomyParams& paramsAnatomy, 
                               TranspirationParams& paramsTranspiration, 
                               WaterStorageParams& paramsWaterStorage,
                               std::vector<double>& VCroot_kmax, std::vector<double>& VGrhizo_kmax,
                               std::vector<double>& PsiSoil, std::vector<double>& VG_n, std::vector<double>& VG_alpha,
                               ControlParameters& control, double sapFluidityDay) {
  
  //Params
  initSureauParams_inner_c(network.params, c, internalWater, 
                           paramsTranspiration, paramsWaterStorage,
                           VCroot_kmax, VGrhizo_kmax,
                           control, sapFluidityDay);
  
  //LAI
  network.LAI = LAIphe[c];
  
  //Water potentials
  network.Psi_LApo = internalWater.LeafPsi[c]; 
  network.Psi_LSym = internalWater.LeafSympPsi[c];
  network.Psi_RCApo = internalWater.RootCrownPsi[c];
  
  network.Psi_SApo = internalWater.StemPsi[c]; 
  network.Psi_SSym = internalWater.StemSympPsi[c];
  network.Psi_SApo_cav = std::min(0.0, invPLC_c(internalWater.StemPLC[c]*100.0, paramsTranspiration.VCstem_slope[c], paramsTranspiration.VCstem_P50[c])); //Sureau operates with %
  network.Psi_LApo_cav = std::min(0.0, invPLC_c(internalWater.LeafPLC[c]*100.0, paramsTranspiration.VCleaf_slope[c], paramsTranspiration.VCleaf_P50[c])); //Sureau operates with %
  
  //PLC levels
  network.PLC_Stem = internalWater.StemPLC[c]*100.0; //Sureau operates with %
  network.PLC_Leaf = internalWater.LeafPLC[c]*100.0; //Sureau operates with %
  //Capacitances (mmol m-2 MPa-1)
  network.C_SApo = medfate::NA_DOUBLE; //Capacitance of the stem apoplasm
  network.C_LApo = medfate::NA_DOUBLE; //Capacitance of the leaf apoplasm
  network.C_SSym = medfate::NA_DOUBLE; //Capacitance of the stem symplasm (HOW TO ESTIMATE THEM?)
  network.C_LSym = medfate::NA_DOUBLE; //Capacitance of the leaf symplasm (HOW TO ESTIMATE THEM?)
  //Conductances (mmol m-2 MPa-1 s-1)
  network.k_SLApo = medfate::NA_DOUBLE; //Conductance from trunk apoplasm to leaf apoplasm
  network.k_CSApo = medfate::NA_DOUBLE; //Conductance from root crown to trunk apoplasm
  network.k_SSym = sapFluidityDay*paramsTranspiration.kstem_symp[c]; //Conductance from trunk apoplasm to trunk symplasm (CONTROL PARAMETER?)
  network.k_LSym = sapFluidityDay*paramsTranspiration.kleaf_symp[c]; //Conductance from leaf apoplasm to leaf symplasm
  
  network.k_RSApo = new double[network.params.npools];
  network.k_SoilToStem = new double[network.params.npools];
  network.k_Soil = new double[network.params.npools];
  network.PsiSoil = new double[network.params.npools];
  for(int l=0;l < network.params.npools; l++) {
    network.PsiSoil[l] = PsiSoil[l]; // Copies external soil Psi
    network.k_SoilToStem[l] = medfate::NA_DOUBLE; //Conductance from soil to trunk apoplasm
    network.k_RSApo[l] = medfate::NA_DOUBLE; //Conductance from rhizosphere surface to trunk apoplasm
    network.k_Soil[l] = vanGenuchtenConductance_c(PsiSoil[l],
                                                  VGrhizo_kmax[l], 
                                                  VG_n[l], VG_alpha[l]); 
  }
  network.k_Plant = medfate::NA_DOUBLE; //Whole-plant conductance  = Plant_kmax
  
  //Water content (mmol m-2)
  double l2mmol = 1.0e6/18.0;
  network.Q_SApo_sat_mmol_perLeafArea = paramsWaterStorage.Vsapwood[c]*paramsWaterStorage.StemAF[c]*l2mmol; //Water content in stem apoplasm
  network.Q_LApo_sat_mmol_perLeafArea = paramsWaterStorage.Vleaf[c]*paramsWaterStorage.LeafAF[c]*l2mmol; //Water content in leaf apoplasm
  network.Q_SSym_sat_mmol_perLeafArea = paramsWaterStorage.Vsapwood[c]*(1.0 - paramsWaterStorage.StemAF[c])*l2mmol; //Water content in stem symplasm
  network.Q_LSym_sat_mmol_perLeafArea = paramsWaterStorage.Vleaf[c]*(1.0 - paramsWaterStorage.LeafAF[c])*l2mmol; //Water content in leaf symplasm
  
  //Flows (mmol m-2 s-1)
  network.Einst = internalWater.Einst[c]; //Total transpiration
  network.Einst_SL = medfate::NA_DOUBLE; //Total transpiration (sunlit leaves)
  network.Einst_SH = medfate::NA_DOUBLE; //Total transpiration (shade leaves)
  network.Elim = internalWater.Elim[c]; //Stomatal transpiration
  network.Elim_SL = medfate::NA_DOUBLE; //Stomatal transpiration (sunlit leaves)
  network.Elim_SH = medfate::NA_DOUBLE; //Stomatal transpiration (shade leaves)
  network.Emin_L = internalWater.Emin_L[c]; //Leaf cuticular transpiration
  network.Emin_L_SL = medfate::NA_DOUBLE; //Leaf cuticular transpiration (sunlit leaves)
  network.Emin_L_SH = medfate::NA_DOUBLE; //Leaf cuticular transpiration (shade leaves)
  network.Emin_S = internalWater.Emin_S[c]; //Stem cuticular transpiration
  
  //Diagnostics
  network.Diag_nwhile_cavit = medfate::NA_INTEGER;
  network.Diag_deltaRegulMax = medfate::NA_DOUBLE;
  network.Diag_deltaPLCMax = medfate::NA_DOUBLE;
  network.Diag_timeStepInSeconds = medfate::NA_DOUBLE;
  
  // Update plant conductances and capacitances according to network status
  update_conductances_c(network);
  update_capacitances_c(network);
}

