#define STRICT_R_HEADERS
#include <RcppArmadillo.h>
#include "inner_sureau_c.h"
#include "photosynthesis.h"
#include "biophysicsutils_c.h"
#include "hydraulics.h"
#include "hydraulics_c.h"
#include "soil.h"
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

void semi_implicit_integration_inner_c(SureauNetwork& network,
                                       double dt, 
                                       const SureauOpt& opt, 
                                       const std::string& stemCavitationRecovery = "annual", 
                                       const std::string& leafCavitationRecovery = "total") {
  
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
