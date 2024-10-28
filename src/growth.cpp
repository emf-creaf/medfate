// [[Rcpp::interfaces(r,cpp)]]
#define STRICT_R_HEADERS
#include <Rcpp.h>
#include <numeric>
#include "communication_structures.h"
#include "lightextinction_basic.h"
#include "lightextinction_advanced.h"
#include "phenology.h"
#include "biophysicsutils.h"
#include "forestutils.h"
#include "fireseverity.h"
#include "firebehaviour.h"
#include "tissuemoisture.h"
#include "hydraulics.h"
#include "hydrology.h"
#include "carbon.h"
#include "root.h"
#include "woodformation.h"
#include "soil.h"
#include "spwb.h"
#include "spwb_day.h"
#include <meteoland.h>
using namespace Rcpp;

//' Mortality
//' 
//' A simple sigmoid function to determine a daily mortality likelihood according 
//' to the value of a stress variable.
//'
//' @param stressValue Current value of the stress variable (0 to 1, 
//'                    with higher values indicate stronger stress).
//' @param stressThreshold Threshold to indicate 50\% annual mortality probability.
//' 
//' @return Returns a probability (between 0 and 1)
//' 
//' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
//' 
//' @seealso \code{\link{growth}}, \code{\link{regeneration}}
//' 
//' @keywords internal
// [[Rcpp::export("mortality_dailyProbability")]]
double dailyMortalityProbability(double stressValue, double stressThreshold) {
  double exponent = 40.0;
  double y = (stressValue - stressThreshold);
  double P_annual = std::min(1.0, 1.0 - exp(exponent*y)/(1.0 + exp(exponent*y)));
  double P_daily = 1.0 - exp(log(1.0 - P_annual)/356.0);
  return(P_daily);
}

/**
 * phloem flow (Holtta et al. 2017)
 *  psiUpstream, psiDownstream - water potential upstream (leaves)  and downstream
 *  concUpstream, concDownstream - sugar concentration upstream (leaves) and downstream (stem)
 *  k_f - phloem conductance per leaf area basis (l*m-2*MPa-1*s-1)
 *  
 *  out mol*s-1*m-2 (flow per leaf area basis)
 */
double phloemFlow(double psiUpstream, double psiDownstream,
                 double concUpstream, double concDownstream,
                 double temp, double k_f, double nonSugarConc) {
  double turgor_up = turgor(psiUpstream, concUpstream, temp, nonSugarConc);
  double turgor_down = turgor(psiDownstream, concDownstream, temp, nonSugarConc);
  if(temp < 0.0) k_f = 0.0; // No phloem flow if temperature below zero
  double relVisc = relativeSapViscosity((concUpstream+concDownstream)/2.0, temp);
  if(turgor_up>turgor_down) {
    return(k_f*concUpstream*(turgor_up - turgor_down)/relVisc);
  } else {
    return(k_f*concDownstream*(turgor_up - turgor_down)/relVisc);
  }
}
// // [[Rcpp::export("growth_dailyphloemFlow")]]
// NumericMatrix dailyPhloemFlow(List x, List spwbOut, 
//                              NumericVector concLeaf, NumericVector concSapwood) {
//   DataFrame paramStorage =  Rcpp::as<Rcpp::DataFrame>(x["paramsWaterStorage"]);
//   NumericVector Vleaf = Rcpp::as<Rcpp::NumericVector>(paramStorage["Vleaf"]);
//   NumericVector stemPI0 = Rcpp::as<Rcpp::NumericVector>(paramStorage["StemPI0"]);
//   NumericVector leafPI0 = Rcpp::as<Rcpp::NumericVector>(paramStorage["LeafPI0"]);
//   NumericVector stemEPS = Rcpp::as<Rcpp::NumericVector>(paramStorage["StemEPS"]);
//   NumericVector leafEPS = Rcpp::as<Rcpp::NumericVector>(paramStorage["LeafEPS"]);
//   List plantsInst = spwbOut["PlantsInst"];  
//   NumericMatrix rwcStem =  Rcpp::as<Rcpp::NumericMatrix>(plantsInst["RWCstem"]);
//   NumericMatrix rwcLeaf =  Rcpp::as<Rcpp::NumericMatrix>(plantsInst["RWCleaf"]);
//   List eb = spwbOut["EnergyBalance"];  
//   DataFrame tempDF =  Rcpp::as<Rcpp::DataFrame>(eb["Temperature"]);
//   NumericVector Tcan = Rcpp::as<Rcpp::NumericVector>(tempDF["Tcan"]);
//   
//   int numCohorts = stemPI0.length();
//   int numSteps = Tcan.length();
//   NumericMatrix ff(numCohorts, numSteps);
//   for(int c=0;c<numCohorts;c++) {
//     for(int s=0;s<numSteps;s++) {
//       // double leafPI = osmoticWaterPotential(concLeaf[c], Tcan[s]);
//       // double sapwoodPI = osmoticWaterPotential(concs[c], Tcan[s]);
//       double psiUp = symplasticWaterPotential(rwcLeaf(c,s), leafPI0[c], leafEPS[c]);
//       double psiDown = symplasticWaterPotential(rwcStem(c,s), stemPI0[c], stemEPS[c]);
//       ff(c,s) = phloemFlow(psiUp, psiDown, concLeaf[c], concSapwood[c], Tcan[s], k_phloem, nonSugarConc)*3600.0; //flow as mol per hour and leaf area basis
//     }
//   }
//   ff.attr("dimnames") = rwcStem.attr("dimnames");
//   return(ff);
// }

double qResp(double Tmean) {
  // Tjoelker, M. G., J. Oleksyn, and P. B. Reich. 2001. Modelling respiration of vegetation: Evidence for a general temperature-dependent Q10. Global Change Biology 7:223–230.
  double Q10_resp = 3.22 - 0.046 * Tmean; 
  return(pow(Q10_resp,(Tmean-20.0)/10.0));
}

// double storageTransferRelativeRate(double fastCstorage, double fastCstoragemax) {
//   double f = ((2.0/(1.0+exp(-5.0*((fastCstorage/fastCstoragemax)-tlpConcSapwood)/tlpConcSapwood)))-1.0);
//   return(f);
// }
// double carbonGrowthFactor(double conc, double threshold) {
//   double k =10.0;
//   return(std::max(0.0,(1.0 - exp(k*(threshold-conc)))/(1.0 - exp(k*(-conc)))));
// }


void fillInitialPlantBiomassBalance(DataFrame plantBiomassBalance, DataFrame ccIni, DataFrame above) {

  NumericVector N = above["N"];
  int numCohorts = N.size();

  //Initial Biomass compartments
  NumericVector SapwoodBiomass = Rcpp::as<Rcpp::NumericVector>(ccIni["SapwoodStructuralBiomass"]);
  NumericVector TotalBiomass = Rcpp::as<Rcpp::NumericVector>(ccIni["TotalBiomass"]);
  NumericVector TotalLivingBiomass= Rcpp::as<Rcpp::NumericVector>(ccIni["TotalLivingBiomass"]);
  NumericVector LabileBiomass = Rcpp::as<Rcpp::NumericVector>(ccIni["LabileBiomass"]);
  NumericVector StructuralBiomass = Rcpp::as<Rcpp::NumericVector>(ccIni["StructuralBiomass"]);
  
  NumericVector Nprev = Rcpp::as<Rcpp::NumericVector>(plantBiomassBalance["InitialDensity"]);
  NumericVector InitialSapwoodBiomass = Rcpp::as<Rcpp::NumericVector>(plantBiomassBalance["InitialSapwoodBiomass"]);
  NumericVector InitialStructuralBiomass = Rcpp::as<Rcpp::NumericVector>(plantBiomassBalance["InitialStructuralBiomass"]);
  NumericVector InitialLabileBiomass = Rcpp::as<Rcpp::NumericVector>(plantBiomassBalance["InitialLabileBiomass"]);
  NumericVector InitialPlantBiomass = Rcpp::as<Rcpp::NumericVector>(plantBiomassBalance["InitialPlantBiomass"]);
  NumericVector InitialLivingPlantBiomass = Rcpp::as<Rcpp::NumericVector>(plantBiomassBalance["InitialLivingPlantBiomass"]);
  NumericVector InitialCohortBiomass = Rcpp::as<Rcpp::NumericVector>(plantBiomassBalance["InitialCohortBiomass"]);

  for(int c=0; c < numCohorts;c++) {
    Nprev[c] = N[c];
    InitialSapwoodBiomass[c] = SapwoodBiomass[c];
    InitialStructuralBiomass[c] = StructuralBiomass[c];
    InitialLabileBiomass[c] = LabileBiomass[c];
    InitialPlantBiomass[c] = TotalBiomass[c];
    InitialCohortBiomass[c] = TotalBiomass[c]*(Nprev[c]/10000.0);
    InitialLivingPlantBiomass[c] = TotalLivingBiomass[c];
    InitialStructuralBiomass[c] = StructuralBiomass[c];
  }
}


void closePlantBiomassBalance(List initialFinalCC, DataFrame plantBiomassBalance, List x,
                         NumericVector LabileCarbonBalance,
                         NumericVector LeafBiomassBalance,
                         NumericVector FineRootBiomassBalance) {
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  int numCohorts = above.nrow();
  
  NumericVector Nfinal = above["N"];

  DataFrame ccFin = as<DataFrame>(initialFinalCC["ccFin_g_ind"]);
  fillCarbonCompartments(ccFin, x, "g_ind");

  NumericVector finalSapwoodBiomass_ind= Rcpp::as<Rcpp::NumericVector>(ccFin["SapwoodStructuralBiomass"]);
  NumericVector plantFinalBiomass_ind = Rcpp::as<Rcpp::NumericVector>(ccFin["TotalBiomass"]);
  NumericVector cohortFinalBiomass_m2 = Rcpp::as<Rcpp::NumericVector>(ccFin["TotalBiomass"]);
  NumericVector labileFinalBiomass_ind = Rcpp::as<Rcpp::NumericVector>(ccFin["LabileBiomass"]);
  NumericVector structuralFinalBiomass_ind = Rcpp::as<Rcpp::NumericVector>(ccFin["StructuralBiomass"]);

  NumericVector Nprev = Rcpp::as<Rcpp::NumericVector>(plantBiomassBalance["InitialDensity"]);
  
  NumericVector InitialSapwoodBiomass = Rcpp::as<Rcpp::NumericVector>(plantBiomassBalance["InitialSapwoodBiomass"]);
  NumericVector StructuralBiomassBalance = Rcpp::as<Rcpp::NumericVector>(plantBiomassBalance["StructuralBiomassBalance"]);
  NumericVector InitialStructuralBiomass = Rcpp::as<Rcpp::NumericVector>(plantBiomassBalance["InitialStructuralBiomass"]);
  NumericVector StructuralBiomassChange = Rcpp::as<Rcpp::NumericVector>(plantBiomassBalance["StructuralBiomassChange"]);
  
  NumericVector LabileBiomassBalance = Rcpp::as<Rcpp::NumericVector>(plantBiomassBalance["LabileBiomassBalance"]);
  NumericVector InitialLabileBiomass = Rcpp::as<Rcpp::NumericVector>(plantBiomassBalance["InitialLabileBiomass"]);
  NumericVector LabileBiomassChange = Rcpp::as<Rcpp::NumericVector>(plantBiomassBalance["LabileBiomassChange"]);

  NumericVector PlantBiomassBalance = Rcpp::as<Rcpp::NumericVector>(plantBiomassBalance["PlantBiomassBalance"]);
  NumericVector InitialPlantBiomass = Rcpp::as<Rcpp::NumericVector>(plantBiomassBalance["InitialPlantBiomass"]);
  NumericVector InitialLivingPlantBiomass = Rcpp::as<Rcpp::NumericVector>(plantBiomassBalance["InitialLivingPlantBiomass"]);
  NumericVector PlantBiomassChange = Rcpp::as<Rcpp::NumericVector>(plantBiomassBalance["PlantBiomassChange"]);

  NumericVector MortalityBiomassLoss = Rcpp::as<Rcpp::NumericVector>(plantBiomassBalance["MortalityBiomassLoss"]);
  
  NumericVector CohortBiomassBalance = Rcpp::as<Rcpp::NumericVector>(plantBiomassBalance["CohortBiomassBalance"]);
  NumericVector InitialCohortBiomass = Rcpp::as<Rcpp::NumericVector>(plantBiomassBalance["InitialCohortBiomass"]);
  NumericVector CohortBiomassChange = Rcpp::as<Rcpp::NumericVector>(plantBiomassBalance["CohortBiomassChange"]);
  
  //PLANT BIOMASS balance (g_ind)
  for(int j=0; j<numCohorts;j++) {
    double sapwoodBiomassBalance = finalSapwoodBiomass_ind[j] - InitialSapwoodBiomass[j];
    StructuralBiomassChange[j] = structuralFinalBiomass_ind[j] - InitialStructuralBiomass[j];
    LabileBiomassChange[j] = labileFinalBiomass_ind[j] - InitialLabileBiomass[j];
    PlantBiomassChange[j] = plantFinalBiomass_ind[j] - InitialPlantBiomass[j];
    
    StructuralBiomassBalance[j] = LeafBiomassBalance[j] + sapwoodBiomassBalance + FineRootBiomassBalance[j];
    LabileBiomassBalance[j] = LabileCarbonBalance[j]*InitialLivingPlantBiomass[j];
    
    PlantBiomassBalance[j] = LabileBiomassBalance[j] + StructuralBiomassBalance[j];
    
    //Biomass loss as the decrease in density multiplied by the total biomass after including individual biomass changes (g/m2)
    MortalityBiomassLoss[j] = (Nprev[j] - Nfinal[j])*(InitialPlantBiomass[j]+PlantBiomassBalance[j])/(10000.0);
    
    //Change units to g/m2    
    StructuralBiomassBalance[j] = StructuralBiomassBalance[j]*(Nprev[j]/10000.0);
    LabileBiomassBalance[j] = LabileBiomassBalance[j]*(Nprev[j]/10000.0);
    PlantBiomassBalance[j] = PlantBiomassBalance[j]*(Nprev[j]/10000.0);
    
    //COHORT BIOMASS balance (g/m2) 
    CohortBiomassBalance[j] = PlantBiomassBalance[j] - MortalityBiomassLoss[j];
    CohortBiomassChange[j] = cohortFinalBiomass_m2[j]*(Nfinal[j]/10000.0) - InitialCohortBiomass[j];
  }
  
}

NumericVector standLevelBiomassBalance(DataFrame biomassBalance) {
  return(NumericVector::create(
      _["StructuralBalance"] = sum(Rcpp::as<Rcpp::NumericVector>(biomassBalance["StructuralBiomassBalance"])),
      _["LabileBalance"] = sum(Rcpp::as<Rcpp::NumericVector>(biomassBalance["LabileBiomassBalance"])),
      _["PlantBalance"] = sum(Rcpp::as<Rcpp::NumericVector>(biomassBalance["PlantBiomassBalance"])),
      _["MortalityLoss"] = sum(Rcpp::as<Rcpp::NumericVector>(biomassBalance["MortalityBiomassLoss"])),
      _["CohortBalance"] = sum(Rcpp::as<Rcpp::NumericVector>(biomassBalance["CohortBiomassBalance"]))
  ));
}

void updateStructuralVariables(List x, NumericVector deltaSAgrowth) {
  
  //Control params
  List control = x["control"];  
  bool shrubDynamics = control["shrubDynamics"];
  bool herbDynamics = control["herbDynamics"];
  
  //Cohort info
  DataFrame cohorts = Rcpp::as<Rcpp::DataFrame>(x["cohorts"]);
  IntegerVector SP = Rcpp::as<Rcpp::IntegerVector>(cohorts["SP"]);
  int numCohorts = SP.size();
  
  //Aboveground parameters  
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  NumericVector DBH = above["DBH"];
  NumericVector Cover = above["Cover"];
  NumericVector H = above["H"];
  NumericVector N = above["N"];
  NumericVector CR = above["CR"];
  NumericVector LAI_live = above["LAI_live"];
  NumericVector LAI_expanded = above["LAI_expanded"];
  NumericVector LAI_dead = above["LAI_dead"];
  NumericVector LAI_nocomp = above["LAI_nocomp"];
  NumericVector SA = above["SA"];
  
  //Allometric parameters
  DataFrame paramsAllometries = Rcpp::as<Rcpp::DataFrame>(x["paramsAllometries"]);
  NumericVector Afbt  = paramsAllometries["Afbt"];
  NumericVector Bfbt  = paramsAllometries["Bfbt"];
  NumericVector Cfbt  = paramsAllometries["Cfbt"];
  NumericVector Aash  = paramsAllometries["Aash"];
  NumericVector Bash  = paramsAllometries["Bash"];
  NumericVector Absh  = paramsAllometries["Absh"];
  NumericVector Bbsh  = paramsAllometries["Bbsh"];
  NumericVector Acw  = paramsAllometries["Acw"];
  NumericVector Bcw  = paramsAllometries["Bcw"];
  NumericVector Acr  = paramsAllometries["Acr"];
  NumericVector B1cr  = paramsAllometries["B1cr"];
  NumericVector B2cr  = paramsAllometries["B2cr"];
  NumericVector B3cr  = paramsAllometries["B3cr"];
  NumericVector C1cr  = paramsAllometries["C1cr"];
  NumericVector C2cr  = paramsAllometries["C2cr"];
  
  //Phenology parameters
  DataFrame internalPhenology = Rcpp::as<Rcpp::DataFrame>(x["internalPhenology"]);
  LogicalVector budFormation = internalPhenology["budFormation"];
  
  //Allometric parameters
  DataFrame paramsGrowth = Rcpp::as<Rcpp::DataFrame>(x["paramsGrowth"]);
  NumericVector fHDmin= paramsGrowth["fHDmin"];
  NumericVector fHDmax= paramsGrowth["fHDmax"];
  
  //Inteception parameters
  DataFrame paramsInterception = Rcpp::as<Rcpp::DataFrame>(x["paramsInterception"]);
  NumericVector kPAR = paramsInterception["kPAR"];
  
  //Anatomy parameters
  DataFrame paramsAnatomy = Rcpp::as<Rcpp::DataFrame>(x["paramsAnatomy"]);
  NumericVector Hmax  = paramsAnatomy["Hmax"];
  NumericVector SLA = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["SLA"]);
  NumericVector Al2As = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["Al2As"]);
  NumericVector r635  = paramsAnatomy["r635"];
  
  DataFrame internalAllocation = Rcpp::as<Rcpp::DataFrame>(x["internalAllocation"]);
  NumericVector leafAreaTarget = internalAllocation["leafAreaTarget"];
  NumericVector sapwoodAreaTarget = internalAllocation["sapwoodAreaTarget"];
  
  //Update DBH
  NumericVector deltaDBH(numCohorts, 0.0);
  for(int j=0;j<numCohorts; j++) {
    if(!NumericVector::is_na(DBH[j])) {
      deltaDBH[j] = 2.0*sqrt(pow(DBH[j]/2.0,2.0)+(deltaSAgrowth[j]/M_PI)) - DBH[j];
      DBH[j] = DBH[j] + deltaDBH[j];
    }
  }
  //Update height
  NumericVector L = parcohortC(H, LAI_live, LAI_dead, kPAR, CR);
  for(int j=0;j<numCohorts; j++) {
    if(!NumericVector::is_na(DBH[j]) && N[j]>0.0) {
      double fHmod = std::max(0.0,std::min(1.0,(1.0-((H[j]-137.0)/(Hmax[j]-137.0)))));
      double fHD = (fHDmin[j]*(L[j]/100.0) + fHDmax[j]*(1.0-(L[j]/100.0)))*fHmod;
      H[j] = H[j] + fHD*deltaDBH[j];
    }
  }
  //Update crown ratio
  NumericVector crNew = treeCrownRatioAllometric(N, DBH, H, Acw, Bcw, Acr, B1cr, B2cr, B3cr, C1cr, C2cr);
  for(int j=0;j<numCohorts; j++) {
    if(!NumericVector::is_na(DBH[j]) && N[j]>0.0) {
      CR[j] = crNew[j];
    }
  }
  //Update tree leaf area target
  NumericVector ltba = largerTreeBasalArea(N, DBH, 1.0); //Allometries were calibrated including the target cohort
  for(int j=0;j<numCohorts;j++) {
    if(!NumericVector::is_na(DBH[j]) && N[j]>0.0) {
      double leafBiomassNoComp = Afbt[j]*pow(std::min(100.0,DBH[j]), Bfbt[j])*exp(-0.0001*N[j]);//Correct for high density packing
      LAI_nocomp[j] = SLA[j]*leafBiomassNoComp*N[j]/10000.0; //LAI without competition effect
      if(budFormation[j]) { //Update target if buds are active
        // Rcout <<j<< " "<< ltba[j]<< " "<<leafAreaTarget[j];
        leafAreaTarget[j] = SLA[j]*leafBiomassNoComp*exp(Cfbt[j]*ltba[j]); //Include competition effect in leaf biomass estimation
        LAI_live[j] = leafAreaTarget[j]*N[j]/10000.0;
        // Rcout << " "<< leafAreaTarget[j]<<"\n";
      }
    }
  }
  //Shrub variables
  if(shrubDynamics) {
    double treeLAI = 0.0;
    for(int j=0;j<numCohorts;j++) {
      if(!NumericVector::is_na(DBH[j])) treeLAI +=LAI_live[j];
    }
    for(int j=0;j<numCohorts; j++) {
      if(NumericVector::is_na(DBH[j]) && N[j]>0.0) {
        if(budFormation[j]) {
          leafAreaTarget[j] = (Al2As[j]*SA[j])/10000.0; // Set leaf area target according to current sapwood area
          double Wleaves = leafAreaTarget[j]/SLA[j];  //Calculates the biomass (kg dry weight) of leaves
          Wleaves = Wleaves/exp(-0.235*treeLAI); //Correct depending on tree leaf area
          double PV = pow(Wleaves*r635[j]/Absh[j], 1.0/Bbsh[j]); //Calculates phytovolume (in m3/ind)
          H[j] = pow(1e6*PV/Aash[j], 1.0/(1.0+Bash[j])); //Updates shrub height
          // Rcout<< Wleaves << " " << PV << " " << H[j]<<"\n";
          if(H[j]> Hmax[j]) { //Limit height (and update the former variables)
            H[j] = Hmax[j];
            PV = (Aash[j]/1e6)*pow(H[j], (1.0+Bash[j])); //recalculate phytovolume from H
            Wleaves = (Absh[j]/r635[j])*pow(PV, Bbsh[j]); //recalculate Wleaves from phytovolume
            Wleaves = Wleaves*exp(-0.235*treeLAI); //Correct depending on tree leaf area
            leafAreaTarget[j] = Wleaves * SLA[j]; //recalculate leaf area target from Wleaves
            sapwoodAreaTarget[j] = 10000.0*leafAreaTarget[j]/Al2As[j]; //Set target sapwood area (may generate sapwood senescence)
          }
          Cover[j] = std::min(100.0, N[j]*Aash[j]*pow(H[j],Bash[j])/1e6); //Updates shrub cover
          LAI_live[j] = leafAreaTarget[j]*N[j]/10000.0;
        }
        //Update LAI without tree competition effect
        LAI_nocomp[j] = LAI_live[j]/exp(-0.235*treeLAI);
      }
    }
  }
  //Herb variables
  if(herbDynamics) {
    double woodyLAI = sum(LAI_live);
    double herbLAImax = x["herbLAImax"];
    x["herbLAI"] = herbLAImax*exp(-0.235*woodyLAI);
  }
}
void growthDay_private(List internalCommunication, List x, NumericVector meteovec, 
                    double latitude, double elevation, double slope, double aspect,
                    double solarConstant, double delta, 
                    double runon = 0.0, Nullable<NumericVector> lateralFlows = R_NilValue, double waterTableDepth = NA_REAL, 
                    bool verbose = false) {
  
  
  //Get previous PLC so that defoliation occurs only when PLC increases
  //Internal state variables
  DataFrame internalWater = Rcpp::as<Rcpp::DataFrame>(x["internalWater"]);
  NumericVector StemPLCprev = clone(Rcpp::as<Rcpp::NumericVector>(internalWater["StemPLC"]));
  
  //Control params
  List control = x["control"];  
  String transpirationMode = control["transpirationMode"];
  List modelOutput;
  List spwbOutput;
  
  //Soil-plant water balance (this creates communication structures as well)
  if(transpirationMode=="Granier") {
    spwbOutput = internalCommunication["basicSPWBOutput"];
    modelOutput = internalCommunication["basicGROWTHOutput"];
    spwbDay_basic(internalCommunication, x, meteovec, 
                  elevation, slope, aspect,
                  runon, lateralFlows, waterTableDepth,
                  verbose); 
  } else {
    spwbOutput = internalCommunication["advancedSPWBOutput"];
    modelOutput = internalCommunication["advancedGROWTHOutput"];
    spwbDay_advanced(internalCommunication, x, meteovec, 
                     latitude, elevation, slope, aspect,
                     solarConstant, delta, 
                     runon, lateralFlows, waterTableDepth,
                     verbose);
  }

  
  String soilFunctions = control["soilFunctions"];
  String mortalityMode = control["mortalityMode"];
  int ntimesteps = control["ndailysteps"];
  bool subdailyCarbonBalance = control["subdailyCarbonBalance"];
  double mortalityRelativeSugarThreshold= control["mortalityRelativeSugarThreshold"];
  double mortalityRWCThreshold= control["mortalityRWCThreshold"];
  bool allowDessication = control["allowDessication"];
  bool allowStarvation = control["allowStarvation"];
  bool sinkLimitation = control["sinkLimitation"];
  bool shrubDynamics = control["shrubDynamics"];
  String allocationStrategy = control["allocationStrategy"];
  if(transpirationMode=="Granier") allocationStrategy = "Al2As";
  String rhizosphereOverlap = control["rhizosphereOverlap"];
  bool plantWaterPools = (rhizosphereOverlap!="total");
  bool taper = control["taper"];
  double averageFracRhizosphereResistance = control["averageFracRhizosphereResistance"];
  double phloemConductanceFactor = control["phloemConductanceFactor"];
  double nonSugarConcentration = control["nonSugarConcentration"];
  List equilibriumOsmoticConcentration  = control["equilibriumOsmoticConcentration"];
  double equilibriumLeafTotalConc = equilibriumOsmoticConcentration["leaf"];
  double equilibriumSapwoodTotalConc = equilibriumOsmoticConcentration["sapwood"];
  
  //Weather
  double tday = meteovec["tday"];
  double tmax = meteovec["tmax"];
  double Patm = meteovec["Patm"];
  double pfire = meteovec["pfire"];

  bool fireOccurrence = false;
  NumericVector fireBehavior(1,NA_REAL);
  if(R::runif(0.0,1.0) < pfire) {
    fireBehavior = fccsHazard(x, meteovec, spwbOutput, slope);
    fireOccurrence = true;
  }
  
  //Cohort info
  DataFrame cohorts = Rcpp::as<Rcpp::DataFrame>(x["cohorts"]);
  IntegerVector SP = Rcpp::as<Rcpp::IntegerVector>(cohorts["SP"]);
  int numCohorts = SP.size();
  
  //Soil
  DataFrame soil = Rcpp::as<Rcpp::DataFrame>(x["soil"]);
  int nlayers = soil.nrow();
  NumericVector psiSoil = psi(soil,"soilFunctions");
  NumericVector widths = soil["widths"];
  NumericVector rfc = soil["rfc"];
  NumericVector Ksat = soil["Ksat"];
  NumericVector VG_n = soil["VG_n"];
  NumericVector VG_alpha = soil["VG_alpha"];
  NumericVector Tsoil = soil["Temp"];
  
  //Aboveground parameters  
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  NumericVector DBH = above["DBH"];
  NumericVector Cover = above["Cover"];
  NumericVector H = above["H"];
  NumericVector N = above["N"];
  NumericVector CR = above["CR"];
  NumericVector LAI_live = above["LAI_live"];
  NumericVector LAI_expanded = above["LAI_expanded"];
  NumericVector LAI_dead = above["LAI_dead"];
  NumericVector LAI_nocomp = above["LAI_nocomp"];
  NumericVector SA = above["SA"];

  //Belowground parameters  
  DataFrame belowdf = Rcpp::as<Rcpp::DataFrame>(x["below"]);
  NumericVector Z95 = Rcpp::as<Rcpp::NumericVector>(belowdf["Z95"]);
  NumericVector Z50 = Rcpp::as<Rcpp::NumericVector>(belowdf["Z50"]);
  NumericVector fineRootBiomass = Rcpp::as<Rcpp::NumericVector>(belowdf["fineRootBiomass"]);
  NumericVector CRSV = Rcpp::as<Rcpp::NumericVector>(belowdf["coarseRootSoilVolume"]);
  List belowLayers = Rcpp::as<Rcpp::List>(x["belowLayers"]);
  List RHOP;
  if(plantWaterPools) RHOP = belowLayers["RHOP"];
  NumericMatrix V = Rcpp::as<Rcpp::NumericMatrix>(belowLayers["V"]);
  NumericMatrix L = Rcpp::as<Rcpp::NumericMatrix>(belowLayers["L"]);
  NumericMatrix RhizoPsi, VCroot_kmax, VGrhizo_kmax;
  if((transpirationMode=="Sperry") || (transpirationMode=="Sureau")) {
    RhizoPsi = Rcpp::as<Rcpp::NumericMatrix>(belowLayers["RhizoPsi"]);
    VCroot_kmax = Rcpp::as<Rcpp::NumericMatrix>(belowLayers["VCroot_kmax"]);
    VGrhizo_kmax = Rcpp::as<Rcpp::NumericMatrix>(belowLayers["VGrhizo_kmax"]);
  }
  
  //Internal state variables
  internalWater = Rcpp::as<Rcpp::DataFrame>(x["internalWater"]);

  //Values at the end of the day (after calling spwb)
  NumericVector StemPLC = Rcpp::as<Rcpp::NumericVector>(internalWater["StemPLC"]);
  NumericVector PlantPsi, psiApoLeaf, psiApoStem, psiSympLeaf, psiSympStem, psiRootCrown;
  NumericVector LeafPLC = Rcpp::as<Rcpp::NumericVector>(internalWater["LeafPLC"]);
  if(transpirationMode=="Granier") {
    PlantPsi  = Rcpp::as<Rcpp::NumericVector>(internalWater["PlantPsi"]);
    psiApoLeaf = Rcpp::as<Rcpp::NumericVector>(internalWater["PlantPsi"]);
  } else {
    psiApoLeaf = Rcpp::as<Rcpp::NumericVector>(internalWater["LeafPsi"]);
    if(transpirationMode == "Sperry") {
      psiApoStem = Rcpp::as<Rcpp::NumericVector>(internalWater["StemPsi"]);
      psiSympLeaf = Rcpp::as<Rcpp::NumericVector>(internalWater["LeafSympPsi"]);
    } else {
      psiApoStem = Rcpp::as<Rcpp::NumericVector>(internalWater["StemPsi"]);
      psiSympLeaf = Rcpp::as<Rcpp::NumericVector>(internalWater["LeafPsi"]);
    }
    psiRootCrown = Rcpp::as<Rcpp::NumericVector>(internalWater["RootCrownPsi"]);
    psiSympStem = Rcpp::as<Rcpp::NumericVector>(internalWater["StemSympPsi"]);
  }
  
  DataFrame internalCarbon = Rcpp::as<Rcpp::DataFrame>(x["internalCarbon"]);
  NumericVector sugarLeaf = internalCarbon["sugarLeaf"]; //Concentrations assuming RWC = 1
  NumericVector starchLeaf = internalCarbon["starchLeaf"];
  NumericVector sugarSapwood = internalCarbon["sugarSapwood"];
  NumericVector starchSapwood = internalCarbon["starchSapwood"];

  DataFrame internalMortality = Rcpp::as<Rcpp::DataFrame>(x["internalMortality"]);
  NumericVector N_dead = internalMortality["N_dead"];
  NumericVector N_starvation = internalMortality["N_starvation"];
  NumericVector N_dessication = internalMortality["N_dessication"];
  NumericVector N_burnt = internalMortality["N_burnt"];
  NumericVector Cover_dead = internalMortality["Cover_dead"];
  NumericVector Cover_starvation = internalMortality["Cover_starvation"];
  NumericVector Cover_dessication = internalMortality["Cover_dessication"];
  NumericVector Cover_burnt = internalMortality["Cover_burnt"];
  
  DataFrame internalAllocation = Rcpp::as<Rcpp::DataFrame>(x["internalAllocation"]);
  NumericVector allocationTarget = internalAllocation["allocationTarget"];
  NumericVector sapwoodAreaTarget = internalAllocation["sapwoodAreaTarget"];
  NumericVector leafAreaTarget = internalAllocation["leafAreaTarget"];
  NumericVector fineRootBiomassTarget = internalAllocation["fineRootBiomassTarget"];
  NumericVector crownBudPercent = internalAllocation["crownBudPercent"];
  
  DataFrame internalPhenology = Rcpp::as<Rcpp::DataFrame>(x["internalPhenology"]);
  LogicalVector leafUnfolding = internalPhenology["leafUnfolding"];
  LogicalVector budFormation = internalPhenology["budFormation"];
  LogicalVector leafSenescence = internalPhenology["leafSenescence"];
  
  DataFrame Plants = Rcpp::as<Rcpp::DataFrame>(spwbOutput["Plants"]);
  List PlantsInst;
  NumericVector Ag = Plants["GrossPhotosynthesis"];
  NumericVector LFMC = Plants["LFMC"];
  NumericVector PARcohort= Plants["FPAR"];
  NumericMatrix AgStep, AnStep;
  int numSteps = 1;
  if(transpirationMode!="Granier") {
    PlantsInst = spwbOutput["PlantsInst"];
    AgStep  =  Rcpp::as<Rcpp::NumericMatrix>(PlantsInst["Ag"]);
    AnStep  =  Rcpp::as<Rcpp::NumericMatrix>(PlantsInst["An"]);
    numSteps = AgStep.ncol();
  }

  //Data from spwb
  NumericVector Tcan(ntimesteps);
  NumericMatrix StemSympPsiInst, LeafSympPsiInst;
  List eb;
  double tcan_day = NA_REAL;
  if((transpirationMode=="Sperry") || (transpirationMode=="Sureau")) {
    StemSympPsiInst =  Rcpp::as<Rcpp::NumericMatrix>(PlantsInst["StemSympPsi"]);
    LeafSympPsiInst =  Rcpp::as<Rcpp::NumericMatrix>(PlantsInst["LeafSympPsi"]);

    eb = spwbOutput["EnergyBalance"];  
    DataFrame tempDF =  Rcpp::as<Rcpp::DataFrame>(eb["Temperature"]);
    NumericVector TcanIN = Rcpp::as<Rcpp::NumericVector>(tempDF["Tcan"]);
    for(int n=0;n<ntimesteps;n++) Tcan[n] = TcanIN[n];
    tcan_day = meteoland::utils_averageDaylightTemperature(min(Tcan), max(Tcan));
  }
  

  //Anatomy parameters
  DataFrame paramsAnatomy = Rcpp::as<Rcpp::DataFrame>(x["paramsAnatomy"]);
  NumericVector Hmed = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["Hmed"]);
  NumericVector SLA = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["SLA"]);
  NumericVector Al2As = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["Al2As"]);
  NumericVector Ar2Al;
  NumericVector RLD;
  if(transpirationMode=="Granier") Ar2Al = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["Ar2Al"]);
  else RLD = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["RLD"]);
  NumericVector WoodDensity = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["WoodDensity"]);
  NumericVector LeafDensity = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["LeafDensity"]);
  NumericVector FineRootDensity = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["FineRootDensity"]);
  NumericVector SRL = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["SRL"]);
  NumericVector conduit2sapwood = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["conduit2sapwood"]);
  
  //Growth parameters
  DataFrame paramsGrowth = Rcpp::as<Rcpp::DataFrame>(x["paramsGrowth"]);
  NumericVector WoodC = Rcpp::as<Rcpp::NumericVector>(paramsGrowth["WoodC"]);
  NumericVector RERleaf = Rcpp::as<Rcpp::NumericVector>(paramsGrowth["RERleaf"]);
  NumericVector RERsapwood = Rcpp::as<Rcpp::NumericVector>(paramsGrowth["RERsapwood"]);
  NumericVector RERfineroot = Rcpp::as<Rcpp::NumericVector>(paramsGrowth["RERfineroot"]);
  NumericVector CCleaf = Rcpp::as<Rcpp::NumericVector>(paramsGrowth["CCleaf"]);
  NumericVector CCsapwood = Rcpp::as<Rcpp::NumericVector>(paramsGrowth["CCsapwood"]);
  NumericVector CCfineroot = Rcpp::as<Rcpp::NumericVector>(paramsGrowth["CCfineroot"]);
  NumericVector RGRleafmax = Rcpp::as<Rcpp::NumericVector>(paramsGrowth["RGRleafmax"]);
  NumericVector RGRcambiummax = Rcpp::as<Rcpp::NumericVector>(paramsGrowth["RGRcambiummax"]);
  NumericVector RGRsapwoodmax = Rcpp::as<Rcpp::NumericVector>(paramsGrowth["RGRsapwoodmax"]);
  NumericVector RGRfinerootmax = Rcpp::as<Rcpp::NumericVector>(paramsGrowth["RGRfinerootmax"]);
  NumericVector SRsapwood = Rcpp::as<Rcpp::NumericVector>(paramsGrowth["SRsapwood"]);
  NumericVector SRfineroot = Rcpp::as<Rcpp::NumericVector>(paramsGrowth["SRfineroot"]);
  NumericVector RSSG = Rcpp::as<Rcpp::NumericVector>(paramsGrowth["RSSG"]);

  //Mortality/regeneration parameters
  DataFrame paramsMortalityRegeneration = Rcpp::as<Rcpp::DataFrame>(x["paramsMortalityRegeneration"]);
  NumericVector MortalityBaselineRate = Rcpp::as<Rcpp::NumericVector>(paramsMortalityRegeneration["MortalityBaselineRate"]);
  NumericVector SurvivalModelStep = Rcpp::as<Rcpp::NumericVector>(paramsMortalityRegeneration["SurvivalModelStep"]);
  NumericVector SurvivalB0 = Rcpp::as<Rcpp::NumericVector>(paramsMortalityRegeneration["SurvivalB0"]);
  NumericVector SurvivalB1 = Rcpp::as<Rcpp::NumericVector>(paramsMortalityRegeneration["SurvivalB1"]);
  NumericVector RecrTreeDensity = Rcpp::as<Rcpp::NumericVector>(paramsMortalityRegeneration["RecrTreeDensity"]);
  NumericVector RecrTreeDBH = Rcpp::as<Rcpp::NumericVector>(paramsMortalityRegeneration["RecrTreeDBH"]);
  NumericVector IngrowthTreeDensity = Rcpp::as<Rcpp::NumericVector>(paramsMortalityRegeneration["IngrowthTreeDensity"]);
  NumericVector IngrowthTreeDBH = Rcpp::as<Rcpp::NumericVector>(paramsMortalityRegeneration["IngrowthTreeDBH"]);
  
  //Phenology parameters
  DataFrame paramsPhenology = Rcpp::as<Rcpp::DataFrame>(x["paramsPhenology"]);
  CharacterVector phenoType = Rcpp::as<Rcpp::CharacterVector>(paramsPhenology["PhenologyType"]);
  NumericVector leafDuration = Rcpp::as<Rcpp::NumericVector>(paramsPhenology["LeafDuration"]);
  
  //Transpiration parameters
  DataFrame paramsTransp = Rcpp::as<Rcpp::DataFrame>(x["paramsTranspiration"]);
  NumericVector Psi_Extract, Kmax_stemxylem, Plant_kmax, VCleaf_kmax, VCleaf_c, VCleaf_d;
  NumericVector VCstem_kmax, VCstem_c, VCstem_d, VCroot_kmaxVEC, VCroot_c, VCroot_d, VGrhizo_kmaxVEC;
  NumericVector WUE_par(numCohorts, 0.3643);
  VCleaf_c = Rcpp::as<Rcpp::NumericVector>(paramsTransp["VCleaf_c"]);
  VCleaf_d = Rcpp::as<Rcpp::NumericVector>(paramsTransp["VCleaf_d"]);
  VCstem_c = Rcpp::as<Rcpp::NumericVector>(paramsTransp["VCstem_c"]);
  VCstem_d = Rcpp::as<Rcpp::NumericVector>(paramsTransp["VCstem_d"]);
  if(transpirationMode=="Granier") {
    Psi_Extract = Rcpp::as<Rcpp::NumericVector>(paramsTransp["Psi_Extract"]);
    if(paramsTransp.containsElementNamed("WUE_par")) {
      WUE_par = Rcpp::as<Rcpp::NumericVector>(paramsTransp["WUE_par"]);
    }
    if(paramsTransp.containsElementNamed("WUE_decay")) { //For compatibility with previous versions (2.7.5)
      WUE_par = Rcpp::as<Rcpp::NumericVector>(paramsTransp["WUE_decay"]);
    }
  } else {
    Kmax_stemxylem = paramsTransp["Kmax_stemxylem"];
    Plant_kmax= paramsTransp["Plant_kmax"];
    VCleaf_kmax = Rcpp::as<Rcpp::NumericVector>(paramsTransp["VCleaf_kmax"]);
    VCstem_kmax = paramsTransp["VCstem_kmax"];
    VCroot_kmaxVEC= paramsTransp["VCroot_kmax"];
    VCroot_c = paramsTransp["VCroot_c"];
    VCroot_d = paramsTransp["VCroot_d"];
    VGrhizo_kmaxVEC= paramsTransp["VGrhizo_kmax"];
  }

  //Water storage parameters
  DataFrame paramsWaterStorage = Rcpp::as<Rcpp::DataFrame>(x["paramsWaterStorage"]);
  NumericVector StemPI0 = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["StemPI0"]);
  NumericVector StemEPS = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["StemEPS"]);
  NumericVector StemAF = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["StemAF"]);
  NumericVector Vsapwood = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["Vsapwood"]); //l·m-2 = mm
  NumericVector LeafAF = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["LeafAF"]);
  NumericVector Vleaf = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["Vleaf"]); //l·m-2 = mm
  
  //Allometry parameters
  DataFrame paramsAllometries = Rcpp::as<Rcpp::DataFrame>(x["paramsAllometries"]);
  NumericVector Abt  = paramsAllometries["Abt"];
  NumericVector Bbt  = paramsAllometries["Bbt"];
  NumericVector BTsh  = paramsAllometries["BTsh"];
  
  //Ring of forming vessels
  // List ringList = as<Rcpp::List>(x["internalRings"]);

  //Daily output vectors
  DataFrame labileCarbonBalance = as<DataFrame>(modelOutput["LabileCarbonBalance"]);
  NumericVector GrossPhotosynthesis = labileCarbonBalance["GrossPhotosynthesis"];
  NumericVector MaintenanceRespiration = labileCarbonBalance["MaintenanceRespiration"];
  NumericVector GrowthCosts = labileCarbonBalance["GrowthCosts"];
  NumericVector RootExudation = labileCarbonBalance["RootExudation"];
  NumericVector LabileCarbonBalance = labileCarbonBalance["LabileCarbonBalance"];
  NumericVector PlantSugarLeaf = labileCarbonBalance["SugarLeaf"];
  NumericVector PlantStarchLeaf = labileCarbonBalance["StarchLeaf"];
  NumericVector PlantSugarSapwood = labileCarbonBalance["SugarSapwood"];
  NumericVector PlantStarchSapwood = labileCarbonBalance["StarchSapwood"];
  NumericVector PlantSugarTransport = labileCarbonBalance["SugarTransport"];
  for(int c=0;c<numCohorts;c++) {
    GrossPhotosynthesis[c] = 0.0;
    MaintenanceRespiration[c] = 0.0;
    GrowthCosts[c] = 0.0;
    RootExudation[c] = 0.0;
    LabileCarbonBalance[c] = 0.0;
    PlantSugarLeaf[c] = 0.0;
    PlantStarchLeaf[c] = 0.0;
    PlantSugarSapwood[c] = 0.0;
    PlantStarchSapwood[c] = 0.0;
    PlantSugarTransport[c] = 0.0;
  }
  
  DataFrame plantStructure = as<DataFrame>(modelOutput["PlantStructure"]);
  NumericVector LeafBiomass = plantStructure["LeafBiomass"];
  NumericVector SapwoodBiomass = plantStructure["SapwoodBiomass"];
  NumericVector FineRootBiomass = plantStructure["FineRootBiomass"];
  NumericVector LeafArea = plantStructure["LeafArea"];
  NumericVector SapwoodArea = plantStructure["SapwoodArea"];
  NumericVector FineRootArea = plantStructure["FineRootArea"];
  NumericVector HuberValue = plantStructure["HuberValue"];
  NumericVector RootAreaLeafArea = plantStructure["RootAreaLeafArea"];
  NumericVector OutputDBH = plantStructure["DBH"];
  NumericVector OutputHeight = plantStructure["Height"];
  
  DataFrame growthMortality = as<DataFrame>(modelOutput["GrowthMortality"]);
  NumericVector SAgrowth = growthMortality["SAgrowth"];
  NumericVector LAgrowth = growthMortality["LAgrowth"];
  NumericVector FRAgrowth = growthMortality["FRAgrowth"];
  NumericVector StarvationRate = growthMortality["StarvationRate"];
  NumericVector DessicationRate = growthMortality["DessicationRate"];
  NumericVector MortalityRate = growthMortality["MortalityRate"];
  
  for(int c=0;c<numCohorts;c++) {
    SAgrowth[c] = 0.0;
    LAgrowth[c] = 0.0;
    FRAgrowth[c] = 0.0;
    StarvationRate[c] = 0.0;
    DessicationRate[c] = 0.0;
    MortalityRate[c] = 0.0;
    LeafBiomass[c] = 0.0;
    SapwoodBiomass[c] = 0.0;
    FineRootBiomass[c] = 0.0;
    LeafArea[c] =  0.0;
    SapwoodArea[c] = 0.0;
    FineRootArea[c] = 0.0;
    HuberValue[c] = 0.0;
    RootAreaLeafArea[c] = 0.0;
  }
  
  NumericVector deltaLAgrowth(numCohorts,0.0);
  NumericVector deltaSAgrowth(numCohorts,0.0);
  
  //Stand carbon balance
  double standGrossPrimaryProduction = 0.0;
  double standMaintenanceRespiration = 0.0;
  double standSynthesisRespiration = 0.0;
  
  double equilibriumLeafSugarConc = equilibriumLeafTotalConc - nonSugarConcentration;
  double equilibriumSapwoodSugarConc = equilibriumSapwoodTotalConc - nonSugarConcentration;

  double rcellmax = relative_expansion_rate(0.0 ,30.0, -1.0, 0.5, 0.05, 5.0);
  
  //Initial Biomass balance
  List initialFinalCC = internalCommunication["initialFinalCC"];
  DataFrame ccIni = as<DataFrame>(initialFinalCC["ccIni_g_ind"]);
  fillCarbonCompartments(ccIni, x, "g_ind");
  NumericVector LeafBiomassBalance(numCohorts,0.0), FineRootBiomassBalance(numCohorts,0.0);
  
  DataFrame plantBiomassBalance = as<DataFrame>(modelOutput["PlantBiomassBalance"]);
  fillInitialPlantBiomassBalance(plantBiomassBalance, ccIni, above);
  
  NumericVector Volume_leaves = Rcpp::as<Rcpp::NumericVector>(ccIni["LeafStorageVolume"]);
  NumericVector Volume_sapwood = Rcpp::as<Rcpp::NumericVector>(ccIni["SapwoodStorageVolume"]);
  NumericVector Starch_max_leaves = Rcpp::as<Rcpp::NumericVector>(ccIni["LeafStarchMaximumConcentration"]);
  NumericVector Starch_max_sapwood = Rcpp::as<Rcpp::NumericVector>(ccIni["SapwoodStarchMaximumConcentration"]);
  NumericVector LeafStructBiomass = Rcpp::as<Rcpp::NumericVector>(ccIni["LeafStructuralBiomass"]);
  NumericVector SapwoodLivingStructBiomass = Rcpp::as<Rcpp::NumericVector>(ccIni["SapwoodLivingStructuralBiomass"]);
  NumericVector TotalLivingBiomass = Rcpp::as<Rcpp::NumericVector>(ccIni["TotalLivingBiomass"]);
  NumericVector TotalBiomass = Rcpp::as<Rcpp::NumericVector>(ccIni["TotalBiomass"]);


  //For survival model based on basal area  
  double treeBasalArea = 0.0;
  for(int j=0;j<numCohorts;j++){
    if(!NumericVector::is_na(DBH[j])) treeBasalArea += N[j]*3.141593*pow(DBH[j]/200,2.0);
  }

  //3. Carbon balance, growth, senescence and mortality by cohort
  for(int j=0;j<numCohorts;j++){
    if(N[j] > 0.0) {
      double costPerLA = 1000.0*CCleaf[j]/SLA[j]; // Construction cost in g gluc · m-2 of leaf area
      double costPerSA = CCsapwood[j]*sapwoodStructuralBiomass(1.0, H[j], L(j,_),V(j,_),WoodDensity[j]); // Construction cost in g gluc · cm-2 of sapwood area
      NumericVector deltaFRBgrowth(nlayers, 0.0);
      
      double LAexpanded = leafArea(LAI_expanded[j], N[j]);
      double LAlive = leafArea(LAI_live[j], N[j]);
      double LAdead = leafArea(LAI_dead[j], N[j]);

      double minimumStarchForSecondaryGrowth = Starch_max_sapwood[j]*RSSG[j];
      double minimumStarchForPrimaryGrowth = Starch_max_sapwood[j]*0.1;
        
      double leafAgG = 0.0;
      double leafRespDay = 0.0;
      double sapwoodResp = 0.0;
      double finerootResp = 0.0;
      double synthesisRespLA = 0.0;
      double synthesisRespSA = 0.0;
      double synthesisRespFRB = 0.0;   
      
      //Estimate phloem conductance as a factor of stem conductance
      double k_phloem = NA_REAL;
      if(subdailyCarbonBalance) k_phloem = VCstem_kmax[j]*phloemConductanceFactor*(0.018/1000.0);
      
      //Xylogenesis
      // List ring = ringList[j];
      double rleafcell = NA_REAL, rcambiumcell = NA_REAL;
      NumericVector rfineroot(nlayers);
      double relative_hormone_factor = std::max(0.0, std::min(1.0, LAI_expanded[j]/LAI_nocomp[j]));
      if(transpirationMode=="Granier") {
        // grow_ring(ring, PlantPsi[j] ,tday, 10.0);
        rleafcell = std::min(rcellmax, relative_expansion_rate(PlantPsi[j] ,tday, -1.0, 0.5,0.05,5.0));
        rcambiumcell = std::min(rcellmax, relative_hormone_factor*relative_expansion_rate(PlantPsi[j] ,tday, -1.0, 0.5,0.05,5.0));
        for(int l=0;l<nlayers;l++) rfineroot[l] = std::min(rcellmax, relative_expansion_rate(psiSoil[l] ,tday, -1.0 ,0.5,0.05,5.0));
        // if(j==0) Rcout<<j<< " Psi:"<< PlantPsi[j]<< " r:"<< rcambiumcell<<"\n";
      } else {
        rleafcell = std::min(rcellmax, relative_expansion_rate(psiRootCrown[j] ,tcan_day, -1.0, 0.5,0.05,5.0));
        rcambiumcell = std::min(rcellmax, relative_hormone_factor*relative_expansion_rate(psiRootCrown[j] ,tcan_day, -1.0, 0.5,0.05,5.0));
        for(int l=0;l<nlayers;l++) rfineroot[l] = std::min(rcellmax, relative_expansion_rate(RhizoPsi(j,l) ,Tsoil[l], -1.0, 0.5,0.05,5.0));
        // Rcout<<j<< " rcellmax "<< rcellmax<<" tcan_day "<< tcan_day<< " Psi:"<< psiRootCrown[j]<< " r:"<< rcambiumcell<<"\n";
      }
      if(!subdailyCarbonBalance) {
        //MAINTENANCE RESPIRATION
        //Respiratory biomass (g dw · ind-1)
        double leafSugarMass = sugarLeaf[j]*(Volume_leaves[j]*glucoseMolarMass);
        double sapwoodSugarMass = sugarSapwood[j]*(Volume_sapwood[j]*glucoseMolarMass);
        double B_resp_leaves = LeafStructBiomass[j] + leafSugarMass;
        double B_resp_sapwood = SapwoodLivingStructBiomass[j] + sapwoodSugarMass;
        // Rcout<<j<< " maintenance costs of leaf sugars: "<< (leafSugarMass/B_resp_leaves)<<" sapwood sugars: "<< (sapwoodSugarMass/B_resp_sapwood)<<"\n";
        double B_resp_fineroots = fineRootBiomass[j];
        double QR = qResp(tday);
        if(LAexpanded>0.0) {
          leafRespDay = B_resp_leaves*RERleaf[j]*QR*std::min(1.0, pow(PARcohort[j]/100.0, 0.5)); //Reduction under shade
        }
        sapwoodResp = B_resp_sapwood*RERsapwood[j]*QR;
        finerootResp = B_resp_fineroots*RERfineroot[j]*QR*(LAexpanded/LAlive);
        MaintenanceRespiration[j] += (leafRespDay+sapwoodResp+finerootResp)/TotalLivingBiomass[j]; 

        //PHOTOSYNTHESIS
        if(LAexpanded>0.0) {
          standGrossPrimaryProduction += Ag[j]; //Add GPP in gC · m-2
          
          //gross fotosynthesis
          double leafAgC = Ag[j]/(N[j]/10000.0); //Translate g C · m-2 to g C ·ind-1
          leafAgG = leafAgC*(glucoseMolarMass/(carbonMolarMass*6.0)); // from g C·ind-1 to g gluc · ind-1
          
          //Update output values
          GrossPhotosynthesis[j] = leafAgG/TotalLivingBiomass[j]; //Ag in g gluc · gdry-1 
        }
        
        //GROWTH
        double growthCostLA = 0.0;
        double growthCostSA = 0.0;
        double growthCostFRB = 0.0;   
        
        
        if(leafUnfolding[j]) {
          double deltaLApheno = std::max(leafAreaTarget[j]*(1.0 - StemPLC[j]) - LAexpanded, 0.0);
          double deltaLAsink = std::min(deltaLApheno, (crownBudPercent[j]/100.0)*SA[j]*RGRleafmax[j]*(rleafcell/rcellmax));
          if(!sinkLimitation) deltaLAsink = std::min(deltaLApheno, (crownBudPercent[j]/100.0)*SA[j]*RGRleafmax[j]); //Deactivates temperature and turgor limitation
          double deltaLAavailable = 0.0;
          deltaLAavailable = std::max(0.0, (starchSapwood[j]-minimumStarchForPrimaryGrowth)*(glucoseMolarMass*Volume_sapwood[j])/costPerLA);
          deltaLAgrowth[j] = std::min(deltaLAsink, deltaLAavailable);
          growthCostLA = deltaLAgrowth[j]*costPerLA;
          synthesisRespLA = growthCostLA*(CCleaf[j] - 1.0)/CCleaf[j];
        }
        
        //fine root growth
        if(fineRootBiomass[j] < fineRootBiomassTarget[j]) {
          for(int s = 0;s<nlayers;s++) {
            double deltaFRBpheno = std::max(fineRootBiomassTarget[j] - fineRootBiomass[j], 0.0);
            double deltaFRBsink = (V(j,s)*fineRootBiomass[j])*RGRfinerootmax[j]*(rfineroot[s]/rcellmax);
            if(!sinkLimitation) deltaFRBsink = (V(j,s)*fineRootBiomass[j])*RGRfinerootmax[j]; //Deactivates temperature and turgor limitation
            double deltaFRBavailable = std::max(0.0,(starchSapwood[j]-minimumStarchForPrimaryGrowth)*(glucoseMolarMass*Volume_sapwood[j])/CCfineroot[j]);
            deltaFRBgrowth[s] = std::min(deltaFRBpheno, std::min(deltaFRBsink, deltaFRBavailable));
            growthCostFRB += deltaFRBgrowth[s]*CCfineroot[j];
            synthesisRespFRB += growthCostFRB*(CCfineroot[j] - 1.0)/CCfineroot[j];
          }
        }
        
        if(LAexpanded>0.0) {
          // NumericVector SAring = ring["SA"];
          // double deltaSAring = 0.0;
          // if(SAring.size()==1) deltaSAring = SAring[0];
          // else deltaSAring = SAring[SAring.size()-1] - SAring[SAring.size()-2];
          // double cellfileareamaxincrease = 950.0; //Found empirically with T = 30 degrees and Psi = -0.033
          // double rgrcellfile = (deltaSAring/10.0)/cellfileareamaxincrease;
          double deltaSAsink = NA_REAL;
          if(!NumericVector::is_na(DBH[j])) { //Trees
            deltaSAsink = (3.141592*DBH[j]*RGRcambiummax[j]*(rcambiumcell/rcellmax)); 
            if(!sinkLimitation) deltaSAsink = 3.141592*DBH[j]*RGRcambiummax[j]; //Deactivates temperature and turgor limitation
          } else { // Shrubs
            deltaSAsink = (SA[j]*RGRsapwoodmax[j]*(rcambiumcell/rcellmax)); 
            if(!sinkLimitation) deltaSAsink = SA[j]*RGRsapwoodmax[j]; //Deactivates temperature and turgor limitation
          }
          double deltaSAavailable = std::max(0.0, (starchSapwood[j]-minimumStarchForSecondaryGrowth)*(glucoseMolarMass*Volume_sapwood[j])/costPerSA);
          deltaSAgrowth[j] = std::min(deltaSAsink, deltaSAavailable);
          // Rcout<< rcambiumcell<<" deltaSAsink "<< deltaSAsink << " dSAgrowth "<< deltaSAgrowth<<"\n"; 
          growthCostSA = deltaSAgrowth[j]*costPerSA; //increase cost (may be non zero if leaf growth was charged onto sapwood)
          synthesisRespSA = growthCostSA*(CCsapwood[j] - 1.0)/CCsapwood[j];
        }
        
        // Rcout<<growthCostLA<<" "<<growthCostSA<<" "<<growthCostFRB<< " "<< TotalLivingBiomass[j]<<"\n";
        GrowthCosts[j] +=(growthCostLA + growthCostSA + growthCostFRB)/TotalLivingBiomass[j]; //growth cost in g gluc · gdry-1
        
        //PARTIAL CARBON BALANCE
        double leafSugarMassDelta = leafAgG - leafRespDay;
        double sapwoodSugarMassDelta =  - sapwoodResp - finerootResp; 
        double sapwoodStarchMassDelta =  - growthCostFRB - growthCostLA - growthCostSA;
        
        sugarSapwood[j] += sapwoodSugarMassDelta/(Volume_sapwood[j]*glucoseMolarMass);
        starchSapwood[j] += sapwoodStarchMassDelta/(Volume_sapwood[j]*glucoseMolarMass);
        if(LAexpanded>0.0) sugarLeaf[j] += leafSugarMassDelta/(Volume_leaves[j]*glucoseMolarMass);
        
        //PHLOEM TRANSPORT AND SUGAR-STARCH DYNAMICS     
        if(LAexpanded>0.0) {
          double ff = (sugarLeaf[j]-sugarSapwood[j])/2.0; 
          sugarLeaf[j] -=ff;
          PlantSugarTransport[j] = 1000.0*(ff*Volume_leaves[j])/(3600.0*24.0); //mmol · s-1
          sugarSapwood[j] +=(Volume_leaves[j]/Volume_sapwood[j])*ff;
          double conversionLeaf = std::max(-starchLeaf[j], sugarLeaf[j] - equilibriumLeafSugarConc);
          starchLeaf[j] +=conversionLeaf;
          sugarLeaf[j] -=conversionLeaf;
        }
        double conversionSapwood = std::max(-starchSapwood[j], sugarSapwood[j] - equilibriumSapwoodSugarConc);
        starchSapwood[j] +=conversionSapwood;
        sugarSapwood[j] -=conversionSapwood;
      } else {
        //Get output matrices
        List labileCBInst = modelOutput["LabileCarbonBalanceInst"];
        NumericMatrix LabileCarbonBalanceInst = labileCBInst["LabileCarbonBalance"];
        NumericMatrix GrossPhotosynthesisInst = labileCBInst["GrossPhotosynthesis"];
        NumericMatrix MaintenanceRespirationInst = labileCBInst["MaintenanceRespiration"];
        NumericMatrix GrowthCostsInst = labileCBInst["GrowthCosts"];
        NumericMatrix RootExudationInst = labileCBInst["RootExudation"];
        NumericMatrix PlantSugarTransportInst = labileCBInst["SugarTransport"];
        NumericMatrix PlantSugarLeafInst = labileCBInst["SugarLeaf"];
        NumericMatrix PlantStarchLeafInst = labileCBInst["StarchLeaf"];
        NumericMatrix PlantSugarSapwoodInst = labileCBInst["SugarSapwood"]; 
        NumericMatrix PlantStarchSapwoodInst = labileCBInst["StarchSapwood"];
        
        //3.1 Carbon balance and growth by steps
        for(int s=0;s<numSteps;s++) {
          
          //Transform sugar concentration (mol gluc · l-1) to sugar mass (g gluc)
          // double leafSugarMassStep = sugarLeaf[j]*(Volume_leaves[j]*glucoseMolarMass);
          double sapwoodSugarMassStep = sugarSapwood[j]*(Volume_sapwood[j]*glucoseMolarMass);
          
          //LEAF PHOTOSYNTHESIS and RESPIRATION
          double leafRespStep = 0.0;
          double leafAgStepG = 0.0, leafAnStepG=0.0;
          if(LAexpanded>0.0){
            standGrossPrimaryProduction += AgStep(j,s); //Add GPP in gC · m-2
            
            double leafAgStepC = AgStep(j,s)/(N[j]/10000.0); //Translate g C · m-2 · h-1 to g C · h-1
            double leafAnStepC = AnStep(j,s)/(N[j]/10000.0); //Translate g C · m-2 · h-1 to g C · h-1
            leafAgStepG = leafAgStepC*(glucoseMolarMass/(carbonMolarMass*6.0)); // from g C· h-1 to g gluc · h-1
            leafAnStepG = leafAnStepC*(glucoseMolarMass/(carbonMolarMass*6.0)); // from g C· h-1 to g gluc · h-1
            leafAgG += leafAgStepG;
            leafRespStep = leafAgStepG - leafAnStepG; //Respiration as Ag - An
            GrossPhotosynthesisInst(j,s) = leafAgStepG/TotalLivingBiomass[j]; //Ag in g gluc · gdry-1
            GrossPhotosynthesis[j] += GrossPhotosynthesisInst(j,s); 
          }
          
          //MAINTENANCE RESPIRATION
          // double B_resp_leaves = LeafStructBiomass[j] + leafSugarMassStep;
          double B_resp_sapwood = SapwoodLivingStructBiomass[j] + sapwoodSugarMassStep;
          double B_resp_fineroots = fineRootBiomass[j];
          double QR = qResp(Tcan[s]);
          leafRespDay +=leafRespStep;
          double sapwoodRespStep = B_resp_sapwood*RERsapwood[j]*QR/((double) numSteps);
          sapwoodResp += sapwoodRespStep;
          double finerootRespStep = B_resp_fineroots*RERfineroot[j]*QR*(LAexpanded/LAlive)/((double) numSteps);
          finerootResp += finerootRespStep;
          MaintenanceRespirationInst(j,s) = (leafRespStep+sapwoodRespStep+finerootRespStep)/TotalLivingBiomass[j];//Rm in g gluc· gdry-1
          MaintenanceRespiration[j] += MaintenanceRespirationInst(j,s); 
          
          
          //GROWTH        
          double growthCostLAStep = 0.0;
          double growthCostSAStep = 0.0;
          double growthCostFRBStep = 0.0;   
          
          //Leaf growth
          if(leafUnfolding[j]) {
            double deltaLApheno = std::max(leafAreaTarget[j]*(1.0 - StemPLC[j]) - LAexpanded, 0.0);
            double deltaLAsink = std::min(deltaLApheno, (crownBudPercent[j]/100.0)*(1.0/((double) numSteps))*SA[j]*RGRleafmax[j]*(rleafcell/rcellmax));
            if(!sinkLimitation) deltaLAsink = std::min(deltaLApheno, (crownBudPercent[j]/100.0)*(1.0/((double) numSteps))*SA[j]*RGRleafmax[j]); //Deactivates temperature and turgor limitation
            //Grow at expense of stem sugar
            double deltaLAavailable = std::max(0.0, (starchSapwood[j]-minimumStarchForPrimaryGrowth)*(glucoseMolarMass*Volume_sapwood[j])/costPerLA);
            double deltaLAgrowthStep = std::min(deltaLAsink, deltaLAavailable);
            growthCostLAStep += deltaLAgrowthStep*costPerLA;
            deltaLAgrowth[j] += deltaLAgrowthStep;
            synthesisRespLA += deltaLAgrowthStep*(CCleaf[j] - 1.0)/CCleaf[j];
          }
          //fine root growth
          if(fineRootBiomass[j] < fineRootBiomassTarget[j]) {
            for(int l = 0;l<nlayers;l++) {
              double deltaFRBpheno = std::max(fineRootBiomassTarget[j] - fineRootBiomass[j], 0.0);
              double deltaFRBsink = (1.0/((double) numSteps))*(V(j,l)*fineRootBiomass[j])*RGRfinerootmax[j]*(rfineroot[l]/rcellmax);
              if(!sinkLimitation) deltaFRBsink = (1.0/((double) numSteps))*(V(j,l)*fineRootBiomass[j])*RGRfinerootmax[j]; //Deactivates temperature and turgor limitation
              double deltaFRBavailable = std::max(0.0, (starchSapwood[j]-minimumStarchForPrimaryGrowth)*(glucoseMolarMass*Volume_sapwood[j])/CCfineroot[j]);
              double deltaFRBgrowthStep = std::min(deltaFRBpheno, std::min(deltaFRBsink, deltaFRBavailable));
              growthCostFRBStep += deltaFRBgrowthStep*CCfineroot[j];
              deltaFRBgrowth[l] += deltaFRBgrowthStep;
              synthesisRespFRB += deltaFRBgrowthStep*(CCfineroot[j] - 1.0)/CCfineroot[j];
            }
          }
          //sapwood area growth
          if(LAexpanded>0.0) {
            // List ring = ringList[j];
            // NumericVector SAring = ring["SA"];
            // double deltaSAring = 0.0;
            // if(SAring.size()==1) deltaSAring = SAring[0];
            // else deltaSAring = SAring[SAring.size()-1] - SAring[SAring.size()-2];
            // double cellfileareamaxincrease = 950.0; //Found empirically with T = 30 degrees and Psi = -0.033
            // double rgrcellfile = (deltaSAring/10.0)/cellfileareamaxincrease;
            double deltaSAsink = NA_REAL;
            if(!NumericVector::is_na(DBH[j])) { //Trees
              deltaSAsink = (3.141592*DBH[j]*RGRcambiummax[j]*(rcambiumcell/rcellmax))/((double) numSteps); 
              if(!sinkLimitation) deltaSAsink = 3.141592*DBH[j]*RGRcambiummax[j]/((double) numSteps); //Deactivates temperature and turgor limitation
            } else { // Shrubs
              deltaSAsink = (SA[j]*RGRsapwoodmax[j]*(rcambiumcell/rcellmax))/((double) numSteps); 
              if(!sinkLimitation) deltaSAsink = SA[j]*RGRsapwoodmax[j]/((double) numSteps); //Deactivates temperature and turgor limitation
            }
            double deltaSAavailable = std::max(0.0, (starchSapwood[j] - minimumStarchForSecondaryGrowth)*(glucoseMolarMass*Volume_sapwood[j])/costPerSA);
            double deltaSAgrowthStep = std::min(deltaSAsink, deltaSAavailable);
            if(deltaSAgrowthStep<0.0) {
              Rcout<<deltaSAsink<<" "<< deltaSAavailable<< " "<< starchSapwood[j]<<"\n";
              stop("negative growth!"); 
            }
            growthCostSAStep += deltaSAgrowthStep*costPerSA; //increase cost (may be non zero if leaf growth was charged onto sapwood)
            deltaSAgrowth[j]  +=deltaSAgrowthStep;
            synthesisRespSA += growthCostSAStep*(CCsapwood[j] - 1.0)/CCsapwood[j];
          }
          GrowthCostsInst(j,s) += (growthCostLAStep + growthCostSAStep + growthCostFRBStep)/TotalLivingBiomass[j];
          GrowthCosts[j] +=GrowthCostsInst(j,s); //growth cost in g gluc · gdry-1
          
          //PHLOEM TRANSPORT AND SUGAR-STARCH DYNAMICS (INCLUDING PARTIAL MASS BALANCE)
          //sugar mass balance
          double leafSugarMassDeltaStep = leafAgStepG - leafRespStep;
          double sapwoodSugarMassDeltaStep = - finerootRespStep - sapwoodRespStep;
          double sapwoodStarchMassDeltaStep = - growthCostFRBStep - growthCostLAStep - growthCostSAStep;
          double ff = 0.0;
          double ctl = 3600.0*Volume_leaves[j]*glucoseMolarMass;
          double cts = 3600.0*Volume_sapwood[j]*glucoseMolarMass;
          for(int t=0;t<3600;t++) {
            sugarSapwood[j] += sapwoodSugarMassDeltaStep/cts;
            starchSapwood[j] += sapwoodStarchMassDeltaStep/cts;
            // double conversionSapwood = sugarStarchDynamicsStem(sugarSapwood[j]/StemSympRWCInst(j,s), starchSapwood[j]/StemSympRWCInst(j,s), equilibriumSapwoodSugarConc);
            double conversionSapwood = sugarStarchDynamicsStem(sugarSapwood[j], starchSapwood[j], equilibriumSapwoodSugarConc);
            // if(j==2) Rcout<<" coh:"<<j<< " s:"<<s<< " Lsugar: "<< sugarSapwood[j] << " Lstarch: "<< starchSapwood[j]<<" starch formation: "<<conversionSapwood<< "\n";
            // double starchSapwoodIncrease = conversionSapwood*StemSympRWCInst(j,s);
            double starchSapwoodIncrease = conversionSapwood;
            
            starchSapwood[j] += starchSapwoodIncrease;
            
            if(LAexpanded>0.0) {
              sugarLeaf[j] += leafSugarMassDeltaStep/ctl;
              double ft = phloemFlow(LeafSympPsiInst(j,s), StemSympPsiInst(j,s), sugarLeaf[j], sugarSapwood[j], Tcan[s], k_phloem, nonSugarConcentration)*LAlive; //flow as mol glucose per s
              // sugar-starch dynamics
              double conversionLeaf = sugarStarchDynamicsLeaf(sugarLeaf[j], starchLeaf[j], equilibriumLeafSugarConc);
              double starchLeafIncrease = conversionLeaf;
              starchLeaf[j]  += starchLeafIncrease;
              // Rcout<<" coh:"<<j<< " s:"<<s<< " Ssugar: "<< sugarSapwood[j] << " Sstarch: "<< starchSapwood[j]<<" starch formation: "<<conversionSapwood<< "\n";
              //Apply phloem transport (mol gluc) to sugar concentrations (mol gluc· l-1)
              sugarLeaf[j]  +=  (-ft/Volume_leaves[j]) - starchLeafIncrease;
              sugarSapwood[j] +=  (ft/Volume_sapwood[j]) - starchSapwoodIncrease;
              ff +=ft;
            } else {
              sugarSapwood[j] += (- starchSapwoodIncrease);
            }
          }
          //Divert to root exudation if sapwood starch is over maximum capacity
          if(starchSapwood[j] > Starch_max_sapwood[j]) {
            RootExudationInst(j,s) += ((starchSapwood[j] - Starch_max_sapwood[j])*(Volume_sapwood[j]*glucoseMolarMass)/TotalLivingBiomass[j]);
            starchSapwood[j] = Starch_max_sapwood[j];
          }
          

          //Add instantaneous root exudation to daily root exudation
          RootExudation[j] += RootExudationInst(j,s);
          
          //Instantaneous carbon balance
          LabileCarbonBalanceInst(j,s) = GrossPhotosynthesisInst(j,s) - MaintenanceRespirationInst(j,s) - GrowthCostsInst(j,s) - RootExudationInst(j,s);
          LabileCarbonBalance[j] +=LabileCarbonBalanceInst(j,s);
          
          PlantSugarLeafInst(j,s) = sugarLeaf[j];
          PlantSugarSapwoodInst(j,s) = sugarSapwood[j];
          PlantStarchLeafInst(j,s) = starchLeaf[j];
          PlantStarchSapwoodInst(j,s) = starchSapwood[j];
          PlantSugarTransportInst(j,s) = 1000.0*ff/(3600.0); //mmol·s-1
          PlantSugarTransport[j] += ff/((double) numSteps); //Average daily rate To calculate daily phloem balance (positive means towards stem)
          // Rcout<<" coh:"<<j<< " s:"<<s<< " conc leaf: "<< sugarLeaf[j] << " conc sap: "<< sugarSapwood[j]<<" ff: "<<ff<< "\n";
        }
      }

      double leafBiomassIncrement = deltaLAgrowth[j]*(1000.0/SLA[j]);
      double finerootBiomassIncrement = sum(deltaFRBgrowth);
      
      //add maintenance and synthesis respiration as g C·m-2
      standMaintenanceRespiration += (leafRespDay+sapwoodResp+finerootResp)*(N[j]/10000.0)*((carbonMolarMass*6.0)/glucoseMolarMass);
      standSynthesisRespiration += (synthesisRespSA+synthesisRespLA+synthesisRespFRB)*(N[j]/10000.0)*((carbonMolarMass*6.0)/glucoseMolarMass);
      

      //SENESCENCE
      //Leaf senescence
      double propLeafSenescence = 0.0;
      //Leaf senescence due to age (Ca+ accumulation) only in evergreen species
      if(phenoType[j] == "progressive-evergreen") {
        propLeafSenescence = std::min(1.0,(LAexpanded/(365.25*LAlive*leafDuration[j])));
      }
      else if((phenoType[j] == "oneflush-evergreen") && (leafSenescence[j])) {
        propLeafSenescence = std::min(1.0,(LAexpanded/(LAlive*leafDuration[j]))); // Fraction of old leaves that die
        leafSenescence[j] = false; //To prevent further loss
      }
      else if(((phenoType[j] == "winter-deciduous") || (phenoType[j] == "winter-semideciduous")) && leafSenescence[j]) {
        propLeafSenescence = 1.0;
        leafSenescence[j] = false; //To prevent further loss
      }

      //Leaf senescence and bud senescence due to drought (only when PLC increases)
      double PLCinc = (StemPLC[j]-StemPLCprev[j]);
      if(PLCinc>0.0) {
        double LeafPDEF = proportionDefoliationWeibull(psiApoLeaf[j], VCleaf_c[j], VCleaf_d[j]);
        double LApdef = std::min(LAexpanded, (1.0 - LeafPDEF)*leafAreaTarget[j]);
        if(LApdef<LAexpanded) {
          propLeafSenescence = std::max((LAexpanded-LApdef)/LAexpanded, propLeafSenescence);
        }
        double budplc = 100.0*(1.0-StemPLC[j]);
        if(budplc < crownBudPercent[j]) {
          crownBudPercent[j] = budplc;
        }
      }

      double deltaLAsenescence = std::min(LAexpanded, LAexpanded*propLeafSenescence);
      double senescenceLeafLoss = deltaLAsenescence*(1000.0/SLA[j]);

      //Define sapwood senescence as maximum of turnover and sapwood exceeding the target
      double propSASenescence = SRsapwood[j]*std::max(0.0,(tday-5.0)/20.0)/(1.0+15.0*exp(-0.01*H[j]));
      double deltaSASenescence = std::max(0.0, SA[j] - sapwoodAreaTarget[j]);
      propSASenescence = std::max(propSASenescence, deltaSASenescence/SA[j]);
        
      //FRB SENESCENCE
      NumericVector deltaFRBsenescence(nlayers, 0.0);
      for(int l=0;l<nlayers;l++) {
        double daySenescence = NA_REAL;
        if(transpirationMode=="Granier") daySenescence = SRfineroot[j]*std::max(0.0,(tday-5.0)/20.0);
        else daySenescence = SRfineroot[j]*std::max(0.0,(Tsoil[l]-5.0)/20.0);
        deltaFRBsenescence[l] = fineRootBiomass[j]*V(j,l)*daySenescence;
      }
      double senescenceFinerootLoss = sum(deltaFRBsenescence);
      
      // if(j==(numCohorts-1)) Rcout<< j << " before translocation "<< sugarLeaf[j]<< " "<< starchLeaf[j]<<"\n";
      

      //TRANSLOCATION (in mol gluc) of labile carbon
      // Rcout<<"-translocation";
      double translocationSugarLeaf = propLeafSenescence*Volume_leaves[j]*sugarLeaf[j];
      double translocationStarchLeaf = propLeafSenescence*Volume_leaves[j]*starchLeaf[j];
      double translocationSugarSapwood = propSASenescence*Volume_sapwood[j]*sugarSapwood[j];
      if(Volume_leaves[j]>0) {
        if(starchLeaf[j] > Starch_max_leaves[j]) { // Add excess leaf starch to translocation
          translocationStarchLeaf += ((starchLeaf[j] - Starch_max_leaves[j])*Volume_leaves[j]);
        }
        sugarLeaf[j] = ((sugarLeaf[j]*Volume_leaves[j]) - translocationSugarLeaf)/Volume_leaves[j]; 
        starchLeaf[j] = ((starchLeaf[j]*Volume_leaves[j]) - translocationStarchLeaf)/Volume_leaves[j]; 
      }
      sugarSapwood[j] = ((sugarSapwood[j]*Volume_sapwood[j]) - translocationSugarSapwood)/Volume_sapwood[j]; 
      starchSapwood[j] = ((starchSapwood[j]*Volume_sapwood[j]) + translocationSugarLeaf + translocationStarchLeaf + translocationSugarSapwood)/Volume_sapwood[j]; 

      //ROOT EXUDATION and close carbon balance (non-subdaily carbon balance)
      if(!subdailyCarbonBalance) {
        //Excess sapwood starch carbon is lost as root exudation
        if(starchSapwood[j] > Starch_max_sapwood[j]) {
          RootExudation[j] += ((starchSapwood[j] - Starch_max_sapwood[j])*(Volume_sapwood[j]*glucoseMolarMass)/TotalLivingBiomass[j]);
          starchSapwood[j] = Starch_max_sapwood[j];
        }
        //Labile CARBON balance
        LabileCarbonBalance[j] = GrossPhotosynthesis[j] - MaintenanceRespiration[j] - GrowthCosts[j] - RootExudation[j];
      }
      
        
      //UPDATE LEAF AREA, SAPWOOD AREA, FINE ROOT BIOMASS AND CONCENTRATION IN LABILE POOLS
      // Rcout<<"-update";
      // double LAprev = LAexpanded;
      LAexpanded += deltaLAgrowth[j] - deltaLAsenescence;
      if(LAexpanded < 0.0) {
        deltaLAsenescence -= LAexpanded;
        LAexpanded = 0.0;
      }
      LAdead += deltaLAsenescence;
      LAI_dead[j] = LAdead*N[j]/10000.0;
      LAI_expanded[j] = LAexpanded*N[j]/10000.0;
      // double SAprev = SA[j];
      SA[j] = SA[j] + deltaSAgrowth[j] - deltaSASenescence; 
      NumericVector newFRB(nlayers,0.0);
      for(int s=0;s<nlayers;s++) {
        newFRB[s] = fineRootBiomass[j]*V(j,s) + deltaFRBgrowth[s] - deltaFRBsenescence[s];
      }
      fineRootBiomass[j] = sum(newFRB);
      for(int s=0;s<nlayers;s++) { 
        V(j,s) = newFRB[s]/fineRootBiomass[j];
      }

      
      //UPDATE DERIVED QUANTITIES (individual level)   
      // Rcout<<"-updatederived";
      if(transpirationMode=="Granier") {
        //Update Huber value
        // if(LAlive>0.0) {
        //   Al2As[j] = (LAlive)/(SA[j]/10000.0);
        // }
        for(int c=0;c<numCohorts;c++){
          L(c,_) = coarseRootLengths(V(c,_), widths, 0.5); //Arbitrary ratio (to revise some day)
          CRSV[c] = coarseRootSoilVolume(V(c,_), widths, 0.5);
        }
      } else { //SPERRY/Sureau
        if(LAlive>0.0) {
          //Update Huber value, stem and root hydraulic conductance
          double oldstemR = 1.0/VCstem_kmax[j];
          double oldrootR = 1.0/VCroot_kmaxVEC[j];
          double oldrootprop = oldrootR/(oldrootR+oldstemR);
          
          // Al2As[j] = (LAlive)/(SA[j]/10000.0);
          VCstem_kmax[j]=maximumStemHydraulicConductance(Kmax_stemxylem[j], Hmed[j], Al2As[j] ,H[j], taper); 

          //Update rhizosphere maximum conductance
          NumericVector VGrhizo_new = rhizosphereMaximumConductance(Ksat, newFRB, LAI_live[j], N[j],
                                                                    SRL[j], FineRootDensity[j], RLD[j]);
          for(int s=0;s<nlayers;s++) { 
            VGrhizo_kmax(j,s) = VGrhizo_new[s];
          }
          VGrhizo_kmaxVEC[j] = sum(VGrhizo_kmax(j,_));
          
          //Update root maximum conductance so that it keeps the same resistance proportion with stem conductance
          double newstemR = 1.0/VCstem_kmax[j];
          double newrootR = oldrootprop*newstemR/(1.0-oldrootprop);
          VCroot_kmaxVEC[j] = 1.0/newrootR;
          //Update coarse root soil volume
          CRSV[j] = coarseRootSoilVolumeFromConductance(Kmax_stemxylem[j], VCroot_kmaxVEC[j], Al2As[j],
                                                        V(j,_), widths, rfc);
          //Update coarse root length and root maximum conductance
          L(j,_) = coarseRootLengthsFromVolume(CRSV[j], V(j,_), widths, rfc);
          NumericVector xp = rootxylemConductanceProportions(L(j,_), V(j,_));
          VCroot_kmax(j,_) = VCroot_kmaxVEC[j]*xp;
          //Update Plant_kmax
          Plant_kmax[j] = 1.0/((1.0/VCleaf_kmax[j])+(1.0/VCstem_kmax[j])+(1.0/VCroot_kmaxVEC[j]));
        }
      }
      //Decrease PLC due to new SA growth
      StemPLC[j] = std::max(0.0, StemPLC[j] - (deltaSAgrowth[j]/SA[j])); 
      LeafPLC[j] = std::max(0.0, LeafPLC[j] - (deltaLAgrowth[j]/LAexpanded)); 
      
      //Increase crown buds to new SA growth
      crownBudPercent[j] = std::min(100.0, crownBudPercent[j] + 100.0*(deltaSAgrowth[j]/SA[j]));
      
      
      //LEAF/FINE ROOT BIOMASS balance (g_ind)
      LeafBiomassBalance[j] = leafBiomassIncrement - senescenceLeafLoss;
      FineRootBiomassBalance[j] = finerootBiomassIncrement - senescenceFinerootLoss;
      
      //MORTALITY Death by carbon starvation or dessication
      double Ndead_day = 0.0;
      bool dynamicCohort = true;
      bool isShrub = !NumericVector::is_na(Cover[j]);
      if((!shrubDynamics) && isShrub) dynamicCohort = false;
      double stemSympRWC = NA_REAL;
      if(transpirationMode=="Granier") stemSympRWC = symplasticRelativeWaterContent(PlantPsi[j], StemPI0[j], StemEPS[j]);
      else stemSympRWC = symplasticRelativeWaterContent(psiSympStem[j], StemPI0[j], StemEPS[j]);
      if(NumericVector::is_na(stemSympRWC)) stop("Missing value for stem symp RWC");

      //Sapwood sugar relative to equilibrium, indicator of starvation
      double relativeSugarSapwood = (sugarSapwood[j]/equilibriumSapwoodSugarConc);
      if(dynamicCohort) {
        String cause = "undertermined";
        //Determine fire severity if fire occurred
        double LAI_burnt_change = 0.0;
        bool abovegroundFireSurvival = true;
        if(fireOccurrence) {
          double rho_air = meteoland::utils_airDensity(tmax, Patm);
          double foliar_factor = leafThermalFactor(SLA[j]);
          //Determine foliage/bud burn
          double Ib_surf = fireBehavior["I_b_surface [kW/m]"];
          double t_res_surf = fireBehavior["t_r_surface [s]"];
          double t_r_crown = fireBehavior["t_r_crown [s]"];
          double fm_dead = fireBehavior["DFMC [%]"];
          double Ic_ratio = fireBehavior["Ic_ratio"];
          double Hn_leaves = 100.0*necrosisHeight(Ib_surf, t_res_surf, foliar_factor, tmax, rho_air); //Necrosis height (cm)
          double Hn_buds = 100.0*necrosisHeight(Ib_surf, t_res_surf, 0.130, tmax, rho_air); //Bud necrosis height (cm)
          double cbh = H[j]*(1.0 - CR[j]);
          double burnRatioLeaves = leafAreaProportion(0.0, Hn_leaves, cbh, H[j]);
          double burnRatioBuds = leafAreaProportion(0.0, Hn_buds, cbh, H[j]);
          // Rcout << " foliar_factor "<< foliar_factor << " Hn_leaves "<<Hn_leaves << " br_leaves "<< burnRatioLeaves<< " Hn_buds "<<Hn_buds << " br_buds "<< burnRatioBuds<<"\n";
          //Determine crown fire or torching effects
          double canopyFMC = (LFMC[j]*(1.0 - StemPLC[j]) + fm_dead*StemPLC[j]);
          double Ib_crit = criticalFirelineIntensity(cbh/100.0, canopyFMC);
          // Rcout << "Ic_ratio "<< Ic_ratio <<" Ib_crit "<<Ib_crit<< " Ib_surf "<< Ib_surf<<"\n";
          if((Ic_ratio > 1.0) || (Ib_surf > Ib_crit)) {
            burnRatioLeaves = 1.0;
            double Tc = necrosisCriticalTemperature(t_r_crown, 0.130 , tmax, rho_air);
            if(Tc < 900.0) burnRatioBuds = 1.0;
          }
          // Rcout << "br_leaves "<< burnRatioLeaves<< " br_buds "<< burnRatioBuds<<"\n";
          //Surface fire effects on cambium
          double bark_diff = barkThermalDiffusivity(fm_dead, 500.0, tmax);
          double xn = radialBoleNecrosis(Ib_surf, t_res_surf, bark_diff, tmax, rho_air);
          double bark_thickness = 1.0;
          if(isShrub) {
            bark_thickness = BTsh[j]*0.1; // from mm to cm
          } else {
            bark_thickness = Abt[j]*pow(DBH[j],Bbt[j])*0.1;// from mm to cm
          } 
          // Rcout << "xn "<< xn<< " xa "<<bark_thickness<<"\n";
          
          //Effects
          double LAburned = LAexpanded * burnRatioLeaves;
          LAI_burnt_change = LAburned*N[j]/10000.0;
          crownBudPercent[j] = crownBudPercent[j]*(1.0 - burnRatioBuds);
          abovegroundFireSurvival = (burnRatioBuds < 1.0) && (bark_thickness > xn);
          // Rcout << "abovegroundSurvival "<< abovegroundFireSurvival<< " burnRatioLeaves "<< burnRatioLeaves<< " LAI_burnt_change "<<LAI_burnt_change<<"\n";
        }
        if(abovegroundFireSurvival) {
          if(mortalityMode=="whole-cohort/deterministic") {
            if((relativeSugarSapwood < mortalityRelativeSugarThreshold) && allowStarvation) {
              Ndead_day = N[j];
              if(verbose) Rcout<<" [Cohort "<< j<<" died from starvation] ";
              cause = "starvation";
            } else if( (stemSympRWC < mortalityRWCThreshold) && allowDessication) {
              Ndead_day = N[j];
              if(verbose) Rcout<<" [Cohort "<< j<<" died from dessication] ";
              cause = "dessication";
            }
          } else {
            //Daily basal mortality rate based on constant year probability
            double basalMortalityRate = 1.0 - exp(log(1.0 - MortalityBaselineRate[j])/356.0);
            //If survival model is available, replace basal mortality value
            if(!NumericVector::is_na(SurvivalModelStep[j]) && !NumericVector::is_na(SurvivalB0[j]) && !NumericVector::is_na(SurvivalB1[j])) {
              //Probability of dying in model step years
              double lp  = SurvivalB0[j] + SurvivalB1[j]*sqrt(treeBasalArea);
              double Pmodel = 1.0 - exp(lp)/(1.0 + exp(lp));
              //Probability of dying in 1 year
              double Pmodel1yr = (1.0- exp(log(1.0-Pmodel)/SurvivalModelStep[j]));
              //Daily basal mortality rate based on model
              basalMortalityRate = 1.0 - exp(log(1.0 - Pmodel1yr)/356.0);
            }
            if(allowStarvation) StarvationRate[j] = dailyMortalityProbability(relativeSugarSapwood, mortalityRelativeSugarThreshold);
            if(allowDessication) DessicationRate[j] = dailyMortalityProbability((stemSympRWC + (1.0 - StemPLC[j]))/2.0, mortalityRWCThreshold);
            MortalityRate[j] = max(NumericVector::create(basalMortalityRate, DessicationRate[j],  StarvationRate[j]));
            if((DessicationRate[j] > basalMortalityRate) && (DessicationRate[j] > StarvationRate[j])) {
              cause = "dessication";
            } else if((StarvationRate[j] > basalMortalityRate) && (StarvationRate[j] > DessicationRate[j])) {
              cause = "starvation";
            }
            if(NumericVector::is_na(MortalityRate[j])) {
              Rcout<< " Basal mortality rate " << basalMortalityRate << " Dessication rate " << DessicationRate[j] << " starvation rate "<< StarvationRate[j]<<"\n"; 
              stop("Missing value for mortality rate");
            }
            // Rcout<< j << " "<< stemSympRWC<< " "<< DessicationRate[j]<<"\n";
            if(mortalityMode =="density/deterministic") {
              Ndead_day = N[j]*MortalityRate[j];
            } else if(mortalityMode =="whole-cohort/stochastic") {
              if(R::runif(0.0,1.0) < MortalityRate[j]) {
                Ndead_day = N[j];
                Rcout<<" [Cohort "<< j<<" died from " << cause.get_cstring() << "] ";
              }
            } else if(mortalityMode == "density/stochastic") {
              Ndead_day = R::rbinom(round(N[j]), MortalityRate[j]);
            }
          }
        } else { //Cohort burned
          Ndead_day = N[j];
          cause = "burnt";
        }
        // Update density and increase the number of dead plants
        if(NumericVector::is_na(Ndead_day)) stop("Missing value for Ndead_day");
        Ndead_day = std::min(Ndead_day, N[j]);
        double Cdead_day = Cover[j]*(Ndead_day/N[j]);
        if(cause == "starvation") {
          N_starvation[j] = N_starvation[j] + Ndead_day;
        } else if(cause == "burnt") {
          N_burnt[j] = N_burnt[j] + Ndead_day;
        } else if(cause == "dessication") {
          N_dessication[j] = N_dessication[j] + Ndead_day;
        } else if(!isShrub) { // Self-thinning occurring in tree cohorts
          if(DBH[j] < IngrowthTreeDBH[j]) {
            double b_st = log(RecrTreeDensity[j]/IngrowthTreeDensity[j])/log(RecrTreeDBH[j]/IngrowthTreeDBH[j]);
            double a_st = IngrowthTreeDensity[j]/pow(IngrowthTreeDBH[j], b_st);
            double N_st = a_st*pow(DBH[j], b_st);
            double N_dead_selfthinning = N[j] - std::min(N[j], N_st);
            // Rcout<< b_st<< " "<< a_st<< " "<< N_st<< " "<< N_dead_selfthinning<<"\n";
            Ndead_day = Ndead_day + N_dead_selfthinning;
          }
        }
        N[j] = N[j] - Ndead_day;
        N_dead[j] = N_dead[j] + Ndead_day;
        if(isShrub) {
          Cover[j] = std::max(0.0, Cover[j] - Cdead_day);
          Cover_dead[j] = Cover_dead[j] + Cdead_day;
          if(cause == "starvation") {
            Cover_starvation[j] = Cover_starvation[j] + Cdead_day;
          } else if(cause == "dessication") {
            Cover_dessication[j] = Cover_dessication[j] + Cdead_day;
          } else if(cause == "burnt") {
            Cover_burnt[j] = Cover_burnt[j] + Cdead_day;
          }
        }
        //Update LAI dead and LAI expanded as a result of density decrease
        if(abovegroundFireSurvival) {
          double LAI_change = LAexpanded*Ndead_day/10000.0;
          LAI_dead[j] = LAI_dead[j] + LAI_change;
          LAI_expanded[j] = std::max(0.0, LAI_expanded[j] - LAI_change - LAI_burnt_change);
        } else {
          fineRootBiomassTarget[j] = 0.0;
          sapwoodAreaTarget[j] = 0.0;
          leafAreaTarget[j] = 0.0;
          LAI_live[j] = 0.0;
          LAI_dead[j] = 0.0;
          LAI_expanded[j] = 0.0;
        }
      }
      

      //UPDATE TARGETS
      // Rcout<<"-updatetargets";
      //Set target leaf area if bud formation is allowed
      // if(budFormation[j]) {
      //   if(allocationStrategy == "Plant_kmax") {
      //     leafAreaTarget[j] = LAlive*(Plant_kmax[j]/allocationTarget[j]);
      //   } else if(allocationStrategy =="Al2As") {
      //     leafAreaTarget[j] = (SA[j]/10000.0)*allocationTarget[j];
      //   }
      //   LAI_live[j] = leafAreaTarget[j]*N[j]/10000.0;
      // }
      //Update fine root biomass target     
      if((LAI_live[j]>0.0) && (N[j]>0.0)) {
        if(transpirationMode=="Granier") {
          sapwoodAreaTarget[j] = 10000.0*leafAreaTarget[j]/Al2As[j];
          fineRootBiomassTarget[j] = (Ar2Al[j]*leafAreaTarget[j])/(specificRootSurfaceArea(SRL[j], FineRootDensity[j])*1e-4);
        } else {
          if(allocationStrategy == "Plant_kmax") {
            sapwoodAreaTarget[j] = 10000.0*(leafAreaTarget[j]/Al2As[j])*(allocationTarget[j]/Plant_kmax[j]);
          } else if(allocationStrategy =="Al2As") {
            sapwoodAreaTarget[j] = 10000.0*leafAreaTarget[j]/Al2As[j];
          }
          NumericVector VGrhizo_target(nlayers,0.0);
          for(int s=0;s<nlayers;s++) {
            // Rcout<<VCroot_kmaxVEC[j]<< " "<<VCstem_kmax[j]<<"\n";
            VGrhizo_target[s] = V(j,s)*findRhizosphereMaximumConductance(averageFracRhizosphereResistance*100.0,
                                                   VG_n[s], VG_alpha[s],
                                                   VCroot_kmaxVEC[j], VCroot_c[j], VCroot_d[j],
                                                   VCstem_kmax[j], VCstem_c[j], VCstem_d[j],
                                                   VCleaf_kmax[j], VCleaf_c[j], VCleaf_d[j],
                                                   log(VGrhizo_kmax(j,s)));
          }
          fineRootBiomassTarget[j] = fineRootBiomassPerIndividual(Ksat, VGrhizo_target, LAI_live[j], N[j],
                                                                  SRL[j], FineRootDensity[j], RLD[j]);
          
        }
      }
      
      //Output variables
      // Rcout<<"-output";
      FineRootBiomass[j] = fineRootBiomass[j];
      SapwoodArea[j] = SA[j];
      FineRootArea[j] = fineRootBiomass[j]*specificRootSurfaceArea(SRL[j], FineRootDensity[j])*1e-4;
      SapwoodBiomass[j] = sapwoodStructuralBiomass(SA[j], H[j], L(j,_), V(j,_),WoodDensity[j]);
      LeafBiomass[j] = leafStructuralBiomass(LAI_expanded[j],N[j],SLA[j]);
      SAgrowth[j] += deltaSAgrowth[j]; //Store sapwood area growth rate (cm2/day)
      LAgrowth[j] += deltaLAgrowth[j];//Store Leaf area growth rate (m2/day)
      LeafArea[j] = LAexpanded;
      HuberValue[j] = SA[j]/leafAreaTarget[j]; 
      RootAreaLeafArea[j] = FineRootArea[j]/leafAreaTarget[j]; 
      FRAgrowth[j] = sum(deltaFRBgrowth)*specificRootSurfaceArea(SRL[j], FineRootDensity[j])*1e-4;//Store fine root area growth rate (m2·d-1)
    }
  }
  
  //UPDATE STRUCTURAL VARIABLES
  updateStructuralVariables(x, deltaSAgrowth);
  
  for(int j=0;j<numCohorts;j++){
    //RECALCULATE storage concentrations (SA, LA and H may have changed)
    double newVolumeSapwood = sapwoodStorageVolume(SA[j], H[j], L(j,_),V(j,_),WoodDensity[j], conduit2sapwood[j]);
    double newVolumeLeaves = leafStorageVolume(LAI_expanded[j],  N[j], SLA[j], LeafDensity[j]);
    if(newVolumeLeaves > 0.0) {
      sugarLeaf[j] = std::max(0.0, sugarLeaf[j]*(Volume_leaves[j]/newVolumeLeaves));
      starchLeaf[j] = std::max(0.0, starchLeaf[j]*(Volume_leaves[j]/newVolumeLeaves)); 
    } else {
      sugarLeaf[j] = 0.0;
      starchLeaf[j] = 0.0;
    }
    sugarSapwood[j] = std::max(0.0, sugarSapwood[j]*(Volume_sapwood[j]/newVolumeSapwood));
    starchSapwood[j] = std::max(0.0, starchSapwood[j]*(Volume_sapwood[j]/newVolumeSapwood)); 
    // if(j==(numCohorts-1)) Rcout<< j << " after recalculation "<< sugarLeaf[j]<< " "<< starchLeaf[j]<<"\n";
    
    //OUTPUT VARIABLES
    OutputDBH[j] = DBH[j];
    OutputHeight[j] = H[j];
    PlantSugarLeaf[j] = sugarLeaf[j];
    PlantStarchLeaf[j] = starchLeaf[j];
    PlantSugarSapwood[j] = sugarSapwood[j];
    PlantStarchSapwood[j] = starchSapwood[j];
  }

  //CLOSE BIOMASS BALANCE
  closePlantBiomassBalance(initialFinalCC, plantBiomassBalance, x,
                           LabileCarbonBalance, LeafBiomassBalance, FineRootBiomassBalance);

  //Update pool proportions and rhizosphere overlap
  if(plantWaterPools) {
    NumericVector poolProportions = Rcpp::as<Rcpp::NumericVector>(belowdf["poolProportions"]);
    for(int j=0;j<numCohorts;j++) poolProportions[j] = LAI_live[j]/sum(LAI_live);
    //Update RHOP
    List newRHOP;
    if(rhizosphereOverlap=="none") newRHOP = nonoverlapHorizontalProportions(V);
    else newRHOP = horizontalProportions(poolProportions, CRSV, N, V, widths, rfc);
    for(int j=0;j<numCohorts;j++) RHOP[j] = newRHOP[j];
  }
  
  
  NumericVector standCB = modelOutput["CarbonBalance"];
  standCB["GrossPrimaryProduction"] = standGrossPrimaryProduction;
  standCB["MaintenanceRespiration"] = standMaintenanceRespiration;
  standCB["SynthesisRespiration"] = standSynthesisRespiration;
  standCB["NetPrimaryProduction"] = standGrossPrimaryProduction - standMaintenanceRespiration - standSynthesisRespiration;

}


//' @rdname communication
//' @keywords internal
// [[Rcpp::export("growth_day_inner")]]
void growthDay_inner(List internalCommunication, List x, CharacterVector date, NumericVector meteovec, 
                double latitude, double elevation, double slope = NA_REAL, double aspect = NA_REAL,  
                double runon = 0.0, Nullable<NumericVector> lateralFlows = R_NilValue, double waterTableDepth = NA_REAL, 
                bool modifyInput = true) {
   
   double tmin = meteovec["MinTemperature"];
   double tmax = meteovec["MaxTemperature"];
   if(tmin > tmax) {
     warning("tmin > tmax. Swapping values.");
     double swap = tmin;
     tmin = tmax;
     tmax = swap;
   }
   double rhmin = meteovec["MinRelativeHumidity"];
   double rhmax = meteovec["MaxRelativeHumidity"];
   if(NumericVector::is_na(rhmax)) {
     rhmax = 100.0;
   }
   if(NumericVector::is_na(rhmin)) {
     double vp_tmin = meteoland::utils_saturationVP(tmin);
     double vp_tmax = meteoland::utils_saturationVP(tmax);
     rhmin = std::min(rhmax, 100.0*(vp_tmin/vp_tmax));
   }
   if(rhmin > rhmax) {
     warning("rhmin > rhmax. Swapping values.");
     double swap = rhmin;
     rhmin = rhmax;
     rhmax = swap;
   }
   double rad = meteovec["Radiation"];
   double prec = meteovec["Precipitation"];
   double wind = NA_REAL;
   if(meteovec.containsElementNamed("WindSpeed")) wind = meteovec["WindSpeed"];
   double Catm = NA_REAL; 
   if(meteovec.containsElementNamed("CO2")) Catm = meteovec["CO2"];
   double Patm = NA_REAL; 
   if(meteovec.containsElementNamed("Patm")) Patm = meteovec["Patm"];
   double Rint = NA_REAL; 
   if(meteovec.containsElementNamed("RainfallIntensity")) Rint = meteovec["RainfallIntensity"];
   double pfire = 0.0; 
   if(meteovec.containsElementNamed("FireProbability")) pfire = meteovec["FireProbability"];
   
   //Control parameters
   List control = x["control"];
   bool verbose = control["verbose"];
   
   bool leafPhenology = control["leafPhenology"];
   String transpirationMode = control["transpirationMode"];
   if(NumericVector::is_na(Catm)) Catm = control["defaultCO2"];
   
   //Will not modify input x 
   if(!modifyInput) {
     x = clone(x);
   }
   
   std::string c = as<std::string>(date[0]);
   int month = std::atoi(c.substr(5,2).c_str());
   int J = meteoland::radiation_julianDay(std::atoi(c.substr(0, 4).c_str()),std::atoi(c.substr(5,2).c_str()),std::atoi(c.substr(8,2).c_str()));
   double delta = meteoland::radiation_solarDeclination(J);
   double solarConstant = meteoland::radiation_solarConstant(J);
   double latrad = latitude * (M_PI/180.0);
   if(NumericVector::is_na(aspect)) aspect = 0.0;
   if(NumericVector::is_na(slope)) slope = 0.0;
   double asprad = aspect * (M_PI/180.0);
   double slorad = slope * (M_PI/180.0);
   double photoperiod = meteoland::radiation_daylength(latrad, 0.0, 0.0, delta);
   double tday = meteoland::utils_averageDaylightTemperature(tmin, tmax);
   double pet = meteoland::penman(latrad, elevation, slorad, asprad, J, tmin, tmax, rhmin, rhmax, rad, wind);
   //Derive doy from date  
   int J0101 = meteoland::radiation_julianDay(std::atoi(c.substr(0, 4).c_str()),1,1);
   int doy = J - J0101+1;
   if(NumericVector::is_na(wind)) wind = control["defaultWindSpeed"]; 
   if(wind<0.1) wind = 0.1; //Minimum windspeed abovecanopy
   
   NumericVector defaultRainfallIntensityPerMonth = control["defaultRainfallIntensityPerMonth"];
   if(NumericVector::is_na(Rint)) Rint = rainfallIntensity(month, prec, defaultRainfallIntensityPerMonth);
   
   //Update phenology
   if(leafPhenology) {
     updatePhenology(x, doy, photoperiod, tday);
     updateLeaves(x, wind, true);
   }
   
   NumericVector meteovec_inner = NumericVector::create(
     Named("tday") = tday,
     Named("tmin") = tmin, 
     Named("tmax") = tmax,
     Named("tminPrev") = tmin, 
     Named("tmaxPrev") = tmax, 
     Named("tminNext") = tmin, 
     Named("prec") = prec,
     Named("rhmin") = rhmin, 
     Named("rhmax") = rhmax, 
     Named("rad") = rad, 
     Named("wind") = wind, 
     Named("Catm") = Catm,
     Named("Patm") = Patm,
     Named("pet") = pet,
     Named("rint") = Rint,
     Named("pfire") = pfire);
   growthDay_private(internalCommunication, x, meteovec_inner, 
                     latitude, elevation, slope, aspect,
                     solarConstant, delta, 
                     runon, lateralFlows, waterTableDepth,
                     verbose);
 }

//' @rdname spwb_day
// [[Rcpp::export("growth_day")]]
List growthDay(List x, CharacterVector date, NumericVector meteovec, 
               double latitude, double elevation, double slope = NA_REAL, double aspect = NA_REAL,  
               double runon = 0.0, Nullable<NumericVector> lateralFlows = R_NilValue, double waterTableDepth = NA_REAL, 
               bool modifyInput = true) {
  //Instance communication structures
  List internalCommunication = instanceCommunicationStructures(x, "growth");
  
  growthDay_inner(internalCommunication, x, date, meteovec,
                                     latitude, elevation, slope, aspect,
                                     runon, lateralFlows, waterTableDepth,
                                     modifyInput);
  
  List modelOutput = copyModelOutput(internalCommunication, x, "growth");
  return(modelOutput);
}


void checkgrowthInput(List x, String transpirationMode, String soilFunctions) {
  
  DataFrame soil = Rcpp::as<Rcpp::DataFrame>(x["soil"]);
  if(!x.containsElementNamed("above")) stop("above missing in growthInput");
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  if(!above.containsElementNamed("LAI_live")) stop("LAI_live missing in growthInput$above");
  if(!above.containsElementNamed("LAI_expanded")) stop("LAI_expanded missing in growthInput$above");
  if(!above.containsElementNamed("LAI_dead")) stop("LAI_dead missing in growthInput$above");
  if(!above.containsElementNamed("SA")) stop("SA missing in growthInput$above");
  if(!above.containsElementNamed("CR")) stop("CR missing in growthInput$above");
  if(!above.containsElementNamed("H")) stop("H missing in growthInput$above");
  if(!above.containsElementNamed("N")) stop("N missing in growthInput$above");
  if(!above.containsElementNamed("DBH")) stop("DBH missing in growthInput$above");
  
  if(!x.containsElementNamed("below")) stop("below missing in growthInput");
  DataFrame below = Rcpp::as<Rcpp::DataFrame>(x["below"]);
  if(!below.containsElementNamed("Z50")) stop("Z50 missing in growthInput$below");
  if(!below.containsElementNamed("Z95")) stop("Z95 missing in growthInput$below");
  if(!x.containsElementNamed("belowLayers")) stop("belowLayers missing in growthInput");
  List belowLayers = Rcpp::as<Rcpp::List>(x["belowLayers"]);
  if(!belowLayers.containsElementNamed("V")) stop("V missing in growthInput$belowLayers");
  if((transpirationMode=="Sperry") || (transpirationMode=="Sureau")) {
    if(!belowLayers.containsElementNamed("VGrhizo_kmax")) stop("VGrhizo_kmax missing in growthInput$below");
    if(!belowLayers.containsElementNamed("VCroot_kmax")) stop("VCroot_kmax missing in growthInput$below");
  }  
  
  if(!x.containsElementNamed("paramsPhenology")) stop("paramsPhenology missing in growthInput");
  DataFrame paramsPhenology = Rcpp::as<Rcpp::DataFrame>(x["paramsPhenology"]);
  if(!paramsPhenology.containsElementNamed("Sgdd")) stop("Sgdd missing in paramsPhenology");
  
  if(!x.containsElementNamed("paramsInterception")) stop("paramsInterception missing in growthInput");
  DataFrame paramsInterception = Rcpp::as<Rcpp::DataFrame>(x["paramsInterception"]);
  if(!paramsInterception.containsElementNamed("kPAR")) stop("kPAR missing in growthInput$paramsInterception");
  if(!paramsInterception.containsElementNamed("g")) stop("g missing in growthInput$paramsInterception");
  
  if(!x.containsElementNamed("paramsGrowth")) stop("paramsGrowth missing in growthInput");
  DataFrame paramsGrowth = Rcpp::as<Rcpp::DataFrame>(x["paramsGrowth"]);
  if(!paramsGrowth.containsElementNamed("WoodC")) stop("WoodC missing in growthInput$paramsGrowth");
  if(!paramsGrowth.containsElementNamed("RGRleafmax")) stop("RGRleafmax missing in growthInput$paramsGrowth");
  if(!paramsGrowth.containsElementNamed("RGRsapwoodmax")) stop("RGRsapwoodmax missing in growthInput$paramsGrowth");
  if(!paramsGrowth.containsElementNamed("RGRfinerootmax")) stop("RGRfinerootmax missing in growthInput$paramsGrowth");
  if(!paramsGrowth.containsElementNamed("RERleaf")) stop("RERleaf missing in growthInput$paramsGrowth");
  if(!paramsGrowth.containsElementNamed("RERsapwood")) stop("RERsapwood missing in growthInput$paramsGrowth");
  if(!paramsGrowth.containsElementNamed("RERfineroot")) stop("RERfineroot missing in growthInput$paramsGrowth");
  
  if(!x.containsElementNamed("paramsAnatomy")) stop("paramsAnatomy missing in growthInput");
  DataFrame paramsAnatomy = Rcpp::as<Rcpp::DataFrame>(x["paramsAnatomy"]);
  if(!paramsAnatomy.containsElementNamed("SLA")) stop("SLA missing in paramsAnatomy$paramsGrowth");
  if(!paramsAnatomy.containsElementNamed("Al2As")) stop("Al2As missing in paramsAnatomy$paramsGrowth");
  if(!paramsAnatomy.containsElementNamed("WoodDensity")) stop("WoodDensity missing in paramsAnatomy$paramsGrowth");
  
  if(!x.containsElementNamed("paramsTranspiration")) stop("paramsTranspiration missing in growthInput");
  DataFrame paramsTranspiration = Rcpp::as<Rcpp::DataFrame>(x["paramsTranspiration"]);
  if(transpirationMode=="Granier") {
    if(!paramsTranspiration.containsElementNamed("Psi_Extract")) stop("Psi_Extract missing in growthInput$paramsTransp");
    if(!paramsTranspiration.containsElementNamed("WUE")) stop("WUE missing in growthInput$paramsTransp");
  } else if((transpirationMode=="Sperry") || (transpirationMode=="Sureau")) {
    if(!soil.containsElementNamed("VG_n")) stop("VG_n missing in soil");
    if(!soil.containsElementNamed("VG_alpha")) stop("VG_alpha missing in soil");
    
    if(!paramsTranspiration.containsElementNamed("VCstem_kmax")) stop("VCstem_kmax missing in growthInput");
    if(!paramsTranspiration.containsElementNamed("VCstem_c")) stop("VCstem_c missing in growthInput$paramsTransp");
    if(!paramsTranspiration.containsElementNamed("VCstem_d")) stop("VCstem_d missing in growthInput$paramsTransp");
    if(!paramsTranspiration.containsElementNamed("VCroot_c")) stop("VCroot_c missing in growthInput$paramsTransp");
    if(!paramsTranspiration.containsElementNamed("VCroot_d")) stop("VCroot_d missing in growthInput$paramsTransp");
    if(!paramsTranspiration.containsElementNamed("Gswmax")) stop("Gswmax missing in growthInput$paramsTransp");
    if(!paramsTranspiration.containsElementNamed("Vmax298")) stop("Vmax298 missing in growthInput$paramsTransp");
    if(!paramsTranspiration.containsElementNamed("Jmax298")) stop("Jmax298 missing in growthInput$paramsTransp");
  }
  if(!soil.containsElementNamed("W")) stop("W missing in soil");
  if(!soil.containsElementNamed("widths")) stop("widths missing in soil");
  if(!soil.containsElementNamed("macro")) stop("macro missing in soil");
  if(soilFunctions=="SX") {
    if(!soil.containsElementNamed("clay")) stop("clay missing in soil");
    if(!soil.containsElementNamed("sand")) stop("sand missing in soil");
  }
  if(soilFunctions=="VG") {
    if(!soil.containsElementNamed("VG_n")) stop("VG_n missing in soil");
    if(!soil.containsElementNamed("VG_alpha")) stop("VG_alpha missing in soil");
    if(!soil.containsElementNamed("VG_theta_res")) stop("VG_theta_res missing in soil");
    if(!soil.containsElementNamed("VG_theta_sat")) stop("VG_theta_sat missing in soil");
  }
}

// [[Rcpp::export(".defineGrowthDailyOutput")]]
List defineGrowthDailyOutput(double latitude, double elevation, double slope, double aspect, 
                             CharacterVector dateStrings, List x) {
  
  NumericVector topo = NumericVector::create(elevation, slope, aspect);
  topo.attr("names") = CharacterVector::create("elevation", "slope", "aspect");
  
  List growthInput = clone(x);
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  int numCohorts = above.nrow();
  int numDays = dateStrings.size();
  List control = x["control"];
  DataFrame soil = Rcpp::as<Rcpp::DataFrame>(x["soil"]);
  String transpirationMode = control["transpirationMode"];
  
  DataFrame DWB = defineWaterBalanceDailyOutput(dateStrings, transpirationMode);
  List Soil = defineSoilDailyOutput(dateStrings, soil, true);
  DataFrame Snow = defineSnowDailyOutput(dateStrings);
  List sunlitDO = defineSunlitShadeLeavesDailyOutput(dateStrings, above);
  List shadeDO = defineSunlitShadeLeavesDailyOutput(dateStrings, above);
  List plantDWOL = definePlantWaterDailyOutput(dateStrings, above, soil, control);
  DataFrame Stand = defineStandDailyOutput(dateStrings);
  
  //Detailed subday results
  List subdailyRes(numDays);
  subdailyRes.attr("names") = dateStrings ;
  
  //Plant carbon output variables
  NumericMatrix LabileCarbonBalance(numDays, numCohorts);
  NumericMatrix MaintenanceRespiration(numDays, numCohorts);
  NumericMatrix GrowthCosts(numDays, numCohorts);
  NumericMatrix PlantSugarLeaf(numDays, numCohorts);
  NumericMatrix PlantStarchLeaf(numDays, numCohorts);
  NumericMatrix PlantSugarSapwood(numDays, numCohorts);
  NumericMatrix PlantStarchSapwood(numDays, numCohorts);
  NumericMatrix PlantSugarTransport(numDays, numCohorts);
  NumericMatrix SapwoodBiomass(numDays, numCohorts);
  NumericMatrix LeafBiomass(numDays, numCohorts);
  NumericMatrix SapwoodArea(numDays, numCohorts);
  NumericMatrix LeafArea(numDays, numCohorts);
  NumericMatrix FineRootArea(numDays, numCohorts);
  NumericMatrix FineRootBiomass(numDays, numCohorts);
  NumericMatrix HuberValue(numDays, numCohorts);
  NumericMatrix RootAreaLeafArea(numDays, numCohorts);
  NumericMatrix DBH(numDays, numCohorts);
  NumericMatrix Height(numDays, numCohorts);
  NumericMatrix LabileBiomass(numDays, numCohorts);
  NumericMatrix TotalBiomass(numDays, numCohorts);
  NumericMatrix SAgrowth(numDays, numCohorts), LAgrowth(numDays, numCohorts), FRAgrowth(numDays, numCohorts);
  NumericMatrix starvationRate(numDays, numCohorts), dessicationRate(numDays, numCohorts), mortalityRate(numDays, numCohorts);
  NumericMatrix GrossPhotosynthesis(numDays, numCohorts);
  NumericMatrix PlantLAIexpanded(numDays, numCohorts), PlantLAIdead(numDays, numCohorts), PlantLAIlive(numDays, numCohorts);
  NumericMatrix RootExudation(numDays, numCohorts);
  NumericMatrix StructuralBiomassBalance(numDays, numCohorts);
  NumericMatrix LabileBiomassBalance(numDays, numCohorts);
  NumericMatrix PlantBiomassBalance(numDays, numCohorts);
  NumericMatrix MortalityBiomassLoss(numDays, numCohorts);
  NumericMatrix CohortBiomassBalance(numDays, numCohorts);
  NumericMatrix StandBiomassBalance(numDays, 5);
  NumericMatrix StandCarbonBalance(numDays, 4);
  
  
  //Add matrix dimnames
  LabileCarbonBalance.attr("dimnames") = List::create(dateStrings, above.attr("row.names"));
  GrossPhotosynthesis.attr("dimnames") = List::create(dateStrings, above.attr("row.names"));
  MaintenanceRespiration.attr("dimnames") = List::create(dateStrings, above.attr("row.names"));
  GrowthCosts.attr("dimnames") = List::create(dateStrings, above.attr("row.names"));
  PlantSugarLeaf.attr("dimnames") = List::create(dateStrings, above.attr("row.names"));
  PlantStarchLeaf.attr("dimnames") = List::create(dateStrings, above.attr("row.names")) ;
  PlantSugarSapwood.attr("dimnames") = List::create(dateStrings, above.attr("row.names")) ;
  PlantStarchSapwood.attr("dimnames") = List::create(dateStrings, above.attr("row.names")) ;
  PlantSugarTransport.attr("dimnames") = List::create(dateStrings, above.attr("row.names")) ;
  SapwoodBiomass.attr("dimnames") = List::create(dateStrings, above.attr("row.names")) ;
  LeafBiomass.attr("dimnames") = List::create(dateStrings, above.attr("row.names")) ;
  FineRootArea.attr("dimnames") = List::create(dateStrings, above.attr("row.names")) ;
  SapwoodArea.attr("dimnames") = List::create(dateStrings, above.attr("row.names")) ;
  LeafArea.attr("dimnames") = List::create(dateStrings, above.attr("row.names")) ;
  HuberValue.attr("dimnames") = List::create(dateStrings, above.attr("row.names")) ;
  RootAreaLeafArea.attr("dimnames") = List::create(dateStrings, above.attr("row.names")) ;
  FineRootBiomass.attr("dimnames") = List::create(dateStrings, above.attr("row.names"));
  DBH.attr("dimnames") = List::create(dateStrings, above.attr("row.names"));
  Height.attr("dimnames") = List::create(dateStrings, above.attr("row.names"));
  LAgrowth.attr("dimnames") = List::create(dateStrings, above.attr("row.names")) ;
  SAgrowth.attr("dimnames") = List::create(dateStrings, above.attr("row.names")) ;
  FRAgrowth.attr("dimnames") = List::create(dateStrings, above.attr("row.names"));
  dessicationRate.attr("dimnames") = List::create(dateStrings, above.attr("row.names")) ;
  mortalityRate.attr("dimnames") = List::create(dateStrings, above.attr("row.names")) ;
  starvationRate.attr("dimnames") = List::create(dateStrings, above.attr("row.names"));
  RootExudation.attr("dimnames") = List::create(dateStrings, above.attr("row.names"));
  StructuralBiomassBalance.attr("dimnames") = List::create(dateStrings, above.attr("row.names"));
  LabileBiomassBalance.attr("dimnames") = List::create(dateStrings, above.attr("row.names"));
  PlantBiomassBalance.attr("dimnames") = List::create(dateStrings, above.attr("row.names"));
  MortalityBiomassLoss.attr("dimnames") = List::create(dateStrings, above.attr("row.names"));
  CohortBiomassBalance.attr("dimnames") = List::create(dateStrings, above.attr("row.names"));
  
  StandBiomassBalance.attr("dimnames") = List::create(dateStrings, 
                           CharacterVector::create("StructuralBalance", "LabileBalance", "PlantBalance", "MortalityLoss", "CohortBalance"));
  StandCarbonBalance.attr("dimnames") = List::create(dateStrings, 
                          CharacterVector::create("GrossPrimaryProduction", "MaintenanceRespiration", "SynthesisRespiration", "NetPrimaryProduction"));
  
  // Assemble output
  List labileCarbonBalance = List::create(
    Named("GrossPhotosynthesis") = GrossPhotosynthesis,
    Named("MaintenanceRespiration") = MaintenanceRespiration,
    Named("GrowthCosts") = GrowthCosts,
    Named("RootExudation") = RootExudation,
    Named("LabileCarbonBalance") = LabileCarbonBalance,
    Named("SugarLeaf") = PlantSugarLeaf,
    Named("StarchLeaf") = PlantStarchLeaf,
    Named("SugarSapwood") = PlantSugarSapwood,
    Named("StarchSapwood") = PlantStarchSapwood,
    Named("SugarTransport") = PlantSugarTransport
  );
  List plantBiomassBalance = List::create(_["StructuralBiomassBalance"] = StructuralBiomassBalance,
                                          _["LabileBiomassBalance"] = LabileBiomassBalance,
                                          _["PlantBiomassBalance"] = PlantBiomassBalance,
                                          _["MortalityBiomassLoss"] = MortalityBiomassLoss,
                                          _["CohortBiomassBalance"] = CohortBiomassBalance);

  
  List growthMortality, plantStructure;
  plantStructure = List::create(Named("LeafBiomass")=LeafBiomass,
                                Named("SapwoodBiomass") = SapwoodBiomass,
                                Named("FineRootBiomass") = FineRootBiomass,
                                Named("LeafArea") = LeafArea,
                                Named("SapwoodArea")=SapwoodArea,
                                Named("FineRootArea") = FineRootArea,
                                Named("HuberValue") = HuberValue,
                                Named("RootAreaLeafArea") = RootAreaLeafArea,
                                Named("DBH") = DBH,
                                Named("Height") = Height);
  growthMortality = List::create(Named("LAgrowth") = LAgrowth,
                                 Named("SAgrowth") = SAgrowth,
                                 Named("FRAgrowth") = FRAgrowth,
                                 Named("StarvationRate") = starvationRate,
                                 Named("DessicationRate") = dessicationRate,
                                 Named("MortalityRate") = mortalityRate);

  List l;
  if(transpirationMode=="Granier") {
    l = List::create(Named("latitude") = latitude,
                     Named("topography") = topo,
                     Named("weather") = NA_REAL,
                     Named("growthInput") = growthInput,
                     Named("growthOutput") = x,
                     Named("WaterBalance")= DWB, 
                     Named("CarbonBalance")=StandCarbonBalance, 
                     Named("BiomassBalance") = StandBiomassBalance);
    if(control["soilResults"]) l.push_back(Soil, "Soil");
    if(control["snowResults"]) l.push_back(Snow, "Snow");
    if(control["standResults"]) l.push_back(Stand, "Stand");
    if(control["plantResults"]) l.push_back(plantDWOL, "Plants");
    if(control["labileCarbonBalanceResults"]) l.push_back(labileCarbonBalance, "LabileCarbonBalance");
    l.push_back(plantBiomassBalance, "PlantBiomassBalance");
    if(control["plantStructureResults"]) l.push_back(plantStructure, "PlantStructure");
    if(control["growthMortalityResults"]) l.push_back(growthMortality, "GrowthMortality");
    if(control["fireHazardResults"]) {
      DataFrame fireHazard = defineFireHazardOutput(dateStrings);
      l.push_back(fireHazard, "FireHazard");
    }
  } else {
    
    DataFrame canopy = Rcpp::as<Rcpp::DataFrame>(x["canopy"]);
    DataFrame DEB = defineEnergyBalanceDailyOutput(dateStrings);
    DataFrame DT = defineTemperatureDailyOutput(dateStrings);
    NumericMatrix DLT =  defineTemperatureLayersDailyOutput(dateStrings, canopy);
    l = List::create(Named("latitude") = latitude,
                     Named("topography") = topo,
                     Named("weather") = NA_REAL,
                     Named("growthInput") = growthInput,
                     Named("growthOutput") = x,
                     Named("WaterBalance")= DWB, 
                     Named("CarbonBalance")=StandCarbonBalance, 
                     Named("EnergyBalance") = DEB,
                     Named("BiomassBalance") = StandBiomassBalance);
    if(control["temperatureResults"]) {
      l.push_back(DT, "Temperature");
      if(control["multiLayerBalance"]) l.push_back(DLT,"TemperatureLayers");
    }
    if(control["soilResults"]) l.push_back(Soil, "Soil");
    if(control["snowResults"]) l.push_back(Snow, "Snow");
    if(control["standResults"]) l.push_back(Stand, "Stand");
    if(control["plantResults"]) l.push_back(plantDWOL, "Plants");
    if(control["labileCarbonBalanceResults"]) l.push_back(labileCarbonBalance, "LabileCarbonBalance");
    l.push_back(plantBiomassBalance, "PlantBiomassBalance");
    if(control["plantStructureResults"]) l.push_back(plantStructure, "PlantStructure");
    if(control["growthMortalityResults"]) l.push_back(growthMortality, "GrowthMortality");
    if(control["leafResults"]) {
      l.push_back(sunlitDO, "SunlitLeaves");
      l.push_back(shadeDO, "ShadeLeaves");
    }
    if(control["fireHazardResults"]) {
      DataFrame fireHazard = defineFireHazardOutput(dateStrings);
      l.push_back(fireHazard, "FireHazard");
    }
  }
  if(control["subdailyResults"]) l.push_back(subdailyRes,"subdaily");
  
  l.attr("class") = CharacterVector::create("growth","list");
  return(l);
}


// [[Rcpp::export(".fillGrowthDailyOutput")]]
void fillGrowthDailyOutput(List l, List x, List sDay, int iday) {
  
  List control = x["control"];
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  int numCohorts = above.nrow();
  int ntimesteps = control["ndailysteps"];
  DataFrame canopyParams = Rcpp::as<Rcpp::DataFrame>(x["canopy"]);
  int ncanlayers = canopyParams.nrow(); //Number of canopy layers
  DataFrame soil = Rcpp::as<Rcpp::DataFrame>(x["soil"]);
  String transpirationMode = control["transpirationMode"];
  
  DataFrame DWB = Rcpp::as<Rcpp::DataFrame>(l["WaterBalance"]);
  int numDays = DWB.nrow();
  fillWaterBalanceDailyOutput(DWB, sDay, iday, transpirationMode);
  
  if(control["soilResults"]) {
    String soilFunctions = control["soilFunctions"];
    List Soil = Rcpp::as<Rcpp::List>(l["Soil"]);
    fillSoilDailyOutput(Soil, soil, sDay, 
                        iday, numDays, soilFunctions,
                        true);
  }
  if(control["snowResults"]) {
    DataFrame Snow = Rcpp::as<Rcpp::DataFrame>(l["Snow"]);
    fillSnowDailyOutput(Snow, x, iday);
  }
  if(control["standResults"]) {
    DataFrame Stand = Rcpp::as<Rcpp::DataFrame>(l["Stand"]);
    fillStandDailyOutput(Stand, sDay,iday); 
  }
  if(control["plantResults"]) {
    List plantDWOL = l["Plants"];
    fillPlantWaterDailyOutput(plantDWOL, sDay, iday, transpirationMode); 
    if(transpirationMode!= "Granier") {
      List sunlitDO = l["SunlitLeaves"];
      List shadeDO = l["ShadeLeaves"];
      fillSunlitShadeLeavesDailyOutput(sunlitDO, shadeDO, sDay, iday, numCohorts); 
    } 
  }
  if(transpirationMode!= "Granier") {
    List DEB = l["EnergyBalance"];
    fillEnergyBalanceDailyOutput(DEB,sDay, iday, ntimesteps);
    if(control["temperatureResults"]) {
      List DT = l["Temperature"];
      fillTemperatureDailyOutput(DT,sDay, iday, ntimesteps);
      if(control["multiLayerBalance"]) {
        NumericMatrix DLT = l["TemperatureLayers"];
        fillTemperatureLayersDailyOutput(DLT,sDay, iday, ncanlayers, ntimesteps);
      }
    }
  } 
  if(control["fireHazardResults"]) {
    DataFrame fireHazard = Rcpp::as<Rcpp::DataFrame>(l["FireHazard"]);
    fillFireHazardOutput(fireHazard, sDay, iday);
  }
  
  if(control["subdailyResults"]) {
    List subdailyRes = Rcpp::as<Rcpp::List>(l["subdaily"]);
    subdailyRes[iday] = copyAdvancedGROWTHOutput(sDay, x); //Clones subdaily results because they are communication structures
  }
  
  List sb = sDay["Soil"];
  List db = sDay["WaterBalance"];
  List Plants = sDay["Plants"];
  DataFrame gm = Rcpp::as<Rcpp::DataFrame>(sDay["GrowthMortality"]);
  DataFrame cb = Rcpp::as<Rcpp::DataFrame>(sDay["LabileCarbonBalance"]);
  DataFrame bb = Rcpp::as<Rcpp::DataFrame>(sDay["PlantBiomassBalance"]);
  DataFrame ps = Rcpp::as<Rcpp::DataFrame>(sDay["PlantStructure"]);
  

  if(control["labileCarbonBalanceResults"]) {
    List labileCarbonBalance = Rcpp::as<Rcpp::List>(l["LabileCarbonBalance"]);
    NumericMatrix LabileCarbonBalance = Rcpp::as<Rcpp::NumericMatrix>(labileCarbonBalance["LabileCarbonBalance"]);
    NumericMatrix GrossPhotosynthesis = Rcpp::as<Rcpp::NumericMatrix>(labileCarbonBalance["GrossPhotosynthesis"]);
    NumericMatrix GrowthCosts = Rcpp::as<Rcpp::NumericMatrix>(labileCarbonBalance["GrowthCosts"]);
    NumericMatrix RootExudation = Rcpp::as<Rcpp::NumericMatrix>(labileCarbonBalance["RootExudation"]);
    NumericMatrix MaintenanceRespiration = Rcpp::as<Rcpp::NumericMatrix>(labileCarbonBalance["MaintenanceRespiration"]);
    NumericMatrix PlantSugarLeaf = Rcpp::as<Rcpp::NumericMatrix>(labileCarbonBalance["SugarLeaf"]);
    NumericMatrix PlantStarchLeaf = Rcpp::as<Rcpp::NumericMatrix>(labileCarbonBalance["StarchLeaf"]);
    NumericMatrix PlantSugarSapwood = Rcpp::as<Rcpp::NumericMatrix>(labileCarbonBalance["SugarSapwood"]);
    NumericMatrix PlantStarchSapwood = Rcpp::as<Rcpp::NumericMatrix>(labileCarbonBalance["StarchSapwood"]);
    NumericMatrix PlantSugarTransport = Rcpp::as<Rcpp::NumericMatrix>(labileCarbonBalance["SugarTransport"]);
    NumericVector LabileCarbonBalanceIN = Rcpp::as<Rcpp::NumericVector>(cb["LabileCarbonBalance"]);
    NumericVector MaintenanceRespirationIN = Rcpp::as<Rcpp::NumericVector>(cb["MaintenanceRespiration"]);
    NumericVector GrowthCostsIN = Rcpp::as<Rcpp::NumericVector>(cb["GrowthCosts"]);
    NumericVector GrossPhotosynthesisIN = Rcpp::as<Rcpp::NumericVector>(cb["GrossPhotosynthesis"]);
    NumericVector PlantSugarLeafIN = Rcpp::as<Rcpp::NumericVector>(cb["SugarLeaf"]);
    NumericVector PlantStarchLeafIN = Rcpp::as<Rcpp::NumericVector>(cb["StarchLeaf"]);
    NumericVector PlantSugarSapwoodIN = Rcpp::as<Rcpp::NumericVector>(cb["SugarSapwood"]);
    NumericVector PlantStarchSapwoodIN = Rcpp::as<Rcpp::NumericVector>(cb["StarchSapwood"]);
    NumericVector PlantSugarTransportIN = Rcpp::as<Rcpp::NumericVector>(cb["SugarTransport"]);
    NumericVector RootExudationIN = Rcpp::as<Rcpp::NumericVector>(cb["RootExudation"]);
    for(int i =0;i<numCohorts;i++) {
      LabileCarbonBalance(iday,i) = LabileCarbonBalanceIN[i];
      MaintenanceRespiration(iday,i) = MaintenanceRespirationIN[i];
      GrowthCosts(iday,i) = GrowthCostsIN[i];
      GrossPhotosynthesis(iday,i) = GrossPhotosynthesisIN[i];
      PlantSugarLeaf(iday,i) = PlantSugarLeafIN[i];
      PlantStarchLeaf(iday,i) = PlantStarchLeafIN[i];
      PlantSugarSapwood(iday,i) = PlantSugarSapwoodIN[i];
      PlantStarchSapwood(iday,i) = PlantStarchSapwoodIN[i];
      PlantSugarTransport(iday,i) = PlantSugarTransportIN[i];
      RootExudation(iday,i) = RootExudationIN[i];
    }
  }
  
  if(control["plantStructureResults"]) {
    List plantStructure = Rcpp::as<Rcpp::List>(l["PlantStructure"]);
    NumericMatrix LeafBiomass = Rcpp::as<Rcpp::NumericMatrix>(plantStructure["LeafBiomass"]);
    NumericMatrix SapwoodBiomass = Rcpp::as<Rcpp::NumericMatrix>(plantStructure["SapwoodBiomass"]);
    NumericMatrix FineRootBiomass = Rcpp::as<Rcpp::NumericMatrix>(plantStructure["FineRootBiomass"]);
    NumericMatrix LeafArea = Rcpp::as<Rcpp::NumericMatrix>(plantStructure["LeafArea"]);
    NumericMatrix SapwoodArea = Rcpp::as<Rcpp::NumericMatrix>(plantStructure["SapwoodArea"]);
    NumericMatrix FineRootArea = Rcpp::as<Rcpp::NumericMatrix>(plantStructure["FineRootArea"]);
    NumericMatrix HuberValue = Rcpp::as<Rcpp::NumericMatrix>(plantStructure["HuberValue"]);
    NumericMatrix RootAreaLeafArea = Rcpp::as<Rcpp::NumericMatrix>(plantStructure["RootAreaLeafArea"]);
    NumericMatrix DBH = Rcpp::as<Rcpp::NumericMatrix>(plantStructure["DBH"]);
    NumericMatrix Height = Rcpp::as<Rcpp::NumericMatrix>(plantStructure["Height"]);
    
    NumericVector SapwoodBiomassIN = Rcpp::as<Rcpp::NumericVector>(ps["SapwoodBiomass"]);
    NumericVector LeafBiomassIN = Rcpp::as<Rcpp::NumericVector>(ps["LeafBiomass"]);
    NumericVector FineRootBiomassIN = Rcpp::as<Rcpp::NumericVector>(ps["FineRootBiomass"]);
    NumericVector SapwoodAreaIN = Rcpp::as<Rcpp::NumericVector>(ps["SapwoodArea"]);
    NumericVector LeafAreaIN = Rcpp::as<Rcpp::NumericVector>(ps["LeafArea"]);
    NumericVector FineRootAreaIN = Rcpp::as<Rcpp::NumericVector>(ps["FineRootArea"]);
    NumericVector HuberValueIN = Rcpp::as<Rcpp::NumericVector>(ps["HuberValue"]);
    NumericVector RootAreaLeafAreaIN = Rcpp::as<Rcpp::NumericVector>(ps["RootAreaLeafArea"]);
    NumericVector DBHIN = Rcpp::as<Rcpp::NumericVector>(ps["DBH"]);
    NumericVector HeightIN = Rcpp::as<Rcpp::NumericVector>(ps["Height"]);
    
    for(int i =0;i<numCohorts;i++) {
      SapwoodBiomass(iday,i) = SapwoodBiomassIN[i];
      LeafBiomass(iday,i) = LeafBiomassIN[i];
      FineRootBiomass(iday,i) = FineRootBiomassIN[i];
      SapwoodArea(iday,i) = SapwoodAreaIN[i];
      LeafArea(iday,i) = LeafAreaIN[i];
      FineRootArea(iday,i) = FineRootAreaIN[i];
      HuberValue(iday,i) = HuberValueIN[i];
      RootAreaLeafArea(iday,i) = RootAreaLeafAreaIN[i];
      DBH(iday,i) = DBHIN[i];
      Height(iday,i) = HeightIN[i];
    }
  }
  
  List plantBiomassBalance = Rcpp::as<Rcpp::List>(l["PlantBiomassBalance"]);
  NumericMatrix StructuralBiomassBalance = Rcpp::as<Rcpp::NumericMatrix>(plantBiomassBalance["StructuralBiomassBalance"]);
  NumericMatrix LabileBiomassBalance = Rcpp::as<Rcpp::NumericMatrix>(plantBiomassBalance["LabileBiomassBalance"]);
  NumericMatrix PlantBiomassBalance = Rcpp::as<Rcpp::NumericMatrix>(plantBiomassBalance["PlantBiomassBalance"]);
  NumericMatrix MortalityBiomassLoss = Rcpp::as<Rcpp::NumericMatrix>(plantBiomassBalance["MortalityBiomassLoss"]);
  NumericMatrix CohortBiomassBalance = Rcpp::as<Rcpp::NumericMatrix>(plantBiomassBalance["CohortBiomassBalance"]);
  
  NumericVector StructuralBiomassBalanceIN = Rcpp::as<Rcpp::NumericVector>(bb["StructuralBiomassBalance"]);
  NumericVector LabileBiomassBalanceIN = Rcpp::as<Rcpp::NumericVector>(bb["LabileBiomassBalance"]);
  NumericVector PlantBiomassBalanceIN = Rcpp::as<Rcpp::NumericVector>(bb["PlantBiomassBalance"]);
  NumericVector MortalityBiomassLossIN = Rcpp::as<Rcpp::NumericVector>(bb["MortalityBiomassLoss"]);
  NumericVector CohortBiomassBalanceIN = Rcpp::as<Rcpp::NumericVector>(bb["CohortBiomassBalance"]);

  for(int i =0;i<numCohorts;i++) {
    StructuralBiomassBalance(iday,i) = StructuralBiomassBalanceIN[i];
    LabileBiomassBalance(iday,i) = LabileBiomassBalanceIN[i];
    PlantBiomassBalance(iday,i) = PlantBiomassBalanceIN[i];
    MortalityBiomassLoss(iday,i) = MortalityBiomassLossIN[i];
    CohortBiomassBalance(iday,i) = CohortBiomassBalanceIN[i];
  }
  
  if(control["growthMortalityResults"]) {
    List growthMortality = Rcpp::as<Rcpp::List>(l["GrowthMortality"]);
    NumericMatrix LAgrowth = Rcpp::as<Rcpp::NumericMatrix>(growthMortality["LAgrowth"]);
    NumericMatrix SAgrowth = Rcpp::as<Rcpp::NumericMatrix>(growthMortality["SAgrowth"]);
    NumericMatrix FRAgrowth = Rcpp::as<Rcpp::NumericMatrix>(growthMortality["FRAgrowth"]);
    NumericMatrix starvationRate = Rcpp::as<Rcpp::NumericMatrix>(growthMortality["StarvationRate"]);
    NumericMatrix dessicationRate = Rcpp::as<Rcpp::NumericMatrix>(growthMortality["DessicationRate"]);
    NumericMatrix mortalityRate = Rcpp::as<Rcpp::NumericMatrix>(growthMortality["MortalityRate"]);
    NumericVector LAgrowthIN = Rcpp::as<Rcpp::NumericVector>(gm["LAgrowth"]);
    NumericVector SAgrowthIN = Rcpp::as<Rcpp::NumericVector>(gm["SAgrowth"]);
    NumericVector FRAgrowthIN = Rcpp::as<Rcpp::NumericVector>(gm["FRAgrowth"]);
    NumericVector starvationRateIN = Rcpp::as<Rcpp::NumericVector>(gm["StarvationRate"]);
    NumericVector dessicationRateIN = Rcpp::as<Rcpp::NumericVector>(gm["DessicationRate"]);
    NumericVector mortalityRateIN = Rcpp::as<Rcpp::NumericVector>(gm["MortalityRate"]);
    
    for(int i =0;i<numCohorts;i++) {
      LAgrowth(iday,i) = LAgrowthIN[i];
      SAgrowth(iday,i) = SAgrowthIN[i];
      FRAgrowth(iday,i) = FRAgrowthIN[i];
      starvationRate(iday,i) = starvationRateIN[i];
      dessicationRate(iday,i) = dessicationRateIN[i];
      mortalityRate(iday,i) = mortalityRateIN[i];
    }
  }
  
  NumericMatrix StandBiomassBalance = Rcpp::as<Rcpp::NumericMatrix>(l["BiomassBalance"]);
  StandBiomassBalance(iday,_) = standLevelBiomassBalance(bb);
  
  NumericMatrix StandCarbonBalance = Rcpp::as<Rcpp::NumericMatrix>(l["CarbonBalance"]);
  StandCarbonBalance(iday,_) = Rcpp::as<Rcpp::NumericVector>(sDay["CarbonBalance"]);
  
}

//' Forest growth
//' 
//' Function \code{growth} is a process-based model that performs energy, water and carbon balances; 
//' and determines changes in water/carbon pools, functional variables (leaf area, sapwood area, root area) 
//' and structural ones (tree diameter, tree height, shrub cover) for woody plant cohorts in a given forest stand 
//' during a period specified in the input climatic data. 
//' 
//' @param x An object of class \code{\link{growthInput}}.
//' @param meteo A data frame with daily meteorological data series (see \code{\link{spwb}}).
//' @param latitude Latitude (in degrees).
//' @param elevation,slope,aspect Elevation above sea level (in m), slope (in degrees) and aspect (in degrees from North). 
//' @param CO2ByYear A named numeric vector with years as names and atmospheric CO2 concentration (in ppm) as values. Used to specify annual changes in CO2 concentration along the simulation (as an alternative to specifying daily values in \code{meteo}).
//' @param waterTableDepth Water table depth (in mm). When not missing, capillarity rise will be allowed if lower than total soil depth.
//' 
//' @details
//' Detailed model description is available in the medfate book. 
//' Simulations using the 'Sperry' or 'Sureau' transpiration modes are computationally much more expensive 
//' than those using the 'Granier' transpiration mode. 
//' 
//' @return
//' A list of class 'growth' with the following elements:
//' \itemize{
//'   \item{\code{"latitude"}: Latitude (in degrees) given as input.} 
//'   \item{\code{"topography"}: Vector with elevation, slope and aspect given as input.} 
//'   \item{\code{"weather"}: A copy of the input weather data frame.}
//'   \item{\code{"growthInput"}: A copy of the object \code{x} of class \code{\link{growthInput}} given as input.}
//'   \item{\code{"growthOutput"}: An copy of the final state of the object \code{x} of class \code{\link{growthInput}}.}
//'   \item{\code{"WaterBalance"}: A data frame where different water balance variables (see \code{\link{spwb}}).}
//'   \item{\code{"EnergyBalance"}: A data frame with the daily values of energy balance components for the soil and the canopy (only for \code{transpirationMode = "Sperry"} or \code{transpirationMode = "Sureau"}; see \code{\link{spwb}}).}
//'   \item{\code{"CarbonBalance"}: A data frame where different stand-level carbon balance components (gross primary production, maintenance respiration, synthesis respiration and net primary production), all in g C · m-2.}
//'   \item{\code{"BiomassBalance"}: A data frame with the daily values of stand biomass balance components (in g dry · m-2.}
//'   \item{\code{"Temperature"}: A data frame with the daily values of minimum/mean/maximum temperatures for the atmosphere (input), canopy and soil (only for \code{transpirationMode = "Sperry"} or \code{transpirationMode = "Sureau"}; see \code{\link{spwb}}).}
//'   \item{\code{"Soil"}: A data frame where different soil variables  (see \code{\link{spwb}}).}
//'   \item{\code{"Stand"}: A data frame where different stand-level variables (see \code{\link{spwb}}).}
//'   \item{\code{"Plants"}: A list of daily results for plant cohorts (see \code{\link{spwb}}).}
//'   \item{\code{"SunlitLeaves"} and \code{"ShadeLeaves"}: A list with daily results for sunlit and shade leaves (only for \code{transpirationMode = "Sperry"} or \code{transpirationMode = "Sureau"}; see \code{\link{spwb}}).}
//'   \item{\code{"LabileCarbonBalance"}: A list of daily labile carbon balance results for plant cohorts, with elements:}
//'   \itemize{
//'     \item{\code{"GrossPhotosynthesis"}: Daily gross photosynthesis per dry weight of living biomass (g gluc · g dry-1).}
//'     \item{\code{"MaintentanceRespiration"}: Daily maintenance respiration per dry weight of living biomass (g gluc · g dry-1).}
//'     \item{\code{"GrowthCosts"}: Daily growth costs per dry weight of living biomass (g gluc · g dry-1).}
//'     \item{\code{"RootExudation"}: Root exudation per dry weight of living biomass (g gluc · g dry-1).}    
//'     \item{\code{"LabileCarbonBalance"}: Daily labile carbon balance (photosynthesis - maintenance respiration - growth costs - root exudation) per dry weight of living biomass (g gluc · g dry-1).}
//'     \item{\code{"SugarLeaf"}: Sugar concentration (mol·l-1) in leaves.}
//'     \item{\code{"StarchLeaf"}: Starch concentration (mol·l-1) in leaves.}
//'     \item{\code{"SugarSapwood"}: Sugar concentration (mol·l-1) in sapwood.}
//'     \item{\code{"StarchSapwood"}: Starch concentration (mol·l-1) in sapwood.}
//'     \item{\code{"SugarTransport"}:  Average instantaneous rate of carbon transferred between leaves and stem compartments via floem (mol gluc·s-1).}
//'   }
//'   \item{\code{"PlantBiomassBalance"}: A list of daily plant biomass balance results for plant cohorts, with elements:}
//'   \itemize{
//'     \item{\code{"StructuralBiomassBalance"}: Daily structural biomass balance (g dry · m-2).}
//'     \item{\code{"LabileBiomassBalance"}: Daily labile biomass balance (g dry · m-2).}
//'     \item{\code{"PlantBiomassBalance"}: Daily plant biomass balance, i.e. labile change + structural change (g dry · m-2).}
//'     \item{\code{"MortalityBiomassLoss"}: Biomass loss due to mortality (g dry · m-2).}    
//'     \item{\code{"CohortBiomassBalance"}: Daily cohort biomass balance (including mortality) (g dry · m-2).}
//'   }
//'   \item{\code{"PlantStructure"}: A list of daily area and biomass values for compartments of plant cohorts, with elements:}
//'   \itemize{
//'     \item{\code{"LeafBiomass"}: Daily amount of leaf structural biomass (in g dry) for an average individual of each plant cohort.}
//'     \item{\code{"SapwoodBiomass"}: Daily amount of sapwood structural biomass (in g dry) for an average individual of each plant cohort.}
//'     \item{\code{"FineRootBiomass"}: Daily amount of fine root biomass (in g dry) for an average individual of each plant cohort.}
//'     \item{\code{"LeafArea"}: Daily amount of leaf area (in m2) for an average individual of each plant cohort.}
//'     \item{\code{"SapwoodArea"}: Daily amount of sapwood area (in cm2) for an average individual of each plant cohort.}
//'     \item{\code{"FineRootArea"}: Daily amount of fine root area (in m2) for an average individual of each plant cohort.}
//'     \item{\code{"HuberValue"}: The ratio of sapwood area to (target) leaf area (in cm2/m2).}
//'     \item{\code{"RootAreaLeafArea"}: The ratio of fine root area to (target) leaf area (in m2/m2).}
//'     \item{\code{"DBH"}: Diameter at breast height (in cm) for an average individual of each plant cohort.}
//'     \item{\code{"Height"}: Height (in cm) for an average individual of each plant cohort.}
//'   }
//'   \item{\code{"GrowthMortality"}: A list of daily growth and mortality rates for plant cohorts, with elements:}
//'   \itemize{
//'     \item{\code{"LAgrowth"}: Leaf area growth (in m2·day-1) for an average individual of each plant cohort.}
//'     \item{\code{"SAgrowth"}: Sapwood area growth rate (in cm2·day-1) for an average individual of each plant cohort.}
//'     \item{\code{"FRAgrowth"}: Fine root area growth (in m2·day-1) for an average individual of each plant cohort.}
//'     \item{\code{"StarvationRate"}: Daily mortality rate from starvation (ind/d-1).}
//'     \item{\code{"DessicationRate"}: Daily mortality rate from dessication (ind/d-1).}
//'     \item{\code{"MortalityRate"}: Daily mortality rate (any cause) (ind/d-1).}
//'   }
//'   \item{\code{"subdaily"}: A list of objects of class \code{\link{growth_day}}, one per day simulated (only if required in \code{control} parameters, see \code{\link{defaultControl}}).}
//' }
//' 
//' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
//' 
//' @seealso \code{\link{growthInput}}, \code{\link{growth_day}}, \code{\link{plot.growth}}
//' 
//' @references
//' De Cáceres M, Molowny-Horas R, Cabon A, Martínez-Vilalta J, Mencuccini M, García-Valdés R, Nadal-Sala D, Sabaté S, 
//' Martin-StPaul N, Morin X, D'Adamo F, Batllori E, Améztegui A (2023) MEDFATE 2.9.3: A trait-enabled model to simulate 
//' Mediterranean forest function and dynamics at regional scales. 
//' Geoscientific Model Development 16: 3165-3201 (https://doi.org/10.5194/gmd-16-3165-2023).
//' 
//' @examples
//' \donttest{
//' #Load example daily meteorological data
//' data(examplemeteo)
//' 
//' #Load example plot plant data
//' data(exampleforest)
//' 
//' #Default species parameterization
//' data(SpParamsMED)
//'   
//' #Initialize control parameters
//' control <- defaultControl("Granier")
//'   
//' #Initialize soil with default soil params (4 layers)
//' examplesoil <- defaultSoilParams(4)
//' 
//' #Initialize model input
//' x1 <- growthInput(exampleforest, examplesoil, SpParamsMED, control)
//' 
//' #Call simulation function
//' G1 <- growth(x1, examplemeteo, latitude = 41.82592, elevation = 100)
//'  
//' #Switch to 'Sperry' transpiration mode
//' control <- defaultControl("Sperry")
//' 
//' #Initialize model input
//' x2 <- growthInput(exampleforest,examplesoil, SpParamsMED, control)
//' 
//' #Call simulation function
//' G2 <-growth(x2, examplemeteo, latitude = 41.82592, elevation = 100)
//' 
//' #Switch to 'Sureau' transpiration mode
//' control <- defaultControl("Sureau")
//' 
//' #Initialize model input
//' x3 <- growthInput(exampleforest,examplesoil, SpParamsMED, control)
//' 
//' #Call simulation function
//' G3 <-growth(x3, examplemeteo, latitude = 41.82592, elevation = 100)
//' }
//'       
// [[Rcpp::export("growth")]]
List growth(List x, DataFrame meteo, double latitude, 
            double elevation, double slope = NA_REAL, double aspect = NA_REAL,
            NumericVector CO2ByYear = NumericVector(0), double waterTableDepth = NA_REAL) {

  //Clone input
  x = clone(x);
  
  //Control params 
  List control =x["control"];  
  String transpirationMode = control["transpirationMode"];
  String soilFunctions = control["soilFunctions"];
  bool verbose = control["verbose"];
  bool leafPhenology = control["leafPhenology"];
  bool unlimitedSoilWater = control["unlimitedSoilWater"];
  NumericVector defaultRainfallIntensityPerMonth = control["defaultRainfallIntensityPerMonth"];
  checkgrowthInput(x, transpirationMode, soilFunctions);


  //Cohort info
  DataFrame cohorts = Rcpp::as<Rcpp::DataFrame>(x["cohorts"]);
  IntegerVector SP = Rcpp::as<Rcpp::IntegerVector>(cohorts["SP"]);
  IntegerVector SPunique = uniqueSpp(SP);

  if(NumericVector::is_na(latitude)) stop("Value for 'latitude' should not be missing.");
  double latrad = latitude * (M_PI/180.0);
  if(NumericVector::is_na(aspect)) aspect = 0.0;
  if(NumericVector::is_na(slope)) slope = 0.0;
  double asprad = aspect * (M_PI/180.0);
  double slorad = slope * (M_PI/180.0);
  
  //Meteorological input    
  NumericVector MinTemperature, MaxTemperature;
  NumericVector MinRelativeHumidity, MaxRelativeHumidity;
  NumericVector Radiation;
  if(!meteo.containsElementNamed("Precipitation")) stop("Please include variable 'Precipitation' in weather input.");
  NumericVector Precipitation = meteo["Precipitation"];
  int numDays = Precipitation.size();
  if(NumericVector::is_na(elevation)) stop("Value for 'elevation' should not be missing.");
  if(!meteo.containsElementNamed("MinTemperature")) stop("Please include variable 'MinTemperature' in weather input.");
  MinTemperature = meteo["MinTemperature"];
  if(!meteo.containsElementNamed("MaxTemperature")) stop("Please include variable 'MaxTemperature' in weather input.");
  MaxTemperature = meteo["MaxTemperature"];
  if(!meteo.containsElementNamed("MinRelativeHumidity")) stop("Please include variable 'MinRelativeHumidity' in weather input.");
  MinRelativeHumidity = meteo["MinRelativeHumidity"];
  if(!meteo.containsElementNamed("MaxRelativeHumidity")) stop("Please include variable 'MaxRelativeHumidity' in weather input.");
  MaxRelativeHumidity = meteo["MaxRelativeHumidity"];
  if(!meteo.containsElementNamed("Radiation")) stop("Please include variable 'Radiation' in weather input.");
  Radiation = meteo["Radiation"];
  
  if(any(is_na(Precipitation))) stop("Missing values in 'Precipitation'");
  if(any(is_na(MinTemperature))) stop("Missing values in 'MinTemperature'");
  if(any(is_na(MaxTemperature))) stop("Missing values in 'MaxTemperature'");
  if(any(is_na(MinRelativeHumidity))) warning("Missing values in 'MinRelativeHumidity' were estimated from temperature range");
  if(any(is_na(MaxRelativeHumidity))) warning("Missing values in 'MaxRelativeHumidity' were assumed to be 100");
  if(any(is_na(Radiation))) warning("Missing values in 'Radiation' were estimated");
  
  NumericVector WindSpeed(numDays, NA_REAL);
  if(meteo.containsElementNamed("WindSpeed")) WindSpeed = meteo["WindSpeed"];
  
  NumericVector FireProbability(numDays, 0.0);
  if(meteo.containsElementNamed("FireProbability")) FireProbability = meteo["FireProbability"];
  
  NumericVector PET = NumericVector(numDays, NA_REAL);
  
  NumericVector CO2(Precipitation.length(), NA_REAL);
  if(meteo.containsElementNamed("CO2")) {
    CO2 = meteo["CO2"];
    if(verbose) {
      Rcout<<"CO2 taken from input column 'CO2'\n";
    }
    if(any(is_na(CO2))) stop("Missing values in 'CO2'");
  }
  NumericVector Patm(Precipitation.length(), NA_REAL);
  if(meteo.containsElementNamed("Patm")) {
    Patm = meteo["Patm"];
    if(verbose) {
      Rcout<<"Patm taken from input column 'Patm'\n";
    }
  }
  NumericVector RainfallIntensity(Precipitation.length(), NA_REAL);
  if(meteo.containsElementNamed("RainfallIntensity")) {
    RainfallIntensity = meteo["RainfallIntensity"];
    if(verbose) {
      Rcout<<"Rainfall intensity taken from input column 'RainfallIntensity'\n";
    }
  }
  
  IntegerVector DOY, JulianDay;
  NumericVector Photoperiod;
  bool doy_input = false, photoperiod_input = false, julianday_input = false;
  if(meteo.containsElementNamed("DOY")) {
    DOY = meteo["DOY"];
    doy_input = true;
    if(verbose) {
      Rcout<<"DOY taken from input column 'DOY'\n";
    }
  }
  if(meteo.containsElementNamed("Photoperiod")) {
    Photoperiod = meteo["Photoperiod"];
    photoperiod_input = true;
    if(verbose) {
      Rcout<<"Photoperiod taken from input column 'Photoperiod'\n";
    }
  }
  if(meteo.containsElementNamed("JulianDay")) {
    JulianDay = meteo["JulianDay"];
    julianday_input = true;
    if(verbose) {
      Rcout<<"Julian day taken from input column 'JulianDay'\n";
    }
  }
  
  // Dates
  CharacterVector dateStrings = getWeatherDates(meteo);
  if(!doy_input) DOY = date2doy(dateStrings);
  if(!photoperiod_input) Photoperiod = date2photoperiod(dateStrings, latrad);
  
  //Soil params 
  DataFrame soil = Rcpp::as<Rcpp::DataFrame>(x["soil"]);
  
  //Output list
  List outputList = defineGrowthDailyOutput(latitude, elevation, slope, aspect,
                                            dateStrings, x);
  outputList["weather"] = clone(meteo);

  //Count years (times structural variables will be updated)
  int numYears = 0;
  for(int i=0;i<numDays;i++) {
    if(((DOY[i]==1) && (i>0)) || ((i==(numDays-1)) && (DOY[i]>=365))) numYears = numYears + 1;
  }

  NumericVector initialSoilContent = water(soil, soilFunctions);
  NumericVector initialPlantContent = plantWaterContent(x);
  double initialSnowContent = x["snowpack"];
  
  
  DataFrame ccIni_m2 = carbonCompartments(x, "g_m2");
  double cohortBiomassBalanceSum = 0.0;
  double initialCohortBiomass = sum(Rcpp::as<Rcpp::NumericVector>(ccIni_m2["TotalBiomass"]));
  
  if(verbose) {
    Rcout<<"Initial plant cohort biomass (g/m2): "<<initialCohortBiomass<<"\n";
    Rcout<<"Initial plant water content (mm): "<< sum(initialPlantContent)<<"\n";
    Rcout<<"Initial soil water content (mm): "<< sum(initialSoilContent)<<"\n";
    Rcout<<"Initial snowpack content (mm): "<< initialSnowContent<<"\n";
  }
  
  //Instance communication structures
  List internalCommunication = instanceCommunicationStructures(x, "growth");

  bool error_occurence = false;
  if(verbose) Rcout << "Performing daily simulations\n";
  List s;
  std::string yearString;
  for(int i=0;i<numDays;i++) {
    std::string c = as<std::string>(dateStrings[i]);
    yearString = c.substr(0, 4);
    if(verbose) {
      if(DOY[i]==1 || i==0) {
        Rcout<<"\n Year "<< yearString<< ":";
      } 
      else if(i%30 == 0) Rcout<<".";//<<i;
    } 
    
    double wind = WindSpeed[i];
    if(NumericVector::is_na(wind)) wind = control["defaultWindSpeed"]; //Default 1 m/s -> 10% of fall every day
    if(wind<0.1) wind = 0.1; //Minimum windspeed abovecanopy
    
    double Catm = CO2[i];
    //If missing, use
    if(NumericVector::is_na(Catm)) {
      if(CO2ByYear.attr("names") != R_NilValue) Catm = CO2ByYear[yearString];
    }
    //If still missing, use default control value
    if(NumericVector::is_na(Catm)) {
      Catm = control["defaultCO2"];
    }
    
    double Rint = RainfallIntensity[i];
    if(NumericVector::is_na(Rint)) {
      int month = std::atoi(c.substr(5,2).c_str());
      Rint = rainfallIntensity(month, Precipitation[i], defaultRainfallIntensityPerMonth);
    }
    
    if(unlimitedSoilWater) {
      NumericVector W = soil["W"];
      for(int h=0;h<W.size();h++) W[h] = 1.0;
    }
    

    //Julian day from either input column or date
    int J = NA_INTEGER;
    if(julianday_input) J = JulianDay[i];
    if(IntegerVector::is_na(J)){
      std::string c = as<std::string>(dateStrings[i]);
      J = meteoland::radiation_julianDay(std::atoi(c.substr(0, 4).c_str()),std::atoi(c.substr(5,2).c_str()),std::atoi(c.substr(8,2).c_str())); 
    }
    double delta = meteoland::radiation_solarDeclination(J);
    double solarConstant = meteoland::radiation_solarConstant(J);
    
    double tmin = MinTemperature[i];
    double tmax = MaxTemperature[i];
    if(tmin > tmax) {
      warning("tmin > tmax. Swapping values.");
      double swap = tmin;
      tmin = tmax;
      tmax = swap;
    }
    double prec = Precipitation[i];
    double tday = meteoland::utils_averageDaylightTemperature(tmin, tmax);
    double rhmin = MinRelativeHumidity[i];
    double rhmax = MaxRelativeHumidity[i];
    double rad = Radiation[i];
    if(NumericVector::is_na(rhmax)) {
      rhmax = 100.0;
    }
    if(NumericVector::is_na(rhmin)) {
      double vp_tmin = meteoland::utils_saturationVP(tmin);
      double vp_tmax = meteoland::utils_saturationVP(tmax);
      rhmin = std::min(rhmax, 100.0*(vp_tmin/vp_tmax));
    }
    if(rhmin > rhmax) {
      warning("rhmin > rhmax. Swapping values.");
      double swap = rhmin;
      rhmin = rhmax;
      rhmax = swap;
    }
    if(NumericVector::is_na(rad)) {
      double vpa = meteoland::utils_averageDailyVP(tmin, tmax, rhmin, rhmax);
      rad = meteoland::radiation_solarRadiation(solarConstant, latrad, elevation,
                                                slorad, asprad, delta, tmax -tmin, tmax-tmin,
                                                vpa, prec);
    }
    PET[i] = meteoland::penman(latrad, elevation, slorad, asprad, J, 
                               tmin, tmax, rhmin, rhmax, rad, wind);
    
    //1. Phenology (only leaf fall)
    if(leafPhenology) {
      updatePhenology(x, DOY[i], Photoperiod[i], tday);
      updateLeaves(x, wind, true);
    }
    
    //2. Water balance and photosynthesis
    if(transpirationMode=="Granier") {
      NumericVector meteovec_inner = NumericVector::create(
        Named("tday") = tday, Named("tmax") = tmax, Named("tmin") = tmin,
        Named("prec") = prec, 
        Named("rhmin") = rhmin, Named("rhmax") = rhmax,
        Named("rad") = rad, 
        Named("wind") = wind, 
        Named("pet") = PET[i],
        Named("Catm") = Catm);
      meteovec_inner.push_back(Patm[i], "Patm");
      meteovec_inner.push_back(Rint, "rint");
      meteovec_inner.push_back(FireProbability[i], "pfire"); 
      try{
        growthDay_private(internalCommunication, x, meteovec_inner,  
                           latitude, elevation, slope, aspect,
                           solarConstant, delta, 
                           0.0, R_NilValue, waterTableDepth,
                           false); //No Runon in simulations for a single cell
        fillGrowthDailyOutput(outputList, x, internalCommunication["basicGROWTHOutput"], i);
      } catch(std::exception& ex) {
        Rcerr<< "c++ error: "<< ex.what() <<"\n";
        error_occurence = true;
      }
    } else {
      double tmaxPrev = tmax;
      double tminPrev = tmin;
      double tminNext = tmin;
      if(i>0) {
        tmaxPrev = MaxTemperature[i-1];
        tminPrev = MinTemperature[i-1];
      }
      if(i<(numDays-1)) tminNext = MinTemperature[i+1]; 
      NumericVector meteovec = NumericVector::create(
        Named("tday") = tday,
        Named("tmin") = tmin, 
        Named("tmax") = tmax,
        Named("tminPrev") = tminPrev, 
        Named("tmaxPrev") = tmaxPrev, 
        Named("tminNext") = tminNext, 
        Named("prec") = Precipitation[i],
        Named("rhmin") = rhmin, 
        Named("rhmax") = rhmax, 
        Named("rad") = rad, 
        Named("wind") = wind, 
        Named("Catm") = Catm,
        Named("Patm") = Patm[i],
        Named("pet") = PET[i],
        Named("rint") = Rint);
      meteovec.push_back(FireProbability[i], "pfire"); 
      try{
        growthDay_private(internalCommunication, x, meteovec, 
                           latitude, elevation, slope, aspect,
                           solarConstant, delta, 
                           0.0, R_NilValue, waterTableDepth,
                           verbose);
        fillGrowthDailyOutput(outputList, x, internalCommunication["advancedGROWTHOutput"], i);
      } catch(std::exception& ex) {
        Rcerr<< "c++ error: "<< ex.what() <<"\n";
        error_occurence = true;
      }
    }    
    
    //Add cohort biomass sum
    List plantBiomassBalance = Rcpp::as<Rcpp::List>(outputList["PlantBiomassBalance"]);
    NumericMatrix CohortBiomassBalance = Rcpp::as<Rcpp::NumericMatrix>(plantBiomassBalance["CohortBiomassBalance"]);
    cohortBiomassBalanceSum += sum(CohortBiomassBalance(i,_));
  }
  if(verbose) Rcout << "\n\n";
  
  // Check biomass balance
  DataFrame ccFin_m2 = carbonCompartments(x, "g_m2");
  double finalCohortBiomass = sum(Rcpp::as<Rcpp::NumericVector>(ccFin_m2["TotalBiomass"]));
  if(verbose) {
    NumericMatrix StandBiomassBalance = Rcpp::as<Rcpp::NumericMatrix>(outputList["BiomassBalance"]);
    
    Rcout<<"Final plant biomass (g/m2): "<<finalCohortBiomass<<"\n";
    Rcout<<"Change in plant biomass (g/m2): " << finalCohortBiomass - initialCohortBiomass <<"\n";
    Rcout<<"Plant biomass balance result (g/m2): " <<  cohortBiomassBalanceSum<<"\n";
    Rcout<<"Plant biomass balance components:\n";
    
    Rcout<<"  Structural balance (g/m2) "  <<round(sum(StandBiomassBalance(_,0)))<<" Labile balance (g/m2) "  <<round(sum(StandBiomassBalance(_,1))) <<"\n";
    Rcout<<"  Plant individual balance (g/m2) "  <<round(sum(StandBiomassBalance(_,2)))<<" Mortality loss (g/m2) "  <<round(sum(StandBiomassBalance(_,3))) <<"\n";
    
    printWaterBalanceResult(outputList, x,
                            initialPlantContent, initialSoilContent, initialSnowContent,
                            transpirationMode);
    
    if(error_occurence) {
      Rcout<< " ERROR: Calculations stopped because of numerical error: Revise parameters\n";
    }
  }
  
  return(outputList);
}
