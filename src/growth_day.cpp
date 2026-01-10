// [[Rcpp::interfaces(r,cpp)]]
#define STRICT_R_HEADERS
#include <Rcpp.h>
#include <numeric>
#include "biophysicsutils.h"
#include "carbon.h"
#include "communication_structures.h"
#include "decomposition.h"
#include "forestutils.h"
#include "fireseverity.h"
#include "firebehaviour.h"
#include "hydraulics.h"
#include "hydrology.h"
#include "lightextinction_basic.h"
#include "lightextinction_advanced.h"
#include "modelInput.h"
#include "phenology.h"
#include "root.h"
#include "soil.h"
#include "spwb.h"
#include "spwb_day.h"
#include "tissuemoisture.h"
#include "woodformation.h"
#include <meteoland.h>
using namespace Rcpp;

//' Mortality
//' 
//' A simple sigmoid function to determine a daily mortality likelihood according 
//' to the value of a stress variable.
//'
//' @param stressValue Current value of the stress variable (0 to 1, 
//'                    with higher values indicate stronger stress).
//' @param stressThreshold Threshold to indicate 50% annual mortality probability.
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

void updateStructuralVariables(List x, NumericVector deltaSAgrowth) {
  
  //Control params
  List control = x["control"];  
  bool shrubDynamics = control["shrubDynamics"];
  bool herbDynamics = control["herbDynamics"];
  
  //Cohort info
  DataFrame cohorts = Rcpp::as<Rcpp::DataFrame>(x["cohorts"]);
  CharacterVector ctype = cohortType(cohorts.attr("row.names"));
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
      if(ctype[j] == "tree") treeLAI +=LAI_live[j];
    }
    for(int j=0;j<numCohorts; j++) {
      if((ctype[j]=="shrub") && N[j]>0.0) {
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

  ///////////////////////////////////////////////////////////////////////////////////
  ///// A. WATER-ENERGY BALANCE (this creates communication structures as well) /////
  ///////////////////////////////////////////////////////////////////////////////////
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

  ///////////////////////////////////////////
  ///// B. FIRE OCCURRENCE AND BEHAVIOR /////
  ///////////////////////////////////////////
  bool fireOccurrence = false;
  double pfire = meteovec["pfire"];
  NumericVector fireBehavior = communicationFireHazard();
  if(R::runif(0.0,1.0) < pfire) {
    fccsHazard(fireBehavior, x, meteovec, spwbOutput, slope);
    fireOccurrence = true;
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
  double tmin = meteovec["tmin"];
  double tmax = meteovec["tmax"];
  double rhmin = meteovec["rhmin"];
  double rhmax = meteovec["rhmax"];
  double Patm = meteovec["Patm"];
  
  //Atmospheric pressure (if missing)
  if(NumericVector::is_na(Patm)) Patm = meteoland::utils_atmosphericPressure(elevation);
  
  //Cohort info
  DataFrame cohorts = Rcpp::as<Rcpp::DataFrame>(x["cohorts"]);
  CharacterVector speciesNames = cohorts["Name"];
  CharacterVector ctype = cohortType(cohorts.attr("row.names"));
  IntegerVector SP = Rcpp::as<Rcpp::IntegerVector>(cohorts["SP"]);
  int numCohorts = SP.size();
  
  //Soil
  DataFrame soil = Rcpp::as<Rcpp::DataFrame>(x["soil"]);
  int nlayers = soil.nrow();
  NumericVector psiSoil = psi(soil, soilFunctions);
  NumericVector widths = soil["widths"];
  NumericVector sand = soil["sand"];
  NumericVector clay = soil["clay"];
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
  NumericVector N_resprouting_stumps = internalMortality["N_resprouting_stumps"];
  NumericVector Cover_dead = internalMortality["Cover_dead"];
  NumericVector Cover_starvation = internalMortality["Cover_starvation"];
  NumericVector Cover_dessication = internalMortality["Cover_dessication"];
  NumericVector Cover_burnt = internalMortality["Cover_burnt"];
  NumericVector Cover_resprouting_stumps = internalMortality["Cover_resprouting_stumps"];
  NumericVector Snag_smallbranches = internalMortality["Snag_smallbranches"];
  NumericVector Snag_largewood = internalMortality["Snag_largewood"];
  
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
  NumericVector r635 = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["r635"]);
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
  NumericVector RespFire = Rcpp::as<Rcpp::NumericVector>(paramsMortalityRegeneration["RespFire"]);
  NumericVector RespDist = Rcpp::as<Rcpp::NumericVector>(paramsMortalityRegeneration["RespDist"]);
  
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
  
  //Litter
  DataFrame internalLitter, internalSnags;
  DataFrame paramsLitterDecomposition;
  NumericVector internalSOC;
  bool exportLitter = false;
  if(x.containsElementNamed("internalLitter") && x.containsElementNamed("internalSnags") && x.containsElementNamed("internalSOC") && x.containsElementNamed("paramsLitterDecomposition")) {
    paramsLitterDecomposition = Rcpp::as<Rcpp::DataFrame>(x["paramsLitterDecomposition"]);
    internalSnags = Rcpp::as<Rcpp::DataFrame>(x["internalSnags"]);
    internalLitter = Rcpp::as<Rcpp::DataFrame>(x["internalLitter"]);
    internalSOC = Rcpp::as<Rcpp::NumericVector>(x["internalSOC"]);
    exportLitter = true;
  }
  
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
  
  // Fire combustion (in C/m2) for stand carbon balance
  double fireCombustion = 0.0;
  
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
  NumericVector AbovegroundWoodBiomass = Rcpp::as<Rcpp::NumericVector>(ccIni["AbovegroundWoodBiomass"]);
  NumericVector BelowgroundWoodBiomass = Rcpp::as<Rcpp::NumericVector>(ccIni["BelowgroundWoodBiomass"]);
  NumericVector SapwoodLivingStructBiomass = Rcpp::as<Rcpp::NumericVector>(ccIni["SapwoodLivingStructuralBiomass"]);
  NumericVector TwigLivingStructuralBiomass = Rcpp::as<Rcpp::NumericVector>(ccIni["TwigLivingStructuralBiomass"]);
  NumericVector TotalLivingBiomass = Rcpp::as<Rcpp::NumericVector>(ccIni["TotalLivingBiomass"]);


  //For survival model based on basal area  
  double treeBasalArea = 0.0;
  for(int j=0;j<numCohorts;j++){
    if(!NumericVector::is_na(DBH[j])) treeBasalArea += N[j]*3.141593*pow(DBH[j]/200,2.0);
  }

  /////////////////////////////////////////////////////////////////////
  ///// C. LABILE CARBON BALANCE, GROWTH AND SENESCENCE BY COHORT /////
  /////////////////////////////////////////////////////////////////////
  for(int j=0;j<numCohorts;j++){
    if(N[j] > 0.0) {
      
      ///// C1. COSTS AND DERIVED VARIABLES /////
      double costPerLA = 1000.0*CCleaf[j]/SLA[j]; // Construction cost in g gluc · m-2 of leaf area
      double twigCostPerLA = 1000.0*CCsapwood[j]*(r635[j] - 1.0)/SLA[j];// Include cost of twigs in g gluc · m-2 of leaf area
      double costPerSA = CCsapwood[j]*sapwoodStructuralBiomass(1.0, H[j], L(j,_),V(j,_),WoodDensity[j]); // Construction cost in g gluc · cm-2 of sapwood area
      NumericVector deltaFRBgrowth(nlayers, 0.0);
      
      double LAexpanded = leafArea(LAI_expanded[j], N[j]);
      double LAlive = leafArea(LAI_live[j], N[j]);
      double LAdead = leafArea(LAI_dead[j], N[j]);

      double minimumStarchForSecondaryGrowth = Starch_max_sapwood[j]*RSSG[j];
      double minimumStarchForPrimaryGrowth = Starch_max_sapwood[j]*0.1;
        
      //Photosynthesis
      double leafAgG = 0.0;
      //Maintenance respiration
      double leafRespDay = 0.0;
      double twigResp = 0.0;
      double sapwoodResp = 0.0;
      double finerootResp = 0.0;
      //Growth costs
      double growthCostLA = 0.0;
      double twigGrowthCostLA = 0.0;
      double growthCostFRB = 0.0;   
      //Synthesis respiration
      double synthesisRespLA = 0.0;
      double twigSynthesisRespLA = 0.0;
      double synthesisRespSA = 0.0;
      double synthesisRespFRB = 0.0;   
      double growthCostSA = 0.0;
      //Initial labile (g gluc)
      double initialLabile =  glucoseMolarMass*((sugarLeaf[j] + starchLeaf[j])*Volume_leaves[j] + (sugarSapwood[j] + starchSapwood[j])*Volume_sapwood[j]);
      
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
        
        ////// C2. PHOTOSYNTHESIS //////
        if(LAexpanded>0.0) {
          standGrossPrimaryProduction += Ag[j]; //Add GPP in gC · m-2
          
          //gross fotosynthesis
          double leafAgC = Ag[j]/(N[j]/10000.0); //Translate g C · m-2 to g C ·ind-1
          leafAgG = leafAgC*(glucoseMolarMass/(carbonMolarMass*6.0)); // from g C·ind-1 to g gluc · ind-1
          
          //Update output values
          GrossPhotosynthesis[j] = leafAgG/TotalLivingBiomass[j]; //Ag in g gluc · gdry-1 
        }
        
        ///// C3. MAINTENANCE RESPIRATION ///////
        //Respiratory biomass (g dw · ind-1)
        double leafSugarMass = sugarLeaf[j]*(Volume_leaves[j]*glucoseMolarMass);
        double sapwoodSugarMass = sugarSapwood[j]*(Volume_sapwood[j]*glucoseMolarMass);
        double B_resp_leaves = LeafStructBiomass[j] + leafSugarMass;
        double B_resp_sapwood = SapwoodLivingStructBiomass[j]  + sapwoodSugarMass;
        double B_resp_twig = TwigLivingStructuralBiomass[j];
        // Rcout<<j<< " maintenance costs of leaf sugars: "<< (leafSugarMass/B_resp_leaves)<<" sapwood sugars: "<< (sapwoodSugarMass/B_resp_sapwood)<<"\n";
        double B_resp_fineroots = fineRootBiomass[j];
        double QR = qResp(tday);
        if(LAexpanded>0.0) {
          leafRespDay = B_resp_leaves*RERleaf[j]*QR*std::min(1.0, pow(PARcohort[j]/100.0, 0.5)); //Reduction under shade
        }
        twigResp = B_resp_twig*RERsapwood[j]*QR;
        sapwoodResp = B_resp_sapwood*RERsapwood[j]*QR;
        finerootResp = B_resp_fineroots*RERfineroot[j]*QR*(LAexpanded/LAlive);
        MaintenanceRespiration[j] = (leafRespDay+twigResp + sapwoodResp+finerootResp)/TotalLivingBiomass[j]; 

        ///// C4. LEAF/TWIG GROWTH /////
        if(leafUnfolding[j]) {
          double deltaLApheno = std::max(leafAreaTarget[j]*(1.0 - StemPLC[j]) - LAexpanded, 0.0);
          double deltaLAsink = std::min(deltaLApheno, (crownBudPercent[j]/100.0)*SA[j]*RGRleafmax[j]*(rleafcell/rcellmax));
          if(!sinkLimitation) deltaLAsink = std::min(deltaLApheno, (crownBudPercent[j]/100.0)*SA[j]*RGRleafmax[j]); //Deactivates temperature and turgor limitation
          double deltaLAavailable = 0.0;
          deltaLAavailable = std::max(0.0, (starchSapwood[j]-minimumStarchForPrimaryGrowth)*(glucoseMolarMass*Volume_sapwood[j])/costPerLA);
          deltaLAgrowth[j] = std::min(deltaLAsink, deltaLAavailable);
          growthCostLA = deltaLAgrowth[j]*costPerLA;
          synthesisRespLA = growthCostLA*(CCleaf[j] - 1.0)/CCleaf[j];
          twigGrowthCostLA = deltaLAgrowth[j]*twigCostPerLA;
          twigSynthesisRespLA = twigGrowthCostLA*(CCsapwood[j] - 1.0)/CCsapwood[j];
        }
       
        ///// C5. SAPWOOD GROWTH /////
        if(LAexpanded>0.0) {
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
        
        
        ///// C6. FINE ROOT GROWTH /////
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
        
        // Rcout<<growthCostLA<<" "<<growthCostSA<<" "<<growthCostFRB<< " "<< TotalLivingBiomass[j]<<"\n";
        GrowthCosts[j] +=(growthCostLA + twigGrowthCostLA + growthCostSA + growthCostFRB)/TotalLivingBiomass[j]; //growth cost in g gluc · gdry-1
        
        ///// C7a PARTIAL CARBON BALANCE: photosynthesis, maintenance respiration and growth /////
        double leafSugarMassDelta = leafAgG - leafRespDay;
        double sapwoodSugarMassDelta =  - twigResp - sapwoodResp - finerootResp; 
        double sapwoodStarchMassDelta =  - growthCostFRB - growthCostLA - twigGrowthCostLA - growthCostSA;
        sugarSapwood[j] += sapwoodSugarMassDelta/(Volume_sapwood[j]*glucoseMolarMass);
        starchSapwood[j] += sapwoodStarchMassDelta/(Volume_sapwood[j]*glucoseMolarMass);
        if(LAexpanded>0.0) sugarLeaf[j] += leafSugarMassDelta/(Volume_leaves[j]*glucoseMolarMass);
    
        ///// C7b PHLOEM TRANSPORT AND SUGAR-STARCH DYNAMICS /////  
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
          double B_resp_sapwood = SapwoodLivingStructBiomass[j] + TwigLivingStructuralBiomass[j] + sapwoodSugarMassStep;
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
      
      //add MR and SR to stand-level maintenance and synthesis respiration as g C·m-2
      standMaintenanceRespiration += (leafRespDay+twigResp+sapwoodResp+finerootResp)*(N[j]/10000.0)*((carbonMolarMass*6.0)/glucoseMolarMass);
      standSynthesisRespiration += (synthesisRespLA+ twigSynthesisRespLA +synthesisRespSA+synthesisRespFRB)*(N[j]/10000.0)*((carbonMolarMass*6.0)/glucoseMolarMass);
      
      ///// C8. LEAF SENESCENCE /////
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

      ///// C9. SAPWOOD AREA SENESCENCE /////
      //Define sapwood senescence as maximum of turnover and sapwood exceeding the target
      double propSASenescence = SRsapwood[j]*std::max(0.0,(tday-5.0)/20.0)/(1.0+15.0*exp(-0.01*H[j]));
      double deltaSASenescence = std::max(0.0, SA[j] - sapwoodAreaTarget[j]);
      propSASenescence = std::max(propSASenescence, deltaSASenescence/SA[j]);
      //For shrubs, convert SA senescence into mass of standing dead branches (g C/m2)
      if(ctype[j]=="shrub") {
        double sh_wood_biomass = (deltaSASenescence/SA[j])*(N[j]/10000.0)*AbovegroundWoodBiomass[j]*WoodC[j]; //Aboveground necromass
        Snag_smallbranches[j] += sh_wood_biomass;
      }
      
      ///// C10. FINE ROOT BIOMASS SENESCENCE /////
      NumericVector deltaFRBsenescence(nlayers, 0.0);
      for(int l=0;l<nlayers;l++) {
        double daySenescence = NA_REAL;
        if(transpirationMode=="Granier") daySenescence = SRfineroot[j]*std::max(0.0,(tday-5.0)/20.0);
        else daySenescence = SRfineroot[j]*std::max(0.0,(Tsoil[l]-5.0)/20.0);
        deltaFRBsenescence[l] = fineRootBiomass[j]*V(j,l)*daySenescence;
      }
      double senescenceFinerootLoss = sum(deltaFRBsenescence);
      if(exportLitter) {
        //From g/ind to g C/m2
        double fineroot_litter = rootCperDry*senescenceFinerootLoss*(N[j]/10000.0); 
        addFineRootLitter(speciesNames[j], fineroot_litter,
                          internalLitter, paramsLitterDecomposition,
                          internalSOC);
      }
      // if(j==(numCohorts-1)) Rcout<< j << " before translocation "<< sugarLeaf[j]<< " "<< starchLeaf[j]<<"\n";
      

      ///// 11. TRANSLOCATION (in mol gluc) of labile carbon due to tissue senescence /////
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

      ///// C12. ROOT EXUDATION and close labile carbon balance (non-subdaily carbon balance)
      if(!subdailyCarbonBalance) {
        //Excess sapwood starch carbon is lost as root exudation
        if(starchSapwood[j] > Starch_max_sapwood[j]) {
          RootExudation[j] = ((starchSapwood[j] - Starch_max_sapwood[j])*(Volume_sapwood[j]*glucoseMolarMass)/TotalLivingBiomass[j]);
          if(exportLitter) {
            //From g gluc/g ind to g C/m2
            double metabolicSoil = (6.0*carbonMolarMass/glucoseMolarMass) * RootExudation[j]*TotalLivingBiomass[j]*(N[j]/10000.0); 
            internalSOC["SoilMetabolic"] = internalSOC["SoilMetabolic"] + metabolicSoil;
          }
          starchSapwood[j] = Starch_max_sapwood[j];
        }
        //Labile carbon balance
        LabileCarbonBalance[j] = GrossPhotosynthesis[j] - MaintenanceRespiration[j] - GrowthCosts[j] - RootExudation[j];
      }
      
      //CHECK LABILE BIOMASS BALANCE (g gluc)
      //Final labile (g gluc)
      double finalLabile =  glucoseMolarMass*((sugarLeaf[j] + starchLeaf[j])*Volume_leaves[j] + (sugarSapwood[j] + starchSapwood[j])*Volume_sapwood[j]);
      double labileChange = finalLabile - initialLabile;
      double labileBalance_ggluc = LabileCarbonBalance[j]*TotalLivingBiomass[j];
      if(std::abs(labileChange - labileBalance_ggluc) > 0.1) {
        Rcout << j << " Initial labile " << initialLabile << " Final Labile " << finalLabile << " Labile change " << (finalLabile-initialLabile) << " A - MR - GC - RE " <<  LabileCarbonBalance[j]*TotalLivingBiomass[j]<<"\n";
        stop("Wrong labile biomass balance");
      }
        
      ///// C13. UPDATE INDIVIDUAL LEAF AREA, DEAD LEAF AREA, SAPWOOD AREA, FINE ROOT BIOMASS AND CONCENTRATION IN LABILE POOLS /////
      // Rcout<<"-update";
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
      //Decrease PLC due to new SA growth
      StemPLC[j] = std::max(0.0, StemPLC[j] - (deltaSAgrowth[j]/SA[j])); 
      LeafPLC[j] = std::max(0.0, LeafPLC[j] - (deltaLAgrowth[j]/LAexpanded)); 
      //Increase crown buds to new SA growth
      crownBudPercent[j] = std::min(100.0, crownBudPercent[j] + 100.0*(deltaSAgrowth[j]/SA[j]));
      
      ///// C14. UPDATE DERIVED HYDRAULIC PARAMETERS ///// 
      if(transpirationMode=="Granier") {
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
      
      
      ///// C15. CLOSE LEAF/FINE ROOT BIOMASS balance (g_ind)
      LeafBiomassBalance[j] = leafBiomassIncrement - senescenceLeafLoss;
      FineRootBiomassBalance[j] = finerootBiomassIncrement - senescenceFinerootLoss;
      
      //INDIVIDUAL-LEVEL OUTPUT
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
  
  /////////////////////////////////////////////////////////
  ///// D. PLANT MORTALITY AND FIRE EFFECTS BY COHORT /////
  /////////////////////////////////////////////////////////
  for(int j=0;j<numCohorts;j++){
    if(N[j] > 0.0) {
      
      //MORTALITY Death by carbon starvation or dessication
      double Ndead_day = 0.0, N_resprouting_stump_day = 0.0, Cover_resprouting_stump_day = 0.0;
      bool dynamicCohort = true;
      if((ctype[j] == "shrub") && (!shrubDynamics)) dynamicCohort = false;
      else if(ctype[j] =="herb") dynamicCohort = false;
      double stemSympRWC = NA_REAL;
      if(transpirationMode=="Granier") stemSympRWC = symplasticRelativeWaterContent(PlantPsi[j], StemPI0[j], StemEPS[j]);
      else stemSympRWC = symplasticRelativeWaterContent(psiSympStem[j], StemPI0[j], StemEPS[j]);
      if(NumericVector::is_na(stemSympRWC)) stop("Missing value for stem symp RWC");
      
      
      //Sapwood sugar relative to equilibrium, indicator of starvation
      double relativeSugarSapwood = (sugarSapwood[j]/equilibriumSapwoodSugarConc);
      if(dynamicCohort) {
        String cause = "undertermined";
        //Determine fire severity if fire occurred
        double burnRatioLeaves = 0.0, burnRatioBuds = 0.0;
        bool abovegroundFireSurvival = true;
        if(fireOccurrence) {
          double rho_air = meteoland::utils_airDensity(tmax, Patm);
          double foliar_factor = leafThermalFactor(SLA[j]);
          double Ib_surf = fireBehavior["I_b_surface [kW/m]"];
          double t_res_surf = fireBehavior["t_r_surface [s]"];
          double t_r_crown = fireBehavior["t_r_crown [s]"];
          double fm_dead = fireBehavior["DFMC [%]"];
          double Ic_ratio = fireBehavior["Ic_ratio"];
          double cbh = H[j]*(1.0 - CR[j]);
          double bark_thickness = 1.0;
          double Hn_leaves = 0.0, Hn_buds = 0.0;
          double xn = 0.0;
          if(ctype[j] == "shrub") {
            bark_thickness = BTsh[j]*0.1; // from mm to cm
          } else {
            bark_thickness = Abt[j]*pow(DBH[j],Bbt[j])*0.1;// from mm to cm
          } 
          //Determine foliage/bud burn
          if(!NumericVector::is_na(Ib_surf) && !NumericVector::is_na(t_res_surf)) {
            Hn_leaves =100.0*necrosisHeight(Ib_surf, t_res_surf, foliar_factor, tmax, rho_air); //Necrosis height (cm)
            Hn_buds = 100.0*necrosisHeight(Ib_surf, t_res_surf, 0.130, tmax, rho_air); //Bud necrosis height (cm)
            burnRatioLeaves = leafAreaProportion(0.0, Hn_leaves, cbh, H[j]);
            burnRatioBuds = leafAreaProportion(0.0, Hn_buds, cbh, H[j]);
            // Rcout << " tmax " << tmax << " rho_air " << rho_air <<" foliar_factor "<< foliar_factor << " Ib_surf "<< Ib_surf << " t_res_surf " << t_res_surf<< " foliar_factor "<< foliar_factor << " Hn_leaves "<<Hn_leaves << " br_leaves "<< burnRatioLeaves<< " Hn_buds "<<Hn_buds << " br_buds "<< burnRatioBuds<<"\n";
          }
          //Determine crown fire or torching effects
          if(!NumericVector::is_na(Ib_surf)) {
            double canopyFMC = (LFMC[j]*(1.0 - StemPLC[j]) + fm_dead*StemPLC[j]);
            double Ib_crit = criticalFirelineIntensity(cbh/100.0, canopyFMC);
            // Rcout << "Ic_ratio "<< Ic_ratio <<" Ib_crit "<<Ib_crit<< " Ib_surf "<< Ib_surf<<"\n";
            if((Ic_ratio > 1.0) || (Ib_surf > Ib_crit)) {
              burnRatioLeaves = 1.0;
              double Tc = necrosisCriticalTemperature(t_r_crown, 0.130 , tmax, rho_air);
              if(Tc < 900.0) burnRatioBuds = 1.0;
            }
            // Rcout << "br_leaves "<< burnRatioLeaves<< " br_buds "<< burnRatioBuds<<"\n";
          }
          //Surface fire effects on cambium
          if(!NumericVector::is_na(Ib_surf) && !NumericVector::is_na(t_res_surf)) {
            double bark_diff = barkThermalDiffusivity(fm_dead, 500.0, tmax);
            xn =  radialBoleNecrosis(Ib_surf, t_res_surf, bark_diff, tmax, rho_air);
            // Rcout << "xn "<< xn<< " xa "<<bark_thickness<<"\n";
          }
          
          crownBudPercent[j] = crownBudPercent[j]*(1.0 - burnRatioBuds);
          abovegroundFireSurvival = (burnRatioBuds < 1.0) && (bark_thickness > xn);
          
          //If does not survive set ratios to 1.0
          if(!abovegroundFireSurvival) {
            burnRatioLeaves = 1.0;
            burnRatioBuds = 1.0;
          }
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
          N_resprouting_stump_day = Ndead_day*RespFire[j];
        } else if(cause == "dessication") {
          N_dessication[j] = N_dessication[j] + Ndead_day;
          N_resprouting_stump_day = Ndead_day*RespDist[j];
        } else if(ctype[j]=="tree") { // Self-thinning occurring in tree cohorts
          if(DBH[j] < IngrowthTreeDBH[j]) {
            double b_st = log(RecrTreeDensity[j]/IngrowthTreeDensity[j])/log(RecrTreeDBH[j]/IngrowthTreeDBH[j]);
            double a_st = IngrowthTreeDensity[j]/pow(IngrowthTreeDBH[j], b_st);
            double N_st = a_st*pow(DBH[j], b_st);
            double N_dead_selfthinning = N[j] - std::min(N[j], N_st);
            // Rcout<< b_st<< " "<< a_st<< " "<< N_st<< " "<< N_dead_selfthinning<<"\n";
            Ndead_day = Ndead_day + N_dead_selfthinning;
          }
        }
        N_resprouting_stumps[j] = N_resprouting_stumps[j] + N_resprouting_stump_day;
        // Rcout << Ndead_day <<"\n";
        N[j] = N[j] - Ndead_day;
        N_dead[j] = N_dead[j] + Ndead_day;
        if(ctype[j]=="shrub") {
          Cover[j] = std::max(0.0, Cover[j] - Cdead_day);
          Cover_dead[j] = Cover_dead[j] + Cdead_day;
          if(cause == "starvation") {
            Cover_starvation[j] = Cover_starvation[j] + Cdead_day;
          } else if(cause == "dessication") {
            Cover_dessication[j] = Cover_dessication[j] + Cdead_day;
            Cover_resprouting_stump_day =  Cdead_day*RespDist[j];
          } else if(cause == "burnt") {
            Cover_burnt[j] = Cover_burnt[j] + Cdead_day;
            Cover_resprouting_stump_day =  Cdead_day*RespFire[j];
          }
        }
        Cover_resprouting_stumps[j] = Cover_resprouting_stumps[j] + Cover_resprouting_stump_day;
        

        //Update LAI live, LAI dead and LAI expanded as a result of density decrease and foliage burn
        double LAI_senescence = LAI_expanded[j]*Ndead_day/10000.0;
        double LAI_senescence_burned = LAI_senescence*burnRatioLeaves;
        double LAI_senescence_notburned = LAI_senescence*(1.0 - burnRatioLeaves);
        double LAI_live_burned = (LAI_expanded[j] - LAI_senescence) * burnRatioLeaves;
        double LAI_dead_burned = LAI_dead[j]*burnRatioLeaves;
        LAI_dead[j] = LAI_dead[j] + LAI_senescence_notburned;
        LAI_expanded[j] = LAI_expanded[j] - LAI_senescence - LAI_live_burned;
        if(!abovegroundFireSurvival) {
          LAI_live[j] = 0.0;
        }
        
        // FOLIAGE BURNING EMISSIONS (in g C / m2) 
        // from m2/m2 to g C/m2
        double leaf_burnt_dry = 1000.0*(LAI_live_burned + LAI_dead_burned + LAI_senescence_burned)/SLA[j];
        double leaf_burnt_C = leafCperDry*leaf_burnt_dry;
        double twig_burnt_C = WoodC[j]*leaf_burnt_dry/(r635[j] - 1.0);
        fireCombustion +=leaf_burnt_C + twig_burnt_C;
        
        //ADD STANDING DEADWOOD (in g C / m2) according to non-resprouting
        double snag_wood_biomass = (Ndead_day/10000.0)*AbovegroundWoodBiomass[j]*WoodC[j]; //Aboveground necromass (includes burning and dessication)
        if(cause != "burnt") {
          if(ctype[j]=="tree") { // 20% goes to small branches (diameter < 7.5) SHOULD depend on size allometry
            Snag_smallbranches[j] += 0.2*snag_wood_biomass; 
            Snag_largewood[j] += 0.8*snag_wood_biomass; 
          } else if(ctype[j] =="shrub") { // For shrubs 100% goes to small branches
            Snag_smallbranches[j] += snag_wood_biomass;
          }
        } else { // burnt
          if(ctype[j]=="tree") { 
            fireCombustion += 0.2*snag_wood_biomass; //20% goes to fire emissions
            Snag_largewood[j] += 0.8*snag_wood_biomass; 
          } else if(ctype[j] =="shrub") { // For shrubs 100% goes to fire emissions
            fireCombustion += snag_wood_biomass;
          }
        }
        //Fine root and coarse root litter corresponding to non-resprouting trees
        if(exportLitter) {
          double coarseroot_litter = ((Ndead_day - N_resprouting_stump_day)/10000.0)*BelowgroundWoodBiomass[j]*WoodC[j]; //Only adds non-resprouting death biomass
          addCoarseRootLitter(speciesNames[j], coarseroot_litter, internalLitter);
          double fineroot_litter = rootCperDry*((Ndead_day - N_resprouting_stump_day)/10000.0)*fineRootBiomass[j]; 
          addFineRootLitter(speciesNames[j], fineroot_litter,
                            internalLitter, paramsLitterDecomposition,
                            internalSOC);
        }
      }
    }
  }
  
  ///////////////////////////////////////////
  ///// E1. UPDATE STRUCTURAL VARIABLES /////
  ///////////////////////////////////////////
  updateStructuralVariables(x, deltaSAgrowth);
  
  ///////////////////////////////////////////////////////////////////////////////////////////////
  ///// E2. UPDATE SAPWOOD AREA AND FINEROOT BIOMASS TARGETS AND RECALCULATE CONCENTRATIONS /////
  ///////////////////////////////////////////////////////////////////////////////////////////////
  for(int j=0;j<numCohorts;j++){
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

    //OUTPUT VARIABLES
    OutputDBH[j] = DBH[j];
    OutputHeight[j] = H[j];
    PlantSugarLeaf[j] = sugarLeaf[j];
    PlantStarchLeaf[j] = starchLeaf[j];
    PlantSugarSapwood[j] = sugarSapwood[j];
    PlantStarchSapwood[j] = starchSapwood[j];
  }
  
  ////////////////////////////////////////
  ///// E3. UPDATE PLANT WATER POOLS /////
  ////////////////////////////////////////
  if(plantWaterPools) {
    NumericVector poolProportions = Rcpp::as<Rcpp::NumericVector>(belowdf["poolProportions"]);
    for(int j=0;j<numCohorts;j++) poolProportions[j] = LAI_live[j]/sum(LAI_live);
    //Update RHOP
    List newRHOP;
    if(rhizosphereOverlap=="none") newRHOP = nonoverlapHorizontalProportions(V);
    else newRHOP =equaloverlapHorizontalProportions(poolProportions, V);
    // else newRHOP = horizontalProportions(poolProportions, CRSV, N, V, widths, rfc);
    for(int j=0;j<numCohorts;j++) RHOP[j] = newRHOP[j];
  }
  
  ////////////////////////////////////////////
  ///// E4. CLOSE PLANT BIOMASS BALANCE //////
  ////////////////////////////////////////////
  closePlantBiomassBalance(initialFinalCC, plantBiomassBalance, x,
                           LabileCarbonBalance, LeafBiomassBalance, FineRootBiomassBalance);
  
  ///////////////////////////////////
  ///// F. CARBON DECOMPOSITION /////
  ///////////////////////////////////
  double heterotrophicRespiration = 0.0;
  if(exportLitter) {
    NumericVector baseAnnualRates = control["decompositionAnnualBaseRates"];
    double annualTurnoverRate = control["decompositionAnnualTurnoverRate"];
    List commDecomp = internalCommunication["decomposition"];
    double soilPH = 7.0;
    if(soil.containsElementNamed("pH")) {
      NumericVector pH = soil["pH"]; 
      if(!NumericVector::is_na(pH[0])) soilPH = pH[0];
    }
    NumericVector soilTheta = theta(soil, soilFunctions);
    NumericVector soilThetaSAT = thetaSAT(soil, soilFunctions);
    double vpatm = meteoland::utils_averageDailyVP(tmin, tmax, rhmin, rhmax);
    double tday = meteoland::utils_averageDaylightTemperature(tmin, tmax);
    double vpsat = meteoland::utils_saturationVP(tday);
    double rhmean = 100.0*std::max(1.0, vpatm/vpsat);
    heterotrophicRespiration = DAYCENTInner(commDecomp, internalSnags, internalLitter, internalSOC,
                                            paramsLitterDecomposition,
                                            baseAnnualRates, annualTurnoverRate,
                                            tday, rhmean, 
                                            sand[0], clay[0], Tsoil[0], soilTheta[0]/soilThetaSAT[0], soilPH);
  }
  
  
  ///////////////////////////////////////////////
  ///// G. FIRE EFFECTS ON SNAGS AND LITTER /////
  ///////////////////////////////////////////////
  if(fireOccurrence) {
    NumericVector litter_leaves = internalLitter["Leaves"];
    NumericVector litter_twigs = internalLitter["Twigs"];
    NumericVector litter_smallbranches = internalLitter["SmallBranches"];
    fireCombustion += sum(litter_leaves) + sum(litter_twigs) + sum(litter_smallbranches);
    for(int i=0;i<litter_leaves.size();i++) {
      litter_leaves[i] = 0.0;
      litter_twigs[i] = 0.0;
      litter_smallbranches[i] = 0.0;
    }
  }
  
  ///////////////////////////////////////////////
  ///// H. CLOSE STAND-LEVEL CARBON BALANCE /////
  ///////////////////////////////////////////////
  NumericVector standCB = modelOutput["CarbonBalance"];
  standCB["GrossPrimaryProduction"] = standGrossPrimaryProduction;
  standCB["MaintenanceRespiration"] = standMaintenanceRespiration;
  standCB["SynthesisRespiration"] = standSynthesisRespiration;
  standCB["NetPrimaryProduction"] = standGrossPrimaryProduction - standMaintenanceRespiration - standSynthesisRespiration;
  standCB["FireCombustion"] = fireCombustion;
  standCB["HeterotrophicRespiration"] = heterotrophicRespiration;
  standCB["NetEcosystemProduction"] = standCB["NetPrimaryProduction"] - fireCombustion - heterotrophicRespiration; 
  
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

//' Single-day forest growth
//'
//' Function \code{growth_day} performs water and carbon balance for a single day.
//' 
//' @param x An object of class \code{\link{growthInput}}.
//' @param date Date as string "yyyy-mm-dd".
//' @param meteovec A named numerical vector with weather data. See variable names in parameter \code{meteo} of \code{\link{spwb}}.
//' @param latitude Latitude (in degrees).
//' @param elevation,slope,aspect Elevation above sea level (in m), slope (in degrees) and aspect (in degrees from North). 
//' @param runon Surface water amount running on the target area from upslope (in mm).
//' @param lateralFlows Lateral source/sink terms for each soil layer (interflow/to from adjacent locations) as mm/day.
//' @param waterTableDepth Water table depth (in mm). When not missing, capillarity rise will be allowed if lower than total soil depth.
//' @param modifyInput Boolean flag to indicate that the input \code{x} object is allowed to be modified during the simulation.
//' 
//' @details
//' Detailed model description is available in the medfate book. 
//' 
//' Forest growth simulations allow using different sub-models for bulk soil water flows and different sub-models of transpiration and photosynthesis (see details in \code{\link{spwb_day}}). 
//' 
//' @return
//' Function \code{growth_day()} returns a list of class \code{growth_day} with the 
//' same elements as \code{\link{spwb_day}} and the following:
//' \itemize{
//'   \item{\code{"CarbonBalance"}: A vector of different stand-level carbon balance components (gross primary production, maintenance respiration, synthesis respiration, net primary production, heterotrophic respiration and net ecosystem exchange), all in g C · m-2.}
//'   \item{\code{"LabileCarbonBalance"}: A data frame with labile carbon balance results for plant cohorts, with elements:}
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
//'   \item{\code{"PlantBiomassBalance"}: A data frame with plant biomass balance results for plant cohorts, with elements:}
//'   \itemize{
//'     \item{\code{"StructuralBiomassBalance"}: Daily structural biomass balance (g dry · m-2).}
//'     \item{\code{"LabileBiomassBalance"}: Daily labile biomass balance (g dry · m-2).}
//'     \item{\code{"PlantBiomassBalance"}: Daily plant biomass balance, i.e. labile change + structural change (g dry · m-2).}
//'     \item{\code{"MortalityBiomassLoss"}: Biomass loss due to mortality (g dry · m-2).}    
//'     \item{\code{"CohortBiomassBalance"}: Daily cohort biomass balance (including mortality) (g dry · m-2).}
//'   }
//'   \item{\code{"PlantStructure"}: A data frame with area and biomass values for compartments of plant cohorts, with elements:}
//'   \itemize{
//'     \item{\code{"LeafBiomass"}: Leaf structural biomass (in g dry) for an average individual of each plant cohort.}
//'     \item{\code{"SapwoodBiomass"}: Sapwood structural biomass (in g dry) for an average individual of each plant cohort.}
//'     \item{\code{"FineRootBiomass"}: Fine root biomass (in g dry) for an average individual of each plant cohort.}
//'     \item{\code{"LeafArea"}: Leaf area (in m2) for an average individual of each plant cohort.}
//'     \item{\code{"SapwoodArea"}: Sapwood area (in cm2) for an average individual of each plant cohort.}
//'     \item{\code{"FineRootArea"}: Fine root area (in m2) for an average individual of each plant cohort.}
//'     \item{\code{"HuberValue"}: Sapwood area to (target) leaf area (in cm2/m2).}
//'     \item{\code{"RootAreaLeafArea"}: The ratio of fine root area to (target) leaf area (in m2/m2).}
//'     \item{\code{"DBH"}: Diameter at breast height (in cm) for an average individual of each plant cohort.}
//'     \item{\code{"Height"}: Height (in cm) for an average individual of each plant cohort.}
//'   }
//'   \item{\code{"GrowthMortality"}: A data frame with growth and mortality rates for plant cohorts, with elements:}
//'   \itemize{
//'     \item{\code{"LAgrowth"}: Leaf area growth (in m2·day-1) for an average individual of each plant cohort.}
//'     \item{\code{"SAgrowth"}: Sapwood area growth rate (in cm2·day-1) for an average individual of each plant cohort.}
//'     \item{\code{"FRAgrowth"}: Fine root area growth (in m2·day-1) for an average individual of each plant cohort.}
//'     \item{\code{"StarvationRate"}: Mortality rate from starvation (ind/d-1).}
//'     \item{\code{"DessicationRate"}: Mortality rate from dessication (ind/d-1).}
//'     \item{\code{"MortalityRate"}: Mortality rate (any cause) (ind/d-1).}
//'   }
//' }
//'   
//' @references
//' De \enc{Cáceres}{Caceres} M, \enc{Martínez}{Martinez}-Vilalta J, Coll L, Llorens P, Casals P, Poyatos R, Pausas JG, Brotons L. (2015) Coupling a water balance model with forest inventory data to predict drought stress: the role of forest structural changes vs. climate changes. Agricultural and Forest Meteorology 213: 77-90 (doi:10.1016/j.agrformet.2015.06.012).
//' 
//' De \enc{Cáceres}{Caceres} M, Mencuccini M, Martin-StPaul N, Limousin JM, Coll L, Poyatos R, Cabon A, Granda V, Forner A, Valladares F, \enc{Martínez}{Martinez}-Vilalta J (2021) Unravelling the effect of species mixing on water use and drought stress in holm oak forests: a modelling approach. Agricultural and Forest Meteorology 296 (doi:10.1016/j.agrformet.2020.108233).
//' 
//' Granier A, \enc{Bréda}{Breda} N, Biron P, Villette S (1999) A lumped water balance model to evaluate duration and intensity of drought constraints in forest stands. Ecol Modell 116:269–283. https://doi.org/10.1016/S0304-3800(98)00205-1.
//' 
//' Ruffault J, Pimont F, Cochard H, Dupuy JL, Martin-StPaul N (2022) 
//' SurEau-Ecos v2.0: a trait-based plant hydraulics model for simulations of plant water status and drought-induced mortality at the ecosystem level.
//' Geoscientific Model Development 15, 5593-5626 (doi:10.5194/gmd-15-5593-2022).
//' 
//' Sperry, J. S., M. D. Venturas, W. R. L. Anderegg, M. Mencuccini, D. S. Mackay, Y. Wang, and D. M. Love. 2017. Predicting stomatal responses to the environment from the optimization of photosynthetic gain and hydraulic cost. Plant Cell and Environment 40, 816-830 (doi: 10.1111/pce.12852).
//' 
//' @author
//' \itemize{
//'   \item{Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF}
//'   \item{Nicolas Martin-StPaul, URFM-INRAE}
//' }
//' 
//' @seealso
//' \code{\link{spwb_day}}, \code{\link{growthInput}}, \code{\link{growth}},
//' \code{\link{plot.growth_day}}  
//' 
//' @examples
//' #Load example daily meteorological data
//' data(examplemeteo)
//' 
//' #Load example plot plant data
//' data(exampleforest)
//' 
//' #Default species parameterization
//' data(SpParamsMED)
//' 
//' #Define soil parameters
//' examplesoil <- defaultSoilParams(4)
//' 
//' # Day to be simulated
//' d <- 100
//' meteovec <- unlist(examplemeteo[d,-1])
//' date <- as.character(examplemeteo$dates[d])
//' 
//' #Simulate water and carbon balance for one day only (Granier mode)
//' control <- defaultControl("Granier")
//' x4  <- growthInput(exampleforest,examplesoil, SpParamsMED, control)
//' sd4 <- growth_day(x4, date, meteovec,
//'                 latitude = 41.82592, elevation = 100, slope=0, aspect=0)
//' 
//' #Simulate water and carbon balance for one day only (Sperry mode)
//' control <- defaultControl("Sperry")
//' x5  <- growthInput(exampleforest,examplesoil, SpParamsMED, control)
//' sd5 <- growth_day(x5, date, meteovec,
//'                 latitude = 41.82592, elevation = 100, slope=0, aspect=0)
//' 
//' #Simulate water and carbon balance for one day only (Sureau mode)
//' control <- defaultControl("Sureau")
//' x6  <- growthInput(exampleforest,examplesoil, SpParamsMED, control)
//' sd6 <- growth_day(x6, date, meteovec,
//'                 latitude = 41.82592, elevation = 100, slope=0, aspect=0)
//' 
//' @name growth_day
// [[Rcpp::export("growth_day")]]
List growthDay(List x, CharacterVector date, NumericVector meteovec, 
               double latitude, double elevation, double slope = NA_REAL, double aspect = NA_REAL,  
               double runon = 0.0, Nullable<NumericVector> lateralFlows = R_NilValue, double waterTableDepth = NA_REAL, 
               bool modifyInput = true) {
  
  //Check if input version is lower than current medfate version. If so, try to complete fields
  if(isLowerVersion(x)) growthInputVersionUpdate(x);
  
  //Instance communication structures
  List internalCommunication = instanceCommunicationStructures(x, "growth");
  growthDay_inner(internalCommunication, x, date, meteovec,
                                     latitude, elevation, slope, aspect,
                                     runon, lateralFlows, waterTableDepth,
                                     modifyInput);

  List modelOutput = copyModelOutput(internalCommunication, x, "growth");
  return(modelOutput);
}
