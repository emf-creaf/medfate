// [[Rcpp::interfaces(r,cpp)]]
#define STRICT_R_HEADERS
#include <Rcpp.h>
#include <numeric>
#include "lightextinction.h"
#include "phenology.h"
#include "biophysicsutils.h"
#include "forestutils.h"
#include "tissuemoisture.h"
#include "hydraulics.h"
#include "hydrology.h"
#include "carbon.h"
#include "root.h"
#include "woodformation.h"
#include "soil.h"
#include "spwb.h"
#include <meteoland.h>
using namespace Rcpp;

const double Q10_resp = 2.0;

// [[Rcpp::export("mortality_dailyProbability")]]
double dailyMortalityProbability(double mortalityBaselineRate, 
                                 double stressValue, double stressThreshold, 
                                 bool allowStress = true,
                                 double minValue = 0.0, double slope = 1.0) {
  double P_day = 1.0 - exp(log(1.0 - mortalityBaselineRate)/356.0);
  if(allowStress) {
    double P_stress = (1.0-exp(slope*(stressValue - stressThreshold)))/(1.0-exp(slope*(minValue-stressThreshold)));
    P_day = std::max(P_day, P_stress);
  }
  return(P_day);
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
// // [[Rcpp::export(".growth_defoliationFraction")]]
// double defoliationFraction(double conc, double threshold) {
//   double k =-10.0;
//   return(std::max(0.0,(exp(k*conc)-exp(k*threshold))/(1.0-exp(k*threshold))));
// }

DataFrame initPlantBiomassBalance(DataFrame ccIni, DataFrame above) {

  NumericVector Nprev = above["N"];
  int numCohorts = Nprev.length();
  
  //Initial Biomass compartments
  NumericVector SapwoodBiomass = Rcpp::as<Rcpp::NumericVector>(ccIni["SapwoodStructuralBiomass"]);
  NumericVector TotalBiomass = Rcpp::as<Rcpp::NumericVector>(ccIni["TotalBiomass"]);
  NumericVector TotalLivingBiomass= Rcpp::as<Rcpp::NumericVector>(ccIni["TotalLivingBiomass"]);
  NumericVector CohortBiomass = TotalBiomass*(Nprev/10000.0);
  NumericVector LabileBiomass = Rcpp::as<Rcpp::NumericVector>(ccIni["LabileBiomass"]);
  NumericVector StructuralBiomass = Rcpp::as<Rcpp::NumericVector>(ccIni["StructuralBiomass"]);
  
  NumericVector LabileBiomassBalance(numCohorts,0.0), StructuralBiomassBalance(numCohorts,0.0), PlantBiomassBalance(numCohorts,0.0), MortalityBiomassLoss(numCohorts,0.0), CohortBiomassBalance(numCohorts,0.0);
  NumericVector StructuralBiomassChange(numCohorts,0.0), LabileBiomassChange(numCohorts,0.0), PlantBiomassChange(numCohorts,0.0), CohortBiomassChange(numCohorts,0.0);
  
  DataFrame plantBiomassBalance = DataFrame::create(_["InitialDensity"] = clone(Nprev),
                                               _["InitialSapwoodBiomass"] = SapwoodBiomass,
                                               _["InitialStructuralBiomass"] = StructuralBiomass,
                                               _["StructuralBiomassBalance"] = StructuralBiomassBalance,
                                               _["StructuralBiomassChange"] = StructuralBiomassChange,
                                               _["InitialLabileBiomass"] = LabileBiomass,
                                               _["LabileBiomassBalance"] = LabileBiomassBalance,
                                               _["LabileBiomassChange"] = LabileBiomassChange,
                                               _["InitialLivingPlantBiomass"] = TotalLivingBiomass,
                                               _["InitialPlantBiomass"] = TotalBiomass,
                                               _["PlantBiomassBalance"] = PlantBiomassBalance,
                                               _["PlantBiomassChange"] = PlantBiomassChange,
                                               _["MortalityBiomassLoss"] = MortalityBiomassLoss,
                                               _["InitialCohortBiomass"] = CohortBiomass,
                                               _["CohortBiomassBalance"] = CohortBiomassBalance,
                                               _["CohortBiomassChange"] = CohortBiomassChange);
  plantBiomassBalance.attr("row.names") = above.attr("row.names");
  return(plantBiomassBalance);  
}


void closePlantBiomassBalance(DataFrame plantBiomassBalance, List x,
                         NumericVector LabileCarbonBalance,
                         NumericVector LeafBiomassBalance,
                         NumericVector FineRootBiomassBalance) {
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  NumericVector Nfinal = above["N"];
  int numCohorts = Nfinal.length();
  
  DataFrame ccFin = carbonCompartments(x, "g_ind");
  
  NumericVector finalSapwoodBiomass_ind= Rcpp::as<Rcpp::NumericVector>(ccFin["SapwoodStructuralBiomass"]);
  NumericVector plantFinalBiomass_ind = Rcpp::as<Rcpp::NumericVector>(ccFin["TotalBiomass"]);
  NumericVector cohortFinalBiomass_m2 = Rcpp::as<Rcpp::NumericVector>(ccFin["TotalBiomass"])*(Nfinal/10000.0);
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
  for(int j; j<numCohorts;j++) {
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
    CohortBiomassChange[j] = cohortFinalBiomass_m2[j] - InitialCohortBiomass[j];
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
  NumericVector SA = above["SA"];
  
  //Allometric parameters
  DataFrame paramsAllometries = Rcpp::as<Rcpp::DataFrame>(x["paramsAllometries"]);
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
  NumericVector r635  = paramsAnatomy["r635"];
  
  DataFrame internalAllocation = Rcpp::as<Rcpp::DataFrame>(x["internalAllocation"]);
  NumericVector leafAreaTarget = internalAllocation["leafAreaTarget"];
  
  NumericVector deltaDBH(numCohorts, 0.0);
  for(int j=0;j<numCohorts; j++) {
    if(!NumericVector::is_na(DBH[j])) {
      deltaDBH[j] = 2.0*sqrt(pow(DBH[j]/2.0,2.0)+(deltaSAgrowth[j]/M_PI)) - DBH[j];
      DBH[j] = DBH[j] + deltaDBH[j];
    }
  }
  NumericVector L = parcohortC(H, LAI_live, LAI_dead, kPAR, CR);
  for(int j=0;j<numCohorts; j++) {
    if(!NumericVector::is_na(DBH[j]) && N[j]>0.0) {
      double fHmod = std::max(0.0,std::min(1.0,(1.0-((H[j]-137.0)/(Hmax[j]-137.0)))));
      double fHD = (fHDmin[j]*(L[j]/100.0) + fHDmax[j]*(1.0-(L[j]/100.0)))*fHmod;
      H[j] = H[j] + fHD*deltaDBH[j];
    }
  }
  NumericVector crNew = treeCrownRatioMED(N, DBH, H, Acw, Bcw, Acr, B1cr, B2cr, B3cr, C1cr, C2cr);
  for(int j=0;j<numCohorts; j++) {
    if(!NumericVector::is_na(DBH[j]) && N[j]>0.0) {
      CR[j] = crNew[j];
    }
  }
  
  //Shrub variables
  if(shrubDynamics) {
    for(int j=0;j<numCohorts; j++) {
      if(NumericVector::is_na(DBH[j]) && N[j]>0.0) {
        double Wleaves = leafAreaTarget[j]/SLA[j];  //Calculates the biomass (kg dry weight) of leaves
        double PV = pow(Wleaves*r635[j]/Absh[j], 1.0/Bbsh[j]); //Calculates phytovolume (in m3/ind)
        H[j] = pow(1e6*PV/Aash[j], 1.0/(1.0+Bash[j])); //Updates shrub height
        // Rcout<< Wleaves << " " << PV << " " << H[j]<<"\n";
        if(H[j]> Hmax[j]) { //Limit height (and update the former variables)
          H[j] = Hmax[j];
        }
        Cover[j] = (N[j]*Aash[j]*pow(H[j],Bash[j])/1e6); //Updates shrub cover
      }
    }
  }
}
List growthDay1(List x, NumericVector meteovec, 
                double elevation = NA_REAL, 
                double runon=0.0, bool verbose = false) {
  
  
  //Soil-plant water balance
  List spwbOut = spwbDay1(x, meteovec, 
                          elevation, 
                          runon, verbose);
  
  //Control params
  List control = x["control"];  
  
  String mortalityMode = control["mortalityMode"];
  double mortalityBaselineRate = control["mortalityBaselineRate"];
  double mortalityRelativeSugarThreshold= control["mortalityRelativeSugarThreshold"];
  double mortalityRWCThreshold= control["mortalityRWCThreshold"];
  bool allowDessication = control["allowDessication"];
  bool allowStarvation = control["allowStarvation"];
  bool allowDefoliation = control["allowDefoliation"];
  bool sinkLimitation = control["sinkLimitation"];
  bool shrubDynamics = control["shrubDynamics"];
  String allocationStrategy = control["allocationStrategy"];
  String cavitationRefill = control["cavitationRefill"];
  bool plantWaterPools = control["plantWaterPools"];
  double nonSugarConcentration = control["nonSugarConcentration"];
  List equilibriumOsmoticConcentration  = control["equilibriumOsmoticConcentration"];
  double equilibriumLeafTotalConc = equilibriumOsmoticConcentration["leaf"];
  double equilibriumSapwoodTotalConc = equilibriumOsmoticConcentration["sapwood"];
  double minimumRelativeSugarForGrowth = control["minimumRelativeSugarForGrowth"];
  List constructionCosts = control["constructionCosts"];
  double leaf_CC = constructionCosts["leaf"];
  double sapwood_CC = constructionCosts["sapwood"];
  double fineroot_CC = constructionCosts["fineroot"];

  //Cohort info
  DataFrame cohorts = Rcpp::as<Rcpp::DataFrame>(x["cohorts"]);
  IntegerVector SP = Rcpp::as<Rcpp::IntegerVector>(cohorts["SP"]);
  int numCohorts = SP.size();
  
  //Soil

  //Weather
  double tday = meteovec["tday"];
  
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
  NumericVector SA = above["SA"];

  //Belowground parameters  
  DataFrame belowdf = Rcpp::as<Rcpp::DataFrame>(x["below"]);
  NumericVector Z95 = Rcpp::as<Rcpp::NumericVector>(belowdf["Z95"]);
  NumericVector Z50 = Rcpp::as<Rcpp::NumericVector>(belowdf["Z50"]);
  NumericVector fineRootBiomass = Rcpp::as<Rcpp::NumericVector>(belowdf["fineRootBiomass"]);
  List belowLayers = Rcpp::as<Rcpp::List>(x["belowLayers"]);
  NumericMatrix V = belowLayers["V"];
  NumericMatrix L = Rcpp::as<Rcpp::NumericMatrix>(belowLayers["L"]);
  
  //Internal state variables
  DataFrame internalWater = Rcpp::as<Rcpp::DataFrame>(x["internalWater"]);
  //Values at the end of the day (after calling spwb)
  NumericVector StemPLC = Rcpp::as<Rcpp::NumericVector>(internalWater["StemPLC"]);
  NumericVector PlantPsi = Rcpp::as<Rcpp::NumericVector>(internalWater["PlantPsi"]);
  
  DataFrame internalCarbon = Rcpp::as<Rcpp::DataFrame>(x["internalCarbon"]);
  NumericVector sugarLeaf = internalCarbon["sugarLeaf"]; //Concentrations assuming RWC = 1
  NumericVector starchLeaf = internalCarbon["starchLeaf"];
  NumericVector sugarSapwood = internalCarbon["sugarSapwood"];
  NumericVector starchSapwood = internalCarbon["starchSapwood"];

  DataFrame internalMortality = Rcpp::as<Rcpp::DataFrame>(x["internalMortality"]);
  NumericVector N_dead = internalMortality["N_dead"];
  NumericVector Cover_dead = internalMortality["Cover_dead"];
  
  DataFrame internalAllocation = Rcpp::as<Rcpp::DataFrame>(x["internalAllocation"]);
  NumericVector allocationTarget = internalAllocation["allocationTarget"];
  NumericVector leafAreaTarget = internalAllocation["leafAreaTarget"];
  
  DataFrame internalPhenology = Rcpp::as<Rcpp::DataFrame>(x["internalPhenology"]);
  LogicalVector leafUnfolding = internalPhenology["leafUnfolding"];
  LogicalVector budFormation = internalPhenology["budFormation"];
  LogicalVector leafSenescence = internalPhenology["leafSenescence"];
  
  List stand = spwbOut["Stand"];
  DataFrame Plants = Rcpp::as<Rcpp::DataFrame>(spwbOut["Plants"]);
  NumericVector Ag = Plants["GrossPhotosynthesis"];

  //Data from spwb
  // NumericVector LeafRWC = Plants["LeafRWC"];
  // NumericVector StemRWC = Plants["StemRWC"];

  //Anatomy parameters
  DataFrame paramsAnatomy = Rcpp::as<Rcpp::DataFrame>(x["paramsAnatomy"]);
  NumericVector Hmed = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["Hmed"]);
  NumericVector SLA = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["SLA"]);
  NumericVector Al2As = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["Al2As"]);
  NumericVector WoodDensity = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["WoodDensity"]);
  NumericVector LeafDensity = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["LeafDensity"]);
  NumericVector conduit2sapwood = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["conduit2sapwood"]);
  
  //Growth parameters
  DataFrame paramsGrowth = Rcpp::as<Rcpp::DataFrame>(x["paramsGrowth"]);
  NumericVector WoodC = Rcpp::as<Rcpp::NumericVector>(paramsGrowth["WoodC"]);
  NumericVector RERleaf = Rcpp::as<Rcpp::NumericVector>(paramsGrowth["RERleaf"]);
  NumericVector RERsapwood = Rcpp::as<Rcpp::NumericVector>(paramsGrowth["RERsapwood"]);
  NumericVector RERfineroot = Rcpp::as<Rcpp::NumericVector>(paramsGrowth["RERfineroot"]);
  NumericVector RGRleafmax = Rcpp::as<Rcpp::NumericVector>(paramsGrowth["RGRleafmax"]);
  NumericVector RGRsapwoodmax = Rcpp::as<Rcpp::NumericVector>(paramsGrowth["RGRsapwoodmax"]);
  NumericVector SRsapwood = Rcpp::as<Rcpp::NumericVector>(paramsGrowth["SRsapwood"]);

  //Phenology parameters
  DataFrame paramsPhenology = Rcpp::as<Rcpp::DataFrame>(x["paramsPhenology"]);
  CharacterVector phenoType = Rcpp::as<Rcpp::CharacterVector>(paramsPhenology["PhenologyType"]);
  NumericVector leafDuration = Rcpp::as<Rcpp::NumericVector>(paramsPhenology["LeafDuration"]);
  
  //Transpiration parameters
  DataFrame paramsTransp = Rcpp::as<Rcpp::DataFrame>(x["paramsTranspiration"]);
  NumericVector Psi_Extract = Rcpp::as<Rcpp::NumericVector>(paramsTransp["Psi_Extract"]);
  

  //Water storage parameters
  DataFrame paramsWaterStorage = Rcpp::as<Rcpp::DataFrame>(x["paramsWaterStorage"]);
  NumericVector StemPI0 = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["StemPI0"]);
  NumericVector StemEPS = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["StemEPS"]);
  NumericVector StemAF = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["StemAF"]);
  NumericVector Vsapwood = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["Vsapwood"]); //l·m-2 = mm
  NumericVector LeafPI0 = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["LeafPI0"]);
  NumericVector LeafEPS = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["LeafEPS"]);
  NumericVector LeafAF = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["LeafAF"]);
  NumericVector Vleaf = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["Vleaf"]); //l·m-2 = mm
  

  //Ring of forming vessels
  List ringList = as<Rcpp::List>(x["internalRings"]);
  
  
  //Daily output vectors
  NumericVector LabileCarbonBalance(numCohorts,0.0);
  NumericVector MaintenanceRespiration(numCohorts,0.0);
  NumericVector GrowthCosts(numCohorts,0.0);
  NumericVector PlantSugarTransport(numCohorts,0.0), PlantSugarLeaf(numCohorts,0.0), PlantStarchLeaf(numCohorts,0.0);
  NumericVector PlantSugarSapwood(numCohorts,0.0), PlantStarchSapwood(numCohorts,0.0);
  NumericVector SapwoodArea(numCohorts,0.0), LeafArea(numCohorts,0.0);
  NumericVector SAgrowth(numCohorts,0.0);
  NumericVector LAgrowth(numCohorts,0.0);
  NumericVector GrossPhotosynthesis(numCohorts,0.0);
  NumericVector RootExudation(numCohorts,0.0);

  NumericVector deltaLAgrowth(numCohorts,0.0);
  NumericVector deltaSAgrowth(numCohorts,0.0);
  
  double equilibriumLeafSugarConc = equilibriumLeafTotalConc - nonSugarConcentration;
  double equilibriumSapwoodSugarConc = equilibriumSapwoodTotalConc - nonSugarConcentration;
  double minimumSugarForGrowth = equilibriumSapwoodSugarConc*minimumRelativeSugarForGrowth;
  double mortalitySugarThreshold = equilibriumSapwoodSugarConc*mortalityRelativeSugarThreshold;
    
  double rleafcellmax = relative_expansion_rate(0.0 ,30.0, -2.0,0.5,0.05,5.0);
  
  //Initial Biomass balance
  NumericVector LeafBiomassBalance(numCohorts,0.0), FineRootBiomassBalance(numCohorts,0.0);
  DataFrame ccIni = carbonCompartments(x, "g_ind");
  DataFrame plantBiomassBalance = initPlantBiomassBalance(ccIni, above);
  NumericVector Volume_leaves = Rcpp::as<Rcpp::NumericVector>(ccIni["LeafStorageVolume"]);
  NumericVector Volume_sapwood = Rcpp::as<Rcpp::NumericVector>(ccIni["SapwoodStorageVolume"]);
  NumericVector Starch_max_leaves = Rcpp::as<Rcpp::NumericVector>(ccIni["LeafStarchMaximumConcentration"]);
  NumericVector Starch_max_sapwood = Rcpp::as<Rcpp::NumericVector>(ccIni["SapwoodStarchMaximumConcentration"]);
  NumericVector LeafStructBiomass = Rcpp::as<Rcpp::NumericVector>(ccIni["LeafStructuralBiomass"]);
  NumericVector SapwoodLivingStructBiomass = Rcpp::as<Rcpp::NumericVector>(ccIni["SapwoodLivingStructuralBiomass"]);
  NumericVector TotalLivingBiomass = Rcpp::as<Rcpp::NumericVector>(ccIni["TotalLivingBiomass"]);
  NumericVector TotalBiomass = Rcpp::as<Rcpp::NumericVector>(ccIni["TotalBiomass"]);

  //3. Carbon balance and growth
  for(int j=0;j<numCohorts;j++){
    if(N[j] > 0.0) {
      double costPerLA = 1000.0*leaf_CC/SLA[j]; // Construction cost in g gluc · m-2 of leaf area
      double costPerSA = sapwood_CC*sapwoodStructuralBiomass(1.0, H[j], L(j,_),V(j,_),WoodDensity[j]); // Construction cost in g gluc · cm-2 of sapwood area
      
      double LAexpanded = leafArea(LAI_expanded[j], N[j]);
      double LAlive = leafArea(LAI_live[j], N[j]);
      double LAdead = leafArea(LAI_dead[j], N[j]);

      double leafRespDay = 0.0;

      //MAINTENANCE RESPIRATION
      //Respiratory biomass (g dw · ind-1)
      double leafSugarMass = sugarLeaf[j]*(Volume_leaves[j]*glucoseMolarMass);
      double sapwoodSugarMass = sugarSapwood[j]*(Volume_sapwood[j]*glucoseMolarMass);
      double B_resp_leaves = LeafStructBiomass[j] + leafSugarMass;
      double B_resp_sapwood = SapwoodLivingStructBiomass[j] + sapwoodSugarMass;
      // Rcout<<j<< " maintenance costs of leaf sugars: "<< (leafSugarMass/B_resp_leaves)<<" sapwood sugars: "<< (sapwoodSugarMass/B_resp_sapwood)<<"\n";
      double B_resp_fineroots = fineRootBiomass[j];
      double QR = qResp(tday);
      if(LAexpanded>0.0) leafRespDay = B_resp_leaves*RERleaf[j]*QR;
      double sapwoodResp = B_resp_sapwood*RERsapwood[j]*QR;
      double finerootResp = B_resp_fineroots*RERfineroot[j]*QR;
      MaintenanceRespiration[j] += (leafRespDay+sapwoodResp+finerootResp)/TotalLivingBiomass[j]; 
      
      //PHOTOSYNTHESIS
      double leafAgG = 0.0;
      if(LAexpanded>0.0) {
        //gross fotosynthesis
        double leafAgC = Ag[j]/(N[j]/10000.0); //Translate g C · m-2 · h-1 to g C · h-1
        leafAgG = leafAgC*(glucoseMolarMass/(carbonMolarMass*6.0)); // from g C· h-1 to g gluc · h-1
        
        //Update output values
        GrossPhotosynthesis[j] = leafAgG/TotalLivingBiomass[j]; //Ag in g gluc · gdry-1 
      }
      
      //GROWTH
      double growthCostLA = 0.0;
      double growthCostSA = 0.0;
      List ring = ringList[j];
      grow_ring(ring, PlantPsi[j] ,tday, 10.0);
      double rleafcell = std::min(rleafcellmax, relative_expansion_rate(PlantPsi[j] ,tday, LeafPI0[j],0.5,0.05,5.0));
      
      if(leafUnfolding[j]) {
        double deltaLApheno = std::max(leafAreaTarget[j] - LAexpanded, 0.0);
        double deltaLAsink = std::min(deltaLApheno, SA[j]*RGRleafmax[j]*(rleafcell/rleafcellmax));
        if(!sinkLimitation) deltaLAsink = std::min(deltaLApheno, SA[j]*RGRleafmax[j]); //Deactivates temperature and turgor limitation
        double deltaLAavailable = 0.0;
        deltaLAavailable = std::max(0.0,((sugarSapwood[j] - minimumSugarForGrowth)*(glucoseMolarMass*Volume_sapwood[j]))/costPerLA);
        deltaLAgrowth[j] = std::min(deltaLAsink, deltaLAavailable);
        growthCostLA = deltaLAgrowth[j]*costPerLA;
      }
      double leafBiomassIncrement = deltaLAgrowth[j]*(1000.0/SLA[j]);
      
      //Assumes increment biomass of fine roots is half leaf biomass increment
      double finerootBiomassIncrement = leafBiomassIncrement/2.0;
      double growthCostFR = (growthCostLA/2.0)*(fineroot_CC/leaf_CC);   
        
      if(LAexpanded>0.0) {
        NumericVector SAring = ring["SA"];
        double deltaSAring = 0.0;
        if(SAring.size()==1) deltaSAring = SAring[0];
        else deltaSAring = SAring[SAring.size()-1] - SAring[SAring.size()-2];
        double cellfileareamaxincrease = 950.0; //Found empirically with T = 30 degrees and Psi = -0.033
        double rgrcellfile = (deltaSAring/10.0)/cellfileareamaxincrease;
        double deltaSAsink = (SA[j]*RGRsapwoodmax[j]*rgrcellfile); 
        if(!sinkLimitation) deltaSAsink = SA[j]*RGRsapwoodmax[j]; //Deactivates temperature and turgor limitation
        double deltaSAavailable = 0.0;
        if(sugarSapwood[j] > minimumSugarForGrowth) {
          deltaSAavailable = (starchSapwood[j]*(glucoseMolarMass*Volume_sapwood[j]))/costPerSA;
        }
        deltaSAgrowth[j] = std::min(deltaSAsink, deltaSAavailable);
        // Rcout<< SAring.size()<<" " <<j<< " "<< PlantPsi[j]<< " "<< LeafPI0[j]<<" dSAring "<<deltaSAring<< " dSAsink "<< deltaSAsink<<" dSAgrowth "<< deltaSAgrowth<<" rgrcellfile"<< rgrcellfile<<"\n";
        growthCostSA = deltaSAgrowth[j]*costPerSA; //increase cost (may be non zero if leaf growth was charged onto sapwood)
      }
      
      GrowthCosts[j] +=(growthCostLA + growthCostSA + growthCostFR)/TotalLivingBiomass[j]; //growth cost in g gluc · gdry-1
      
      
      //PARTIAL CARBON BALANCE
      double leafSugarMassDelta = leafAgG - leafRespDay;
      double sapwoodSugarMassDelta =  - sapwoodResp - finerootResp; 
      double sapwoodStarchMassDelta =  - growthCostFR - growthCostLA - growthCostSA;
      
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


      //SENESCENCE
      //Leaf senescence
      double propLeafSenescence = 0.0;
      //Leaf senescence due to age (Ca+ accumulation) only in evergreen species
      if(phenoType[j] == "progressive-evergreen") {
        propLeafSenescence = (1.0/(365.25*leafDuration[j]));
      }
      else if((phenoType[j] == "oneflush-evergreen") & (leafSenescence[j])) {
        propLeafSenescence = (1.0/leafDuration[j]); // Fraction of old leaves that die
        leafSenescence[j] = false; //To prevent further loss
      }
      else if(((phenoType[j] == "winter-deciduous") || (phenoType[j] == "winter-semideciduous")) & leafSenescence[j]) {
        propLeafSenescence = 1.0;
        leafSenescence[j] = false; //To prevent further loss
      }
      //Leaf senescence due to drought 
      double LAplc = std::min(LAexpanded, (1.0 - StemPLC[j])*leafAreaTarget[j]);
      if(LAplc<LAexpanded) {
        propLeafSenescence = std::max((LAexpanded-LAplc)/LAexpanded, propLeafSenescence); 
      }
      double deltaLAsenescence = std::min(LAexpanded, LAexpanded*propLeafSenescence);
      double senescenceLeafLoss = deltaLAsenescence*(1000.0/SLA[j]);
      double senescenceFinerootLoss = (senescenceLeafLoss/2.0);
      
      //Define sapwood senescense
      double propSASenescence = SRsapwood[j]/(1.0+15.0*exp(-0.01*H[j]));
      double deltaSASenescence = propSASenescence*SA[j];
      

      //TRANSLOCATION (in mol gluc) of labile carbon
      double translocationSugarLeaf = propLeafSenescence*Volume_leaves[j]*sugarLeaf[j];
      double translocationStarchLeaf = propLeafSenescence*Volume_leaves[j]*starchLeaf[j];
      double translocationSugarSapwood = propSASenescence*Volume_sapwood[j]*starchSapwood[j];
      if(Volume_leaves[j]>0) {
        sugarLeaf[j] = ((sugarLeaf[j]*Volume_leaves[j]) - translocationSugarLeaf)/Volume_leaves[j]; 
        starchLeaf[j] = ((starchLeaf[j]*Volume_leaves[j]) - translocationStarchLeaf)/Volume_leaves[j]; 
      }
      sugarSapwood[j] = ((sugarSapwood[j]*Volume_sapwood[j]) - translocationSugarSapwood)/Volume_sapwood[j]; 
      starchSapwood[j] = ((starchSapwood[j]*Volume_sapwood[j]) + translocationSugarLeaf + translocationStarchLeaf + translocationSugarSapwood)/Volume_sapwood[j]; 
      
      //ROOT EXUDATION
      //Excess starch carbon is lost as root exudation
      if(starchLeaf[j] > Starch_max_leaves[j]) {
        RootExudation[j] += ((starchLeaf[j] - Starch_max_leaves[j])*(Volume_leaves[j]*glucoseMolarMass)/TotalLivingBiomass[j]);
        starchLeaf[j] = Starch_max_leaves[j];
      }
      if(starchSapwood[j] > Starch_max_sapwood[j]) {
        RootExudation[j] += ((starchSapwood[j] - Starch_max_sapwood[j])*(Volume_sapwood[j]*glucoseMolarMass)/TotalLivingBiomass[j]);
        starchSapwood[j] = Starch_max_sapwood[j];
      }
      
      //Labile CARBON balance
      LabileCarbonBalance[j] = GrossPhotosynthesis[j] - MaintenanceRespiration[j] - GrowthCosts[j] - RootExudation[j];
        
      //UPDATE LEAF AREA, SAPWOOD AREA, FINE ROOT BIOMASS AND CONCENTRATION IN LABILE POOLS
      double LAprev = LAexpanded;
      LAexpanded += deltaLAgrowth[j] - deltaLAsenescence;
      LAdead += deltaLAsenescence;
      LAI_dead[j] = LAdead*N[j]/10000.0;
      LAI_expanded[j] = LAexpanded*N[j]/10000.0;
      double SAprev = SA[j];
      SA[j] = SA[j] + deltaSAgrowth[j] - deltaSASenescence; 
      fineRootBiomass[j] = fineRootBiomass[j] + finerootBiomassIncrement - senescenceFinerootLoss;
      
      //UPDATE DERIVED QUANTITIES (individual level)      
      //Update Huber value
      if(LAlive>0.0) {
        Al2As[j] = (LAlive)/(SA[j]/10000.0);
      }
      //Decrease PLC due to new SA growth
      if(cavitationRefill=="growth") StemPLC[j] = std::max(0.0, StemPLC[j] - (deltaSAgrowth[j]/SA[j]));
      
      
      
      //LEAF/FINE ROOT BIOMASS balance (g_ind)
      LeafBiomassBalance[j] = leafBiomassIncrement - senescenceLeafLoss;
      FineRootBiomassBalance[j] = finerootBiomassIncrement - senescenceFinerootLoss;
      
      //MORTALITY Death by carbon starvation or dessication
      double Nprev = N[j]; //Store initial density (for biomass balance)
      double Ndead_day = 0.0;
      bool dynamicCohort = true;
      bool isShrub = !NumericVector::is_na(Cover[j]);
      if((!shrubDynamics) & isShrub) dynamicCohort = false;
      double stemSympRWC = symplasticRelativeWaterContent(PlantPsi[j], StemPI0[j], StemEPS[j]);
      if(dynamicCohort) {
        if(mortalityMode=="whole-cohort/deterministic") {
          if((sugarSapwood[j]<mortalitySugarThreshold) & allowStarvation) {
            Ndead_day = N[j];
            if(verbose) Rcout<<" [Cohort "<< j<<" died from starvation] ";
          } else if( (stemSympRWC < mortalityRWCThreshold) & allowDessication) {
            Ndead_day = N[j];
            if(verbose) Rcout<<" [Cohort "<< j<<" died from dessication] ";
          }
        } else {
          double P_starv = dailyMortalityProbability(mortalityBaselineRate, sugarSapwood[j], 
                                                     mortalitySugarThreshold, allowStarvation,
                                                     0.0, 2.0);
          double P_dess = dailyMortalityProbability(mortalityBaselineRate, stemSympRWC, 
                                                    mortalityRWCThreshold, allowDessication,
                                                    0.0, 2.0);
          double P_day = std::max(P_dess, P_starv);
          if(mortalityMode =="density/deterministic") {
            Ndead_day = N[j]*P_day;
          } else if(mortalityMode =="whole-cohort/stochastic") {
            if(R::runif(0.0,1.0) < P_day) {
              Ndead_day = N[j];
              if(P_dess>P_starv) {
                Rcout<<" [Cohort "<< j<<" died from dessication] ";
              } else {
                Rcout<<" [Cohort "<< j<<" died from starvation] ";
              }
            }
          } else if(mortalityMode == "density/stochastic") {
            Ndead_day = R::rbinom(round(N[j]), P_day);
            // Rcout<< j<< " "<< P_day<< " "<< N[j]<< " "<< Ndead_day<< " "<< R::rbinom(750, 2.82309e-05)<< "\n";
          }
        }
        // Update density and increase the number of dead plants
        Ndead_day = std::min(Ndead_day, N[j]);
        double Cdead_day = Cover[j]*(Ndead_day/N[j]);
        N[j] = N[j] - Ndead_day;
        N_dead[j] = N_dead[j] + Ndead_day;
        if(isShrub) {
          Cover[j] = std::max(0.0, Cover[j] - Cdead_day);
          Cover_dead[j] = Cover_dead[j] + Cdead_day;
        }
        //Update LAI dead and LAI expanded as a result of density decrease
        double LAI_change = LAexpanded*Ndead_day/10000.0;
        LAI_dead[j] = LAI_dead[j] + LAI_change;
        LAI_expanded[j] = LAI_expanded[j] - LAI_change;
      }
      

      //UPDATE TARGETS
      //Set target leaf area if bud formation is allowed
      if(budFormation[j]) {
        leafAreaTarget[j] = (SA[j]/10000.0)*allocationTarget[j];
        LAI_live[j] = leafAreaTarget[j]*N[j]/10000.0;
      }
      
      //Output variables
      SapwoodArea[j] = SA[j];
      SAgrowth[j] += deltaSAgrowth[j]/SAprev; //Store sapwood area growth rate (cm2/cm2)
      LAgrowth[j] += deltaLAgrowth[j]/SAprev;//Store Leaf area growth rate in relation to sapwood area (m2/cm2)
      LeafArea[j] = LAexpanded;
    }
  }
  
  //UPDATE STRUCTURAL VARIABLES
  updateStructuralVariables(x, deltaSAgrowth);
  
  //RECALCULATE storage concentrations (SA, LA and H may have changed)
  for(int j=0;j<numCohorts;j++){
    double newVolumeSapwood = sapwoodStorageVolume(SA[j], H[j], L(j,_),V(j,_),WoodDensity[j], conduit2sapwood[j]);
    double newVolumeLeaves = leafStorageVolume(LAI_expanded[j],  N[j], SLA[j], LeafDensity[j]);
    if(newVolumeLeaves > 0.0) {
      sugarLeaf[j] = sugarLeaf[j]*(Volume_leaves[j]/newVolumeLeaves);
      starchLeaf[j] = starchLeaf[j]*(Volume_leaves[j]/newVolumeLeaves); 
    } else {
      sugarLeaf[j] = 0.0;
      starchLeaf[j] = 0.0;
    }
    sugarSapwood[j] = sugarSapwood[j]*(Volume_sapwood[j]/newVolumeSapwood);
    starchSapwood[j] = starchSapwood[j]*(Volume_sapwood[j]/newVolumeSapwood); 
    
    //OUTPUT VARIABLES
    PlantSugarLeaf[j] = sugarLeaf[j];
    PlantStarchLeaf[j] = starchLeaf[j];
    PlantSugarSapwood[j] = sugarSapwood[j];
    PlantStarchSapwood[j] = starchSapwood[j];
  }
  
  //CLOSE BIOMASS BALANCE
  closePlantBiomassBalance(plantBiomassBalance, x,
                      LabileCarbonBalance, LeafBiomassBalance, FineRootBiomassBalance);
  
  //Update pool proportions??
  if(plantWaterPools) {
    NumericVector poolProportions = Rcpp::as<Rcpp::NumericVector>(belowdf["poolProportions"]);
  }
  

  DataFrame labileCarbonBalance = DataFrame::create(_["GrossPhotosynthesis"] = GrossPhotosynthesis,
                                                   _["MaintenanceRespiration"] = MaintenanceRespiration,
                                                   _["GrowthCosts"] = GrowthCosts,
                                                   _["RootExudation"] = RootExudation,
                                                   _["LabileCarbonBalance"] = LabileCarbonBalance,
                                                   _["SugarLeaf"] = PlantSugarLeaf,
                                                   _["StarchLeaf"] = PlantStarchLeaf,
                                                   _["SugarSapwood"] = PlantSugarSapwood,
                                                   _["StarchSapwood"] = PlantStarchSapwood,
                                                   _["SugarTransport"] = PlantSugarTransport,
                                                   _["StemPI0"] = clone(StemPI0), //Store a copy of the current osmotic potential at full turgor
                                                   _["LeafPI0"] = clone(LeafPI0));
  labileCarbonBalance.attr("row.names") = above.attr("row.names");
  

  //Final Biomass compartments
  DataFrame plantStructure = DataFrame::create(
    _["LeafArea"] = LeafArea,
    _["SapwoodArea"] = SapwoodArea,
    _["FineRootBiomass"] = clone(fineRootBiomass),
    _["DBH"] = clone(DBH),
    _["Height"] = clone(H)
  );
  
  DataFrame plantGrowth = DataFrame::create(
    _["SAgrowth"] = SAgrowth,
    _["LAgrowth"] = LAgrowth
  );
  plantGrowth.attr("row.names") = above.attr("row.names");
  
  List l = List::create(_["cohorts"] = clone(cohorts),
                        _["WaterBalance"] = spwbOut["WaterBalance"], 
                        _["Soil"] = spwbOut["Soil"], 
                        _["Stand"] = spwbOut["Stand"], 
                        _["Plants"] = spwbOut["Plants"],
                        _["LabileCarbonBalance"] = labileCarbonBalance,
                        _["PlantBiomassBalance"] = plantBiomassBalance,
                        _["PlantStructure"] = plantStructure,
                        _["PlantGrowth"] = plantGrowth);
  l.attr("class") = CharacterVector::create("growth_day","list");
  return(l);
}




List growthDay2(List x, NumericVector meteovec, 
                double latitude, double elevation, double slope, double aspect,
                double solarConstant, double delta, 
                double runon=0.0, bool verbose = false) {
  
  //1. Soil-plant water balance
  List spwbOut = spwbDay2(x, meteovec, 
                          latitude, elevation, slope, aspect,
                          solarConstant, delta, 
                          runon, verbose);
  

  //2. Retrieve state
  
  //Control params
  List control = x["control"];  

  //Meteo input
  double tmin = meteovec["tmin"];
  double tmax = meteovec["tmax"];
  double tminPrev = meteovec["tminPrev"];
  double tmaxPrev = meteovec["tmaxPrev"];
  double tminNext = meteovec["tminNext"];
  double rhmin = meteovec["rhmin"];
  double rhmax = meteovec["rhmax"];
  double rad = meteovec["rad"];
  double wind = meteovec["wind"];
  double Catm = meteovec["Catm"];
  
  double tday = meteoland::utils_averageDaylightTemperature(tmin, tmax);
  
  String mortalityMode = control["mortalityMode"];
  double mortalityBaselineRate = control["mortalityBaselineRate"];
  double mortalityRelativeSugarThreshold= control["mortalityRelativeSugarThreshold"];
  double mortalityRWCThreshold= control["mortalityRWCThreshold"];
  
  bool allowDessication = control["allowDessication"];
  bool allowStarvation = control["allowStarvation"];
  bool allowDefoliation = control["allowDefoliation"];
  bool sinkLimitation = control["sinkLimitation"];
  bool shrubDynamics = control["shrubDynamics"];
  String allocationStrategy = control["allocationStrategy"];
  String cavitationRefill = control["cavitationRefill"];
  bool plantWaterPools = control["plantWaterPools"];
  bool taper = control["taper"];
  bool nonStomatalPhotosynthesisLimitation = control["nonStomatalPhotosynthesisLimitation"];
  double averageFracRhizosphereResistance = control["averageFracRhizosphereResistance"];
  double phloemConductanceFactor = control["phloemConductanceFactor"];
  double nonSugarConcentration = control["nonSugarConcentration"];
  List equilibriumOsmoticConcentration  = control["equilibriumOsmoticConcentration"];
  double equilibriumLeafTotalConc = equilibriumOsmoticConcentration["leaf"];
  double equilibriumSapwoodTotalConc = equilibriumOsmoticConcentration["sapwood"];
  double minimumRelativeSugarForGrowth = control["minimumRelativeSugarForGrowth"];
  List constructionCosts = control["constructionCosts"];
  double leaf_CC = constructionCosts["leaf"];
  double sapwood_CC = constructionCosts["sapwood"];
  double fineroot_CC = constructionCosts["fineroot"];

  //Soil params
  List soil  = x["soil"];
  NumericVector Ksat = soil["Ksat"];
  NumericVector dVec = soil["dVec"];
  NumericVector rfc = soil["rfc"];
  NumericVector VG_n = soil["VG_n"];
  NumericVector VG_alpha = soil["VG_alpha"];
  NumericVector Tsoil = soil["Temp"];
  
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
  NumericVector SA = above["SA"];

  
  
  //Belowground parameters  
  DataFrame belowdf = Rcpp::as<Rcpp::DataFrame>(x["below"]);
  NumericVector Z95 = Rcpp::as<Rcpp::NumericVector>(belowdf["Z95"]);
  NumericVector Z50 = Rcpp::as<Rcpp::NumericVector>(belowdf["Z50"]);
  NumericVector fineRootBiomass = Rcpp::as<Rcpp::NumericVector>(belowdf["fineRootBiomass"]);
  NumericVector CRSV = Rcpp::as<Rcpp::NumericVector>(belowdf["coarseRootSoilVolume"]);
  List belowLayers = Rcpp::as<Rcpp::List>(x["belowLayers"]);
  NumericMatrix V = Rcpp::as<Rcpp::NumericMatrix>(belowLayers["V"]);
  NumericMatrix L = Rcpp::as<Rcpp::NumericMatrix>(belowLayers["L"]);
  NumericMatrix RhizoPsi = Rcpp::as<Rcpp::NumericMatrix>(belowLayers["RhizoPsi"]);
  NumericMatrix VCroot_kmax = Rcpp::as<Rcpp::NumericMatrix>(belowLayers["VCroot_kmax"]);
  NumericMatrix VGrhizo_kmax = Rcpp::as<Rcpp::NumericMatrix>(belowLayers["VGrhizo_kmax"]);
  int numLayers = VCroot_kmax.ncol();
  
  //Internal state variables
  DataFrame internalWater = Rcpp::as<Rcpp::DataFrame>(x["internalWater"]);
  NumericVector NSPL = Rcpp::as<Rcpp::NumericVector>(internalWater["NSPL"]);

  //Values at the end of the day (after calling spwb)
  NumericVector psiApoLeaf = Rcpp::as<Rcpp::NumericVector>(internalWater["LeafPsi"]);
  NumericVector psiApoStem = Rcpp::as<Rcpp::NumericVector>(internalWater["Stem1Psi"]);
  NumericVector psiSympLeaf = Rcpp::as<Rcpp::NumericVector>(internalWater["LeafSympPsi"]);
  NumericVector psiSympStem = Rcpp::as<Rcpp::NumericVector>(internalWater["StemSympPsi"]);
  NumericVector StemPLC = Rcpp::as<Rcpp::NumericVector>(internalWater["StemPLC"]);
  
  DataFrame internalCarbon = Rcpp::as<Rcpp::DataFrame>(x["internalCarbon"]);
  NumericVector sugarLeaf = internalCarbon["sugarLeaf"]; //Concentrations assuming RWC = 1
  NumericVector starchLeaf = internalCarbon["starchLeaf"];
  NumericVector sugarSapwood = internalCarbon["sugarSapwood"];
  NumericVector starchSapwood = internalCarbon["starchSapwood"];
  
  DataFrame internalMortality = Rcpp::as<Rcpp::DataFrame>(x["internalMortality"]);
  NumericVector N_dead = internalMortality["N_dead"];
  NumericVector Cover_dead = internalMortality["Cover_dead"];
  
  DataFrame internalAllocation = Rcpp::as<Rcpp::DataFrame>(x["internalAllocation"]);
  NumericVector allocationTarget = internalAllocation["allocationTarget"];
  NumericVector leafAreaTarget = internalAllocation["leafAreaTarget"];
  NumericVector fineRootBiomassTarget = internalAllocation["fineRootBiomassTarget"];

  DataFrame internalPhenology = Rcpp::as<Rcpp::DataFrame>(x["internalPhenology"]);
  LogicalVector leafUnfolding = internalPhenology["leafUnfolding"];
  LogicalVector budFormation = internalPhenology["budFormation"];
  LogicalVector leafSenescence = internalPhenology["leafSenescence"];
  
  List stand = spwbOut["Stand"];
  DataFrame Plants = Rcpp::as<Rcpp::DataFrame>(spwbOut["Plants"]);
  List PlantsInst = spwbOut["PlantsInst"];
  
  //Recover module-communication state variables
  NumericMatrix AgStep  =  Rcpp::as<Rcpp::NumericMatrix>(PlantsInst["Ag"]);
  int numSteps = AgStep.ncol();
  
  //Data from spwb
  NumericVector LeafSympRWC = Plants["LeafSympRWC"];
  NumericVector StemSympRWC = Plants["StemSympRWC"];
  NumericMatrix StemSympPsiInst =  Rcpp::as<Rcpp::NumericMatrix>(PlantsInst["StemSympPsi"]);
  NumericMatrix LeafSympPsiInst =  Rcpp::as<Rcpp::NumericMatrix>(PlantsInst["LeafSympPsi"]);
  NumericMatrix StemSympRWCInst =  Rcpp::as<Rcpp::NumericMatrix>(PlantsInst["StemSympRWC"]);
  NumericMatrix LeafSympRWCInst =  Rcpp::as<Rcpp::NumericMatrix>(PlantsInst["LeafSympRWC"]);
  
  List eb = spwbOut["EnergyBalance"];  
  DataFrame tempDF =  Rcpp::as<Rcpp::DataFrame>(eb["Temperature"]);
  NumericVector Tcan = Rcpp::as<Rcpp::NumericVector>(tempDF["Tcan"]);
  
  //Anatomy parameters
  DataFrame paramsAnatomy = Rcpp::as<Rcpp::DataFrame>(x["paramsAnatomy"]);
  NumericVector Hmed = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["Hmed"]);
  NumericVector SLA = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["SLA"]);
  NumericVector Al2As = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["Al2As"]);
  NumericVector WoodDensity = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["WoodDensity"]);
  NumericVector LeafDensity = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["LeafDensity"]);
  NumericVector FineRootDensity = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["FineRootDensity"]);
  NumericVector SRL = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["SRL"]);
  NumericVector RLD = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["RLD"]);
  NumericVector conduit2sapwood = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["conduit2sapwood"]);
  
  //Growth parameters
  DataFrame paramsGrowth = Rcpp::as<Rcpp::DataFrame>(x["paramsGrowth"]);
  NumericVector WoodC = Rcpp::as<Rcpp::NumericVector>(paramsGrowth["WoodC"]);
  NumericVector RERleaf = Rcpp::as<Rcpp::NumericVector>(paramsGrowth["RERleaf"]);
  NumericVector RERsapwood = Rcpp::as<Rcpp::NumericVector>(paramsGrowth["RERsapwood"]);
  NumericVector RERfineroot = Rcpp::as<Rcpp::NumericVector>(paramsGrowth["RERfineroot"]);
  NumericVector RGRleafmax = Rcpp::as<Rcpp::NumericVector>(paramsGrowth["RGRleafmax"]);
  NumericVector RGRsapwoodmax = Rcpp::as<Rcpp::NumericVector>(paramsGrowth["RGRsapwoodmax"]);
  NumericVector RGRfinerootmax = Rcpp::as<Rcpp::NumericVector>(paramsGrowth["RGRfinerootmax"]);
  NumericVector SRsapwood = Rcpp::as<Rcpp::NumericVector>(paramsGrowth["SRsapwood"]);
  NumericVector SRfineroot = Rcpp::as<Rcpp::NumericVector>(paramsGrowth["SRfineroot"]);
  
  
  //Phenology parameters
  DataFrame paramsPhenology = Rcpp::as<Rcpp::DataFrame>(x["paramsPhenology"]);
  CharacterVector phenoType = Rcpp::as<Rcpp::CharacterVector>(paramsPhenology["PhenologyType"]);
  NumericVector leafDuration = Rcpp::as<Rcpp::NumericVector>(paramsPhenology["LeafDuration"]);
  
  // NumericVector Cstoragepmax= Rcpp::as<Rcpp::NumericVector>(paramsGrowth["Cstoragepmax"]);
  // NumericVector slowCstorage_max(numCohorts), fastCstorage_max(numCohorts);
  //Transpiration parameters
  DataFrame paramsTransp = Rcpp::as<Rcpp::DataFrame>(x["paramsTranspiration"]);
  NumericVector Kmax_stemxylem = paramsTransp["Kmax_stemxylem"];
  NumericVector Plant_kmax= paramsTransp["Plant_kmax"];
  NumericVector VCleaf_kmax = Rcpp::as<Rcpp::NumericVector>(paramsTransp["VCleaf_kmax"]);
  NumericVector VCleaf_c = Rcpp::as<Rcpp::NumericVector>(paramsTransp["VCleaf_c"]);
  NumericVector VCleaf_d = Rcpp::as<Rcpp::NumericVector>(paramsTransp["VCleaf_d"]);
  NumericVector VCstem_kmax = paramsTransp["VCstem_kmax"];
  NumericVector VCstem_c = Rcpp::as<Rcpp::NumericVector>(paramsTransp["VCstem_c"]);
  NumericVector VCstem_d = Rcpp::as<Rcpp::NumericVector>(paramsTransp["VCstem_d"]);
  NumericVector VCroot_kmaxVEC= paramsTransp["VCroot_kmax"];
  NumericVector VCroot_c = paramsTransp["VCroot_c"];
  NumericVector VCroot_d = paramsTransp["VCroot_d"];
  NumericVector VGrhizo_kmaxVEC= paramsTransp["VGrhizo_kmax"];
  
  //Water storage parameters
  DataFrame paramsWaterStorage = Rcpp::as<Rcpp::DataFrame>(x["paramsWaterStorage"]);
  NumericVector StemPI0 = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["StemPI0"]);
  NumericVector StemEPS = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["StemEPS"]);
  NumericVector StemAF = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["StemAF"]);
  NumericVector Vsapwood = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["Vsapwood"]); //l·m-2 = mm
  NumericVector LeafPI0 = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["LeafPI0"]);
  NumericVector LeafEPS = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["LeafEPS"]);
  NumericVector LeafAF = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["LeafAF"]);
  NumericVector Vleaf = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["Vleaf"]); //l·m-2 = mm
  
  
  //Ring of forming vessels
  List ringList = as<Rcpp::List>(x["internalRings"]);
  
  //Subdaily output matrices
  NumericMatrix LabileCarbonBalanceInst(numCohorts, numSteps);  
  NumericMatrix GrossPhotosynthesisInst(numCohorts, numSteps);  
  NumericMatrix MaintenanceRespirationInst(numCohorts, numSteps);  
  NumericMatrix GrowthCostsInst(numCohorts, numSteps);  
  NumericMatrix RootExudationInst(numCohorts, numSteps);  
  NumericMatrix PlantSugarTransportInst(numCohorts, numSteps);
  NumericMatrix PlantSugarLeafInst(numCohorts, numSteps), PlantStarchLeafInst(numCohorts, numSteps);
  NumericMatrix PlantSugarSapwoodInst(numCohorts, numSteps), PlantStarchSapwoodInst(numCohorts, numSteps);
  
  //Daily output vectors
  NumericVector LabileCarbonBalance(numCohorts,0.0);
  NumericVector MaintenanceRespiration(numCohorts,0.0);
  NumericVector GrowthCosts(numCohorts,0.0);
  NumericVector PlantSugarTransport(numCohorts,0.0), PlantSugarLeaf(numCohorts,0.0), PlantStarchLeaf(numCohorts,0.0);
  NumericVector PlantSugarSapwood(numCohorts,0.0), PlantStarchSapwood(numCohorts,0.0);
  NumericVector SapwoodArea(numCohorts,0.0), LeafArea(numCohorts,0.0), FineRootArea(numCohorts, 0.0);
  NumericVector SAgrowth(numCohorts,0.0), LAgrowth(numCohorts,0.0), FRAgrowth(numCohorts,0.0);
  NumericVector GrossPhotosynthesis(numCohorts,0.0);
  NumericVector RootExudation(numCohorts,0.0);

  NumericVector deltaLAgrowth(numCohorts,0.0);
  NumericVector deltaSAgrowth(numCohorts,0.0);
  
  double equilibriumLeafSugarConc = equilibriumLeafTotalConc - nonSugarConcentration;
  double equilibriumSapwoodSugarConc = equilibriumSapwoodTotalConc - nonSugarConcentration;
  double minimumSugarForGrowth = equilibriumSapwoodSugarConc*minimumRelativeSugarForGrowth;
  double mortalitySugarThreshold = equilibriumSapwoodSugarConc*mortalityRelativeSugarThreshold;
  
  double rleafcellmax = relative_expansion_rate(0.0,30.0, -2.0,0.5,0.05,5.0);

  //Initial Biomass balance
  NumericVector LeafBiomassBalance(numCohorts,0.0), FineRootBiomassBalance(numCohorts,0.0);
  DataFrame ccIni = carbonCompartments(x, "g_ind");
  DataFrame plantBiomassBalance = initPlantBiomassBalance(ccIni, above);
  NumericVector Volume_leaves = Rcpp::as<Rcpp::NumericVector>(ccIni["LeafStorageVolume"]);
  NumericVector Volume_sapwood = Rcpp::as<Rcpp::NumericVector>(ccIni["SapwoodStorageVolume"]);
  NumericVector Starch_max_leaves = Rcpp::as<Rcpp::NumericVector>(ccIni["LeafStarchMaximumConcentration"]);
  NumericVector Starch_max_sapwood = Rcpp::as<Rcpp::NumericVector>(ccIni["SapwoodStarchMaximumConcentration"]);
  NumericVector LeafStructBiomass = Rcpp::as<Rcpp::NumericVector>(ccIni["LeafStructuralBiomass"]);
  NumericVector SapwoodLivingStructBiomass = Rcpp::as<Rcpp::NumericVector>(ccIni["SapwoodLivingStructuralBiomass"]);
  NumericVector TotalLivingBiomass = Rcpp::as<Rcpp::NumericVector>(ccIni["TotalLivingBiomass"]);
  NumericVector TotalBiomass = Rcpp::as<Rcpp::NumericVector>(ccIni["TotalBiomass"]);
  
  //3. Carbon balance, growth and senescence by cohort
  for(int j=0;j<numCohorts;j++){
    if(N[j]>0.0){
      double LAexpanded = leafArea(LAI_expanded[j], N[j]);
      double LAlive = leafArea(LAI_live[j], N[j]);
      double LAdead = leafArea(LAI_dead[j], N[j]);

      
      double costPerLA = 1000.0*leaf_CC/SLA[j]; // Construction cost in g gluc · m-2 of leaf area
      double costPerSA = sapwood_CC*sapwoodStructuralBiomass(1.0, H[j], L(j,_),V(j,_),WoodDensity[j]); // Construction cost in g gluc · cm-2 of sapwood area
      NumericVector deltaFRBgrowth(numLayers, 0.0);
        

      double leafRespDay = 0.0;
      // double sfrRespDay = 0.0;
      
      //Estimate phloem conductance as a factor of stem conductance
      double k_phloem = VCstem_kmax[j]*phloemConductanceFactor*(0.018/1000.0);
        
      //3.0 Xylogenesis
      grow_ring(ringList[j], psiSympStem[j] ,tday, 10.0);
      double rleafcell = std::min(rleafcellmax, relative_expansion_rate(psiSympLeaf[j] ,tday, LeafPI0[j],0.5,0.05,5.0));
      NumericVector rfineroot(numLayers);
      for(int s=0;s<numLayers;s++) rfineroot[s] = relative_expansion_rate(RhizoPsi(j,s) ,tday, StemPI0[j],0.5,0.05,5.0);
      
      //3.1 Carbon balance and growth by steps
      for(int s=0;s<numSteps;s++) {
        
        // minimum concentration (mol gluc·l-1) to avoid turgor loss
        // double leafTLP = turgorLossPoint(LeafPI0[j], LeafEPS[j]);
        // double stemTLP = turgorLossPoint(StemPI0[j], StemEPS[j]);
        // double rwcLeafTLP = symplasticRelativeWaterContent(leafTLP, LeafPI0[j], LeafEPS[j]);
        // double rwcStemTLP = symplasticRelativeWaterContent(stemTLP, StemPI0[j], StemEPS[j]);
        // tlpConcLeaf = sugarConcentration(leafTLP,Tcan[s], nonSugarConc)*(rwcLeafTLP/rwcLeaf(j,s)); 
        // tlpConcSapwood = sugarConcentration(stemTLP,Tcan[s], nonSugarConc)*(rwcStemTLP/rwcStem(j,s)); 

        //Transform sugar concentration (mol gluc · l-1) to sugar mass (g gluc)
        // double lstvol = 0.001*(starchLeaf[j]/starchDensity);
        // double sstvol = 0.001*(starchSapwood[j]/starchDensity);
        double leafSugarMassStep = sugarLeaf[j]*(Volume_leaves[j]*glucoseMolarMass);
        double sapwoodSugarMassStep = sugarSapwood[j]*(Volume_sapwood[j]*glucoseMolarMass);
        
        //MAINTENANCE RESPIRATION
        double B_resp_leaves = LeafStructBiomass[j] + leafSugarMassStep;
        double B_resp_sapwood = SapwoodLivingStructBiomass[j] + sapwoodSugarMassStep;
        double B_resp_fineroots = fineRootBiomass[j];
        double QR = qResp(Tcan[s]);
        double leafRespStep = 0.0;
        if(LAexpanded>0.0) leafRespStep = B_resp_leaves*RERleaf[j]*QR/((double) numSteps);
        double sapwoodRespStep = B_resp_sapwood*RERsapwood[j]*QR/((double) numSteps);
        double finerootRespStep = B_resp_fineroots*RERfineroot[j]*QR/((double) numSteps);
        leafRespDay +=leafRespStep;
        MaintenanceRespirationInst(j,s) = (leafRespStep+sapwoodRespStep+finerootRespStep)/TotalLivingBiomass[j];//Rm in g gluc· gdry-1
        MaintenanceRespiration[j] += MaintenanceRespirationInst(j,s); 
        
        //PHOTOSYNTHESIS
        double leafAgStepG = 0.0;
        if(LAexpanded>0.0) {
          double leafAgStepC = AgStep(j,s)/(N[j]/10000.0); //Translate g C · m-2 · h-1 to g C · h-1
          leafAgStepG = leafAgStepC*(glucoseMolarMass/(carbonMolarMass*6.0)); // from g C· h-1 to g gluc · h-1
          GrossPhotosynthesisInst(j,s) = leafAgStepG/TotalLivingBiomass[j]; //Ag in g gluc · gdry-1
          GrossPhotosynthesis[j] += GrossPhotosynthesisInst(j,s); 
        }

        //GROWTH        
        double growthCostLAStep = 0.0;
        double growthCostSAStep = 0.0;
        double growthCostFRBStep = 0.0;   

        //Leaf growth
        if(leafUnfolding[j]) {
          double deltaLApheno = std::max(leafAreaTarget[j] - LAexpanded, 0.0);
          double deltaLAsink = std::min(deltaLApheno, (1.0/((double) numSteps))*SA[j]*RGRleafmax[j]*(rleafcell/rleafcellmax));
          if(!sinkLimitation) deltaLAsink = std::min(deltaLApheno, (1.0/((double) numSteps))*SA[j]*RGRleafmax[j]); //Deactivates temperature and turgor limitation
          //Grow at expense of stem sugar
          double deltaLAavailable = std::max(0.0,((sugarSapwood[j] - minimumSugarForGrowth)*(glucoseMolarMass*Volume_sapwood[j]))/costPerLA);
          double deltaLAgrowthStep = std::min(deltaLAsink, deltaLAavailable);
          growthCostLAStep += deltaLAgrowthStep*costPerLA;
          deltaLAgrowth[j] += deltaLAgrowthStep;
        }
        //fine root growth
        if(fineRootBiomass[j] < fineRootBiomassTarget[j]) {
          for(int s = 0;s<numLayers;s++) {
            double deltaFRBpheno = std::max(fineRootBiomassTarget[j] - fineRootBiomass[j], 0.0);
            double deltaFRBsink = (1.0/((double) numSteps))*(V(j,s)*fineRootBiomass[j])*RGRfinerootmax[j]*(rfineroot[s]/rleafcellmax);
            if(!sinkLimitation) deltaFRBsink = (1.0/((double) numSteps))*(V(j,s)*fineRootBiomass[j])*RGRfinerootmax[j]; //Deactivates temperature and turgor limitation
            double deltaFRBavailable = std::max(0.0,((sugarSapwood[j] - minimumSugarForGrowth)*(glucoseMolarMass*Volume_sapwood[j]))/fineroot_CC);
            double deltaFRBgrowthStep = std::min(deltaFRBpheno, std::min(deltaFRBsink, deltaFRBavailable));
            growthCostFRBStep += deltaFRBgrowthStep*fineroot_CC;
            deltaFRBgrowth[s] += deltaFRBgrowthStep;
          }
        }
        //sapwood area growth
        if(LAexpanded>0.0) {
          List ring = ringList[j];
          NumericVector SAring = ring["SA"];
          double deltaSAring = 0.0;
          if(SAring.size()==1) deltaSAring = SAring[0];
          else deltaSAring = SAring[SAring.size()-1] - SAring[SAring.size()-2];
          double cellfileareamaxincrease = 950.0; //Found empirically with T = 30 degrees and Psi = -0.033
          double rgrcellfile = (deltaSAring/10.0)/cellfileareamaxincrease;
          double deltaSAsink = (SA[j]*RGRsapwoodmax[j]*rgrcellfile)/((double) numSteps); 
          if(!sinkLimitation) deltaSAsink = SA[j]*RGRsapwoodmax[j]/((double) numSteps); //Deactivates temperature and turgor limitation
          double deltaSAavailable = 0.0;
          if(sugarSapwood[j] > minimumSugarForGrowth) {
            deltaSAavailable = std::max(0.0, starchSapwood[j]*(glucoseMolarMass*Volume_sapwood[j])/costPerSA);
          }
          double deltaSAgrowthStep = std::min(deltaSAsink, deltaSAavailable);
          if(deltaSAgrowthStep<0.0) {
            Rcout<<deltaSAsink<<" "<< deltaSAavailable<< " "<< starchSapwood[j]<<"\n";
            stop("negative growth!"); 
          }
          growthCostSAStep += deltaSAgrowthStep*costPerSA; //increase cost (may be non zero if leaf growth was charged onto sapwood)
          deltaSAgrowth[j]  +=deltaSAgrowthStep;
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
            // double ft = phloemFlow(LeafSympPsiInst(j,s), StemSympPsiInst(j,s), sugarLeaf[j]/LeafSympRWCInst(j,s), sugarSapwood[j]/StemSympRWCInst(j,s), Tcan[s], k_phloem, nonSugarConcentration)*LAlive; //flow as mol glucose per s
            // sugar-starch dynamics
            // double conversionLeaf = sugarStarchDynamicsLeaf(sugarLeaf[j]/LeafSympRWCInst(j,s), starchLeaf[j]/LeafSympRWCInst(j,s), equilibriumLeafSugarConc);
            double conversionLeaf = sugarStarchDynamicsLeaf(sugarLeaf[j], starchLeaf[j], equilibriumLeafSugarConc);
            // double starchLeafIncrease = conversionLeaf*LeafSympRWCInst(j,s);
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
        //Divert to root exudation if starch is over maximum capacity
        if(starchLeaf[j] > Starch_max_leaves[j]) {
          RootExudationInst(j,s) += ((starchLeaf[j] - Starch_max_leaves[j])*(Volume_leaves[j]*glucoseMolarMass)/TotalLivingBiomass[j]);
          starchLeaf[j] = Starch_max_leaves[j];
        }
        //Divert to root exudation if starch is over maximum capacity
        if(starchSapwood[j] > Starch_max_sapwood[j]) {
          RootExudationInst(j,s) += ((starchSapwood[j] - Starch_max_sapwood[j])*(Volume_sapwood[j]*glucoseMolarMass)/TotalLivingBiomass[j]);
          starchSapwood[j] = Starch_max_sapwood[j];
        }
        
        // if(starchSapwoodIncrease > Starch_max_sapwood[j] - starchSapwood[j]) {
        //   RootExudationInst(j,s) += ((starchSapwood[j] + starchSapwoodIncrease - Starch_max_sapwood[j])*(Volume_sapwood[j]*glucoseMolarMass)/TotalLivingBiomass[j]);
        //   starchSapwoodIncrease = Starch_max_sapwood[j] - starchSapwood[j];
        // }
        // if(starchLeafIncrease > Starch_max_leaves[j] - starchLeaf[j]) {
        //   RootExudationInst(j,s) += ((starchLeaf[j] + starchLeafIncrease - Starch_max_leaves[j])*(Volume_leaves[j]*glucoseMolarMass)/TotalLivingBiomass[j]);
        //   starchLeafIncrease = Starch_max_leaves[j] - starchLeaf[j];
        // }
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
        
        // Rcout<<j<<" LeafTLP "<< turgorLossPoint(LeafPI0[j], LeafEPS[j])<< " Leaf PI "<< osmoticWaterPotential(sugarLeaf[j], tday)<< " Conc "<< sugarLeaf[j]<< " TLPconc"<< tlpConcLeaf<<"\n";
      }
      double leafBiomassIncrement = deltaLAgrowth[j]*(1000.0/SLA[j]);
      double finerootBiomassIncrement = sum(deltaFRBgrowth);
      

      //SENESCENCE
      //Leaf senescence
      double propLeafSenescence = 0.0;
      //Leaf senescence due to age (Ca+ accumulation) only in evergreen species
      if(phenoType[j] == "progressive-evergreen") {
        propLeafSenescence = (1.0/(365.25*leafDuration[j]));
      }
      else if((phenoType[j] == "oneflush-evergreen") & (leafSenescence[j])) {
        propLeafSenescence = (1.0/leafDuration[j]); // Fraction of old leaves that die
        leafSenescence[j] = false; //To prevent further loss
      }
      else if(((phenoType[j] == "winter-deciduous") || (phenoType[j] == "winter-semideciduous")) & leafSenescence[j]) {
        propLeafSenescence = 1.0;
        leafSenescence[j] = false; //To prevent further loss
      }
      //Leaf senescence due to drought 
      double LAplc = std::min(LAexpanded, (1.0 - StemPLC[j])*leafAreaTarget[j]);
      if(LAplc<LAexpanded) {
        propLeafSenescence = std::max((LAexpanded-LAplc)/LAexpanded, propLeafSenescence); 
      }
      //Complete defoliation if RWCsymp < 0.5
      if(LAexpanded > 0.0 && LeafSympRWC[j]<0.5){
        if(allowDefoliation) {
          propLeafSenescence = 1.0;
          if(verbose) Rcout<<" [Cohort "<< j<<" defoliated ] ";
        }
      }
      double deltaLAsenescence = LAexpanded*propLeafSenescence;
      //SA senescence
      double propSASenescence = SRsapwood[j]/(1.0+15.0*exp(-0.01*H[j]));
      double deltaSASenescence = propSASenescence*SA[j];
      //FRB SENESCENCE
      NumericVector deltaFRBsenescence(numLayers, 0.0);
      for(int s=0;s<numLayers;s++) {
        double daySenescence = SRfineroot[j]*std::max(0.0,(Tsoil[s]-5.0)/20.0);
        deltaFRBsenescence[s] = fineRootBiomass[j]*V(j,s)*daySenescence;
      }
      double senescenceLeafLoss = deltaLAsenescence*(1000.0/SLA[j]);
      double senescenceFinerootLoss = sum(deltaFRBsenescence);
      
      
      //TRANSLOCATION (in mol gluc) of labile carbon
      double translocationSugarLeaf = propLeafSenescence*Volume_leaves[j]*sugarLeaf[j];
      double translocationStarchLeaf = propLeafSenescence*Volume_leaves[j]*starchLeaf[j];
      double translocationSugarSapwood = propSASenescence*Volume_sapwood[j]*starchSapwood[j];
      if(Volume_leaves[j]>0) {
        sugarLeaf[j] = ((sugarLeaf[j]*Volume_leaves[j]) - translocationSugarLeaf)/Volume_leaves[j];
        starchLeaf[j] = ((starchLeaf[j]*Volume_leaves[j]) - translocationStarchLeaf)/Volume_leaves[j];
      }
      sugarSapwood[j] = ((sugarSapwood[j]*Volume_sapwood[j]) - translocationSugarSapwood)/Volume_sapwood[j];
      starchSapwood[j] = ((starchSapwood[j]*Volume_sapwood[j]) + translocationSugarLeaf + translocationStarchLeaf + translocationSugarSapwood)/Volume_sapwood[j];

      
      //UPDATE LEAF AREA, SAPWOOD AREA, FINE ROOT BIOMASS/DISTRIBUTION
      double initialSapwoodBiomass = sapwoodStructuralBiomass(SA[j], H[j], L(j,_), V(j,_),WoodDensity[j]);
      double LAprev = LAexpanded;
      LAexpanded +=deltaLAgrowth[j] - deltaLAsenescence;
      if(LAexpanded < 0.0) {
        deltaLAsenescence -= LAexpanded;
        LAexpanded = 0.0;
      }
      LAdead += deltaLAsenescence;
      LAI_dead[j] = LAdead*N[j]/10000.0;
      LAI_expanded[j] = LAexpanded*N[j]/10000.0;
      double SAprev = SA[j];
      SA[j] = SA[j] + deltaSAgrowth[j] - deltaSASenescence; 
      NumericVector newFRB(numLayers,0.0);
      for(int s=0;s<numLayers;s++) {
        newFRB[s] = fineRootBiomass[j]*V(j,s) + deltaFRBgrowth[s] - deltaFRBsenescence[s];
      }
      fineRootBiomass[j] = sum(newFRB);
      for(int s=0;s<numLayers;s++) { 
        V(j,s) = newFRB[s]/fineRootBiomass[j];
      }
      
      //UPDATE DERIVED QUANTITIES   
      if(LAlive>0.0) {
        //Update Huber value, stem and root hydraulic conductance
        double oldstemR = 1.0/VCstem_kmax[j];
        double oldrootR = 1.0/VCroot_kmaxVEC[j];
        double oldrootprop = oldrootR/(oldrootR+oldstemR);
        
        Al2As[j] = (LAlive)/(SA[j]/10000.0);
        VCstem_kmax[j]=maximumStemHydraulicConductance(Kmax_stemxylem[j], Hmed[j], Al2As[j] ,H[j], taper); 
        
        //Update rhizosphere maximum conductance
        NumericVector VGrhizo_new = rhizosphereMaximumConductance(Ksat, newFRB, LAI_live[j], N[j],
                                                                  SRL[j], FineRootDensity[j], RLD[j]);
        for(int s=0;s<numLayers;s++) { 
          VGrhizo_kmax(j,s) = VGrhizo_new[s];
        }
        VGrhizo_kmaxVEC[j] = sum(VGrhizo_kmax(j,_));
        
        //Update root maximum conductance so that it keeps the same resistance proportion with stem conductance
        double newstemR = 1.0/VCstem_kmax[j];
        double newrootR = oldrootprop*newstemR/(1.0-oldrootprop);
        VCroot_kmaxVEC[j] = 1.0/newrootR;
        //Update coarse root soil volume
        CRSV[j] = coarseRootSoilVolumeFromConductance(Kmax_stemxylem[j], VCroot_kmaxVEC[j], Al2As[j],
                                                      V(j,_), dVec, rfc);
        //Update coarse root length and root maximum conductance
        L(j,_) = coarseRootLengthsFromVolume(CRSV[j], V(j,_), dVec, rfc);
        NumericVector xp = rootxylemConductanceProportions(L(j,_), V(j,_));
        VCroot_kmax(j,_) = VCroot_kmaxVEC[j]*xp;
        //Update Plant_kmax
        Plant_kmax[j] = 1.0/((1.0/VCleaf_kmax[j])+(1.0/VCstem_kmax[j])+(1.0/VCroot_kmaxVEC[j]));
        //Update leaf and stem osmotic water potential at full turgor
        LeafPI0[j] = osmoticWaterPotential(sugarLeaf[j], 20.0, nonSugarConcentration); //Osmotic potential at full turgor assuming RWC = 1 and 20ºC
        StemPI0[j] = osmoticWaterPotential(sugarSapwood[j], 20.0, nonSugarConcentration);
        //Update non-stomatal photosynthesis limitations
        if(nonStomatalPhotosynthesisLimitation) NSPL[j] = 1.0 - std::max(0.0, std::min(1.0, sugarLeaf[j] - 0.5)); //photosynthesis limited when conc > 0.5 and zero when conc > 1.5 mol·l-1
        else NSPL[j] = 1.0;
      }
      //Decrease PLC due to new SA growth
      if(cavitationRefill=="growth") StemPLC[j] = std::max(0.0, StemPLC[j] - (deltaSAgrowth[j]/SA[j]));
      
      
      //LEAF/FINE ROOT BIOMASS balance (g_ind)
      LeafBiomassBalance[j] = leafBiomassIncrement - senescenceLeafLoss;
      FineRootBiomassBalance[j] = finerootBiomassIncrement - senescenceFinerootLoss;

      //MORTALITY Death by carbon starvation or dessication
      double Nprev = N[j]; //Store initial density (for biomass balance)
      double Ndead_day = 0.0;
      bool dynamicCohort = true;
      bool isShrub = !NumericVector::is_na(Cover[j]);
      if((!shrubDynamics) & isShrub) dynamicCohort = false;
      if(dynamicCohort) {
        if(mortalityMode=="whole-cohort/deterministic") {
          if((sugarSapwood[j]<mortalitySugarThreshold) & allowStarvation) {
            Ndead_day = N[j];
            if(verbose) Rcout<<" [Cohort "<< j<<" died from starvation] ";
          } else if( (StemSympRWC[j] < mortalityRWCThreshold) & allowDessication) {
            Ndead_day = N[j];
            if(verbose) Rcout<<" [Cohort "<< j<<" died from dessication] ";
          }
        } else {
          double P_starv = dailyMortalityProbability(mortalityBaselineRate, sugarSapwood[j], 
                                                     mortalitySugarThreshold, allowStarvation,
                                                     0.0, 2.0);
          double P_dess = dailyMortalityProbability(mortalityBaselineRate, StemSympRWC[j], 
                                                    mortalityRWCThreshold, allowDessication,
                                                    0.0, 2.0);
          double P_day = std::max(P_dess, P_starv);
          if(mortalityMode =="density/deterministic") {
            Ndead_day = N[j]*P_day;
          } else if(mortalityMode =="whole-cohort/stochastic") {
            if(R::runif(0.0,1.0) < P_day) {
              Ndead_day = N[j];
              if(P_dess>P_starv) {
                Rcout<<" [Cohort "<< j<<" died from dessication] ";
              } else {
                Rcout<<" [Cohort "<< j<<" died from starvation] ";
              }
            }
          } else if(mortalityMode =="density/stochastic") {
            // NumericVector nv = Rcpp::runif(round(N[j]), 0.0,1.0);
            // for(int k=0;k<nv.size();k++) if(nv[k]< P_day) Ndead_day += 1.0;
            Ndead_day = R::rbinom(round(N[j]), P_day);
            // Rcout<< j<< " "<< P_day<< " "<< N[j]<< " "<< Ndead_day<< " "<< R::rbinom(750, 2.82309e-05) << "\n";
          }
          // if(Ndead_day>0.0) Rcout<< j<< " "<< P_day<< " "<< Ndead_day<< "\n";
        }
        // Update density and increase the number of dead plants
        Ndead_day = std::min(Ndead_day, N[j]);
        double Cdead_day = Cover[j]*(Ndead_day/N[j]);
        N[j] = N[j] - Ndead_day;
        N_dead[j] = N_dead[j] + Ndead_day;
        if(isShrub) {
          Cover[j] = std::max(0.0, Cover[j] - Cdead_day);
          Cover_dead[j] = Cover_dead[j] + Cdead_day;
        }
        //Update LAI dead and LAI expanded as a result of density decrease
        double LAI_change = LAexpanded*Ndead_day/10000.0;
        LAI_dead[j] = LAI_dead[j] + LAI_change;
        LAI_expanded[j] = LAI_expanded[j] - LAI_change;
      }


      //UPDATE TARGETS
      //Set leaf area target if bud formation is allowed
      if(budFormation[j]) {
        if(allocationStrategy == "Plant_kmax") {
          leafAreaTarget[j] = LAlive*(Plant_kmax[j]/allocationTarget[j]);
        } else if(allocationStrategy =="Al2As") {
          leafAreaTarget[j] = (SA[j]/10000.0)*allocationTarget[j];
        }
        LAI_live[j] =  leafAreaTarget[j]*N[j]/10000.0;
      }
      //Update fine root biomass target     
      if(LAI_live[j]>0) {
        NumericVector VGrhizo_target(numLayers,0.0);
        for(int s=0;s<numLayers;s++) {
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
      
      
      //OUTPUT VARIABLES
      SapwoodArea[j] = SA[j];
      LeafArea[j] = LAexpanded;
      FineRootArea[j] = fineRootBiomass[j]*specificRootSurfaceArea(SRL[j], FineRootDensity[j])*1e-4;
      LAgrowth[j] += deltaLAgrowth[j]/SAprev;//Store Leaf area growth rate in relation to sapwood area (m2·cm-2·d-1)
      SAgrowth[j] += deltaSAgrowth[j]/SAprev; //Store sapwood area growth rate (cm2·cm-2·d-1)
      FRAgrowth[j] = sum(deltaFRBgrowth)*specificRootSurfaceArea(SRL[j], FineRootDensity[j])*1e-4/SA[j];//Store fine root area growth rate (m2·cm-2·d-1)

    }
  }
  //UPDATE STRUCTURAL VARIABLES
  updateStructuralVariables(x, deltaSAgrowth);
  
  //RECALCULATE storage concentrations (SA, LA and H may have changed)
  for(int j=0;j<numCohorts;j++){
    double newVolumeSapwood = sapwoodStorageVolume(SA[j], H[j], L(j,_),V(j,_),WoodDensity[j], conduit2sapwood[j]);
    double newVolumeLeaves = leafStorageVolume(LAI_expanded[j],  N[j], SLA[j], LeafDensity[j]);
    if(newVolumeLeaves > 0.0) {
      sugarLeaf[j] = sugarLeaf[j]*(Volume_leaves[j]/newVolumeLeaves);
      starchLeaf[j] = starchLeaf[j]*(Volume_leaves[j]/newVolumeLeaves); 
    } else {
      sugarLeaf[j] = 0.0;
      starchLeaf[j] = 0.0;
    }
    sugarSapwood[j] = sugarSapwood[j]*(Volume_sapwood[j]/newVolumeSapwood);
    starchSapwood[j] = starchSapwood[j]*(Volume_sapwood[j]/newVolumeSapwood); 
    
    //OUTPUT VARIABLES
    PlantSugarLeaf[j] = sugarLeaf[j];
    PlantStarchLeaf[j] = starchLeaf[j];
    PlantSugarSapwood[j] = sugarSapwood[j];
    PlantStarchSapwood[j] = starchSapwood[j];
  }
  
  //CLOSE BIOMASS BALANCE
  closePlantBiomassBalance(plantBiomassBalance, x, 
                      LabileCarbonBalance, LeafBiomassBalance, FineRootBiomassBalance);
  
  
  //Update pool proportions??
  if(plantWaterPools) {
    NumericVector poolProportions = Rcpp::as<Rcpp::NumericVector>(belowdf["poolProportions"]);
    // for(int j=0;j<numCohorts;j++) poolProportions[j] = LAI_live[j]/sum(LAI_live);
    //Update RHOP
    List RHOP = belowLayers["RHOP"];
    List newRHOP = horizontalProportions(poolProportions, CRSV, N, V,dVec, rfc);
    for(int j=0;j<numCohorts;j++) RHOP[j] = newRHOP[j];
  }
  
  GrossPhotosynthesisInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,numSteps));
  MaintenanceRespirationInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,numSteps));
  LabileCarbonBalanceInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,numSteps));
  GrowthCostsInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,numSteps));
  PlantSugarLeafInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,numSteps));
  PlantStarchLeafInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,numSteps));
  PlantSugarSapwoodInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,numSteps));
  PlantStarchSapwoodInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,numSteps));
  PlantSugarTransportInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,numSteps));
  RootExudationInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,numSteps));
  List labileCBInst = List::create(
    _["GrossPhotosynthesis"] = GrossPhotosynthesisInst,
    _["MaintenanceRespiration"] = MaintenanceRespirationInst,
    _["GrowthCosts"] = GrowthCostsInst,
    _["RootExudation"] = RootExudationInst,
    _["LabileCarbonBalance"] = LabileCarbonBalanceInst,
    _["SugarLeaf"] = PlantSugarLeafInst,
    _["StarchLeaf"] = PlantStarchLeafInst,
    _["SugarSapwood"] = PlantSugarSapwoodInst,
    _["StarchSapwood"] = PlantStarchSapwoodInst,
    _["SugarTransport"] = PlantSugarTransportInst
  );
  
  DataFrame labileCarbonBalance = DataFrame::create(_["GrossPhotosynthesis"] = GrossPhotosynthesis,
                                   _["MaintenanceRespiration"] = MaintenanceRespiration,
                                   _["GrowthCosts"] = GrowthCosts,
                                   _["RootExudation"] = RootExudation,
                                   _["LabileCarbonBalance"] = LabileCarbonBalance,
                                   _["SugarLeaf"] = PlantSugarLeaf,
                                   _["StarchLeaf"] = PlantStarchLeaf,
                                   _["SugarSapwood"] = PlantSugarSapwood,
                                   _["StarchSapwood"] = PlantStarchSapwood,
                                   _["SugarTransport"] = PlantSugarTransport,
                                   _["StemPI0"] = clone(StemPI0), //Store a copy of the current osmotic potential at full turgor
                                   _["LeafPI0"] = clone(LeafPI0));
  labileCarbonBalance.attr("row.names") = above.attr("row.names");
  
  //Final Biomass compartments
  DataFrame plantStructure = DataFrame::create(
    _["LeafArea"] = LeafArea,
    _["SapwoodArea"] = SapwoodArea,
    _["FineRootArea"] = FineRootArea,
    _["FineRootBiomass"] = clone(fineRootBiomass),
    _["DBH"] = clone(DBH),
    _["Height"] = clone(H)
  );

  DataFrame plantGrowth = DataFrame::create(
    _["SAgrowth"] = SAgrowth,
    _["LAgrowth"] = LAgrowth,
    _["FRAgrowth"] = FRAgrowth
  );
  plantGrowth.attr("row.names") = above.attr("row.names");
  
  List l = List::create(_["cohorts"] = clone(cohorts),
                        _["WaterBalance"] = spwbOut["WaterBalance"], 
                        _["EnergyBalance"] = spwbOut["EnergyBalance"],
                        _["Soil"] = spwbOut["Soil"], 
                        _["Stand"] = spwbOut["Stand"], 
                        _["Plants"] = spwbOut["Plants"],
                        _["LabileCarbonBalance"] = labileCarbonBalance,
                        _["PlantBiomassBalance"] = plantBiomassBalance,
                        _["PlantStructure"] = plantStructure,                        
                        _["PlantGrowth"] = plantGrowth,
                        _["RhizoPsi"] = spwbOut["RhizoPsi"],
                        _["SunlitLeaves"] = spwbOut["SunlitLeaves"],
                        _["ShadeLeaves"] = spwbOut["ShadeLeaves"],
                        _["ExtractionInst"] = spwbOut["ExtractionInst"],
                        _["PlantsInst"] = spwbOut["PlantsInst"],
                        _["SunlitLeavesInst"] = spwbOut["SunlitLeavesInst"],
                        _["ShadeLeavesInst"] = spwbOut["ShadeLeavesInst"],
                        _["LabileCarbonBalanceInst"] = labileCBInst,
                        _["LightExtinction"] = spwbOut["LightExtinction"],
                        _["CanopyTurbulence"] = spwbOut["CanopyTurbulence"]);
  l.attr("class") = CharacterVector::create("growth_day","list");
  return(l);
}


// [[Rcpp::export("growth_day")]]
List growthDay(List x, CharacterVector date, double tmin, double tmax, double rhmin, 
               double rhmax, double rad, double wind, 
               double latitude, double elevation, double slope, double aspect,  
               double prec, double runon=0.0) {
  //Control parameters
  List control = x["control"];
  bool verbose = control["verbose"];
  bool modifyInput = control["modifyInput"];
  bool leafPhenology = control["leafPhenology"];
  String transpirationMode = control["transpirationMode"];
  double Catm = control["Catm"];
  
  //Will not modify input x 
  if(!modifyInput) {
    x = clone(x);
  }
  
  std::string c = as<std::string>(date[0]);
  int J = meteoland::radiation_julianDay(std::atoi(c.substr(0, 4).c_str()),std::atoi(c.substr(5,2).c_str()),std::atoi(c.substr(8,2).c_str()));
  double delta = meteoland::radiation_solarDeclination(J);
  double solarConstant = meteoland::radiation_solarConstant(J);
  double tday = meteoland::utils_averageDaylightTemperature(tmin, tmax);
  double latrad = latitude * (M_PI/180.0);
  double asprad = aspect * (M_PI/180.0);
  double slorad = slope * (M_PI/180.0);
  double photoperiod = meteoland::radiation_daylength(latrad, 0.0, 0.0, delta);
  double pet = meteoland::penman(latrad, elevation, slorad, asprad, J, tmin, tmax, rhmin, rhmax, rad, wind);
  //Derive doy from date  
  int J0101 = meteoland::radiation_julianDay(std::atoi(c.substr(0, 4).c_str()),1,1);
  int doy = J - J0101+1;
  if(NumericVector::is_na(wind)) wind = control["defaultWindSpeed"]; 
  if(wind<0.1) wind = 0.1; //Minimum windspeed abovecanopy
  
  //Update phenology
  if(leafPhenology) {
    updatePhenology(x, doy, photoperiod, tday);
    updateLeaves(x, wind, true);
  }
  
  double er = erFactor(doy, pet, prec);
  List s;
  if(transpirationMode=="Granier") {
    NumericVector meteovec = NumericVector::create(
      Named("tday") = tday, 
      Named("prec") = prec,
      Named("rad") = rad, 
      Named("pet") = pet,
      Named("er") = er);
    s = growthDay1(x, meteovec, 
                   elevation, 
                   runon, verbose);
  } else {
    NumericVector meteovec = NumericVector::create(
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
      Named("pet") = pet,
      Named("er") = er);
    s = growthDay2(x, meteovec, 
                 latitude, elevation, slope, aspect,
                 solarConstant, delta, 
                 runon, verbose);
  }
  // Rcout<<"hola4\n";
  return(s);
}


void checkgrowthInput(List x, String transpirationMode, String soilFunctions) {
  
  List soil = x["soil"];
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
  if(transpirationMode=="Sperry"){
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
  // if(!paramsTransp.containsElementNamed("pRootDisc")) stop("pRootDisc missing in growthInput$paramsTransp");
  if(transpirationMode=="Granier") {
    if(!paramsTranspiration.containsElementNamed("Psi_Extract")) stop("Psi_Extract missing in growthInput$paramsTransp");
    if(!paramsTranspiration.containsElementNamed("WUE")) stop("WUE missing in growthInput$paramsTransp");
  } else if(transpirationMode=="Sperry") {
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
  if(!soil.containsElementNamed("dVec")) stop("dVec missing in soil");
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

// [[Rcpp::export("growth")]]
List growth(List x, DataFrame meteo, double latitude, double elevation = NA_REAL, double slope = NA_REAL, double aspect = NA_REAL) {

  //Control params 
  List control =x["control"];  
  String transpirationMode = control["transpirationMode"];
  String soilFunctions = control["soilFunctions"];
  bool verbose = control["verbose"];
  bool modifyInput = control["modifyInput"];
  bool subdailyResults = control["subdailyResults"];
  bool leafPhenology = control["leafPhenology"];
  bool unlimitedSoilWater = control["unlimitedSoilWater"];
  bool multiLayerBalance = control["multiLayerBalance"];
  bool shrubDynamics = control["shrubDynamics"];
  checkgrowthInput(x, transpirationMode, soilFunctions);

  //Store input
  List growthInput;
  if(modifyInput) {
    growthInput = clone(x); // Will modify x and return the unmodified object
  } else {
    x = clone(x);
    growthInput = x; //Will not modified input x and return the final modified object
  }
  
  //Soil params 
  List soil = x["soil"];
  
  //Cohort info
  DataFrame cohorts = Rcpp::as<Rcpp::DataFrame>(x["cohorts"]);
  IntegerVector SP = Rcpp::as<Rcpp::IntegerVector>(cohorts["SP"]);
  IntegerVector SPunique = uniqueSpp(SP);
  int numSpecies = SPunique.size();
  
  if(NumericVector::is_na(latitude)) stop("Value for 'latitude' should not be missing.");
  double latrad = latitude * (M_PI/180.0);
  
  //Meteorological input    
  NumericVector MinTemperature, MaxTemperature;
  NumericVector MinRelativeHumidity, MaxRelativeHumidity;
  NumericVector Radiation;
  if(!meteo.containsElementNamed("Precipitation")) stop("Please include variable 'Precipitation' in weather input.");
  NumericVector Precipitation = meteo["Precipitation"];
  int numDays = Precipitation.size();
  if(!meteo.containsElementNamed("MeanTemperature")) stop("Please include variable 'MeanTemperature' in weather input.");
  NumericVector MeanTemperature = meteo["MeanTemperature"];
  NumericVector WindSpeed(numDays, NA_REAL);
  if(meteo.containsElementNamed("WindSpeed")) WindSpeed = meteo["WindSpeed"];
  NumericVector PET = NumericVector(numDays,0.0);
  NumericVector CO2(Precipitation.length(), NA_REAL);
  IntegerVector DOY, JulianDay;
  NumericVector Photoperiod;
  bool doy_input = false, photoperiod_input = false, julianday_input = false;
  if(transpirationMode=="Granier") {
    if(!meteo.containsElementNamed("PET")) stop("Please include variable 'PET' in weather input.");
    PET = meteo["PET"];
    if(control["snowpack"]) {
      if(!meteo.containsElementNamed("Radiation")) stop("If 'snowpack = TRUE', variable 'Radiation' must be provided.");
      else Radiation = meteo["Radiation"];
    }
  } else if(transpirationMode=="Sperry") {
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
    if(meteo.containsElementNamed("CO2")) {
      CO2 = meteo["CO2"];
      if(verbose) {
        Rcout<<"CO2 taken from input column 'CO2'\n";
      }
    }
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
  }
  CharacterVector dateStrings = meteo.attr("row.names");
  
  if(!doy_input) DOY = date2doy(dateStrings);
  if(!photoperiod_input) Photoperiod = date2photoperiod(dateStrings, latrad);
  
  
  //Canopy scalars
  DataFrame canopy = Rcpp::as<Rcpp::DataFrame>(x["canopy"]);
  
  //Aboveground parameters  
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  int numCohorts = SP.size();

  //Belowground state variables  
  List below = Rcpp::as<Rcpp::List>(x["below"]);
  NumericVector Z95 = Rcpp::as<Rcpp::NumericVector>(below["Z95"]);


  //Base parameters
  DataFrame paramsPhenology = Rcpp::as<Rcpp::DataFrame>(x["paramsPhenology"]);
  NumericVector Sgdd = Rcpp::as<Rcpp::NumericVector>(paramsPhenology["Sgdd"]);
  
  
  // Rcout<<"3";


  //Soil input
  NumericVector Water_FC = waterFC(soil, soilFunctions);

  
  
  
  //Detailed subday results
  List subdailyRes(numDays);
  
  //EnergyBalance output variables
  DataFrame DEB = defineEnergyBalanceDailyOutput(meteo);
  DataFrame DT = defineTemperatureDailyOutput(meteo);
  NumericMatrix DLT;
  if(transpirationMode=="Sperry") DLT =  defineTemperatureLayersDailyOutput(meteo, canopy);
  
  //Plant carbon output variables
  NumericMatrix LabileCarbonBalance(numDays, numCohorts);
  NumericMatrix MaintenanceRespiration(numDays, numCohorts);
  NumericMatrix GrowthCosts(numDays, numCohorts);
  NumericMatrix PlantSugarLeaf(numDays, numCohorts);
  NumericMatrix PlantStarchLeaf(numDays, numCohorts);
  NumericMatrix PlantSugarSapwood(numDays, numCohorts);
  NumericMatrix PlantStarchSapwood(numDays, numCohorts);
  NumericMatrix PlantSugarTransport(numDays, numCohorts);
  NumericMatrix SapwoodArea(numDays, numCohorts);
  NumericMatrix LeafArea(numDays, numCohorts);
  NumericMatrix FineRootArea(numDays, numCohorts);
  NumericMatrix FineRootBiomass(numDays, numCohorts);
  NumericMatrix DBH(numDays, numCohorts);
  NumericMatrix Height(numDays, numCohorts);
  NumericMatrix LabileBiomass(numDays, numCohorts);
  NumericMatrix TotalBiomass(numDays, numCohorts);
  NumericMatrix SAgrowth(numDays, numCohorts);
  NumericMatrix LAgrowth(numDays, numCohorts);
  NumericMatrix FRAgrowth(numDays, numCohorts);
  NumericMatrix GrossPhotosynthesis(numDays, numCohorts);
  NumericMatrix PlantLAIexpanded(numDays, numCohorts), PlantLAIdead(numDays, numCohorts), PlantLAIlive(numDays, numCohorts);
  NumericMatrix StemPI0(numDays, numCohorts), LeafPI0(numDays, numCohorts);
  NumericMatrix RootExudation(numDays, numCohorts);
  NumericMatrix StructuralBiomassBalance(numDays, numCohorts);
  NumericMatrix LabileBiomassBalance(numDays, numCohorts);
  NumericMatrix PlantBiomassBalance(numDays, numCohorts);
  NumericMatrix MortalityBiomassLoss(numDays, numCohorts);
  NumericMatrix CohortBiomassBalance(numDays, numCohorts);
  NumericMatrix StandBiomassBalance(numDays, 5);
  
  //Water balance output variables
  DataFrame DWB = defineWaterBalanceDailyOutput(meteo, PET, transpirationMode);
  DataFrame SWB = defineSoilWaterBalanceDailyOutput(meteo, soil, transpirationMode);
  
  
  NumericVector LAI(numDays), LAIlive(numDays), LAIexpanded(numDays), LAIdead(numDays);
  NumericVector Cm(numDays);
  NumericVector LgroundPAR(numDays);
  NumericVector LgroundSWR(numDays);

  //Plant water output variables
  List sunlitDO = defineSunlitShadeLeavesDailyOutput(meteo, above);
  List shadeDO = defineSunlitShadeLeavesDailyOutput(meteo, above);
  List plantDWOL = definePlantWaterDailyOutput(meteo, above, soil, control);
  NumericVector EplantCohTot(numCohorts, 0.0);

  
  //Count years (times structural variables will be updated)
  int numYears = 0;
  for(int i=0;i<numDays;i++) {
    if(((DOY[i]==1) && (i>0)) | ((i==(numDays-1)) && (DOY[i]>=365))) numYears = numYears + 1;
  }

  NumericVector initialContent = water(soil, soilFunctions);
  double initialSnowContent = soil["SWE"];
  DataFrame ccIni_m2 = carbonCompartments(x, "g_m2");
  double cohortBiomassBalanceSum = 0.0;
  double initialCohortBiomass = sum(Rcpp::as<Rcpp::NumericVector>(ccIni_m2["TotalBiomass"]));
  
  if(verbose) {
    Rcout<<"Initial plant cohort biomass (g/m2): "<<initialCohortBiomass<<"\n";
    Rcout<<"Initial soil water content (mm): "<< sum(initialContent)<<"\n";
    Rcout<<"Initial snowpack content (mm): "<< initialSnowContent<<"\n";
  }
  
  
  if(verbose) Rcout << "Performing daily simulations\n";
  List s;
  int iyear = 0;
  for(int i=0;i<numDays;i++) {
    if(verbose) {
      if(DOY[i]==1 || i==0) {
        std::string c = as<std::string>(dateStrings[i]);
        Rcout<<"\n [Year "<< c.substr(0, 4)<< "]: ";
      } 
      else if(i%10 == 0) Rcout<<".";//<<i;
    } 
    
    double wind = WindSpeed[i];
    if(NumericVector::is_na(wind)) wind = control["defaultWindSpeed"]; //Default 1 m/s -> 10% of fall every day
    if(wind<0.1) wind = 0.1; //Minimum windspeed abovecanopy
    
    if(unlimitedSoilWater) {
      NumericVector W = soil["W"];
      for(int h=0;h<W.size();h++) W[h] = 1.0;
    }
    
    //1. Phenology (only leaf fall)
    if(leafPhenology) {
      updatePhenology(x, DOY[i], Photoperiod[i], MeanTemperature[i]);
      updateLeaves(x, wind, true);
    }

    //2. Water balance and photosynthesis
    if(transpirationMode=="Granier") {
      NumericVector meteovec = NumericVector::create(
        Named("tday") = MeanTemperature[i], 
        Named("prec") = Precipitation[i],
        Named("rad") = Radiation[i], 
        Named("pet") = PET[i],
        Named("er") = erFactor(DOY[i], PET[i], Precipitation[i]));
      s = growthDay1(x, meteovec,  
                     elevation, 
                     0.0, false); //No Runon in simulations for a single cell
    } else if(transpirationMode=="Sperry") {
      //Julian day from either input column or date
      int J = NA_INTEGER;
      if(julianday_input) J = JulianDay[i];
      if(IntegerVector::is_na(J)){
        std::string c = as<std::string>(dateStrings[i]);
        J = meteoland::radiation_julianDay(std::atoi(c.substr(0, 4).c_str()),std::atoi(c.substr(5,2).c_str()),std::atoi(c.substr(8,2).c_str())); 
      }
      double delta = meteoland::radiation_solarDeclination(J);
      double solarConstant = meteoland::radiation_solarConstant(J);
      double latrad = latitude * (M_PI/180.0);
      if(NumericVector::is_na(aspect)) aspect = 0.0;
      if(NumericVector::is_na(slope)) slope = 0.0;
      double asprad = aspect * (M_PI/180.0);
      double slorad = slope * (M_PI/180.0);
      double tmin = MinTemperature[i];
      double tmax = MaxTemperature[i];
      double tmaxPrev = tmax;
      double tminPrev = tmin;
      double tminNext = tmin;
      if(i>0) {
        tmaxPrev = MaxTemperature[i-1];
        tminPrev = MinTemperature[i-1];
      }
      if(i<(numDays-1)) tminNext = MinTemperature[i+1]; 
      double rhmin = MinRelativeHumidity[i];
      double rhmax = MaxRelativeHumidity[i];
      double rad = Radiation[i];
      double Catm = CO2[i];
      if(NumericVector::is_na(Catm)) Catm = control["Catm"];
      PET[i] = meteoland::penman(latrad, elevation, slorad, asprad, J, tmin, tmax, rhmin, rhmax, rad, wind);
      NumericVector meteovec = NumericVector::create(
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
        Named("pet") = PET[i],
        Named("er") = erFactor(DOY[i], PET[i], Precipitation[i]));
      s = growthDay2(x, meteovec, 
                   latitude, elevation, slope, aspect,
                   solarConstant, delta, 
                   0.0, verbose);

      fillEnergyBalanceTemperatureDailyOutput(DEB,DT,DLT,s,i, multiLayerBalance);
    }    
    
    fillPlantWaterDailyOutput(plantDWOL, sunlitDO, shadeDO, s, i, transpirationMode);
    fillWaterBalanceDailyOutput(DWB, s,i, transpirationMode);
    fillSoilWaterBalanceDailyOutput(SWB, soil, s,
                                    i, numDays, transpirationMode, soilFunctions);
    
    List stand = s["Stand"];
    LgroundPAR[i] = stand["LgroundPAR"];
    LgroundSWR[i] = stand["LgroundSWR"];
    LAI[i] = stand["LAI"];
    LAIlive[i] = stand["LAIlive"];
    LAIexpanded[i] = stand["LAIexpanded"];
    LAIdead[i] = stand["LAIdead"];
    Cm[i] = stand["Cm"];
    
    List sb = s["Soil"];
    List db = s["WaterBalance"];
    List Plants = s["Plants"];
    DataFrame cb = Rcpp::as<Rcpp::DataFrame>(s["LabileCarbonBalance"]);
    DataFrame bb = Rcpp::as<Rcpp::DataFrame>(s["PlantBiomassBalance"]);
    DataFrame pg = Rcpp::as<Rcpp::DataFrame>(s["PlantGrowth"]);
    DataFrame ps = Rcpp::as<Rcpp::DataFrame>(s["PlantStructure"]);
    
    
    //4. Assemble output
    LabileCarbonBalance(i,_) = Rcpp::as<Rcpp::NumericVector>(cb["LabileCarbonBalance"]);
    MaintenanceRespiration(i,_) = Rcpp::as<Rcpp::NumericVector>(cb["MaintenanceRespiration"]);
    GrowthCosts(i,_) = Rcpp::as<Rcpp::NumericVector>(cb["GrowthCosts"]);
    GrossPhotosynthesis(i,_) = Rcpp::as<Rcpp::NumericVector>(cb["GrossPhotosynthesis"]);
    PlantSugarLeaf(i,_) = Rcpp::as<Rcpp::NumericVector>(cb["SugarLeaf"]);
    PlantStarchLeaf(i,_) = Rcpp::as<Rcpp::NumericVector>(cb["StarchLeaf"]);
    PlantSugarSapwood(i,_) = Rcpp::as<Rcpp::NumericVector>(cb["SugarSapwood"]);
    PlantStarchSapwood(i,_) = Rcpp::as<Rcpp::NumericVector>(cb["StarchSapwood"]);
    PlantSugarTransport(i,_) = Rcpp::as<Rcpp::NumericVector>(cb["SugarTransport"]);
    StemPI0(i,_) = Rcpp::as<Rcpp::NumericVector>(cb["StemPI0"]); 
    LeafPI0(i,_) = Rcpp::as<Rcpp::NumericVector>(cb["LeafPI0"]); 
    RootExudation(i,_) = Rcpp::as<Rcpp::NumericVector>(cb["RootExudation"]);
    
    SapwoodArea(i,_) = Rcpp::as<Rcpp::NumericVector>(ps["SapwoodArea"]);
    LeafArea(i,_) = Rcpp::as<Rcpp::NumericVector>(ps["LeafArea"]);
    FineRootBiomass(i,_) = Rcpp::as<Rcpp::NumericVector>(ps["FineRootBiomass"]);
    if(transpirationMode=="Sperry") {
      FineRootArea(i,_) = Rcpp::as<Rcpp::NumericVector>(ps["FineRootArea"]);
    }
    DBH(i,_) = Rcpp::as<Rcpp::NumericVector>(ps["DBH"]);
    Height(i,_) = Rcpp::as<Rcpp::NumericVector>(ps["Height"]);
    
    StructuralBiomassBalance(i,_) = Rcpp::as<Rcpp::NumericVector>(bb["StructuralBiomassBalance"]);
    LabileBiomassBalance(i,_) = Rcpp::as<Rcpp::NumericVector>(bb["LabileBiomassBalance"]);
    PlantBiomassBalance(i,_) = Rcpp::as<Rcpp::NumericVector>(bb["PlantBiomassBalance"]);
    MortalityBiomassLoss(i,_) = Rcpp::as<Rcpp::NumericVector>(bb["MortalityBiomassLoss"]);
    CohortBiomassBalance(i,_) = Rcpp::as<Rcpp::NumericVector>(bb["CohortBiomassBalance"]);
    StandBiomassBalance(i,_) = standLevelBiomassBalance(bb);
    cohortBiomassBalanceSum += sum(CohortBiomassBalance(i,_));

    LAgrowth(i,_) = Rcpp::as<Rcpp::NumericVector>(pg["LAgrowth"]);
    SAgrowth(i,_) = Rcpp::as<Rcpp::NumericVector>(pg["SAgrowth"]);
    if(transpirationMode=="Sperry") {
      FRAgrowth(i,_) = Rcpp::as<Rcpp::NumericVector>(pg["FRAgrowth"]);
    }
      

    //5 Update structural variables
    if(((i<(numDays-1)) && (DOY[i+1]==1)) || (i==(numDays-1))) { 
      //reset ring structures
      List ringList = x["internalRings"];
      for(int j=0;j<numCohorts; j++) ringList[j] = initialize_ring();
    }

    if(subdailyResults) {
      subdailyRes[i] = clone(s);
    }
  }
  if(verbose) Rcout << "\n\n";
  
  // Check biomass balance
  DataFrame ccFin_m2 = carbonCompartments(x, "g_m2");
  double finalCohortBiomass = sum(Rcpp::as<Rcpp::NumericVector>(ccFin_m2["TotalBiomass"]));
  if(verbose) {
    Rcout<<"Final plant biomass (g/m2): "<<finalCohortBiomass<<"\n";
    Rcout<<"Change in plant biomass (g/m2): " << finalCohortBiomass - initialCohortBiomass <<"\n";
    Rcout<<"Plant biomass balance result (g/m2): " <<  cohortBiomassBalanceSum<<"\n";
    Rcout<<"Plant biomass balance components:\n";
    
    Rcout<<"  Structural balance (g/m2) "  <<round(sum(StandBiomassBalance(_,0)))<<" Labile balance (g/m2) "  <<round(sum(StandBiomassBalance(_,1))) <<"\n";
    Rcout<<"  Plant individual balance (g/m2) "  <<round(sum(StandBiomassBalance(_,2)))<<" Mortality loss (g/m2) "  <<round(sum(StandBiomassBalance(_,3))) <<"\n";
    
    printWaterBalanceResult(DWB, plantDWOL, soil, soilFunctions,
                            initialContent, initialSnowContent,
                            transpirationMode);
  }
  
  
  //Add matrix dimnames
  LabileCarbonBalance.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names"));
  GrossPhotosynthesis.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names"));
  MaintenanceRespiration.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names"));
  GrowthCosts.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names"));
  PlantSugarLeaf.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names"));
  PlantStarchLeaf.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names")) ;
  PlantSugarSapwood.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names")) ;
  PlantStarchSapwood.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names")) ;
  PlantSugarTransport.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names")) ;
  SapwoodArea.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names")) ;
  LeafArea.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names")) ;
  FineRootArea.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names")) ;
  FineRootBiomass.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names"));
  DBH.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names"));
  Height.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names"));
  LAgrowth.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names")) ;
  SAgrowth.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names")) ;
  FRAgrowth.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names"));
  StemPI0.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names"));
  LeafPI0.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names"));
  RootExudation.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names"));
  StructuralBiomassBalance.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names"));
  LabileBiomassBalance.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names"));
  PlantBiomassBalance.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names"));
  MortalityBiomassLoss.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names"));
  CohortBiomassBalance.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names"));

  StandBiomassBalance.attr("dimnames") = List::create(meteo.attr("row.names"), 
                           CharacterVector::create("StructuralBalance", "LabileBalance", "PlantBalance", "MortalityLoss", "CohortBalance"));
  
  subdailyRes.attr("names") = meteo.attr("row.names") ;
  
  NumericVector topo = NumericVector::create(elevation, slope, aspect);
  topo.attr("names") = CharacterVector::create("elevation", "slope", "aspect");
  
  Rcpp::DataFrame Stand = DataFrame::create(_["LAI"]=LAI, _["LAIlive"]=LAIlive, _["LAIexpanded"]=LAIexpanded,_["LAIdead"]=LAIdead,
                                            _["Cm"]=Cm, 
                                            _["LgroundPAR"] = LgroundPAR, _["LgroundSWR"] = LgroundSWR);
  Stand.attr("row.names") = meteo.attr("row.names");

  
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
    Named("SugarTransport") = PlantSugarTransport,
    Named("LeafPI0") = LeafPI0,
    Named("StemPI0") = StemPI0
  );
  List plantBiomassBalance = List::create(_["StructuralBiomassBalance"] = StructuralBiomassBalance,
                                     _["LabileBiomassBalance"] = LabileBiomassBalance,
                                     _["PlantBiomassBalance"] = PlantBiomassBalance,
                                     _["MortalityBiomassLoss"] = MortalityBiomassLoss,
                                     _["CohortBiomassBalance"] = CohortBiomassBalance);
  
  List plantGrowth, plantStructure;
  
  List l;
  if(transpirationMode=="Granier") {
    plantStructure = List::create(Named("LeafArea") = LeafArea,
                                  Named("SapwoodArea")=SapwoodArea,
                                  Named("FineRootBiomass")=FineRootBiomass,
                                  Named("DBH") = DBH,
                                  Named("Height") = Height);
    plantGrowth = List::create(Named("LAgrowth") = LAgrowth,
                               Named("SAgrowth") = SAgrowth);
    l = List::create(Named("latitude") = latitude,
                     Named("topography") = topo,
                     Named("growthInput") = growthInput,
                     Named("WaterBalance")=DWB, 
                     Named("BiomassBalance") = StandBiomassBalance,
                     Named("Soil")=SWB,
                     Named("Stand")=Stand,
                     Named("Plants") = plantDWOL,
                     Named("LabileCarbonBalance") = labileCarbonBalance,
                     Named("PlantBiomassBalance") = plantBiomassBalance,
                     Named("PlantStructure") = plantStructure,
                     Named("PlantGrowth") = plantGrowth,
                     Named("subdaily") =  subdailyRes);
    
  } else {
    plantStructure = List::create(Named("LeafArea") = LeafArea,
                                  Named("SapwoodArea")=SapwoodArea,
                                  Named("FineRootArea") = FineRootArea,
                                  Named("FineRootBiomass") = FineRootBiomass,
                                  Named("DBH") = DBH,
                                  Named("Height") = Height);
    plantGrowth = List::create(Named("LAgrowth") = LAgrowth,
                               Named("SAgrowth") = SAgrowth,
                               Named("FRAgrowth") = FRAgrowth);
    l = List::create(Named("latitude") = latitude,
                   Named("topography") = topo,
                   Named("growthInput") = growthInput,
                   Named("WaterBalance")=DWB, 
                   Named("EnergyBalance") = DEB,
                   Named("BiomassBalance") = StandBiomassBalance,
                   Named("Temperature") = DT,
                   Named("TemperatureLayers") = NA_REAL,
                   Named("Soil")=SWB,
                   Named("Stand")=Stand,
                   Named("Plants") = plantDWOL,
                   Named("SunlitLeaves") = sunlitDO,
                   Named("ShadeLeaves") = shadeDO,
                   Named("LabileCarbonBalance") = labileCarbonBalance,
                   Named("PlantBiomassBalance") = plantBiomassBalance,
                   Named("PlantStructure") = plantStructure,
                   Named("PlantGrowth") = plantGrowth,
                   Named("subdaily") =  subdailyRes);
    if(multiLayerBalance) l["TemperatureLayers"] = DLT;
  }
  l.attr("class") = CharacterVector::create("growth","list");
  return(l);
}
