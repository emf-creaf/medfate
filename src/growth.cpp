// [[Rcpp::interfaces(r,cpp)]]
#define STRICT_R_HEADERS
#include <Rcpp.h>
#include <numeric>
#include "lightextinction.h"
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
// [[Rcpp::export("mortality_dailyProbability")]]
double dailyMortalityProbability(double stressValue, double stressThreshold) {
  double exponent = 40.0;
  double y = (stressValue - stressThreshold);
  double P_annual = 1.0 - exp(exponent*y)/(1.0 + exp(exponent*y));
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
      if(budFormation[j]) {
        // Rcout <<j<< " "<< ltba[j]<< " "<<leafAreaTarget[j];
        leafAreaTarget[j] = SLA[j]*(Afbt[j]*pow(std::min(50.0,DBH[j]), Bfbt[j])*exp(Cfbt[j]*ltba[j]));
        leafAreaTarget[j] = leafAreaTarget[j] * exp(-0.0001*N[j]);//Correct for high density packing
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
List growthDayInner(List x, NumericVector meteovec, 
                    double latitude, double elevation, double slope = NA_REAL, double aspect = NA_REAL,
                    double solarConstant = NA_REAL, double delta = NA_REAL, 
                    double runon=0.0, bool verbose = false) {
  
  
  //Get previous PLC so that defoliation occurs only when PLC increases
  //Internal state variables
  DataFrame internalWater = Rcpp::as<Rcpp::DataFrame>(x["internalWater"]);
  NumericVector StemPLCprev = clone(Rcpp::as<Rcpp::NumericVector>(internalWater["StemPLC"]));
  
  //Control params
  List control = x["control"];  
  
  String transpirationMode = control["transpirationMode"];
  String soilFunctions = control["soilFunctions"];
  String mortalityMode = control["mortalityMode"];
  bool subdailyCarbonBalance = control["subdailyCarbonBalance"];
  double mortalityRelativeSugarThreshold= control["mortalityRelativeSugarThreshold"];
  double mortalityRWCThreshold= control["mortalityRWCThreshold"];
  bool allowDessication = control["allowDessication"];
  bool allowStarvation = control["allowStarvation"];
  bool sinkLimitation = control["sinkLimitation"];
  bool shrubDynamics = control["shrubDynamics"];
  String allocationStrategy = control["allocationStrategy"];
  if(transpirationMode=="Granier") allocationStrategy = "Al2As";
  String cavitationRefill = control["cavitationRefill"];
  String rhizosphereOverlap = control["rhizosphereOverlap"];
  bool plantWaterPools = (rhizosphereOverlap!="total");
  bool taper = control["taper"];
  double averageFracRhizosphereResistance = control["averageFracRhizosphereResistance"];
  double phloemConductanceFactor = control["phloemConductanceFactor"];
  double nonSugarConcentration = control["nonSugarConcentration"];
  List equilibriumOsmoticConcentration  = control["equilibriumOsmoticConcentration"];
  double equilibriumLeafTotalConc = equilibriumOsmoticConcentration["leaf"];
  double equilibriumSapwoodTotalConc = equilibriumOsmoticConcentration["sapwood"];
  double maximumStemConductance = control["maximumStemConductance"];
  
  //Soil-plant water balance
  List spwbOut;
  if(transpirationMode=="Granier") {
    spwbOut = spwbDay_basic(x, meteovec, 
                       elevation, slope, aspect,
                       runon, verbose); 
  } else {
    spwbOut = spwbDay_advanced(x, meteovec, 
                       latitude, elevation, slope, aspect,
                       solarConstant, delta, 
                       runon, verbose);
  }
  //Weather
  double tday = meteovec["tday"];
  double tmax = meteovec["tmax"];
  double Patm = meteovec["Patm"];
  double pfire = meteovec["pfire"];

  bool fireOccurrence = false;
  NumericVector fireBehavior(1,NA_REAL);
  if(R::runif(0.0,1.0) < pfire) {
    fireBehavior = fccsHazard(x, meteovec, spwbOut, slope);
    fireOccurrence = true;
  }
  
  //Cohort info
  DataFrame cohorts = Rcpp::as<Rcpp::DataFrame>(x["cohorts"]);
  IntegerVector SP = Rcpp::as<Rcpp::IntegerVector>(cohorts["SP"]);
  int numCohorts = SP.size();
  
  //Soil
  List soil  = x["soil"];
  NumericVector psiSoil = psi(soil,"soilFunctions");
  NumericVector dVec = soil["dVec"];
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
  if((transpirationMode=="Sperry") || (transpirationMode=="Cochard")) {
    RhizoPsi = Rcpp::as<Rcpp::NumericMatrix>(belowLayers["RhizoPsi"]);
    VCroot_kmax = Rcpp::as<Rcpp::NumericMatrix>(belowLayers["VCroot_kmax"]);
    VGrhizo_kmax = Rcpp::as<Rcpp::NumericMatrix>(belowLayers["VGrhizo_kmax"]);
  }
  
  int numLayers = V.ncol();
  
  //Internal state variables
  internalWater = Rcpp::as<Rcpp::DataFrame>(x["internalWater"]);

  //Values at the end of the day (after calling spwb)
  NumericVector StemPLC = Rcpp::as<Rcpp::NumericVector>(internalWater["StemPLC"]);
  NumericVector PlantPsi, psiApoLeaf, psiApoStem, psiSympLeaf, psiSympStem;
  if(transpirationMode=="Granier") {
    PlantPsi  = Rcpp::as<Rcpp::NumericVector>(internalWater["PlantPsi"]);
  } else {
    psiApoLeaf = Rcpp::as<Rcpp::NumericVector>(internalWater["LeafPsi"]);
    if(transpirationMode == "Sperry") {
      psiApoStem = Rcpp::as<Rcpp::NumericVector>(internalWater["Stem1Psi"]);
      psiSympLeaf = Rcpp::as<Rcpp::NumericVector>(internalWater["LeafSympPsi"]);
    } else {
      psiApoStem = Rcpp::as<Rcpp::NumericVector>(internalWater["StemPsi"]);
      psiSympLeaf = Rcpp::as<Rcpp::NumericVector>(internalWater["LeafPsi"]);
    }
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
  
  DataFrame Plants = Rcpp::as<Rcpp::DataFrame>(spwbOut["Plants"]);
  List PlantsInst;
  NumericVector Ag = Plants["GrossPhotosynthesis"];
  NumericVector LFMC = Plants["LFMC"];
  NumericVector PARcohort;
  NumericMatrix AgStep, AnStep;
  int numSteps = 1;
  if(transpirationMode=="Granier") {
    PARcohort= Plants["FPAR"];
  } else {
    PlantsInst = spwbOut["PlantsInst"];
    AgStep  =  Rcpp::as<Rcpp::NumericMatrix>(PlantsInst["Ag"]);
    AnStep  =  Rcpp::as<Rcpp::NumericMatrix>(PlantsInst["An"]);
    numSteps = AgStep.ncol();
  }

  //Data from spwb
  NumericVector Tcan;
  NumericMatrix StemSympPsiInst, LeafSympPsiInst, StemSympRWCInst, LeafSympRWCInst;
  List eb;
  double tcan_day = NA_REAL;
  if((transpirationMode=="Sperry") || (transpirationMode=="Cochard")) {
    StemSympPsiInst =  Rcpp::as<Rcpp::NumericMatrix>(PlantsInst["StemSympPsi"]);
    LeafSympPsiInst =  Rcpp::as<Rcpp::NumericMatrix>(PlantsInst["LeafSympPsi"]);
    StemSympRWCInst =  Rcpp::as<Rcpp::NumericMatrix>(PlantsInst["StemSympRWC"]);
    LeafSympRWCInst =  Rcpp::as<Rcpp::NumericMatrix>(PlantsInst["LeafSympRWC"]);

    eb = spwbOut["EnergyBalance"];  
    DataFrame tempDF =  Rcpp::as<Rcpp::DataFrame>(eb["Temperature"]);
    Tcan = Rcpp::as<Rcpp::NumericVector>(tempDF["Tcan"]);
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
    VCleaf_c = Rcpp::as<Rcpp::NumericVector>(paramsTransp["VCleaf_c"]);
    VCleaf_d = Rcpp::as<Rcpp::NumericVector>(paramsTransp["VCleaf_d"]);
    VCstem_kmax = paramsTransp["VCstem_kmax"];
    VCstem_c = Rcpp::as<Rcpp::NumericVector>(paramsTransp["VCstem_c"]);
    VCstem_d = Rcpp::as<Rcpp::NumericVector>(paramsTransp["VCstem_d"]);
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
  
  //Subdaily output matrices (Sperry/Cochard)
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
  NumericVector LeafBiomass(numCohorts,0.0), SapwoodBiomass(numCohorts, 0.0), SapwoodArea(numCohorts,0.0), LeafArea(numCohorts,0.0), FineRootArea(numCohorts, 0.0), HuberValue(numCohorts, 0.0), RootAreaLeafArea(numCohorts, 0.0);
  NumericVector SAgrowth(numCohorts,0.0), LAgrowth(numCohorts,0.0), FRAgrowth(numCohorts,0.0), starvationRate(numCohorts,0.0), dessicationRate(numCohorts,0.0), mortalityRate(numCohorts,0.0);
  NumericVector GrossPhotosynthesis(numCohorts,0.0);
  NumericVector RootExudation(numCohorts,0.0);

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
      NumericVector deltaFRBgrowth(numLayers, 0.0);
      
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
      NumericVector rfineroot(numLayers);
      if(transpirationMode=="Granier") {
        // grow_ring(ring, PlantPsi[j] ,tday, 10.0);
        rleafcell = std::min(rcellmax, relative_expansion_rate(PlantPsi[j] ,tday, -1.0, 0.5,0.05,5.0));
        rcambiumcell = std::min(rcellmax, relative_expansion_rate(PlantPsi[j] ,tday, -1.0, 0.5,0.05,5.0));
        for(int l=0;l<numLayers;l++) rfineroot[l] = std::min(rcellmax, relative_expansion_rate(psiSoil[l] ,tday, -1.0 ,0.5,0.05,5.0));
        // if(j==0) Rcout<<j<< " Psi:"<< PlantPsi[j]<< " r:"<< rcambiumcell<<"\n";
      } else {
        rleafcell = std::min(rcellmax, relative_expansion_rate(psiSympLeaf[j] ,tcan_day, -1.0, 0.5,0.05,5.0));
        rcambiumcell = std::min(rcellmax, relative_expansion_rate(psiSympStem[j] ,tcan_day, -1.0, 0.5,0.05,5.0));
        for(int l=0;l<numLayers;l++) rfineroot[l] = std::min(rcellmax, relative_expansion_rate(RhizoPsi(j,l) ,Tsoil[l], -1.0, 0.5,0.05,5.0));
        // if(j==0) Rcout<<j<< " Psi:"<< psiSympStem[j]<< " pi0:"<< " r:"<< rcambiumcell<<"\n";
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
          leafRespDay = B_resp_leaves*RERleaf[j]*QR*std::min(1.0, pow(PARcohort[j]/100.0,WUE_par[j]));
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
          for(int s = 0;s<numLayers;s++) {
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
          // Rcout<< SAring.size()<<" " <<j<< " "<< PlantPsi[j]<< " "<<" dSAring "<<deltaSAring<< " dSAsink "<< deltaSAsink<<" dSAgrowth "<< deltaSAgrowth<<" rgrcellfile"<< rgrcellfile<<"\n";
          growthCostSA = deltaSAgrowth[j]*costPerSA; //increase cost (may be non zero if leaf growth was charged onto sapwood)
          synthesisRespSA = growthCostSA*(CCsapwood[j] - 1.0)/CCsapwood[j];
        }
        
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
            for(int l = 0;l<numLayers;l++) {
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
        double LAplc = std::min(LAexpanded, (1.0-StemPLC[j])*leafAreaTarget[j]);
        if(LAplc<LAexpanded) {
          propLeafSenescence = std::max((LAexpanded-LAplc)/LAexpanded, propLeafSenescence);
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
      NumericVector deltaFRBsenescence(numLayers, 0.0);
      for(int l=0;l<numLayers;l++) {
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
      NumericVector newFRB(numLayers,0.0);
      for(int s=0;s<numLayers;s++) {
        newFRB[s] = fineRootBiomass[j]*V(j,s) + deltaFRBgrowth[s] - deltaFRBsenescence[s];
      }
      fineRootBiomass[j] = sum(newFRB);
      for(int s=0;s<numLayers;s++) { 
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
          L(c,_) = coarseRootLengths(V(c,_), dVec, 0.5); //Arbitrary ratio (to revise some day)
          CRSV[c] = coarseRootSoilVolume(V(c,_), dVec, 0.5);
        }
      } else { //SPERRY/COCHARD
        if(LAlive>0.0) {
          //Update Huber value, stem and root hydraulic conductance
          double oldstemR = 1.0/VCstem_kmax[j];
          double oldrootR = 1.0/VCroot_kmaxVEC[j];
          double oldrootprop = oldrootR/(oldrootR+oldstemR);
          
          // Al2As[j] = (LAlive)/(SA[j]/10000.0);
          VCstem_kmax[j]=maximumStemHydraulicConductance(Kmax_stemxylem[j], Hmed[j], Al2As[j] ,H[j], taper); 
          VCstem_kmax[j]=std::min(VCstem_kmax[j], maximumStemConductance);
          
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
        }
      }
      //Decrease PLC due to new SA growth
      if(cavitationRefill=="growth") StemPLC[j] = std::max(0.0, StemPLC[j] - (deltaSAgrowth[j]/SA[j]));
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
      else stemSympRWC = sum(StemSympRWCInst(j,_))/((double) numSteps);
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
            if(allowStarvation) starvationRate[j] = dailyMortalityProbability(relativeSugarSapwood, mortalityRelativeSugarThreshold);
            if(allowDessication) dessicationRate[j] = dailyMortalityProbability(std::max(stemSympRWC, 1.0 - StemPLC[j]), mortalityRWCThreshold);
            mortalityRate[j] = max(NumericVector::create(basalMortalityRate, dessicationRate[j],  starvationRate[j]));
            if((dessicationRate[j] > basalMortalityRate) && (dessicationRate[j] > starvationRate[j])) {
              cause = "dessication";
            } else if((starvationRate[j] > basalMortalityRate) && (starvationRate[j] > dessicationRate[j])) {
              cause = "starvation";
            }
            // Rcout<< j << " "<< stemSympRWC<< " "<< dessicationRate[j]<<"\n";
            if(mortalityMode =="density/deterministic") {
              Ndead_day = N[j]*mortalityRate[j];
            } else if(mortalityMode =="whole-cohort/stochastic") {
              if(R::runif(0.0,1.0) < mortalityRate[j]) {
                Ndead_day = N[j];
                Rcout<<" [Cohort "<< j<<" died from " << cause.get_cstring() << "] ";
              }
            } else if(mortalityMode == "density/stochastic") {
              Ndead_day = R::rbinom(round(N[j]), mortalityRate[j]);
            }
          }
        } else { //Cohort burned
          Ndead_day = N[j];
          cause = "burnt";
        }
        // Update density and increase the number of dead plants
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
          NumericVector VGrhizo_target(numLayers,0.0);
          for(int s=0;s<numLayers;s++) {
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
    PlantSugarLeaf[j] = sugarLeaf[j];
    PlantStarchLeaf[j] = starchLeaf[j];
    PlantSugarSapwood[j] = sugarSapwood[j];
    PlantStarchSapwood[j] = starchSapwood[j];
  }
  
  //CLOSE BIOMASS BALANCE
  closePlantBiomassBalance(plantBiomassBalance, x,
                      LabileCarbonBalance, LeafBiomassBalance, FineRootBiomassBalance);
  
  //Update pool proportions and rhizosphere overlap
  if(plantWaterPools) {
    NumericVector poolProportions = Rcpp::as<Rcpp::NumericVector>(belowdf["poolProportions"]);
    for(int j=0;j<numCohorts;j++) poolProportions[j] = LAI_live[j]/sum(LAI_live);
    //Update RHOP
    List newRHOP;
    if(rhizosphereOverlap=="none") newRHOP = nonoverlapHorizontalProportions(V);
    else newRHOP = horizontalProportions(poolProportions, CRSV, N, V, dVec, rfc);
    for(int j=0;j<numCohorts;j++) RHOP[j] = newRHOP[j];
  }
  

  List labileCBInst;
  if(subdailyCarbonBalance) {
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
    labileCBInst = List::create(
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
                                                   _["SugarTransport"] = PlantSugarTransport);
  labileCarbonBalance.attr("row.names") = above.attr("row.names");
  
  NumericVector standCB = {standGrossPrimaryProduction, standMaintenanceRespiration, standSynthesisRespiration,
                           standGrossPrimaryProduction - standMaintenanceRespiration - standSynthesisRespiration};
  standCB.attr("names") = CharacterVector({"GrossPrimaryProduction", "MaintenanceRespiration", "SynthesisRespiration", "NetPrimaryProduction"});
    
  //Final Biomass compartments
  DataFrame plantStructure = DataFrame::create(
    _["LeafBiomass"] = LeafBiomass,
    _["SapwoodBiomass"] = SapwoodBiomass,
    _["FineRootBiomass"] = clone(fineRootBiomass),
    _["LeafArea"] = LeafArea,
    _["SapwoodArea"] = SapwoodArea,
    _["FineRootArea"] = FineRootArea,
    _["HuberValue"] = HuberValue,
    _["RootAreaLeafArea"] = RootAreaLeafArea,
    _["DBH"] = clone(DBH),
    _["Height"] = clone(H)
  );
  
  DataFrame growthMortality = DataFrame::create(
    _["SAgrowth"] = SAgrowth,
    _["LAgrowth"] = LAgrowth,
    _["FRAgrowth"] = FRAgrowth,
    _["StarvationRate"] = starvationRate,
    _["DessicationRate"] = dessicationRate,
    _["MortalityRate"] = mortalityRate
  );
  growthMortality.attr("row.names") = above.attr("row.names");
  
  List l;
  if(transpirationMode=="Granier"){
       l = List::create(_["cohorts"] = clone(cohorts),
                        _["topography"] = spwbOut["topography"],
                        _["weather"] = spwbOut["weather"],
                        _["WaterBalance"] = spwbOut["WaterBalance"], 
                        _["CarbonBalance"] = standCB,
                        _["Soil"] = spwbOut["Soil"], 
                        _["Stand"] = spwbOut["Stand"], 
                        _["Plants"] = spwbOut["Plants"],
                        _["LabileCarbonBalance"] = labileCarbonBalance,
                        _["PlantBiomassBalance"] = plantBiomassBalance,
                        _["PlantStructure"] = plantStructure,
                        _["GrowthMortality"] = growthMortality);
  } else {
         l = List::create(_["cohorts"] = clone(cohorts),
                          _["topography"] = spwbOut["topography"],
                          _["weather"] = spwbOut["weather"],
                          _["WaterBalance"] = spwbOut["WaterBalance"], 
                          _["EnergyBalance"] = spwbOut["EnergyBalance"],
                          _["CarbonBalance"] = standCB,
                          _["Soil"] = spwbOut["Soil"], 
                          _["Stand"] = spwbOut["Stand"],
                          _["Plants"] = spwbOut["Plants"],
                          _["LabileCarbonBalance"] = labileCarbonBalance,
                          _["PlantBiomassBalance"] = plantBiomassBalance,
                          _["PlantStructure"] = plantStructure,
                          _["GrowthMortality"] = growthMortality,
                          _["RhizoPsi"] = spwbOut["RhizoPsi"],
                          _["SunlitLeaves"] = spwbOut["SunlitLeaves"],
                          _["ShadeLeaves"] = spwbOut["ShadeLeaves"],
                          _["ExtractionInst"] = spwbOut["ExtractionInst"],
                          _["PlantsInst"] = spwbOut["PlantsInst"],
                          _["SunlitLeavesInst"] = spwbOut["SunlitLeavesInst"],
                          _["ShadeLeavesInst"] = spwbOut["ShadeLeavesInst"]);
    if(subdailyCarbonBalance) l.push_back(labileCBInst,"LabileCarbonBalanceInst");
    l.push_back(spwbOut["LightExtinction"], "LightExtinction");
    l.push_back(spwbOut["CanopyTurbulence"], "CanopyTurbulence");
  }
  if(control["fireHazardResults"]) l.push_back(spwbOut["FireHazard"], "FireHazard");
  l.attr("class") = CharacterVector::create("growth_day","list");
  return(l);
}




//' @rdname spwb_day
// [[Rcpp::export("growth_day")]]
List growthDay(List x, CharacterVector date, NumericVector meteovec, 
               double latitude, double elevation, double slope, double aspect,  
               double runon=0.0, bool modifyInput = true) {
  
  double tmin = meteovec["MinTemperature"];
  double tmax = meteovec["MaxTemperature"];
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
  double rad = meteovec["Radiation"];
  double prec = meteovec["Precipitation"];
  double wind = NA_REAL;
  if(meteovec.containsElementNamed("WindSpeed")) wind = meteovec["WindSpeed"];
  double Catm = NA_REAL; 
  if(meteovec.containsElementNamed("CO2")) Catm = meteovec["CO2"];
  double Patm = NA_REAL; 
  if(meteovec.containsElementNamed("Patm")) Patm = meteovec["Patm"];
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
  
  //Update phenology
  if(leafPhenology) {
    updatePhenology(x, doy, photoperiod, tday);
    updateLeaves(x, wind, true);
  }
  
  double er = erFactor(doy, pet, prec);
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
    Named("er") = er,
    Named("pfire") = pfire);
  List s = growthDayInner(x, meteovec_inner, 
                     latitude, elevation, slope, aspect,
                     solarConstant, delta, 
                     runon, verbose);
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
  if((transpirationMode=="Sperry") || (transpirationMode=="Cochard")) {
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
  } else if((transpirationMode=="Sperry") || (transpirationMode=="Cochard")) {
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
//' 
//' @details
//' Detailed model description is available in the medfate book. 
//' Simulations using the 'Sperry' or 'Cochard' transpiration modes are computationally much more expensive 
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
//'   \item{\code{"EnergyBalance"}: A data frame with the daily values of energy balance components for the soil and the canopy (only for \code{transpirationMode = "Sperry"} or \code{transpirationMode = "Cochard"}; see \code{\link{spwb}}).}
//'   \item{\code{"CarbonBalance"}: A data frame where different stand-level carbon balance components (gross primary production, maintenance respiration, synthesis respiration and net primary production), all in g C · m-2.}
//'   \item{\code{"BiomassBalance"}: A data frame with the daily values of stand biomass balance components (in g dry · m-2.}
//'   \item{\code{"Temperature"}: A data frame with the daily values of minimum/mean/maximum temperatures for the atmosphere (input), canopy and soil (only for \code{transpirationMode = "Sperry"} or \code{transpirationMode = "Cochard"}; see \code{\link{spwb}}).}
//'   \item{\code{"Soil"}: A data frame where different soil variables  (see \code{\link{spwb}}).}
//'   \item{\code{"Stand"}: A data frame where different stand-level variables (see \code{\link{spwb}}).}
//'   \item{\code{"Plants"}: A list of daily results for plant cohorts (see \code{\link{spwb}}).}
//'   \item{\code{"SunlitLeaves"} and \code{"ShadeLeaves"}: A list with daily results for sunlit and shade leaves (only for \code{transpirationMode = "Sperry"} or \code{transpirationMode = "Cochard"}; see \code{\link{spwb}}).}
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
//' #Load example daily meteorological data
//' data(examplemeteo)
//' 
//' #Load example plot plant data
//' data(exampleforestMED)
//' 
//' #Default species parameterization
//' data(SpParamsMED)
//'   
//' #Initialize control parameters
//' control <- defaultControl("Granier")
//'   
//' #Initialize soil with default soil params (4 layers)
//' examplesoil <- soil(defaultSoilParams(4))
//' 
//' #Initialize vegetation input
//' x1 <- forest2growthInput(exampleforestMED, examplesoil, SpParamsMED, control)
//' 
//' #Call simulation function
//' G1 <- growth(x1, examplemeteo, latitude = 41.82592, elevation = 100)
//'  
//' \donttest{
//' #Switch to 'Sperry' transpiration mode
//' control <- defaultControl("Sperry")
//' 
//' #Initialize vegetation input
//' x2 <- forest2growthInput(exampleforestMED,examplesoil, SpParamsMED, control)
//' 
//' #Call simulation function
//' G2 <-growth(x2, examplemeteo, latitude = 41.82592, elevation = 100)
//' 
//' #Switch to 'Cochard' transpiration mode
//' control <- defaultControl("Cochard")
//' 
//' #Initialize vegetation input
//' x3 <- forest2growthInput(exampleforestMED,examplesoil, SpParamsMED, control)
//' 
//' #Call simulation function
//' G3 <-growth(x3, examplemeteo, latitude = 41.82592, elevation = 100)
//' }
//'       
// [[Rcpp::export("growth")]]
List growth(List x, DataFrame meteo, double latitude, 
            double elevation = NA_REAL, double slope = NA_REAL, double aspect = NA_REAL,
            NumericVector CO2ByYear = NumericVector(0)) {

  //Control params 
  List control =x["control"];  
  String transpirationMode = control["transpirationMode"];
  String soilFunctions = control["soilFunctions"];
  bool verbose = control["verbose"];
  bool leafPhenology = control["leafPhenology"];
  bool unlimitedSoilWater = control["unlimitedSoilWater"];
  bool multiLayerBalance = control["multiLayerBalance"];
  checkgrowthInput(x, transpirationMode, soilFunctions);

  //Store input
  List growthInput = x; // Store initial object
  x = clone(x); //Ensure a copy will be modified

  //Soil params 
  List soil = x["soil"];
  
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
  DataFrame DEB = defineEnergyBalanceDailyOutput(dateStrings);
  DataFrame DT = defineTemperatureDailyOutput(dateStrings);
  NumericMatrix DLT;
  if((transpirationMode=="Sperry") || (transpirationMode == "Cochard")) DLT =  defineTemperatureLayersDailyOutput(dateStrings, canopy);
  
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
  
  //Water balance output variables
  DataFrame DWB = defineWaterBalanceDailyOutput(dateStrings, PET, transpirationMode);
  DataFrame SWB = defineSoilWaterBalanceDailyOutput(dateStrings, soil, transpirationMode);
  
  
  NumericVector LAI(numDays), LAIherb(numDays) , LAIlive(numDays), LAIexpanded(numDays), LAIdead(numDays);
  NumericVector Cm(numDays);
  NumericVector LgroundPAR(numDays);
  NumericVector LgroundSWR(numDays);

  //Plant water output variables
  List sunlitDO = defineSunlitShadeLeavesDailyOutput(dateStrings, above);
  List shadeDO = defineSunlitShadeLeavesDailyOutput(dateStrings, above);
  List plantDWOL = definePlantWaterDailyOutput(dateStrings, above, soil, control);
  NumericVector EplantCohTot(numCohorts, 0.0);

  //Fire hazard output variables
  DataFrame fireHazard;
  if(control["fireHazardResults"]) fireHazard = defineFireHazardOutput(dateStrings);
  
  //Count years (times structural variables will be updated)
  int numYears = 0;
  for(int i=0;i<numDays;i++) {
    if(((DOY[i]==1) && (i>0)) || ((i==(numDays-1)) && (DOY[i]>=365))) numYears = numYears + 1;
  }

  NumericVector initialSoilContent = water(soil, soilFunctions);
  NumericVector initialPlantContent = plantWaterContent(x);
  double initialSnowContent = soil["SWE"];
  DataFrame ccIni_m2 = carbonCompartments(x, "g_m2");
  double cohortBiomassBalanceSum = 0.0;
  double initialCohortBiomass = sum(Rcpp::as<Rcpp::NumericVector>(ccIni_m2["TotalBiomass"]));
  
  if(verbose) {
    Rcout<<"Initial plant cohort biomass (g/m2): "<<initialCohortBiomass<<"\n";
    Rcout<<"Initial plant water content (mm): "<< sum(initialPlantContent)<<"\n";
    Rcout<<"Initial soil water content (mm): "<< sum(initialSoilContent)<<"\n";
    Rcout<<"Initial snowpack content (mm): "<< initialSnowContent<<"\n";
  }
  
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
      else if(i%10 == 0) Rcout<<".";//<<i;
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
      meteovec_inner.push_back(erFactor(DOY[i], PET[i], prec), "er");
      meteovec_inner.push_back(FireProbability[i], "pfire"); 
      try{
        s = growthDayInner(x, meteovec_inner,  
                           latitude, elevation, slope, aspect,
                           solarConstant, delta, 
                           0.0, false); //No Runon in simulations for a single cell
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
        Named("er") = erFactor(DOY[i], PET[i], Precipitation[i]));
      meteovec.push_back(FireProbability[i], "pfire"); 
      try{
        s = growthDayInner(x, meteovec, 
                           latitude, elevation, slope, aspect,
                           solarConstant, delta, 
                           0.0, verbose);
      } catch(std::exception& ex) {
        Rcerr<< "c++ error: "<< ex.what() <<"\n";
        error_occurence = true;
      }
      fillEnergyBalanceTemperatureDailyOutput(DEB,DT,DLT,s,i, multiLayerBalance);
    }    
    
    fillPlantWaterDailyOutput(plantDWOL, sunlitDO, shadeDO, s, i, transpirationMode);
    fillWaterBalanceDailyOutput(DWB, s,i, transpirationMode);
    fillSoilWaterBalanceDailyOutput(SWB, soil, s,
                                    i, numDays, transpirationMode, soilFunctions);
    if(control["fireHazardResults"]) fillFireHazardOutput(fireHazard, s, i);
    
    List stand = s["Stand"];
    LgroundPAR[i] = stand["LgroundPAR"];
    LgroundSWR[i] = stand["LgroundSWR"];
    LAI[i] = stand["LAI"];
    LAIherb[i] = stand["LAIherb"];
    LAIlive[i] = stand["LAIlive"];
    LAIexpanded[i] = stand["LAIexpanded"];
    LAIdead[i] = stand["LAIdead"];
    Cm[i] = stand["Cm"];
    
    List sb = s["Soil"];
    List db = s["WaterBalance"];
    List Plants = s["Plants"];
    DataFrame cb = Rcpp::as<Rcpp::DataFrame>(s["LabileCarbonBalance"]);
    DataFrame bb = Rcpp::as<Rcpp::DataFrame>(s["PlantBiomassBalance"]);
    DataFrame ps = Rcpp::as<Rcpp::DataFrame>(s["PlantStructure"]);
    DataFrame gm = Rcpp::as<Rcpp::DataFrame>(s["GrowthMortality"]);
    
    
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
    RootExudation(i,_) = Rcpp::as<Rcpp::NumericVector>(cb["RootExudation"]);
    
    SapwoodBiomass(i,_) = Rcpp::as<Rcpp::NumericVector>(ps["SapwoodBiomass"]);
    LeafBiomass(i,_) = Rcpp::as<Rcpp::NumericVector>(ps["LeafBiomass"]);
    FineRootBiomass(i,_) = Rcpp::as<Rcpp::NumericVector>(ps["FineRootBiomass"]);
    SapwoodArea(i,_) = Rcpp::as<Rcpp::NumericVector>(ps["SapwoodArea"]);
    LeafArea(i,_) = Rcpp::as<Rcpp::NumericVector>(ps["LeafArea"]);
    FineRootArea(i,_) = Rcpp::as<Rcpp::NumericVector>(ps["FineRootArea"]);
    HuberValue(i,_) = Rcpp::as<Rcpp::NumericVector>(ps["HuberValue"]);
    RootAreaLeafArea(i,_) = Rcpp::as<Rcpp::NumericVector>(ps["RootAreaLeafArea"]);
    DBH(i,_) = Rcpp::as<Rcpp::NumericVector>(ps["DBH"]);
    Height(i,_) = Rcpp::as<Rcpp::NumericVector>(ps["Height"]);
    
    StructuralBiomassBalance(i,_) = Rcpp::as<Rcpp::NumericVector>(bb["StructuralBiomassBalance"]);
    LabileBiomassBalance(i,_) = Rcpp::as<Rcpp::NumericVector>(bb["LabileBiomassBalance"]);
    PlantBiomassBalance(i,_) = Rcpp::as<Rcpp::NumericVector>(bb["PlantBiomassBalance"]);
    MortalityBiomassLoss(i,_) = Rcpp::as<Rcpp::NumericVector>(bb["MortalityBiomassLoss"]);
    CohortBiomassBalance(i,_) = Rcpp::as<Rcpp::NumericVector>(bb["CohortBiomassBalance"]);
    StandBiomassBalance(i,_) = standLevelBiomassBalance(bb);
    cohortBiomassBalanceSum += sum(CohortBiomassBalance(i,_));

    LAgrowth(i,_) = Rcpp::as<Rcpp::NumericVector>(gm["LAgrowth"]);
    SAgrowth(i,_) = Rcpp::as<Rcpp::NumericVector>(gm["SAgrowth"]);
    FRAgrowth(i,_) = Rcpp::as<Rcpp::NumericVector>(gm["FRAgrowth"]);
    
    starvationRate(i,_) = Rcpp::as<Rcpp::NumericVector>(gm["StarvationRate"]);
    dessicationRate(i,_) = Rcpp::as<Rcpp::NumericVector>(gm["DessicationRate"]);
    mortalityRate(i,_) = Rcpp::as<Rcpp::NumericVector>(gm["MortalityRate"]);
    
    StandCarbonBalance(i,_) = Rcpp::as<Rcpp::NumericVector>(s["CarbonBalance"]);
    
    //5 Update structural variables
    // if(((i<(numDays-1)) && (DOY[i+1]==1)) || (i==(numDays-1))) { 
    //   //reset ring structures
    //   List ringList = x["internalRings"];
    //   for(int j=0;j<numCohorts; j++) ringList[j] = initialize_ring();
    // }

    if(control["subdailyResults"]) {
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
    
    printWaterBalanceResult(DWB, plantDWOL, x,
                            initialPlantContent, initialSoilContent, initialSnowContent,
                            transpirationMode);
    
    if(error_occurence) {
      Rcout<< " ERROR: Calculations stopped because of numerical error: Revise parameters\n";
    }
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
  SapwoodBiomass.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names")) ;
  LeafBiomass.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names")) ;
  FineRootArea.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names")) ;
  SapwoodArea.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names")) ;
  LeafArea.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names")) ;
  HuberValue.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names")) ;
  RootAreaLeafArea.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names")) ;
  FineRootBiomass.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names"));
  DBH.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names"));
  Height.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names"));
  LAgrowth.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names")) ;
  SAgrowth.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names")) ;
  FRAgrowth.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names"));
  dessicationRate.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names")) ;
  mortalityRate.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names")) ;
  starvationRate.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names"));
  RootExudation.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names"));
  StructuralBiomassBalance.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names"));
  LabileBiomassBalance.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names"));
  PlantBiomassBalance.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names"));
  MortalityBiomassLoss.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names"));
  CohortBiomassBalance.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names"));

  StandBiomassBalance.attr("dimnames") = List::create(meteo.attr("row.names"), 
                           CharacterVector::create("StructuralBalance", "LabileBalance", "PlantBalance", "MortalityLoss", "CohortBalance"));
  StandCarbonBalance.attr("dimnames") = List::create(meteo.attr("row.names"), 
                          CharacterVector::create("GrossPrimaryProduction", "MaintenanceRespiration", "SynthesisRespiration", "NetPrimaryProduction"));
  subdailyRes.attr("names") = meteo.attr("row.names") ;
  
  NumericVector topo = NumericVector::create(elevation, slope, aspect);
  topo.attr("names") = CharacterVector::create("elevation", "slope", "aspect");
  
  Rcpp::DataFrame Stand = DataFrame::create(_["LAI"]=LAI, _["LAIherb"]=LAIherb, 
                                            _["LAIlive"]=LAIlive, _["LAIexpanded"]=LAIexpanded,_["LAIdead"]=LAIdead,
                                            _["Cm"]=Cm, _["LgroundPAR"] = LgroundPAR, _["LgroundSWR"] = LgroundSWR);
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
    Named("SugarTransport") = PlantSugarTransport
  );
  List plantBiomassBalance = List::create(_["StructuralBiomassBalance"] = StructuralBiomassBalance,
                                     _["LabileBiomassBalance"] = LabileBiomassBalance,
                                     _["PlantBiomassBalance"] = PlantBiomassBalance,
                                     _["MortalityBiomassLoss"] = MortalityBiomassLoss,
                                     _["CohortBiomassBalance"] = CohortBiomassBalance);
  
  List growthMortality, plantStructure;
  
  List l;
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
  if(transpirationMode=="Granier") {
    l = List::create(Named("latitude") = latitude,
                     Named("topography") = topo,
                     Named("weather") = clone(meteo),
                     Named("growthInput") = growthInput,
                     Named("growthOutput") = clone(x),
                     Named("WaterBalance")=DWB, 
                     Named("CarbonBalance")=StandCarbonBalance, 
                     Named("BiomassBalance") = StandBiomassBalance);
    if(control["soilResults"]) l.push_back(SWB, "Soil");
    if(control["standResults"]) l.push_back(Stand, "Stand");
    if(control["plantResults"]) {
      l.push_back(plantDWOL, "Plants");
      l.push_back(labileCarbonBalance, "LabileCarbonBalance");
      l.push_back(plantBiomassBalance, "PlantBiomassBalance");
      l.push_back(plantStructure, "PlantStructure");
      l.push_back(growthMortality, "GrowthMortality");
    }
    if(control["fireHazardResults"]) l.push_back(fireHazard, "FireHazard");
  } else {
    l = List::create(Named("latitude") = latitude,
                   Named("topography") = topo,
                   Named("weather") = clone(meteo),
                   Named("growthInput") = growthInput,
                   Named("growthOutput") = clone(x),
                   Named("WaterBalance")=DWB, 
                   Named("CarbonBalance")=StandCarbonBalance, 
                   Named("EnergyBalance") = DEB,
                   Named("BiomassBalance") = StandBiomassBalance);
    if(control["temperatureResults"]) {
      l.push_back(DT, "Temperature");
      if(multiLayerBalance) l.push_back(DLT,"TemperatureLayers");
    }
    if(control["soilResults"]) l.push_back(SWB, "Soil");
    if(control["standResults"]) l.push_back(Stand, "Stand");
    if(control["plantResults"]) l.push_back(plantDWOL, "Plants");
    if(control["leafResults"]) {
      l.push_back(sunlitDO, "SunlitLeaves");
      l.push_back(shadeDO, "ShadeLeaves");
    }
    if(control["plantResults"]){
      l.push_back(labileCarbonBalance, "LabileCarbonBalance");
      l.push_back(plantBiomassBalance, "PlantBiomassBalance");
      l.push_back(plantStructure, "PlantStructure");
      l.push_back(growthMortality, "GrowthMortality");
    }
    if(control["fireHazardResults"]) l.push_back(fireHazard, "FireHazard");
  }
  if(control["subdailyResults"]) l.push_back(subdailyRes,"subdaily");
  l.attr("class") = CharacterVector::create("growth","list");
  return(l);
}
