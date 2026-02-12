#include <RcppArmadillo.h>
#include "biophysicsutils_c.h"
#include "carbon_c.h"
#include "modelInput_c.h"
#include "phenology_c.h"
#include "spwb_day_c.h"
#include "growth_day_c.h"
#include "forestutils_c.h"
#include "lowlevel_structures_c.h"
#include "hydraulics_c.h"
#include "root_c.h"
#include "woodformation.h"
#include "meteoland/pet_c.hpp"

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
double dailyMortalityProbability_c(double stressValue, double stressThreshold) {
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
double phloemFlow_c(double psiUpstream, double psiDownstream,
                    double concUpstream, double concDownstream,
                    double temp, double k_f, double nonSugarConc) {
  double turgor_up = turgor_c(psiUpstream, concUpstream, temp, nonSugarConc);
  double turgor_down = turgor_c(psiDownstream, concDownstream, temp, nonSugarConc);
  if(temp < 0.0) k_f = 0.0; // No phloem flow if temperature below zero
  double relVisc = relativeSapViscosity_c((concUpstream+concDownstream)/2.0, temp);
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
//       double psiUp = symplasticWaterPotential_c(rwcLeaf(c,s), leafPI0[c], leafEPS[c]);
//       double psiDown = symplasticWaterPotential_c(rwcStem(c,s), stemPI0[c], stemEPS[c]);
//       ff(c,s) = phloemFlow(psiUp, psiDown, concLeaf[c], concSapwood[c], Tcan[s], k_phloem, nonSugarConc)*3600.0; //flow as mol per hour and leaf area basis
//     }
//   }
//   ff.attr("dimnames") = rwcStem.attr("dimnames");
//   return(ff);
// }

double qResp_c(double Tmean) {
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
Rcpp::DataFrame copyLabileCarbonBalanceResult_c(const LabileCarbonBalance_RESULT& LCBres, ModelInput& x) {
  Rcpp::DataFrame DF = Rcpp::DataFrame::create(
    Rcpp::Named("GrossPhotosynthesis") = Rcpp::wrap(LCBres.GrossPhotosynthesis),
    Rcpp::Named("MaintenanceRespiration") = Rcpp::wrap(LCBres.MaintenanceRespiration),
    Rcpp::Named("GrowthCosts") = Rcpp::wrap(LCBres.GrowthCosts),
    Rcpp::Named("RootExudation") = Rcpp::wrap(LCBres.RootExudation),
    Rcpp::Named("LabileCarbonBalance") = Rcpp::wrap(LCBres.LabileCarbonBalance),
    Rcpp::Named("SugarLeaf") = Rcpp::wrap(LCBres.SugarLeaf),
    Rcpp::Named("StarchLeaf") = Rcpp::wrap(LCBres.StarchLeaf),
    Rcpp::Named("SugarSapwood") = Rcpp::wrap(LCBres.SugarSapwood),
    Rcpp::Named("StarchSapwood") = Rcpp::wrap(LCBres.StarchSapwood),
    Rcpp::Named("SugarTransport") = Rcpp::wrap(LCBres.SugarTransport)
  );
  DF.attr("row.names") = x.cohorts.CohortCode;
  return(DF);
}
Rcpp::DataFrame copyPlantBiomassBalanceResult_c(const PlantBiomassBalance_RESULT& PBBres, ModelInput& x) {
  Rcpp::DataFrame DF = Rcpp::DataFrame::create(
    Rcpp::Named("InitialDensity") = Rcpp::wrap(PBBres.InitialDensity),
    Rcpp::Named("InitialSapwoodBiomass") = Rcpp::wrap(PBBres.InitialSapwoodBiomass),
    Rcpp::Named("InitialStructuralBiomass") = Rcpp::wrap(PBBres.InitialStructuralBiomass),
    Rcpp::Named("StructuralBiomassBalance") = Rcpp::wrap(PBBres.StructuralBiomassBalance),
    Rcpp::Named("StructuralBiomassChange") = Rcpp::wrap(PBBres.StructuralBiomassChange),
    Rcpp::Named("InitialLabileBiomass") = Rcpp::wrap(PBBres.InitialLabileBiomass),
    Rcpp::Named("LabileBiomassBalance") = Rcpp::wrap(PBBres.LabileBiomassBalance),
    Rcpp::Named("LabileBiomassChange") = Rcpp::wrap(PBBres.LabileBiomassChange),
    Rcpp::Named("InitialLivingPlantBiomass") = Rcpp::wrap(PBBres.InitialLivingPlantBiomass),
    Rcpp::Named("InitialPlantBiomass") = Rcpp::wrap(PBBres.InitialPlantBiomass),
    Rcpp::Named("PlantBiomassBalance") = Rcpp::wrap(PBBres.PlantBiomassBalance),
    Rcpp::Named("PlantBiomassChange") = Rcpp::wrap(PBBres.PlantBiomassChange),
    Rcpp::Named("MortalityBiomassLoss") = Rcpp::wrap(PBBres.MortalityBiomassLoss),
    Rcpp::Named("InitialCohortBiomass") = Rcpp::wrap(PBBres.InitialCohortBiomass),
    Rcpp::Named("CohortBiomassBalance") = Rcpp::wrap(PBBres.CohortBiomassBalance),
    Rcpp::Named("CohortBiomassChange") = Rcpp::wrap(PBBres.CohortBiomassChange)
  );
  DF.attr("row.names") = x.cohorts.CohortCode;
  return(DF);
}
Rcpp::DataFrame copyPlantStructureResult_c(const PlantStructure_RESULT& PSres, ModelInput& x) {
  Rcpp::DataFrame DF = Rcpp::DataFrame::create(
    Rcpp::Named("LeafBiomass") = Rcpp::wrap(PSres.LeafBiomass),
    Rcpp::Named("SapwoodBiomass") = Rcpp::wrap(PSres.SapwoodBiomass),
    Rcpp::Named("FineRootBiomass") = Rcpp::wrap(PSres.FineRootBiomass),
    Rcpp::Named("LeafArea") = Rcpp::wrap(PSres.LeafArea),
    Rcpp::Named("SapwoodArea") = Rcpp::wrap(PSres.SapwoodArea),
    Rcpp::Named("FineRootArea") = Rcpp::wrap(PSres.FineRootArea),
    Rcpp::Named("HuberValue") = Rcpp::wrap(PSres.HuberValue),
    Rcpp::Named("RootAreaLeafArea") = Rcpp::wrap(PSres.RootAreaLeafArea),
    Rcpp::Named("DBH") = Rcpp::wrap(PSres.DBH),
    Rcpp::Named("Height") = Rcpp::wrap(PSres.Height)
  );
  DF.attr("row.names") = x.cohorts.CohortCode;
  return(DF);
}
Rcpp::DataFrame copyGrowthMortalityResult_c(const GrowthMortality_RESULT& GMres, ModelInput& x) {
  Rcpp::DataFrame DF = Rcpp::DataFrame::create(
    Rcpp::Named("SAgrowth") = Rcpp::wrap(GMres.SAgrowth),
    Rcpp::Named("LAgrowth") = Rcpp::wrap(GMres.LAgrowth),
    Rcpp::Named("FRAgrowth") = Rcpp::wrap(GMres.FRAgrowth),
    Rcpp::Named("StarvationRate") = Rcpp::wrap(GMres.StarvationRate),
    Rcpp::Named("DessicationRate") = Rcpp::wrap(GMres.DessicationRate),
    Rcpp::Named("MortalityRate") = Rcpp::wrap(GMres.MortalityRate)
  );
  DF.attr("row.names") = x.cohorts.CohortCode;
  return(DF);
}
Rcpp::List copyBasicGROWTHResult_c(BasicGROWTH_RESULT& GROWTHres, ModelInput& x) {
  Rcpp::List l = copyBasicSPWBResult_c(GROWTHres.BSPWBres, x);
  l.push_back(copyCarbonBalanceResult_c(GROWTHres.standCB), "CarbonBalance");
  l.push_back(copyLabileCarbonBalanceResult_c(GROWTHres.LCBres, x), "LabileCarbonBalance");
  l.push_back(copyPlantBiomassBalanceResult_c(GROWTHres.PBBres, x), "PlantBiomassBalance");
  l.push_back(copyPlantStructureResult_c(GROWTHres.PSres, x), "PlantStructure");
  l.push_back(copyGrowthMortalityResult_c(GROWTHres.GMres, x), "GrowthMortality");
  l.attr("class") = Rcpp::CharacterVector::create("growth_day","list");
  return(l);
}

Rcpp::List copyAdvancedGROWTHResult_c(AdvancedGROWTH_RESULT& GROWTHres, ModelInput& x) {
  Rcpp::List l = copyAdvancedSPWBResult_c(GROWTHres.ASPWBres, x);
  l.push_back(copyCarbonBalanceResult_c(GROWTHres.standCB), "CarbonBalance");
  l.push_back(copyLabileCarbonBalanceResult_c(GROWTHres.LCBres, x), "LabileCarbonBalance");
  l.push_back(copyPlantBiomassBalanceResult_c(GROWTHres.PBBres, x), "PlantBiomassBalance");
  l.push_back(copyPlantStructureResult_c(GROWTHres.PSres, x), "PlantStructure");
  l.push_back(copyGrowthMortalityResult_c(GROWTHres.GMres, x), "GrowthMortality");
  l.attr("class") = Rcpp::CharacterVector::create("growth_day","list");
  return(l);
}

Rcpp::List copyGROWTHResult_c(GROWTH_RESULT& GROWTHres, ModelInput& x) {
  Rcpp::List l;
  if(x.control.transpirationMode=="Granier") {
    try {
      auto& BGROWTHres = dynamic_cast<BasicGROWTH_RESULT&>(GROWTHres);
      l = copyBasicGROWTHResult_c(BGROWTHres, x);
    } catch(const std::bad_cast&) {
      throw medfate::MedfateInternalError("Control transpiration mode set to basic but result object is not basic");
    }
  } else {
    try {
      auto& AGROWTHres = dynamic_cast<AdvancedGROWTH_RESULT&>(GROWTHres);
      l = copyAdvancedGROWTHResult_c(AGROWTHres, x);
    } catch(const std::bad_cast&) {
      throw medfate::MedfateInternalError("Control transpiration mode set to advanced but result object is not advanced");
    }
  }
  return(l);
}
void growthDay_private_c(GROWTH_RESULT& GROWTHres, GROWTHCommunicationStructures& GROWTHcomm, ModelInput& x, 
                         WeatherInputVector meteovec, 
                         double latitude, double elevation, double slope, double aspect,
                         double solarConstant, double delta, 
                         const double runon, 
                         const std::vector<double>& lateralFlows, const double waterTableDepth) {


  int ntimesteps = x.control.advancedWB.ndailysteps;
  int nlayers = x.soil.getNlayers();
  int numCohorts = x.internalWater.StemPLC.size();
  
  
  ///////////////////////////////////////////////////////////////////////////////////
  ///// A. WATER-ENERGY BALANCE (this creates communication structures as well) /////
  ///////////////////////////////////////////////////////////////////////////////////
  std::vector<double> Tcan(ntimesteps);
  std::vector<double> Ag(numCohorts);
  std::vector<double> PARcohort(numCohorts);
  
  double tcan_day = medfate::NA_DOUBLE;

  if(x.control.transpirationMode=="Granier") {
    try {
      auto& BGROWTHres = dynamic_cast<BasicGROWTH_RESULT&>(GROWTHres);
      spwbDay_basic_c(BGROWTHres.BSPWBres, GROWTHcomm.SPWBcomm.BSPWBcomm, x, 
                      meteovec,
                      elevation, slope, aspect,
                      runon, 
                      lateralFlows, waterTableDepth);
      
      for(int i=0;i<numCohorts;i++) {
        Ag[i] = BGROWTHres.BSPWBres.BTres.plants.GrossPhotosynthesis[i];
        PARcohort[i] = BGROWTHres.BSPWBres.BTres.plants.FPAR[i];
      }
    } catch(const std::bad_cast&) {
      throw medfate::MedfateInternalError("Control transpiration mode set to basic but result object is not basic");
    }
  } else {
    try {
      auto& AGROWTHres = dynamic_cast<AdvancedGROWTH_RESULT&>(GROWTHres);
      spwbDay_advanced_c(AGROWTHres.ASPWBres, GROWTHcomm.SPWBcomm.ASPWBcomm, x, 
                         meteovec,
                         latitude, elevation, slope, aspect,
                         solarConstant, delta,
                         runon, 
                         lateralFlows, waterTableDepth);
      for(int n=0;n<ntimesteps;n++) Tcan[n] = AGROWTHres.ASPWBres.ATres.energy.Tcan[n];
      tcan_day = averageDaylightTemperature_c(*std::min_element(Tcan.begin(), Tcan.end()), 
                                              *std::max_element(Tcan.begin(), Tcan.end()));
      for(int i=0;i<numCohorts;i++) {
        Ag[i] = AGROWTHres.ASPWBres.ATres.plants.GrossPhotosynthesis[i];
        PARcohort[i] = AGROWTHres.ASPWBres.ATres.plants.FPAR[i];
      }
      
    } catch(const std::bad_cast&) {
      throw medfate::MedfateInternalError("Control transpiration mode set to advanced but result object is not advanced");
    }
  }
  

  // String soilFunctions = control["soilFunctions"];
  // String mortalityMode = control["mortalityMode"];
  bool subdailyCarbonBalance = x.control.growth.subdailyCarbonBalance;
  // double mortalityRelativeSugarThreshold= control["mortalityRelativeSugarThreshold"];
  // double mortalityRWCThreshold= control["mortalityRWCThreshold"];
  // bool allowDessication = control["allowDessication"];
  // bool allowStarvation = control["allowStarvation"];
  bool sinkLimitation = x.control.growth.sinkLimitation;
  // bool shrubDynamics = control["shrubDynamics"];
  std::string allocationStrategy = x.control.growth.allocationStrategy;
  if(x.control.transpirationMode=="Granier")  {
    allocationStrategy = "Al2As";
    subdailyCarbonBalance = false;
  }
  // String rhizosphereOverlap = control["rhizosphereOverlap"];
  // bool plantWaterPools = (rhizosphereOverlap!="total");
  // bool taper = control["taper"];
  // double averageFracRhizosphereResistance = control["averageFracRhizosphereResistance"];
  // double phloemConductanceFactor = control["phloemConductanceFactor"];
  double nonSugarConcentration = x.control.growth.nonSugarConcentration;
  double equilibriumLeafTotalConc = x.control.growth.equilibriumOsmoticConcentration.leaf;
  double equilibriumSapwoodTotalConc = x.control.growth.equilibriumOsmoticConcentration.sapwood;
  // 
  //Weather
  double tday = meteovec.tday;
  double tmin = meteovec.tmin;
  double tmax = meteovec.tmax;
  double rhmin = meteovec.rhmin;
  double rhmax = meteovec.rhmax;
  double Patm = meteovec.Patm;

  //Atmospheric pressure (if missing)
  if(std::isnan(Patm)) Patm = atmosphericPressure_c(elevation);

  
  //Get previous PLC so that defoliation occurs only when PLC increases
  std::vector<double> StemPLCprev(numCohorts);
  for(int i=0; i< numCohorts;i++) StemPLCprev[i] = x.internalWater.StemPLC[i];
  std::vector<std::string> ctype = cohortType_c(x.cohorts.CohortCode);
  
  // //Cohort info
  // DataFrame cohorts = Rcpp::as<Rcpp::DataFrame>(x["cohorts"]);
  // CharacterVector speciesNames = cohorts["Name"];
  // CharacterVector ctype = cohortType(cohorts.attr("row.names"));
  // IntegerVector SP = Rcpp::as<Rcpp::IntegerVector>(cohorts["SP"]);
  // 
  // //Soil
  // DataFrame soil = Rcpp::as<Rcpp::DataFrame>(x["soil"]);
  // int nlayers = soil.nrow();
  // NumericVector psiSoil = psi(soil, soilFunctions);
  // NumericVector widths = soil["widths"];
  // NumericVector sand = soil["sand"];
  // NumericVector clay = soil["clay"];
  // NumericVector rfc = soil["rfc"];
  // NumericVector Ksat = soil["Ksat"];
  // NumericVector VG_n = soil["VG_n"];
  // NumericVector VG_alpha = soil["VG_alpha"];
  // NumericVector Tsoil = soil["Temp"];

  //Aboveground parameters
  std::vector<double>& DBH = x.above.DBH;
  std::vector<double>& Cover = x.above.Cover;
  std::vector<double>& H = x.above.H;
  std::vector<double>& N = x.above.N;
  std::vector<double>& CR = x.above.CR;
  std::vector<double>& LAI_live = x.above.LAI_live;
  std::vector<double>& LAI_expanded = x.above.LAI_expanded;
  std::vector<double>& LAI_dead = x.above.LAI_dead;
  std::vector<double>& LAI_nocomp = x.above.LAI_nocomp;
  std::vector<double>& SA = x.above.SA;

  // //Belowground parameters  
  // DataFrame belowdf = Rcpp::as<Rcpp::DataFrame>(x["below"]);
  // NumericVector Z95 = Rcpp::as<Rcpp::NumericVector>(belowdf["Z95"]);.
  // NumericVector Z50 = Rcpp::as<Rcpp::NumericVector>(belowdf["Z50"]);
  std::vector<double>& fineRootBiomass = x.below.fineRootBiomass;
  // NumericVector CRSV = Rcpp::as<Rcpp::NumericVector>(belowdf["coarseRootSoilVolume"]);
  // List belowLayers = Rcpp::as<Rcpp::List>(x["belowLayers"]);
  // List RHOP;
  // if(plantWaterPools) RHOP = belowLayers["RHOP"];
  // NumericMatrix V = Rcpp::as<Rcpp::NumericMatrix>(belowLayers["V"]);
  // NumericMatrix L = Rcpp::as<Rcpp::NumericMatrix>(belowLayers["L"]);
  // NumericMatrix RhizoPsi, VCroot_kmax, VGrhizo_kmax;
  // if((transpirationMode=="Sperry") || (transpirationMode=="Sureau")) {
  //   RhizoPsi = Rcpp::as<Rcpp::NumericMatrix>(belowLayers["RhizoPsi"]);
  //   VCroot_kmax = Rcpp::as<Rcpp::NumericMatrix>(belowLayers["VCroot_kmax"]);
  //   VGrhizo_kmax = Rcpp::as<Rcpp::NumericMatrix>(belowLayers["VGrhizo_kmax"]);
  // }
  // 
  // //Internal state variables
  // internalWater = Rcpp::as<Rcpp::DataFrame>(x["internalWater"]);
  // 
  
  //Values at the end of the day (after calling spwb)
  std::vector<double>& LeafPLC = x.internalWater.LeafPLC;
  std::vector<double>& StemPLC = x.internalWater.StemPLC;
  std::vector<double> psiApoLeaf(numCohorts), psiApoStem(numCohorts), psiSympLeaf(numCohorts), psiSympStem(numCohorts), psiRootCrown(numCohorts);
  for(int i=0;i<numCohorts;i++) {
    if(x.control.transpirationMode=="Granier") {
      psiApoLeaf[i] = x.internalWater.PlantPsi[i];
    } else {
      psiApoLeaf[i] = x.internalWater.LeafPsi[i];
      if(x.control.transpirationMode == "Sperry") {
        psiApoStem[i] = x.internalWater.StemPsi[i];
        psiSympLeaf[i] = x.internalWater.LeafSympPsi[i];
      } else {
        psiApoStem[i] = x.internalWater.StemPsi[i];
        psiSympLeaf[i] = x.internalWater.LeafPsi[i];
      }
      psiRootCrown[i] = x.internalWater.RootCrownPsi[i];
      psiSympStem[i] = x.internalWater.StemSympPsi[i];
    }
  }
  
  std::vector<double>& sugarLeaf = x.internalCarbon.sugarLeaf; //Concentrations assuming RWC = 1
  std::vector<double>& starchLeaf = x.internalCarbon.starchLeaf;
  std::vector<double>& sugarSapwood = x.internalCarbon.sugarSapwood;
  std::vector<double>& starchSapwood = x.internalCarbon.starchSapwood;
  // 
  // DataFrame internalMortality = Rcpp::as<Rcpp::DataFrame>(x["internalMortality"]);
  // NumericVector N_dead = internalMortality["N_dead"];
  // NumericVector N_starvation = internalMortality["N_starvation"];
  // NumericVector N_dessication = internalMortality["N_dessication"];
  // NumericVector N_burnt = internalMortality["N_burnt"];
  // NumericVector N_resprouting_stumps = internalMortality["N_resprouting_stumps"];
  // NumericVector Cover_dead = internalMortality["Cover_dead"];
  // NumericVector Cover_starvation = internalMortality["Cover_starvation"];
  // NumericVector Cover_dessication = internalMortality["Cover_dessication"];
  // NumericVector Cover_burnt = internalMortality["Cover_burnt"];
  // NumericVector Cover_resprouting_stumps = internalMortality["Cover_resprouting_stumps"];
  // NumericVector Snag_smallbranches = internalMortality["Snag_smallbranches"];
  // NumericVector Snag_largewood = internalMortality["Snag_largewood"];
  // 
  std::vector<double>& allocationTarget = x.internalAllocation.allocationTarget;
  std::vector<double>& sapwoodAreaTarget = x.internalAllocation.sapwoodAreaTarget;
  std::vector<double>& leafAreaTarget = x.internalAllocation.leafAreaTarget;
  std::vector<double>& fineRootBiomassTarget = x.internalAllocation.fineRootBiomassTarget;
  std::vector<double>& crownBudPercent = x.internalAllocation.crownBudPercent;
  // 
  // DataFrame internalPhenology = Rcpp::as<Rcpp::DataFrame>(x["internalPhenology"]);
  // LogicalVector leafUnfolding = internalPhenology["leafUnfolding"];
  // LogicalVector budFormation = internalPhenology["budFormation"];
  // LogicalVector leafSenescence = internalPhenology["leafSenescence"];
  // 
  // DataFrame Plants = Rcpp::as<Rcpp::DataFrame>(spwbOutput["Plants"]);
  // List PlantsInst;
  // NumericVector Ag = Plants["GrossPhotosynthesis"];
  // NumericVector LFMC = Plants["LFMC"];
  // NumericVector PARcohort= Plants["FPAR"];
  // NumericMatrix AgStep, AnStep;
  // int numSteps = 1;
  // if(transpirationMode!="Granier") {
  //   PlantsInst = spwbOutput["PlantsInst"];
  //   AgStep  =  Rcpp::as<Rcpp::NumericMatrix>(PlantsInst["Ag"]);
  //   AnStep  =  Rcpp::as<Rcpp::NumericMatrix>(PlantsInst["An"]);
  //   numSteps = AgStep.ncol();
  // }


  // 
  // //Anatomy parameters
  // DataFrame paramsAnatomy = Rcpp::as<Rcpp::DataFrame>(x["paramsAnatomy"]);
  // NumericVector Hmed = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["Hmed"]);
  std::vector<double>& SLA = x.paramsAnatomy.SLA;
  std::vector<double>& r635 = x.paramsAnatomy.r635;
  // NumericVector Al2As = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["Al2As"]);
  // NumericVector Ar2Al;
  // NumericVector RLD;
  // if(transpirationMode=="Granier") Ar2Al = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["Ar2Al"]);
  // else RLD = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["RLD"]);
  std::vector<double>& WoodDensity = x.paramsAnatomy.WoodDensity;
  std::vector<double>& LeafDensity =x.paramsAnatomy.LeafDensity;
  // NumericVector FineRootDensity = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["FineRootDensity"]);
  // NumericVector SRL = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["SRL"]);
  // NumericVector conduit2sapwood = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["conduit2sapwood"]);
  // 
  // //Growth parameters
  // DataFrame paramsGrowth = Rcpp::as<Rcpp::DataFrame>(x["paramsGrowth"]);
  // NumericVector fHDmin= paramsGrowth["fHDmin"];
  // NumericVector fHDmax= paramsGrowth["fHDmax"];
  // NumericVector WoodC = Rcpp::as<Rcpp::NumericVector>(paramsGrowth["WoodC"]);
  // NumericVector RERleaf = Rcpp::as<Rcpp::NumericVector>(paramsGrowth["RERleaf"]);
  // NumericVector RERsapwood = Rcpp::as<Rcpp::NumericVector>(paramsGrowth["RERsapwood"]);
  // NumericVector RERfineroot = Rcpp::as<Rcpp::NumericVector>(paramsGrowth["RERfineroot"]);
  // NumericVector CCleaf = Rcpp::as<Rcpp::NumericVector>(paramsGrowth["CCleaf"]);
  // NumericVector CCsapwood = Rcpp::as<Rcpp::NumericVector>(paramsGrowth["CCsapwood"]);
  // NumericVector CCfineroot = Rcpp::as<Rcpp::NumericVector>(paramsGrowth["CCfineroot"]);
  // NumericVector RGRleafmax = Rcpp::as<Rcpp::NumericVector>(paramsGrowth["RGRleafmax"]);
  // NumericVector RGRcambiummax = Rcpp::as<Rcpp::NumericVector>(paramsGrowth["RGRcambiummax"]);
  // NumericVector RGRsapwoodmax = Rcpp::as<Rcpp::NumericVector>(paramsGrowth["RGRsapwoodmax"]);
  // NumericVector RGRfinerootmax = Rcpp::as<Rcpp::NumericVector>(paramsGrowth["RGRfinerootmax"]);
  // NumericVector SRsapwood = Rcpp::as<Rcpp::NumericVector>(paramsGrowth["SRsapwood"]);
  // NumericVector SRfineroot = Rcpp::as<Rcpp::NumericVector>(paramsGrowth["SRfineroot"]);
  // NumericVector RSSG = Rcpp::as<Rcpp::NumericVector>(paramsGrowth["RSSG"]);
  // 
  // //Mortality/regeneration parameters
  // DataFrame paramsMortalityRegeneration = Rcpp::as<Rcpp::DataFrame>(x["paramsMortalityRegeneration"]);
  // NumericVector MortalityBaselineRate = Rcpp::as<Rcpp::NumericVector>(paramsMortalityRegeneration["MortalityBaselineRate"]);
  // NumericVector SurvivalModelStep = Rcpp::as<Rcpp::NumericVector>(paramsMortalityRegeneration["SurvivalModelStep"]);
  // NumericVector SurvivalB0 = Rcpp::as<Rcpp::NumericVector>(paramsMortalityRegeneration["SurvivalB0"]);
  // NumericVector SurvivalB1 = Rcpp::as<Rcpp::NumericVector>(paramsMortalityRegeneration["SurvivalB1"]);
  // NumericVector RecrTreeDensity = Rcpp::as<Rcpp::NumericVector>(paramsMortalityRegeneration["RecrTreeDensity"]);
  // NumericVector RecrTreeDBH = Rcpp::as<Rcpp::NumericVector>(paramsMortalityRegeneration["RecrTreeDBH"]);
  // NumericVector IngrowthTreeDensity = Rcpp::as<Rcpp::NumericVector>(paramsMortalityRegeneration["IngrowthTreeDensity"]);
  // NumericVector IngrowthTreeDBH = Rcpp::as<Rcpp::NumericVector>(paramsMortalityRegeneration["IngrowthTreeDBH"]);
  // NumericVector RespFire = Rcpp::as<Rcpp::NumericVector>(paramsMortalityRegeneration["RespFire"]);
  // NumericVector RespDist = Rcpp::as<Rcpp::NumericVector>(paramsMortalityRegeneration["RespDist"]);
  // 
  // //Phenology parameters
  // DataFrame paramsPhenology = Rcpp::as<Rcpp::DataFrame>(x["paramsPhenology"]);
  // CharacterVector phenoType = Rcpp::as<Rcpp::CharacterVector>(paramsPhenology["PhenologyType"]);
  // NumericVector leafDuration = Rcpp::as<Rcpp::NumericVector>(paramsPhenology["LeafDuration"]);
  // 
  // //Transpiration parameters
  // DataFrame paramsTransp = Rcpp::as<Rcpp::DataFrame>(x["paramsTranspiration"]);
  // NumericVector Psi_Extract, Kmax_stemxylem, Plant_kmax, VCleaf_kmax, VCleaf_c, VCleaf_d;
  // NumericVector VCstem_kmax, VCstem_c, VCstem_d, VCroot_kmaxVEC, VCroot_c, VCroot_d, VGrhizo_kmaxVEC;
  // NumericVector WUE_par(numCohorts, 0.3643);
  // VCleaf_c = Rcpp::as<Rcpp::NumericVector>(paramsTransp["VCleaf_c"]);
  // VCleaf_d = Rcpp::as<Rcpp::NumericVector>(paramsTransp["VCleaf_d"]);
  // VCstem_c = Rcpp::as<Rcpp::NumericVector>(paramsTransp["VCstem_c"]);
  // VCstem_d = Rcpp::as<Rcpp::NumericVector>(paramsTransp["VCstem_d"]);
  // if(transpirationMode=="Granier") {
  //   Psi_Extract = Rcpp::as<Rcpp::NumericVector>(paramsTransp["Psi_Extract"]);
  //   if(paramsTransp.containsElementNamed("WUE_par")) {
  //     WUE_par = Rcpp::as<Rcpp::NumericVector>(paramsTransp["WUE_par"]);
  //   }
  //   if(paramsTransp.containsElementNamed("WUE_decay")) { //For compatibility with previous versions (2.7.5)
  //     WUE_par = Rcpp::as<Rcpp::NumericVector>(paramsTransp["WUE_decay"]);
  //   }
  // } else {
  //   Kmax_stemxylem = paramsTransp["Kmax_stemxylem"];
  //   Plant_kmax= paramsTransp["Plant_kmax"];
  //   VCleaf_kmax = Rcpp::as<Rcpp::NumericVector>(paramsTransp["VCleaf_kmax"]);
  //   VCstem_kmax = paramsTransp["VCstem_kmax"];
  //   VCroot_kmaxVEC= paramsTransp["VCroot_kmax"];
  //   VCroot_c = paramsTransp["VCroot_c"];
  //   VCroot_d = paramsTransp["VCroot_d"];
  //   VGrhizo_kmaxVEC= paramsTransp["VGrhizo_kmax"];
  // }
  // 
  // //Water storage parameters
  // DataFrame paramsWaterStorage = Rcpp::as<Rcpp::DataFrame>(x["paramsWaterStorage"]);
  // NumericVector StemPI0 = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["StemPI0"]);
  // NumericVector StemEPS = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["StemEPS"]);
  // NumericVector StemAF = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["StemAF"]);
  // NumericVector Vsapwood = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["Vsapwood"]); //l·m-2 = mm
  // NumericVector LeafAF = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["LeafAF"]);
  // NumericVector Vleaf = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["Vleaf"]); //l·m-2 = mm
  // 
  // //Allometry parameters
  // DataFrame paramsAllometries = Rcpp::as<Rcpp::DataFrame>(x["paramsAllometries"]);
  // NumericVector Hmax  = paramsAnatomy["Hmax"];
  // NumericVector Abt  = paramsAllometries["Abt"];
  // NumericVector Bbt  = paramsAllometries["Bbt"];
  // NumericVector BTsh  = paramsAllometries["BTsh"];
  // 
  // //Inteception parameters
  // DataFrame paramsInterception = Rcpp::as<Rcpp::DataFrame>(x["paramsInterception"]);
  // NumericVector kPAR = paramsInterception["kPAR"];
  // 
  // //Litter
  // DataFrame internalLitter, internalSnags;
  // DataFrame paramsLitterDecomposition;
  // NumericVector internalSOC(0);
  // bool exportLitter = false;
  // if(x.containsElementNamed("internalLitter") && x.containsElementNamed("internalSnags") && x.containsElementNamed("internalSOC") && x.containsElementNamed("paramsLitterDecomposition")) {
  //   paramsLitterDecomposition = Rcpp::as<Rcpp::DataFrame>(x["paramsLitterDecomposition"]);
  //   internalSnags = Rcpp::as<Rcpp::DataFrame>(x["internalSnags"]);
  //   internalLitter = Rcpp::as<Rcpp::DataFrame>(x["internalLitter"]);
  //   internalSOC = Rcpp::as<Rcpp::NumericVector>(x["internalSOC"]);
  //   exportLitter = true;
  // }
  // 
  // //Ring of forming vessels
  // // List ringList = as<Rcpp::List>(x["internalRings"]);
  // 
  // //Daily output vectors
  // DataFrame labileCarbonBalance = as<DataFrame>(modelOutput["LabileCarbonBalance"]);
  // NumericVector GrossPhotosynthesis = labileCarbonBalance["GrossPhotosynthesis"];
  // NumericVector MaintenanceRespiration = labileCarbonBalance["MaintenanceRespiration"];
  // NumericVector GrowthCosts = labileCarbonBalance["GrowthCosts"];
  // NumericVector RootExudation = labileCarbonBalance["RootExudation"];
  // NumericVector LabileCarbonBalance = labileCarbonBalance["LabileCarbonBalance"];
  // NumericVector PlantSugarLeaf = labileCarbonBalance["SugarLeaf"];
  // NumericVector PlantStarchLeaf = labileCarbonBalance["StarchLeaf"];
  // NumericVector PlantSugarSapwood = labileCarbonBalance["SugarSapwood"];
  // NumericVector PlantStarchSapwood = labileCarbonBalance["StarchSapwood"];
  // NumericVector PlantSugarTransport = labileCarbonBalance["SugarTransport"];
  for(int c=0;c<numCohorts;c++) {
    GROWTHres.LCBres.GrossPhotosynthesis[c] = 0.0;
    GROWTHres.LCBres.MaintenanceRespiration[c] = 0.0;
    GROWTHres.LCBres.GrowthCosts[c] = 0.0;
    GROWTHres.LCBres.RootExudation[c] = 0.0;
    GROWTHres.LCBres.LabileCarbonBalance[c] = 0.0;
    GROWTHres.LCBres.SugarLeaf[c] = 0.0;
    GROWTHres.LCBres.StarchLeaf[c] = 0.0;
    GROWTHres.LCBres.SugarSapwood[c] = 0.0;
    GROWTHres.LCBres.StarchSapwood[c] = 0.0;
    GROWTHres.LCBres.SugarTransport[c] = 0.0;
    GROWTHres.GMres.SAgrowth[c] = 0.0;
    GROWTHres.GMres.LAgrowth[c] = 0.0;
    GROWTHres.GMres.FRAgrowth[c] = 0.0;
    GROWTHres.GMres.StarvationRate[c] = 0.0;
    GROWTHres.GMres.DessicationRate[c] = 0.0;
    GROWTHres.GMres.MortalityRate[c] = 0.0;
    GROWTHres.PSres.LeafBiomass[c] = 0.0;
    GROWTHres.PSres.SapwoodBiomass[c] = 0.0;
    GROWTHres.PSres.FineRootBiomass[c] = 0.0;
    GROWTHres.PSres.LeafArea[c] =  0.0;
    GROWTHres.PSres.SapwoodArea[c] = 0.0;
    GROWTHres.PSres.FineRootArea[c] = 0.0;
    GROWTHres.PSres.HuberValue[c] = 0.0;
    GROWTHres.PSres.RootAreaLeafArea[c] = 0.0;
  }
  // 
  // DataFrame plantStructure = as<DataFrame>(modelOutput["PlantStructure"]);
  // NumericVector LeafBiomass = plantStructure["LeafBiomass"];
  // NumericVector SapwoodBiomass = plantStructure["SapwoodBiomass"];
  // NumericVector FineRootBiomass = plantStructure["FineRootBiomass"];
  // NumericVector LeafArea = plantStructure["LeafArea"];
  // NumericVector SapwoodArea = plantStructure["SapwoodArea"];
  // NumericVector FineRootArea = plantStructure["FineRootArea"];
  // NumericVector HuberValue = plantStructure["HuberValue"];
  // NumericVector RootAreaLeafArea = plantStructure["RootAreaLeafArea"];
  // NumericVector OutputDBH = plantStructure["DBH"];
  // NumericVector OutputHeight = plantStructure["Height"];
  // 
  // DataFrame growthMortality = as<DataFrame>(modelOutput["GrowthMortality"]);
  // NumericVector SAgrowth = growthMortality["SAgrowth"];
  // NumericVector LAgrowth = growthMortality["LAgrowth"];
  // NumericVector FRAgrowth = growthMortality["FRAgrowth"];
  // NumericVector StarvationRate = growthMortality["StarvationRate"];
  // NumericVector DessicationRate = growthMortality["DessicationRate"];
  // NumericVector MortalityRate = growthMortality["MortalityRate"];
  // 
  // 
  // 
  std::vector<double> deltaLAgrowth(numCohorts,0.0);
  std::vector<double> deltaSAgrowth(numCohorts,0.0);

  //Stand carbon balance
  double standGrossPrimaryProduction = 0.0;
  double standMaintenanceRespiration = 0.0;
  double standSynthesisRespiration = 0.0;

  double equilibriumLeafSugarConc = equilibriumLeafTotalConc - nonSugarConcentration;
  double equilibriumSapwoodSugarConc = equilibriumSapwoodTotalConc - nonSugarConcentration;

  double rcellmax = relative_expansion_rate(0.0 ,30.0, -1.0, 0.5, 0.05, 5.0);



  // Initial Biomass balance
  std::string units_g_ind = "g_ind";
  std::string units_gC_m2 = "gC_m2";
  fillCarbonCompartments_c(GROWTHcomm.initialFinalCC.ccIni_g_ind, x, units_g_ind);
  fillCarbonCompartments_c(GROWTHcomm.initialFinalCC.ccIni_gC_m2, x, units_gC_m2);
  
  // NumericVector initialTotalCohortBiomass_gC_m2 = ccIni_gC_m2["TotalBiomass"];
  // double initialStand_gC_m2 = sum(initialTotalCohortBiomass_gC_m2);
  // // Rcout<< "initial cohort gC_m2 = "<<sum(initialTotalCohortBiomass_gC_m2)<<"\n";
  // if(exportLitter) {
  //   NumericVector sn_smallBranches = internalSnags["SmallBranches"];
  //   NumericVector sn_largeWood = internalSnags["LargeWood"];
  //   initialStand_gC_m2 += (sum(sn_smallBranches) + sum(sn_largeWood));
  //   NumericVector li_leaves = internalLitter["Leaves"];
  //   NumericVector li_twigs = internalLitter["Twigs"];
  //   NumericVector li_smallBranches = internalLitter["SmallBranches"];
  //   NumericVector li_largeWood = internalLitter["LargeWood"];
  //   NumericVector li_coarseRoots = internalLitter["CoarseRoots"];
  //   NumericVector li_fineRoots = internalLitter["FineRoots"];
  //   initialStand_gC_m2 += (sum(li_leaves) + sum(li_twigs) + sum(li_smallBranches) + sum(li_largeWood) + sum(li_coarseRoots) + sum(li_fineRoots));
  //   initialStand_gC_m2 += sum(internalSOC);
  // }
  // initialStand_gC_m2 += (sum(Snag_smallbranches) + sum(Snag_largewood));
  // 
  std::vector<double> LeafBiomassBalance(numCohorts,0.0), FineRootBiomassBalance(numCohorts,0.0);
  // 
  // DataFrame plantBiomassBalance = as<DataFrame>(modelOutput["PlantBiomassBalance"]);
  // fillInitialPlantBiomassBalance(plantBiomassBalance, ccIni_g_ind, above);
  // 
  std::vector<double>& Volume_leaves = GROWTHcomm.initialFinalCC.ccIni_g_ind.LeafStorageVolume;
  std::vector<double>& Volume_sapwood = GROWTHcomm.initialFinalCC.ccIni_g_ind.SapwoodStorageVolume;
  std::vector<double>& Starch_max_leaves = GROWTHcomm.initialFinalCC.ccIni_g_ind.LeafStarchMaximumConcentration;
  std::vector<double>& Starch_max_sapwood = GROWTHcomm.initialFinalCC.ccIni_g_ind.SapwoodStarchMaximumConcentration;
  std::vector<double>& LeafStructBiomass = GROWTHcomm.initialFinalCC.ccIni_g_ind.LeafStructuralBiomass;
  std::vector<double>& AbovegroundWoodBiomass = GROWTHcomm.initialFinalCC.ccIni_g_ind.AbovegroundWoodBiomass;
  std::vector<double>& BelowgroundWoodBiomass = GROWTHcomm.initialFinalCC.ccIni_g_ind.BelowgroundWoodBiomass;
  std::vector<double>& SapwoodLivingStructBiomass = GROWTHcomm.initialFinalCC.ccIni_g_ind.SapwoodLivingStructuralBiomass;
  std::vector<double>& TwigLivingStructuralBiomass = GROWTHcomm.initialFinalCC.ccIni_g_ind.TwigLivingStructuralBiomass;
  std::vector<double>& TotalLivingBiomass = GROWTHcomm.initialFinalCC.ccIni_g_ind.TotalLivingBiomass;


  //For survival model based on basal area
  double treeBasalArea = 0.0;
  for(int j=0;j<numCohorts;j++){
    if(!std::isnan(DBH[j])) treeBasalArea += N[j]*3.141593*pow(DBH[j]/200,2.0);
  }

  /////////////////////////////////////////////////////////////////////
  ///// B. LABILE CARBON BALANCE, GROWTH AND SENESCENCE BY COHORT /////
  /////////////////////////////////////////////////////////////////////
  for(int j=0;j<numCohorts;j++){
    if(N[j] > 0.0) {
      
      std::vector<double> Vj(nlayers);
      std::vector<double> Lj(nlayers);
      for(int l=0;l<nlayers;l++) {
        Vj[l] = x.belowLayers.V(j,l);
        Lj[l] = x.belowLayers.L(j,l);
      }
      
      ///// B1. COSTS AND DERIVED VARIABLES /////
      double costPerLA = 1000.0*x.paramsGrowth.CCleaf[j]/SLA[j]; // Construction cost in g gluc · m-2 of leaf area
      double twigCostPerLA = 1000.0*x.paramsGrowth.CCsapwood[j]*(r635[j] - 1.0)/SLA[j];// Include cost of twigs i<n g gluc · m-2 of leaf area
      double costPerSA = x.paramsGrowth.CCsapwood[j]*sapwoodStructuralBiomass_c(1.0, H[j], Lj, Vj,WoodDensity[j]); // Construction cost in g gluc · cm-2 of sapwood area
      if(ctype[j]=="tree") {
        double section_cm2 = pow(DBH[j]/2.0, 2.0)*M_PI;
        double deltaDBH = 2.0*sqrt(pow(DBH[j]/2.0,2.0)+(1.0/M_PI)) - DBH[j];
        double fHmod = std::max(0.0,std::min(1.0,(1.0-((H[j]-137.0)/(x.paramsAnatomy.Hmax[j]-137.0)))));
        double fHD = (x.paramsGrowth.fHDmin[j]*(PARcohort[j]/100.0) + x.paramsGrowth.fHDmax[j]*(1.0-(PARcohort[j]/100.0)))*fHmod;
        double deltaH = fHD*deltaDBH; //Increase in height
        costPerSA *= (1.0 + deltaH); //Cost of increase in sapwood due to height increase
        costPerSA += x.paramsGrowth.CCsapwood[j]*(section_cm2 - SA[j])*deltaH*WoodDensity[j]; //Cost of increase in heartwood due to height increase
      }

      std::vector<double> deltaFRBgrowth(nlayers, 0.0);

      double LAexpanded = leafArea_c(LAI_expanded[j], N[j]);
      double LAlive = leafArea_c(LAI_live[j], N[j]);
      double LAdead = leafArea_c(LAI_dead[j], N[j]);

      double minimumStarchForSecondaryGrowth = Starch_max_sapwood[j]*x.paramsGrowth.RSSG[j];
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
      // double k_phloem = NA_REAL;
      // if(subdailyCarbonBalance) k_phloem = VCstem_kmax[j]*phloemConductanceFactor*(0.018/1000.0);

      //Xylogenesis
      // List ring = ringList[j];
      double rleafcell = medfate::NA_DOUBLE, rcambiumcell = medfate::NA_DOUBLE;
      std::vector<double> rfineroot(nlayers);
      double relative_hormone_factor = std::max(0.0, std::min(1.0, LAI_expanded[j]/LAI_nocomp[j]));
      if(x.control.transpirationMode=="Granier") {
        // grow_ring(ring, PlantPsi[j] ,tday, 10.0);
        rleafcell = std::min(rcellmax, relative_expansion_rate(x.internalWater.PlantPsi[j] ,tday, -1.0, 0.5,0.05,5.0));
        rcambiumcell = std::min(rcellmax, relative_hormone_factor*relative_expansion_rate(x.internalWater.PlantPsi[j] ,tday, -1.0, 0.5,0.05,5.0));
        for(int l=0;l<nlayers;l++) rfineroot[l] = std::min(rcellmax, relative_expansion_rate(x.soil.getPsi(l) ,tday, -1.0 ,0.5,0.05,5.0));
        // if(j==0) Rcout<<j<< " Psi:"<< PlantPsi[j]<< " r:"<< rcambiumcell<<"\n";
      } else {
        rleafcell = std::min(rcellmax, relative_expansion_rate(psiRootCrown[j] ,tcan_day, -1.0, 0.5,0.05,5.0));
        rcambiumcell = std::min(rcellmax, relative_hormone_factor*relative_expansion_rate(psiRootCrown[j] ,tcan_day, -1.0, 0.5,0.05,5.0));
        for(int l=0;l<nlayers;l++) rfineroot[l] = std::min(rcellmax, relative_expansion_rate(x.belowLayers.RhizoPsi(j,l), x.soil.getTemp(l), -1.0, 0.5,0.05,5.0));
        // Rcout<<j<< " rcellmax "<< rcellmax<<" tcan_day "<< tcan_day<< " Psi:"<< psiRootCrown[j]<< " r:"<< rcambiumcell<<"\n";
      }

      //Respiratory biomass (g dw · ind-1)
      double B_resp_leaves = LeafStructBiomass[j] + sugarLeaf[j]*(Volume_leaves[j]*glucoseMolarMass);
      double B_resp_sapwood = SapwoodLivingStructBiomass[j]  + sugarSapwood[j]*(Volume_sapwood[j]*glucoseMolarMass);
      double B_resp_twig = TwigLivingStructuralBiomass[j];
      double B_resp_fineroots = fineRootBiomass[j];

      if(!subdailyCarbonBalance) {

        ////// B2. PHOTOSYNTHESIS //////
        if(LAexpanded>0.0) {
          standGrossPrimaryProduction += Ag[j]; //Add GPP in gC · m-2
          //gross fotosynthesis
          double leafAgC = Ag[j]/(N[j]/10000.0); //Translate g C · m-2 to g C ·ind-1
          leafAgG = leafAgC*(glucoseMolarMass/(carbonMolarMass*6.0)); // from g C·ind-1 to g gluc · ind-1
        }

        ///// B3. MAINTENANCE RESPIRATION ///////
        double QR = qResp_c(tday);
        if(LAexpanded>0.0) {
          leafRespDay = B_resp_leaves*x.paramsGrowth.RERleaf[j]*QR*std::min(1.0, pow(PARcohort[j]/100.0, 0.5)); //Reduction under shade
        }
        twigResp = B_resp_twig*x.paramsGrowth.RERsapwood[j]*QR;
        sapwoodResp = B_resp_sapwood*x.paramsGrowth.RERsapwood[j]*QR;
        finerootResp = B_resp_fineroots*x.paramsGrowth.RERfineroot[j]*QR*(LAexpanded/LAlive);

        ///// B4. LEAF/TWIG GROWTH /////
        if(x.internalPhenology.leafUnfolding[j]) {
          double deltaLApheno = std::max(leafAreaTarget[j]*(1.0 - StemPLC[j]) - LAexpanded, 0.0);
          double deltaLAsink = std::min(deltaLApheno, (crownBudPercent[j]/100.0)*SA[j]*x.paramsGrowth.RGRleafmax[j]*(rleafcell/rcellmax));
          if(!sinkLimitation) deltaLAsink = std::min(deltaLApheno, (crownBudPercent[j]/100.0)*SA[j]*x.paramsGrowth.RGRleafmax[j]); //Deactivates temperature and turgor limitation
          double deltaLAavailable = 0.0;
          deltaLAavailable = std::max(0.0, (starchSapwood[j]-minimumStarchForPrimaryGrowth)*(glucoseMolarMass*Volume_sapwood[j])/costPerLA);
          deltaLAgrowth[j] = std::min(deltaLAsink, deltaLAavailable);
          growthCostLA = deltaLAgrowth[j]*costPerLA;
          synthesisRespLA = growthCostLA*(x.paramsGrowth.CCleaf[j] - 1.0)/x.paramsGrowth.CCleaf[j];
          twigGrowthCostLA = deltaLAgrowth[j]*twigCostPerLA;
          twigSynthesisRespLA = twigGrowthCostLA*(x.paramsGrowth.CCsapwood[j] - 1.0)/x.paramsGrowth.CCsapwood[j];
        }

        ///// B5. SAPWOOD GROWTH /////
        if(LAexpanded>0.0) {
          double deltaSAsink = NA_REAL;
          if(!NumericVector::is_na(DBH[j])) { //Trees
            deltaSAsink = (3.141592*DBH[j]*x.paramsGrowth.RGRcambiummax[j]*(rcambiumcell/rcellmax));
            if(!sinkLimitation) deltaSAsink = 3.141592*DBH[j]*x.paramsGrowth.RGRcambiummax[j]; //Deactivates temperature and turgor limitation
          } else { // Shrubs
            deltaSAsink = (SA[j]*x.paramsGrowth.RGRsapwoodmax[j]*(rcambiumcell/rcellmax));
            if(!sinkLimitation) deltaSAsink = SA[j]*x.paramsGrowth.RGRsapwoodmax[j]; //Deactivates temperature and turgor limitation
          }
          double deltaSAavailable = std::max(0.0, (starchSapwood[j]-minimumStarchForSecondaryGrowth)*(glucoseMolarMass*Volume_sapwood[j])/costPerSA);
          deltaSAgrowth[j] = std::min(deltaSAsink, deltaSAavailable);
          // Rcout<< rcambiumcell<<" deltaSAsink "<< deltaSAsink << " dSAgrowth "<< deltaSAgrowth<<"\n";
          growthCostSA = deltaSAgrowth[j]*costPerSA; //increase cost (may be non zero if leaf growth was charged onto sapwood)
          synthesisRespSA = growthCostSA*(x.paramsGrowth.CCsapwood[j] - 1.0)/x.paramsGrowth.CCsapwood[j];
        }


        ///// B6. FINE ROOT GROWTH /////
        if(fineRootBiomass[j] < fineRootBiomassTarget[j]) {
          for(int s = 0;s<nlayers;s++) {
            double deltaFRBpheno = std::max(fineRootBiomassTarget[j] - fineRootBiomass[j], 0.0);
            double deltaFRBsink = (Vj[s]*fineRootBiomass[j])*x.paramsGrowth.RGRfinerootmax[j]*(rfineroot[s]/rcellmax);
            if(!sinkLimitation) deltaFRBsink = (Vj[s]*fineRootBiomass[j])*x.paramsGrowth.RGRfinerootmax[j]; //Deactivates temperature and turgor limitation
            double deltaFRBavailable = std::max(0.0,(starchSapwood[j]-minimumStarchForPrimaryGrowth)*(glucoseMolarMass*Volume_sapwood[j])/x.paramsGrowth.CCfineroot[j]);
            deltaFRBgrowth[s] = std::min(deltaFRBpheno, std::min(deltaFRBsink, deltaFRBavailable));
            growthCostFRB += deltaFRBgrowth[s]*x.paramsGrowth.CCfineroot[j];
            synthesisRespFRB += deltaFRBgrowth[s]*(x.paramsGrowth.CCfineroot[j] - 1.0);
          }
        }


        ///// B7a PARTIAL CARBON BALANCE: photosynthesis, maintenance respiration and growth /////
        double leafSugarMassDelta = leafAgG - leafRespDay;
        double sapwoodSugarMassDelta =  - twigResp - sapwoodResp - finerootResp;
        double sapwoodStarchMassDelta =  - growthCostFRB - growthCostLA - twigGrowthCostLA - growthCostSA;
        sugarSapwood[j] += sapwoodSugarMassDelta/(Volume_sapwood[j]*glucoseMolarMass);
        starchSapwood[j] += sapwoodStarchMassDelta/(Volume_sapwood[j]*glucoseMolarMass);
        if(LAexpanded>0.0) sugarLeaf[j] += leafSugarMassDelta/(Volume_leaves[j]*glucoseMolarMass);

        ///// B7b PHLOEM TRANSPORT AND SUGAR-STARCH DYNAMICS /////
        if(LAexpanded>0.0) {
          double ff = (sugarLeaf[j]-sugarSapwood[j])/2.0;
          sugarLeaf[j] -=ff;
          GROWTHres.LCBres.SugarTransport[j] = 1000.0*(ff*Volume_leaves[j])/(3600.0*24.0); //mmol · s-1
          sugarSapwood[j] +=(Volume_leaves[j]/Volume_sapwood[j])*ff;
          double conversionLeaf = std::max(-starchLeaf[j], sugarLeaf[j] - equilibriumLeafSugarConc);
          starchLeaf[j] +=conversionLeaf;
          sugarLeaf[j] -=conversionLeaf;
        }
        double conversionSapwood = std::max(-starchSapwood[j], sugarSapwood[j] - equilibriumSapwoodSugarConc);
        starchSapwood[j] +=conversionSapwood;
        sugarSapwood[j] -=conversionSapwood;

        // LABILE BALANCE OUTPUT VARIABLES
        GROWTHres.LCBres.GrossPhotosynthesis[j] = leafAgG/TotalLivingBiomass[j]; //Ag in g gluc · gdry-1
        GROWTHres.LCBres.GrowthCosts[j] = (growthCostLA + twigGrowthCostLA + growthCostSA + growthCostFRB)/TotalLivingBiomass[j]; //growth cost in g gluc · gdry-1
        GROWTHres.LCBres.MaintenanceRespiration[j] = (leafRespDay+twigResp + sapwoodResp+finerootResp)/TotalLivingBiomass[j];

      } else {
  //       //Get output matrices
  //       List labileCBInst = modelOutput["LabileCarbonBalanceInst"];
  //       NumericMatrix LabileCarbonBalanceInst = labileCBInst["LabileCarbonBalance"];
  //       NumericMatrix GrossPhotosynthesisInst = labileCBInst["GrossPhotosynthesis"];
  //       NumericMatrix MaintenanceRespirationInst = labileCBInst["MaintenanceRespiration"];
  //       NumericMatrix GrowthCostsInst = labileCBInst["GrowthCosts"];
  //       NumericMatrix RootExudationInst = labileCBInst["RootExudation"];
  //       NumericMatrix PlantSugarTransportInst = labileCBInst["SugarTransport"];
  //       NumericMatrix PlantSugarLeafInst = labileCBInst["SugarLeaf"];
  //       NumericMatrix PlantStarchLeafInst = labileCBInst["StarchLeaf"];
  //       NumericMatrix PlantSugarSapwoodInst = labileCBInst["SugarSapwood"]; 
  //       NumericMatrix PlantStarchSapwoodInst = labileCBInst["StarchSapwood"];
  //       
  //       //Carbon balance and growth by steps
  //       for(int s=0;s<numSteps;s++) {
  //         
  //         double leafRespStep = 0.0;
  //         double leafAgStepG = 0.0, leafAnStepG=0.0;
  //         double growthCostLAStep = 0.0;
  //         double twigGrowthCostLAStep = 0.0;
  //         double growthCostSAStep = 0.0;
  //         double growthCostFRBStep = 0.0;   
  //         
  //         //B2. LEAF PHOTOSYNTHESIS and MAINTENANCE RESPIRATION
  //         if(LAexpanded>0.0){
  //           standGrossPrimaryProduction += AgStep(j,s); //Add GPP in gC · m-2
  //           
  //           double leafAgStepC = AgStep(j,s)/(N[j]/10000.0); //Translate g C · m-2 · h-1 to g C · h-1
  //           double leafAnStepC = AnStep(j,s)/(N[j]/10000.0); //Translate g C · m-2 · h-1 to g C · h-1
  //           leafAgStepG = leafAgStepC*(glucoseMolarMass/(carbonMolarMass*6.0)); // from g C· h-1 to g gluc · h-1
  //           leafAnStepG = leafAnStepC*(glucoseMolarMass/(carbonMolarMass*6.0)); // from g C· h-1 to g gluc · h-1
  //           leafAgG += leafAgStepG;
  //           leafRespStep = leafAgStepG - leafAnStepG; //Respiration as Ag - An
  //           GrossPhotosynthesisInst(j,s) = leafAgStepG/TotalLivingBiomass[j]; //Ag in g gluc · gdry-1
  //           GrossPhotosynthesis[j] += GrossPhotosynthesisInst(j,s); 
  //         }
  //         
  //         //B3. MAINTENANCE RESPIRATION
  //         double QR = qResp_c(Tcan[s]);
  //         leafRespDay +=leafRespStep;
  //         double sapwoodRespStep = B_resp_sapwood*RERsapwood[j]*QR/((double) numSteps);
  //         sapwoodResp += sapwoodRespStep;
  //         double twigRespStep = B_resp_twig*RERsapwood[j]*QR/((double) numSteps);
  //         twigResp += twigRespStep;
  //         double finerootRespStep = B_resp_fineroots*RERfineroot[j]*QR*(LAexpanded/LAlive)/((double) numSteps);
  //         finerootResp += finerootRespStep;
  //         MaintenanceRespirationInst(j,s) = (leafRespStep+twigRespStep+sapwoodRespStep+finerootRespStep)/TotalLivingBiomass[j];//Rm in g gluc· gdry-1
  //         MaintenanceRespiration[j] += MaintenanceRespirationInst(j,s); 
  //         
  //         //B.4 Leaf growth
  //         if(leafUnfolding[j]) {
  //           double deltaLApheno = std::max(leafAreaTarget[j]*(1.0 - StemPLC[j]) - LAexpanded, 0.0);
  //           double deltaLAsink = std::min(deltaLApheno, (crownBudPercent[j]/100.0)*(1.0/((double) numSteps))*SA[j]*RGRleafmax[j]*(rleafcell/rcellmax));
  //           if(!sinkLimitation) deltaLAsink = std::min(deltaLApheno, (crownBudPercent[j]/100.0)*(1.0/((double) numSteps))*SA[j]*RGRleafmax[j]); //Deactivates temperature and turgor limitation
  //           //Grow at expense of stem sugar
  //           double deltaLAavailable = std::max(0.0, (starchSapwood[j]-minimumStarchForPrimaryGrowth)*(glucoseMolarMass*Volume_sapwood[j])/costPerLA);
  //           double deltaLAgrowthStep = std::min(deltaLAsink, deltaLAavailable);
  //           growthCostLAStep = deltaLAgrowthStep*costPerLA;
  //           twigGrowthCostLAStep = deltaLAgrowthStep*twigCostPerLA;
  //           deltaLAgrowth[j] += deltaLAgrowthStep;
  //           synthesisRespLA += growthCostLAStep*(CCleaf[j] - 1.0)/CCleaf[j];
  //           twigSynthesisRespLA += twigGrowthCostLAStep*(CCsapwood[j] - 1.0)/CCsapwood[j];
  //         }
  //         
  //         //B5. sapwood area growth
  //         if(LAexpanded>0.0) {
  //           double deltaSAsink = NA_REAL;
  //           if(!NumericVector::is_na(DBH[j])) { //Trees
  //             deltaSAsink = (3.141592*DBH[j]*RGRcambiummax[j]*(rcambiumcell/rcellmax))/((double) numSteps); 
  //             if(!sinkLimitation) deltaSAsink = 3.141592*DBH[j]*RGRcambiummax[j]/((double) numSteps); //Deactivates temperature and turgor limitation
  //           } else { // Shrubs
  //             deltaSAsink = (SA[j]*RGRsapwoodmax[j]*(rcambiumcell/rcellmax))/((double) numSteps); 
  //             if(!sinkLimitation) deltaSAsink = SA[j]*RGRsapwoodmax[j]/((double) numSteps); //Deactivates temperature and turgor limitation
  //           }
  //           double deltaSAavailable = std::max(0.0, (starchSapwood[j] - minimumStarchForSecondaryGrowth)*(glucoseMolarMass*Volume_sapwood[j])/costPerSA);
  //           double deltaSAgrowthStep = std::min(deltaSAsink, deltaSAavailable);
  //           if(deltaSAgrowthStep<0.0) {
  //             Rcout<<deltaSAsink<<" "<< deltaSAavailable<< " "<< starchSapwood[j]<<"\n";
  //             stop("negative growth!"); 
  //           }
  //           growthCostSAStep = deltaSAgrowthStep*costPerSA; //increase cost (may be non zero if leaf growth was charged onto sapwood)
  //           growthCostSA += growthCostSAStep;
  //           deltaSAgrowth[j]  +=deltaSAgrowthStep;
  //           synthesisRespSA += growthCostSAStep*(CCsapwood[j] - 1.0)/CCsapwood[j];
  //         }
  //         
  //         //B6. fine root growth
  //         if(fineRootBiomass[j] < fineRootBiomassTarget[j]) {
  //           for(int l = 0;l<nlayers;l++) {
  //             double deltaFRBpheno = std::max(fineRootBiomassTarget[j] - fineRootBiomass[j], 0.0);
  //             double deltaFRBsink = (1.0/((double) numSteps))*(V(j,l)*fineRootBiomass[j])*RGRfinerootmax[j]*(rfineroot[l]/rcellmax);
  //             if(!sinkLimitation) deltaFRBsink = (1.0/((double) numSteps))*(V(j,l)*fineRootBiomass[j])*RGRfinerootmax[j]; //Deactivates temperature and turgor limitation
  //             double deltaFRBavailable = std::max(0.0, (starchSapwood[j]-minimumStarchForPrimaryGrowth)*(glucoseMolarMass*Volume_sapwood[j])/CCfineroot[j]);
  //             double deltaFRBgrowthStep = std::min(deltaFRBpheno, std::min(deltaFRBsink, deltaFRBavailable));
  //             growthCostFRBStep += deltaFRBgrowthStep*CCfineroot[j];
  //             deltaFRBgrowth[l] += deltaFRBgrowthStep;
  //             synthesisRespFRB += deltaFRBgrowthStep*(CCfineroot[j] - 1.0);
  //           }
  //           growthCostFRB += growthCostFRBStep;
  //         }
  //         
  //         GrowthCostsInst(j,s) = (growthCostLAStep + twigGrowthCostLAStep + growthCostSAStep + growthCostFRBStep)/TotalLivingBiomass[j];
  //         GrowthCosts[j] +=GrowthCostsInst(j,s); //growth cost in g gluc · gdry-1
  //         
  //         ///// B7a PARTIAL CARBON BALANCE: photosynthesis, maintenance respiration and growth /////
  //         double leafSugarMassDeltaStep = leafAgStepG - leafRespStep;
  //         double sapwoodSugarMassDeltaStep =  - twigRespStep - sapwoodRespStep - finerootRespStep; 
  //         double sapwoodStarchMassDeltaStep =  - growthCostFRBStep - growthCostLAStep - twigGrowthCostLAStep - growthCostSAStep;
  //         sugarSapwood[j] += sapwoodSugarMassDeltaStep/(Volume_sapwood[j]*glucoseMolarMass);
  //         starchSapwood[j] += sapwoodStarchMassDeltaStep/(Volume_sapwood[j]*glucoseMolarMass);
  //         if(LAexpanded>0.0) sugarLeaf[j] += leafSugarMassDeltaStep/(Volume_leaves[j]*glucoseMolarMass);
  //         
  //         ///// B7b PHLOEM TRANSPORT AND SUGAR-STARCH DYNAMICS /////  
  //         if(LAexpanded>0.0) {
  //           double ff = (sugarLeaf[j]-sugarSapwood[j])/2.0; 
  //           sugarLeaf[j] -=ff;
  //           PlantSugarTransportInst(j,s) = 1000.0*(ff*Volume_leaves[j])/(3600.0*((double) numSteps)); //mmol · s-1
  //           PlantSugarTransport[j] += PlantSugarTransportInst(j,s)/((double) numSteps); //Average daily rate To calculate daily phloem balance (positive means towards stem)
  //           sugarSapwood[j] +=(Volume_leaves[j]/Volume_sapwood[j])*ff;
  //           double conversionLeaf = std::max(-starchLeaf[j], sugarLeaf[j] - equilibriumLeafSugarConc);
  //           starchLeaf[j] +=conversionLeaf;
  //           sugarLeaf[j] -=conversionLeaf;
  //         }
  //         double conversionSapwood = std::max(-starchSapwood[j], sugarSapwood[j] - equilibriumSapwoodSugarConc);
  //         starchSapwood[j] +=conversionSapwood;
  //         sugarSapwood[j] -=conversionSapwood;
  //         
  //         //Divert to root exudation if sapwood starch is over maximum capacity
  //         if(starchSapwood[j] > Starch_max_sapwood[j]) {
  //           RootExudationInst(j,s) = ((starchSapwood[j] - Starch_max_sapwood[j])*(Volume_sapwood[j]*glucoseMolarMass)/TotalLivingBiomass[j]);
  //           RootExudation[j] += RootExudationInst(j,s);//root exudation in g gluc · gdry-1
  //           starchSapwood[j] = Starch_max_sapwood[j];
  //         }
  //         
  //         //Instantaneous carbon balance
  //         LabileCarbonBalanceInst(j,s) = GrossPhotosynthesisInst(j,s) - MaintenanceRespirationInst(j,s) - GrowthCostsInst(j,s) - RootExudationInst(j,s);
  //         PlantSugarLeafInst(j,s) = sugarLeaf[j];
  //         PlantSugarSapwoodInst(j,s) = sugarSapwood[j];
  //         PlantStarchLeafInst(j,s) = starchLeaf[j];
  //         PlantStarchSapwoodInst(j,s) = starchSapwood[j];
  //         // Rcout<<" coh:"<<j<< " s:"<<s<< " conc leaf: "<< sugarLeaf[j] << " conc sap: "<< sugarSapwood[j]<<" ff: "<<ff<< "\n";
  //       }
      }

      double leafBiomassIncrement = deltaLAgrowth[j]*(1000.0/SLA[j]);
      double finerootBiomassIncrement = std::accumulate(deltaFRBgrowth.begin(), deltaFRBgrowth.end(), 0.0);

      //add MR and SR to stand-level maintenance and synthesis respiration as g C·m-2
      standMaintenanceRespiration += (leafRespDay+twigResp+sapwoodResp+finerootResp)*(N[j]/10000.0)*((carbonMolarMass*6.0)/glucoseMolarMass);
      standSynthesisRespiration += (synthesisRespLA+ twigSynthesisRespLA +synthesisRespSA+synthesisRespFRB)*(N[j]/10000.0)*((carbonMolarMass*6.0)/glucoseMolarMass);

      ///// B8. LEAF SENESCENCE /////
      double propLeafSenescence = 0.0;
      //Leaf senescence due to age (Ca+ accumulation) only in evergreen species
      if(x.paramsPhenology.phenoType[j] == "progressive-evergreen") {
        propLeafSenescence = std::min(1.0,(LAexpanded/(365.25*LAlive*x.paramsPhenology.leafDuration[j])));
      } else if((x.paramsPhenology.phenoType[j] == "oneflush-evergreen") && (x.internalPhenology.leafSenescence[j])) {
        propLeafSenescence = std::min(1.0,(LAexpanded/(LAlive*x.paramsPhenology.leafDuration[j]))); // Fraction of old leaves that die
        x.internalPhenology.leafSenescence[j] = false; //To prevent further loss
      } else if(((x.paramsPhenology.phenoType[j] == "winter-deciduous") || (x.paramsPhenology.phenoType[j] == "winter-semideciduous")) && x.internalPhenology.leafSenescence[j]) {
        propLeafSenescence = 1.0;
        x.internalPhenology.leafSenescence[j] = false; //To prevent further loss
      }
      //Leaf senescence and bud senescence due to drought (only when PLC increases)
      double PLCinc = (StemPLC[j]-StemPLCprev[j]);
      if(PLCinc>0.0) {
        double LeafPDEF = proportionDefoliationWeibull_c(psiApoLeaf[j], x.paramsTranspiration.VCleaf_c[j], x.paramsTranspiration.VCleaf_d[j], 0.88, 10);
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

      ///// B9. SAPWOOD AREA SENESCENCE /////
      //Define sapwood senescence as maximum of turnover and sapwood exceeding the target
      double propSASenescence = x.paramsGrowth.SRsapwood[j]*std::max(0.0,(tday-5.0)/20.0)/(1.0+15.0*exp(-0.01*H[j]));
      double deltaSASenescence = std::max(0.0, SA[j] - sapwoodAreaTarget[j]);
      propSASenescence = std::max(propSASenescence, deltaSASenescence/SA[j]);

      ///// B10. FINE ROOT BIOMASS SENESCENCE /////
      std::vector<double> deltaFRBsenescence(nlayers, 0.0);
      for(int l=0;l<nlayers;l++) {
        double daySenescence = NA_REAL;
        if(x.control.transpirationMode=="Granier") daySenescence = x.paramsGrowth.SRfineroot[j]*std::max(0.0,(tday-5.0)/20.0);
        else daySenescence = x.paramsGrowth.SRfineroot[j]*std::max(0.0,(x.soil.getTemp(l)-5.0)/20.0);
        deltaFRBsenescence[l] = fineRootBiomass[j]*Vj[l]*daySenescence;
      }
      double senescenceFinerootLoss = std::accumulate(deltaFRBsenescence.begin(), deltaFRBsenescence.end(),0.0);
      // if(j==(numCohorts-1)) Rcout<< j << " before translocation "<< sugarLeaf[j]<< " "<< starchLeaf[j]<<"\n";


      ///// 11. TRANSLOCATION (in mol gluc) of labile carbon due to tissue senescence /////
      double translocationSugarLeaf = propLeafSenescence*Volume_leaves[j]*sugarLeaf[j];
      double translocationStarchLeaf = propLeafSenescence*Volume_leaves[j]*starchLeaf[j];
      double translocationSugarSapwood = propSASenescence*Volume_sapwood[j]*sugarSapwood[j];
      if(Volume_leaves[j]>0.0) {
        if(starchLeaf[j] > Starch_max_leaves[j]) { // Add excess leaf starch to translocation
          translocationStarchLeaf += ((starchLeaf[j] - Starch_max_leaves[j])*Volume_leaves[j]);
        }
        sugarLeaf[j] = ((sugarLeaf[j]*Volume_leaves[j]) - translocationSugarLeaf)/Volume_leaves[j];
        starchLeaf[j] = ((starchLeaf[j]*Volume_leaves[j]) - translocationStarchLeaf)/Volume_leaves[j];
      }
      sugarSapwood[j] = ((sugarSapwood[j]*Volume_sapwood[j]) - translocationSugarSapwood)/Volume_sapwood[j];
      starchSapwood[j] = ((starchSapwood[j]*Volume_sapwood[j]) + translocationSugarLeaf + translocationStarchLeaf + translocationSugarSapwood)/Volume_sapwood[j];

      ///// C12. ROOT EXUDATION and close labile carbon balance (non-subdaily carbon balance)
      //Excess sapwood starch carbon is lost as root exudation
      if(starchSapwood[j] > Starch_max_sapwood[j]) {
        GROWTHres.LCBres.RootExudation[j] += ((starchSapwood[j] - Starch_max_sapwood[j])*(Volume_sapwood[j]*glucoseMolarMass)/TotalLivingBiomass[j]);
        starchSapwood[j] = Starch_max_sapwood[j];
      }

      GROWTHres.LCBres.LabileCarbonBalance[j] = GROWTHres.LCBres.GrossPhotosynthesis[j] - GROWTHres.LCBres.MaintenanceRespiration[j] - GROWTHres.LCBres.GrowthCosts[j] - GROWTHres.LCBres.RootExudation[j];

      //CHECK LABILE BIOMASS BALANCE (g gluc)
      //Final labile (g gluc)
      double finalLabile =  glucoseMolarMass*((sugarLeaf[j] + starchLeaf[j])*Volume_leaves[j] + (sugarSapwood[j] + starchSapwood[j])*Volume_sapwood[j]);
      double labileChange = finalLabile - initialLabile;
      double labileBalance_ggluc = GROWTHres.LCBres.LabileCarbonBalance[j]*TotalLivingBiomass[j];
      if(std::abs(labileChange - labileBalance_ggluc) > 0.1) {
        throw medfate::MedfateInternalError(" Labile biomass balance for cohort not closed!");
        // Rcout << "   Initial labile " << initialLabile << " Final Labile " << finalLabile << " Labile change " << (finalLabile-initialLabile);
        // Rcout << " A " << GROWTHres.LCBres.GrossPhotosynthesis[j]*TotalLivingBiomass[j] << " MR " << GROWTHres.LCBres.MaintenanceRespiration[j]*TotalLivingBiomass[j];
        // Rcout << " GC "<< GROWTHres.LCBres.GrowthCosts[j]*TotalLivingBiomass[j] << " RE " << GROWTHres.LCBres.RootExudation[j]*TotalLivingBiomass[j] << " DeltaS " <<  GROWTHres.LCBres.LabileCarbonBalance[j]*TotalLivingBiomass[j]<<"\n";
      }

  //     ///// C13. EXPORT ROOT EXUDATION AND FINEROOT LITTER TO DECOMPOSING POOLS
  //     if(exportLitter) {
  //       //From g/ind to g C/m2
  //       double fineroot_litter = rootCperDry*senescenceFinerootLoss*(N[j]/10000.0); 
  //       addFineRootLitter(speciesNames[j], fineroot_litter,
  //                         internalLitter, paramsLitterDecomposition,
  //                         internalSOC);
  //       //From g gluc/g ind to g C/m2
  //       double metabolicSoil = (6.0*carbonMolarMass/glucoseMolarMass) * RootExudation[j]*TotalLivingBiomass[j]*(N[j]/10000.0); 
  //       internalSOC["SoilMetabolic"] = internalSOC["SoilMetabolic"] + metabolicSoil;
  //     }
  //     //For shrubs, convert SA senescence into mass of standing dead branches (g C/m2)
  //     if(ctype[j]=="shrub") {
  //       double sh_wood_biomass = (deltaSASenescence/SA[j])*(N[j]/10000.0)*AbovegroundWoodBiomass[j]*WoodC[j]; //Aboveground necromass
  //       Snag_smallbranches[j] += sh_wood_biomass;
  //     }
  //     
  //     
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
      SA[j] = SA[j] + deltaSAgrowth[j] - deltaSASenescence;
      // Rcout << j << " deltaSA-G " << deltaSAgrowth[j] << " deltaSA-S "<< deltaSASenescence <<"\n";
      std::vector<double> newFRB(nlayers,0.0);
      for(int s=0;s<nlayers;s++) {
        newFRB[s] = fineRootBiomass[j]*Vj[s] + deltaFRBgrowth[s] - deltaFRBsenescence[s];
      }
      fineRootBiomass[j] = std::accumulate(newFRB.begin(), newFRB.end(), 0.0);
      for(int s=0;s<nlayers;s++) {
        x.belowLayers.V(j,s) = newFRB[s]/fineRootBiomass[j];
        Vj[s] = x.belowLayers.V(j,s);
      }
      //Decrease PLC due to new SA growth
      StemPLC[j] = std::max(0.0, StemPLC[j] - (deltaSAgrowth[j]/SA[j]));
      LeafPLC[j] = std::max(0.0, LeafPLC[j] - (deltaLAgrowth[j]/LAexpanded));
      //Increase crown buds to new SA growth
      crownBudPercent[j] = std::min(100.0, crownBudPercent[j] + 100.0*(deltaSAgrowth[j]/SA[j]));

      ///// B14a. UPDATE DERIVED HYDRAULIC PARAMETERS /////
      if(x.control.transpirationMode=="Granier") {
        Lj = coarseRootLengths_c(Vj, x.soil.getWidths(), 0.5); //Arbitrary ratio (to revise some day)
        for(int l =0;l<nlayers; l++) x.belowLayers.L(j,l) = Lj[l]; 
        x.below.coarseRootSoilVolume[j] = coarseRootSoilVolume_c(Vj, x.soil.getWidths(), 0.5);
      } else { //SPERRY/Sureau
        if(LAlive>0.0) {
  //         //Update Huber value, stem and root hydraulic conductance
  //         double oldstemR = 1.0/VCstem_kmax[j];
  //         double oldrootR = 1.0/VCroot_kmaxVEC[j];
  //         double oldrootprop = oldrootR/(oldrootR+oldstemR);
  //         
  //         // Al2As[j] = (LAlive)/(SA[j]/10000.0);
  //         VCstem_kmax[j]=maximumStemHydraulicConductance(Kmax_stemxylem[j], Hmed[j], Al2As[j] ,H[j], taper); 
  //         
  //         //Update rhizosphere maximum conductance
  //         NumericVector VGrhizo_new = rhizosphereMaximumConductance(Ksat, newFRB, LAI_live[j], N[j],
  //                                                                   SRL[j], FineRootDensity[j], RLD[j]);
  //         for(int s=0;s<nlayers;s++) { 
  //           VGrhizo_kmax(j,s) = VGrhizo_new[s];
  //         }
  //         VGrhizo_kmaxVEC[j] = sum(VGrhizo_kmax(j,_));
  //         
  //         //Update root maximum conductance so that it keeps the same resistance proportion with stem conductance
  //         double newstemR = 1.0/VCstem_kmax[j];
  //         double newrootR = oldrootprop*newstemR/(1.0-oldrootprop);
  //         VCroot_kmaxVEC[j] = 1.0/newrootR;
  //         //Update coarse root soil volume
  //         CRSV[j] = coarseRootSoilVolumeFromConductance(Kmax_stemxylem[j], VCroot_kmaxVEC[j], Al2As[j],
  //                                                       V(j,_), widths, rfc);
  //         //Update coarse root length and root maximum conductance
  //         L(j,_) = coarseRootLengthsFromVolume(CRSV[j], V(j,_), widths, rfc);
  //         NumericVector xp = rootxylemConductanceProportions(L(j,_), V(j,_));
  //         VCroot_kmax(j,_) = VCroot_kmaxVEC[j]*xp;
  //         //Update Plant_kmax
  //         Plant_kmax[j] = 1.0/((1.0/VCleaf_kmax[j])+(1.0/VCstem_kmax[j])+(1.0/VCroot_kmaxVEC[j]));
        }
      }


      ///// B14b. ESTIMATE LEAF/FINE ROOT BIOMASS balance (g_ind)
      LeafBiomassBalance[j] = leafBiomassIncrement - senescenceLeafLoss;
      FineRootBiomassBalance[j] = finerootBiomassIncrement - senescenceFinerootLoss;

      //INDIVIDUAL-LEVEL OUTPUT
      GROWTHres.PSres.FineRootBiomass[j] = fineRootBiomass[j];
      GROWTHres.PSres.LeafArea[j] = LAexpanded;
      GROWTHres.PSres.SapwoodArea[j] = SA[j];
      GROWTHres.PSres.FineRootArea[j] = fineRootBiomass[j]*specificRootSurfaceArea_c(x.paramsAnatomy.SRL[j], x.paramsAnatomy.FineRootDensity[j])*1e-4;
      GROWTHres.PSres.SapwoodBiomass[j] = sapwoodStructuralBiomass_c(SA[j], H[j], Lj, Vj,WoodDensity[j]);
      GROWTHres.PSres.LeafBiomass[j] = leafStructuralBiomass_c(LAI_expanded[j],N[j],SLA[j]);
      GROWTHres.PSres.RootAreaLeafArea[j] = GROWTHres.PSres.FineRootArea[j]/leafAreaTarget[j];
      GROWTHres.PSres.HuberValue[j] = SA[j]/leafAreaTarget[j];
      GROWTHres.GMres.SAgrowth[j] += deltaSAgrowth[j]; //Store sapwood area growth rate (cm2/day)
      GROWTHres.GMres.LAgrowth[j] += deltaLAgrowth[j];//Store Leaf area growth rate (m2/day)
      GROWTHres.GMres.FRAgrowth[j] = std::accumulate(deltaFRBgrowth.begin(), deltaFRBgrowth.end(), 0.0)*specificRootSurfaceArea_c(x.paramsAnatomy.SRL[j], x.paramsAnatomy.FineRootDensity[j])*1e-4;//Store fine root area growth rate (m2·d-1)
    }
  }
  ///// B15. UPDATE STRUCTURAL VARIABLES /////
  // updateStructuralVariables(x, deltaSAgrowth);

  ///// B16. UPDATE SAPWOOD AREA AND FINEROOT BIOMASS TARGETS AND RECALCULATE CONCENTRATIONS /////
  for(int j=0;j<numCohorts;j++){
    //Update fine root biomass target
    if((LAI_live[j]>0.0) && (N[j]>0.0)) {
      if(x.control.transpirationMode=="Granier") {
        sapwoodAreaTarget[j] = 10000.0*leafAreaTarget[j]/x.paramsAnatomy.Al2As[j];
        fineRootBiomassTarget[j] = (x.paramsAnatomy.Ar2Al[j]*leafAreaTarget[j])/(specificRootSurfaceArea_c(x.paramsAnatomy.SRL[j], x.paramsAnatomy.FineRootDensity[j])*1e-4);
      } else {
        // if(allocationStrategy == "Plant_kmax") {
        //   sapwoodAreaTarget[j] = 10000.0*(leafAreaTarget[j]/x.paramsAnatomy.Al2As[j])*(allocationTarget[j]/x.paramsTranspiration.Plant_kmax[j]);
        // } else if(allocationStrategy =="Al2As") {
        //   sapwoodAreaTarget[j] = 10000.0*leafAreaTarget[j]/x.paramsAnatomy.Al2As[j];
        // }
        // NumericVector VGrhizo_target(nlayers,0.0);
        // for(int s=0;s<nlayers;s++) {
        //   VGrhizo_target[s] = V(j,s)*findRhizosphereMaximumConductance(averageFracRhizosphereResistance*100.0,
        //                         VG_n[s], VG_alpha[s],
        //                                          VCroot_kmaxVEC[j], VCroot_c[j], VCroot_d[j],
        //                                                                                  VCstem_kmax[j], VCstem_c[j], VCstem_d[j],
        //                                                                                                                       VCleaf_kmax[j], VCleaf_c[j], VCleaf_d[j],
        //                                                                                                                                                            log(VGrhizo_kmax(j,s)));
        // }
        // fineRootBiomassTarget[j] = fineRootBiomassPerIndividual(Ksat, VGrhizo_target, LAI_live[j], N[j],
        //                                                         SRL[j], FineRootDensity[j], RLD[j]);
      }
    }

    //RECALCULATE storage concentrations (SA, LA and H may have changed)
    std::vector<double> Vj(nlayers);
    std::vector<double> Lj(nlayers);
    for(int l=0;l<nlayers;l++) {
      Vj[l] = x.belowLayers.V(j,l);
      Lj[l] = x.belowLayers.L(j,l);
    }
    double newVolumeSapwood = sapwoodStorageVolume_c(SA[j], H[j], Lj, Vj, WoodDensity[j], x.paramsAnatomy.conduit2sapwood[j]);
    double newVolumeLeaves = leafStorageVolume_c(LAI_expanded[j],  N[j], SLA[j], LeafDensity[j]);
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
    GROWTHres.PSres.DBH[j] = DBH[j];
    GROWTHres.PSres.Height[j] = H[j];
    GROWTHres.LCBres.SugarLeaf[j] = sugarLeaf[j];
    GROWTHres.LCBres.StarchLeaf[j] = starchLeaf[j];
    GROWTHres.LCBres.SugarSapwood[j] = sugarSapwood[j];
    GROWTHres.LCBres.StarchSapwood[j] = starchSapwood[j];
  }

  ///////////////////////////////////////////
  ///// C. FIRE OCCURRENCE AND BEHAVIOR /////
  ///////////////////////////////////////////
  double fireCombustion = 0.0;// Fire combustion (in C/m2) for stand carbon balance
  // bool fireOccurrence = false;
  // double pfire = meteovec["pfire"];
  // NumericVector fireBehavior = communicationFireHazard();
  // if(R::runif(0.0,1.0) < pfire) {
  //   fccsHazard(fireBehavior, x, meteovec, spwbOutput, slope);
  //   fireOccurrence = true;
  // }
  // 
  // /////////////////////////////////////////////////////////
  // ///// D. PLANT MORTALITY AND FIRE EFFECTS BY COHORT /////
  // /////////////////////////////////////////////////////////
  // for(int j=0;j<numCohorts;j++){
  //   if(N[j] > 0.0) {
  //     
  //     //MORTALITY Death by carbon starvation or dessication
  //     double Ndead_day = 0.0, N_resprouting_stump_day = 0.0, Cover_resprouting_stump_day = 0.0;
  //     bool dynamicCohort = true;
  //     if((ctype[j] == "shrub") && (!shrubDynamics)) dynamicCohort = false;
  //     else if(ctype[j] =="herb") dynamicCohort = false;
  //     double stemSympRWC = NA_REAL;
  //     if(transpirationMode=="Granier") stemSympRWC = symplasticRelativeWaterContent_c(PlantPsi[j], StemPI0[j], StemEPS[j]);
  //     else stemSympRWC = symplasticRelativeWaterContent_c(psiSympStem[j], StemPI0[j], StemEPS[j]);
  //     if(NumericVector::is_na(stemSympRWC)) stop("Missing value for stem symp RWC");
  //     
  //     
  //     //Sapwood sugar relative to equilibrium, indicator of starvation
  //     double relativeSugarSapwood = (sugarSapwood[j]/equilibriumSapwoodSugarConc);
  //     if(dynamicCohort) {
  //       String cause = "undertermined";
  //       //Determine fire severity if fire occurred
  //       double burnRatioLeaves = 0.0, burnRatioBuds = 0.0;
  //       bool abovegroundFireSurvival = true;
  //       if(fireOccurrence) {
  //         double rho_air = airDensity_c(tmax, Patm);
  //         double foliar_factor = leafThermalFactor_c(SLA[j], 130.0, 2500.0);
  //         double Ib_surf = fireBehavior["I_b_surface [kW/m]"];
  //         double t_res_surf = fireBehavior["t_r_surface [s]"];
  //         double t_r_crown = fireBehavior["t_r_crown [s]"];
  //         double fm_dead = fireBehavior["DFMC [%]"];
  //         double Ic_ratio = fireBehavior["Ic_ratio"];
  //         double cbh = H[j]*(1.0 - CR[j]);
  //         double bark_thickness = 1.0;
  //         double Hn_leaves = 0.0, Hn_buds = 0.0;
  //         double xn = 0.0;
  //         if(ctype[j] == "shrub") {
  //           bark_thickness = BTsh[j]*0.1; // from mm to cm
  //         } else {
  //           bark_thickness = Abt[j]*pow(DBH[j],Bbt[j])*0.1;// from mm to cm
  //         }
  //         //Determine foliage/bud burn
  //         if(!NumericVector::is_na(Ib_surf) && !NumericVector::is_na(t_res_surf)) {
  //           Hn_leaves =100.0*necrosisHeight_c(Ib_surf, t_res_surf, foliar_factor, tmax, rho_air); //Necrosis height (cm)
  //           Hn_buds = 100.0*necrosisHeight_c(Ib_surf, t_res_surf, 0.130, tmax, rho_air); //Bud necrosis height (cm)
  //           burnRatioLeaves = leafAreaProportion_c(0.0, Hn_leaves, cbh, H[j]);
  //           burnRatioBuds = leafAreaProportion_c(0.0, Hn_buds, cbh, H[j]);
  //           // Rcout << " tmax " << tmax << " rho_air " << rho_air <<" foliar_factor "<< foliar_factor << " Ib_surf "<< Ib_surf << " t_res_surf " << t_res_surf<< " foliar_factor "<< foliar_factor << " Hn_leaves "<<Hn_leaves << " br_leaves "<< burnRatioLeaves<< " Hn_buds "<<Hn_buds << " br_buds "<< burnRatioBuds<<"\n";
  //         }
  //         //Determine crown fire or torching effects
  //         if(!NumericVector::is_na(Ib_surf)) {
  //           double canopyFMC = (LFMC[j]*(1.0 - StemPLC[j]) + fm_dead*StemPLC[j]);
  //           double Ib_crit = criticalFirelineIntensity_c(cbh/100.0, canopyFMC);
  //           // Rcout << "Ic_ratio "<< Ic_ratio <<" Ib_crit "<<Ib_crit<< " Ib_surf "<< Ib_surf<<"\n";
  //           if((Ic_ratio > 1.0) || (Ib_surf > Ib_crit)) {
  //             burnRatioLeaves = 1.0;
  //             double Tc = necrosisCriticalTemperature_c(t_r_crown, 0.130 , tmax);
  //             if(Tc < 900.0) burnRatioBuds = 1.0;
  //           }
  //           // Rcout << "br_leaves "<< burnRatioLeaves<< " br_buds "<< burnRatioBuds<<"\n";
  //         }
  //         //Surface fire effects on cambium
  //         if(!NumericVector::is_na(Ib_surf) && !NumericVector::is_na(t_res_surf)) {
  //           double bark_diff = barkThermalDiffusivity_c(fm_dead, 500.0, tmax);
  //           xn =  radialBoleNecrosis_c(Ib_surf, t_res_surf, bark_diff, tmax, rho_air);
  //           // Rcout << "xn "<< xn<< " xa "<<bark_thickness<<"\n";
  //         }
  //         
  //         crownBudPercent[j] = crownBudPercent[j]*(1.0 - burnRatioBuds);
  //         abovegroundFireSurvival = (burnRatioBuds < 1.0) && (bark_thickness > xn);
  //         
  //         //If does not survive set ratios to 1.0
  //         if(!abovegroundFireSurvival) {
  //           burnRatioLeaves = 1.0;
  //           burnRatioBuds = 1.0;
  //         }
  //         // Rcout << "abovegroundSurvival "<< abovegroundFireSurvival<< " burnRatioLeaves "<< burnRatioLeaves<< " LAI_burnt_change "<<LAI_burnt_change<<"\n";
  //       }
  //       
  //       
  //       if(abovegroundFireSurvival) {
  //         if(mortalityMode=="whole-cohort/deterministic") {
  //           if((relativeSugarSapwood < mortalityRelativeSugarThreshold) && allowStarvation) {
  //             Ndead_day = N[j];
  //             if(verbose) Rcout<<" [Cohort "<< j<<" died from starvation] ";
  //             cause = "starvation";
  //           } else if( (stemSympRWC < mortalityRWCThreshold) && allowDessication) {
  //             Ndead_day = N[j];
  //             if(verbose) Rcout<<" [Cohort "<< j<<" died from dessication] ";
  //             cause = "dessication";
  //           }
  //         } else {
  //           //Daily basal mortality rate based on constant year probability
  //           double basalMortalityRate = 1.0 - exp(log(1.0 - MortalityBaselineRate[j])/356.0);
  //           //If survival model is available, replace basal mortality value
  //           if(!NumericVector::is_na(SurvivalModelStep[j]) && !NumericVector::is_na(SurvivalB0[j]) && !NumericVector::is_na(SurvivalB1[j])) {
  //             //Probability of dying in model step years
  //             double lp  = SurvivalB0[j] + SurvivalB1[j]*sqrt(treeBasalArea);
  //             double Pmodel = 1.0 - exp(lp)/(1.0 + exp(lp));
  //             //Probability of dying in 1 year
  //             double Pmodel1yr = (1.0- exp(log(1.0-Pmodel)/SurvivalModelStep[j]));
  //             //Daily basal mortality rate based on model
  //             basalMortalityRate = 1.0 - exp(log(1.0 - Pmodel1yr)/356.0);
  //           }
  //           if(allowStarvation) StarvationRate[j] = dailyMortalityProbability_c(relativeSugarSapwood, mortalityRelativeSugarThreshold);
  //           if(allowDessication) DessicationRate[j] = dailyMortalityProbability_c((stemSympRWC + (1.0 - StemPLC[j]))/2.0, mortalityRWCThreshold);
  //           MortalityRate[j] = max(NumericVector::create(basalMortalityRate, DessicationRate[j],  StarvationRate[j]));
  //           if((DessicationRate[j] > basalMortalityRate) && (DessicationRate[j] > StarvationRate[j])) {
  //             cause = "dessication";
  //           } else if((StarvationRate[j] > basalMortalityRate) && (StarvationRate[j] > DessicationRate[j])) {
  //             cause = "starvation";
  //           }
  //           if(NumericVector::is_na(MortalityRate[j])) {
  //             Rcout<< " Basal mortality rate " << basalMortalityRate << " Dessication rate " << DessicationRate[j] << " starvation rate "<< StarvationRate[j]<<"\n";
  //             stop("Missing value for mortality rate");
  //           }
  //           // Rcout<< j << " "<< stemSympRWC<< " "<< DessicationRate[j]<<"\n";
  //           if(mortalityMode =="density/deterministic") {
  //             Ndead_day = N[j]*MortalityRate[j];
  //           } else if(mortalityMode =="whole-cohort/stochastic") {
  //             if(R::runif(0.0,1.0) < MortalityRate[j]) {
  //               Ndead_day = N[j];
  //               Rcout<<" [Cohort "<< j<<" died from " << cause.get_cstring() << "] ";
  //             }
  //           } else if(mortalityMode == "density/stochastic") {
  //             Ndead_day = R::rbinom(round(N[j]), MortalityRate[j]);
  //           }
  //         }
  //       } else { //Cohort burned
  //         Ndead_day = N[j];
  //         cause = "burnt";
  //       }
  //       // Update density and increase the number of dead plants
  //       if(NumericVector::is_na(Ndead_day)) stop("Missing value for Ndead_day");
  //       Ndead_day = std::min(Ndead_day, N[j]);
  //       double Cdead_day = Cover[j]*(Ndead_day/N[j]);
  //       if(cause == "starvation") {
  //         N_starvation[j] = N_starvation[j] + Ndead_day;
  //       } else if(cause == "burnt") {
  //         N_burnt[j] = N_burnt[j] + Ndead_day;
  //         N_resprouting_stump_day = Ndead_day*RespFire[j];
  //       } else if(cause == "dessication") {
  //         N_dessication[j] = N_dessication[j] + Ndead_day;
  //         N_resprouting_stump_day = Ndead_day*RespDist[j];
  //       } else if(ctype[j]=="tree") { // Self-thinning occurring in tree cohorts
  //         if(DBH[j] < IngrowthTreeDBH[j]) {
  //           double b_st = log(RecrTreeDensity[j]/IngrowthTreeDensity[j])/log(RecrTreeDBH[j]/IngrowthTreeDBH[j]);
  //           double a_st = IngrowthTreeDensity[j]/pow(IngrowthTreeDBH[j], b_st);
  //           double N_st = a_st*pow(DBH[j], b_st);
  //           double N_dead_selfthinning = N[j] - std::min(N[j], N_st);
  //           // Rcout<< b_st<< " "<< a_st<< " "<< N_st<< " "<< N_dead_selfthinning<<"\n";
  //           Ndead_day = Ndead_day + N_dead_selfthinning;
  //         }
  //       }
  //       N_resprouting_stumps[j] = N_resprouting_stumps[j] + N_resprouting_stump_day;
  //       // Rcout << Ndead_day <<"\n";
  //       N[j] = N[j] - Ndead_day;
  //       N_dead[j] = N_dead[j] + Ndead_day;
  //       if(ctype[j]=="shrub") {
  //         Cover[j] = std::max(0.0, Cover[j] - Cdead_day);
  //         Cover_dead[j] = Cover_dead[j] + Cdead_day;
  //         if(cause == "starvation") {
  //           Cover_starvation[j] = Cover_starvation[j] + Cdead_day;
  //         } else if(cause == "dessication") {
  //           Cover_dessication[j] = Cover_dessication[j] + Cdead_day;
  //           Cover_resprouting_stump_day =  Cdead_day*RespDist[j];
  //         } else if(cause == "burnt") {
  //           Cover_burnt[j] = Cover_burnt[j] + Cdead_day;
  //           Cover_resprouting_stump_day =  Cdead_day*RespFire[j];
  //         }
  //       }
  //       Cover_resprouting_stumps[j] = Cover_resprouting_stumps[j] + Cover_resprouting_stump_day;
  //       
  //       
  //       //Update LAI live, LAI dead and LAI expanded as a result of density decrease and foliage burn
  //       double LAI_senescence = LAI_expanded[j]*Ndead_day/10000.0;
  //       double LAI_senescence_burned = LAI_senescence*burnRatioLeaves;
  //       double LAI_senescence_notburned = LAI_senescence*(1.0 - burnRatioLeaves);
  //       double LAI_live_burned = (LAI_expanded[j] - LAI_senescence) * burnRatioLeaves;
  //       double LAI_dead_burned = LAI_dead[j]*burnRatioLeaves;
  //       LAI_dead[j] = LAI_dead[j] + LAI_senescence_notburned;
  //       LAI_expanded[j] = LAI_expanded[j] - LAI_senescence - LAI_live_burned;
  //       if(!abovegroundFireSurvival) {
  //         LAI_live[j] = 0.0;
  //       }
  //       
  //       // FOLIAGE BURNING EMISSIONS (in g C / m2)
  //       // from m2/m2 to g C/m2
  //       double leaf_burnt_dry = 1000.0*(LAI_live_burned + LAI_dead_burned + LAI_senescence_burned)/SLA[j];
  //       double leaf_burnt_C = leafCperDry*leaf_burnt_dry;
  //       double twig_burnt_C = WoodC[j]*leaf_burnt_dry/(r635[j] - 1.0);
  //       fireCombustion +=leaf_burnt_C + twig_burnt_C;
  //       
  //       //ADD STANDING DEADWOOD (in g C / m2) according to non-resprouting
  //       double snag_wood_biomass = (Ndead_day/10000.0)*AbovegroundWoodBiomass[j]*WoodC[j]; //Aboveground necromass (includes burning and dessication)
  //       if(cause != "burnt") {
  //         if(ctype[j]=="tree") { // 20% goes to small branches (diameter < 7.5) SHOULD depend on size allometry
  //           Snag_smallbranches[j] += 0.2*snag_wood_biomass;
  //           Snag_largewood[j] += 0.8*snag_wood_biomass;
  //         } else if(ctype[j] =="shrub") { // For shrubs 100% goes to small branches
  //           Snag_smallbranches[j] += snag_wood_biomass;
  //         }
  //       } else { // burnt
  //         if(ctype[j]=="tree") {
  //           fireCombustion += 0.2*snag_wood_biomass; //20% goes to fire emissions
  //           Snag_largewood[j] += 0.8*snag_wood_biomass;
  //         } else if(ctype[j] =="shrub") { // For shrubs 100% goes to fire emissions
  //           fireCombustion += snag_wood_biomass;
  //         }
  //       }
  //       //Fine root and coarse root litter corresponding to non-resprouting trees
  //       if(exportLitter) {
  //         double coarseroot_litter = ((Ndead_day - N_resprouting_stump_day)/10000.0)*BelowgroundWoodBiomass[j]*WoodC[j]; //Only adds non-resprouting death biomass
  //         addCoarseRootLitter(speciesNames[j], coarseroot_litter, internalLitter);
  //         double fineroot_litter = rootCperDry*((Ndead_day - N_resprouting_stump_day)/10000.0)*fineRootBiomass[j];
  //         addFineRootLitter(speciesNames[j], fineroot_litter,
  //                           internalLitter, paramsLitterDecomposition,
  //                           internalSOC);
  //       }
  //     }
  //   }
  // }
  // 
  // 
  // 
  // ////////////////////////////////////////
  // ///// E3. UPDATE PLANT WATER POOLS /////
  // ////////////////////////////////////////
  // if(plantWaterPools) {
  //   NumericVector poolProportions = Rcpp::as<Rcpp::NumericVector>(belowdf["poolProportions"]);
  //   for(int j=0;j<numCohorts;j++) poolProportions[j] = LAI_live[j]/sum(LAI_live);
  //   //Update RHOP
  //   List newRHOP;
  //   if(rhizosphereOverlap=="none") newRHOP = nonoverlapHorizontalProportions(V);
  //   else newRHOP =equaloverlapHorizontalProportions(poolProportions, V);
  //   // else newRHOP = horizontalProportions(poolProportions, CRSV, N, V, widths, rfc);
  //   for(int j=0;j<numCohorts;j++) RHOP[j] = newRHOP[j];
  // }

  ///////////////////////////////////
  ///// F. CARBON DECOMPOSITION /////
  ///////////////////////////////////
  double heterotrophicRespiration = 0.0;
  // if(exportLitter) {
  //   NumericVector baseAnnualRates = control["decompositionAnnualBaseRates"];
  //   double annualTurnoverRate = control["decompositionAnnualTurnoverRate"];
  //   List commDecomp = internalCommunication["decomposition"];
  //   double soilPH = 7.0;
  //   if(soil.containsElementNamed("pH")) {
  //     NumericVector pH = soil["pH"];
  //     if(!NumericVector::is_na(pH[0])) soilPH = pH[0];
  //   }
  //   NumericVector soilTheta = theta(soil, soilFunctions);
  //   NumericVector soilThetaSAT = thetaSAT(soil, soilFunctions);
  //   double vpatm = averageDailyVapourPressure_c(tmin, tmax, rhmin, rhmax);
  //   double tday = averageDaylightTemperature_c(tmin, tmax);
  //   double vpsat = saturationVapourPressure_c(tday);
  //   double rhmean = 100.0*std::max(1.0, vpatm/vpsat);
  //   heterotrophicRespiration = DAYCENTInner(commDecomp, internalSnags, internalLitter, internalSOC,
  //                                           paramsLitterDecomposition,
  //                                           baseAnnualRates, annualTurnoverRate,
  //                                           tday, rhmean,
  //                                           sand[0], clay[0], Tsoil[0], soilTheta[0]/soilThetaSAT[0], soilPH);
  // }
  // 
  // 
  // ///////////////////////////////////////////////
  // ///// G. FIRE EFFECTS ON SNAGS AND LITTER /////
  // ///////////////////////////////////////////////
  // if(fireOccurrence) {
  //   NumericVector li_leaves = internalLitter["Leaves"];
  //   NumericVector li_twigs = internalLitter["Twigs"];
  //   NumericVector li_smallbranches = internalLitter["SmallBranches"];
  //   fireCombustion += sum(li_leaves) + sum(li_twigs) + sum(li_smallbranches);
  //   for(int i=0;i<li_leaves.size();i++) {
  //     li_leaves[i] = 0.0;
  //     li_twigs[i] = 0.0;
  //     li_smallbranches[i] = 0.0;
  //   }
  //   NumericVector sn_smallbranches = internalSnags["SmallBranches"];
  //   fireCombustion += sum(sn_smallbranches);
  //   for(int i=0;i<sn_smallbranches.size();i++) {
  //     sn_smallbranches[i] = 0.0;
  //   }
  // }
  // 
  // ////////////////////////////////////////////
  // ///// E4. CLOSE PLANT BIOMASS BALANCE //////
  // ////////////////////////////////////////////
  // closePlantBiomassBalance(initialFinalCC, plantBiomassBalance, x,
  //                          LabileCarbonBalance, LeafBiomassBalance, FineRootBiomassBalance);
  // DataFrame ccFin_gC_m2 = as<DataFrame>(initialFinalCC["ccFin_gC_m2"]);
  // fillCarbonCompartments(ccFin_gC_m2, x, "gC_m2");
  // NumericVector finalTotalCohortBiomass_gC_m2 = ccFin_gC_m2["TotalBiomass"];
  // // Rcout<< "final cohort gC_m2 = "<< sum(finalTotalCohortBiomass_gC_m2)<<"\n";
  // double finalStand_gC_m2 = sum(finalTotalCohortBiomass_gC_m2);
  // if(exportLitter) {
  //   NumericVector sn_smallBranches = internalSnags["SmallBranches"];
  //   NumericVector sn_largeWood = internalSnags["LargeWood"];
  //   finalStand_gC_m2 += (sum(sn_smallBranches) + sum(sn_largeWood));
  //   NumericVector li_leaves = internalLitter["Leaves"];
  //   NumericVector li_twigs = internalLitter["Twigs"];
  //   NumericVector li_smallBranches = internalLitter["SmallBranches"];
  //   NumericVector li_largeWood = internalLitter["LargeWood"];
  //   NumericVector li_coarseRoots = internalLitter["CoarseRoots"];
  //   NumericVector li_fineRoots = internalLitter["FineRoots"];
  //   finalStand_gC_m2 += (sum(li_leaves) + sum(li_twigs) + sum(li_smallBranches) + sum(li_largeWood) + sum(li_coarseRoots) + sum(li_fineRoots));
  //   finalStand_gC_m2 += sum(internalSOC);
  // }
  // finalStand_gC_m2 += (sum(Snag_smallbranches) + sum(Snag_largewood));
  // double changeStand_gC_m2 = finalStand_gC_m2 - initialStand_gC_m2;
  // 
  ///////////////////////////////////////////////
  ///// H. CLOSE STAND-LEVEL CARBON BALANCE /////
  ///////////////////////////////////////////////
  GROWTHres.standCB.GrossPrimaryProduction = standGrossPrimaryProduction;
  GROWTHres.standCB.MaintenanceRespiration = standMaintenanceRespiration;
  GROWTHres.standCB.SynthesisRespiration = standSynthesisRespiration;
  GROWTHres.standCB.NetPrimaryProduction = standGrossPrimaryProduction - standMaintenanceRespiration - standSynthesisRespiration;
  GROWTHres.standCB.FireCombustion = fireCombustion;
  GROWTHres.standCB.HeterotrophicRespiration = heterotrophicRespiration;
  GROWTHres.standCB.NetEcosystemProduction = GROWTHres.standCB.NetPrimaryProduction - fireCombustion - heterotrophicRespiration; 

  double NEP_gC_m2 = GROWTHres.standCB.NetEcosystemProduction;
  // if(std::abs(changeStand_gC_m2 - NEP_gC_m2) > 0.001) {
  //   // Rcout << " Stand biomass balance not closed!\n";
  //   // Rcout<< "GPP " << standGrossPrimaryProduction<<" MR "<<standMaintenanceRespiration<< " SR "<<standSynthesisRespiration<< " HR "<< heterotrophicRespiration << " FC "<< fireCombustion<<"\n";
  //   // Rcout<< "NEP: " << NEP_gC_m2<< " Mass change: " << changeStand_gC_m2<< " dif: "<< (changeStand_gC_m2 - NEP_gC_m2)<< "\n";
  // }
}

void growthDay_inner_c(GROWTH_RESULT& GROWTHres, GROWTHCommunicationStructures& GROWTHcomm, ModelInput& x, 
                       std::string date,
                       WeatherInputVector meteovec, 
                       double latitude, double elevation, double slope, double aspect,
                       const double runon, 
                       const std::vector<double>& lateralFlows, const double waterTableDepth) {
  
  if(std::isnan(meteovec.prec)) throw medfate::MedfateInternalError("Missing precipitation value");
  if(std::isnan(meteovec.tmin)) throw medfate::MedfateInternalError("Missing minimum temperature value");
  if(std::isnan(meteovec.tmax)) throw medfate::MedfateInternalError("Missing maximum temperature value");
  if(meteovec.tmin > meteovec.tmax) {
    double swap = meteovec.tmin;
    meteovec.tmin = meteovec.tmax;
    meteovec.tmax = swap;
  }
  if(std::isnan(meteovec.rhmax)) {
    meteovec.rhmax = 100.0;
  }
  if(std::isnan(meteovec.rhmin)) {
    double vp_tmin = saturationVapourPressure_c(meteovec.tmin);
    double vp_tmax = saturationVapourPressure_c(meteovec.tmax);
    meteovec.rhmin = std::min(meteovec.rhmax, 100.0*(vp_tmin/vp_tmax));
  }
  if(meteovec.rhmin > meteovec.rhmax) {
    // warning("rhmin > rhmax. Swapping values.");
    double swap = meteovec.rhmin;
    meteovec.rhmin = meteovec.rhmax;
    meteovec.rhmax = swap;
  }
  if(std::isnan(meteovec.wind)) meteovec.wind = x.control.weather.defaultWindSpeed; 
  if(meteovec.wind<0.1) meteovec.wind = 0.1; //Minimum windspeed abovecanopy
  
  if(std::isnan(meteovec.Catm)) meteovec.Catm = x.control.weather.defaultCO2;
  
  int month = std::atoi(date.substr(5,2).c_str());
  int J = julianDay_c(std::atoi(date.substr(0, 4).c_str()),std::atoi(date.substr(5,2).c_str()),std::atoi(date.substr(8,2).c_str()));
  double delta = solarDeclination_c(J);
  double solarConstant = solarConstant_c(J);
  double latrad = latitude * (M_PI/180.0);
  if(std::isnan(aspect)) aspect = 0.0;
  if(std::isnan(slope)) slope = 0.0;
  double asprad = aspect * (M_PI/180.0);
  double slorad = slope * (M_PI/180.0);
  double photoperiod = daylength_c(latrad, 0.0, 0.0, delta);
  if(std::isnan(meteovec.tday)) {
    meteovec.tday = averageDaylightTemperature_c(meteovec.tmin, meteovec.tmax);
  }
  if(std::isnan(meteovec.rad)) {
    // warning("Estimating solar radiation");
    double vpa = averageDailyVapourPressure_c(meteovec.tmin, meteovec.tmax, meteovec.rhmin, meteovec.rhmax);
    meteovec.rad = RDay_c(solarConstant, latrad, elevation,
                          slorad, asprad, delta, meteovec.tmax -meteovec.tmin, meteovec.tmax-meteovec.tmin,
                          vpa, meteovec.prec);
  }
  if(std::isnan(meteovec.pet)) {
    meteovec.pet = PenmanPET_c(latrad, elevation, slorad, asprad, J, 
                               meteovec.tmin, meteovec.tmax, meteovec.rhmin, meteovec.rhmax, meteovec.rad, meteovec.wind);
  }
  
  //Derive doy from date  
  int J0101 = julianDay_c(std::atoi(date.substr(0, 4).c_str()),1,1);
  int doy = J - J0101+1;
  
  std::vector<double> defaultRainfallIntensityPerMonth = x.control.weather.defaultRainfallIntensityPerMonth;
  if(std::isnan(meteovec.rint)) meteovec.rint = rainfallIntensity_c(month, meteovec.prec, defaultRainfallIntensityPerMonth);
  
  bool leafPhenology = x.control.phenology.leafPhenology;
  
  //Update phenology
  if(leafPhenology) {
    updatePhenology_c(x, doy, photoperiod, meteovec.tday);
    updateLeaves_c(x, meteovec.wind, false);
  }
  
  meteovec.tminPrev = meteovec.tmin;
  meteovec.tmaxPrev = meteovec.tmax;
  meteovec.tminNext = meteovec.tmin;
  
  // Rcpp::Rcout << "about to enter growth day_private_c\n";
  growthDay_private_c(GROWTHres, GROWTHcomm, x, 
                      meteovec, 
                      latitude, elevation, slope, aspect,
                      solarConstant, delta,
                      runon, 
                      lateralFlows, waterTableDepth);
  
}
