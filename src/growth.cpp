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
#include <Rcpp.h>
#include <meteoland.h>
using namespace Rcpp;

const double Q10_resp = 2.0;


/**
 * floem flow (Holtta et al. 2017)
 *  psiUpstream, psiDownstream - water potential upstream (leaves)  and downstream
 *  concUpstream, concDownstream - sugar concentration upstream (leaves) and downstream (stem)
 *  k_f - floem conductance per leaf area basis (l*m-2*MPa-1*s-1)
 *  
 *  out mol*s-1*m-2 (flow per leaf area basis)
 */
double floemFlow(double psiUpstream, double psiDownstream,
                 double concUpstream, double concDownstream,
                 double temp, double k_f, double nonSugarConc) {
  double turgor_up = turgor(psiUpstream, concUpstream, temp, nonSugarConc);
  double turgor_down = turgor(psiDownstream, concDownstream, temp, nonSugarConc);
  if(temp < 0.0) k_f = 0.0; // No floem flow if temperature below zero
  double relVisc = relativeSapViscosity((concUpstream+concDownstream)/2.0, temp);
  if(turgor_up>turgor_down) {
    return(k_f*concUpstream*(turgor_up - turgor_down)/relVisc);
  } else {
    return(k_f*concDownstream*(turgor_up - turgor_down)/relVisc);
  }
}
// // [[Rcpp::export("growth_dailyFloemFlow")]]
// NumericMatrix dailyFloemFlow(List x, List spwbOut, 
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
//       ff(c,s) = floemFlow(psiUp, psiDown, concLeaf[c], concSapwood[c], Tcan[s], k_floem, nonSugarConc)*3600.0; //flow as mol per hour and leaf area basis
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


List growthDay1(List x, List soil, double tday, double pet, double prec, double er, double runon=0.0, 
              double rad = NA_REAL, double elevation = NA_REAL, bool verbose = false) {
  
  
  //Soil-plant water balance
  List spwbOut = spwbDay1(x, soil, tday, pet, prec, er,runon,
                          rad, elevation, verbose);
  
  //Control params
  List control = x["control"];  
  
  String transpirationMode = control["transpirationMode"];
  String allocationStrategy = control["allocationStrategy"];
  String cavitationRefill = control["cavitationRefill"];
  bool plantWaterPools = control["plantWaterPools"];
  double nonSugarConcentration = control["nonSugarConcentration"];
  NumericVector equilibriumOsmoticConcentration  = control["equilibriumOsmoticConcentration"];
  double equilibriumLeafTotalConc = equilibriumOsmoticConcentration["leaf"];
  double equilibriumSapwoodTotalConc = equilibriumOsmoticConcentration["sapwood"];
  NumericVector minimumSugarForGrowth = control["minimumSugarForGrowth"];
  double minimumSugarGrowthLeaves = minimumSugarForGrowth["leaf"];
  double minimumStarchGrowthSapwood = minimumSugarForGrowth["sapwood"];
  double minimumSugarGrowthFineRoots = minimumSugarForGrowth["fineroot"];
  NumericVector respirationRates = control["respirationRates"];
  double leaf_RR = respirationRates["leaf"];
  double sapwood_RR = respirationRates["sapwood"];
  double fineroot_RR = respirationRates["fineroot"];
  NumericVector turnoverRates = control["turnoverRates"];
  double dailySapwoodTurnoverProportion = turnoverRates["sapwood"];
  double dailyFineRootTurnoverProportion = turnoverRates["fineroot"];
  NumericVector constructionCosts = control["constructionCosts"];
  double leaf_CC = constructionCosts["leaf"];
  double sapwood_CC = constructionCosts["sapwood"];
  NumericVector maximumRGR = control["maximumRelativeGrowthRates"];
  double RGRleafmax = maximumRGR["leaf"];
  double RGRfinerootmax = maximumRGR["fineroot"];
  
  //Cohort info
  DataFrame cohorts = Rcpp::as<Rcpp::DataFrame>(x["cohorts"]);
  NumericVector SP = cohorts["SP"];
  int numCohorts = SP.size();
  
  //Soil
  NumericVector dVec = soil["dVec"];
  int nlayers = dVec.size();
  
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
  StringVector Status = Rcpp::as<Rcpp::StringVector>(above["Status"]);
  
  //Belowground parameters  
  DataFrame belowdf = Rcpp::as<Rcpp::DataFrame>(x["below"]);
  NumericVector Z95 = Rcpp::as<Rcpp::NumericVector>(belowdf["Z95"]);
  NumericVector Z50 = Rcpp::as<Rcpp::NumericVector>(belowdf["Z50"]);
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
  
  DataFrame internalAllocation = Rcpp::as<Rcpp::DataFrame>(x["internalAllocation"]);
  NumericVector allocationTarget = internalAllocation["allocationTarget"];
  NumericVector leafAreaTarget = internalAllocation["leafAreaTarget"];
  
  DataFrame internalPhenology = Rcpp::as<Rcpp::DataFrame>(x["internalPhenology"]);
  LogicalVector leafUnfolding = internalPhenology["leafUnfolding"];
  LogicalVector budFormation = internalPhenology["budFormation"];
  
  List stand = spwbOut["Stand"];
  DataFrame Plants = Rcpp::as<Rcpp::DataFrame>(spwbOut["Plants"]);
  NumericVector Ag = Plants["Photosynthesis"];

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
  //Growth parameters
  DataFrame paramsGrowth = Rcpp::as<Rcpp::DataFrame>(x["paramsGrowth"]);
  NumericVector WoodC = Rcpp::as<Rcpp::NumericVector>(paramsGrowth["WoodC"]);
  NumericVector RGRsapwoodmax = Rcpp::as<Rcpp::NumericVector>(paramsGrowth["RGRsapwoodmax"]);
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
  NumericVector CarbonBalance(numCohorts,0.0);
  NumericVector MaintenanceRespiration(numCohorts,0.0);
  NumericVector GrowthRespiration(numCohorts,0.0);
  NumericVector LabileMassLeaf(numCohorts,0.0), LabileMassSapwood(numCohorts,0.0);
  NumericVector PlantSugarTransport(numCohorts,0.0), PlantSugarLeaf(numCohorts,0.0), PlantStarchLeaf(numCohorts,0.0);
  NumericVector PlantSugarSapwood(numCohorts,0.0), PlantStarchSapwood(numCohorts,0.0);
  NumericVector SapwoodArea(numCohorts,0.0), LeafArea(numCohorts,0.0);
  NumericVector SAgrowth(numCohorts,0.0);
  NumericVector LAgrowth(numCohorts,0.0);
  NumericVector GrossPhotosynthesis(numCohorts,0.0);
  NumericVector PlantLAIdead(numCohorts,0.0), PlantLAIlive(numCohorts,0.0),PlantLAIexpanded(numCohorts,0.0);
  
  //Storage volume and maximum starch capacity for leaves and sapwood  
  NumericVector Volume_leaves(numCohorts,0.0);
  NumericVector Volume_sapwood(numCohorts,0.0);
  NumericVector Starch_max_leaves(numCohorts,0.0);
  NumericVector Starch_max_sapwood(numCohorts,0.0);
  NumericVector B_struct_leaves(numCohorts,0.0);
  NumericVector B_struct_sapwood(numCohorts,0.0);
  NumericVector B_struct_fineroots(numCohorts,0.0);
  
  double minimumLeafSugarConc = equilibriumLeafTotalConc - nonSugarConcentration;
  double minimumSapwoodSugarConc = equilibriumSapwoodTotalConc - nonSugarConcentration;
  
  double rleafcellmax = relative_expansion_rate(0.0 ,25, -2.0,0.5,0.05,5.0);
  
  //3. Carbon balance and growth
  for(int j=0;j<numCohorts;j++){
    if(Status[j]=="alive") {
      double costPerLA = 1000.0*leaf_CC/SLA[j]; // Construction cost in g gluc · m-2 of leaf area
      double costPerSA = sapwood_CC*(H[j]+(Z95[j]/10.0))*WoodDensity[j];  //Construction cost in g gluc ·cm-2 of sapwood
      
      Volume_leaves[j] = leafStorageVolume(LAI_expanded[j],  N[j], SLA[j], LeafDensity[j]);
      Volume_sapwood[j] = sapwoodStorageVolume(SA[j], H[j], L(j,_), V(j,_),WoodDensity[j], 0.5);
      Starch_max_leaves[j] = leafStarchCapacity(LAI_expanded[j],  N[j], SLA[j], 0.3)/Volume_leaves[j];
      Starch_max_sapwood[j] = sapwoodStarchCapacity(SA[j], H[j], L(j,_), V(j,_),WoodDensity[j], 0.2)/Volume_sapwood[j];
      B_struct_leaves[j] = leafStructuralBiomass(LAI_expanded[j],N[j],SLA[j]);
      B_struct_sapwood[j] = sapwoodStructuralLivingBiomass(SA[j], H[j], L(j,_), V(j,_), WoodDensity[j], 0.5);
      B_struct_fineroots[j] = B_struct_leaves[j]/2.0; //TO BE CHANGED
      
      double labileMassLeafIni = (sugarLeaf[j]+starchLeaf[j])*(glucoseMolarMass*Volume_leaves[j]);
      double labileMassSapwoodIni = (sugarSapwood[j]+starchSapwood[j])*(glucoseMolarMass*Volume_sapwood[j]);
      
      double B_total = (B_struct_leaves[j] + B_struct_sapwood[j] + B_struct_fineroots[j]+labileMassSapwoodIni+labileMassLeafIni);
      
      // Rcout << j << " Lvol: "<< Volume_leaves[j] << " Svol: "<<Volume_sapwood[j]<< " LStarchMax: "<<Starch_max_leaves[j]
            // << " SStarchMax: "<<Starch_max_sapwood[j]<< " Bleaf "<< B_struct_leaves[j]<< " Bsap "<< B_struct_sapwood[j]<< " Bfr "<< B_struct_fineroots[j]<<"\n";
      
      double LAexpanded = leafArea(LAI_expanded[j], N[j]);
      double LAlive = leafArea(LAI_live[j], N[j]);
      double LAdead = leafArea(LAI_dead[j], N[j]);
      double LAlive_ini = LAlive;
      
      double leafRespDay = 0.0;
      // double sfrRespDay = 0.0;
      
      //Set target leaf area if bud formation is allowed
      if(budFormation[j]) {
        leafAreaTarget[j] = (SA[j]/10000.0)*allocationTarget[j];
      }
      
      //Carbon balance for labile carbon of leaves and stems
      double leafSugarMass = sugarLeaf[j]*(Volume_leaves[j]*glucoseMolarMass);
      double sapwoodSugarMass = sugarSapwood[j]*(Volume_sapwood[j]*glucoseMolarMass);

      //Respiratory biomass (g dw · ind-1)
      double B_resp_leaves = B_struct_leaves[j] + leafSugarMass;
      double B_resp_sapwood = B_struct_sapwood[j] + sapwoodSugarMass;
      double B_resp_fineroots = B_struct_fineroots[j];
      double QR = qResp(tday);
      if(LAlive>0.0) leafRespDay = B_resp_leaves*leaf_RR*QR;
      double sapwoodResp = B_resp_sapwood*sapwood_RR*QR;
      double finerootResp = B_resp_fineroots*fineroot_RR*QR;

      double leafAgG = 0.0;
      if(LAlive>0.0) {
        //gross fotosynthesis
        double leafAgC = Ag[j]/(N[j]/10000.0); //Translate g C · m-2 · h-1 to g C · h-1
        leafAgG = leafAgC*(glucoseMolarMass/(carbonMolarMass*6.0)); // from g C· h-1 to g gluc · h-1
        
        //Update output values
        GrossPhotosynthesis[j] = leafAgG/B_total; //Ag in g gluc · gdry-1 
      }
      MaintenanceRespiration[j] += (leafRespDay+sapwoodResp+finerootResp)/B_total; 
      
      
      //Growth
      double growthCostLA = 0.0;
      double growthCostSA = 0.0;
      double deltaLAgrowth = 0.0;
      double deltaSAgrowth = 0.0;
      List ring = ringList[j];
      grow_ring(ring, PlantPsi[j] ,tday, 10.0);
      double rleafcell = relative_expansion_rate(PlantPsi[j] ,tday, LeafPI0[j],0.5,0.05,5.0);
      
      
      if(leafUnfolding[j]) {
        double deltaLApheno = std::max(leafAreaTarget[j] - LAlive, 0.0);
        double deltaLAsink = std::min(deltaLApheno, SA[j]*RGRleafmax*(rleafcell/rleafcellmax));
        if(LAlive>0.0) {
          double deltaLAavailable = std::max(0.0,((sugarLeaf[j] - minimumSugarGrowthLeaves)*(glucoseMolarMass*Volume_leaves[j]))/costPerLA);
          deltaLAgrowth += std::min(deltaLAsink, deltaLAavailable);
          growthCostLA += deltaLAgrowth*costPerLA;
        } else { //Grow at expense of stem sugar
          double deltaLAavailable = std::max(0.0,((sugarSapwood[j] - minimumSugarGrowthLeaves)*(glucoseMolarMass*Volume_sapwood[j]))/costPerLA);
          deltaLAgrowth += std::min(deltaLAsink, deltaLAavailable);
          growthCostSA += deltaLAgrowth*costPerLA;
        }
      }

      if(LAlive > 0.0) {
        NumericVector SAring = ring["SA"];
        double deltaSAring = 0.0;
        if(SAring.size()==1) deltaSAring = SAring[0];
        else deltaSAring = SAring[SAring.size()-1] - SAring[SAring.size()-2];
        double RGRcellmax = (2e-8/SA[j]);
        double deltaSAsink = (1e-8*(deltaSAring/10.0))*(RGRsapwoodmax[j]/RGRcellmax); //Correction for the difference in the number of cells
        double deltaSAavailable = std::max(0.0,((starchSapwood[j]- minimumStarchGrowthSapwood)*(glucoseMolarMass*Volume_sapwood[j]))/costPerSA);
        deltaSAgrowth = std::min(deltaSAsink, deltaSAavailable);
        // Rcout<< SAring.size()<<" " <<j<< " "<< PlantPsi[j]<< " "<< LeafPI0[j]<<" dSAring "<<deltaSAring<< " dSAsink "<< deltaSAsink<<" dSAgrowth "<< deltaSAgrowth<<"\n";
        growthCostSA = deltaSAgrowth*costPerSA; //increase cost (may be non zero if leaf growth was charged onto sapwood)
      }      
      
      GrowthRespiration[j] +=(growthCostLA + growthCostSA)/B_total; //growth cost in g gluc · gdry-1
      //Instantaneous carbon balance
      CarbonBalance[j] +=GrossPhotosynthesis[j] - MaintenanceRespiration[j] - GrowthRespiration[j];
      
      //sugar mass balance
      double leafSugarMassDelta = leafAgG - leafRespDay - growthCostLA;
      double sapwoodSugarMassDelta =  - finerootResp;
      double sapwoodStarchMassDelta = - sapwoodResp - growthCostSA;
        
      sugarSapwood[j] += sapwoodSugarMassDelta/(Volume_sapwood[j]*glucoseMolarMass);
      starchSapwood[j] += sapwoodStarchMassDelta/(Volume_sapwood[j]*glucoseMolarMass);
      
      if(LAlive>0.0) {
        sugarLeaf[j] += leafSugarMassDelta/(Volume_leaves[j]*glucoseMolarMass);
        //floem transport to make sugar concentrations equal     
        double ff = (sugarLeaf[j]-sugarSapwood[j])/2.0; 
        sugarLeaf[j] -=ff;
        PlantSugarTransport[j] = (ff*Volume_leaves[j])/(3600.0*24.0); //mol · s-1
        sugarSapwood[j] +=(Volume_leaves[j]/Volume_sapwood[j])*ff;
        double conversionLeaf = std::max(-starchLeaf[j], sugarLeaf[j] - minimumLeafSugarConc);
        starchLeaf[j] +=conversionLeaf;
        sugarLeaf[j] -=conversionLeaf;
      }
      double conversionSapwood = std::max(-starchSapwood[j], sugarSapwood[j] - minimumSapwoodSugarConc);
      starchSapwood[j] +=conversionSapwood;
      sugarSapwood[j] -=conversionSapwood;

      //Excess is lost
      starchLeaf[j] = std::min(starchLeaf[j], Starch_max_leaves[j]);
      starchSapwood[j] = std::min(starchSapwood[j], Starch_max_sapwood[j]);
      
      if(sugarLeaf[j] < 0.0) { //Leaf senescense due to C starvation
        double respirationExcess = -sugarLeaf[j]*(Volume_leaves[j]*glucoseMolarMass); //g gluc
        double propExcess = respirationExcess/leafRespDay; //day
        // Rcout<< j <<" Excess respiration: " << respirationExcess << " Prop:"<< propExcess<< " LAlive " << LAlive << " LAlivenew "<< LAlive*(1.0 - propExcess) <<"\n";
        LAdead = LAdead + LAexpanded*propExcess;
        LAexpanded = LAexpanded*(1.0 - propExcess);
        LAlive = LAlive*(1.0 - propExcess);
        if(LAlive<0.0001) {
          LAlive = 0.0;
          LAexpanded = 0.0;
        }
        sugarLeaf[j] = 0.0;
      }
      
      //Leaf growth
      LAlive += deltaLAgrowth; //Update leaf area
      LAexpanded +=deltaLAgrowth;
      // Rcout<< j<<" "<< growthCostLA << " "<<deltaLAgrowth<< " "<<LAlive<<"\n";
      //Leaf senescence
      double propLeafSenescence = 0.0;
      //Leaf senescence due to age (Ca+ accumulation) only in evergreen species
      if(phenoType[j] == "oneflush-evergreen" || phenoType[j] == "progressive-evergreen") {
        propLeafSenescence = (1.0/(365.25*leafDuration[j]));
      }
      //Leaf senescence due to drought 
      double LAplc = std::min(LAlive, (1.0 - StemPLC[j])*leafAreaTarget[j]);
      if(LAplc<LAlive) {
        // Rcout<<j<< " "<< LAplc<< " "<< LAlive<<"\n";
        propLeafSenescence = std::max((LAlive-LAplc)/LAlive, propLeafSenescence); 
      }
      // if(LeafRWC[j] < 0.5) {
      //   double k = -5.0;
      //   propLeafSenescence = std::min(propLeafSenescence,
      //                                 std::max(0.0,(exp(k*LeafRWC[j])-exp(k*0.5))/(1.0-exp(k*0.5))));
      // }
      double LA_exp_prev= LAexpanded; //Store previous value
      LAdead += LAexpanded*propLeafSenescence;
      LAexpanded = LAexpanded*(1.0 - propLeafSenescence); //Update expanded leaf area
      LAlive = LAlive*(1.0 - propLeafSenescence); //Update expanded leaf area
      
      //SA growth senescense
      double SAprev = SA[j];
      double deltaSAturnover = (dailySapwoodTurnoverProportion/(1.0+15.0*exp(-0.01*H[j])))*SA[j];
      SA[j] = SA[j] - deltaSAturnover; //Update sapwood area
      //SA growth     
      SA[j] += deltaSAgrowth; //Update sapwood area
      //Decrease PLC due to new SA growth
      if(cavitationRefill=="growth") StemPLC[j] = std::max(0.0, StemPLC[j] - (deltaSAgrowth/SA[j]));
      
      //Death by carbon starvation or dessication
      if((sugarSapwood[j]<0.0) || (StemPLC[j] > 0.5)) {
        LAdead = LAlive;
        LAlive = 0.0;
        LAexpanded = 0.0;
        if(sugarSapwood[j]<0.0) Status(j) = "starvation";
        else if(StemPLC[j]>0.5) Status(j) = "dessication";
        Rcout<<" [Cohort "<< j<<" died from " << Status(j)<<"] ";
      }
      
      
      //Update LAI
      LAI_live[j] = LAlive*N[j]/10000.0;
      LAI_expanded[j] = LAexpanded*N[j]/10000.0;
      LAI_dead[j] = LAdead*N[j]/10000.0;
      
      //Update Huber value
      if(LAlive>0.0) {
        Al2As[j] = (LAlive)/(SA[j]/10000.0);
      }
      //Update leaf and stem osmotic water potential at full turgor
      // LeafPI0[j] = osmoticWaterPotential(sugarLeaf[j], 20.0, nonSugarConcentration); //Osmotic potential at full turgor assuming RWC = 1 and 20ºC
      // StemPI0[j] = osmoticWaterPotential(sugarSapwood[j], 20.0, nonSugarConcentration);
      // 
      //Output variables
      PlantSugarLeaf[j] = sugarLeaf[j];
      PlantStarchLeaf[j] = starchLeaf[j];
      PlantSugarSapwood[j] = sugarSapwood[j];
      PlantStarchSapwood[j] = starchSapwood[j];
      SapwoodArea[j] = SA[j];
      SAgrowth[j] += deltaSAgrowth/SA[j]; //Store sapwood area growth rate (cm2/cm2)
      LAgrowth[j] += deltaLAgrowth/SA[j];//Store Leaf area growth rate in relation to sapwood area (m2/cm2)
      LeafArea[j] = LAexpanded;
      PlantLAIlive[j] = LAI_live[j];
      PlantLAIexpanded[j] = LAI_expanded[j];
      PlantLAIdead[j] = LAI_dead[j];
      LabileMassLeaf[j] = (sugarLeaf[j]+starchLeaf[j])*(glucoseMolarMass*Volume_leaves[j]);
      LabileMassSapwood[j] = (sugarSapwood[j]+starchSapwood[j])*(glucoseMolarMass*Volume_sapwood[j]);
      
      //Carbon balance check
      // double sugarTransportMass = PlantSugarTransport[j]*glucoseMolarMass;
      // double sumLeaf = GrossPhotosynthesis[j] - LeafMaintenanceRespiration[j] - LeafGrowthRespiration[j] - sugarTransportMass;
      // double sumSapwood = sugarTransportMass - SapwoodMaintenanceRespiration[j] - FineRootMaintenanceRespiration[j] - SapwoodGrowthRespiration[j] - FineRootGrowthRespiration[j];
      // Rcout<<j<<" CBLeaf "<< sumLeaf << " ChLabLeaf: "<< (LabileMassLeaf[j] - labileMassLeafIni);
      // Rcout<<" CBSapwood "<< sumSapwood << " ChLabSapwood: "<< (LabileMassSapwood[j] - labileMassSapwoodIni)<<"\n";
    }
  }
  //Update pool proportions??
  if(plantWaterPools) {
    NumericVector poolProportions = Rcpp::as<Rcpp::NumericVector>(belowdf["poolProportions"]);
    // for(int j=0;j<numCohorts;j++) poolProportions[j] = LAI_live[j]/sum(LAI_live);
  }
  
  //Needed with string vectors
  above["Status"] = Status;
  

  DataFrame plantCarbonBalance = DataFrame::create(_["GrossPhotosynthesis"] = GrossPhotosynthesis,
                                                   _["MaintenanceRespiration"] = MaintenanceRespiration,
                                                   _["GrowthRespiration"] = GrowthRespiration,
                                                   _["CarbonBalance"] = CarbonBalance,
                                                   _["SugarLeaf"] = PlantSugarLeaf,
                                                   _["StarchLeaf"] = PlantStarchLeaf,
                                                   _["SugarSapwood"] = PlantSugarSapwood,
                                                   _["StarchSapwood"] = PlantStarchSapwood,
                                                   _["SugarTransport"] = PlantSugarTransport,
                                                   _["LabileMassLeaf"] = LabileMassLeaf,
                                                   _["LabileMassSapwood"] = LabileMassSapwood,
                                                   _["StemPI0"] = clone(StemPI0), //Store a copy of the current osmotic potential at full turgor
                                                   _["LeafPI0"] = clone(LeafPI0));
  plantCarbonBalance.attr("row.names") = above.attr("row.names");
  
  DataFrame plantGrowth = List::create(
    _["SapwoodArea"] = SapwoodArea,
    _["LeafArea"] = LeafArea,
    _["SAgrowth"] = SAgrowth,
    _["LAgrowth"] = LAgrowth
  );
  plantGrowth.attr("row.names") = above.attr("row.names");
  
  List l = List::create(_["cohorts"] = clone(cohorts),
                        _["WaterBalance"] = spwbOut["WaterBalance"], 
                        _["Soil"] = spwbOut["Soil"], 
                        _["Stand"] = spwbOut["Stand"], 
                        _["Plants"] = spwbOut["Plants"],
                        _["PlantCarbonBalance"] = plantCarbonBalance,
                        _["PlantGrowth"] = plantGrowth);
  l.attr("class") = CharacterVector::create("growth_day","list");
  return(l);
}




List growthDay2(List x, List soil, double tmin, double tmax, double tminPrev, double tmaxPrev, double tminNext, 
                double rhmin, double rhmax, double rad, double wind, 
                double latitude, double elevation, double slope, double aspect,
                double solarConstant, double delta, 
                double prec, double pet, double er, double runon=0.0, bool verbose = false) {
  
  //1. Soil-plant water balance
  List spwbOut = spwbDay2(x, soil, tmin, tmax, tminPrev, tmaxPrev, tminNext,
                          rhmin, rhmax, rad, wind, 
                          latitude, elevation, slope, aspect,
                          solarConstant, delta, 
                          prec, pet, er, runon, verbose);
  

  //2. Retrieve state
  
  //Control params
  List control = x["control"];  
  
  double tday = meteoland::utils_averageDaylightTemperature(tmin, tmax);
  
  String transpirationMode = control["transpirationMode"];
  String allocationStrategy = control["allocationStrategy"];
  String cavitationRefill = control["cavitationRefill"];
  bool plantWaterPools = control["plantWaterPools"];
  bool taper = control["taper"];
  bool nonStomatalPhotosynthesisLimitation = control["nonStomatalPhotosynthesisLimitation"];
  double averageFracRhizosphereResistance = control["averageFracRhizosphereResistance"];
  double floemConductanceFactor = control["floemConductanceFactor"];
  double nonSugarConcentration = control["nonSugarConcentration"];
  NumericVector equilibriumOsmoticConcentration  = control["equilibriumOsmoticConcentration"];
  double equilibriumLeafTotalConc = equilibriumOsmoticConcentration["leaf"];
  double equilibriumSapwoodTotalConc = equilibriumOsmoticConcentration["sapwood"];
  NumericVector minimumSugarForGrowth = control["minimumSugarForGrowth"];
  double minimumSugarGrowthLeaves = minimumSugarForGrowth["leaf"];
  double minimumStarchGrowthSapwood = minimumSugarForGrowth["sapwood"];
  double minimumSugarGrowthFineRoots = minimumSugarForGrowth["fineroot"];
  NumericVector respirationRates = control["respirationRates"];
  double leaf_RR = respirationRates["leaf"];
  double sapwood_RR = respirationRates["sapwood"];
  double fineroot_RR = respirationRates["fineroot"];
  NumericVector turnoverRates = control["turnoverRates"];
  double dailySapwoodTurnoverProportion = turnoverRates["sapwood"];
  double dailyFineRootTurnoverProportion = turnoverRates["fineroot"];
  NumericVector constructionCosts = control["constructionCosts"];
  double leaf_CC = constructionCosts["leaf"];
  double sapwood_CC = constructionCosts["sapwood"];
  double fineroot_CC = constructionCosts["fineroot"];
  NumericVector maximumRGR = control["maximumRelativeGrowthRates"];
  double RGRleafmax = maximumRGR["leaf"];
  double RGRfinerootmax = maximumRGR["fineroot"];
  
  //Soil info
  NumericVector Ksat = soil["Ksat"];
  NumericVector dVec = soil["dVec"];
  NumericVector rfc = soil["rfc"];
  NumericVector VG_n = soil["VG_n"];
  NumericVector VG_alpha = soil["VG_alpha"];
  NumericVector Tsoil = soil["Temp"];
  
  //Cohort info
  DataFrame cohorts = Rcpp::as<Rcpp::DataFrame>(x["cohorts"]);
  NumericVector SP = cohorts["SP"];
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
  StringVector Status = Rcpp::as<Rcpp::StringVector>(above["Status"]);
  
  
  
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
  
  DataFrame internalAllocation = Rcpp::as<Rcpp::DataFrame>(x["internalAllocation"]);
  NumericVector allocationTarget = internalAllocation["allocationTarget"];
  NumericVector leafAreaTarget = internalAllocation["leafAreaTarget"];
  NumericVector fineRootBiomassTarget = internalAllocation["fineRootBiomassTarget"];

  DataFrame internalPhenology = Rcpp::as<Rcpp::DataFrame>(x["internalPhenology"]);
  LogicalVector leafUnfolding = internalPhenology["leafUnfolding"];
  LogicalVector budFormation = internalPhenology["budFormation"];
  
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
  
  //Growth parameters
  DataFrame paramsGrowth = Rcpp::as<Rcpp::DataFrame>(x["paramsGrowth"]);
  NumericVector WoodC = Rcpp::as<Rcpp::NumericVector>(paramsGrowth["WoodC"]);
  NumericVector RGRsapwoodmax = Rcpp::as<Rcpp::NumericVector>(paramsGrowth["RGRsapwoodmax"]);
  
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
  NumericMatrix CarbonBalanceInst(numCohorts, numSteps);  
  NumericMatrix GrossPhotosynthesisInst(numCohorts, numSteps);  
  NumericMatrix MaintenanceRespirationInst(numCohorts, numSteps);  
  NumericMatrix GrowthRespirationInst(numCohorts, numSteps);  
  NumericMatrix PlantSugarTransportInst(numCohorts, numSteps);
  NumericMatrix PlantSugarLeafInst(numCohorts, numSteps), PlantStarchLeafInst(numCohorts, numSteps);
  NumericMatrix PlantSugarSapwoodInst(numCohorts, numSteps), PlantStarchSapwoodInst(numCohorts, numSteps);
  
  //Daily output vectors
  NumericVector CarbonBalance(numCohorts,0.0);
  NumericVector MaintenanceRespiration(numCohorts,0.0);
  NumericVector GrowthRespiration(numCohorts,0.0);
  NumericVector LabileMassLeaf(numCohorts,0.0), LabileMassSapwood(numCohorts,0.0);
  NumericVector PlantSugarTransport(numCohorts,0.0), PlantSugarLeaf(numCohorts,0.0), PlantStarchLeaf(numCohorts,0.0);
  NumericVector PlantSugarSapwood(numCohorts,0.0), PlantStarchSapwood(numCohorts,0.0);
  NumericVector SapwoodArea(numCohorts,0.0), LeafArea(numCohorts,0.0), FineRootArea(numCohorts, 0.0);
  NumericVector SAgrowth(numCohorts,0.0), LAgrowth(numCohorts,0.0), FRAgrowth(numCohorts,0.0);
  NumericVector GrossPhotosynthesis(numCohorts,0.0);
  NumericVector PlantLAIdead(numCohorts,0.0), PlantLAIlive(numCohorts,0.0),PlantLAIexpanded(numCohorts,0.0);
  
  //Storage volume and maximum starch capacity for leaves and sapwood  
  NumericVector Volume_leaves(numCohorts,0.0);
  NumericVector Volume_sapwood(numCohorts,0.0);
  NumericVector Starch_max_leaves(numCohorts,0.0);
  NumericVector Starch_max_sapwood(numCohorts,0.0);
  NumericVector B_struct_leaves(numCohorts,0.0);
  NumericVector B_struct_sapwood(numCohorts,0.0);

  double minimumLeafSugarConc = equilibriumLeafTotalConc - nonSugarConcentration;
  double minimumSapwoodSugarConc = equilibriumSapwoodTotalConc - nonSugarConcentration;

  double rleafcellmax = relative_expansion_rate(0.0 ,25, -2.0,0.5,0.05,5.0);

  //3. Carbon balance, growth and senescence by cohort
  for(int j=0;j<numCohorts;j++){
    if(Status[j]=="alive") {
      double LAexpanded = leafArea(LAI_expanded[j], N[j]);
      double LAlive = leafArea(LAI_live[j], N[j]);
      double LAdead = leafArea(LAI_dead[j], N[j]);
      double LAlive_ini = LAlive;
      
      
      double costPerLA = 1000.0*leaf_CC/SLA[j]; // Construction cost in g gluc · m-2 of leaf area
      double costPerSA = sapwood_CC*(H[j]+(Z95[j]/10.0))*WoodDensity[j];  //Construction cost in g gluc ·cm-2 of sapwood
      double deltaLAgrowth = 0.0;
      double deltaSAgrowth = 0.0;
      NumericVector deltaFRBgrowth(numLayers, 0.0);
        
      Volume_leaves[j] = leafStorageVolume(LAI_expanded[j],  N[j], SLA[j], LeafDensity[j]);
      Volume_sapwood[j] = sapwoodStorageVolume(SA[j], H[j], L(j,_),V(j,_),WoodDensity[j], 0.5);
      Starch_max_leaves[j] = leafStarchCapacity(LAI_expanded[j],  N[j], SLA[j], 0.3)/Volume_leaves[j];
      Starch_max_sapwood[j] = sapwoodStarchCapacity(SA[j], H[j],L(j,_),V(j,_),WoodDensity[j], 0.2)/Volume_sapwood[j];
      B_struct_leaves[j] = leafStructuralBiomass(LAI_expanded[j],N[j],SLA[j]);
      B_struct_sapwood[j] = sapwoodStructuralLivingBiomass(SA[j], H[j], L(j,_),V(j,_), WoodDensity[j], 0.5);

      double labileMassLeafIni = (sugarLeaf[j]+starchLeaf[j])*(glucoseMolarMass*Volume_leaves[j]);
      double labileMassSapwoodIni = (sugarSapwood[j]+starchSapwood[j])*(glucoseMolarMass*Volume_sapwood[j]);
      
      double B_total = (B_struct_leaves[j] + B_struct_sapwood[j] + fineRootBiomass[j]+labileMassSapwoodIni+labileMassLeafIni);
      
      // Rcout << j << " Lvol: "<< Volume_leaves[j] << " Svol: "<<Volume_sapwood[j]<< " LStarchMax: "<<Starch_max_leaves[j]
      //       << " SStarchMax: "<<Starch_max_sapwood[j]<< " Bleaf "<< B_struct_leaves[j]<< " Bsap "<< B_struct_sapwood[j]<< " Bfr "<< B_struct_fineroots[j]<<"\n";
      
      double leafRespDay = 0.0;
      // double sfrRespDay = 0.0;
      
      //Estimate floem conductance as a factor of stem conductance
      double k_floem = VCstem_kmax[j]*floemConductanceFactor*(0.018/1000.0);
        
      //3.0 Xylogenesis
      grow_ring(ringList[j], psiSympStem[j] ,tday, 10.0);
      double rleafcell = relative_expansion_rate(psiSympLeaf[j] ,tday, LeafPI0[j],0.5,0.05,5.0);
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
        
        //Respiratory biomass (g dw · ind-1)
        double B_resp_leaves = B_struct_leaves[j] + leafSugarMassStep;
        double B_resp_sapwood = B_struct_sapwood[j] + sapwoodSugarMassStep;
        double B_resp_fineroots = fineRootBiomass[j];
        double QR = qResp(Tcan[s]);
        double leafRespStep = 0.0;
        if(LAlive>0.0) leafRespStep = B_resp_leaves*leaf_RR*QR/((double) numSteps);
        double sapwoodRespStep = B_resp_sapwood*sapwood_RR*QR/((double) numSteps);
        double finerootRespStep = B_resp_fineroots*fineroot_RR*QR/((double) numSteps);
        leafRespDay +=leafRespStep;
        
        double leafAgStepG = 0.0;
        if(LAlive>0.0) {
          //gross fotosynthesis
          double leafAgStepC = AgStep(j,s)/(N[j]/10000.0); //Translate g C · m-2 · h-1 to g C · h-1
          leafAgStepG = leafAgStepC*(glucoseMolarMass/(carbonMolarMass*6.0)); // from g C· h-1 to g gluc · h-1
          
          //Update output values
          GrossPhotosynthesisInst(j,s) = leafAgStepG/B_total; //Ag in g gluc · gdry-1
          GrossPhotosynthesis[j] += GrossPhotosynthesisInst(j,s); 
        }
        MaintenanceRespirationInst(j,s) = (leafRespStep+sapwoodRespStep+finerootRespStep)/B_total;//Rm in g gluc· gdry-1
        MaintenanceRespiration[j] += MaintenanceRespirationInst(j,s); 
        
        //Leaf growth
        double growthCostLAStep = 0.0;
        double growthCostSAStep = 0.0;
        double growthCostFRBStep = 0.0;        
        
        if(leafUnfolding[j]) {
          double deltaLApheno = std::max(leafAreaTarget[j] - LAlive, 0.0);
          double deltaLAsink = std::min(deltaLApheno, (1.0/((double) numSteps))*SA[j]*RGRleafmax*(rleafcell/rleafcellmax));
          if(LAlive>0.0) {
            double deltaLAavailable = std::max(0.0,((sugarLeaf[j] - minimumSugarGrowthLeaves)*(glucoseMolarMass*Volume_leaves[j]))/costPerLA);
            double deltaLAgrowthStep = std::min(deltaLAsink, deltaLAavailable);
            growthCostLAStep += deltaLAgrowthStep*costPerLA;
            deltaLAgrowth += deltaLAgrowthStep;
          } else { //Grow at expense of stem sugar
            double deltaLAavailable = std::max(0.0,((sugarSapwood[j] - minimumSugarGrowthLeaves)*(glucoseMolarMass*Volume_sapwood[j]))/costPerLA);
            double deltaLAgrowthStep = std::min(deltaLAsink, deltaLAavailable);
            // Rcout<<"hola"<< deltaLAavailable<< " "<< deltaLAsink<< " "<< deltaLAgrowthStep<<"\n";
            growthCostSAStep += deltaLAgrowthStep*costPerLA;
            deltaLAgrowth += deltaLAgrowthStep;
          }
        }
        
        //fine root growth
        if(fineRootBiomass[j] < fineRootBiomassTarget[j]) {
          for(int s = 0;s<numLayers;s++) {
            double deltaFRBsink = (1.0/((double) numSteps))*(V(j,s)*fineRootBiomass[j])*RGRfinerootmax*(rfineroot[s]/rleafcellmax);
            double deltaFRBavailable = std::max(0.0,((sugarSapwood[j] - minimumSugarGrowthFineRoots)*(glucoseMolarMass*Volume_sapwood[j]))/fineroot_CC);
            double deltaFRBgrowthStep = std::min(deltaFRBsink, deltaFRBavailable);
            growthCostFRBStep += deltaFRBgrowthStep*fineroot_CC;
            deltaFRBgrowth[s] += deltaFRBgrowthStep;
            // if(deltaFRBgrowthStep>0) Rcout << deltaFRBgrowthStep<<"\n";
          }
        }
        
        if(LAlive > 0.0) {
          List ring = ringList[j];
          NumericVector SAring = ring["SA"];
          double deltaSAring = 0.0;
          if(SAring.size()==1) deltaSAring = SAring[0];
          else deltaSAring = SAring[SAring.size()-1] - SAring[SAring.size()-2];
          double RGRcellmax = (2e-8/SA[j]);
          double deltaSAsink = (1e-8*(deltaSAring/10.0))*(RGRsapwoodmax[j]/RGRcellmax)/((double) numSteps); //Correction for the difference in the number of cells
          double deltaSAavailable = std::max(0.0,((starchSapwood[j]- minimumStarchGrowthSapwood)*(glucoseMolarMass*Volume_sapwood[j]))/costPerSA);
          double deltaSAgrowthStep = std::min(deltaSAsink, deltaSAavailable);
          growthCostSAStep += deltaSAgrowthStep*costPerSA; //increase cost (may be non zero if leaf growth was charged onto sapwood)
          deltaSAgrowth  +=deltaSAgrowthStep;
          // Rcout<< SAring.size()<<" " <<j<< " "<< psiSympStem[j]<< " "<< StemPI0[j]<<" dSAring "<<deltaSAring<< " dSAsink "<< deltaSAsink<<" dSAgrowth "<< deltaSAgrowthStep<<"\n";
        }      
        
        GrowthRespirationInst(j,s) += (growthCostLAStep + growthCostSAStep + growthCostFRBStep)/B_total;
        GrowthRespiration[j] +=GrowthRespirationInst(j,s); //growth cost in g gluc · gdry-1
        //Instantaneous carbon balance
        CarbonBalanceInst(j,s) = GrossPhotosynthesisInst(j,s) - MaintenanceRespirationInst(j,s) - GrowthRespirationInst(j,s);
        CarbonBalance[j] +=CarbonBalanceInst(j,s);
        
        //sugar mass balance
        double leafSugarMassDeltaStep = leafAgStepG - leafRespStep - growthCostLAStep;
        double sapwoodSugarMassDeltaStep = - finerootRespStep  - growthCostFRBStep;
        double sapwoodStarchMassDeltaStep = - growthCostSAStep - sapwoodRespStep;
        // Rcout<<" coh:"<<j<< " s:"<<s<<" dS: "<< leafSugarMassDeltaStep<<" sugar mass leaf: "<< leafSugarMassStep << " dS:"<< sapwoodSugarMassDeltaStep<< " sugar mass sap: "<< sapwoodSugarMassStep<<"\n";
        
        
        //floem transport      
        
        double ff = 0.0;
        double ctl = 3600.0*Volume_leaves[j]*glucoseMolarMass;
        double cts = 3600.0*Volume_sapwood[j]*glucoseMolarMass;
        for(int t=0;t<3600;t++) {
          sugarSapwood[j] += sapwoodSugarMassDeltaStep/cts;
          starchSapwood[j] += sapwoodStarchMassDeltaStep/cts;
          
          double conversionSapwood = sugarStarchDynamicsStem(sugarSapwood[j]/StemSympRWCInst(j,s), starchSapwood[j]/StemSympRWCInst(j,s), minimumSapwoodSugarConc);
          // Rcout<<" coh:"<<j<< " s:"<<s<< " Lsugar: "<< sugarLeaf[j] << " Lstarch: "<< sugarSapwood[j]<<" starch formation: "<<conversionLeaf<< "\n";
          double starchSapwoodIncrease = conversionSapwood*StemSympRWCInst(j,s);
          starchSapwoodIncrease = std::min(starchSapwoodIncrease, Starch_max_sapwood[j] - starchSapwood[j]);
          starchSapwood[j] += starchSapwoodIncrease;
          
          if(LAlive>0.0) {
            sugarLeaf[j] += leafSugarMassDeltaStep/ctl;
            double ft = floemFlow(LeafSympPsiInst(j,s), StemSympPsiInst(j,s), sugarLeaf[j]/LeafSympRWCInst(j,s), sugarSapwood[j]/StemSympRWCInst(j,s), Tcan[s], k_floem, nonSugarConcentration)*LAlive; //flow as mol glucose per s
            // sugar-starch dynamics
            double conversionLeaf = sugarStarchDynamicsLeaf(sugarLeaf[j]/LeafSympRWCInst(j,s), starchLeaf[j]/LeafSympRWCInst(j,s), minimumLeafSugarConc);
            double starchLeafIncrease = conversionLeaf*LeafSympRWCInst(j,s);
            starchLeafIncrease = std::min(starchLeafIncrease, Starch_max_leaves[j] - starchLeaf[j]);
            starchLeaf[j]  += starchLeafIncrease;
            // Rcout<<" coh:"<<j<< " s:"<<s<< " Ssugar: "<< sugarSapwood[j] << " Sstarch: "<< starchSapwood[j]<<" starch formation: "<<conversionSapwood<< "\n";
            //Apply floem transport (mol gluc) to sugar concentrations (mol gluc· l-1)
            sugarLeaf[j]  +=  (-ft/Volume_leaves[j]) - starchLeafIncrease;
            sugarSapwood[j] +=  (ft/Volume_sapwood[j]) - starchSapwoodIncrease;
            ff +=ft;
          } else {
            sugarSapwood[j] += - starchSapwoodIncrease;
          }
        }
        PlantSugarLeafInst(j,s) = sugarLeaf[j];
        PlantSugarSapwoodInst(j,s) = sugarSapwood[j];
        PlantStarchLeafInst(j,s) = starchLeaf[j];
        PlantStarchSapwoodInst(j,s) = starchSapwood[j];
        PlantSugarTransportInst(j,s) = 1000.0*ff/(3600.0); //mmol·s-1
        PlantSugarTransport[j] += ff; //To calculate daily floem balance (positive means towards stem)
        // Rcout<<" coh:"<<j<< " s:"<<s<< " conc leaf: "<< sugarLeaf[j] << " conc sap: "<< sugarSapwood[j]<<" ff: "<<ff<< "\n";
        
        // Rcout<<j<<" LeafTLP "<< turgorLossPoint(LeafPI0[j], LeafEPS[j])<< " Leaf PI "<< osmoticWaterPotential(sugarLeaf[j], tday)<< " Conc "<< sugarLeaf[j]<< " TLPconc"<< tlpConcLeaf<<"\n";
      }
  
      //3.2 Apply growth and senescence
      //3.2.1 Leaf growth and senescence
      //Leaf senescense due to C starvation
      if(sugarLeaf[j] < 0.0) { 
        double respirationExcess = -sugarLeaf[j]*(Volume_leaves[j]*glucoseMolarMass); //g gluc
        double propExcess = respirationExcess/leafRespDay; //day
        // Rcout<< j <<" Excess respiration: " << respirationExcess << " Prop:"<< propExcess<< " LAlive " << LAlive << " LAlivenew "<< LAlive*(1.0 - propExcess) <<"\n";
        LAdead = LAdead + LAexpanded*propExcess;
        LAexpanded = LAexpanded*(1.0 - propExcess);
        LAlive = LAlive*(1.0 - propExcess);
        if(LAlive<0.0001) {
          LAlive = 0.0;
          LAexpanded = 0.0;
        }
        sugarLeaf[j] = 0.0;
      }
      //Leaf growth
      LAlive += deltaLAgrowth; //Update leaf area
      LAexpanded +=deltaLAgrowth;
      //Leaf senescence
      double propLeafSenescence = 0.0;
      //Leaf senescence due to age (Ca+ accumulation) only in evergreen species
      if(phenoType[j] == "oneflush-evergreen" || phenoType[j] == "progressive-evergreen") {
        propLeafSenescence = (1.0/(365.25*leafDuration[j]));
      }
      //Leaf senescence due to drought 
      double LAplc = std::min(LAlive, (1.0 - StemPLC[j])*leafAreaTarget[j]);
      if(LAplc<LAlive) {
        // Rcout<<j<< " "<< LAplc<< " "<< LAlive<<"\n";
        propLeafSenescence = std::max((LAlive-LAplc)/LAlive, propLeafSenescence); 
      }
      //Complete defoliation if RWCsymp < 0.5
      if(LAlive > 0.0 && LeafSympRWC[j]<0.5){
        propLeafSenescence = 1.0;
        Rcout<<" [Cohort "<< j<<" defoliated ] ";
      }
      //Apply senescence due to age/drought
      double LA_exp_prev= LAexpanded; //Store previous value
      LAdead += LAexpanded*propLeafSenescence;
      LAexpanded = LAexpanded*(1.0 - propLeafSenescence); //Update expanded leaf area
      LAlive = LAlive*(1.0 - propLeafSenescence); //Update expanded leaf area
      
      //3.2.2 FRB growth and senescence
      NumericVector newFRB(numLayers,0.0);
      for(int s=0;s<numLayers;s++) {
        double initialFRB = fineRootBiomass[j]*V(j,s);
        double dayTurnover = dailyFineRootTurnoverProportion*std::max(0.0,(Tsoil[s]-5.0)/20.0);
        newFRB[s] = std::max(0.0,initialFRB*(1.0 - dayTurnover) + deltaFRBgrowth[s]);
      }
      fineRootBiomass[j] = sum(newFRB);
      //Update vertical fine root distribution
      for(int s=0;s<numLayers;s++) { 
        V(j,s) = newFRB[s]/fineRootBiomass[j];
      }

      //3.2.3 SA growth and senescense
      double SAprev = SA[j];
      double deltaSAturnover = (dailySapwoodTurnoverProportion/(1.0+15.0*exp(-0.01*H[j])))*SA[j];
      SA[j] = SA[j] - deltaSAturnover; //Update sapwood area
      //SA growth     
      SA[j] += deltaSAgrowth; //Update sapwood area
      //Decrease PLC due to new SA growth
      if(cavitationRefill=="growth") StemPLC[j] = std::max(0.0, StemPLC[j] - (deltaSAgrowth/SA[j]));
      
      //3.3. Determine plant death by carbon starvation or dessication
      if((sugarSapwood[j]<0.0) || (StemSympRWC[j] <0.5)) {
        LAdead = LAlive;
        LAlive = 0.0;
        LAexpanded = 0.0;
        if(sugarSapwood[j]<0.0) Status(j) = "starvation";
        else if(StemSympRWC[j]<0.5) Status(j) = "dessication";
        Rcout<<" [Cohort "<< j<<" died from " << Status(j)<<"] ";
      }
      
      
      //3.4 Update functional variables (feed-back to soil-plant-water balance)
      LAI_live[j] = LAlive*N[j]/10000.0;
      LAI_expanded[j] = LAexpanded*N[j]/10000.0;
      LAI_dead[j] = LAdead*N[j]/10000.0;
      
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
        // StemPI0[j] = osmoticWaterPotential(sugarSapwood[j], 20.0, nonSugarConc);
        //Update non-stomatal photosynthesis limitations
        if(nonStomatalPhotosynthesisLimitation) NSPL[j] = 1.0 - std::max(0.0, std::min(1.0, sugarLeaf[j] - 0.5)); //photosynthesis limited when conc > 0.5 and zero when conc > 1.5 mol·l-1
        else NSPL[j] = 1.0;

        //3.5 Update allocation (leaf area and fine root biomass) targets
        //Set leaf area target if bud formation is allowed
        if(budFormation[j]) {
          if(allocationStrategy == "Plant_kmax") {
            leafAreaTarget[j] = LAlive*(Plant_kmax[j]/allocationTarget[j]);
          } else if(allocationStrategy =="Al2As") {
            leafAreaTarget[j] = (SA[j]/10000.0)*allocationTarget[j];
          }
        }
        //Update fine root biomass target      
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
      

      //3.6 Output variables by cohort
      PlantSugarLeaf[j] = sugarLeaf[j];
      PlantStarchLeaf[j] = starchLeaf[j];
      PlantSugarSapwood[j] = sugarSapwood[j];
      PlantStarchSapwood[j] = starchSapwood[j];
      SapwoodArea[j] = SA[j];
      LeafArea[j] = LAexpanded;
      FineRootArea[j] = fineRootBiomass[j]*specificRootSurfaceArea(SRL[j], FineRootDensity[j])*1e-4;
      LAgrowth[j] += deltaLAgrowth/SA[j];//Store Leaf area growth rate in relation to sapwood area (m2·cm-2·d-1)
      SAgrowth[j] += deltaSAgrowth/SA[j]; //Store sapwood area growth rate (cm2·cm-2·d-1)
      FRAgrowth[j] = sum(deltaFRBgrowth)*specificRootSurfaceArea(SRL[j], FineRootDensity[j])*1e-4/SA[j];//Store fine root area growth rate (m2·cm-2·d-1)
      PlantLAIlive[j] = LAI_live[j];
      PlantLAIexpanded[j] = LAI_expanded[j];
      PlantLAIdead[j] = LAI_dead[j];
      LabileMassLeaf[j] = (sugarLeaf[j]+starchLeaf[j])*(glucoseMolarMass*Volume_leaves[j]);
      LabileMassSapwood[j] = (sugarSapwood[j]+starchSapwood[j])*(glucoseMolarMass*Volume_sapwood[j]);
      
      //Carbon balance check
      // double sugarTransportMass = PlantSugarTransport[j]*glucoseMolarMass;
      // double sumLeaf = GrossPhotosynthesis[j] - LeafMaintenanceRespiration[j] - LeafGrowthRespiration[j] - sugarTransportMass;
      // double sumSapwood = sugarTransportMass - SapwoodMaintenanceRespiration[j] - FineRootMaintenanceRespiration[j] - SapwoodGrowthRespiration[j] - FineRootGrowthRespiration[j];
      // Rcout<<j<<" CBLeaf "<< sumLeaf << " ChLabLeaf: "<< (LabileMassLeaf[j] - labileMassLeafIni);
      // Rcout<<" CBSapwood "<< sumSapwood << " ChLabSapwood: "<< (LabileMassSapwood[j] - labileMassSapwoodIni)<<"\n";
    }
  }
  //Update pool proportions??
  if(plantWaterPools) {
    NumericVector poolProportions = Rcpp::as<Rcpp::NumericVector>(belowdf["poolProportions"]);
    // for(int j=0;j<numCohorts;j++) poolProportions[j] = LAI_live[j]/sum(LAI_live);
    //Update RHOP
    List RHOP = belowLayers["RHOP"];
    List newRHOP = horizontalProportions(poolProportions, CRSV, N, V,dVec, rfc);
    for(int j=0;j<numCohorts;j++) RHOP[j] = newRHOP[j];
  }
  
  //Needed with string vectors
  above["Status"] = Status;
  
  GrossPhotosynthesisInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,numSteps));
  MaintenanceRespirationInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,numSteps));
  CarbonBalanceInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,numSteps));
  GrowthRespirationInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,numSteps));
  PlantSugarLeafInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,numSteps));
  PlantStarchLeafInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,numSteps));
  PlantSugarSapwoodInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,numSteps));
  PlantStarchSapwoodInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,numSteps));
  PlantSugarTransportInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,numSteps));
  List plantCBInst = List::create(
    _["GrossPhotosynthesis"] = GrossPhotosynthesisInst,
    _["MaintenanceRespiration"] = MaintenanceRespirationInst,
    _["GrowthRespiration"] = GrowthRespirationInst,
    _["CarbonBalance"] = CarbonBalanceInst,
    _["SugarLeaf"] = PlantSugarLeafInst,
    _["StarchLeaf"] = PlantStarchLeafInst,
    _["SugarSapwood"] = PlantSugarSapwoodInst,
    _["StarchSapwood"] = PlantStarchSapwoodInst,
    _["SugarTransport"] = PlantSugarTransportInst
  );
  
  DataFrame plantCarbonBalance = DataFrame::create(_["GrossPhotosynthesis"] = GrossPhotosynthesis,
                                   _["MaintenanceRespiration"] = MaintenanceRespiration,
                                   _["GrowthRespiration"] = GrowthRespiration,
                                   _["CarbonBalance"] = CarbonBalance,
                                   _["SugarLeaf"] = PlantSugarLeaf,
                                   _["StarchLeaf"] = PlantStarchLeaf,
                                   _["SugarSapwood"] = PlantSugarSapwood,
                                   _["StarchSapwood"] = PlantStarchSapwood,
                                   _["SugarTransport"] = PlantSugarTransport,
                                   _["LabileMassLeaf"] = LabileMassLeaf,
                                   _["LabileMassSapwood"] = LabileMassSapwood,
                                   _["StemPI0"] = clone(StemPI0), //Store a copy of the current osmotic potential at full turgor
                                   _["LeafPI0"] = clone(LeafPI0));
  plantCarbonBalance.attr("row.names") = above.attr("row.names");
  
  DataFrame plantGrowth = List::create(
    _["SapwoodArea"] = SapwoodArea,
    _["LeafArea"] = LeafArea,
    _["FineRootArea"] = FineRootArea,
    _["SAgrowth"] = SAgrowth,
    _["LAgrowth"] = LAgrowth,
    _["FRAgrowth"] = FRAgrowth
    // _["LAIlive"] = PlantLAIlive,
    // _["LAIexpanded"] = PlantLAIexpanded,
    // _["LAIdead"] = PlantLAIdead,
  );
  plantGrowth.attr("row.names") = above.attr("row.names");
  
  List l = List::create(_["cohorts"] = clone(cohorts),
                        _["WaterBalance"] = spwbOut["WaterBalance"], 
                        _["EnergyBalance"] = spwbOut["EnergyBalance"],
                        _["Soil"] = spwbOut["Soil"], 
                        _["Stand"] = spwbOut["Stand"], 
                        _["Plants"] = spwbOut["Plants"],
                        _["PlantCarbonBalance"] = plantCarbonBalance,
                        _["PlantGrowth"] = plantGrowth,
                        _["RhizoPsi"] = spwbOut["RhizoPsi"],
                        _["SunlitLeaves"] = spwbOut["SunlitLeaves"],
                        _["ShadeLeaves"] = spwbOut["ShadeLeaves"],
                        _["ExtractionInst"] = spwbOut["ExtractionInst"],
                        _["PlantsInst"] = spwbOut["PlantsInst"],
                        _["SunlitLeavesInst"] = spwbOut["SunlitLeavesInst"],
                        _["ShadeLeavesInst"] = spwbOut["ShadeLeavesInst"],
                        _["PlantCBInst"] = plantCBInst,
                        _["LightExtinction"] = spwbOut["LightExtinction"],
                        _["WindExtinction"] = spwbOut["WindExtinction"]);
  l.attr("class") = CharacterVector::create("growth_day","list");
  return(l);
}


// [[Rcpp::export("growth_day")]]
List growthDay(List x, List soil, CharacterVector date, double tmin, double tmax, double rhmin, 
               double rhmax, double rad, double wind, 
               double latitude, double elevation, double slope, double aspect,  
               double prec, double runon=0.0) {
  //Control parameters
  List control = x["control"];
  bool verbose = control["verbose"];
  bool leafPhenology = control["leafPhenology"];
  String transpirationMode = control["transpirationMode"];
  std::string c = as<std::string>(date[0]);
  int J = meteoland::radiation_julianDay(std::atoi(c.substr(0, 4).c_str()),std::atoi(c.substr(5,2).c_str()),std::atoi(c.substr(8,2).c_str()));
  double delta = meteoland::radiation_solarDeclination(J);
  double solarConstant = meteoland::radiation_solarConstant(J);
  double tday = meteoland::utils_averageDaylightTemperature(tmin, tmax);
  double latrad = latitude * (PI/180.0);
  double asprad = aspect * (PI/180.0);
  double slorad = slope * (PI/180.0);
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
    s = growthDay1(x,soil, tday, pet, prec, er, runon, rad, elevation, verbose);
  } else {
    s = growthDay2(x,soil, tmin, tmax, tmin, tmax, tmin, rhmin, rhmax, rad, wind, 
                 latitude, elevation, slope, aspect,
                 solarConstant, delta, prec, pet, er, runon, verbose);
  }
  // Rcout<<"hola4\n";
  return(s);
}


void checkgrowthInput(List x, List soil, String transpirationMode, String soilFunctions) {
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
  if(!paramsGrowth.containsElementNamed("RGRsapwoodmax")) stop("RGRsapwoodmax missing in growthInput$paramsGrowth");
  
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
    if(!paramsTranspiration.containsElementNamed("Gwmax")) stop("Gwmax missing in growthInput$paramsTransp");
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

void recordStandSummary(DataFrame standSummary, DataFrame above, int pos) {
  
  NumericVector DBH = above["DBH"];
  NumericVector Cover = above["Cover"];
  NumericVector H = above["H"];
  NumericVector N = above["N"];
  NumericVector LAI_live = above["LAI_live"];
  StringVector Status = above["Status"];
  
  NumericVector SLAI = as<Rcpp::NumericVector>(standSummary["LeafAreaIndex"]);
  SLAI[pos] = 0.0;
  NumericVector TBAL = as<Rcpp::NumericVector>(standSummary["TreeBasalAreaLive"]);
  TBAL[pos] = 0.0;
  NumericVector TBAD = as<Rcpp::NumericVector>(standSummary["TreeBasalAreaDead"]);
  TBAD[pos] = 0.0;
  NumericVector SCoverL = as<Rcpp::NumericVector>(standSummary["ShrubCoverLive"]);
  SCoverL[pos] = 0.0;
  NumericVector SCoverD = as<Rcpp::NumericVector>(standSummary["ShrubCoverDead"]);
  SCoverD[pos] = 0.0;
  NumericVector MaxHeight = as<Rcpp::NumericVector>(standSummary["MaxHeight"]);
  MaxHeight[pos] = 0.0;
  int numCohorts = N.length();
  NumericVector treeBA = treeBasalArea(N, DBH);
  for(int i=0;i<numCohorts;i++) {
    SLAI[pos] += LAI_live[i];
    if(!NumericVector::is_na(treeBA[i])) {
      if(Status[i]=="alive") TBAL[pos] += treeBA[i];
      else TBAD[pos] += treeBA[i];
      MaxHeight[pos] = std::max(MaxHeight[pos], H[i]);
    } else {
      if(Status[i]=="alive") SCoverL[pos] +=Cover[i];
      else SCoverD[pos] +=Cover[i];
    }
  }
}


// [[Rcpp::export("growth")]]
List growth(List x, List soil, DataFrame meteo, double latitude, double elevation = NA_REAL, double slope = NA_REAL, double aspect = NA_REAL) {
  
  //Control params 
  List control =x["control"];  
  //Store input
  List growthInput = clone(x);
  List soilInput = clone(soil);
  
    
  // Rcout<<"1";
  
  //Cohort info
  DataFrame cohorts = Rcpp::as<Rcpp::DataFrame>(x["cohorts"]);
  NumericVector SP = cohorts["SP"];

  String transpirationMode = control["transpirationMode"];
  String soilFunctions = control["soilFunctions"];
  
  bool verbose = control["verbose"];
  bool snowpack = control["snowpack"];
  bool subdailyResults = control["subdailyResults"];
  bool leafPhenology = control["leafPhenology"];
  bool unlimitedSoilWater = control["unlimitedSoilWater"];
  checkgrowthInput(x, soil, transpirationMode, soilFunctions);
  
  if(NumericVector::is_na(latitude)) stop("Value for 'latitude' should not be missing.");
  double latrad = latitude * (PI/180.0);
  
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
  }
  CharacterVector dateStrings = meteo.attr("row.names");
  
  IntegerVector DOY = date2doy(dateStrings);
  NumericVector Photoperiod = date2photoperiod(dateStrings, latrad);
  
  //Canpopy parameters
  List canopyParams = x["canopy"];
  
  //Aboveground parameters  
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  int numCohorts = SP.size();

  //Belowground state variables  
  List below = Rcpp::as<Rcpp::List>(x["below"]);
  NumericVector Z95 = Rcpp::as<Rcpp::NumericVector>(below["Z95"]);

  //Internal state variables
  DataFrame internalWater = Rcpp::as<Rcpp::List>(x["internalWater"]);
  DataFrame internalCarbon = Rcpp::as<Rcpp::List>(x["internalCarbon"]);
  

  //Base parameters
  DataFrame paramsPhenology = Rcpp::as<Rcpp::DataFrame>(x["paramsPhenology"]);
  NumericVector Sgdd = Rcpp::as<Rcpp::NumericVector>(paramsPhenology["Sgdd"]);
  
  DataFrame paramsInterception = Rcpp::as<Rcpp::DataFrame>(x["paramsInterception"]);
  NumericVector kPAR = paramsInterception["kPAR"];
  
  // Rcout<<"3";


  //Soil input
  NumericVector Water_FC = waterFC(soil, soilFunctions);
  int nlayers = Water_FC.size();
  
  //Anatomy parameters
  DataFrame paramsAnatomy = Rcpp::as<Rcpp::DataFrame>(x["paramsAnatomy"]);
  NumericVector Hmax  = paramsAnatomy["Hmax"];
  NumericVector SLA = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["SLA"]);
  NumericVector r635  = paramsAnatomy["r635"];
  
  // NumericVector Al2As = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["Al2As"]);
  // NumericVector WoodDensity = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["WoodDensity"]);

  
  //Allometric parameters
  DataFrame paramsAllometries = Rcpp::as<Rcpp::DataFrame>(x["paramsAllometries"]);
  NumericVector Aash  = paramsAllometries["Aash"];
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
  
  
  
  //Detailed subday results
  List subdailyRes(numDays);
  
  //EnergyBalance output variables
  DataFrame DEB = defineEnergyBalanceDailyOutput(meteo);
  DataFrame DT = defineTemperatureDailyOutput(meteo);
  
  //Plant carbon output variables
  NumericMatrix CarbonBalance(numDays, numCohorts);
  NumericMatrix MaintenanceRespiration(numDays, numCohorts);
  NumericMatrix GrowthRespiration(numDays, numCohorts);
  // NumericMatrix LabileMassLeaf(numDays, numCohorts);
  // NumericMatrix LabileMassSapwood(numDays, numCohorts);
  NumericMatrix PlantSugarLeaf(numDays, numCohorts);
  NumericMatrix PlantStarchLeaf(numDays, numCohorts);
  NumericMatrix PlantSugarSapwood(numDays, numCohorts);
  NumericMatrix PlantStarchSapwood(numDays, numCohorts);
  NumericMatrix PlantSugarTransport(numDays, numCohorts);
  NumericMatrix SapwoodArea(numDays, numCohorts);
  NumericMatrix LeafArea(numDays, numCohorts);
  NumericMatrix FineRootArea(numDays, numCohorts);
  NumericMatrix SAgrowth(numDays, numCohorts);
  NumericMatrix LAgrowth(numDays, numCohorts);
  NumericMatrix FRAgrowth(numDays, numCohorts);
  NumericMatrix GrossPhotosynthesis(numDays, numCohorts);
  NumericMatrix PlantLAIexpanded(numDays, numCohorts), PlantLAIdead(numDays, numCohorts), PlantLAIlive(numDays, numCohorts);
  NumericVector SAgrowthcum(numCohorts, 0.0);
  NumericMatrix StemPI0(numDays, numCohorts), LeafPI0(numDays, numCohorts);
  
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
  DataFrame standSummary = DataFrame::create(
    _["LeafAreaIndex"] = NumericVector(numYears+1,0.0),
    _["TreeBasalAreaLive"] = NumericVector(numYears+1,0.0),
    _["TreeBasalAreaDead"] = NumericVector(numYears+1,0.0),
    _["ShrubCoverLive"] = NumericVector(numYears+1,0.0),
    _["ShrubCoverDead"] = NumericVector(numYears+1,0.0),
    _["MaxHeight"] = NumericVector(numYears+1,0.0)
  );
  List standStructures(numYears+1);
  CharacterVector nss(numYears+1);
  for(int i=0;i<(numYears+1);i++) {
    if(i==0) {
      nss[i] = "Initial"; 
    } else {
      char Result[16];
      sprintf(Result, "Year_%d", i);
      nss[i] = Result;
    }
  }
  standSummary.attr("row.names") = nss;
  standStructures.attr("names") = nss;
  standStructures[0] = clone(above);
  recordStandSummary(standSummary, above, 0);
  
  NumericVector initialContent = water(soil, soilFunctions);
  double initialSnowContent = soil["SWE"];
  if(verbose) {
    Rcout<<"Initial soil water content (mm): "<< sum(initialContent)<<"\n";
    Rcout<<"Initial snowpack content (mm): "<< initialSnowContent<<"\n";
  }
  
  if(verbose) Rcout << "Performing daily simulations ";
  List s;
  int iyear = 0;
  for(int i=0;i<numDays;i++) {
    if(verbose && (i%10 == 0)) Rcout<<".";//<<i;
    
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
      double er = erFactor(DOY[i], PET[i], Precipitation[i]);
      s = growthDay1(x, soil, MeanTemperature[i], PET[i], Precipitation[i], er, 0.0, 
                     Radiation[i], elevation, false); //No Runon in simulations for a single cell
    } else if(transpirationMode=="Sperry") {
      int ntimesteps = control["ndailysteps"];
      double tstep = 86400.0/((double) ntimesteps);
      std::string c = as<std::string>(dateStrings[i]);
      int J = meteoland::radiation_julianDay(std::atoi(c.substr(0, 4).c_str()),std::atoi(c.substr(5,2).c_str()),std::atoi(c.substr(8,2).c_str()));
      double delta = meteoland::radiation_solarDeclination(J);
      double solarConstant = meteoland::radiation_solarConstant(J);
      double latrad = latitude * (PI/180.0);
      if(NumericVector::is_na(aspect)) aspect = 0.0;
      if(NumericVector::is_na(slope)) slope = 0.0;
      double asprad = aspect * (PI/180.0);
      double slorad = slope * (PI/180.0);
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
      PET[i] = meteoland::penman(latrad, elevation, slorad, asprad, J, tmin, tmax, rhmin, rhmax, rad, wind);
      double er = erFactor(DOY[i], PET[i], Precipitation[i]);
      s = growthDay2(x, soil, tmin, tmax, tminPrev, tmaxPrev, tminNext,
                   rhmin, rhmax, rad, wind, 
                   latitude, elevation, slope, aspect,
                   solarConstant, delta, Precipitation[i], PET[i], 
                   er, 0.0, verbose);
      // Rcout<<" coh 1: "<< Status[1]<<"\n";
      
      fillEnergyBalanceTemperatureDailyOutput(DEB,DT,s,i);
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
    DataFrame cb = Rcpp::as<Rcpp::DataFrame>(s["PlantCarbonBalance"]);
    DataFrame pg = Rcpp::as<Rcpp::DataFrame>(s["PlantGrowth"]);
    
    
    //4. Assemble output
    CarbonBalance(i,_) = Rcpp::as<Rcpp::NumericVector>(cb["CarbonBalance"]);
    MaintenanceRespiration(i,_) = Rcpp::as<Rcpp::NumericVector>(cb["MaintenanceRespiration"]);
    GrowthRespiration(i,_) = Rcpp::as<Rcpp::NumericVector>(cb["GrowthRespiration"]);
    GrossPhotosynthesis(i,_) = Rcpp::as<Rcpp::NumericVector>(cb["GrossPhotosynthesis"]);
    // LabileMassLeaf(i,_) = Rcpp::as<Rcpp::NumericVector>(cb["LabileMassLeaf"]);
    // LabileMassSapwood(i,_) = Rcpp::as<Rcpp::NumericVector>(cb["LabileMassSapwood"]);
    PlantSugarLeaf(i,_) = Rcpp::as<Rcpp::NumericVector>(cb["SugarLeaf"]);
    PlantStarchLeaf(i,_) = Rcpp::as<Rcpp::NumericVector>(cb["StarchLeaf"]);
    PlantSugarSapwood(i,_) = Rcpp::as<Rcpp::NumericVector>(cb["SugarSapwood"]);
    PlantStarchSapwood(i,_) = Rcpp::as<Rcpp::NumericVector>(cb["StarchSapwood"]);
    PlantSugarTransport(i,_) = Rcpp::as<Rcpp::NumericVector>(cb["SugarTransport"]);
    StemPI0(i,_) = Rcpp::as<Rcpp::NumericVector>(cb["StemPI0"]); 
    LeafPI0(i,_) = Rcpp::as<Rcpp::NumericVector>(cb["LeafPI0"]); 
    
    SapwoodArea(i,_) = Rcpp::as<Rcpp::NumericVector>(pg["SapwoodArea"]);
    LeafArea(i,_) = Rcpp::as<Rcpp::NumericVector>(pg["LeafArea"]);
    LAgrowth(i,_) = Rcpp::as<Rcpp::NumericVector>(pg["LAgrowth"]);
    SAgrowth(i,_) = Rcpp::as<Rcpp::NumericVector>(pg["SAgrowth"]);
    if(transpirationMode=="Sperry") {
      FRAgrowth(i,_) = Rcpp::as<Rcpp::NumericVector>(pg["FRAgrowth"]);
      FineRootArea(i,_) = Rcpp::as<Rcpp::NumericVector>(pg["FineRootArea"]);
    }
      
    for(int j=0;j<numCohorts;j++){
      SAgrowthcum[j] += SAgrowth(i,j)*SapwoodArea(i,j); //Store cumulative SA growth (for structural variable update)
    }
    
    //4 Update structural variables
    if(((DOY[i]==1) && (i>0)) | ((i==(numDays-1)) && (DOY[i]>=365))) { 
      if(verbose) Rcout<<" [update structural variables] ";
      iyear++;
      
      NumericVector DBH = above["DBH"];
      NumericVector Cover = above["Cover"];
      NumericVector H = above["H"];
      NumericVector N = above["N"];
      NumericVector CR = above["CR"];
      NumericVector LAI_live = above["LAI_live"];
      NumericVector LAI_expanded = above["LAI_expanded"];
      NumericVector LAI_dead = above["LAI_dead"];
      NumericVector SA = above["SA"];
      StringVector Status = above["Status"];
      
      DataFrame internalAllocation = Rcpp::as<Rcpp::DataFrame>(x["internalAllocation"]);
      NumericVector allocationTarget = internalAllocation["allocationTarget"];
      NumericVector leafAreaTarget = internalAllocation["leafAreaTarget"];
      
      NumericVector deltaDBH(numCohorts, 0.0);
      for(int j=0;j<numCohorts; j++) {
        if(!NumericVector::is_na(DBH[j])) {
          deltaDBH[j] = 2.0*sqrt(pow(DBH[j]/2.0,2.0)+(SAgrowthcum[j]/PI)) - DBH[j];
          DBH[j] = DBH[j] + deltaDBH[j];
        } 
        SAgrowthcum[j] = 0.0; //Reset cumulative growth
      }

      NumericVector L = parcohortC(H, LAI_live, LAI_dead, kPAR, CR);
      for(int j=0;j<numCohorts; j++) {
        if(!NumericVector::is_na(DBH[j]) && Status[j]=="alive") {
          double fHmod = std::max(0.0,std::min(1.0,(1.0-((H[j]-137.0)/(Hmax[j]-137.0)))));
          double fHD = (fHDmin[j]*(L[j]/100.0) + fHDmax[j]*(1.0-(L[j]/100.0)))*fHmod;
          // Rcout << fHmod<<" "<< fHD<<" "<< L[j]<<"\n";
          H[j] = H[j] + fHD*deltaDBH[j];
        }
      }
      NumericVector crNew = treeCrownRatioMED(N, DBH, H, Acw, Bcw, Acr, B1cr, B2cr, B3cr, C1cr, C2cr);
      for(int j=0;j<numCohorts; j++) {
        if(!NumericVector::is_na(DBH[j]) && Status[j]=="alive") {
          CR[j] = crNew[j];
        }
      }

      //Shrub variables
      for(int j=0;j<numCohorts; j++) {
        if(NumericVector::is_na(DBH[j]) && Status[j]=="alive") {
          double Wleaves = leafAreaTarget[j]/SLA[j];  //Calculates the biomass (kg dry weight) of leaves
          double PV = pow(Wleaves*r635[j]/Absh[j], 1.0/Bbsh[j]); //Calculates crown phytovolume (in m3)
          H[j] = pow(1e6*PV/(Aash[j]*CR[j]), 1.0/3.0); //Updates shrub height
          if(H[j]> Hmax[j]) { //Limit height (and update the former variables)
            H[j] = Hmax[j];
            // PV = (Aash[j]*pow(H[j],2.0)/10000.0)*(H[j]/100.0)*CR[j];
            // Wleaves = Absh[j]*pow(PV, Bbsh[j])/r635[j];
            // double prevLive = LAI_live[j];
            // LAI_live[j] = Wleaves*((N[j]/10000)*SLA[j]); //Update LAI_live to the maximum
            // LAI_dead[j] += prevLive - LAI_live[j]; //Increment dead LAI with the difference
          }
          Cover[j] = (N[j]*Aash[j]*pow(H[j],2.0)/1e6); //Updates shrub cover
        }
      }
      
      //reset ring structures
      List ringList = x["internalRings"];
      for(int j=0;j<numCohorts; j++) ringList[j] = initialize_ring();
      // Store stand structure
      standStructures[iyear] = clone(above);
      recordStandSummary(standSummary, above,iyear);
    }

    if(subdailyResults) {
      subdailyRes[i] = clone(s);
    }
  }
  if(verbose) Rcout << "done.\n";
  
  if(verbose) {
    printWaterBalanceResult(DWB, plantDWOL, soil, soilFunctions,
                            initialContent, initialSnowContent,
                            transpirationMode);
  }
  
  //Add matrix dimnames
  CarbonBalance.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names"));
  GrossPhotosynthesis.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names"));
  MaintenanceRespiration.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names"));
  GrowthRespiration.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names"));
  // LabileMassLeaf.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names"));
  // LabileMassSapwood.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names"));
  PlantSugarLeaf.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names"));
  PlantStarchLeaf.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names")) ;
  PlantSugarSapwood.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names")) ;
  PlantStarchSapwood.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names")) ;
  PlantSugarTransport.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names")) ;
  SapwoodArea.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names")) ;
  LeafArea.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names")) ;
  FineRootArea.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names")) ;
  LAgrowth.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names")) ;
  SAgrowth.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names")) ;
  FRAgrowth.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names")) ;
  StemPI0.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  LeafPI0.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  // PlantLAIdead.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names")) ;
  // PlantLAIlive.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names")) ;
  // PlantLAIexpanded.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names")) ;
  
  subdailyRes.attr("names") = meteo.attr("row.names") ;
  
  NumericVector topo = NumericVector::create(elevation, slope, aspect);
  topo.attr("names") = CharacterVector::create("elevation", "slope", "aspect");
  
  Rcpp::DataFrame Stand = DataFrame::create(_["LAI"]=LAI, _["LAIlive"]=LAIlive, _["LAIexpanded"]=LAIexpanded,_["LAIdead"]=LAIdead,
                                            _["Cm"]=Cm, 
                                            _["LgroundPAR"] = LgroundPAR, _["LgroundSWR"] = LgroundSWR);
  Stand.attr("row.names") = meteo.attr("row.names");

  
  // Assemble output
  List plantCarbonBalance = List::create(
    Named("GrossPhotosynthesis") = GrossPhotosynthesis,
    Named("MaintenanceRespiration") = MaintenanceRespiration,
    Named("GrowthRespiration") = GrowthRespiration,
    Named("CarbonBalance") = CarbonBalance,
    Named("SugarLeaf") = PlantSugarLeaf,
    Named("StarchLeaf") = PlantStarchLeaf,
    Named("SugarSapwood") = PlantSugarSapwood,
    Named("StarchSapwood") = PlantStarchSapwood,
    Named("SugarTransport") = PlantSugarTransport,
    // Named("LabileMassLeaf") = LabileMassLeaf,
    // Named("LabileMassSapwood") = LabileMassSapwood,
    Named("LeafPI0") = LeafPI0,
    Named("StemPI0") = StemPI0
  );

  List plantGrowth;
  
  List l;
  if(transpirationMode=="Granier") {
    plantGrowth = List::create(Named("SapwoodArea")=SapwoodArea,
                               Named("LeafArea") = LeafArea,
                               Named("LAgrowth") = LAgrowth,
                               Named("SAgrowth") = SAgrowth);
    l = List::create(Named("latitude") = latitude,
                     Named("topography") = topo,
                     Named("growthInput") = growthInput,
                     Named("soilInput") = soilInput,
                     Named("WaterBalance")=DWB, 
                     Named("Soil")=SWB,
                     Named("Stand")=Stand,
                     Named("Plants") = plantDWOL,
                     Named("PlantCarbonBalance") = plantCarbonBalance,
                     Named("PlantGrowth") = plantGrowth,
                     Named("StandStructures") = standStructures,
                     Named("StandSummary") = standSummary,
                     Named("subdaily") =  subdailyRes);
    
  } else {
    plantGrowth = List::create(Named("SapwoodArea")=SapwoodArea,
                               Named("LeafArea") = LeafArea,
                               Named("FineRootArea") = FineRootArea,
                               Named("LAgrowth") = LAgrowth,
                               Named("SAgrowth") = SAgrowth,
                               Named("FRAgrowth") = FRAgrowth);
    l = List::create(Named("latitude") = latitude,
                   Named("topography") = topo,
                   Named("growthInput") = growthInput,
                   Named("soilInput") = soilInput,
                   Named("WaterBalance")=DWB, 
                   Named("EnergyBalance") = DEB,
                   Named("Temperature") = DT,
                   Named("Soil")=SWB,
                   Named("Stand")=Stand,
                   Named("Plants") = plantDWOL,
                   Named("SunlitLeaves") = sunlitDO,
                   Named("ShadeLeaves") = shadeDO,
                   Named("PlantCarbonBalance") = plantCarbonBalance,
                   Named("PlantGrowth") = plantGrowth,
                   Named("StandStructures") = standStructures,
                   Named("StandSummary") = standSummary,
                   Named("subdaily") =  subdailyRes);
  
  }
  l.attr("class") = CharacterVector::create("growth","list");
  return(l);
}
