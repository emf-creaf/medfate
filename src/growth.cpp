#include <numeric>
#include "lightextinction.h"
#include "phenology.h"
#include "biophysicsutils.h"
#include "forestutils.h"
#include "tissuemoisture.h"
#include "hydraulics.h"
#include "hydrology.h"
#include "carbon.h"
#include "soil.h"
#include "spwb.h"
#include <Rcpp.h>
#include <meteoland.h>
using namespace Rcpp;

//Ogle & Pacala 2010
//Tree Physiology 29, 587–605
const double leaf_RR = 0.00260274; // g gluc · g dw -1 · day -1
const double sapwood_RR = 6.849315e-05; // g gluc · g dw -1 · day -1
const double fineroot_RR = 0.002054795; // g gluc · g dw -1 · day -1
// const double leaf_RR = 0.0260274; // g gluc · g dw -1 · day -1
// const double sapwood_RR = 6.849315e-04; // g gluc · g dw -1 · day -1
// const double fineroot_RR = 0.02054795; // g gluc · g dw -1 · day -1
const double Q10_resp = 2.0;
//construction costs
const double leaf_CC = 1.5; // g gluc · g dw -1
const double sapwood_CC = 1.47; // g gluc · g dw -1
const double fineroots_CC = 1.30; // g gluc · g dw -1

// const double leaf_RR = 0.05; // g gluc · g dw -1 · day -1
// const double sapwood_RR = 0.005; // g gluc · g dw -1 · day -1
// const double fineroot_RR = 0.05; // g gluc · g dw -1 · day -1

//Maximum relative growth rate of leaves (1%/day)
const double RGRleafmax = 0.01; // m2·m-2·day-1

// minimum concentration (mol gluc·l-1) to ensure metabolism
const double equilibriumConcentrationLeaves = 0.1; 
const double equilibriumConcentrationSapwood = 0.1; 

//Ogle & Pacala 2010
//Tree Physiology 29, 587–605
const double dailySAturnoverProportion = 0.0001261398; //day-1 Equivalent to annual 4.5% 1-(1-0.045)^(1.0/365)
const double dailyFineRootTurnoverProportion = 0.001897231; //day-1 Equivalent to annual 4.5% 1-(1-0.5)^(1.0/365)



/**
 * floem flow (Holtta et al. 2017)
 *  psiUpstream, psiDownstream - water potential upstream (leaves)  and downstream
 *  concUpstream, concDownstream - sugar concentration upstream (leaves) and downstream (stem)
 *  k_f - floem conductance per leaf area basis (l*m-2*MPa-1*s-1)
 *  
 *  out mol*s-1*m-2 (flow per leaf area basis)
 */
// [[Rcpp::export("growth_floemFlow")]]
double floemFlow(double psiUpstream, double psiDownstream,
                 double concUpstream, double concDownstream,
                 double temp, double k_f = 3.0e-5) {
  double turgor_up = turgor(psiUpstream, concUpstream, temp);
  double turgor_down = turgor(psiDownstream, concDownstream, temp);
  if(temp < 0.0) k_f = 0.0; // No floem flow if temperature below zero
  if(turgor_up>turgor_down) {
    double relVisc = relativeSapViscosity(concUpstream, temp);
    return(k_f*concUpstream*(turgor_up - turgor_down)/relVisc);
  } else {
    double relVisc = relativeSapViscosity(concDownstream, temp);
    return(k_f*concDownstream*(turgor_up - turgor_down)/relVisc);
  }
}
// [[Rcpp::export("growth_dailyFloemFlow")]]
NumericMatrix dailyFloemFlow(List x, List spwbOut, 
                             NumericVector concLeaf, NumericVector concSapwood) {
  DataFrame paramStorage =  Rcpp::as<Rcpp::DataFrame>(x["paramsWaterStorage"]);
  NumericVector Vleaf = Rcpp::as<Rcpp::NumericVector>(paramStorage["Vleaf"]);
  NumericVector stemPI0 = Rcpp::as<Rcpp::NumericVector>(paramStorage["StemPI0"]);
  NumericVector leafPI0 = Rcpp::as<Rcpp::NumericVector>(paramStorage["LeafPI0"]);
  NumericVector stemEPS = Rcpp::as<Rcpp::NumericVector>(paramStorage["StemEPS"]);
  NumericVector leafEPS = Rcpp::as<Rcpp::NumericVector>(paramStorage["LeafEPS"]);
  List plantsInst = spwbOut["PlantsInst"];  
  NumericMatrix rwcStem =  Rcpp::as<Rcpp::NumericMatrix>(plantsInst["RWCstem"]);
  NumericMatrix rwcLeaf =  Rcpp::as<Rcpp::NumericMatrix>(plantsInst["RWCleaf"]);
  List eb = spwbOut["EnergyBalance"];  
  DataFrame tempDF =  Rcpp::as<Rcpp::DataFrame>(eb["Temperature"]);
  NumericVector Tcan = Rcpp::as<Rcpp::NumericVector>(tempDF["Tcan"]);
  
  int numCohorts = stemPI0.length();
  int numSteps = Tcan.length();
  NumericMatrix ff(numCohorts, numSteps);
  for(int c=0;c<numCohorts;c++) {
    for(int s=0;s<numSteps;s++) {
      // double leafPI = osmoticWaterPotential(concLeaf[c], Tcan[s]);
      // double sapwoodPI = osmoticWaterPotential(concs[c], Tcan[s]);
      double psiUp = symplasticWaterPotential(rwcLeaf(c,s), leafPI0[c], leafEPS[c]);
      double psiDown = symplasticWaterPotential(rwcStem(c,s), stemPI0[c], stemEPS[c]);
      ff(c,s) = floemFlow(psiUp, psiDown, concLeaf[c], concSapwood[c], Tcan[s])*3600.0; //flow as mol per hour and leaf area basis
    }
  }
  ff.attr("dimnames") = rwcStem.attr("dimnames");
  return(ff);
}

double qResp(double Tmean) {
  return(pow(Q10_resp,(Tmean-20.0)/10.0));
}

double storageTransferRelativeRate(double fastCstorage, double fastCstoragemax) {
  double f = ((2.0/(1.0+exp(-5.0*((fastCstorage/fastCstoragemax)-equilibriumConcentrationSapwood)/equilibriumConcentrationSapwood)))-1.0);
  return(f);
}
double temperatureGrowthFactor(double Tmean) {
  double Tlow = 5.0;
  double Thigh = 40.0;
  double Topt = 25.0;
  double f = ((Tmean-Tlow)*(Thigh-Tmean))/((Topt-Tlow)*(Thigh-Topt));
  return(std::min(std::max(f,0.0),1.0));
}
double turgorGrowthFactor(double psi, double psi_tlp) {
  return(std::max(0.0,1.0 - pow(exp((psi/psi_tlp)-1.0),5.0)));
}
double carbonGrowthFactor(double conc, double threshold) {
  double k =10.0;
  return(std::max(0.0,(1.0 - exp(k*(threshold-conc)))/(1.0 - exp(k*(-conc)))));
}
// [[Rcpp::export(".growth_defoliationFraction")]]
double defoliationFraction(double conc, double threshold) {
  double k =-10.0;
  return(std::max(0.0,(exp(k*conc)-exp(k*threshold))/(1.0-exp(k*threshold))));
}


/**
* Daily plant carbon balance
*/
DataFrame growthDay(List x, List spwbOut, double tday) {
  //Control params
  List control = x["control"];  
  
  String transpirationMode = control["transpirationMode"];
  String cavitationRefill = control["cavitationRefill"];
  bool taper = control["taper"];
  
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
  
  
  //Belowground parameters  
  List below = Rcpp::as<Rcpp::List>(x["below"]);
  NumericVector Z = Rcpp::as<Rcpp::NumericVector>(below["Z"]);
  
  //Internal state variables
  List internalWater = Rcpp::as<Rcpp::List>(x["internalWater"]);
  NumericVector NSPL = Rcpp::as<Rcpp::NumericVector>(internalWater["NSPL"]);
  NumericVector psiSympLeaf = Rcpp::as<Rcpp::NumericVector>(internalWater["psiSympLeaf"]);
  NumericVector psiSympStem = Rcpp::as<Rcpp::NumericVector>(internalWater["psiSympStem"]);
  
  List internalCarbon = Rcpp::as<Rcpp::List>(x["internalCarbon"]);
  NumericVector sugarLeaf = internalCarbon["sugarLeaf"];
  NumericVector starchLeaf = internalCarbon["starchLeaf"];
  NumericVector sugarSapwood = internalCarbon["sugarSapwood"];
  NumericVector starchSapwood = internalCarbon["starchSapwood"];
  NumericVector longtermStorage = internalCarbon["longtermStorage"];
  
  
  List stand = spwbOut["Stand"];
  List Plants = spwbOut["Plants"];
  List PlantsInst = spwbOut["PlantsInst"];
  
  //Recover module-communication state variables
  NumericMatrix AgStep  =  Rcpp::as<Rcpp::NumericMatrix>(PlantsInst["Ag"]);
  int numSteps = AgStep.ncol();
  
  NumericMatrix rwcStem =  Rcpp::as<Rcpp::NumericMatrix>(PlantsInst["RWCstem"]);
  NumericMatrix rwcLeaf =  Rcpp::as<Rcpp::NumericMatrix>(PlantsInst["RWCleaf"]);
  
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
  //Growth parameters
  DataFrame paramsGrowth = Rcpp::as<Rcpp::DataFrame>(x["paramsGrowth"]);
  NumericVector WoodC = Rcpp::as<Rcpp::NumericVector>(paramsGrowth["WoodC"]);
  NumericVector RGRmax = Rcpp::as<Rcpp::NumericVector>(paramsGrowth["RGRmax"]);
  NumericVector leafDuration = Rcpp::as<Rcpp::NumericVector>(paramsGrowth["leafDuration"]);
  // NumericVector Cstoragepmax= Rcpp::as<Rcpp::NumericVector>(paramsGrowth["Cstoragepmax"]);
  // NumericVector slowCstorage_max(numCohorts), fastCstorage_max(numCohorts);
  //Transpiration parameters
  DataFrame paramsTransp = Rcpp::as<Rcpp::DataFrame>(x["paramsTransp"]);
  NumericVector Kmax_stemxylem, VCstem_kmax, Psi_Extract, VCstem_c, VCstem_d;
  NumericMatrix VGrhizo_kmax, VCroot_kmax;
  Kmax_stemxylem = paramsTransp["Kmax_stemxylem"];
  VCstem_kmax = paramsTransp["VCstem_kmax"];
  VGrhizo_kmax = Rcpp::as<Rcpp::NumericMatrix>(below["VGrhizo_kmax"]);
  VCroot_kmax = Rcpp::as<Rcpp::NumericMatrix>(below["VCroot_kmax"]);
  int numLayers = VCroot_kmax.ncol();
  
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
  
  //Output vectors
  NumericVector MaintenanceRespiration(numCohorts,0.0);
  NumericVector GrowthRespiration(numCohorts,0.0);
  NumericVector PlantSugarTransport(numCohorts,0.0);
  NumericVector PlantSugarLeaf(numCohorts,0.0);
  NumericVector PlantStarchLeaf(numCohorts,0.0);
  NumericVector PlantSugarSapwood(numCohorts,0.0);
  NumericVector PlantStarchSapwood(numCohorts,0.0);
  NumericVector PlantSA(numCohorts,0.0);
  NumericVector PlantSAgrowth(numCohorts,0.0);
  NumericVector GrossPhotosynthesis(numCohorts,0.0);
  NumericVector PlantLAIdead(numCohorts,0.0);
  NumericVector PlantLAIlive(numCohorts,0.0);
  

  //Storage volume and maximum starch capacity for leaves and sapwood  
  NumericVector Volume_leaves(numCohorts,0.0);
  NumericVector Volume_sapwood(numCohorts,0.0);
  NumericVector Starch_max_leaves(numCohorts,0.0);
  NumericVector Starch_max_sapwood(numCohorts,0.0);
  NumericVector B_struct_leaves(numCohorts,0.0);
  NumericVector B_struct_sapwood(numCohorts,0.0);
  NumericVector B_struct_fineroots(numCohorts,0.0);

  //3. Carbon balance and growth
  for(int j=0;j<numCohorts;j++){

    Volume_leaves[j] = leafStorageVolume(LAI_expanded[j],  N[j], SLA[j], LeafDensity[j]);
    Volume_sapwood[j] = sapwoodStorageVolume(SA[j], H[j],Z[j],WoodDensity[j], 0.5);
    Starch_max_leaves[j] = leafStarchCapacity(LAI_expanded[j],  N[j], SLA[j], 0.3)/Volume_leaves[j];
    Starch_max_sapwood[j] = sapwoodStarchCapacity(SA[j], H[j],Z[j],WoodDensity[j], 0.2)/Volume_sapwood[j];
    B_struct_leaves[j] = leafStructuralBiomass(LAI_expanded[j],N[j],SLA[j]);
    B_struct_sapwood[j] = sapwoodStructuralLivingBiomass(SA[j], H[j], Z[j], WoodDensity[j], 0.5);
    B_struct_fineroots[j] = B_struct_leaves[j]/2.0; //TO BE CHANGED
    
    // Rcout << j << " Lvol: "<< Volume_leaves[j] << " Svol: "<<Volume_sapwood[j]<< " LStarchMax: "<<Starch_max_leaves[j]
    //       << " SStarchMax: "<<Starch_max_sapwood[j]<< " Bleaf "<< B_struct_leaves[j]<< " Bsap "<< B_struct_sapwood[j]<< " Bfr "<< B_struct_fineroots[j]<<"\n";
    
    double LAexpanded = leafArea(LAI_expanded[j], N[j]);
    double LAlive = leafArea(LAI_live[j], N[j]);
    double LAdead = leafArea(LAI_dead[j], N[j]);
    double LAlive_ini = LAlive;
    
    double leafRespDay = 0.0;
    double sfrRespDay = 0.0;
    //Carbon balance for labile carbon of leaves and stems
    for(int s=0;s<numSteps;s++) {
      
      //transform concentrations relative to RWC = 1 to concentrations relative to RWC of the step
      double sugarLeafStep = sugarLeaf[j]/rwcLeaf(j,s);
      double sugarSapwoodStep = sugarSapwood[j]/rwcStem(j,s);
      double starchLeafStep = starchLeaf[j]/rwcLeaf(j,s);
      double starchSapwoodStep = starchSapwood[j]/rwcStem(j,s);
      
      //Transform sugar concentration (mol gluc · l-1) to sugar mass (g gluc)
      // double lstvol = 0.001*(starchLeaf[j]/starchDensity);
      // double sstvol = 0.001*(starchSapwood[j]/starchDensity);
      double leafSugarMassStep = sugarLeafStep*(Volume_leaves[j]*glucoseMolarMass);
      double sapwoodSugarMassStep = sugarSapwoodStep*(Volume_sapwood[j]*glucoseMolarMass);
      
      //Respiratory biomass (g dw · ind-1)
      double B_resp_leaves = B_struct_leaves[j] + leafSugarMassStep;
      double B_resp_sapwood = B_struct_sapwood[j] + sapwoodSugarMassStep;
      double B_resp_fineroots = B_struct_fineroots[j];
      double QR = qResp(Tcan[s]);
      double leafRespStep = B_resp_leaves*leaf_RR*QR/((double) numSteps);
      double sfrRespStep = (B_resp_sapwood*sapwood_RR + B_resp_fineroots*fineroot_RR)*QR/((double) numSteps);
      leafRespDay +=leafRespStep;
      sfrRespDay += sfrRespStep;
      
      //gross fotosynthesis
      double leafAgStepC = AgStep(j,s)/(N[j]/10000.0); //Translate g C · m-2 · h-1 to g C · h-1
      double leafAgStepG = leafAgStepC*(glucoseMolarMass/(carbonMolarMass*6.0)); // from g C· h-1 to g gluc · h-1
      
      GrossPhotosynthesis[j] += leafAgStepG; //Ag in g gluc
      MaintenanceRespiration[j] += (sfrRespStep + leafRespStep); //Rm in g gluc
      
      //sugar mass balance
      leafSugarMassStep = leafSugarMassStep + leafAgStepG - leafRespStep;
      sapwoodSugarMassStep = sapwoodSugarMassStep - sfrRespStep;
      // Rcout<<" coh:"<<j<< " s:"<<s<<" A: "<< leafAgStepG << "Rl: " << leafRespStep<<" sugar mass leaf: "<< leafSugarMassStep << " Rs"<<  sfrRespStep<< " sugar mass sap: "<< sapwoodSugarMassStep<<"\n";
      
      //re-transform to concentration (mol gluc · l-1)
      sugarLeafStep = leafSugarMassStep/(Volume_leaves[j]*glucoseMolarMass);
      sugarSapwoodStep = sapwoodSugarMassStep/(Volume_sapwood[j]*glucoseMolarMass);
      
      //floem transport      
      double psiUp = symplasticWaterPotential(rwcLeaf(j,s), LeafPI0[j], LeafEPS[j]);
      double psiDown = symplasticWaterPotential(rwcStem(j,s), StemPI0[j], StemEPS[j]);
      
      double ff = 0.0;
      for(int t=0;t<3600;t++) {
        double ft = floemFlow(psiUp, psiDown, sugarLeafStep, sugarSapwoodStep, Tcan[s])*LAlive; //flow as mol glucose per s
        // sugar-starch dynamics
        double conversionLeaf = sugarStarchDynamicsLeaf(sugarLeafStep, starchLeafStep, equilibriumConcentrationLeaves);
        double conversionSapwood = sugarStarchDynamicsStem(sugarSapwoodStep, starchSapwoodStep, equilibriumConcentrationSapwood);
        // Rcout<<" coh:"<<j<< " s:"<<s<< " Lsugar: "<< sugarLeafStep << " Lstarch: "<< sugarSapwoodStep<<" starch formation: "<<conversionLeaf<< "\n";
        starchLeafStep  += conversionLeaf;
        starchSapwoodStep += conversionSapwood;
        // Rcout<<" coh:"<<j<< " s:"<<s<< " Ssugar: "<< sugarSapwoodStep << " Sstarch: "<< starchSapwoodStep<<" starch formation: "<<conversionSapwood<< "\n";
        //Apply floem transport (mol gluc) to sugar concentrations (mol gluc· l-1)
        sugarLeafStep = ((sugarLeafStep*Volume_leaves[j]) - ft)/Volume_leaves[j] - conversionLeaf;
        sugarSapwoodStep = ((sugarSapwoodStep*Volume_sapwood[j]) + ft)/Volume_sapwood[j] - conversionSapwood;
        // sugarLeafStep = ((sugarLeafStep*Volume_leaves[j]) - ft)/Volume_leaves[j];
        // sugarSapwoodStep = ((sugarSapwoodStep*Volume_sapwood[j]) + ft)/Volume_sapwood[j];
        
        ff +=ft;
      }
      PlantSugarTransport[j] += ff; //To calculate daily floem balance (positive means towards stem)
      // Rcout<<" coh:"<<j<< " s:"<<s<< " conc leaf: "<< sugarLeafStep << " conc sap: "<< sugarSapwoodStep<<" ff: "<<ff<< "\n";
      
      
      sugarLeaf[j] = sugarLeafStep*rwcLeaf(j,s);
      starchLeaf[j] = starchLeafStep*rwcLeaf(j,s);
      sugarSapwood[j] = sugarSapwoodStep*rwcStem(j,s);
      starchSapwood[j] = starchSapwoodStep*rwcStem(j,s);
    }
    //Check starch limits
    // if(starchLeaf[j]> Starch_max_leaves[j]) {
    //   longtermStorage[j] += (starchLeaf[j]-Starch_max_leaves[j])*Volume_leaves[j];
    //   starchLeaf[j] = Starch_max_leaves[j];
    // }
    // if(starchSapwood[j]> Starch_max_sapwood[j]) {
    //   longtermStorage[j] += (starchSapwood[j]-Starch_max_sapwood[j])*Volume_sapwood[j];
    //   starchSapwood[j] = Starch_max_sapwood[j];
    // }
    
    if(sugarLeaf[j] < 0.0) { //Leaf senescense due to C starvation
      double respirationExcess = -sugarLeaf[j]*(Volume_leaves[j]*glucoseMolarMass); //g gluc
      double propExcess = respirationExcess/leafRespDay; //day
      Rcout<< j <<" Excess respiration: " << respirationExcess << " Prop:"<< propExcess<< " LAlive " << LAlive << " LAlivenew "<< LAlive*(1.0 - propExcess) <<"\n";
      LAdead = LAdead + LAexpanded*propExcess;
      LAexpanded = LAexpanded*(1.0 - propExcess);
      LAlive = LAlive*(1.0 - propExcess);
      sugarLeaf[j] = 0.0;
    } else { //Leaf senescence due to age (Ca+ accumulation) (should change with better phenology modelling)
      double propAged = (1.0/(365.25*leafDuration[j]));
      LAdead = LAdead + LAexpanded*propAged;
      LAexpanded = LAexpanded*(1.0 - propAged); 
      LAlive = LAlive*(1.0 - propAged); 
    }
    
    //Leaf growth
    double f_temp = temperatureGrowthFactor(tday);
    double fLA_turgor = turgorGrowthFactor(psiSympLeaf[j],turgorLossPoint(LeafPI0[j], LeafEPS[j]));
    if(fLA_turgor>0.0) {
      double costLA = leaf_CC/(1000.0*SLA[j]); // g gluc · m-2 of leaf area
      double fLA_sinkC = carbonGrowthFactor(sugarLeaf[j], equilibriumConcentrationLeaves);
      double deltaLAsink = RGRleafmax*LAlive*f_temp*fLA_turgor*fLA_sinkC;
      double deltaLAavailable = (sugarLeaf[j]*(glucoseMolarMass*Volume_leaves[j]))/costLA;
      double deltaLAgrowth = std::min(deltaLAsink, deltaLAavailable);
      double growthCost = deltaLAgrowth*costLA;
      sugarLeaf[j] = sugarLeaf[j] - (growthCost/(glucoseMolarMass*Volume_leaves[j]));
      LAlive += deltaLAgrowth; //Update leaf area
      LAexpanded +=deltaLAgrowth;
      GrowthRespiration[j] +=growthCost; //growth cost in g gluc
      // Rcout<< j << " fLAturgor: "<< fLA_turgor<< " fLA_sinkC "<< fLA_sinkC<< "f_temp"<< f_temp<<"\n";
    }

    //SA growth senescense
    double deltaSAturnover = (dailySAturnoverProportion/(1.0+15.0*exp(-0.01*H[j])))*SA[j];
    SA[j] = SA[j] - deltaSAturnover; //Update sapwood area
    
    //SA growth     
         
    
    //Update LAI
    LAI_live[j] = LAlive*N[j]/10000.0;
    LAI_expanded[j] = LAexpanded*N[j]/10000.0;
    LAI_dead[j] = LAdead*N[j]/10000.0;
    //Update Huber value, stem and root hydraulic conductance
    Al2As[j] = (LAlive)/(SA[j]/10000.0);
    VCstem_kmax[j]=maximumStemHydraulicConductance(Kmax_stemxylem[j], Hmed[j], Al2As[j] ,H[j], taper);
    for(int s=0;s<numLayers;s++) {
      VCroot_kmax(j,s) = VCroot_kmax(j,s)*(LAlive_ini/LAlive);
    }     
    
    //Update leaf and stem osmotic water potential at maximum turgor
    LeafPI0[j] = osmoticWaterPotential(sugarLeaf[j], tday)/rwcLeaf(j, numSteps-1);
    StemPI0[j] = osmoticWaterPotential(sugarSapwood[j], tday)/rwcStem(j, numSteps-1);
    
    //Update non-stomatal photosynthesis limitations
    NSPL[j] = 1.0 - std::max(0.0, std::min(1.0, sugarLeaf[j] - 0.5)); //photosynthesis limited when conc > 0.5 and zero when conc > 1.5 mol·l-1
      
    //Output variables
    PlantSugarLeaf[j] = sugarLeaf[j];
    PlantStarchLeaf[j] = starchLeaf[j];
    PlantSugarSapwood[j] = sugarSapwood[j];
    PlantStarchSapwood[j] = starchSapwood[j];
    PlantSA[j] = SA[j];
    PlantLAIlive[j] = LAI_live[j];
    PlantLAIdead[j] = LAI_dead[j];
  }
  
  DataFrame df = DataFrame::create(_["GrossPhotosynthesis"] = GrossPhotosynthesis,
                                   _["MaintenanceRespiration"] = MaintenanceRespiration,
                                   _["GrowthRespiration"] = GrowthRespiration,
                                   _["PlantSugarLeaf"] = PlantSugarLeaf,
                                   _["PlantStarchLeaf"] = PlantStarchLeaf,
                                   _["PlantSugarSapwood"] = PlantSugarSapwood,
                                   _["PlantStarchSapwood"] = PlantStarchSapwood,
                                   _["PlantSugarTransport"] = PlantSugarTransport,
                                   _["PlantSA"] = PlantSA,
                                   _["PlantLAIlive"] = PlantLAIlive,
                                   _["PlantLAIdead"] = PlantLAIdead,
                                   _["PlantSAgrowth"] = PlantSAgrowth);
  
  df.attr("row.names") = cohorts.attr("row.names");
  return(df);
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
  List below = Rcpp::as<Rcpp::List>(x["below"]);
  if(!below.containsElementNamed("Z")) stop("Z missing in growthInput$below");
  if(!below.containsElementNamed("V")) stop("V missing in growthInput$below");
  if(transpirationMode=="Sperry"){
    if(!below.containsElementNamed("VGrhizo_kmax")) stop("VGrhizo_kmax missing in growthInput$below");
    if(!below.containsElementNamed("VCroot_kmax")) stop("VCroot_kmax missing in growthInput$below");
  }  
  
  if(!x.containsElementNamed("paramsBase")) stop("paramsBase missing in growthInput");
  DataFrame paramsBase = Rcpp::as<Rcpp::DataFrame>(x["paramsBase"]);
  if(!paramsBase.containsElementNamed("Sgdd")) stop("Sgdd missing in growthInput$paramsBase");
  if(!paramsBase.containsElementNamed("kPAR")) stop("kPAR missing in growthInput$paramsBase");
  if(!paramsBase.containsElementNamed("g")) stop("g missing in growthInput$paramsBase");
  
  if(!x.containsElementNamed("paramsGrowth")) stop("paramsGrowth missing in growthInput");
  DataFrame paramsGrowth = Rcpp::as<Rcpp::DataFrame>(x["paramsGrowth"]);
  if(!paramsGrowth.containsElementNamed("WoodC")) stop("WoodC missing in growthInput$paramsGrowth");
  if(!paramsGrowth.containsElementNamed("Cstoragepmax")) stop("Cstoragepmax missing in growthInput$paramsGrowth");
  if(!paramsGrowth.containsElementNamed("RGRmax")) stop("RGRmax missing in growthInput$paramsGrowth");
  
  if(!x.containsElementNamed("paramsAnatomy")) stop("paramsAnatomy missing in growthInput");
  DataFrame paramsAnatomy = Rcpp::as<Rcpp::DataFrame>(x["paramsAnatomy"]);
  if(!paramsAnatomy.containsElementNamed("SLA")) stop("SLA missing in paramsAnatomy$paramsGrowth");
  if(!paramsAnatomy.containsElementNamed("Al2As")) stop("Al2As missing in paramsAnatomy$paramsGrowth");
  if(!paramsAnatomy.containsElementNamed("WoodDensity")) stop("WoodDensity missing in paramsAnatomy$paramsGrowth");
  
  if(!x.containsElementNamed("paramsTransp")) stop("paramsTransp missing in growthInput");
  DataFrame paramsTransp = Rcpp::as<Rcpp::DataFrame>(x["paramsTransp"]);
  // if(!paramsTransp.containsElementNamed("pRootDisc")) stop("pRootDisc missing in growthInput$paramsTransp");
  if(transpirationMode=="Granier") {
    if(!paramsTransp.containsElementNamed("Psi_Extract")) stop("Psi_Extract missing in growthInput$paramsTransp");
    if(!paramsTransp.containsElementNamed("WUE")) stop("WUE missing in growthInput$paramsTransp");
  } else if(transpirationMode=="Sperry") {
    if(!soil.containsElementNamed("VG_n")) stop("VG_n missing in soil");
    if(!soil.containsElementNamed("VG_alpha")) stop("VG_alpha missing in soil");
    
    if(!paramsTransp.containsElementNamed("VCstem_kmax")) stop("VCstem_kmax missing in growthInput");
    if(!paramsTransp.containsElementNamed("VCstem_c")) stop("VCstem_c missing in growthInput$paramsTransp");
    if(!paramsTransp.containsElementNamed("VCstem_d")) stop("VCstem_d missing in growthInput$paramsTransp");
    if(!paramsTransp.containsElementNamed("VCroot_c")) stop("VCroot_c missing in growthInput$paramsTransp");
    if(!paramsTransp.containsElementNamed("VCroot_d")) stop("VCroot_d missing in growthInput$paramsTransp");
    if(!paramsTransp.containsElementNamed("Gwmax")) stop("Gwmax missing in growthInput$paramsTransp");
    if(!paramsTransp.containsElementNamed("Vmax298")) stop("Vmax298 missing in growthInput$paramsTransp");
    if(!paramsTransp.containsElementNamed("Jmax298")) stop("Jmax298 missing in growthInput$paramsTransp");
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

void recordStandSummary(DataFrame standSummary, NumericVector LAIlive,
                        NumericVector N, 
                        NumericVector DBH,  NumericVector Cover, 
                        NumericVector H, int pos) {
  NumericVector SLAI = as<Rcpp::NumericVector>(standSummary["LeafAreaIndex"]);
  SLAI[pos] = 0.0;
  NumericVector TBA = as<Rcpp::NumericVector>(standSummary["TreeBasalArea"]);
  TBA[pos] = 0.0;
  NumericVector TDensity = as<Rcpp::NumericVector>(standSummary["TreeDensity"]);
  TDensity[pos] = 0.0;
  NumericVector SCover = as<Rcpp::NumericVector>(standSummary["ShrubCover"]);
  SCover[pos] = 0.0;
  NumericVector MaxHeight = as<Rcpp::NumericVector>(standSummary["MaxHeight"]);
  MaxHeight[pos] = 0.0;
  int numCohorts = N.length();
  NumericVector treeBA = treeBasalArea(N, DBH);
  for(int i=0;i<numCohorts;i++) {
    SLAI[pos] += LAIlive[i];
    if(!NumericVector::is_na(treeBA[i])) {
      TBA[pos] += treeBA[i];
      TDensity[pos] +=N[i];
      MaxHeight[pos] = std::max(MaxHeight[pos], H[i]);
    } else {
      SCover[pos] +=Cover[i];
    }
  }
}

// DataFrame growthPipe(List x, List spwbOut, double tday) {
//   //Control params
//   List control = x["control"];  
//   
//   String transpirationMode = control["transpirationMode"];
//   String cavitationRefill = control["cavitationRefill"];
//   bool taper = control["taper"];
//   
//   //Cohort info
//   DataFrame cohorts = Rcpp::as<Rcpp::DataFrame>(x["cohorts"]);
//   NumericVector SP = cohorts["SP"];
//   int numCohorts = SP.size();
//   
//   //Aboveground state variables  
//   DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
//   NumericVector DBH = above["DBH"];
//   NumericVector Cover = above["Cover"];
//   NumericVector H = above["H"];
//   NumericVector N = above["N"];
//   NumericVector CR = above["CR"];
//   NumericVector LAI_live = above["LAI_live"];
//   NumericVector LAI_expanded = above["LAI_expanded"];
//   NumericVector LAI_dead = above["LAI_dead"];
//   NumericVector SA = above["SA"];
//   
//   //Belowground parameters  
//   List below = Rcpp::as<Rcpp::List>(x["below"]);
//   NumericVector Z = Rcpp::as<Rcpp::NumericVector>(below["Z"]);
//   
//   //Internal state variables
//   List internalWater = Rcpp::as<Rcpp::List>(x["internalWater"]);
//   List internalCarbon = Rcpp::as<Rcpp::List>(x["internalCarbon"]);
//   NumericVector fastCstorage = internalCarbon["fastCstorage"];
//   NumericVector slowCstorage = internalCarbon["slowCstorage"];
//   
//   
//   List stand = spwbOut["Stand"];
//   List Plants = spwbOut["Plants"];
//   
//   //Recover module-communication state variables
//   NumericVector Ag, psiCoh;
//   if(transpirationMode=="Granier") {
//     Ag =  Rcpp::as<Rcpp::NumericVector>(Plants["Photosynthesis"]);
//     psiCoh =  Rcpp::as<Rcpp::NumericVector>(Plants["psi"]);  
//   } else {
//     Ag =  Rcpp::as<Rcpp::NumericVector>(Plants["GrossPhotosynthesis"]);
//     psiCoh =  clone(Rcpp::as<Rcpp::NumericVector>(x["psiLeaf"]));  
//   }
//   
//   //Anatomy parameters
//   DataFrame paramsAnatomy = Rcpp::as<Rcpp::DataFrame>(x["paramsAnatomy"]);
//   NumericVector SLA = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["SLA"]);
//   NumericVector Al2As = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["Al2As"]);
//   NumericVector WoodDensity = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["WoodDensity"]);
//   //Growth parameters
//   DataFrame paramsGrowth = Rcpp::as<Rcpp::DataFrame>(x["paramsGrowth"]);
//   NumericVector WoodC = Rcpp::as<Rcpp::NumericVector>(paramsGrowth["WoodC"]);
//   NumericVector RGRmax = Rcpp::as<Rcpp::NumericVector>(paramsGrowth["RGRmax"]);
//   NumericVector Cstoragepmax= Rcpp::as<Rcpp::NumericVector>(paramsGrowth["Cstoragepmax"]);
//   NumericVector slowCstorage_max(numCohorts), fastCstorage_max(numCohorts);
//   
//   
//   //Transpiration parameters
//   DataFrame paramsTransp = Rcpp::as<Rcpp::DataFrame>(x["paramsTransp"]);
//   NumericVector Kmax_stemxylem, VCstem_kmax, Psi_Extract, VCstem_c, VCstem_d;
//   if(transpirationMode=="Sperry") {
//     Kmax_stemxylem = paramsTransp["Kmax_stemxylem"];
//     VCstem_kmax = paramsTransp["VCstem_kmax"];
//     // VCstem_c = paramsTransp["VCstem_c"];
//     // VCstem_d = paramsTransp["VCstem_d"];
//   } else if(transpirationMode == "Granier"){
//     // Psi_Extract = paramsTransp["Psi_Extract"];
//   }
//   
//   //Output vectors
//   NumericVector PlantRespiration(numCohorts,0.0);
//   NumericVector PlantCstorageFast(numCohorts,0.0);
//   NumericVector PlantCstorageSlow(numCohorts,0.0);
//   NumericVector PlantSA(numCohorts,0.0);
//   NumericVector PlantSAgrowth(numCohorts,0.0);
//   NumericVector PlantGrossPhotosynthesis(numCohorts,0.0);
//   NumericVector PlantLAIdead(numCohorts,0.0);
//   NumericVector PlantLAIlive(numCohorts,0.0);
//   
//   //3. Carbon balance and growth
//   double B_leaf_expanded, B_stem, B_fineroot;
//   for(int j=0;j<numCohorts;j++){
//     //3.1 Live biomass and maximum C pool capacity
//     NumericVector compartments = carbonCompartments(SA[j], LAI_expanded[j], H[j], Z[j], N[j], SLA[j], WoodDensity[j], WoodC[j]);
//     B_leaf_expanded = compartments[0];
//     B_stem = compartments[1];
//     B_fineroot = compartments[2];
//     fastCstorage_max[j] = 0.05*(B_leaf_expanded+B_stem+B_fineroot);
//     slowCstorage_max[j] = std::max(slowCstorage_max[j],(Cstoragepmax[j]-0.05)*(B_leaf_expanded+B_stem+B_fineroot)); //Slow pool capacity cannot decrease
//     
//     //3.2 Respiration and photosynthesis 
//     double Agj = Ag[j]/(N[j]/10000.0); //Translate g C · m-2 to g C · ind-1
//     // double Anj = 0.0;
//     double QR = qResp(tday);
//     double Rj = (B_leaf_expanded*leaf_RR + B_stem*stem_RR + B_fineroot*root_RR)*QR;
//     
//     //3.3. Carbon balance, update of fast C pool and C available for growth
//     double growthAvailableC = 0.0;
//     growthAvailableC = std::max(0.0,fastCstorage[j]+(Agj-Rj));
//     fastCstorage[j] = growthAvailableC;
//     
//     //3.4 Growth in LAI_live and SA
//     double deltaSAturnover = (dailySAturnoverProportion/(1.0+15*exp(-0.01*H[j])))*SA[j];
//     double f_turgor = turgorGrowthFactor(psiCoh[j],-1.5);
//     double deltaSAgrowth = 0.0;
//     // if(f_turgor>0.0) { //Growth is possible
//     double costLA = 0.1*leafCperDry*(Al2As[j]/SLA[j]); //Construction cost in g C·cm-2 of sapwood
//     double costSA = WoodC[j]*(H[j]+Z[j])*WoodDensity[j];  //Construction cost in g C·cm-2 of sapwood
//     double costFR = costLA/2.5;
//     double cost = 1.3*(costLA+costSA+costFR);  //Construction cost in g C·cm-2 of sapwood (including 30% growth respiration)
//     double deltaSAavailable = growthAvailableC/cost;
//     double f_source = 1.0;
//     f_source = carbonGrowthFactor(fastCstorage[j]/fastCstorage_max[j], growthCarbonConcentrationThreshold);
//     double f_temp = temperatureGrowthFactor(tday);
//     double deltaSAsink = RGRmax[j]*SA[j]*f_temp*f_turgor;
//     deltaSAgrowth = std::min(deltaSAsink*f_source, deltaSAavailable);
//     
//     //update pools
//     fastCstorage[j] = fastCstorage[j]-deltaSAgrowth*cost; //Remove construction costs from (fast) C pool
//     if(transpirationMode=="Granier"){
//       if(cavitationRefill!="total") { //If we track cavitation update proportion of embolized conduits
//         NumericVector pEmb =  Rcpp::as<Rcpp::NumericVector>(x["PLC"]);
//         pEmb[j] = pEmb[j]*((SA[j] - deltaSAturnover)/(SA[j] + deltaSAgrowth - deltaSAturnover));
//       }
//     } else {
//       NumericMatrix PLCstem =  Rcpp::as<Rcpp::NumericMatrix>(x["PLCstem"]);
//       int nStemSegments = PLCstem.ncol();
//       for(int s=0;s<nStemSegments;s++) {
//         PLCstem(j,s) = PLCstem(j,s)*((SA[j] - deltaSAturnover)/(SA[j] + deltaSAgrowth - deltaSAturnover));
//       }
//       
//     }
//     SA[j] = SA[j] + deltaSAgrowth - deltaSAturnover; //Update sapwood area
//     
//     double leafDie = std::min((N[j]/10000.0)*(deltaSAturnover/10000.0)*Al2As[j], LAI_live[j]);
//     LAI_dead[j] += leafDie; //Update dead LAI
//     LAI_live[j] += (N[j]/10000.0)*((deltaSAgrowth-deltaSAturnover)/10000.0)*Al2As[j]; //Update live LAI
//     LAI_live[j] = std::max(LAI_live[j], 0.0); //Check negative values do not occur
//     
//     //3.5 transfer between pools and constrain of C pools
//     //Relative transfer rate (maximum 5% of the source pool per day)
//     double reltransferRate = 0.05*storageTransferRelativeRate(fastCstorage[j], fastCstorage_max[j]);
//     if(reltransferRate>0.0) { //Transfer from fast to slow 
//       double transfer = std::min(reltransferRate*fastCstorage[j],(slowCstorage_max[j]-slowCstorage[j])*0.9);
//       fastCstorage[j] -= transfer;
//       slowCstorage[j] += transfer*0.9; //10% cost in respiration (not added to slow pool)
//     } else { //Transfer from slow to fast 
//       double transfer = std::min(-reltransferRate*slowCstorage[j],(fastCstorage_max[j]-fastCstorage[j])*0.9);
//       fastCstorage[j] += transfer*0.9; //10% cost in respiration (removed from what actually reaches fast pool)
//       slowCstorage[j] -= transfer;
//     }
//     //Trim pools to maximum capacity
//     fastCstorage[j] = std::max(0.0,std::min(fastCstorage[j], fastCstorage_max[j]));
//     slowCstorage[j] = std::max(0.0,std::min(slowCstorage[j], slowCstorage_max[j]));
//     
//     //3.8 Calculate defoliation if fast storage is low
//     double def = defoliationFraction(fastCstorage[j]/fastCstorage_max[j], growthCarbonConcentrationThreshold);
//     double maxLAI = (SA[j]*Al2As[j]/10000.0)*(N[j]/10000.0);
//     double defLAI = maxLAI * (1.0-def);
//     // Rcout<<defLAI<<" "<<LAI_live[j]<<"\n";
//     LAI_live[j] = std::min(LAI_live[j], defLAI);
//     
//     //3.7 Update stem conductance (Sperry mode)
//     if(transpirationMode=="Sperry") {
//       double al2as = (LAI_expanded[j]/(N[j]/10000.0))/(SA[j]/10000.0);
//       VCstem_kmax[j]=maximumStemHydraulicConductance(Kmax_stemxylem[j], al2as,H[j], taper);
//       // Rcout<<Al2As[j]<<" "<< al2as<<" "<<VCstem_kmax[j]<<"\n";
//     }
//     
//     //Output variables
//     PlantRespiration[j] = Rj*(N[j]/10000.0); //Scaled to cohort level: Translate g C · ind-1 to g C · m-2
//     PlantCstorageFast[j] = fastCstorage[j]/fastCstorage_max[j];
//     PlantCstorageSlow[j] = slowCstorage[j]/slowCstorage_max[j];
//     PlantSA[j] = SA[j];
//     PlantLAIlive[j] = LAI_live[j];
//     PlantLAIdead[j] = LAI_dead[j];
//     PlantSAgrowth[j] = deltaSAgrowth;
//   }
//   
//   DataFrame df = DataFrame::create(_["PlantRespiration"] = PlantRespiration,
//                                    _["PlantCstorageFast"] = PlantCstorageFast,
//                                    _["PlantCstorageSlow"] = PlantCstorageSlow,
//                                    _["PlantSA"] = PlantSA,
//                                    _["PlantLAIlive"] = PlantLAIlive,
//                                    _["PlantLAIdead"] = PlantLAIdead,
//                                    _["PlantSAgrowth"] = PlantSAgrowth);
//   
//   df.attr("row.names") = cohorts.attr("row.names");
//   return(df);
// }

// [[Rcpp::export("growth")]]
List growth(List x, List soil, DataFrame meteo, double latitude = NA_REAL, double elevation = NA_REAL, double slope = NA_REAL, double aspect = NA_REAL) {
  //Control params
  List control = x["control"];  
  
  //Store input
  List spwbInput = clone(x);
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
  String cavitationRefill = control["cavitationRefill"];
  checkgrowthInput(x, soil, transpirationMode, soilFunctions);
  
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
    if(NumericVector::is_na(latitude)) stop("Value for 'latitude' should not be missing.");
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
  
  //Canpopy parameters
  List canopyParams = x["canopy"];
  
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
  int numCohorts = SP.size();

  //Belowground state variables  
  List below = Rcpp::as<Rcpp::List>(x["below"]);
  NumericVector Z = Rcpp::as<Rcpp::NumericVector>(below["Z"]);

  //Internal state variables
  DataFrame internalWater = Rcpp::as<Rcpp::List>(x["internalWater"]);
  DataFrame internalCarbon = Rcpp::as<Rcpp::List>(x["internalCarbon"]);
  

  //Base parameters
  DataFrame paramsBase = Rcpp::as<Rcpp::DataFrame>(x["paramsBase"]);
  NumericVector Sgdd = paramsBase["Sgdd"];
  NumericVector kPAR = paramsBase["kPAR"];
  
  // Rcout<<"3";

  //Transpiration parameters
  DataFrame paramsTransp = Rcpp::as<Rcpp::DataFrame>(x["paramsTransp"]);
  NumericVector Kmax_stemxylem, VCstem_kmax, Psi_Extract, VCstem_c, VCstem_d;
  if(transpirationMode=="Sperry") {
    Kmax_stemxylem = paramsTransp["Kmax_stemxylem"];
    VCstem_kmax = paramsTransp["VCstem_kmax"];
    VCstem_c = paramsTransp["VCstem_c"];
    VCstem_d = paramsTransp["VCstem_d"];
  } else if(transpirationMode == "Granier"){
    Psi_Extract = paramsTransp["Psi_Extract"];
  }

  //Soil input
  NumericVector Water_FC = waterFC(soil, soilFunctions);
  int nlayers = Water_FC.size();
  
  //Anatomy parameters
  DataFrame paramsAnatomy = Rcpp::as<Rcpp::DataFrame>(x["paramsAnatomy"]);
  NumericVector SLA = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["SLA"]);
  // NumericVector Al2As = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["Al2As"]);
  // NumericVector WoodDensity = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["WoodDensity"]);

  
  //Allometric parameters
  DataFrame paramsAllometries = Rcpp::as<Rcpp::DataFrame>(x["paramsAllometries"]);
  NumericVector Hmax  = paramsAllometries["Hmax"];
  NumericVector Zmax  = paramsAllometries["Zmax"];
  NumericVector Aash  = paramsAllometries["Aash"];
  NumericVector Absh  = paramsAllometries["Absh"];
  NumericVector Bbsh  = paramsAllometries["Bbsh"];
  NumericVector r635  = paramsAllometries["r635"];
  NumericVector Acw  = paramsAllometries["Acw"];
  NumericVector Bcw  = paramsAllometries["Bcw"];
  NumericVector Acr  = paramsAllometries["Acr"];
  NumericVector B1cr  = paramsAllometries["B1cr"];
  NumericVector B2cr  = paramsAllometries["B2cr"];
  NumericVector B3cr  = paramsAllometries["B3cr"];
  NumericVector C1cr  = paramsAllometries["C1cr"];
  NumericVector C2cr  = paramsAllometries["C2cr"];
  NumericVector fHDmin= paramsAllometries["fHDmin"];
  NumericVector fHDmax= paramsAllometries["fHDmax"];
  
  
  
  //Detailed subday results
  List subdailyRes(numDays);
  
  //Water Output variables
  NumericMatrix PlantPsi(numDays, numCohorts);
  NumericMatrix PlantStress(numDays, numCohorts);
  NumericMatrix PlantTranspiration(numDays, numCohorts);
  NumericMatrix MaintenanceRespiration(numDays, numCohorts);
  NumericMatrix GrowthRespiration(numDays, numCohorts);
  NumericMatrix PlantSugarLeaf(numDays, numCohorts);
  NumericMatrix PlantStarchLeaf(numDays, numCohorts);
  NumericMatrix PlantSugarSapwood(numDays, numCohorts);
  NumericMatrix PlantStarchSapwood(numDays, numCohorts);
  NumericMatrix PlantSugarTransport(numDays, numCohorts);
  NumericMatrix PlantSA(numDays, numCohorts);
  NumericMatrix PlantSAgrowth(numDays, numCohorts);
  NumericMatrix GrossPhotosynthesis(numDays, numCohorts);
  NumericMatrix PlantLAIdead(numDays, numCohorts);
  NumericMatrix PlantLAIlive(numDays, numCohorts);
  NumericVector SAgrowthcum(numCohorts, 0.0);
  
  //Water balance output variables
  NumericVector GDD(numDays);
  NumericVector SoilEvaporation(numDays);
  NumericVector LAIcell(numDays);
  NumericVector LAIcelldead(numDays);
  NumericVector Cm(numDays);
  NumericVector LgroundPAR(numDays);
  NumericVector LgroundSWR(numDays);
  NumericVector Runoff(numDays);
  NumericVector NetRain(numDays);
  NumericVector Rain(numDays);
  NumericVector Snow(numDays);
  NumericVector Interception(numDays);
  NumericVector Infiltration(numDays);
  NumericVector DeepDrainage(numDays);
  NumericMatrix Eplantdays(numDays, nlayers);
  NumericMatrix Wdays(numDays, nlayers); //Soil moisture content in relation to field capacity
  NumericMatrix psidays(numDays, nlayers);
  NumericMatrix MLdays(numDays, nlayers);
  NumericVector WaterTable(numDays, NA_REAL);
  NumericVector MLTot(numDays, 0.0);
  NumericVector SWE(numDays, 0.0);
  

  
  //Count years (times structural variables will be updated)
  int numYears = 0;
  for(int i=0;i<numDays;i++) {
    if(((DOY[i]==1) & (i>0)) | ((i==(numDays-1)) & (DOY[i]>=365))) numYears = numYears + 1;
  }
  DataFrame standSummary = DataFrame::create(
    _["LeafAreaIndex"] = NumericVector(numYears+1,0.0),
    _["TreeBasalArea"] = NumericVector(numYears+1,0.0),
    _["TreeDensity"] = NumericVector(numYears+1,0.0),
    _["ShrubCover"] = NumericVector(numYears+1,0.0),
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
  recordStandSummary(standSummary, LAI_live, N, DBH, Cover, H, 0);
  
  NumericVector Wini = soil["W"];
  Wdays(0,_) = Wini;
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
    if(verbose & (i%10 == 0)) Rcout<<".";//<<i;
    
    double wind = WindSpeed[i];
    if(NumericVector::is_na(wind)) wind = control["defaultWindSpeed"]; //Default 1 m/s -> 10% of fall every day
    if(wind<0.1) wind = 0.1; //Minimum windspeed abovecanopy
    
    //If DOY == 1 reset PLC (Growth assumed)
    if(cavitationRefill=="annual") {
      if(DOY[i]==1) {
        DataFrame internalWater = Rcpp::as<Rcpp::DataFrame>(x["internalWater"]);
        if(transpirationMode=="Granier") {
          NumericVector PLC = Rcpp::as<Rcpp::NumericVector>(internalWater["PLC"]);
          for(int j=0;j<PLC.length();j++) PLC[j] = 0.0;
        } else {
          NumericVector StemPLC = Rcpp::as<Rcpp::NumericVector>(internalWater["PLCstem"]);
          for(int j=0;j<StemPLC.length();j++) StemPLC[j] = 0.0;
        }
      }
    }
    
    if(unlimitedSoilWater) {
      NumericVector W = soil["W"];
      for(int h=0;h<W.size();h++) W[h] = 1.0;
    }
    
    //1. Phenology and leaf fall
    if(leafPhenology) updateLeaves(x, DOY[i], MeanTemperature[i], wind);
    
    //Store GDD
    GDD[i] = canopyParams["gdd"];
    
    //2. Water balance and photosynthesis
    if(transpirationMode=="Granier") {
      double er = erFactor(DOY[i], PET[i], Precipitation[i]);
      s = spwbDay1(x, soil, MeanTemperature[i], PET[i], Precipitation[i], er, 0.0, 
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
      double rhmin = MinRelativeHumidity[i];
      double rhmax = MaxRelativeHumidity[i];
      double rad = Radiation[i];
      PET[i] = meteoland::penman(latrad, elevation, slorad, asprad, J, tmin, tmax, rhmin, rhmax, rad, wind);
      double er = erFactor(DOY[i], PET[i], Precipitation[i]);
      s = spwbDay2(x, soil, tmin, tmax, 
                   rhmin, rhmax, rad, wind, 
                   latitude, elevation, slope, aspect,
                   solarConstant, delta, Precipitation[i], PET[i], 
                   er, 0.0, verbose);
    }    
    
    List stand = s["Stand"];
    LgroundPAR[i] = stand["LgroundPAR"];
    LgroundSWR[i] = stand["LgroundSWR"];
    LAIcell[i] = stand["LAIcell"];
    LAIcelldead[i] = stand["LAIcelldead"];
    Cm[i] = stand["Cm"];
    
    List db = s["WaterBalance"];
    DeepDrainage[i] = db["DeepDrainage"];
    Infiltration[i] = db["Infiltration"];
    Runoff[i] = db["Runoff"];
    NetRain[i] = db["NetRain"];
    Rain[i] = db["Rain"];
    Snow[i] = db["Snow"];
    Interception[i] = Precipitation[i]-NetRain[i];
    
    List sb = s["Soil"];
    NumericVector psi = sb["psi"];
    SoilEvaporation[i] = sum(Rcpp::as<Rcpp::NumericVector>(sb["SoilEvaporation"]));
    psidays(i,_) = psi;
    NumericVector EplantVec = sb["PlantExtraction"];
    SWE[i] = soil["SWE"];
    
    List Plants = s["Plants"];
    NumericVector EplantCoh = Plants["Transpiration"];
    Eplantdays(i,_) = EplantVec;
    PlantTranspiration(i,_) = EplantCoh;
    PlantStress(i,_) = Rcpp::as<Rcpp::NumericVector>(Plants["DDS"]);
    NumericVector psiCoh;
    if(transpirationMode=="Granier") {
      psiCoh =  Rcpp::as<Rcpp::NumericVector>(Plants["psi"]);  
    } else {
      psiCoh =  clone(Rcpp::as<Rcpp::NumericVector>(internalWater["psiLeaf"])); 
      NumericVector Agd =  Rcpp::as<Rcpp::NumericVector>(Plants["GrossPhotosynthesis"]);
    }
    PlantPsi(i,_) = psiCoh;

    //3. Carbon balance and growth
    DataFrame growDay = growthDay(x,s,MeanTemperature[i]);  
    
    //4. Assemble output
    NumericVector MaintenanceRespirationDay = Rcpp::as<Rcpp::NumericVector>(growDay["MaintenanceRespiration"]);
    NumericVector GrowthRespirationDay = Rcpp::as<Rcpp::NumericVector>(growDay["GrowthRespiration"]);
    NumericVector GrossPhotosynthesisDay = Rcpp::as<Rcpp::NumericVector>(growDay["GrossPhotosynthesis"]);
    NumericVector PlantSugarLeafDay = Rcpp::as<Rcpp::NumericVector>(growDay["PlantSugarLeaf"]);
    NumericVector PlantStarchLeafDay = Rcpp::as<Rcpp::NumericVector>(growDay["PlantStarchLeaf"]);
    NumericVector PlantSugarSapwoodDay = Rcpp::as<Rcpp::NumericVector>(growDay["PlantSugarSapwood"]);
    NumericVector PlantStarchSapwoodDay = Rcpp::as<Rcpp::NumericVector>(growDay["PlantStarchSapwood"]);
    NumericVector PlantSugarTransportDay = Rcpp::as<Rcpp::NumericVector>(growDay["PlantSugarTransport"]);
    NumericVector PlantSADay = Rcpp::as<Rcpp::NumericVector>(growDay["PlantSA"]);
    NumericVector PlantLAIliveDay = Rcpp::as<Rcpp::NumericVector>(growDay["PlantLAIlive"]);
    NumericVector PlantLAIdeadDay = Rcpp::as<Rcpp::NumericVector>(growDay["PlantLAIdead"]);
    NumericVector PlantSAgrowthDay = Rcpp::as<Rcpp::NumericVector>(growDay["PlantSAgrowth"]);
    for(int j=0;j<numCohorts;j++){
      GrossPhotosynthesis(i,j) = GrossPhotosynthesisDay[j];
      MaintenanceRespiration(i,j) = MaintenanceRespirationDay[j];
      GrowthRespiration(i,j) = GrowthRespirationDay[j];
      PlantSugarLeaf(i,j) = PlantSugarLeafDay[j];
      PlantStarchLeaf(i,j) = PlantStarchLeafDay[j];
      PlantSugarSapwood(i,j) = PlantSugarSapwoodDay[j];
      PlantStarchSapwood(i,j) = PlantStarchSapwoodDay[j];
      PlantSugarTransport(i,j) = PlantSugarTransportDay[j];
      PlantSA(i,j) = PlantSADay[j];
      PlantLAIlive(i,j) = PlantLAIliveDay[j];
      PlantLAIdead(i,j) = PlantLAIdeadDay[j];
      PlantSAgrowth(i,j) = PlantSAgrowthDay[j];
      SAgrowthcum[j] += PlantSAgrowthDay[j]; //Store cumulative SA growth (for structural variable update)
    }
    
    //4 Update structural variables
    if(((DOY[i]==1) & (i>0)) | ((i==(numDays-1)) & (DOY[i]>=365))) { 
      if(verbose) Rcout<<" [update structural variables] ";
      iyear++;
      
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
        if(!NumericVector::is_na(DBH[j])) {
          double fHmod = std::max(0.0,std::min(1.0,(1.0-((H[j]-137.0)/(Hmax[j]-137.0)))));
          double fHD = (fHDmin[j]*(L[j]/100.0) + fHDmax[j]*(1.0-(L[j]/100.0)))*fHmod;
          // Rcout << fHmod<<" "<< fHD<<" "<< L[j]<<"\n";
          H[j] = H[j] + fHD*deltaDBH[j];
        }
      }
      NumericVector crNew = treeCrownRatioMED(N, DBH, H, Acw, Bcw, Acr, B1cr, B2cr, B3cr, C1cr, C2cr);
      for(int j=0;j<numCohorts; j++) {
        if(!NumericVector::is_na(DBH[j])) {
          CR[j] = crNew[j];
        }
      }

      //Shrub variables
      for(int j=0;j<numCohorts; j++) {
        if(NumericVector::is_na(DBH[j])) {
          double Wleaves = LAI_live[j]/((N[j]/10000)*SLA[j]);  //Calculates the biomass (dry weight) of leaves
          double PV = pow(Wleaves*r635[j]/Absh[j], 1.0/Bbsh[j]); //Calculates crown phytovolume (in m3)
          H[j] = pow(1e6*PV/(Aash[j]*CR[j]), 1.0/3.0); //Updates shrub height
          if(H[j]> Hmax[j]) { //Limit height (and update the former variables)
            H[j] = Hmax[j];
            PV = (Aash[j]*pow(H[j],2.0)/10000.0)*(H[j]/100.0)*CR[j];
            Wleaves = Absh[j]*pow(PV, Bbsh[j])/r635[j];
            double prevLive = LAI_live[j];
            LAI_live[j] = Wleaves*((N[j]/10000)*SLA[j]); //Update LAI_live to the maximum
            LAI_dead[j] += prevLive - LAI_live[j]; //Increment dead LAI with the difference
          }
          Cover[j] = (N[j]*Aash[j]*pow(H[j],2.0)/1e6); //Updates shrub cover
        }
      }
      
      // Store stand structure
      standStructures[iyear] = clone(above);
      recordStandSummary(standSummary, LAI_live, N, DBH, Cover, H,iyear);
    }
    if(i<(numDays-1))  Wdays(i+1,_) = as<Rcpp::NumericVector>(soil["W"]);
    WaterTable[i] = waterTableDepth(soil, soilFunctions);
  }
  if(verbose) Rcout << "done.\n";
  
  NumericVector Transpiration(numDays,0.0);
  for(int l=0;l<nlayers;l++) {
    MLdays(_,l) = Wdays(_,l)*Water_FC[l]; 
    MLTot = MLTot + MLdays(_,l);
    Transpiration = Transpiration + Eplantdays(_,l);
  }
  
  NumericVector Evapotranspiration = Transpiration+SoilEvaporation + Interception;
  

  DataFrame SWB = DataFrame::create(_["W"]=Wdays, _["ML"]=MLdays,_["MLTot"]=MLTot,
                                    _["WTD"] = WaterTable,
                                    _["SWE"] = SWE, _["PlantExt"]=Eplantdays, _["psi"]=psidays);
  Rcpp::DataFrame Stand = DataFrame::create(_["GDD"] = GDD,
                                          _["LAIcell"]=LAIcell, _["LAIcelldead"]=LAIcelldead,
                                          _["Cm"]=Cm, _["LgroundPAR"] = LgroundPAR, _["LgroundSWR"] = LgroundSWR);
  Rcpp::DataFrame DWB = DataFrame::create(_["PET"]=PET, 
                                          _["Precipitation"] = Precipitation, _["Rain"] = Rain, _["Snow"] = Snow, 
                                          _["NetRain"]=NetRain,_["Infiltration"]=Infiltration, _["Runoff"]=Runoff, _["DeepDrainage"]=DeepDrainage, 
                                          _["Evapotranspiration"]=Evapotranspiration,_["SoilEvaporation"]=SoilEvaporation,
                                          _["Transpiration"]=Transpiration);
  
  SWB.attr("row.names") = meteo.attr("row.names");
  DWB.attr("row.names") = meteo.attr("row.names");
  Stand.attr("row.names") = meteo.attr("row.names");
  
  PlantTranspiration.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names"));
  GrossPhotosynthesis.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names")) ;
  MaintenanceRespiration.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names")) ;
  GrowthRespiration.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names")) ;
  PlantSugarLeaf.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names")) ;
  PlantStarchLeaf.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names")) ;
  PlantSugarSapwood.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names")) ;
  PlantStarchSapwood.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names")) ;
  PlantSugarTransport.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names")) ;
  PlantSA.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names")) ;
  PlantSAgrowth.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names")) ;
  PlantPsi.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names")) ;
  PlantStress.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names")) ;
  PlantLAIdead.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names")) ;
  PlantLAIlive.attr("dimnames") = List::create(meteo.attr("row.names"), cohorts.attr("row.names")) ;
  
  NumericVector topo = NumericVector::create(elevation, slope, aspect);
  topo.attr("names") = CharacterVector::create("elevation", "slope", "aspect");
  
  List plantStorage = List::create(
    Named("SugarLeaf") = PlantSugarLeaf,
    Named("StarchLeaf") = PlantStarchLeaf,
    Named("SugarSapwood") = PlantSugarSapwood,
    Named("StarchSapwood") = PlantStarchSapwood,
    Named("SugarTransport") = PlantSugarTransport
  );
  List carbonBalance = List::create(
    Named("GrossPhotosynthesis") = GrossPhotosynthesis,
    Named("MaintenanceRespiration") = MaintenanceRespiration,
    Named("GrowthRespiration") = GrowthRespiration
  );
  List l = List::create(Named("latitude") = latitude,
                        Named("topography") = topo,
                        Named("spwbInput") = spwbInput,
                        Named("soilInput") = soilInput,
                        Named("WaterBalance")=DWB, 
                        Named("Soil")=SWB,
                        Named("Stand")=Stand,
                        Named("CarbonBalance") = carbonBalance,
                        Named("PlantStorage") = plantStorage,
                        Named("PlantTranspiration") = PlantTranspiration,
                        Named("PlantSAgrowth") = PlantSAgrowth,
                        Named("PlantSA")=PlantSA,
                        Named("PlantPsi") = PlantPsi, 
                        Named("PlantStress") = PlantStress,
                        Named("PlantLAIdead") = PlantLAIdead,
                        Named("PlantLAIlive") = PlantLAIlive,
                        Named("StandStructures") = standStructures,
                        Named("StandSummary") = standSummary);
  l.attr("class") = CharacterVector::create("growth","list");
  return(l);
}