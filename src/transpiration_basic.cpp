#define STRICT_R_HEADERS
#include <Rcpp.h>
#include <numeric>
#include "communication_structures.h"
#include "lightextinction_basic.h"
#include "windextinction.h"
#include "windKatul.h"
#include "hydraulics.h"
#include "biophysicsutils.h"
#include "phenology.h"
#include "forestutils.h"
#include "tissuemoisture.h"
#include "carbon.h"
#include "photosynthesis.h"
#include "root.h"
#include "soil.h"
#include "spwb.h"
#include <meteoland.h>
using namespace Rcpp;

struct ParamsVolume{
  double leafpi0;
  double leafeps;
  double leafaf;
  double stempi0;
  double stemeps;
  double stem_c;
  double stem_d;
  double stemaf;
  double Vsapwood;
  double Vleaf;
  double LAI;
  double LAIlive;
}; 

//Plant volume in l·m-2 ground = mm
double plantVol(double plantPsi, ParamsVolume pars) {
  
  double leafrwc = tissueRelativeWaterContent(plantPsi, pars.leafpi0, pars.leafeps, 
                                              plantPsi, pars.stem_c, pars.stem_d, 
                                              pars.leafaf);
  double stemrwc = tissueRelativeWaterContent(plantPsi, pars.stempi0, pars.stemeps, 
                                              plantPsi, pars.stem_c, pars.stem_d, 
                                              pars.stemaf);
  return((pars.Vleaf * leafrwc * pars.LAI) + (pars.Vsapwood * stemrwc * pars.LAIlive));
}


double findNewPlantPsiConnected(double flowFromRoots, double plantPsi, double rootCrownPsi,
                                ParamsVolume parsVol){
  if(plantPsi==rootCrownPsi) return(plantPsi);
  double V = plantVol(plantPsi, parsVol);
  //More negative rootCrownPsi causes increased flow due to water being removed
  double psiStep = rootCrownPsi - plantPsi;
  double Vnew = plantVol(plantPsi + psiStep, parsVol);
  while(std::abs(Vnew - V) > flowFromRoots) {
    psiStep = psiStep/2.0;
    Vnew = plantVol(plantPsi + psiStep, parsVol);
  }
  return(plantPsi+psiStep);
}


void transpirationBasic(List transpOutput, List x, NumericVector meteovec,  
                        double elevation, bool modifyInput = true) {
  
  //Will not modify input x 
  if(!modifyInput) {
    x = clone(x);
  }
  // Should have internal communication 
  NumericVector outputStand = as<NumericVector>(transpOutput["Stand"]);
  DataFrame outputPlants = as<DataFrame>(transpOutput["Plants"]);
  NumericMatrix outputExtraction = as<NumericMatrix>(transpOutput["Extraction"]);
  List outputExtractionPools = transpOutput["ExtractionPools"];

  //Control parameters
  List control = x["control"];
  String stemCavitationRecovery = control["stemCavitationRecovery"];
  String leafCavitationRecovery = control["leafCavitationRecovery"];
  double cavitationRecoveryMaximumRate = control["cavitationRecoveryMaximumRate"];
  String soilFunctions = control["soilFunctions"];
  double verticalLayerSize = control["verticalLayerSize"];
  String rhizosphereOverlap = control["rhizosphereOverlap"];
  double fullRhizosphereOverlapConductivity = 0.01; //For backwards compatibility
  if(control.containsElementNamed("fullRhizosphereOverlapConductivity")) fullRhizosphereOverlapConductivity = control["fullRhizosphereOverlapConductivity"];
  bool plantWaterPools = (rhizosphereOverlap!="total");
  double hydraulicRedistributionFraction = control["hydraulicRedistributionFraction"];
  String lfmcComponent = control["lfmcComponent"];

  //Soil water at field capacity
  DataFrame soil = Rcpp::as<Rcpp::DataFrame>(x["soil"]);
  NumericVector Water_FC = waterFC(soil, soilFunctions);
  
  //Meteo input
  double pet = meteovec["pet"];
  double rhmax = meteovec["rhmax"];
  double rhmin = meteovec["rhmin"];
  double tmax = meteovec["tmax"];
  double tmin = meteovec["tmin"];
  double Catm = meteovec["Catm"];
  double Patm = meteovec["Patm"];
  
  //Atmospheric pressure (if missing)
  if(NumericVector::is_na(Patm)) Patm = meteoland::utils_atmosphericPressure(elevation);
  
  //Daily average water vapor pressure at the atmosphere (kPa)
  double vpatm = meteoland::utils_averageDailyVP(tmin, tmax, rhmin, rhmax);
  double vpd = std::max(0.0, meteoland::utils_saturationVP((tmin+tmax)/2.0) - vpatm);
    
    
  // Canopy
  DataFrame canopyParams = Rcpp::as<Rcpp::DataFrame>(x["canopy"]);
  int ncanlayers = canopyParams.nrow();
  
  //Vegetation input
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  DataFrame cohorts = Rcpp::as<Rcpp::DataFrame>(x["cohorts"]);
  NumericVector LAIlive = Rcpp::as<Rcpp::NumericVector>(above["LAI_live"]);
  NumericVector LAIphe = Rcpp::as<Rcpp::NumericVector>(above["LAI_expanded"]);
  NumericVector LAIdead = Rcpp::as<Rcpp::NumericVector>(above["LAI_dead"]);
  NumericVector H = Rcpp::as<Rcpp::NumericVector>(above["H"]);
  NumericVector CR = Rcpp::as<Rcpp::NumericVector>(above["CR"]);
  int numCohorts = LAIphe.size();
  
  //Root distribution input
  DataFrame belowdf = Rcpp::as<Rcpp::DataFrame>(x["below"]);
  List belowLayers = Rcpp::as<Rcpp::List>(x["belowLayers"]);
  NumericMatrix V = Rcpp::as<Rcpp::NumericMatrix>(belowLayers["V"]);
  
  //Water pools
  NumericMatrix Wpool = Rcpp::as<Rcpp::NumericMatrix>(belowLayers["Wpool"]);
  List RHOP;
  NumericVector poolProportions(numCohorts);
  if(plantWaterPools) {
    RHOP = belowLayers["RHOP"];
    poolProportions = belowdf["poolProportions"];
  }
  
  
  
  //Phenology parameters
  DataFrame paramsPhenology = Rcpp::as<Rcpp::DataFrame>(x["paramsPhenology"]);
  CharacterVector phenoType = paramsPhenology["PhenologyType"];
  
  //Parameters  
  DataFrame paramsAnatomy = Rcpp::as<Rcpp::DataFrame>(x["paramsAnatomy"]);
  NumericVector Al2As = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["Al2As"]);
  NumericVector r635 = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["r635"]);
  
  DataFrame paramsInterception = Rcpp::as<Rcpp::DataFrame>(x["paramsInterception"]);
  NumericVector kPAR = Rcpp::as<Rcpp::NumericVector>(paramsInterception["kPAR"]);
  NumericVector kSWR = Rcpp::as<Rcpp::NumericVector>(paramsInterception["kSWR"]);
  
  DataFrame paramsTransp = Rcpp::as<Rcpp::DataFrame>(x["paramsTranspiration"]);
  NumericVector Gswmin = Rcpp::as<Rcpp::NumericVector>(paramsTransp["Gswmin"]);
  NumericVector Psi_Extract = Rcpp::as<Rcpp::NumericVector>(paramsTransp["Psi_Extract"]);
  NumericVector Exp_Extract = Rcpp::as<Rcpp::NumericVector>(paramsTransp["Exp_Extract"]);
  NumericVector VCstem_c = Rcpp::as<Rcpp::NumericVector>(paramsTransp["VCstem_c"]);
  NumericVector VCstem_d = Rcpp::as<Rcpp::NumericVector>(paramsTransp["VCstem_d"]);
  NumericVector VCleaf_c = Rcpp::as<Rcpp::NumericVector>(paramsTransp["VCleaf_c"]);
  NumericVector VCleaf_d = Rcpp::as<Rcpp::NumericVector>(paramsTransp["VCleaf_d"]);
  NumericVector WUE = Rcpp::as<Rcpp::NumericVector>(paramsTransp["WUE"]);
  NumericVector WUE_par(numCohorts, 0.3643);
  NumericVector WUE_co2(numCohorts, 0.002757);
  NumericVector WUE_vpd(numCohorts, -0.4636);
  NumericVector Tmax_LAI(numCohorts, 0.134);
  NumericVector Tmax_LAIsq(numCohorts, -0.006);
  if(paramsTransp.containsElementNamed("Tmax_LAI")) {
    Tmax_LAI = Rcpp::as<Rcpp::NumericVector>(paramsTransp["Tmax_LAI"]);
    Tmax_LAIsq = Rcpp::as<Rcpp::NumericVector>(paramsTransp["Tmax_LAIsq"]);
  }
  if(paramsTransp.containsElementNamed("WUE_par")) {
    WUE_par = Rcpp::as<Rcpp::NumericVector>(paramsTransp["WUE_par"]);
  }
  if(paramsTransp.containsElementNamed("WUE_decay")) { //For compatibility with previous versions (2.7.5)
    WUE_par = Rcpp::as<Rcpp::NumericVector>(paramsTransp["WUE_decay"]);
  }
  if(paramsTransp.containsElementNamed("WUE_co2")) {
    WUE_co2 = Rcpp::as<Rcpp::NumericVector>(paramsTransp["WUE_co2"]);
  }
  if(paramsTransp.containsElementNamed("WUE_vpd")) {
    WUE_vpd = Rcpp::as<Rcpp::NumericVector>(paramsTransp["WUE_vpd"]);
  }
  //Water storage parameters
  DataFrame paramsWaterStorage = Rcpp::as<Rcpp::DataFrame>(x["paramsWaterStorage"]);
  NumericVector maxFMC = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["maxFMC"]);
  NumericVector maxMCstem = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["maxMCstem"]);
  NumericVector maxMCleaf = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["maxMCleaf"]);
  NumericVector StemPI0 = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["StemPI0"]);
  NumericVector StemEPS = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["StemEPS"]);
  NumericVector StemAF = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["StemAF"]);
  NumericVector Vsapwood = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["Vsapwood"]); //l·m-2 = mm
  NumericVector LeafPI0 = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["LeafPI0"]);
  NumericVector LeafEPS = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["LeafEPS"]);
  NumericVector LeafAF = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["LeafAF"]);
  NumericVector Vleaf = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["Vleaf"]); //l·m-2 = mm
  
  //Communication vectors
  //Comunication with outside
  DataFrame internalPhenology = Rcpp::as<Rcpp::DataFrame>(x["internalPhenology"]);
  NumericVector phi = Rcpp::as<Rcpp::NumericVector>(internalPhenology["phi"]);
  DataFrame internalWater = Rcpp::as<Rcpp::DataFrame>(x["internalWater"]);
  NumericVector PlantPsi = Rcpp::as<Rcpp::NumericVector>(internalWater["PlantPsi"]);
  NumericVector StemPLC = Rcpp::as<Rcpp::NumericVector>(internalWater["StemPLC"]);
  NumericVector LeafPLC = Rcpp::as<Rcpp::NumericVector>(internalWater["LeafPLC"]);
  //LAI distribution
  List internalLAIDistribution = x["internalLAIDistribution"];
  NumericMatrix LAIme = internalLAIDistribution["expanded"];
  NumericMatrix LAImd = internalLAIDistribution["dead"];
  NumericVector PrevLAIexpanded = internalLAIDistribution["PrevLAIexpanded"];
  NumericVector PrevLAIdead = internalLAIDistribution["PrevLAIdead"];
  NumericVector PARcohort = internalLAIDistribution["PARcohort"];
  
  //Determine whether leaves are out (phenology) and the adjusted Leaf area
  double s = 0.0, LAIcell = 0.0, LAIcelllive = 0.0, LAIcellexpanded = 0.0,LAIcelldead = 0.0;
  for(int c=0;c<numCohorts;c++) {
    s += (kPAR[c]*(LAIphe[c]+LAIdead[c]));
    LAIcell += LAIphe[c]+LAIdead[c];
    LAIcelldead += LAIdead[c];
    LAIcellexpanded +=LAIphe[c];
    LAIcelllive += LAIlive[c];
  }
  

  if(numCohorts>0) {
    bool recalc_LAI = false;
    if(NumericVector::is_na(PrevLAIexpanded[0]) || NumericVector::is_na(PrevLAIdead[0])) {
      recalc_LAI = true; 
    } else{
      if(sum(abs(LAIphe - PrevLAIexpanded))>0.001) {
        recalc_LAI = true; 
      } else {
        if(sum(abs(LAIdead - PrevLAIdead))>0.001) recalc_LAI = true;
      }
    }
    if(recalc_LAI) {
      NumericVector z(ncanlayers+1,0.0);
      for(int i=1;i<=ncanlayers;i++) z[i] = z[i-1] + verticalLayerSize;
      for(int i=0; i<numCohorts;i++) {
        PARcohort[i] = availableLight(H[i]*(1.0-(1.0-CR[i])/2.0), H, LAIphe, LAIdead, kPAR, CR);
        PrevLAIexpanded[i] = LAIphe[i];
        PrevLAIdead[i] = LAIdead[i];
      }
      //Update LAI distribution if necessary
      updateLAIdistributionVectors(LAIme, z, LAIphe, H, CR);
      updateLAIdistributionVectors(LAImd, z, LAIdead, H, CR);
    }
  }
  NumericVector CohASWRF = cohortAbsorbedSWRFraction(LAIme, LAImd, kSWR);
  CohASWRF = pow(CohASWRF, 0.75);
  
  //Apply fractions to potential evapotranspiration
  //Maximum canopy transpiration
  //    Tmax = PET[i]*(-0.006*pow(LAIcell[i],2.0)+0.134*LAIcell[i]+0.036); //From Granier (1999)
  NumericVector Tmax = pet*(Tmax_LAIsq*(LAIcell*LAIcell)+ Tmax_LAI*LAIcell); //From Granier (1999)
  
  //Fraction of Tmax attributed to each plant cohort
  double pabs = std::accumulate(CohASWRF.begin(),CohASWRF.end(),0.0);
  NumericVector TmaxCoh(numCohorts,0.0);
  if(pabs>0.0) TmaxCoh = Tmax*(CohASWRF/pabs);
  
  //Actual plant transpiration
  int nlayers = Wpool.ncol();

  NumericVector Eplant(numCohorts, 0.0), Agplant(numCohorts, 0.0);
  NumericVector DDS(numCohorts, 0.0), LFMC(numCohorts, 0.0);
  NumericVector PLCm(numCohorts), RWCsm(numCohorts), RWClm(numCohorts),RWCssm(numCohorts), RWClsm(numCohorts);
  NumericVector PWB(numCohorts,0.0);
  NumericVector Kl, epc, Vl;
  
  //Calculate unsaturated conductivity (mmolH20·m-1·s-1·MPa-1)
  NumericVector Kunsat = conductivity(soil, soilFunctions, true);
  //Calculate soil water potential
  NumericVector psiSoil = psi(soil,soilFunctions);

  
  NumericMatrix WaterM(numCohorts, nlayers);
  NumericMatrix KunsatM(numCohorts, nlayers);
  NumericMatrix psiSoilM(numCohorts, nlayers);
  if(plantWaterPools) {
    DataFrame soil_pool = clone(soil);
    NumericVector Ws_pool = soil_pool["W"];
    for(int j = 0; j<numCohorts;j++) {
      //Copy values of soil moisture from pool of cohort j
      for(int l = 0; l<nlayers;l++) Ws_pool[l] = Wpool(j,l);
      //Calculate unsaturated conductivity (mmolH20·m-1·s-1·MPa-1)
      KunsatM(j,_) = conductivity(soil_pool, soilFunctions, true);
      //Calculate soil water potential
      psiSoilM(j,_) = psi(soil_pool, soilFunctions);
      // Rcout<< Ws_pool[3]<< " "<< psiSoilM(j,3)<<"\n";
      //Calculate available water
      WaterM(j,_) = water(soil_pool, soilFunctions);
    }
  }
  // for(int i= 0;i<psiSoil.size();i++) Rcout<< "S "<<i<<" "<<psiSoil[i]<<"\n";
  
  ParamsVolume parsVol;
  
  for(int c=0;c<numCohorts;c++) {
    
    parsVol.stem_c = VCstem_c[c];
    parsVol.stem_d = VCstem_d[c];
    parsVol.leafpi0 = LeafPI0[c];
    parsVol.leafeps = LeafEPS[c];
    parsVol.leafaf = LeafAF[c];
    parsVol.stempi0 = StemPI0[c];
    parsVol.stemeps = StemEPS[c];
    parsVol.stemaf = StemAF[c];
    parsVol.Vsapwood = Vsapwood[c];
    parsVol.Vleaf = Vleaf[c];
    parsVol.LAI = LAIphe[c];
    parsVol.LAIlive = LAIlive[c];
    
    double rootCrownPsi = NA_REAL;
    
    //Cuticular transpiration    
    double lvp_tmax = leafVapourPressure(tmax,  PlantPsi[c]);
    double lvp_tmin = leafVapourPressure(tmin,  PlantPsi[c]);
    double lvpd_tmax = std::max(0.0, lvp_tmax - vpatm);
    double lvpd_tmin = std::max(0.0, lvp_tmin - vpatm);
    double E_gmin = Gswmin[c]*(lvpd_tmin+lvpd_tmax)/(2.0*Patm); // mol·s-1·m-2
    double E_gmin_day = E_gmin*LAIphe[c]*(24.0*3600.0*0.018); //L·d-1·m-2 = mm·d-1
    
    //Extraction from soil (can later be modified if there are changes in plant water content)
    if(!plantWaterPools) {
      NumericVector Klc(nlayers);
      NumericVector Kunlc(nlayers);
      for(int l=0;l<nlayers;l++) {
        Klc[l] = Psi2K(psiSoil[l], Psi_Extract[c], Exp_Extract[c]);
        //Limit Mean Kl due to previous cavitation
        if(stemCavitationRecovery!="total") {
          Klc[l] = std::min(Klc[l], 1.0-StemPLC[c]); 
        }
        Kunlc[l] = std::sqrt(Kunsat[l])*V(c,l);
      }
      double sumKunlc = sum(Kunlc);
      double Klcmean = sum(Klc*V(c,_));
      for(int l=0;l<nlayers;l++) {
        outputExtraction(c,l) = std::max(TmaxCoh[c]*Klcmean, E_gmin_day)*(Kunlc[l]/sumKunlc);
      }
      rootCrownPsi = averagePsi(psiSoil, V(c,_), Exp_Extract[c], Psi_Extract[c]);
    } else {
      NumericMatrix ExtractionPoolsCoh = outputExtractionPools[c]; //this is used to store extraction of a SINGLE plant cohort from all pools
      
      NumericMatrix RHOPcoh = Rcpp::as<Rcpp::NumericMatrix>(RHOP[c]);
      NumericMatrix RHOPcohDyn(numCohorts, nlayers);
      for(int l=0;l<nlayers;l++) {
        RHOPcohDyn(c,l) = RHOPcoh(c,l);
        for(int j=0; j<numCohorts;j++) {
          if(j!=c) {
            double overlapFactor = std::min(1.0, KunsatM(j,l)/(cmdTOmmolm2sMPa*fullRhizosphereOverlapConductivity));
            RHOPcohDyn(j,l) = RHOPcoh(j,l)*overlapFactor;
            RHOPcohDyn(c,l) = RHOPcohDyn(c,l) + (RHOPcoh(j,l) - RHOPcohDyn(j,l));
          } 
        }
      }
      NumericMatrix Klc(numCohorts, nlayers);
      NumericMatrix Kunlc(numCohorts, nlayers);
      NumericMatrix RHOPcohV(numCohorts, nlayers);
      for(int j = 0;j<numCohorts;j++) {
        for(int l=0;l<nlayers;l++) {
          RHOPcohV(j,l) = RHOPcohDyn(j,l)*V(c,l);
          Klc(j,l) = Psi2K(psiSoilM(c,l), Psi_Extract[c], Exp_Extract[c]);
          //Limit Mean Kl due to previous cavitation
          if(stemCavitationRecovery!="total") Klc(j,l) = std::min(Klc(j,l), 1.0-StemPLC[c]); 
          Kunlc(j,l) = std::sqrt(KunsatM(j,l))*RHOPcohV(j,l);
        }
      }
      double sumKunlc = sum(Kunlc);
      double Klcmean = sum(Klc*RHOPcohV);
      for(int l=0;l<nlayers;l++) {
        for(int j = 0;j<numCohorts;j++) {
          ExtractionPoolsCoh(j,l) = std::max(TmaxCoh[c]*Klcmean, E_gmin_day)*(Kunlc(j,l)/sumKunlc);
        }
        outputExtraction(c,l) = sum(ExtractionPoolsCoh(_,l)); // Sum extraction from all pools (layer l)
      }
      
      rootCrownPsi = averagePsiPool(psiSoilM, RHOPcohV, Exp_Extract[c], Psi_Extract[c]);
      // Rcout<< c << " : "<< psiSoilM(c,0) << " " << psiSoilM(c,1) << " " << psiSoilM(c,2) << " " << psiSoilM(c,3) << " " << rootCrownPsi<<"\n";
      // Rcout<< c << " : "<< RHOPcohV(c,0) << " " << RHOPcohV(c,1) << " " << RHOPcohV(c,2) << " " << RHOPcohV(c,3) << " " << rootCrownPsi<<"\n";
    }


    double oldVol = plantVol(PlantPsi[c], parsVol); 
    
    //Transpiration is the maximum of predicted extraction and cuticular transpiration
    double ext_sum = sum(outputExtraction(c,_));
    Eplant[c] = ext_sum;
    // PlantPsi[c] = findNewPlantPsiConnected(Eplant[c], PlantPsi[c], rootCrownPsi, parsVol);
    //For deciduous species, make water potential follow soil during winter
    // if(LAIphe[c]==0.0) PlantPsi[c] = rootCrownPsi;
    PlantPsi[c] = rootCrownPsi;
    // PlantPsi[c] = rootCrownPsi;
    double newVol = plantVol(PlantPsi[c], parsVol);
    
    double volDiff = newVol - oldVol;
    //Plant transpiration and water balance
    PWB[c] = volDiff;
    
    //Photosynthesis
    double fpar = std::min(1.0, pow(PARcohort[c]/100.0,WUE_par[c]));
    double fco2 = (1.0 - exp((-1.0)*WUE_co2[c]*Catm));
    double fvpd = pow(vpd, WUE_vpd[c]);
    if(vpd < 0.25) fvpd = 2.5 - (2.5 - pow(0.25, WUE_vpd[c]))*(vpd/0.25);
    // Rcout<<fpar<<" "<< fco2 << " "<< fvpd<< " "<< WUE[c]*fpar*fco2*fvpd<<"\n";
    Agplant[c] = WUE[c]*Eplant[c]*fpar*fco2*fvpd;
  }
  
  //Plant water status (StemPLC, RWC, DDS)
  for(int c=0;c<numCohorts;c++) {
    if(stemCavitationRecovery!="total") {
      StemPLC[c] = std::max(1.0 - xylemConductance(PlantPsi[c], 1.0, VCstem_c[c], VCstem_d[c]), StemPLC[c]); //Track current embolism if no refill
    } else {
      StemPLC[c] = 1.0 - xylemConductance(PlantPsi[c], 1.0, VCstem_c[c], VCstem_d[c]);
    }
    if(leafCavitationRecovery!="total") {
      LeafPLC[c] = std::max(1.0 - xylemConductance(PlantPsi[c], 1.0, VCleaf_c[c], VCleaf_d[c]), LeafPLC[c]); //Track current embolism if no refill
    } else {
      LeafPLC[c] = 1.0 - xylemConductance(PlantPsi[c], 1.0, VCleaf_c[c], VCleaf_d[c]);
    }
    
    //Relative water content and fuel moisture from plant water potential
    RWClm[c] =  tissueRelativeWaterContent(PlantPsi[c], LeafPI0[c], LeafEPS[c], 
                                           PlantPsi[c], VCstem_c[c], VCstem_d[c], 
                                           LeafAF[c]);
    RWCsm[c] =  tissueRelativeWaterContent(PlantPsi[c], StemPI0[c], StemEPS[c], 
                                           PlantPsi[c], VCstem_c[c], VCstem_d[c], 
                                           StemAF[c]);
    // The fraction of leaves will decrease due to phenology or processes leading to defoliation
    double fleaf = (1.0/r635[c])*(LAIphe[c]/LAIlive[c]);
    if(lfmcComponent=="fine") { //fine fuel moisture
      LFMC[c] = maxMCleaf[c]*RWClm[c]*fleaf + maxMCstem[c]*RWCsm[c]*(1.0 - fleaf);
    } else { //"leaf"
      LFMC[c] = maxFMC[c]*RWClm[c];
    }
    
    //Daily drought stress from plant WP
    DDS[c] = (1.0 - Psi2K(PlantPsi[c],Psi_Extract[c],Exp_Extract[c])); 
    if(phenoType[c] == "winter-deciduous" || phenoType[c] == "winter-semideciduous") DDS[c] = phi[c]*DDS[c];
      
    double SAmax = 10e4/Al2As[c]; //cm2·m-2 of leaf area
    double r = cavitationRecoveryMaximumRate*std::max(0.0, (PlantPsi[c] + 1.5)/1.5);
    if(stemCavitationRecovery=="rate") {
      StemPLC[c] = std::max(0.0, StemPLC[c] - (r/SAmax));
    }
    if(leafCavitationRecovery=="rate") {
      LeafPLC[c] = std::max(0.0, LeafPLC[c] - (r/SAmax));
    }
  }
  
  //Atempt to implement hydraulic redistribution
  if(hydraulicRedistributionFraction > 0.0) {
    if(!plantWaterPools) {
      for(int c=0;c<numCohorts;c++) {
        double redAmount = Eplant[c]*hydraulicRedistributionFraction;
        // Rcout<<c<< "red amount"<< redAmount;
        NumericVector Ws = soil["W"];
        NumericVector SW = water(soil, soilFunctions);
        double soilRWC = sum(Ws*SW)/sum(SW);
        NumericVector WDiff = Ws - soilRWC; 
        NumericVector DonorDiff = pmax(0.0, WDiff);
        NumericVector ReceiverDiff = pmax(0.0, -WDiff);
        NumericVector HD(nlayers,0.0);
        if(sum(DonorDiff)>0.0) {
          for(int l=0;l<nlayers;l++) {
            if(WDiff[l]>0.0) {
              HD[l] = redAmount*DonorDiff[l]/sum(DonorDiff);
            } else{
              HD[l] = -redAmount*ReceiverDiff[l]/sum(ReceiverDiff);
            }
            // Rcout<<" "<<l<<" "<<HD[l];
            outputExtraction(c,l) += HD[l];
          }
        }
        // Rcout<< "\n";
      }
    } else {
      for(int c=0;c<numCohorts;c++) {
        NumericMatrix ExtractionPoolsCoh = outputExtractionPools[c];
        for(int j=0;j<numCohorts;j++) {
          double redAmount = sum(ExtractionPoolsCoh(j,_))*hydraulicRedistributionFraction;
          NumericVector Ws = Wpool(j,_);
          NumericVector SW = WaterM(j,_);
          double soilRWC = sum(Ws*SW)/sum(SW);
          NumericVector WDiff = Ws - soilRWC; 
          NumericVector DonorDiff = pmax(0.0, WDiff);
          NumericVector ReceiverDiff = pmax(0.0, -WDiff);
          NumericVector HD(nlayers,0.0);
          if(sum(DonorDiff)>0.0) {
            for(int l=0;l<nlayers;l++) {
              if(WDiff[l]>0.0) {
                HD[l] = redAmount*DonorDiff[l]/sum(DonorDiff);
              } else{
                HD[l] = -redAmount*ReceiverDiff[l]/sum(ReceiverDiff);
              }
              outputExtraction(c,l) += HD[l];
              ExtractionPoolsCoh(j,l) += HD[l];
            }
          }
        }
      }
    }
  }
  

  // Copy output stand
  outputStand["LAI"] = LAIcell;
  outputStand["LAIlive"] = LAIcelllive;
  outputStand["LAIexpanded"] = LAIcellexpanded;
  outputStand["LAIdead"] = LAIcelldead;

  // Copy output plants
  NumericVector outputLAI = outputPlants["LAI"];
  NumericVector outputLAIlive = outputPlants["LAIlive"];
  NumericVector outputFPAR = outputPlants["FPAR"];
  NumericVector outputAbsorbedSWRFraction = outputPlants["AbsorbedSWRFraction"];
  NumericVector outputExtractionByPlants = outputPlants["Extraction"];
  NumericVector outputTranspirationByPlants = outputPlants["Transpiration"];
  NumericVector outputGrossPhotosynthesis = outputPlants["GrossPhotosynthesis"];
  NumericVector outputPlantPsi = outputPlants["PlantPsi"];
  NumericVector outputDDS = outputPlants["DDS"];
  NumericVector outputStemRWC = outputPlants["StemRWC"];
  NumericVector outputLeafRWC = outputPlants["LeafRWC"];
  NumericVector outputLFMC = outputPlants["LFMC"];
  NumericVector outputStemPLC = outputPlants["StemPLC"];
  NumericVector outputLeafPLC = outputPlants["LeafPLC"];
  NumericVector outputWaterBalance = outputPlants["WaterBalance"];
  for(int c =0;c<numCohorts;c++) {
    outputLAI[c] = LAIphe[c];
    outputLAIlive[c] = LAIlive[c];
    outputFPAR[c] = PARcohort[c];
    outputAbsorbedSWRFraction[c] = 100.0*CohASWRF[c];
    outputExtractionByPlants[c] = sum(outputExtraction(c,_));
    outputTranspirationByPlants[c] = Eplant[c];
    outputGrossPhotosynthesis[c] = Agplant[c];
    outputPlantPsi[c] = PlantPsi[c];
    outputDDS[c] = DDS[c];
    outputStemRWC[c] = RWCsm[c];
    outputLeafRWC[c] = RWClm[c];
    outputLFMC[c] = LFMC[c];
    outputStemPLC[c] = StemPLC[c];
    outputLeafPLC[c] = LeafPLC[c];
    outputWaterBalance[c] = PWB[c];
  }

}

//' Transpiration modes
//' 
//' High-level sub-models representing transpiration, plant hydraulics, photosynthesis and water relations 
//' within plants. 
//' 
//' Three sub-models are available: 
//' \itemize{
//'   \item{Sub-model in function \code{transp_transpirationGranier} was described in De \enc{Cáceres}{Caceres} et al. (2015), 
//'   and implements an approach originally described in Granier et al. (1999).} 
//'   \item{Sub-model in function \code{transp_transpirationSperry} was described in De \enc{Cáceres}{Caceres} et al. (2021), and
//'   implements a modelling approach originally described in Sperry et al. (2017).} 
//'   \item{Sub-model in function \code{transp_transpirationSureau} was described for SurEau-Ecos v2.0 model in Ruffault et al. (2022).} 
//' }
//' 
//' @param x An object of class \code{\link{spwbInput}} or \code{\link{growthInput}}, built using the 'Granier', 'Sperry' or 'Sureau' transpiration modes.
//' @param meteo A data frame with daily meteorological data series (see \code{\link{spwb}}).
//' @param day An integer to identify a day (row) within the \code{meteo} data frame.
//' @param latitude Latitude (in degrees).
//' @param elevation,slope,aspect Elevation above sea level (in m), slope (in degrees) and aspect (in degrees from North).
//' @param modifyInput Boolean flag to indicate that the input \code{x} object is allowed to be modified during the simulation.
//' 
//' @return
//' A list with the following elements:
//' \itemize{
//'   \item{\code{"cohorts"}: A data frame with cohort information, copied from \code{\link{spwbInput}}.}
//'   \item{\code{"Stand"}: A vector of stand-level variables.}
//'   \item{\code{"Plants"}: A data frame of results for each plant cohort. When using \code{transp_transpirationGranier}, element \code{"Plants"} includes:
//'     \itemize{
//'       \item{\code{"LAI"}: Leaf area index of the plant cohort.}
//'       \item{\code{"LAIlive"}: Leaf area index of the plant cohort, assuming all leaves are unfolded.}
//'       \item{\code{"AbsorbedSWRFraction"}: Fraction of SWR absorbed by each cohort.}
//'       \item{\code{"Transpiration"}: Transpirated water (in mm) corresponding to each cohort.}
//'       \item{\code{"GrossPhotosynthesis"}: Gross photosynthesis (in gC/m2) corresponding to each cohort.}
//'       \item{\code{"psi"}: Water potential (in MPa) of the plant cohort (average over soil layers).}
//'       \item{\code{"DDS"}: Daily drought stress \[0-1\] (relative whole-plant conductance).}
//'     }
//'   When using \code{transp_transpirationSperry} or \code{transp_transpirationSureau}, element \code{"Plants"} includes:
//'     \itemize{
//'       \item{\code{"LAI"}: Leaf area index of the plant cohort.}
//'       \item{\code{"LAIlive"}: Leaf area index of the plant cohort, assuming all leaves are unfolded.}
//'       \item{\code{"Extraction"}: Water extracted from the soil (in mm) for each cohort.}
//'       \item{\code{"Transpiration"}: Transpirated water (in mm) corresponding to each cohort.}
//'       \item{\code{"GrossPhotosynthesis"}: Gross photosynthesis (in gC/m2) corresponding to each cohort.}
//'       \item{\code{"NetPhotosynthesis"}: Net photosynthesis (in gC/m2) corresponding to each cohort.}
//'       \item{\code{"RootPsi"}: Minimum water potential (in MPa) at the root collar.}
//'       \item{\code{"StemPsi"}: Minimum water potential (in MPa) at the stem.}
//'       \item{\code{"StemPLC"}: Proportion of conductance loss in stem.}
//'       \item{\code{"LeafPsiMin"}: Minimum (predawn) water potential (in MPa) at the leaf (representing an average leaf).}
//'       \item{\code{"LeafPsiMax"}: Maximum (midday) water potential (in MPa) at the leaf (representing an average leaf).}
//'       \item{\code{"LeafPsiMin_SL"}: Minimum (predawn) water potential (in MPa) at sunlit leaves.}
//'       \item{\code{"LeafPsiMax_SL"}: Maximum (midday) water potential (in MPa) at sunlit leaves.}
//'       \item{\code{"LeafPsiMin_SH"}: Minimum (predawn) water potential (in MPa) at shade leaves.}
//'       \item{\code{"LeafPsiMax_SH"}: Maximum (midday) water potential (in MPa) at shade leaves.}
//'       \item{\code{"dEdP"}: Overall soil-plant conductance (derivative of the supply function).}
//'       \item{\code{"DDS"}: Daily drought stress \[0-1\] (relative whole-plant conductance).}
//'       \item{\code{"StemRWC"}: Relative water content of stem tissue (including symplasm and apoplasm).}
//'       \item{\code{"LeafRWC"}: Relative water content of leaf tissue (including symplasm and apoplasm).}
//'       \item{\code{"LFMC"}: Live fuel moisture content (in percent of dry weight).}
//'       \item{\code{"WaterBalance"}: Plant water balance (extraction - transpiration).}
//'     }
//'   }
//'   \item{\code{"Extraction"}: A data frame with mm of water extracted from each soil layer (in columns) by each cohort (in rows). The sum of a given row is equal to the total extraction of the corresponding plant cohort.}
//'   \item{\code{"ExtractionPools"}: A named list with as many elements as plant cohorts, where each element is a matrix data with mm of water extracted from each layer (in columns) of the water pool of each cohort (in rows). The sum of a given matrix is equal to the total extraction of the corresponding plant cohort.}
//' 
//'   The remaining items are only given by \code{transp_transpirationSperry} or \code{transp_transpirationSureau}:
//'   \item{\code{"EnergyBalance"}: A list with the following elements:
//'     \itemize{
//'       \item{\code{"Temperature"}: A data frame with the temperature of the atmosphere ('Tatm'), canopy ('Tcan') and soil ('Tsoil.1', 'Tsoil.2', ...) for each time step.}
//'       \item{\code{"CanopyEnergyBalance"}: A data frame with the components of the canopy energy balance (in W/m2) for each time step.}
//'       \item{\code{"SoilEnergyBalance"}: A data frame with the components of the soil energy balance (in W/m2) for each time step.}
//'     }  
//'   }
//'   \item{\code{"RhizoPsi"}: Minimum water potential (in MPa) inside roots, after crossing rhizosphere, per cohort and soil layer.}
//'   \item{\code{"Sunlitleaves"} and \code{"ShadeLeaves"}: Data frames for sunlit leaves and shade leaves and the following columns per cohort:
//'     \itemize{
//'       \item{\code{"LAI"}: Cumulative leaf area index of sunlit/shade leaves.}
//'       \item{\code{"Vmax298"}: Average maximum carboxilation rate for sunlit/shade leaves.}
//'       \item{\code{"Jmax298"}: Average maximum electron transport rate for sunlit/shade leaves.}
//'     }  
//'   }
//'   \item{\code{"ExtractionInst"}: Water extracted by each plant cohort during each time step.}
//'   \item{\code{"PlantsInst"}: A list with instantaneous (per time step) results for each plant cohort:
//'     \itemize{
//'       \item{\code{"E"}: A data frame with the cumulative transpiration (mm) for each plant cohort during each time step. }
//'       \item{\code{"Ag"}: A data frame with the cumulative gross photosynthesis (gC/m2) for each plant cohort during each time step. }
//'       \item{\code{"An"}: A data frame with the cumulative net photosynthesis (gC/m2) for each plant cohort during each time step. }
//'       \item{\code{"Sunlitleaves"} and \code{"ShadeLeaves"}: Lists with instantaneous (for each time step) results for sunlit leaves and shade leaves and the following items:
//'         \itemize{
//'           \item{\code{"Abs_SWR"}: A data frame with instantaneous absorbed short-wave radiation (SWR).} 
//'           \item{\code{"Net_LWR"}: A data frame with instantaneous net long-wave radiation (LWR).} 
//'           \item{\code{"An"}: A data frame with instantaneous net photosynthesis (in micromol/m2/s). }
//'           \item{\code{"Ci"}: A data frame with instantaneous intercellular CO2 concentration (in ppm). }
//'           \item{\code{"GW"}: A data frame with instantaneous stomatal conductance (in mol/m2/s). }
//'           \item{\code{"VPD"}: A data frame with instantaneous vapour pressure deficit (in kPa). }
//'           \item{\code{"Temp"}: A data frame with leaf temperature (in degrees Celsius). }
//'           \item{\code{"Psi"}: A data frame with leaf water potential (in MPa). }
//'         }
//'       }
//'       \item{\code{"dEdP"}: A data frame with the slope of the plant supply function (an estimation of whole-plant conductance).}
//'       \item{\code{"RootPsi"}: A data frame with root crown water potential (in MPa) for each plant cohort during each time step.}
//'       \item{\code{"StemPsi"}: A data frame with stem water potential (in MPa) for each plant cohort during each time step.}
//'       \item{\code{"LeafPsi"}: A data frame with leaf (average) water potential (in MPa) for each plant cohort during each time step. }
//'       \item{\code{"StemPLC"}: A data frame with the proportion loss of conductance \[0-1\] for each plant cohort during each time step. }
//'       \item{\code{"StemRWC"}: A data frame with the (average) relative water content of stem tissue \[0-1\] for each plant cohort during each time step. }
//'       \item{\code{"LeafRWC"}: A data frame with the relative water content of leaf tissue \[0-1\] for each plant cohort during each time step. }
//'       \item{\code{"StemSympRWC"}: A data frame with the (average) relative water content of symplastic stem tissue \[0-1\] for each plant cohort during each time step. }
//'       \item{\code{"LeafSympRWC"}: A data frame with the relative water content of symplastic leaf tissue \[0-1\] for each plant cohort during each time step. }
//'       \item{\code{"PWB"}: A data frame with plant water balance (extraction - transpiration).}
//'     }
//'   }
//'   \item{\code{"LightExtinction"}: A list of information regarding radiation balance through the canopy, as returned by function \code{\link{light_instantaneousLightExtinctionAbsortion}}.}
//'   \item{\code{"CanopyTurbulence"}: Canopy turbulence (see \code{\link{wind_canopyTurbulence}}).}
//'   \item{\code{"SupplyFunctions"}: If \code{stepFunctions} is not missing, a list of supply functions, photosynthesis functions and profit maximization functions.}
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
//' @seealso \code{\link{spwb_day}}, \code{\link{plot.spwb_day}}
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
//' #Define soil with default soil params (4 layers)
//' examplesoil <- defaultSoilParams(4)
//' 
//' #Initialize control parameters
//' control <- defaultControl("Granier")
//' 
//' #Initialize input
//' x1 <- spwbInput(exampleforest,examplesoil, SpParamsMED, control)
//' 
//' # Transpiration according to Granier's model, plant water potential 
//' # and plant stress for a given day
//' t1 <- transp_transpirationGranier(x1, examplemeteo, 1, 
//'                                  latitude = 41.82592, elevation = 100, slope = 0, aspect = 0, 
//'                                  modifyInput = FALSE)
//' 
//' #Switch to 'Sperry' transpiration mode
//' control <- defaultControl("Sperry")
//' 
//' #Initialize input
//' x2 <- spwbInput(exampleforest,examplesoil, SpParamsMED, control)
//' 
//' # Transpiration according to Sperry's model
//' t2 <- transp_transpirationSperry(x2, examplemeteo, 1, 
//'                                 latitude = 41.82592, elevation = 100, slope = 0, aspect = 0,
//'                                 modifyInput = FALSE)
//'                                 
//' #Switch to 'Sureau' transpiration mode
//' control <- defaultControl("Sureau")
//' 
//' #Initialize input
//' x3 <- spwbInput(exampleforest,examplesoil, SpParamsMED, control)
//' 
//' # Transpiration according to Sureau model
//' t3 <- transp_transpirationSureau(x3, examplemeteo, 1, 
//'                                   latitude = 41.82592, elevation = 100, slope = 0, aspect = 0,
//'                                   modifyInput = FALSE)
//'                                 
//' @name transp_modes
//' @keywords internal
// [[Rcpp::export("transp_transpirationGranier")]]
List transpirationGranier(List x, DataFrame meteo, int day,
                          double latitude, double elevation, double slope, double aspect, 
                          bool modifyInput = true) {


  List control = x["control"];
  if(!meteo.containsElementNamed("MinTemperature")) stop("Please include variable 'MinTemperature' in weather input.");
  NumericVector MinTemperature = meteo["MinTemperature"];
  if(!meteo.containsElementNamed("MaxTemperature")) stop("Please include variable 'MaxTemperature' in weather input.");
  NumericVector MaxTemperature = meteo["MaxTemperature"];
  if(!meteo.containsElementNamed("MinRelativeHumidity")) stop("Please include variable 'MinRelativeHumidity' in weather input.");
  NumericVector MinRelativeHumidity = meteo["MinRelativeHumidity"];
  if(!meteo.containsElementNamed("MaxRelativeHumidity")) stop("Please include variable 'MaxRelativeHumidity' in weather input.");
  NumericVector MaxRelativeHumidity = meteo["MaxRelativeHumidity"];
  if(!meteo.containsElementNamed("Radiation")) stop("Please include variable 'Radiation' in weather input.");
  NumericVector Radiation = meteo["Radiation"];
  NumericVector WindSpeed(MinTemperature.length(), NA_REAL);
  if(meteo.containsElementNamed("WindSpeed")) WindSpeed = meteo["WindSpeed"];
  NumericVector CO2(MinTemperature.length(), NA_REAL);
  if(meteo.containsElementNamed("CO2")) CO2 = meteo["CO2"];
  NumericVector Patm(MinTemperature.length(), NA_REAL);
  if(meteo.containsElementNamed("Patm")) Patm = meteo["Patm"];
  
  if(NumericVector::is_na(latitude)) stop("Value for 'latitude' should not be missing.");
  double latrad = latitude * (M_PI/180.0);
  if(NumericVector::is_na(aspect)) aspect = 0.0;
  if(NumericVector::is_na(slope)) slope = 0.0;
  double asprad = aspect * (M_PI/180.0);
  double slorad = slope * (M_PI/180.0);
  
  CharacterVector dateStrings = getWeatherDates(meteo);
  std::string c = as<std::string>(dateStrings[day-1]);
  int J = meteoland::radiation_julianDay(std::atoi(c.substr(0, 4).c_str()),std::atoi(c.substr(5,2).c_str()),std::atoi(c.substr(8,2).c_str()));

  double tmin = MinTemperature[day-1];
  double tmax = MaxTemperature[day-1];
  double tday = meteoland::utils_averageDaylightTemperature(tmin, tmax);
  double rhmax = MaxRelativeHumidity[day-1];
  double rhmin = MinRelativeHumidity[day-1];
  double rad = Radiation[day-1];
  double wind = WindSpeed[day-1];
  double Catm = CO2[day-1];

  double pet = meteoland::penman(latrad, elevation, slorad, asprad, J, 
                             tmin, tmax, rhmin, rhmax, rad, wind);
  
  if(NumericVector::is_na(Catm)) Catm = control["defaultCO2"];
  NumericVector meteovec = NumericVector::create(
    Named("tmax") = tmax,
    Named("tmin") = tmin,
    Named("rhmin") = rhmin, 
    Named("rhmax") = rhmax,
    Named("tday") = tday, 
    Named("pet") = pet,
    Named("Catm") = Catm,
    Named("Patm") = Patm[day-1]);
  
  DataFrame cohorts = Rcpp::as<Rcpp::DataFrame>(x["cohorts"]);
  DataFrame soil = Rcpp::as<Rcpp::DataFrame>(x["soil"]);
  int nlayers = soil.nrow();
  int numCohorts = cohorts.nrow();
  List transpOutput = basicTranspirationCommunicationOutput(numCohorts, nlayers);
  transpirationBasic(transpOutput, x, meteovec, elevation, modifyInput);

  List transpBasic = copyBasicTranspirationOutput(transpOutput, x);
    
  return(transpBasic);
} 


