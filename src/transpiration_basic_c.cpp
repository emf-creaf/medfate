#include <RcppArmadillo.h>
#include "biophysicsutils_c.h"
#include "control_c.h"
#include "communication_structures_c.h"
#include "windextinction_c.h"
#include "forestutils_c.h"
#include "modelInput_c.h"
#include "hydraulics_c.h"
#include "tissuemoisture_c.h"
#include "transpiration_basic_c.h"
#include "lightextinction_basic_c.h"
#include "soil_c.h"
#include <meteoland.h>

using namespace Rcpp;


//Plant volume in l·m-2 ground = mm
double plantVol_c(double plantPsi, ParamsVolume pars) {
  
  double leafrwc = tissueRelativeWaterContent_c(plantPsi, pars.leafpi0, pars.leafeps, 
                                              plantPsi, pars.stem_c, pars.stem_d, 
                                              pars.leafaf);
  double stemrwc = tissueRelativeWaterContent_c(plantPsi, pars.stempi0, pars.stemeps, 
                                              plantPsi, pars.stem_c, pars.stem_d, 
                                              pars.stemaf);
  return((pars.Vleaf * leafrwc * pars.LAI) + (pars.Vsapwood * stemrwc * pars.LAIlive));
}


double findNewPlantPsiConnected_c(double flowFromRoots, double plantPsi, double rootCrownPsi,
                                  ParamsVolume parsVol){
  if(plantPsi==rootCrownPsi) return(plantPsi);
  double V = plantVol_c(plantPsi, parsVol);
  //More negative rootCrownPsi causes increased flow due to water being removed
  double psiStep = rootCrownPsi - plantPsi;
  double Vnew = plantVol_c(plantPsi + psiStep, parsVol);
  while(std::abs(Vnew - V) > flowFromRoots) {
    psiStep = psiStep/2.0;
    Vnew = plantVol_c(plantPsi + psiStep, parsVol);
  }
  return(plantPsi+psiStep);
}


void transpirationBasic_c(BasicTranspiration_RESULT& BTres, BasicTranspiration_COMM& BT_comm, ModelInput& x, 
                          const WeatherInputVector& meteovec,  const double elevation) {
  
  
  // Should have internal communication structures for output
  StandBasicTranspiration_RESULT& outputStand = BTres.stand;
  PlantsBasicTranspiration_RESULT& outputPlants = BTres.plants;
  arma::mat& outputExtraction = BTres.extraction;
  std::vector<arma::mat>& outputExtractionPools = BTres.extractionPools;

  //Control parameters
  std::string& rhizosphereOverlap = x.control.rhizosphereOverlap;
  std::string& stemCavitationRecovery = x.control.commonWB.stemCavitationRecovery;
  std::string& leafCavitationRecovery = x.control.commonWB.leafCavitationRecovery;
  double cavitationRecoveryMaximumRate = x.control.commonWB.cavitationRecoveryMaximumRate;
  std::string& soilFunctions = x.control.commonWB.soilFunctions;
  double verticalLayerSize = x.control.commonWB.verticalLayerSize;
  double fullRhizosphereOverlapConductivity = x.control.commonWB.fullRhizosphereOverlapConductivity;
  bool plantWaterPools = (rhizosphereOverlap!="total");
  double hydraulicRedistributionFraction = x.control.basicWB.hydraulicRedistributionFraction;
  std::string& lfmcComponent = x.control.fireHazard.lfmcComponent;

  //Soil
  Soil& soil = x.soil;

  //Meteo input
  double pet = meteovec.pet;
  double rhmax = meteovec.rhmax;
  double rhmin = meteovec.rhmin;
  double tmax = meteovec.tmax;
  double tmin = meteovec.tmin;
  double Catm = meteovec.Catm;
  double Patm = meteovec.Patm;

  //Atmospheric pressure (if missing)
  if(std::isnan(Patm)) Patm = meteoland::utils_atmosphericPressure(elevation);

  //Daily average water vapor pressure at the atmosphere (kPa)
  double vpatm = meteoland::utils_averageDailyVP(tmin, tmax, rhmin, rhmax);
  double vpd = std::max(0.0, meteoland::utils_saturationVP((tmin+tmax)/2.0) - vpatm);


  // Canopy
  int ncanlayers = x.canopy.zlow.size();

  //Vegetation input
  std::vector<double>& LAIlive = x.above.LAI_live;
  std::vector<double>& LAIphe = x.above.LAI_expanded;
  std::vector<double>& LAIdead = x.above.LAI_dead;
  std::vector<double>& H = x.above.H;
  std::vector<double>& CR = x.above.CR;
  int numCohorts = LAIphe.size();

  //Root distribution input
  arma::mat& V = x.belowLayers.V;
  
  //Water pools
  arma::mat& Wpool = x.belowLayers.Wpool;
  std::vector<arma::mat>& RHOP = x.belowLayers.RHOP;
  std::vector<double>& poolProportions = x.below.poolProportions;
  
  //Phenology parameters
  const std::vector<std::string>& phenoType = x.paramsPhenology.phenoType;

  //Parameters
  const std::vector<double>& Al2As = x.paramsAnatomy.Al2As;
  const std::vector<double>& r635 = x.paramsAnatomy.r635;
  
  const std::vector<double>& kPAR = x.paramsInterception.kPAR;
  const std::vector<double>& kSWR = x.paramsInterception.kSWR;
  
  const std::vector<double>& Gswmin = x.paramsTranspiration.Gswmin;
  const std::vector<double>& Psi_Extract = x.paramsTranspiration.Psi_Extract;
  const std::vector<double>& Exp_Extract = x.paramsTranspiration.Exp_Extract;
  const std::vector<double>& VCstem_c = x.paramsTranspiration.VCstem_c;
  const std::vector<double>& VCstem_d = x.paramsTranspiration.VCstem_d;
  const std::vector<double>& VCleaf_c = x.paramsTranspiration.VCleaf_c;
  const std::vector<double>& VCleaf_d = x.paramsTranspiration.VCleaf_d;
  const std::vector<double>& WUE = x.paramsTranspiration.WUE;
  const std::vector<double>& WUE_par = x.paramsTranspiration.WUE_par;
  const std::vector<double>& WUE_co2 = x.paramsTranspiration.WUE_co2;
  const std::vector<double>& WUE_vpd = x.paramsTranspiration.WUE_vpd;
  const std::vector<double>& Tmax_LAI = x.paramsTranspiration.Tmax_LAI;
  const std::vector<double>& Tmax_LAIsq = x.paramsTranspiration.Tmax_LAIsq;
  
  const std::vector<double>& maxFMC = x.paramsWaterStorage.maxFMC;
  const std::vector<double>& maxMCstem = x.paramsWaterStorage.maxMCstem;
  const std::vector<double>& maxMCleaf = x.paramsWaterStorage.maxMCleaf;
  const std::vector<double>& StemPI0 = x.paramsWaterStorage.StemPI0;
  const std::vector<double>& StemEPS = x.paramsWaterStorage.StemEPS;
  const std::vector<double>& StemAF = x.paramsWaterStorage.StemAF;
  const std::vector<double>& Vsapwood = x.paramsWaterStorage.Vsapwood; //l·m-2 = mm
  const std::vector<double>& LeafPI0 = x.paramsWaterStorage.LeafPI0;
  const std::vector<double>& LeafEPS = x.paramsWaterStorage.LeafEPS;
  const std::vector<double>& LeafAF = x.paramsWaterStorage.LeafAF;
  const std::vector<double>& Vleaf = x.paramsWaterStorage.Vleaf; //l·m-2 = mm

  
  //Internal state variables
  std::vector<double>& phi = x.internalPhenology.phi;
  std::vector<double>& PlantPsi = x.internalWater.PlantPsi;
  std::vector<double>& StemPLC = x.internalWater.StemPLC;
  std::vector<double>& LeafPLC = x.internalWater.LeafPLC;
  arma::mat& LAIme = x.internalLAIDistribution.expanded;
  arma::mat& LAImd = x.internalLAIDistribution.dead;
  std::vector<double>& PrevLAIexpanded = x.internalLAIDistribution.PrevLAIexpanded;
  std::vector<double>& PrevLAIdead = x.internalLAIDistribution.PrevLAIdead;
  std::vector<double>& PARcohort = x.internalLAIDistribution.PARcohort;

  //Determine whether leaves are out (phenology) and the adjusted Leaf area
  double s = 0.0, LAIcell = 0.0, LAIcelllive = 0.0, LAIcellexpanded = 0.0,LAIcelldead = 0.0;
  double sum_abs_exp = 0.0, sum_abs_dead = 0.0;
  for(int c=0;c<numCohorts;c++) {
    s += (kPAR[c]*(LAIphe[c]+LAIdead[c]));
    LAIcell += LAIphe[c]+LAIdead[c];
    LAIcelldead += LAIdead[c];
    LAIcellexpanded +=LAIphe[c];
    LAIcelllive += LAIlive[c];
    sum_abs_exp += std::abs(LAIphe[c] - PrevLAIexpanded[c]);
    sum_abs_dead += std::abs(LAIdead[c] - PrevLAIdead[c]);
  }

  if(numCohorts>0) {
    bool recalc_LAI = false;
    if(std::isnan(PrevLAIexpanded[0]) || std::isnan(PrevLAIdead[0])) {
      recalc_LAI = true;
    } else{
      if(sum_abs_exp>0.001) {
        recalc_LAI = true;
      } else {
        if(sum_abs_dead>0.001) recalc_LAI = true;
      }
    }
    if(recalc_LAI) {
      std::vector<double> z(ncanlayers+1,0.0);
      for(int i=1;i<=ncanlayers;i++) z[i] = z[i-1] + verticalLayerSize;
      for(int i=0; i<numCohorts;i++) {
        PARcohort[i] = availableLight_c(H[i]*(1.0-(1.0-CR[i])/2.0), H, LAIphe, LAIdead, kPAR, CR);
        PrevLAIexpanded[i] = LAIphe[i];
        PrevLAIdead[i] = LAIdead[i];
      }
      //Update LAI distribution if necessary
      updateLAIdistributionVectors_c(LAIme, z, LAIphe, H, CR);
      updateLAIdistributionVectors_c(LAImd, z, LAIdead, H, CR);
    }
  }
  // NumericVector CohASWRF = cohortAbsorbedSWRFraction(LAIme, LAImd, kSWR);
  // CohASWRF = pow(CohASWRF, 0.75);
  // 
  // //Apply fractions to potential evapotranspiration
  // //Maximum canopy transpiration
  // //    Tmax = PET[i]*(-0.006*pow(LAIcell[i],2.0)+0.134*LAIcell[i]+0.036); //From Granier (1999)
  // NumericVector Tmax = pet*(Tmax_LAIsq*(LAIcell*LAIcell)+ Tmax_LAI*LAIcell); //From Granier (1999)
  // 
  // //Fraction of Tmax attributed to each plant cohort
  // double pabs = std::accumulate(CohASWRF.begin(),CohASWRF.end(),0.0);
  // NumericVector TmaxCoh(numCohorts,0.0);
  // if(pabs>0.0) TmaxCoh = Tmax*(CohASWRF/pabs);
  // 
  // //Actual plant transpiration
  // int nlayers = Wpool.ncol();
  // 
  // NumericVector Eplant(numCohorts, 0.0), Agplant(numCohorts, 0.0);
  // NumericVector DDS(numCohorts, 0.0), LFMC(numCohorts, 0.0);
  // NumericVector PLCm(numCohorts), RWCsm(numCohorts), RWClm(numCohorts),RWCssm(numCohorts), RWClsm(numCohorts);
  // NumericVector PWB(numCohorts,0.0);
  // NumericVector Kl, epc, Vl;
  // 
  // //Calculate unsaturated conductivity (mmolH20·m-1·s-1·MPa-1)
  // NumericVector Kunsat = conductivity(soil, soilFunctions, true);
  // //Calculate soil water potential
  // NumericVector psiSoil = psi(soil,soilFunctions);
  // 
  // 
  // NumericMatrix WaterM(numCohorts, nlayers);
  // NumericMatrix KunsatM(numCohorts, nlayers);
  // NumericMatrix psiSoilM(numCohorts, nlayers);
  // if(plantWaterPools) {
  //   DataFrame soil_pool = clone(soil);
  //   NumericVector Ws_pool = soil_pool["W"];
  //   for(int j = 0; j<numCohorts;j++) {
  //     //Copy values of soil moisture from pool of cohort j
  //     for(int l = 0; l<nlayers;l++) Ws_pool[l] = Wpool(j,l);
  //     //Calculate unsaturated conductivity (mmolH20·m-1·s-1·MPa-1)
  //     KunsatM(j,_) = conductivity(soil_pool, soilFunctions, true);
  //     //Calculate soil water potential
  //     psiSoilM(j,_) = psi(soil_pool, soilFunctions);
  //     // Rcout<< Ws_pool[3]<< " "<< psiSoilM(j,3)<<"\n";
  //     //Calculate available water
  //     WaterM(j,_) = water(soil_pool, soilFunctions);
  //   }
  // }
  // // for(int i= 0;i<psiSoil.size();i++) Rcout<< "S "<<i<<" "<<psiSoil[i]<<"\n";
  // 
  // ParamsVolume parsVol;
  // 
  // for(int c=0;c<numCohorts;c++) {
  //   
  //   parsVol.stem_c = VCstem_c[c];
  //   parsVol.stem_d = VCstem_d[c];
  //   parsVol.leafpi0 = LeafPI0[c];
  //   parsVol.leafeps = LeafEPS[c];
  //   parsVol.leafaf = LeafAF[c];
  //   parsVol.stempi0 = StemPI0[c];
  //   parsVol.stemeps = StemEPS[c];
  //   parsVol.stemaf = StemAF[c];
  //   parsVol.Vsapwood = Vsapwood[c];
  //   parsVol.Vleaf = Vleaf[c];
  //   parsVol.LAI = LAIphe[c];
  //   parsVol.LAIlive = LAIlive[c];
  //   
  //   double rootCrownPsi = NA_REAL;
  //   
  //   //Cuticular transpiration    
  //   double lvp_tmax = leafVapourPressure_c(tmax,  PlantPsi[c]);
  //   double lvp_tmin = leafVapourPressure_c(tmin,  PlantPsi[c]);
  //   double lvpd_tmax = std::max(0.0, lvp_tmax - vpatm);
  //   double lvpd_tmin = std::max(0.0, lvp_tmin - vpatm);
  //   double E_gmin = Gswmin[c]*(lvpd_tmin+lvpd_tmax)/(2.0*Patm); // mol·s-1·m-2
  //   double E_gmin_day = E_gmin*LAIphe[c]*(24.0*3600.0*0.018); //L·d-1·m-2 = mm·d-1
  //   
  //   //Extraction from soil (can later be modified if there are changes in plant water content)
  //   if(!plantWaterPools) {
  //     NumericVector Klc(nlayers);
  //     NumericVector Kunlc(nlayers);
  //     for(int l=0;l<nlayers;l++) {
  //       Klc[l] = Psi2K_c(psiSoil[l], Psi_Extract[c], Exp_Extract[c]);
  //       //Limit Mean Kl due to previous cavitation
  //       if(stemCavitationRecovery!="total") {
  //         Klc[l] = std::min(Klc[l], 1.0-StemPLC[c]); 
  //       }
  //       Kunlc[l] = std::sqrt(Kunsat[l])*V(c,l);
  //     }
  //     double sumKunlc = sum(Kunlc);
  //     double Klcmean = sum(Klc*V(c,_));
  //     for(int l=0;l<nlayers;l++) {
  //       outputExtraction(c,l) = std::max(TmaxCoh[c]*Klcmean, E_gmin_day)*(Kunlc[l]/sumKunlc);
  //     }
  //     rootCrownPsi = averagePsi(psiSoil, V(c,_), Exp_Extract[c], Psi_Extract[c]);
  //   } else {
  //     NumericMatrix ExtractionPoolsCoh = outputExtractionPools[c]; //this is used to store extraction of a SINGLE plant cohort from all pools
  //     
  //     NumericMatrix RHOPcoh = Rcpp::as<Rcpp::NumericMatrix>(RHOP[c]);
  //     NumericMatrix RHOPcohDyn(numCohorts, nlayers);
  //     for(int l=0;l<nlayers;l++) {
  //       RHOPcohDyn(c,l) = RHOPcoh(c,l);
  //       for(int j=0; j<numCohorts;j++) {
  //         if(j!=c) {
  //           double overlapFactor = std::min(1.0, KunsatM(j,l)/(cmdTOmmolm2sMPa*fullRhizosphereOverlapConductivity));
  //           RHOPcohDyn(j,l) = RHOPcoh(j,l)*overlapFactor;
  //           RHOPcohDyn(c,l) = RHOPcohDyn(c,l) + (RHOPcoh(j,l) - RHOPcohDyn(j,l));
  //         } 
  //       }
  //     }
  //     NumericMatrix Klc(numCohorts, nlayers);
  //     NumericMatrix Kunlc(numCohorts, nlayers);
  //     NumericMatrix RHOPcohV(numCohorts, nlayers);
  //     for(int j = 0;j<numCohorts;j++) {
  //       for(int l=0;l<nlayers;l++) {
  //         RHOPcohV(j,l) = RHOPcohDyn(j,l)*V(c,l);
  //         Klc(j,l) = Psi2K_c(psiSoilM(c,l), Psi_Extract[c], Exp_Extract[c]);
  //         //Limit Mean Kl due to previous cavitation
  //         if(stemCavitationRecovery!="total") Klc(j,l) = std::min(Klc(j,l), 1.0-StemPLC[c]); 
  //         Kunlc(j,l) = std::sqrt(KunsatM(j,l))*RHOPcohV(j,l);
  //       }
  //     }
  //     double sumKunlc = sum(Kunlc);
  //     double Klcmean = sum(Klc*RHOPcohV);
  //     for(int l=0;l<nlayers;l++) {
  //       for(int j = 0;j<numCohorts;j++) {
  //         ExtractionPoolsCoh(j,l) = std::max(TmaxCoh[c]*Klcmean, E_gmin_day)*(Kunlc(j,l)/sumKunlc);
  //       }
  //       outputExtraction(c,l) = sum(ExtractionPoolsCoh(_,l)); // Sum extraction from all pools (layer l)
  //     }
  //     
  //     rootCrownPsi = averagePsiPool(psiSoilM, RHOPcohV, Exp_Extract[c], Psi_Extract[c]);
  //     // Rcout<< c << " : "<< psiSoilM(c,0) << " " << psiSoilM(c,1) << " " << psiSoilM(c,2) << " " << psiSoilM(c,3) << " " << rootCrownPsi<<"\n";
  //     // Rcout<< c << " : "<< RHOPcohV(c,0) << " " << RHOPcohV(c,1) << " " << RHOPcohV(c,2) << " " << RHOPcohV(c,3) << " " << rootCrownPsi<<"\n";
  //   }
  // 
  // 
  //   double oldVol = plantVol(PlantPsi[c], parsVol); 
  //   
  //   //Transpiration is the maximum of predicted extraction and cuticular transpiration
  //   double ext_sum = sum(outputExtraction(c,_));
  //   Eplant[c] = ext_sum;
  //   // PlantPsi[c] = findNewPlantPsiConnected(Eplant[c], PlantPsi[c], rootCrownPsi, parsVol);
  //   //For deciduous species, make water potential follow soil during winter
  //   // if(LAIphe[c]==0.0) PlantPsi[c] = rootCrownPsi;
  //   PlantPsi[c] = rootCrownPsi;
  //   // PlantPsi[c] = rootCrownPsi;
  //   double newVol = plantVol(PlantPsi[c], parsVol);
  //   
  //   double volDiff = newVol - oldVol;
  //   //Plant transpiration and water balance
  //   PWB[c] = volDiff;
  //   
  //   //Photosynthesis
  //   double fpar = std::min(1.0, pow(PARcohort[c]/100.0,WUE_par[c]));
  //   double fco2 = (1.0 - exp((-1.0)*WUE_co2[c]*Catm));
  //   double fvpd = pow(vpd, WUE_vpd[c]);
  //   if(vpd < 0.25) fvpd = 2.5 - (2.5 - pow(0.25, WUE_vpd[c]))*(vpd/0.25);
  //   // Rcout<<fpar<<" "<< fco2 << " "<< fvpd<< " "<< WUE[c]*fpar*fco2*fvpd<<"\n";
  //   Agplant[c] = WUE[c]*Eplant[c]*fpar*fco2*fvpd;
  // }
  // 
  // //Plant water status (StemPLC, RWC, DDS)
  // for(int c=0;c<numCohorts;c++) {
  //   if(stemCavitationRecovery!="total") {
  //     StemPLC[c] = std::max(1.0 - xylemConductance_c(PlantPsi[c], 1.0, VCstem_c[c], VCstem_d[c]), StemPLC[c]); //Track current embolism if no refill
  //   } else {
  //     StemPLC[c] = 1.0 - xylemConductance_c(PlantPsi[c], 1.0, VCstem_c[c], VCstem_d[c]);
  //   }
  //   if(leafCavitationRecovery!="total") {
  //     LeafPLC[c] = std::max(1.0 - xylemConductance_c(PlantPsi[c], 1.0, VCleaf_c[c], VCleaf_d[c]), LeafPLC[c]); //Track current embolism if no refill
  //   } else {
  //     LeafPLC[c] = 1.0 - xylemConductance_c(PlantPsi[c], 1.0, VCleaf_c[c], VCleaf_d[c]);
  //   }
  //   
  //   //Relative water content and fuel moisture from plant water potential
  //   RWClm[c] =  tissueRelativeWaterContent_c(PlantPsi[c], LeafPI0[c], LeafEPS[c], 
  //                                          PlantPsi[c], VCstem_c[c], VCstem_d[c], 
  //                                          LeafAF[c]);
  //   RWCsm[c] =  tissueRelativeWaterContent_c(PlantPsi[c], StemPI0[c], StemEPS[c], 
  //                                          PlantPsi[c], VCstem_c[c], VCstem_d[c], 
  //                                          StemAF[c]);
  //   // The fraction of leaves will decrease due to phenology or processes leading to defoliation
  //   double fleaf = (1.0/r635[c])*(LAIphe[c]/LAIlive[c]);
  //   if(lfmcComponent=="fine") { //fine fuel moisture
  //     LFMC[c] = maxMCleaf[c]*RWClm[c]*fleaf + maxMCstem[c]*RWCsm[c]*(1.0 - fleaf);
  //   } else { //"leaf"
  //     LFMC[c] = maxFMC[c]*RWClm[c];
  //   }
  //   
  //   //Daily drought stress from plant WP
  //   DDS[c] = (1.0 - Psi2K_c(PlantPsi[c],Psi_Extract[c],Exp_Extract[c])); 
  //   if(phenoType[c] == "winter-deciduous" || phenoType[c] == "winter-semideciduous") DDS[c] = phi[c]*DDS[c];
  //     
  //   double SAmax = 10e4/Al2As[c]; //cm2·m-2 of leaf area
  //   double r = cavitationRecoveryMaximumRate*std::max(0.0, (PlantPsi[c] + 1.5)/1.5);
  //   if(stemCavitationRecovery=="rate") {
  //     StemPLC[c] = std::max(0.0, StemPLC[c] - (r/SAmax));
  //   }
  //   if(leafCavitationRecovery=="rate") {
  //     LeafPLC[c] = std::max(0.0, LeafPLC[c] - (r/SAmax));
  //   }
  // }
  // 
  // //Atempt to implement hydraulic redistribution
  // if(hydraulicRedistributionFraction > 0.0) {
  //   if(!plantWaterPools) {
  //     for(int c=0;c<numCohorts;c++) {
  //       double redAmount = Eplant[c]*hydraulicRedistributionFraction;
  //       // Rcout<<c<< "red amount"<< redAmount;
  //       NumericVector Ws = soil["W"];
  //       NumericVector SW = water(soil, soilFunctions);
  //       double soilRWC = sum(Ws*SW)/sum(SW);
  //       NumericVector WDiff = Ws - soilRWC; 
  //       NumericVector DonorDiff = pmax(0.0, WDiff);
  //       NumericVector ReceiverDiff = pmax(0.0, -WDiff);
  //       NumericVector HD(nlayers,0.0);
  //       if(sum(DonorDiff)>0.0) {
  //         for(int l=0;l<nlayers;l++) {
  //           if(WDiff[l]>0.0) {
  //             HD[l] = redAmount*DonorDiff[l]/sum(DonorDiff);
  //           } else{
  //             HD[l] = -redAmount*ReceiverDiff[l]/sum(ReceiverDiff);
  //           }
  //           // Rcout<<" "<<l<<" "<<HD[l];
  //           outputExtraction(c,l) += HD[l];
  //         }
  //       }
  //       // Rcout<< "\n";
  //     }
  //   } else {
  //     for(int c=0;c<numCohorts;c++) {
  //       NumericMatrix ExtractionPoolsCoh = outputExtractionPools[c];
  //       for(int j=0;j<numCohorts;j++) {
  //         double redAmount = sum(ExtractionPoolsCoh(j,_))*hydraulicRedistributionFraction;
  //         NumericVector Ws = Wpool(j,_);
  //         NumericVector SW = WaterM(j,_);
  //         double soilRWC = sum(Ws*SW)/sum(SW);
  //         NumericVector WDiff = Ws - soilRWC; 
  //         NumericVector DonorDiff = pmax(0.0, WDiff);
  //         NumericVector ReceiverDiff = pmax(0.0, -WDiff);
  //         NumericVector HD(nlayers,0.0);
  //         if(sum(DonorDiff)>0.0) {
  //           for(int l=0;l<nlayers;l++) {
  //             if(WDiff[l]>0.0) {
  //               HD[l] = redAmount*DonorDiff[l]/sum(DonorDiff);
  //             } else{
  //               HD[l] = -redAmount*ReceiverDiff[l]/sum(ReceiverDiff);
  //             }
  //             outputExtraction(c,l) += HD[l];
  //             ExtractionPoolsCoh(j,l) += HD[l];
  //           }
  //         }
  //       }
  //     }
  //   }
  // }
  // 
  // 
  // // Copy output stand
  // outputStand["LAI"] = LAIcell;
  // outputStand["LAIlive"] = LAIcelllive;
  // outputStand["LAIexpanded"] = LAIcellexpanded;
  // outputStand["LAIdead"] = LAIcelldead;
  // 
  // // Copy output plants
  // NumericVector outputLAI = outputPlants["LAI"];
  // NumericVector outputLAIlive = outputPlants["LAIlive"];
  // NumericVector outputFPAR = outputPlants["FPAR"];
  // NumericVector outputAbsorbedSWRFraction = outputPlants["AbsorbedSWRFraction"];
  // NumericVector outputExtractionByPlants = outputPlants["Extraction"];
  // NumericVector outputTranspirationByPlants = outputPlants["Transpiration"];
  // NumericVector outputGrossPhotosynthesis = outputPlants["GrossPhotosynthesis"];
  // NumericVector outputPlantPsi = outputPlants["PlantPsi"];
  // NumericVector outputDDS = outputPlants["DDS"];
  // NumericVector outputStemRWC = outputPlants["StemRWC"];
  // NumericVector outputLeafRWC = outputPlants["LeafRWC"];
  // NumericVector outputLFMC = outputPlants["LFMC"];
  // NumericVector outputStemPLC = outputPlants["StemPLC"];
  // NumericVector outputLeafPLC = outputPlants["LeafPLC"];
  // NumericVector outputWaterBalance = outputPlants["WaterBalance"];
  // for(int c =0;c<numCohorts;c++) {
  //   outputLAI[c] = LAIphe[c];
  //   outputLAIlive[c] = LAIlive[c];
  //   outputFPAR[c] = PARcohort[c];
  //   outputAbsorbedSWRFraction[c] = 100.0*CohASWRF[c];
  //   outputExtractionByPlants[c] = sum(outputExtraction(c,_));
  //   outputTranspirationByPlants[c] = Eplant[c];
  //   outputGrossPhotosynthesis[c] = Agplant[c];
  //   outputPlantPsi[c] = PlantPsi[c];
  //   outputDDS[c] = DDS[c];
  //   outputStemRWC[c] = RWCsm[c];
  //   outputLeafRWC[c] = RWClm[c];
  //   outputLFMC[c] = LFMC[c];
  //   outputStemPLC[c] = StemPLC[c];
  //   outputLeafPLC[c] = LeafPLC[c];
  //   outputWaterBalance[c] = PWB[c];
  // }

}


Rcpp::List copyBasicTranspirationOutput_c(const BasicTranspiration_RESULT& btc, ModelInput& x) {
  const std::string& rhizosphereOverlap = x.control.rhizosphereOverlap;
  bool plantWaterPools = (rhizosphereOverlap!="total");
  int nlayers = x.soil.getNlayers();
  int numCohorts = x.cohorts.CohortCode.size();
  
  const arma::mat& extractionComm = btc.extraction;
  Rcpp::NumericMatrix Extraction = copyNumericMatrix_c(extractionComm, numCohorts, nlayers); // this is final extraction of each cohort from each layer
  Extraction.attr("dimnames") = Rcpp::List::create(x.cohorts.CohortCode, Rcpp::seq(1,nlayers));
  
  Rcpp::List ExtractionPools(numCohorts);
  const std::vector< arma::mat>& ExtractionPoolsComm = btc.extractionPools;
  if(plantWaterPools) {
    for(int c=0;c<numCohorts;c++) {
      const arma::mat& extractionPoolsCohComm = ExtractionPoolsComm[c];
      Rcpp::NumericMatrix ExtractionPoolsCohComm_c = copyNumericMatrix_c(extractionPoolsCohComm, numCohorts, nlayers); // this is final extraction of each cohort from each layer
      ExtractionPoolsCohComm_c.attr("dimnames") = Rcpp::List::create(x.cohorts.CohortCode, Rcpp::seq(1,nlayers));
      ExtractionPools[c] = ExtractionPoolsCohComm_c;
    }
    ExtractionPools.attr("names") = x.cohorts.CohortCode;
  }
  
  NumericVector standVEC = Rcpp::NumericVector::create(_["LAI"] = btc.stand.LAI,
                                                 _["LAIlive"] = btc.stand.LAIlive, 
                                                 _["LAIexpanded"] = btc.stand.LAIexpanded, 
                                                 _["LAIdead"] = btc.stand.LAIdead);
  
  DataFrame cohortsDF = DataFrame::create(
    _["SP"] = Rcpp::wrap(x.cohorts.SpeciesIndex),
    _["Name"]= Rcpp::wrap(x.cohorts.SpeciesName)
  );
  cohortsDF.attr("row.names") = x.cohorts.CohortCode;
  
  DataFrame plantsDF = DataFrame::create(
    _["LAI"] = Rcpp::wrap(btc.plants.LAI),
    _["LAIlive"] = Rcpp::wrap(btc.plants.LAIlive),
    _["FPAR"] = Rcpp::wrap(btc.plants.FPAR),
    _["AbsorbedSWRFraction"] = Rcpp::wrap(btc.plants.AbsorbedSWRFraction),
    _["Extraction"] = Rcpp::wrap(btc.plants.Extraction),
    _["Transpiration"] = Rcpp::wrap(btc.plants.Transpiration),
    _["GrossPhotosynthesis"] = Rcpp::wrap(btc.plants.GrossPhotosynthesis),
    _["PlantPsi"] = Rcpp::wrap(btc.plants.PlantPsi),
    _["DDS"] = Rcpp::wrap(btc.plants.DDS),
    _["StemRWC"] = Rcpp::wrap(btc.plants.StemRWC),
    _["LeafRWC"] = Rcpp::wrap(btc.plants.LeafRWC),
    _["LFMC"] = Rcpp::wrap(btc.plants.LFMC),
    _["StemPLC"] = Rcpp::wrap(btc.plants.StemPLC),
    _["LeafPLC"] = Rcpp::wrap(btc.plants.LeafPLC),
    _["WaterBalance"] = Rcpp::wrap(btc.plants.WaterBalance)
  );
  plantsDF.attr("row.names") = x.cohorts.CohortCode;
  
  List l = List::create(_["cohorts"] = cohortsDF,
                        _["Stand"] = standVEC,
                        _["Plants"] = plantsDF,
                        _["Extraction"] = Extraction,
                        _["ExtractionPools"] = ExtractionPools);
  return(l);
}
