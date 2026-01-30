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
  double verticalLayerSize = x.control.commonWB.verticalLayerSize;
  double fullRhizosphereOverlapConductivity = x.control.commonWB.fullRhizosphereOverlapConductivity;
  bool plantWaterPools = (rhizosphereOverlap!="total");
  double hydraulicRedistributionFraction = x.control.basicWB.hydraulicRedistributionFraction;
  std::string& lfmcComponent = x.control.fireHazard.lfmcComponent;

  //Soil
  Soil& soil = x.soil;
  int nlayers = soil.getNlayers();

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

  
  // Internal communication vectors
  std::vector<double>& CohASWRF = BT_comm.CohASWRF;
  std::vector<double>& Tmax = BT_comm.Tmax;
  std::vector<double>& TmaxCoh = BT_comm.TmaxCoh;
  arma::mat& RHOPCohDyn = BT_comm.RHOPCohDyn;
  for(int c=0;c<numCohorts;c++) {
    CohASWRF[c] = 0.0;
    Tmax[c] = 0.0;
    TmaxCoh[c] = 0.0;
  }
  
  //Output vectors
  std::vector<double>&  Eplant = outputPlants.Transpiration;
  std::vector<double>&  Agplant = outputPlants.GrossPhotosynthesis;
  std::vector<double>&  DDS = outputPlants.DDS;
  std::vector<double>&  LFMC = outputPlants.LFMC;
  std::vector<double>&  StemRWC = outputPlants.StemRWC;
  std::vector<double>&  LeafRWC = outputPlants.LeafRWC;
  std::vector<double>&  PWB = outputPlants.WaterBalance;
  std::vector<double>&  Extraction = outputPlants.Extraction;

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
  
  //Cohort absorbed fraction
  cohortAbsorbedSWRFraction_c(CohASWRF, BT_comm.AbSWRcomm, LAIme, LAImd, kSWR);
  for(int c=0;c<numCohorts;c++) {
    CohASWRF[c] = pow(CohASWRF[c], 0.75);
  }
  
  //Apply fractions to potential evapotranspiration
  //Maximum canopy transpiration
  //    Tmax = PET[i]*(-0.006*pow(LAIcell[i],2.0)+0.134*LAIcell[i]+0.036); //From Granier (1999)
  for(int c=0;c<numCohorts;c++) {
    Tmax[c] = pet*(Tmax_LAIsq[c]*(LAIcell*LAIcell)+ Tmax_LAI[c]*LAIcell); //From Granier (1999)
  }
  
  //Fraction of Tmax attributed to each plant cohort
  double pabs = std::accumulate(CohASWRF.begin(),CohASWRF.end(),0.0);
  if(pabs>0.0) {
    for(int c=0;c<numCohorts;c++) {
      TmaxCoh[c] = Tmax[c]*(CohASWRF[c]/pabs); 
    }
  }

  arma::mat WaterM;
  arma::mat KunsatM;
  arma::mat psiSoilM;
  if(plantWaterPools) {
    // Initialize psiSoilM, KunsatM, WaterM
    WaterM = arma::mat(numCohorts, nlayers);
    KunsatM = arma::mat(numCohorts, nlayers);
    psiSoilM = arma::mat(numCohorts, nlayers);
    //Store overall soil moisture in a backup copy
    double* Wbackup = new double[nlayers];
    for(int l = 0; l<nlayers;l++) Wbackup[l] = soil.getW(l);
    //   DataFrame soil_pool = clone(soil);
    //   NumericVector Ws_pool = soil_pool["W"];
    for(int j = 0; j<numCohorts;j++) {
      //Copy values of soil moisture from pool of cohort j to general soil
      for(int l = 0; l<nlayers;l++) {
        soil.setW(l,Wpool(j,l)); // this updates psi, theta, ... 
        psiSoilM(j,l) = soil.getPsi(l);
        KunsatM(j,l) = soil.getConductivity(l, true);
        WaterM(j,l) = soil.getWater(l);
      }
    }
    //Restore soil moisture
    for(int l = 0; l<nlayers;l++) soil.setW(l, Wbackup[l]);
    //Delete backup
    delete[] Wbackup;
  }


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
    double lvp_tmax = leafVapourPressure_c(tmax,  PlantPsi[c]);
    double lvp_tmin = leafVapourPressure_c(tmin,  PlantPsi[c]);
    double lvpd_tmax = std::max(0.0, lvp_tmax - vpatm);
    double lvpd_tmin = std::max(0.0, lvp_tmin - vpatm);
    double E_gmin = Gswmin[c]*(lvpd_tmin+lvpd_tmax)/(2.0*Patm); // mol·s-1·m-2
    double E_gmin_day = E_gmin*LAIphe[c]*(24.0*3600.0*0.018); //L·d-1·m-2 = mm·d-1

    //Extraction from soil (can later be modified if there are changes in plant water content)
    if(!plantWaterPools) {
      std::vector<double> psiSoil(nlayers);
      std::vector<double> V_c(nlayers, 0.0);
      double* Klc = new double[nlayers];
      double* Kunlc = new double[nlayers];
      double Klcmean = 0.0;
      double sumKunlc = 0.0;
      for(int l=0;l<nlayers;l++) {
        psiSoil[l] = soil.getPsi(l);
        V_c[l] = V(c,l);
        Klc[l] = Psi2K_c(psiSoil[l], Psi_Extract[c], Exp_Extract[c]);
        //Limit Mean Kl due to previous cavitation
        if(stemCavitationRecovery!="total") {
          Klc[l] = std::min(Klc[l], 1.0-StemPLC[c]);
        }
        Klcmean += Klc[l]*V_c[l];
        Kunlc[l] = std::sqrt(soil.getConductivity(l,true))*V_c[l];
        sumKunlc += Kunlc[l];
      }
      for(int l=0;l<nlayers;l++) {
        outputExtraction(c,l) = std::max(TmaxCoh[c]*Klcmean, E_gmin_day)*(Kunlc[l]/sumKunlc);
      }
      rootCrownPsi = averagePsi_c(psiSoil, V_c, Exp_Extract[c], Psi_Extract[c]);
      delete[] Klc;
      delete[] Kunlc;
    } else {
      //Calculate dynamic overlap      
      arma::mat& ExtractionPoolsCoh = outputExtractionPools[c]; //this is used to store extraction of a SINGLE plant cohort from all pools
      arma::mat& RHOPcoh = RHOP[c];
      for(int l=0;l<nlayers;l++) {
        RHOPCohDyn(c,l) = RHOPcoh(c,l);
        for(int j=0; j<numCohorts;j++) {
          if(j!=c) {
            double overlapFactor = std::min(1.0, KunsatM(j,l)/(cmdTOmmolm2sMPa*fullRhizosphereOverlapConductivity));
            RHOPCohDyn(j,l) = RHOPcoh(j,l)*overlapFactor;
            RHOPCohDyn(c,l) = RHOPCohDyn(c,l) + (RHOPcoh(j,l) - RHOPCohDyn(j,l));
          }
        }
      }
      arma::mat Klc(numCohorts, nlayers);
      arma::mat Kunlc(numCohorts, nlayers);
      arma::mat RHOPcohV(numCohorts, nlayers);
      double Klcmean = 0.0;
      double sumKunlc = 0.0;
      for(int j = 0;j<numCohorts;j++) {
        for(int l=0;l<nlayers;l++) {
          RHOPcohV(j,l) = RHOPCohDyn(j,l)*V(c,l);
          Klc(j,l) = Psi2K_c(psiSoilM(c,l), Psi_Extract[c], Exp_Extract[c]);
          //Limit Mean Kl due to previous cavitation
          if(stemCavitationRecovery!="total") Klc(j,l) = std::min(Klc(j,l), 1.0-StemPLC[c]);
          Klcmean += Klc(j,l)*RHOPcohV(j,l);
          Kunlc(j,l) = std::sqrt(KunsatM(j,l))*RHOPcohV(j,l);
          sumKunlc += Kunlc(j,l);
        }
      }
      for(int l=0;l<nlayers;l++) {
        outputExtraction(c,l) = 0.0;
        for(int j = 0;j<numCohorts;j++) {
          ExtractionPoolsCoh(j,l) = std::max(TmaxCoh[c]*Klcmean, E_gmin_day)*(Kunlc(j,l)/sumKunlc);
          outputExtraction(c,l) += ExtractionPoolsCoh(j,l); // Sum extraction from all pools (layer l)
        }
      }
      rootCrownPsi = averagePsiPool_c(psiSoilM, RHOPcohV, Exp_Extract[c], Psi_Extract[c]);
      // Rcout<< c << " : "<< psiSoilM(c,0) << " " << psiSoilM(c,1) << " " << psiSoilM(c,2) << " " << psiSoilM(c,3) << " " << rootCrownPsi<<"\n";
      // Rcout<< c << " : "<< RHOPcohV(c,0) << " " << RHOPcohV(c,1) << " " << RHOPcohV(c,2) << " " << RHOPcohV(c,3) << " " << rootCrownPsi<<"\n";
    }

    double oldVol = plantVol_c(PlantPsi[c], parsVol);

    //Transpiration is now equal to extraction
    Extraction[c] = arma::sum(outputExtraction.row(c));
    Eplant[c] = Extraction[c];
    //For deciduous species, make water potential follow soil during winter
    PlantPsi[c] = rootCrownPsi;
    double newVol = plantVol_c(PlantPsi[c], parsVol);

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
      StemPLC[c] = std::max(1.0 - xylemConductance_c(PlantPsi[c], 1.0, VCstem_c[c], VCstem_d[c]), StemPLC[c]); //Track current embolism if no refill
    } else {
      StemPLC[c] = 1.0 - xylemConductance_c(PlantPsi[c], 1.0, VCstem_c[c], VCstem_d[c]);
    }
    if(leafCavitationRecovery!="total") {
      LeafPLC[c] = std::max(1.0 - xylemConductance_c(PlantPsi[c], 1.0, VCleaf_c[c], VCleaf_d[c]), LeafPLC[c]); //Track current embolism if no refill
    } else {
      LeafPLC[c] = 1.0 - xylemConductance_c(PlantPsi[c], 1.0, VCleaf_c[c], VCleaf_d[c]);
    }

    //Relative water content and fuel moisture from plant water potential
    LeafRWC[c] =  tissueRelativeWaterContent_c(PlantPsi[c], LeafPI0[c], LeafEPS[c],
                                               PlantPsi[c], VCstem_c[c], VCstem_d[c],
                                               LeafAF[c]);
    StemRWC[c] =  tissueRelativeWaterContent_c(PlantPsi[c], StemPI0[c], StemEPS[c],
                                               PlantPsi[c], VCstem_c[c], VCstem_d[c],
                                               StemAF[c]);
    // The fraction of leaves will decrease due to phenology or processes leading to defoliation
    double fleaf = (1.0/r635[c])*(LAIphe[c]/LAIlive[c]);
    if(lfmcComponent=="fine") { //fine fuel moisture
      LFMC[c] = maxMCleaf[c]*LeafRWC[c]*fleaf + maxMCstem[c]*StemRWC[c]*(1.0 - fleaf);
    } else { //"leaf"
      LFMC[c] = maxFMC[c]*LeafRWC[c];
    }

    //Daily drought stress from plant WP
    DDS[c] = (1.0 - Psi2K_c(PlantPsi[c],Psi_Extract[c],Exp_Extract[c]));
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
    double* WDiff = new double[nlayers];
    double* DonorDiff = new double[nlayers];
    double* ReceiverDiff = new double[nlayers];
    if(!plantWaterPools) {
      double* W = new double[nlayers];
      double sumWater = 0.0;
      double sumWaterW = 0.0;
      double water_l = 0.0;
      for(int l=0;l<nlayers;l++) {
        water_l = soil.getWater(l);
        W[l] = soil.getW(l);
        sumWater += water_l;
        sumWaterW += water_l*W[l];
      }
      double soilWaverage = sumWaterW/sumWater;
      double sumDonor = 0.0;
      double sumReceiver = 0.0;
      for(int l=0;l<nlayers;l++) {
        WDiff[l] = W[l] - soilWaverage;
        DonorDiff[l] = std::max(0.0, WDiff[l]);
        ReceiverDiff[l] = std::max(0.0, -WDiff[l]);
        sumDonor += DonorDiff[l];
        sumReceiver += ReceiverDiff[l];
      }
      if(sumDonor>0.0) {
        for(int c=0;c<numCohorts;c++) {
          double redAmount = Eplant[c]*hydraulicRedistributionFraction;
          double hrd = 0.0;
          for(int l=0;l<nlayers;l++) {
            if(WDiff[l]>0.0) {
              hrd += redAmount*DonorDiff[l]/sumDonor;
            } else{
              hrd += -redAmount*ReceiverDiff[l]/sumReceiver;
            }
            outputExtraction(c,l) += hrd;
          }
        }
      }
      delete[] W;
    } else {
      for(int c=0;c<numCohorts;c++) {
        arma::mat& ExtractionPoolsCoh = outputExtractionPools[c]; //this is used to store extraction of a SINGLE plant cohort from all pools
        for(int j=0;j<numCohorts;j++) {
          double sumWater = 0.0;
          double sumWaterW = 0.0;
          for(int l=0;l<nlayers;l++) {
            sumWater += WaterM(j,l);
            sumWaterW += WaterM(j,l)*Wpool(j,l);
          }
          double soilWaverage = sumWaterW/sumWater;
          double sumDonor = 0.0;
          double sumReceiver = 0.0;
          for(int l=0;l<nlayers;l++) {
            WDiff[l] = Wpool(j,l) - soilWaverage;
            DonorDiff[l] = std::max(0.0, WDiff[l]);
            ReceiverDiff[l] = std::max(0.0, -WDiff[l]);
            sumDonor += DonorDiff[l];
            sumReceiver += ReceiverDiff[l];
          }
          double redAmount = 0.0;
          for(int l=0;l<nlayers;l++) redAmount += ExtractionPoolsCoh(j,l)*hydraulicRedistributionFraction;
          if(sumDonor>0.0) {
            double hd = 0.0;
            for(int l=0;l<nlayers;l++) {
              if(WDiff[l]>0.0) {
                hd = redAmount*DonorDiff[l]/sumDonor;
              } else{
                hd = -redAmount*ReceiverDiff[l]/sumReceiver;
              }
              outputExtraction(c,l) += hd;
              ExtractionPoolsCoh(j,l) += hd;
            }
          }
        }
      }
    }
    delete[] WDiff;
    delete[] DonorDiff;
    delete[] ReceiverDiff;
  }


  // Copy output stand
  outputStand.LAI = LAIcell;
  outputStand.LAIlive = LAIcelllive;
  outputStand.LAIexpanded = LAIcellexpanded;
  outputStand.LAIdead = LAIcelldead;

  // Copy output plants
  for(int c =0;c<numCohorts;c++) {
    outputPlants.LAI[c] = LAIphe[c];
    outputPlants.LAIlive[c] = LAIlive[c];
    outputPlants.FPAR[c] = PARcohort[c];
    outputPlants.AbsorbedSWRFraction[c] = 100.0*CohASWRF[c];
    outputPlants.PlantPsi[c] = PlantPsi[c];
    outputPlants.LeafPLC[c] = LeafPLC[c];
    outputPlants.StemPLC[c] = StemPLC[c];
  }
}

Rcpp::DataFrame copyPlantBasicTranspirationResult_c(const PlantsBasicTranspiration_RESULT& plants, ModelInput& x){
  DataFrame plantsDF = DataFrame::create(
    _["LAI"] = Rcpp::wrap(plants.LAI),
    _["LAIlive"] = Rcpp::wrap(plants.LAIlive),
    _["FPAR"] = Rcpp::wrap(plants.FPAR),
    _["AbsorbedSWRFraction"] = Rcpp::wrap(plants.AbsorbedSWRFraction),
    _["Extraction"] = Rcpp::wrap(plants.Extraction),
    _["Transpiration"] = Rcpp::wrap(plants.Transpiration),
    _["GrossPhotosynthesis"] = Rcpp::wrap(plants.GrossPhotosynthesis),
    _["PlantPsi"] = Rcpp::wrap(plants.PlantPsi),
    _["DDS"] = Rcpp::wrap(plants.DDS),
    _["StemRWC"] = Rcpp::wrap(plants.StemRWC),
    _["LeafRWC"] = Rcpp::wrap(plants.LeafRWC),
    _["LFMC"] = Rcpp::wrap(plants.LFMC),
    _["StemPLC"] = Rcpp::wrap(plants.StemPLC),
    _["LeafPLC"] = Rcpp::wrap(plants.LeafPLC),
    _["WaterBalance"] = Rcpp::wrap(plants.WaterBalance)
  );
  plantsDF.attr("row.names") = x.cohorts.CohortCode;
  return(plantsDF);
}

Rcpp::List copyBasicTranspirationResult_c(const BasicTranspiration_RESULT& btc, ModelInput& x) {
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
  
  List l = List::create(_["cohorts"] = copyCohorts_c(x.cohorts),
                        _["Stand"] = standVEC,
                        _["Plants"] = copyPlantBasicTranspirationResult_c(btc.plants, x),
                        _["Extraction"] = Extraction,
                        _["ExtractionPools"] = ExtractionPools);
  return(l);
}
