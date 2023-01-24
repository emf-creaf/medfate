#define STRICT_R_HEADERS
#include <Rcpp.h>
#include <numeric>
#include "lightextinction.h"
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
#include <meteoland.h>
using namespace Rcpp;


//Plant volume in l·m-2 ground = mm
double plantVol(double plantPsi, NumericVector pars) {
  
  double leafrwc = tissueRelativeWaterContent(plantPsi, pars["leafpi0"], pars["leafeps"], 
                                              plantPsi, pars["stem_c"], pars["stem_d"], 
                                              pars["leafaf"], 0.0);
  double stemrwc = tissueRelativeWaterContent(plantPsi, pars["stempi0"], pars["stemeps"], 
                                              plantPsi, pars["stem_c"], pars["stem_d"], 
                                              pars["stemaf"], pars["stem_plc"]);
  return(((pars["Vleaf"] * leafrwc)*pars["LAIphe"]) + ((pars["Vsapwood"] * stemrwc)*pars["LAIlive"]));
}


double findNewPlantPsiConnected(double flowFromRoots, double plantPsi, double rootCrownPsi,
                                NumericVector parsVol){
  //More negative rootCrownPsi causes increased flow due to water being removed
  if(rootCrownPsi <= plantPsi) return(rootCrownPsi); 
  else {
    double V = plantVol(plantPsi, parsVol);
    double psiStep = rootCrownPsi - plantPsi;
    double Vnew = plantVol(plantPsi + psiStep, parsVol);
    while((Vnew - V) > flowFromRoots) {
      psiStep = psiStep/2.0;
      Vnew = plantVol(plantPsi + psiStep, parsVol);
    }
    return(plantPsi+psiStep);
  }
  return(plantPsi);
}

double findNewPlantPsiCuticular(double E_cut, double plantPsi, NumericVector parsVol){
  double V = plantVol(plantPsi, parsVol);
  // Rcout<< " V: "<< V<<" "<< E_cut<<" ";
  double psiStep = 0.001;
  double psi = plantPsi - psiStep;
  double Vnew = plantVol(psi, parsVol);
  double Vdecrease = V - Vnew;
  double Etol = 1e-6;
  int cnt = 0;
  while((std::abs(Vdecrease - E_cut) > Etol) && (cnt < 100)) {
    cnt++;
    if(Vdecrease > E_cut) {
      //Go one step behind and reduce step size 
      psi = psi + psiStep;
      psiStep = psiStep/2.0; 
    } else {
      psi = psi - psiStep;
    }
    Vnew = plantVol(psi, parsVol);
    Vdecrease = V - Vnew;
  }
  return(psi);
}

List transpirationGranier(List x, NumericVector meteovec,  
                          double elevation, bool modifyInput = true) {
  //Control parameters
  List control = x["control"];
  String cavitationRefill = control["cavitationRefill"];
  double refillMaximumRate = control["refillMaximumRate"];
  String soilFunctions = control["soilFunctions"];
  double verticalLayerSize = control["verticalLayerSize"];
  String rhizosphereOverlap = control["rhizosphereOverlap"];
  bool plantWaterPools = (rhizosphereOverlap!="total");
  double hydraulicRedistributionFraction = control["hydraulicRedistributionFraction"];
  
  //Soil water at field capacity
  List soil = x["soil"];
  NumericVector Water_FC = waterFC(soil, soilFunctions);
  
  //Meteo input
  double pet = meteovec["pet"];
  double rhmax = meteovec["rhmax"];
  double rhmin = meteovec["rhmin"];
  double tmax = meteovec["tmax"];
  double tmin = meteovec["tmin"];
  double Catm = meteovec["Catm"];
  //Daily average water vapor pressure at the atmosphere (kPa)
  double vpatm = meteoland::utils_averageDailyVP(tmin, tmax, rhmin, rhmax);
  double vpd = std::max(0.0, meteoland::utils_saturationVP((tmin+tmax)/2.0) - vpatm);
    
  //Atmospheric pressure (kPa)
  double Patm = meteoland::utils_atmosphericPressure(elevation);
  
  //Vegetation input
  DataFrame cohorts = Rcpp::as<Rcpp::DataFrame>(x["cohorts"]);
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
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
  
  //Parameters  
  DataFrame paramsAnatomy = Rcpp::as<Rcpp::DataFrame>(x["paramsAnatomy"]);
  NumericVector Al2As = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["Al2As"]);
  NumericVector r635 = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["r635"]);
  
  DataFrame paramsInterception = Rcpp::as<Rcpp::DataFrame>(x["paramsInterception"]);
  NumericVector kPAR = Rcpp::as<Rcpp::NumericVector>(paramsInterception["kPAR"]);
  
  DataFrame paramsTransp = Rcpp::as<Rcpp::DataFrame>(x["paramsTranspiration"]);
  NumericVector Gswmin = Rcpp::as<Rcpp::NumericVector>(paramsTransp["Gswmin"]);
  NumericVector Psi_Extract = Rcpp::as<Rcpp::NumericVector>(paramsTransp["Psi_Extract"]);
  NumericVector Exp_Extract = Rcpp::as<Rcpp::NumericVector>(paramsTransp["Exp_Extract"]);
  NumericVector VCstem_c = Rcpp::as<Rcpp::NumericVector>(paramsTransp["VCstem_c"]);
  NumericVector VCstem_d = Rcpp::as<Rcpp::NumericVector>(paramsTransp["VCstem_d"]);
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
  DataFrame internalWater = Rcpp::as<Rcpp::DataFrame>(x["internalWater"]);
  NumericVector PlantPsi = clone(Rcpp::as<Rcpp::NumericVector>(internalWater["PlantPsi"]));
  NumericVector StemPLC = clone(Rcpp::as<Rcpp::NumericVector>(internalWater["StemPLC"]));
  
  //Determine whether leaves are out (phenology) and the adjusted Leaf area
  NumericVector Phe(numCohorts,0.0);
  double s = 0.0, LAIcell = 0.0, canopyHeight = 0.0, LAIcelllive = 0.0, LAIcellexpanded = 0.0,LAIcelldead = 0.0;
  for(int c=0;c<numCohorts;c++) {
    if(LAIlive[c]>0) Phe[c]=LAIphe[c]/LAIlive[c]; //Phenological status
    else Phe[c]=0.0;
    s += (kPAR[c]*(LAIphe[c]+LAIdead[c]));
    LAIcell += LAIphe[c]+LAIdead[c];
    LAIcelldead += LAIdead[c];
    LAIcellexpanded +=LAIphe[c];
    LAIcelllive += LAIlive[c];
    if(canopyHeight<H[c]) canopyHeight = H[c];
  }
  int nz = ceil(canopyHeight/verticalLayerSize); //Number of vertical layers
  NumericVector z(nz+1,0.0);
  NumericVector zmid(nz);
  for(int i=1;i<=nz;i++) {
    z[i] = z[i-1] + verticalLayerSize;
    zmid[i-1] = (verticalLayerSize/2.0) + verticalLayerSize*((double) (i-1));
  }
  
  NumericVector PARcohort = parcohortC(H, LAIphe,  LAIdead, kPAR, CR);
  NumericVector CohASWRF = cohortAbsorbedSWRFraction(z, LAIphe,  LAIdead, H, CR, kPAR);
  CohASWRF = pow(CohASWRF, 0.75);
  
  //Apply fractions to potential evapotranspiration
  //Maximum canopy transpiration
  //    Tmax = PET[i]*(-0.006*pow(LAIcell[i],2.0)+0.134*LAIcell[i]+0.036); //From Granier (1999)
  NumericVector Tmax = pet*(Tmax_LAIsq*pow(LAIcell,2.0)+ Tmax_LAI*LAIcell); //From Granier (1999)
  
  //Fraction of Tmax attributed to each plant cohort
  double pabs = std::accumulate(CohASWRF.begin(),CohASWRF.end(),0.0);
  NumericVector TmaxCoh(numCohorts,0.0);
  if(pabs>0.0) TmaxCoh = Tmax*(CohASWRF/pabs);
  
  //Actual plant transpiration
  int nlayers = Wpool.ncol();
  NumericMatrix Extraction(numCohorts, nlayers); // this is final extraction of each cohort from each layer
  List ExtractionPools(numCohorts);
  
  NumericVector Eplant(numCohorts, 0.0), Agplant(numCohorts, 0.0);
  NumericVector DDS(numCohorts, 0.0), LFMC(numCohorts, 0.0);
  NumericVector PLCm(numCohorts), RWCsm(numCohorts), RWClm(numCohorts),RWCssm(numCohorts), RWClsm(numCohorts);
  NumericVector PWB(numCohorts,0.0);
  NumericVector Kl, epc, Vl;
  
  //Calculate unsaturated conductivity
  NumericVector Kunsat = conductivity(soil);
  //Calculate soil water potential
  NumericVector psiSoil = psi(soil,soilFunctions);

  NumericMatrix WaterM(numCohorts, nlayers);
  NumericMatrix KunsatM(numCohorts, nlayers);
  NumericMatrix psiSoilM(numCohorts, nlayers);
  if(plantWaterPools) {
    List soil_pool = clone(soil);
    NumericVector Ws_pool = soil_pool["W"];
    for(int j = 0; j<numCohorts;j++) {
      //Copy values of soil moisture from pool of cohort j
      for(int l = 0; l<nlayers;l++) Ws_pool[l] = Wpool(j,l);
      //Calculate unsaturated conductivity
      KunsatM(j,_) = conductivity(soil_pool);
      //Calculate soil water potential
      psiSoilM(j,_) = psi(soil_pool, soilFunctions);
      //Calculate available water
      WaterM(j,_) = water(soil_pool, soilFunctions);
    }
  }
  
  for(int c=0;c<numCohorts;c++) {
    NumericMatrix ExtractionPoolsCoh(numCohorts, nlayers); //this is used to store extraction of a SINGLE plant cohort from all pools
    
    NumericVector parsVol = NumericVector::create(_["stem_c"] = VCstem_c[c], _["stem_d"] = VCstem_d[c], _["stem_plc"] = StemPLC[c],
                                                  _["leafpi0"] = LeafPI0[c], _["leafeps"] = LeafEPS[c],
                                                  _["leafaf"] = LeafAF[c],_["stempi0"] = StemPI0[c],_["stemeps"] = StemEPS[c],
                                                  _["stemaf"] = StemAF[c],_["Vsapwood"] = Vsapwood[c],_["Vleaf"] = Vleaf[c],
                                                  _["LAIphe"] = LAIphe[c],_["LAIlive"] = LAIlive[c]);
    
    double rootCrownPsi = NA_REAL;
    
    //Extraction from soil (can later be modified if there are changes in plant water content)
    if(!plantWaterPools) {
      NumericVector Klc(nlayers);
      NumericVector Kunlc(nlayers);
      for(int l=0;l<nlayers;l++) {
        Klc[l] = Psi2K(psiSoil[l], Psi_Extract[c], Exp_Extract[c]);
        //Limit Mean Kl due to previous cavitation
        if(cavitationRefill!="total") {
          Klc[l] = std::min(Klc[l], 1.0-StemPLC[c]); 
        }
        Kunlc[l] = pow(Kunsat[l],0.5)*V(c,l);
      }
      double sumKunlc = sum(Kunlc);
      double Klcmean = sum(Klc*V(c,_));
      for(int l=0;l<nlayers;l++) {
        Extraction(c,l) = std::max(TmaxCoh[c]*Klcmean*(Kunlc[l]/sumKunlc),0.0);
      }
      rootCrownPsi = averagePsi(psiSoil, V(c,_), Exp_Extract[c], Psi_Extract[c]);
      // Rcout<< c << " : " << rootCrownPsi<<"\n";
    } else {
      NumericMatrix RHOPcoh = Rcpp::as<Rcpp::NumericMatrix>(RHOP[c]);
      NumericMatrix Klc(numCohorts, nlayers);
      NumericMatrix Kunlc(numCohorts, nlayers);
      NumericMatrix RHOPcohV(numCohorts, nlayers);
      for(int j = 0;j<numCohorts;j++) {
        for(int l=0;l<nlayers;l++) {
          RHOPcohV(j,l) = RHOPcoh(j,l)*V(c,l);
          Klc(j,l) = Psi2K(psiSoilM(c,l), Psi_Extract[c], Exp_Extract[c]);
          //Limit Mean Kl due to previous cavitation
          if(cavitationRefill!="total") Klc(j,l) = std::min(Klc(j,l), 1.0-StemPLC[c]); 
          Kunlc(j,l) = pow(KunsatM(j,l),0.5)*RHOPcohV(j,l);
        }
      }
      double sumKunlc = sum(Kunlc);
      double Klcmean = sum(Klc*RHOPcohV);
      for(int l=0;l<nlayers;l++) {
        for(int j = 0;j<numCohorts;j++) {
          ExtractionPoolsCoh(j,l) = std::max(TmaxCoh[c]*Klcmean*(Kunlc(j,l)/sumKunlc),0.0);
        }
        Extraction(c,l) = sum(ExtractionPoolsCoh(_,l)); // Sum extraction from all pools (layer l)
      }
      rootCrownPsi = averagePsiPool(psiSoilM, RHOPcohV, Exp_Extract[c], Psi_Extract[c]);
      // Rcout<< c << " : " << rootCrownPsi<<"\n";
    }
    
    //Cuticular transpiration    
    double lvp_tmax = leafVapourPressure(tmax,  PlantPsi[c]);
    double lvp_tmin = leafVapourPressure(tmin,  PlantPsi[c]);
    double lvpd_tmax = std::max(0.0, lvp_tmax - vpatm);
    double lvpd_tmin = std::max(0.0, lvp_tmin - vpatm);
    double E_gmin = Gswmin[c]*(lvpd_tmin+lvpd_tmax)/(2.0*Patm); // mol·s-1·m-2
    double E_cut = E_gmin*LAIphe[c]*(24.0*3600.0*0.018);

    double oldVol = plantVol(PlantPsi[c], parsVol); 
    
    //Transpiration is the maximum of predicted extraction and cuticular transpiration
    double ext_sum = sum(Extraction(c,_));
    double corr_extraction = 0.0;
    if(E_cut > ext_sum) {
      Eplant[c] = E_cut;
      corr_extraction = E_cut - ext_sum; //Correction to be added to extraction
    } else {
      Eplant[c] = ext_sum;
    }
    PlantPsi[c] = findNewPlantPsiConnected(Eplant[c], PlantPsi[c], rootCrownPsi, parsVol);
    //For deciduous species, make water potential follow soil during winter
    if(LAIphe[c]==0.0) PlantPsi[c] = rootCrownPsi;
    double newVol = plantVol(PlantPsi[c], parsVol);
    
    double volDiff = newVol - oldVol;
    //Plant transpiration and water balance
    PWB[c] = volDiff;
    
    //Divide the difference among soil layers extraction
    for(int l=0;l<nlayers;l++) {
      if(!plantWaterPools) { 
        if(ext_sum>0.0) Extraction(c,l) += (volDiff+corr_extraction)*(Extraction(c,l)/ext_sum);
      } else { // recalculate also extraction from soil pools
        for(int j = 0;j<numCohorts;j++) {
          if(ext_sum>0.0) ExtractionPoolsCoh(j,l) += (volDiff+corr_extraction)*(ExtractionPoolsCoh(j,l)/ext_sum);
        }
        //Recalculate extraction from soil layers
        Extraction(c,l) = sum(ExtractionPoolsCoh(_,l)); // Sum extraction from all pools (layer l)
      }
    }
    //Photosynthesis
    double fpar = std::min(1.0, pow(PARcohort[c]/100.0,WUE_par[c]));
    double fco2 = (1.0 - exp((-1.0)*WUE_co2[c]*Catm));
    double fvpd = pow(vpd, WUE_vpd[c]);
    if(vpd < 0.25) fvpd = 2.5 - (2.5 - pow(0.25, WUE_vpd[c]))*(vpd/0.25);
    // Rcout<<fpar<<" "<< fco2 << " "<< fvpd<< " "<< WUE[c]*fpar*fco2*fvpd<<"\n";
    Agplant[c] = WUE[c]*Eplant[c]*fpar*fco2*fvpd;
    
    //Store extractionPool
    ExtractionPools[c] = ExtractionPoolsCoh;
  }
  
  //Plant water status (StemPLC, RWC, DDS)
  for(int c=0;c<numCohorts;c++) {
    if(cavitationRefill!="total") {
      StemPLC[c] = std::max(1.0 - xylemConductance(PlantPsi[c], 1.0, VCstem_c[c], VCstem_d[c]), StemPLC[c]); //Track current embolism if no refill
    } else {
      StemPLC[c] = 1.0 - xylemConductance(PlantPsi[c], 1.0, VCstem_c[c], VCstem_d[c]);
    }
    
    //Relative water content and fuel moisture from plant water potential
    RWClm[c] =  tissueRelativeWaterContent(PlantPsi[c], LeafPI0[c], LeafEPS[c], 
                                           PlantPsi[c], VCstem_c[c], VCstem_d[c], 
                                           LeafAF[c], StemPLC[c]);
    RWCsm[c] =  tissueRelativeWaterContent(PlantPsi[c], StemPI0[c], StemEPS[c], 
                                           PlantPsi[c], VCstem_c[c], VCstem_d[c], 
                                           StemAF[c], StemPLC[c]);
    LFMC[c] = maxFMC[c]*((1.0/r635[c])*RWClm[c]+(1.0 - (1.0/r635[c]))*RWCsm[c]);
    
    //Daily drought stress from plant WP
    DDS[c] = Phe[c]*(1.0 - Psi2K(PlantPsi[c],Psi_Extract[c],Exp_Extract[c])); 
    
    if(cavitationRefill=="rate") {
      double SAmax = 10e4/Al2As[c]; //cm2·m-2 of leaf area
      double r = refillMaximumRate*std::max(0.0, (PlantPsi[c] + 1.5)/1.5);
      StemPLC[c] = std::max(0.0, StemPLC[c] - (r/SAmax));
    }
  }
  
  
  if(modifyInput) {
    internalWater["StemPLC"] = StemPLC;
    internalWater["PlantPsi"] = PlantPsi;
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
            Extraction(c,l) += HD[l];
          }
        }
        // Rcout<< "\n";
      }
    } else {
      for(int c=0;c<numCohorts;c++) {
        NumericMatrix ExtractionPoolsCoh = ExtractionPools[c];
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
              Extraction(c,l) += HD[l];
              ExtractionPoolsCoh(j,l) += HD[l];
            }
          }
        }
      }
    }
  }
  //Modifies input soil
  if(modifyInput) {
    NumericVector Ws = soil["W"];
    for(int l=0;l<nlayers;l++) Ws[l] = Ws[l] - (sum(Extraction(_,l))/Water_FC[l]); 
    for(int c=0;c<numCohorts;c++) {
      for(int l=0;l<nlayers;l++) {
        if(!plantWaterPools){ //copy soil to the pools of all cohorts
          Wpool(c,l) = Ws[l];
        } else {//Applies pool extraction by each plant cohort
          NumericMatrix ExtractionPoolsCoh = ExtractionPools[c];
          for(int j=0;j<numCohorts;j++) {
            Wpool(c,l) = Wpool(c,l) - (ExtractionPoolsCoh(j,l)/(Water_FC[l]*poolProportions[c])); //Apply extraction from pools
          }
        }
      }
    } 
  }
  
  //Copy LAIexpanded for output
  NumericVector LAIcohort(numCohorts);
  NumericVector SoilExtractCoh(numCohorts,0.0);
  for(int c=0;c<numCohorts;c++) {
    LAIcohort[c]= LAIphe[c];
    SoilExtractCoh[c] =  sum(Extraction(c,_));
  }
  NumericVector Stand = NumericVector::create(_["LAI"] = LAIcell,
                                              _["LAIlive"] = LAIcelllive, 
                                              _["LAIexpanded"] = LAIcellexpanded, 
                                              _["LAIdead"] = LAIcelldead);
  
  DataFrame Plants = DataFrame::create(_["LAI"] = LAIcohort,
                                       _["LAIlive"] = clone(LAIlive),
                                       _["FPAR"] = PARcohort,
                                       _["AbsorbedSWRFraction"] = 100.0*CohASWRF, 
                                       _["Extraction"] = SoilExtractCoh,
                                       _["Transpiration"] = Eplant, 
                                       _["GrossPhotosynthesis"] = Agplant,
                                       _["PlantPsi"] = PlantPsi, 
                                       _["DDS"] = DDS,
                                       _["StemRWC"] = RWCsm,
                                       _["LeafRWC"] = RWClm,
                                       _["LFMC"] = LFMC,
                                       _["StemPLC"] = StemPLC,
                                       _["WaterBalance"] = PWB);
  Plants.attr("row.names") = above.attr("row.names");
  Extraction.attr("dimnames") = List::create(above.attr("row.names"), seq(1,nlayers));
  List l = List::create(_["cohorts"] = clone(cohorts),
                        _["Stand"] = Stand,
                        _["Plants"] = Plants,
                        _["Extraction"] = Extraction);
  return(l);
}

//' Transpiration modes
//' 
//' High-level sub-models to represent transpiration, plant hydraulics and water relations 
//' within plants. The two submodels represent a very different degree of complexity, 
//' and correspond to Granier et al. (1999) or Sperry et al. (2017).
//' 
//' @param x An object of class \code{\link{spwbInput}} or \code{\link{growthInput}}, built using the 'Granier' or 'Sperry' transpiration modes, depending on the function to be called.
//' @param meteo A data frame with daily meteorological data series:
//'   \itemize{
//'     \item{\code{DOY}: Day of the year (Julian day).}
//'     \item{\code{Precipitation}: Precipitation (in mm).}
//'     \item{\code{MinTemperature}: Minimum temperature (in degrees Celsius).}
//'     \item{\code{MaxTemperature}: Maximum temperature (in degrees Celsius).}
//'     \item{\code{MinRelativeHumidity}: Minimum relative humidity (in percent).}
//'     \item{\code{MaxRelativeHumidity}: Maximum relative humidity (in percent).}
//'     \item{\code{Radiation}: Solar radiation (in MJ/m2/day).}
//'     \item{\code{WindSpeed}: Wind speed (in m/s). If not available, this column can be left with \code{NA} values.}
//'   }
//' @param day An integer to identify a day within \code{meteo}.
//' @param latitude Latitude (in degrees).
//' @param elevation,slope,aspect Elevation above sea level (in m), slope (in degrees) and aspect (in degrees from North).
//' @param modifyInput Boolean flag to indicate that the input \code{x} object is allowed to be modified during the simulation.
//' 
//' @return
//' Function \code{transp_transpirationGranier} and \code{transp_transpirationSperry} return a list with the following elements:
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
//'       \item{\code{"DDS"}: Daily drought stress [0-1] (relative whole-plant conductance).}
//'     }
//'   When using \code{transp_transpirationSperry}, element \code{"Plants"} includes:
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
//'       \item{\code{"DDS"}: Daily drought stress [0-1] (relative whole-plant conductance).}
//'       \item{\code{"StemRWC"}: Relative water content of stem tissue (including symplasm and apoplasm).}
//'       \item{\code{"LeafRWC"}: Relative water content of leaf tissue (including symplasm and apoplasm).}
//'       \item{\code{"LFMC"}: Live fuel moisture content (in percent of dry weight).}
//'       \item{\code{"WaterBalance"}: Plant water balance (extraction - transpiration).}
//'     }
//'   }
//'   \item{\code{"Extraction"}: A data frame with mm of water extracted from each soil layer (in columns) by each cohort (in rows).}
//' 
//'   The remaining items are only given by \code{transp_transpirationSperry}:
//'   \item{\code{"EnergyBalance"}: When using the 'Sperry' transpiration mode, the model performs energy balance of the stand and 'EnergyBalance' is a list with the following:
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
//'       \item{\code{"StemPLC"}: A data frame with the proportion loss of conductance [0-1] for each plant cohort during each time step. }
//'       \item{\code{"StemRWC"}: A data frame with the (average) relative water content of stem tissue [0-1] for each plant cohort during each time step. }
//'       \item{\code{"LeafRWC"}: A data frame with the relative water content of leaf tissue [0-1] for each plant cohort during each time step. }
//'       \item{\code{"StemSympRWC"}: A data frame with the (average) relative water content of symplastic stem tissue [0-1] for each plant cohort during each time step. }
//'       \item{\code{"LeafSympRWC"}: A data frame with the relative water content of symplastic leaf tissue [0-1] for each plant cohort during each time step. }
//'       \item{\code{"PWB"}: A data frame with plant water balance (extraction - transpiration).}
//'     }
//'   }
//'   \item{\code{"LightExtinction"}: A list of information regarding radiation balance through the canopy, as returned by function \code{\link{light_instantaneousLightExtinctionAbsortion}}.}
//'   \item{\code{"CanopyTurbulence"}: Canopy turbulence (see \code{\link{wind_canopyTurbulence}}).}
//'   \item{\code{"SupplyFunctions"}: If \code{stepFunctions} is not missing, a list of supply functions, photosynthesis functions and profit maximization functions.}
//' }
//' 
//' @references
//' Granier A, \enc{Bréda}{Breda} N, Biron P, Villette S (1999) A lumped water balance model to evaluate duration and intensity of drought constraints in forest stands. Ecol Modell 116:269–283. https://doi.org/10.1016/S0304-3800(98)00205-1.
//' 
//' Sperry, J. S., M. D. Venturas, W. R. L. Anderegg, M. Mencuccini, D. S. Mackay, Y. Wang, and D. M. Love. 2017. Predicting stomatal responses to the environment from the optimization of photosynthetic gain and hydraulic cost. Plant Cell and Environment 40, 816-830 (doi: 10.1111/pce.12852).
//' 
//' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
//' 
//' @seealso \code{\link{spwb_day}}, \code{\link{plot.spwb_day}}
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
//' #Initialize soil with default soil params (4 layers)
//' examplesoil = soil(defaultSoilParams(4))
//' 
//' #Initialize control parameters
//' control = defaultControl("Granier")
//' 
//' #Initialize input
//' x1 = forest2spwbInput(exampleforestMED,examplesoil, SpParamsMED, control)
//' 
//' # Transpiration according to Granier's model, plant water potential 
//' # and plant stress for a given day
//' t1 = transp_transpirationGranier(x1, examplemeteo, 1, 
//'                                  latitude = 41.82592, elevation = 100, slope = 0, aspect = 0, 
//'                                  modifyInput = FALSE)
//' 
//' #Switch to 'Sperry' transpiration mode
//' control = defaultControl("Sperry")
//' 
//' #Initialize input
//' x2 = forest2spwbInput(exampleforestMED,examplesoil, SpParamsMED, control)
//' 
//' # Transpiration according to Sperry's model
//' t2 = transp_transpirationSperry(x2, examplemeteo, 1, 
//'                                 latitude = 41.82592, elevation = 100, slope = 0, aspect = 0,
//'                                 modifyInput = FALSE)
//'                                 
//' @name transp_modes
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
  
  if(NumericVector::is_na(latitude)) stop("Value for 'latitude' should not be missing.");
  double latrad = latitude * (M_PI/180.0);
  if(NumericVector::is_na(aspect)) aspect = 0.0;
  if(NumericVector::is_na(slope)) slope = 0.0;
  double asprad = aspect * (M_PI/180.0);
  double slorad = slope * (M_PI/180.0);
  
  CharacterVector dateStrings = meteo.attr("row.names");
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
    Named("Catm") = Catm);
  return(transpirationGranier(x, meteovec, elevation, modifyInput));
} 


