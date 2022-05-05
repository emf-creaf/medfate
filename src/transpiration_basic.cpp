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

const double WeibullShape=3.0;


//Plant volume in l·m-2 ground = mm
double plantVol(double plantPsi, NumericVector pars) {
  
  double leafrwc = tissueRelativeWaterContent(plantPsi, pars["leafpi0"], pars["leafeps"], 
                                              plantPsi, WeibullShape, pars["psi_critic"], 
                                                                          pars["leafaf"], 0.0);
  double stemrwc = tissueRelativeWaterContent(plantPsi, pars["stempi0"], pars["stemeps"], 
                                              plantPsi, WeibullShape, pars["psi_critic"], 
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
  
  //Soil water at field capacity
  List soil = x["soil"];
  NumericVector Water_FC = waterFC(soil, soilFunctions);
  
  //Meteo input
  double tday = meteovec["tday"];
  double pet = meteovec["pet"];
  double rhmax = meteovec["rhmax"];
  double rhmin = meteovec["rhmin"];
  double tmax = meteovec["tmax"];
  double tmin = meteovec["tmin"];
  //Daily average water vapor pressure at the atmosphere (kPa)
  double vpatm = meteoland::utils_averageDailyVP(tmin, tmax, rhmin,rhmax);
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
  NumericVector Psi_Critic = Rcpp::as<Rcpp::NumericVector>(paramsTransp["Psi_Critic"]);
  NumericVector WUE = Rcpp::as<Rcpp::NumericVector>(paramsTransp["WUE"]);
  NumericVector WUE_decay(numCohorts, 0.2812);
  NumericVector Tmax_LAI(numCohorts, 0.134);
  NumericVector Tmax_LAIsq(numCohorts, -0.006);
  if(paramsTransp.containsElementNamed("Tmax_LAI")) {
    Tmax_LAI = Rcpp::as<Rcpp::NumericVector>(paramsTransp["Tmax_LAI"]);
    Tmax_LAIsq = Rcpp::as<Rcpp::NumericVector>(paramsTransp["Tmax_LAIsq"]);
  }
  if(paramsTransp.containsElementNamed("WUE_decay")) {
    WUE_decay = Rcpp::as<Rcpp::NumericVector>(paramsTransp["WUE_decay"]);
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
  NumericMatrix ExtractionPoolsCoh(numCohorts, nlayers); //this is used to store extraction of a SINGLE plant cohort from all pools
  
  NumericVector Eplant(numCohorts, 0.0), Agplant(numCohorts, 0.0);
  NumericVector DDS(numCohorts, 0.0), LFMC(numCohorts, 0.0);
  NumericVector PLCm(numCohorts), RWCsm(numCohorts), RWClm(numCohorts),RWCssm(numCohorts), RWClsm(numCohorts);
  NumericVector PWB(numCohorts,0.0);
  NumericVector Kl, epc, Vl;
  
  //Calculate unsaturated conductivity
  NumericVector Kunsat = conductivity(soil);
  //Calculate soil water potential
  NumericVector psiSoil = psi(soil,soilFunctions);

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
    }
  }
  
  for(int c=0;c<numCohorts;c++) {
    NumericVector parsVol = NumericVector::create(_["psi_critic"] = Psi_Critic[c], _["stem_plc"] = StemPLC[c],
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
        Klc[l] = Psi2K(psiSoil[l], Psi_Extract[c], WeibullShape);
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
      rootCrownPsi = averagePsi(psiSoil, V(c,_), WeibullShape, Psi_Extract[c]);
      // Rcout<< c << " : " << rootCrownPsi<<"\n";
    } else {
      NumericMatrix RHOPcoh = Rcpp::as<Rcpp::NumericMatrix>(RHOP[c]);
      NumericMatrix Klc(numCohorts, nlayers);
      NumericMatrix Kunlc(numCohorts, nlayers);
      NumericMatrix RHOPcohV(numCohorts, nlayers);
      for(int j = 0;j<numCohorts;j++) {
        for(int l=0;l<nlayers;l++) {
          RHOPcohV(j,l) = RHOPcoh(j,l)*V(c,l);
          Klc(j,l) = Psi2K(psiSoilM(c,l), Psi_Extract[c], WeibullShape);
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
      rootCrownPsi = averagePsiPool(psiSoilM, RHOPcohV, WeibullShape, Psi_Extract[c]);
      // Rcout<< c << " : " << rootCrownPsi<<"\n";
    }
    
    //Cuticular transpiration    
    double lvp_tmax = leafVapourPressure(tmax,  PlantPsi[c]);
    double lvp_tmin = leafVapourPressure(tmin,  PlantPsi[c]);
    double vpd_tmax = std::max(0.0, lvp_tmax - vpatm);
    double vpd_tmin = std::max(0.0, lvp_tmin - vpatm);
    double E_gmin = Gswmin[c]*(vpd_tmin+vpd_tmax)/(2.0*Patm); // mol·s-1·m-2
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
          if(modifyInput) Wpool(j,l) = Wpool(j,l) - (ExtractionPoolsCoh(j,l)/(Water_FC[l]*poolProportions[j])); //Apply extraction from pools
        }
        //Recalculate extraction from soil layers
        Extraction(c,l) = sum(ExtractionPoolsCoh(_,l)); // Sum extraction from all pools (layer l)
      }
    }
    
    //Photosynthesis
    Agplant[c] = WUE[c]*Eplant[c]*std::min(1.0, pow(PARcohort[c]/100.0,WUE_decay[c]));
  }
  
  //Plant water status (StemPLC, RWC, DDS)
  for(int c=0;c<numCohorts;c++) {
    if(cavitationRefill!="total") {
      StemPLC[c] = std::max(1.0 - Psi2K(PlantPsi[c],Psi_Critic[c],WeibullShape), StemPLC[c]); //Track current embolism if no refill
    } else {
      StemPLC[c] = 1.0 - Psi2K(PlantPsi[c],Psi_Critic[c],WeibullShape);
    }
    
    //Relative water content and fuel moisture from plant water potential
    RWClm[c] =  tissueRelativeWaterContent(PlantPsi[c], LeafPI0[c], LeafEPS[c], 
                                           PlantPsi[c], WeibullShape, Psi_Critic[c], 
                                           LeafAF[c], StemPLC[c]);
    RWCsm[c] =  tissueRelativeWaterContent(PlantPsi[c], StemPI0[c], StemEPS[c], 
                                           PlantPsi[c], WeibullShape, Psi_Critic[c], 
                                           StemAF[c], StemPLC[c]);
    LFMC[c] = maxFMC[c]*((1.0/r635[c])*RWClm[c]+(1.0 - (1.0/r635[c]))*RWCsm[c]);
    
    //Daily drought stress from plant WP
    DDS[c] = Phe[c]*(1.0 - Psi2K(PlantPsi[c],Psi_Extract[c],WeibullShape)); 
    
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
  //Modifies input soil
  if(modifyInput) {
    NumericVector Ws = soil["W"];
    for(int l=0;l<nlayers;l++) Ws[l] = Ws[l] - (sum(Extraction(_,l))/Water_FC[l]); 
    if(!plantWaterPools){ //copy soil to the pools of all cohorts
      for(int j=0;j<numCohorts;j++) {
        for(int l=0;l<nlayers;l++) {
          Wpool(j,l) = Ws[l];
        }
      }
    }
  }
  
  //Copy LAIexpanded for output
  NumericVector LAIcohort(numCohorts);
  for(int c=0;c<numCohorts;c++) LAIcohort[c]= LAIphe[c];
  
  NumericVector Stand = NumericVector::create(_["LAI"] = LAIcell,
                                              _["LAIlive"] = LAIcelllive, 
                                              _["LAIexpanded"] = LAIcellexpanded, 
                                              _["LAIdead"] = LAIcelldead);
  
  DataFrame Plants = DataFrame::create(_["LAI"] = LAIcohort,
                                       _["LAIlive"] = clone(LAIlive),
                                       _["FPAR"] = PARcohort,
                                       _["AbsorbedSWRFraction"] = 100.0*CohASWRF, 
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

// [[Rcpp::export("transp_transpirationGranier")]]
List transpirationGranier(List x, DataFrame meteo, int day,
                          double elevation, bool modifyInput = true) {
  if(!meteo.containsElementNamed("MinTemperature")) stop("Please include variable 'MinTemperature' in weather input.");
  NumericVector MinTemperature = meteo["MinTemperature"];
  if(!meteo.containsElementNamed("MaxTemperature")) stop("Please include variable 'MaxTemperature' in weather input.");
  NumericVector MaxTemperature = meteo["MaxTemperature"];
  if(!meteo.containsElementNamed("MinRelativeHumidity")) stop("Please include variable 'MinRelativeHumidity' in weather input.");
  NumericVector MinRelativeHumidity = meteo["MinRelativeHumidity"];
  if(!meteo.containsElementNamed("MaxRelativeHumidity")) stop("Please include variable 'MaxRelativeHumidity' in weather input.");
  NumericVector MaxRelativeHumidity = meteo["MaxRelativeHumidity"];
  if(!meteo.containsElementNamed("MeanTemperature")) stop("Please include variable 'MeanTemperature' in weather input.");
  NumericVector MeanTemperature = meteo["MeanTemperature"];
  if(!meteo.containsElementNamed("PET")) stop("Please include variable 'PET' in weather input.");
  NumericVector PET = meteo["PET"];
  double pet = PET[day-1];
  double tday = MeanTemperature[day-1];
  double tmin = MinTemperature[day-1];
  double tmax = MaxTemperature[day-1];
  double rhmax = MaxRelativeHumidity[day-1];
  double rhmin = MinRelativeHumidity[day-1];
  NumericVector meteovec = NumericVector::create(
    Named("tmax") = tmax,
    Named("tmin") = tmin,
    Named("rhmin") = rhmin, 
    Named("rhmax") = rhmax,
    Named("tday") = tday, 
    Named("pet") = pet);
  return(transpirationGranier(x, meteovec, elevation, modifyInput));
} 


