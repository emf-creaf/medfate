#include <numeric>
#include "lightextinction.h"
#include "windextinction.h"
#include "hydraulics.h"
#include "biophysicsutils.h"
#include "photosynthesis.h"
#include "soil.h"
#include <Rcpp.h>
#include <meteoland.h>
using namespace Rcpp;

// [[Rcpp::export(".er")]]
NumericVector er(IntegerVector DOY, double ERconv=0.05, double ERsyn = 0.2){
  int nDays = DOY.size();
  NumericVector ER=rep(0.0,nDays);
  for(int i=0;i<nDays;i++){
    if((DOY[i]<=120)|(DOY[i]>=335)) {
      ER[i] = ERsyn;
    } else {
      ER[i] = ERconv;
    }
  }
  return(ER);
  
}
// [[Rcpp::export(".gdd")]]
NumericVector gdd(IntegerVector DOY, NumericVector Temp, double Tbase = 5.0){
  int nDays = Temp.size();
  NumericVector GDD(nDays);
  double cum = 0.0;
  for(int i=0;i<nDays;i++){
    if((Temp[i]-Tbase < 0.0) & (DOY[i]>180)) {
      cum = 0.0;
    } else {
      if(Temp[i]-Tbase>0.0) cum = cum + (Temp[i]-Tbase);
    }
    GDD[i] = cum;
    if(DOY[i] >= 365) cum = 0.0;
  }
  return(GDD);
}

// [[Rcpp::export("swb.SoilEvaporation")]]
double soilevaporation(double DEF,double PETs, double Gsoil){
  double t = pow(DEF/Gsoil, 2.0);
  double Esoil = 0.0;
  Esoil = std::min(Gsoil*(sqrt(t+1)-sqrt(t)), PETs);
  return(Esoil);
}
// [[Rcpp::export(".infiltrationDay")]]
double infiltrationDay(double NetPrec, double Ssoil) {
  double I = 0;
  if(NetPrec>0.2*Ssoil) {
    I = NetPrec-(pow(NetPrec-0.2*Ssoil,2.0)/(NetPrec+0.8*Ssoil));
  } else {
    I = NetPrec;
  }
  return(I);
}


// [[Rcpp::export(".interceptionGashDay")]]
double interceptionGashDay(double Precipitation, double Cm, double p, double ER=0.05) {
    double I = 0.0;
    double PG = (-Cm/(ER*(1.0-p)))*log(1.0-ER); //Precipitation need to saturate the canopy
    if(Cm==0.0 || p==1.0) PG = 0.0; //Avoid NAs
    if(Precipitation>PG) {
      I = (1-p)*PG + (1-p)*ER*(Precipitation-PG);
    } else {
      I = (1-p)*Precipitation;
    }
    return(I);
}


// Soil water balance with simple hydraulic model
// [[Rcpp::export(".swbDay1")]]
List swbDay1(List x, List soil, double gdd, double tday, double pet, double rain, double er, double runon=0.0, 
             bool verbose = false) {

  //Soil input
  NumericVector W = soil["W"];
  NumericVector psi = soil["psi"];
  NumericVector dVec = soil["dVec"];
  NumericVector Theta_FC = soil["Theta_FC"];
  NumericVector Water_FC = soil["Water_FC"];
  NumericVector macro = soil["macro"];
  NumericVector clay = soil["clay"];
  NumericVector sand = soil["sand"];
  int nlayers = W.size();

  //Vegetation input
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  NumericVector LAI = Rcpp::as<Rcpp::NumericVector>(above["LAI_live"]);
  NumericVector H = Rcpp::as<Rcpp::NumericVector>(above["H"]);
  NumericVector CR = Rcpp::as<Rcpp::NumericVector>(above["CR"]);
  NumericVector Sgdd = Rcpp::as<Rcpp::NumericVector>(x["Sgdd"]);
  NumericVector kPAR = Rcpp::as<Rcpp::NumericVector>(x["k"]);
  NumericVector gRainIntercept = Rcpp::as<Rcpp::NumericVector>(x["g"]);
  NumericVector transpiration = Rcpp::as<Rcpp::NumericVector>(x["Transpiration"]);
  NumericVector photosynthesis = Rcpp::as<Rcpp::NumericVector>(x["Photosynthesis"]);
  NumericVector Psi_Extract = Rcpp::as<Rcpp::NumericVector>(x["Psi_Extract"]);
  NumericVector WUE = Rcpp::as<Rcpp::NumericVector>(x["WUE"]);
  int numCohorts = LAI.size();

  //Root distribution input
  NumericMatrix V = Rcpp::as<Rcpp::NumericMatrix>(x["V"]);


  //Determine whether leaves are out (phenology) and the adjusted Leaf area
  NumericVector LAIphe(numCohorts);
  NumericVector Phe = pmin(pmax(gdd/Sgdd,0.0),1.0);
  double s = 0.0, LAIcell = 0.0, Cm = 0.0;
  for(int c=0;c<numCohorts;c++) {
    if(Sgdd[c]==0.0) Phe[c]=1.0;
    LAIphe[c] = LAI[c]*Phe[c]; //LAI modified by phenology
    s += (kPAR[c]*LAIphe[c]);
    LAIcell += LAIphe[c];
    Cm += LAIphe[c]*gRainIntercept[c];
  }
  NumericVector CohASWRF = cohortAbsorbedSWRFraction(LAIphe,  H, CR, kPAR);
  NumericVector CohPAR = parcohortC(H, LAIphe, kPAR, CR)/100.0;
  double LgroundPAR = exp((-1)*s);
  double LgroundSWR = 1.0 - std::accumulate(CohASWRF.begin(),CohASWRF.end(),0.0);
  
  //Hydrologic input
  double NetPrec = 0.0, Infiltration= 0.0, Runoff= 0.0, DeepDrainage= 0.0;
  if(rain>0.0) {
    //Interception
    NetPrec = rain - interceptionGashDay(rain,Cm,LgroundPAR,er);
    //Net Runoff and infiltration
    Infiltration = infiltrationDay(NetPrec+runon, soil["Ssoil"]);
    Runoff = (NetPrec+runon) - Infiltration;
    //Input of the first soil layer is infiltration
    double VI = Infiltration;
    double Wn;
    //Update topsoil layer
    for(int l=0;l<nlayers;l++) {
      if((dVec[l]>0) & (VI>0)) {
        Wn = W[l]*Water_FC[l] + VI*(1.0-macro[l]); //Update water volume
        VI = VI*macro[l] + std::max(Wn - Water_FC[l],0.0); //Update VI, adding the excess to the infiltrating water (saturated flow)
        W[l] = std::min(std::max(0.0,Wn/Water_FC[l]),1.0); //Update theta
      } 
      DeepDrainage = VI; //Reset deep drainage
    }
  }
  for(int l=0;l<nlayers;l++) psi[l] = theta2psi(clay[l], sand[l], W[l]*Theta_FC[l]);


  //Proportion of transpiration that absorbed by each plant cohort (old version)
  // NumericVector PP = CohLight*LAIphe;
  // NumericVector f = PP/std::accumulate(PP.begin(),PP.end(),0.0); 
  // if(LAIcell==0.0) f = rep(0.0,numCohorts); //Avoids NaN values

  //Apply fractions to potential evapotranspiration
  //Maximum canopy transpiration
  //    Tmax = PET[i]*(-0.006*pow(LAIcell[i],2.0)+0.134*LAIcell[i]+0.036); //From Granier (1999)
  double Tmax = pet*(-0.006*pow(LAIcell,2.0)+0.134*LAIcell); //From Granier (1999)
  double PETsoil = pet*LgroundSWR;

  //Fraction of Tmax attributed to each plant cohort
  double pabs = std::accumulate(CohASWRF.begin(),CohASWRF.end(),0.0);
  NumericVector TmaxCoh(numCohorts);
  if(pabs>0.0) TmaxCoh = Tmax*(CohASWRF/pabs);
  
  //Actual plant transpiration
  NumericMatrix EplantCoh(numCohorts, nlayers);
  NumericVector PlantPsi(numCohorts, NA_REAL);
  NumericVector Eplant(numCohorts, 0.0);
  NumericVector DDS(numCohorts, 0.0);
  NumericVector EplantVec(nlayers, 0.0);
  NumericVector Kl, epc, Vl;
  double WeibullShape=3.0;
  for(int l=0;l<nlayers;l++) {
    Kl = Psi2K(psi[l], Psi_Extract, WeibullShape);
    Vl = V(_,l);
    epc = pmax(TmaxCoh*Kl*Vl,0.0);
    EplantCoh(_,l) = epc;
    Eplant = Eplant + epc;
    EplantVec[l] = std::accumulate(epc.begin(),epc.end(),0.0);
    DDS = DDS + Phe*(Vl*(1.0 - Kl));
  }
  for(int c=0;c<numCohorts;c++) {
    PlantPsi[c] = averagePsi(psi, V(c,_), WeibullShape, Psi_Extract[c]);
  }

  //Evaporation from bare soil
  double Gsoil = soil["Gsoil"];
  double Ksoil = soil["Ksoil"];
  double Esoil = soilevaporation((Water_FC[0]*(1.0 - W[0])), PETsoil, Gsoil);
  NumericVector EsoilVec(nlayers,0.0);//Exponential decay to divide bare soil evaporation among layers
  for(int l=0;l<nlayers;l++) {
    double cumAnt = 0.0;
    double cumPost = 0.0;
    for(int l2=0;l2<l;l2++) cumAnt +=dVec[l2];
    cumPost = cumAnt+dVec[l];
    if(l<(nlayers-1)) EsoilVec[l] = Esoil*(exp(-Ksoil*cumAnt)-exp(-Ksoil*cumPost));
    else EsoilVec[l] = Esoil*exp(-Ksoil*cumAnt);
    
  }
  // Rcout<<Esoil<<" "<<  std::accumulate(EsoilVec.begin(),EsoilVec.end(),0.0)<<"\n";

  //Apply decrease in soil layers
  for(int l=0;l<nlayers;l++) W[l] =std::max(W[l]-(EplantVec[l]+EsoilVec[l])/Water_FC[l],0.0);


  double alpha = std::max(std::min(tday/20.0,1.0),0.0);
  //For comunication with growth and landscape water balance
  for(int c=0;c<numCohorts;c++) {
    transpiration[c] = Eplant[c];
    photosynthesis[c] = alpha*WUE[c]*transpiration[c];
  }

  List l = List::create(_["NetPrec"] = NetPrec, _["Runon"] = runon, _["Infiltration"] = Infiltration, _["Runoff"] = Runoff, _["DeepDrainage"] = DeepDrainage,
                        _["LAIcell"] = LAIcell, _["Cm"] = Cm, _["Lground"] = LgroundPAR,
                        _["EsoilVec"] = EsoilVec, _["EplantVec"] = EplantVec, _["psiVec"] = psi,
                        _["EplantCoh"] = Eplant, _["psiCoh"] = PlantPsi, _["DDS"] = DDS);
  return(l);
}


// Soil water balance with Sperry hydraulic and stomatal conductance models
// [[Rcpp::export(".swbDay2")]]
List swbDay2(List x, List soil, double gdd, double tmin, double tmax, double rhmin, double rhmax, double rad, double wind, 
             double latitude, double elevation, double slope, double aspect, double delta, 
             double rain, double er, double runon=0.0, bool verbose = false) {
  
  //Soil input
  NumericVector W = soil["W"];
  NumericVector psi = soil["psi"];
  NumericVector dVec = soil["dVec"];
  NumericVector Theta_FC = soil["Theta_FC"];
  NumericVector Water_FC = soil["Water_FC"];
  NumericVector macro = soil["macro"];
  NumericVector clay = soil["clay"];
  NumericVector sand = soil["sand"];
  int nlayers = W.size();
  
  //Vegetation input
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  NumericVector LAI = Rcpp::as<Rcpp::NumericVector>(above["LAI_live"]);
  NumericVector H = Rcpp::as<Rcpp::NumericVector>(above["H"]);
  NumericVector CR = Rcpp::as<Rcpp::NumericVector>(above["CR"]);
  NumericVector Sgdd = Rcpp::as<Rcpp::NumericVector>(x["Sgdd"]);
  NumericVector kPAR = Rcpp::as<Rcpp::NumericVector>(x["k"]);
  NumericVector gRainIntercept = Rcpp::as<Rcpp::NumericVector>(x["g"]);
  NumericVector transpiration = Rcpp::as<Rcpp::NumericVector>(x["Transpiration"]);
  NumericVector photosynthesis = Rcpp::as<Rcpp::NumericVector>(x["Photosynthesis"]);
  NumericVector Gwmax = Rcpp::as<Rcpp::NumericVector>(x["Gwmax"]);
  NumericVector VG_n = Rcpp::as<Rcpp::NumericVector>(soil["VG_n"]);
  NumericVector VG_alpha = Rcpp::as<Rcpp::NumericVector>(soil["VG_alpha"]);
  NumericVector VCstem_kmax = Rcpp::as<Rcpp::NumericVector>(x["VCstem_kmax"]);
  NumericVector VCstem_c = Rcpp::as<Rcpp::NumericVector>(x["VCstem_c"]);
  NumericVector VCstem_d = Rcpp::as<Rcpp::NumericVector>(x["VCstem_d"]);
  NumericVector VCroot_c = x["VCroot_c"];
  NumericVector VCroot_d = x["VCroot_d"];
  NumericVector Vmax298 = x["Vmax298"];
  NumericVector Jmax298 = x["Jmax298"];
  //Root distribution input
  NumericMatrix V =Rcpp::as<Rcpp::NumericMatrix>(x["V"]);
  NumericMatrix VCroot_kmax= Rcpp::as<Rcpp::NumericMatrix>(x["VCroot_kmax"]);
  NumericMatrix VGrhizo_kmax= Rcpp::as<Rcpp::NumericMatrix>(x["VGrhizo_kmax"]);
  int numCohorts = LAI.size();
  
 
  double latrad = latitude * (PI/180.0);
  double asprad = aspect * (PI/180.0);
  double slorad = slope * (PI/180.0);
  //Day length (latitude in radians), daylight temperature and VPD 
  double tauday = meteoland::radiation_daylengthseconds(latrad,  slorad,asprad, delta); 
  double Patm = meteoland::utils_atmosphericPressure(elevation);
  double Catm = 386.0;

 
  //Leaf phenology and the adjusted leaf area index
  NumericVector LAIphe(numCohorts);
  NumericVector Phe = pmin(pmax(gdd/Sgdd,0.0),1.0);
  double s = 0.0, LAIcell = 0.0, Cm = 0.0, canopyHeight = 0.0;
  for(int c=0;c<numCohorts;c++) {
    if(Sgdd[c]==0.0) Phe[c]=1.0;
    LAIphe[c] = LAI[c]*Phe[c]; //LAI modified by phenology
    s += (kPAR[c]*LAIphe[c]);
    LAIcell += LAIphe[c];
    Cm += LAIphe[c]*gRainIntercept[c];
    if(canopyHeight<H[c]) canopyHeight = H[c];
  }
  
  //Light extinction
  NumericVector CohASWRF = cohortAbsorbedSWRFraction(LAIphe,  H, CR, kPAR);
  NumericVector CohPAR = parcohortC(H, LAIphe, kPAR, CR)/100.0;
  double LgroundPAR = exp((-1)*s);
  double LgroundSWR = 1.0 - std::accumulate(CohASWRF.begin(),CohASWRF.end(),0.0);
  
  //Wind extinction
  NumericVector CohWind = windExtinctionCohort(H,CR, wind, LAIcell, canopyHeight);
  
  //Hydrologic input
  double NetPrec = 0.0, Infiltration= 0.0, Runoff= 0.0, DeepDrainage= 0.0;
  if(rain>0.0) {
    //Interception
    NetPrec = rain - interceptionGashDay(rain,Cm,LgroundPAR,er);
    //Net Runoff and infiltration
    Infiltration = infiltrationDay(NetPrec+runon, soil["Ssoil"]);
    Runoff = (NetPrec+runon) - Infiltration;
    //Input of the first soil layer is infiltration
    double VI = Infiltration;
    double Wn;
    //Update topsoil layer
    for(int l=0;l<nlayers;l++) {
      if((dVec[l]>0) & (VI>0)) {
        Wn = W[l]*Water_FC[l] + VI*(1.0-macro[l]); //Update water volume
        VI = VI*macro[l] + std::max(Wn - Water_FC[l],0.0); //Update VI, adding the excess to the infiltrating water (saturated flow)
        W[l] = std::min(std::max(0.0,Wn/Water_FC[l]),1.0); //Update theta
      } 
      DeepDrainage = VI; //Reset deep drainage
    }
  }
  for(int l=0;l<nlayers;l++) {
    psi[l] = theta2psi(clay[l], sand[l], W[l]*Theta_FC[l]);
    // Rcout<<psi[l]<<" ";
  }
  // Rcout<<"\n";
  
  //Radiation balance (latitude in radians)
  double vpa = meteoland::utils_averageDailyVP(tmin, tmax, rhmin,rhmax);
  double Rn = meteoland::radiation_netRadiation(latrad, elevation, slorad, asprad, delta, 
                                                vpa, tmin, tmax, rad);
  // Rcout<<Rn<<" "<<Rninst<<" "<<Iparinst<<" Q:"<<Q<<"\n";
  
  //Transpiration
  NumericMatrix EplantCoh(numCohorts, nlayers);
  NumericVector PlantPsi(numCohorts);
  NumericMatrix K(numCohorts, nlayers);
  NumericVector Eplant(numCohorts);
  NumericVector EplantVec(nlayers);
  NumericVector DDS(numCohorts);
  NumericVector Vc, VCroot_kmaxc, VGrhizo_kmaxc, psic, VG_nc,VG_alphac;
  for(int c=0;c<numCohorts;c++) {
    int nlayersc = 0;
    while((nlayersc<nlayers) & (V(c,nlayersc)>0)) {
      nlayersc++;
    }
    Vc = NumericVector(nlayersc);
    VCroot_kmaxc = NumericVector(nlayersc);
    VGrhizo_kmaxc = NumericVector(nlayersc);
    psic = NumericVector(nlayersc);
    VG_nc = NumericVector(nlayersc);
    VG_alphac= NumericVector(nlayersc);
    for(int l=0;l<nlayersc;l++) {
      Vc[l] = V(c,l);
      VCroot_kmaxc[l] = VCroot_kmax(c,l);
      VGrhizo_kmaxc[l] = VGrhizo_kmax(c,l);
      psic[l] = psi[l];
      VG_nc[l] = VG_n[l];
      VG_alphac[l] = VG_alpha[l];
    }

    
    List supplyNetwork = supplyFunctionNetwork(psic,
                                               VGrhizo_kmaxc,VG_nc,VG_alphac,
                                               VCroot_kmaxc, VCroot_c[c],VCroot_d[c],
                                               VCstem_kmax[c], VCstem_c[c],VCstem_d[c]);
    NumericVector Etot = supplyNetwork["E"];
    NumericMatrix ElayersMat = supplyNetwork["Elayers"];

    double Tair, propRad;

    double minPsi = 0.0;
    int ntimesteps = 10;
    double tstep = tauday/((double) ntimesteps);
    double t = tstep/2.0;
    NumericVector Ec(nlayers,0.0);
    photosynthesis[c] = 0.0;
    double Gwmin = 0.00001;
    for(int n=0;n<ntimesteps;n++) {
      Tair = temperatureDiurnalPattern(t, tmin, tmax, tauday);
      propRad = std::max(0.0,radiationDiurnalPattern(t, tauday));
      double Rninst = pow(10.0, 6.0)*(Rn*propRad)*CohASWRF[c]/LAIphe[c]; //Instantaneous radiation absorbed
      double Iparinst = pow(10.0, 6.0)*(rad*0.5*propRad)*CohPAR[c]; //Instantaneous incident PAR
      double Qinst = irradianceToPhotonFlux(Iparinst);
      

      List photo = photosynthesisFunction(supplyNetwork,
                                          Catm, Patm,
                                          Tair, vpa,
                                          CohWind[c], Rninst, Qinst, Vmax298[c], Jmax298[c]);
      List PM =  profitMaximization(supplyNetwork, photo, Gwmin, Gwmax[c]);
      int iPM = PM["iMaxProfit"];
      NumericVector PsiVec = supplyNetwork["PsiPlant"];
      NumericVector Ect = ElayersMat(iPM,_);
      for(int l=0;l<nlayersc;l++) Ec[l] += Ect[l]*0.001*0.01802*LAIphe[c]*tstep;
      NumericVector An = photo["NetPhotosynthesis"];
      photosynthesis[c] +=pow(10,-6)*12.01017*An[iPM]*LAIphe[c]*tstep;
      minPsi = std::min(minPsi,PsiVec[iPM]);
      // Rcout<<"[T"<<n<<" Rn: "<<Rninst<<" Qinst: "<< Qinst <<" iPM: "<<iPM<<" E: "<< sum(Ect)<<"]";
      t +=tstep;
    }
    // Rcout<<"\n";
    for(int l=0;l<nlayersc;l++) EplantCoh(c,l) =Ec[l];
    if(nlayersc<nlayers) for(int l=nlayersc;l<nlayers;l++) EplantCoh(c,l)=0.0;
    // Rcout<<"\n";
    EplantVec = EplantVec + EplantCoh(c,_);
    Eplant[c] = sum(EplantCoh(c,_));
    PlantPsi[c] = minPsi;
    DDS[c] = Phe[c]*(1.0-Psi2K(PlantPsi[c], VCstem_d[c], VCstem_c[c]));
    transpiration[c] = Eplant[c]; 
  }
  
  //Potential soil evaporation
  double PETsoil = meteoland::penmanmonteith(200.0, elevation, tmin, tmax, rhmin, rhmax, Rn*LgroundSWR, wind);
  

  
  //Evaporation from bare soil
  double Gsoil = soil["Gsoil"];
  double Ksoil = soil["Ksoil"];
  double Esoil = soilevaporation((Water_FC[0]*(1.0 - W[0])), PETsoil, Gsoil);
  NumericVector EsoilVec(nlayers,0.0);//Exponential decay to divide bare soil evaporation among layers
  for(int l=0;l<nlayers;l++) {
    double cumAnt = 0.0;
    double cumPost = 0.0;
    for(int l2=0;l2<l;l2++) cumAnt +=dVec[l2];
    cumPost = cumAnt+dVec[l];
    if(l<(nlayers-1)) EsoilVec[l] = Esoil*(exp(-Ksoil*cumAnt)-exp(-Ksoil*cumPost));
    else EsoilVec[l] = Esoil*exp(-Ksoil*cumAnt);
  }
  
  //Update layers
  for(int l=0;l<nlayers;l++) W[l] =std::max(W[l]-(EplantVec[l]+EsoilVec[l])/Water_FC[l],0.0);
  

  List l = List::create(_["PET"] = NA_REAL, _["NetPrec"] = NetPrec, _["Runon"] = runon, _["Infiltration"] = Infiltration, _["Runoff"] = Runoff, _["DeepDrainage"] = DeepDrainage,
                        _["LAIcell"] = LAIcell, _["Cm"] = Cm, _["Lground"] = LgroundPAR,
                        _["EsoilVec"] = EsoilVec, _["EplantVec"] = EplantVec, _["psiVec"] = psi,
                        _["EplantCoh"] = Eplant, _["psiCoh"] = PlantPsi, _["DDS"] = DDS);
  return(l);
}

NumericVector getTrackSpeciesTranspiration( NumericVector trackSpecies, NumericVector Eplant, DataFrame x) {
  int nTrackSpecies = trackSpecies.size();
  NumericVector Eplantsp(nTrackSpecies, 0.0);
  NumericVector SP = x["SP"];
  int nCoh = SP.size();
  int ts;
  for(int its =0;its<nTrackSpecies;its++) {
    ts = trackSpecies[its];
    for(int i=0;i<nCoh;i++) {
      if(SP[i]==ts) {
        Eplantsp[its] += Eplant[i];
      }
    }
  }
  return(Eplantsp);
}

NumericVector getTrackSpeciesDDS(NumericVector trackSpecies, NumericVector DDS, DataFrame x) {
  int nTrackSpecies = trackSpecies.size();
  NumericVector DDSsp(nTrackSpecies, 0.0);
  NumericVector LAI = x["LAI"];
  NumericVector SP = x["SP"];
  int nCoh = LAI.size();
  int ts;
  double laiSum;
  for(int its =0;its<nTrackSpecies;its++) {
    ts = trackSpecies[its];
    laiSum = 0.0;
    for(int i=0;i<nCoh;i++) {
      if(SP[i]==ts) {
        DDSsp[its] += DDS[i]*LAI[i];
        laiSum +=LAI[i];
      }
    }
    DDSsp = DDSsp/laiSum;
  }
  return(DDSsp);
}


// [[Rcpp::export(".swbgridDay")]]
List swbgridDay(CharacterVector lct, List xList, List soilList, 
                IntegerVector waterO, List queenNeigh, List waterQ,
                NumericVector gddVec, NumericVector petVec, NumericVector rainVec, 
                NumericVector erVec, NumericVector trackSpecies) {
  int nX = xList.size();
  int nTrackSpecies = trackSpecies.size();
  NumericVector NetPrec(nX,NA_REAL), Runon(nX,0.0), Infiltration(nX,NA_REAL);
  NumericVector Runoff(nX,NA_REAL), DeepDrainage(nX,NA_REAL), W1(nX,NA_REAL), W2(nX,NA_REAL);
  NumericVector W3(nX,NA_REAL), Esoil(nX,NA_REAL), Eplant(nX,NA_REAL);
  NumericMatrix Transpiration(nX, nTrackSpecies), DDS(nX, nTrackSpecies);
  double runoffExport = 0.0;
  for(int i=0;i<nX;i++) {
    //get next cell in order
    int iCell = waterO[i]-1; //Decrease index!!!!
    if((lct[iCell]=="Wildland") || (lct[iCell]=="Agriculture") ) {
      DataFrame x = Rcpp::as<Rcpp::DataFrame>(xList[iCell]);
      List soil = Rcpp::as<Rcpp::List>(soilList[iCell]);
      //Run daily soil water balance for the current cell
      List res = swbDay1(x, soil, gddVec[iCell], petVec[iCell], rainVec[iCell], erVec[iCell], Runon[iCell]);
      NetPrec[iCell] = res["NetPrec"];
      Runon[iCell] = res["Runon"];
      Infiltration[iCell] = res["Infiltration"];
      Runoff[iCell] = res["Runoff"];
      DeepDrainage[iCell] = res["DeepDrainage"];
      Esoil[iCell] = sum(Rcpp::as<Rcpp::NumericVector>(res["EsoilVec"]));
      NumericVector EplantCoh = res["EplantCoh"];
      NumericVector DDScell = res["DDS"];
      Eplant[iCell] = sum(EplantCoh);
      if(nTrackSpecies>0) {
        Transpiration(iCell,_) = getTrackSpeciesTranspiration(trackSpecies, EplantCoh, x);
        DDS(iCell,_) = getTrackSpeciesDDS(trackSpecies, DDScell, x);
      }
      NumericVector W = soil["W"];
      W1[iCell] = W[0];
      W2[iCell] = W[1];
      W3[iCell] = W[2];

      //Assign runoff to runon of neighbours
      double ri =  Runoff[iCell];
      if(ri>0.0) {
        IntegerVector ni = Rcpp::as<Rcpp::IntegerVector>(queenNeigh[iCell]);
        NumericVector qi = Rcpp::as<Rcpp::NumericVector>(waterQ[iCell]);
        if(ni.size()>0) {
          for(int j=0;j<ni.size();j++) Runon[ni[j]-1] += qi[j]*ri; //decrease index
        } else {
          runoffExport += ri; //If no suitable neighbours add ri to landscape export via runoff
        }
      }
    } else if(lct[iCell]=="Rock") {//all Precipitation becomes runoff if cell is rock outcrop
      Runoff[iCell] =  Runon[iCell]+rainVec[iCell];
      double ri = Runoff[iCell];
      if(ri>0.0) {
        IntegerVector ni = Rcpp::as<Rcpp::IntegerVector>(queenNeigh[iCell]);
        NumericVector qi = Rcpp::as<Rcpp::NumericVector>(waterQ[iCell]);
        if(ni.size()>0) {
          for(int j=0;j<ni.size();j++) Runon[ni[j]-1] += qi[j]*ri;//decrease index
        } else {
          runoffExport += ri; //If no suitable neighbours add ri to landscape export via runoff
        }
      }
    } else if(lct[iCell]=="Static") {
      // static cells receive water from other cells or Precipitation
      // but do not export to the atmosphere contribute nor to other cells.
      // Hence, water balance over the landscape is achieved by
      // adding this water to the landscape export via runoff.
      runoffExport += Runon[iCell] + rainVec[iCell];
    }
  }
  DataFrame waterBalance = DataFrame::create(_["NetPrec"] = NetPrec, _["Runon"] = Runon, _["Infiltration"] = Infiltration,
                                   _["Runoff"] = Runoff, _["DeepDrainage"] = DeepDrainage,
                                   _["W1"] = W1, _["W2"] = W2, _["W3"] = W3,
                                   _["Esoil"] = Esoil, _["Eplant"] = Eplant);
  return(List::create(_["WaterBalance"] = waterBalance,
                      _["RunoffExport"] = runoffExport,
                      _["Transpiration"] = Transpiration,
                      _["DDS"] = DDS));
}



void checkswbInput(List x, List soil, String transpirationMode) {
  if(!x.containsElementNamed("above")) stop("above missing in swbInput");
  if(!x.containsElementNamed("V")) stop("V missing in swbInput");
  if(!x.containsElementNamed("Sgdd")) stop("Sgdd missing in swbInput");
  if(!x.containsElementNamed("k")) stop("k missing in swbInput");
  if(!x.containsElementNamed("g")) stop("g missing in swbInput");
  if(transpirationMode=="Simple") {
    if(!x.containsElementNamed("Psi_Extract")) stop("Psi_Extract missing in swbInput");
    if(!x.containsElementNamed("WUE")) stop("WUE missing in swbInput");
  } else if(transpirationMode=="Sperry") {
    if(!x.containsElementNamed("VGrhizo_kmax")) stop("VCstem_kmax missing in swbInput");
    if(!soil.containsElementNamed("VG_n")) stop("VG_n missing in soil");
    if(!soil.containsElementNamed("VG_alpha")) stop("VG_alpha missing in soil");
    if(!x.containsElementNamed("VCstem_kmax")) stop("VCstem_kmax missing in swbInput");
    if(!x.containsElementNamed("VCstem_c")) stop("VCstem_c missing in swbInput");
    if(!x.containsElementNamed("VCstem_d")) stop("VCstem_d missing in swbInput");
    if(!x.containsElementNamed("VCroot_kmax")) stop("VCroot_kmax missing in swbInput");
    if(!x.containsElementNamed("VCroot_c")) stop("VCroot_c missing in swbInput");
    if(!x.containsElementNamed("VCroot_d")) stop("VCroot_d missing in swbInput");
    if(!x.containsElementNamed("Gwmax")) stop("Gwmax missing in swbInput");
    if(!x.containsElementNamed("Vmax298")) stop("Vmax298 missing in swbInput");
    if(!x.containsElementNamed("Jmax298")) stop("Jmax298 missing in swbInput");
  }
  if(!soil.containsElementNamed("W")) stop("W missing in soil");
  if(!soil.containsElementNamed("psi")) stop("psi missing in soil");
  if(!soil.containsElementNamed("dVec")) stop("dVec missing in soil");
  if(!soil.containsElementNamed("Theta_FC")) stop("Theta_FC missing in soil");
  if(!soil.containsElementNamed("Water_FC")) stop("Water_FC missing in soil");
  if(!soil.containsElementNamed("macro")) stop("macro missing in soil");
  if(!soil.containsElementNamed("clay")) stop("clay missing in soil");
  if(!soil.containsElementNamed("sand")) stop("sand missing in soil");
}

// [[Rcpp::export("swb")]]
List swb(List x, List soil, DataFrame meteo, double latitude = NA_REAL, double elevation = NA_REAL, double slope = NA_REAL, double aspect = NA_REAL) {
  String transpirationMode = x["TranspirationMode"];
  bool verbose = x["verbose"];
  
  checkswbInput(x, soil, transpirationMode);
    
  NumericVector Precipitation = meteo["Precipitation"];
  NumericVector Temperature = meteo["MeanTemperature"];
  IntegerVector DOY = meteo["DOY"];
  CharacterVector dateStrings = meteo.attr("row.names");
  
  NumericVector GDD = gdd(DOY, Temperature, 5.0);
  NumericVector ER = er(DOY);
  
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  NumericVector SP = above["SP"];
  int numCohorts = SP.size();
  int numDays = Precipitation.size();
  
  NumericVector Water_FC = soil["Water_FC"];
  NumericVector W = soil["W"];
  
  int nlayers = W.size();
  NumericVector PET(numDays);

  //Water balance variables
  NumericVector Esoil(numDays);
  NumericVector LAIcell(numDays);
  NumericVector Cm(numDays);
  NumericVector Lground(numDays);
  NumericVector Runoff(numDays);
  NumericVector NetPrec(numDays);
  NumericVector Interception(numDays);
  NumericVector Infiltration(numDays);
  NumericVector DeepDrainage(numDays);
  NumericMatrix Eplantdays(numDays, nlayers);
  NumericMatrix Wdays(numDays, nlayers); //Soil moisture content in relation to field capacity
  NumericMatrix psidays(numDays, nlayers);
  NumericMatrix MLdays(numDays, nlayers);
  NumericVector MLTot(numDays, 0.0);
  
  NumericMatrix PlantPsi(numDays, numCohorts);
  NumericMatrix PlantStress(numDays, numCohorts);
  NumericMatrix PlantTranspiration(numDays, numCohorts);
  NumericMatrix PlantPhotosynthesis(numDays, numCohorts);
  NumericVector EplantCohTot(numCohorts, 0.0);
  
  
  Wdays(0,_) = W;

  if(verbose) {
    for(int l=0;l<nlayers;l++) Rcout << "W"<<(l+1)<<"i:"<< round(100*W[l])/100<<" ";
    Rcout<<"\n";
  }

  if(verbose) Rcout << "Daily balance:";
  if(transpirationMode=="Simple") {
    PET = meteo["PET"];
    NumericVector Tday = meteo["MeanTemperature"];
    for(int i=0;i<numDays;i++) {
      if(verbose) Rcout<<".";
      List s = swbDay1(x, soil, GDD[i], Tday[i], PET[i], Precipitation[i], ER[i], 0.0, false); //No Runon in simulations for a single cell

      Lground[i] = s["Lground"];
      Esoil[i] = sum(Rcpp::as<Rcpp::NumericVector>(s["EsoilVec"]));
      LAIcell[i] = s["LAIcell"];
      Cm[i] = s["Cm"];
      DeepDrainage[i] = s["DeepDrainage"];
      Infiltration[i] = s["Infiltration"];
      Runoff[i] = s["Runoff"];
      NetPrec[i] = s["NetPrec"];
      Interception[i] = Precipitation[i]-NetPrec[i];
      
      NumericVector psi = s["psiVec"];
      psidays(i,_) = psi;

      NumericVector EplantCoh = s["EplantCoh"];
      NumericVector EplantVec = s["EplantVec"];
      Eplantdays(i,_) = EplantVec;
      PlantPhotosynthesis(i,_) = Rcpp::as<Rcpp::NumericVector>(x["Photosynthesis"]);
      PlantTranspiration(i,_) = EplantCoh;
      PlantStress(i,_) = Rcpp::as<Rcpp::NumericVector>(s["DDS"]);
      PlantPsi(i,_) = Rcpp::as<Rcpp::NumericVector>(s["psiCoh"]);
      EplantCohTot = EplantCohTot + EplantCoh;
      
      if(i<(numDays-1)) Wdays(i+1,_) = W;
    }
  } else if(transpirationMode=="Sperry") {
    if(NumericVector::is_na(latitude)) stop("Value for 'latitude' should not be missing.");
    if(NumericVector::is_na(elevation)) stop("Value for 'elevation' should not be missing.");
    NumericVector MinTemperature = meteo["MinTemperature"];
    NumericVector MaxTemperature = meteo["MaxTemperature"];
    NumericVector MinRelativeHumidity = meteo["MinRelativeHumidity"];
    NumericVector MaxRelativeHumidity = meteo["MaxRelativeHumidity"];
    NumericVector Radiation = meteo["Radiation"];
    NumericVector WindSpeed = meteo["WindSpeed"];
    for(int i=0;i<numDays;i++) {
      std::string c = as<std::string>(dateStrings[i]);
      if(verbose) Rcout<<".";
      //Julian day
      int J = meteoland::radiation_julianDay(std::atoi(c.substr(0, 4).c_str()),std::atoi(c.substr(5,2).c_str()),std::atoi(c.substr(8,2).c_str()));
      double delta = meteoland::radiation_solarDeclination(J);
      List s = swbDay2(x, soil, GDD[i], MinTemperature[i], MaxTemperature[i], 
                       MinRelativeHumidity[i], MaxRelativeHumidity[i], Radiation[i], WindSpeed[i], 
                       latitude, elevation, slope, aspect, delta, Precipitation[i], ER[i], 0.0, false);
      
      PET[i] = s["PET"];
      Lground[i] = s["Lground"];
      Esoil[i] = sum(Rcpp::as<Rcpp::NumericVector>(s["EsoilVec"]));
      LAIcell[i] = s["LAIcell"];
      Cm[i] = s["Cm"];
      DeepDrainage[i] = s["DeepDrainage"];
      Infiltration[i] = s["Infiltration"];
      Runoff[i] = s["Runoff"];
      NetPrec[i] = s["NetPrec"];
      Interception[i] = Precipitation[i]-NetPrec[i];
      
      NumericVector psi = s["psiVec"];
      psidays(i,_) = psi;
      
      NumericVector EplantCoh = s["EplantCoh"];
      NumericVector EplantVec = s["EplantVec"];
      Eplantdays(i,_) = EplantVec;
      PlantTranspiration(i,_) = EplantCoh;
      PlantPhotosynthesis(i,_) = Rcpp::as<Rcpp::NumericVector>(x["Photosynthesis"]);
      PlantPsi(i,_) = Rcpp::as<Rcpp::NumericVector>(s["psiCoh"]);
      PlantStress(i,_) = Rcpp::as<Rcpp::NumericVector>(s["DDS"]);
      EplantCohTot = EplantCohTot + EplantCoh;
      
      if(i<(numDays-1)) Wdays(i+1,_) = W;
    }
  }
  if(verbose) Rcout << "done\n";
  
  NumericVector Eplanttot(numDays,0.0);
  for(int l=0;l<nlayers;l++) {
    MLdays(_,l) = Wdays(_,l)*Water_FC[l]; 
    MLTot = MLTot + MLdays(_,l);
    Eplanttot = Eplanttot + Eplantdays(_,l);
  }

  NumericVector Etot = Eplanttot+Esoil;
  
  if(verbose) {
    double Precipitationsum = sum(Precipitation);
    double NetPrecsum = sum(NetPrec);
    double Interceptionsum = sum(Interception);
    double Esoilsum = sum(Esoil);
    double Runoffsum  = sum(Runoff);
    double Infiltrationsum  = sum(Infiltration);
    double DeepDrainagesum = sum(DeepDrainage);

    double Eplantsum = sum(Eplanttot);
    
    Rcout<<"Total Precipitation (mm) "  <<round(Precipitationsum) <<"\n";
    Rcout<<"Interception (mm) " << round(Interceptionsum)  <<" Net Prec (mm) " << round(NetPrecsum) <<"\n";
    Rcout<<"Infiltration (mm) " << round(Infiltrationsum)  <<
      " Runoff (mm) " << round(Runoffsum) <<
        " Deep drainage (mm) "  << round(DeepDrainagesum)  <<"\n";
    Rcout<<"Esoil (mm) " << round(Esoilsum) <<
      " Eplant (mm) "  <<round(Eplantsum) <<"\n";
    for(int l=0;l<nlayers;l++) Rcout << "W"<<(l+1)<<"f:"<< round(100*W[l])/100<<" ";
    Rcout<<"\n";
    Rcout<<"Final volume: "<< round(MLTot[numDays-1])<<"\n\n";
    
  }
  if(verbose) Rcout<<"Building SWB and DWB output ...";
  
   Rcpp::DataFrame SWB = DataFrame::create(_["W"]=Wdays, _["ML"]=MLdays,_["MLTot"]=MLTot,_["psi"]=psidays);
   Rcpp::DataFrame DWB = DataFrame::create(_["LAIcell"]=LAIcell, _["Cm"]=Cm, _["Lground"] = Lground, _["PET"]=PET, 
                                           _["Precipitation"] = Precipitation, _["NetPrec"]=NetPrec,_["Infiltration"]=Infiltration, _["Runoff"]=Runoff, _["DeepDrainage"]=DeepDrainage, 
                                           _["Etot"]=Etot,_["Esoil"]=Esoil,
                                           _["Eplanttot"]=Eplanttot,
                                           _["Eplant"]=Eplantdays);
  
   SWB.attr("row.names") = meteo.attr("row.names") ;
   DWB.attr("row.names") = meteo.attr("row.names") ;
  if(verbose) Rcout<<"plant output ...";
  
  PlantTranspiration.attr("dimnames") = List::create(meteo.attr("row.names"), SP.attr("names"));
  PlantStress.attr("dimnames") = List::create(meteo.attr("row.names"), SP.attr("names")) ;
  
  if(verbose) Rcout<<"list ...";
  List l = List::create(Named("TranspirationMode") = transpirationMode,
                        Named("NumSoilLayers") = nlayers,
                        Named("DailyBalance")=DWB, 
                        Named("SoilWaterBalance")=SWB,
                        Named("PlantTranspiration") = PlantTranspiration,
                        Named("PlantPhotosynthesis") = PlantPhotosynthesis,
                        Named("PlantPsi") = PlantPsi, 
                        Named("PlantStress") = PlantStress);
  l.attr("class") = CharacterVector::create("swb","list");
  if(verbose) Rcout<<"done.\n";
  return(l);
}
