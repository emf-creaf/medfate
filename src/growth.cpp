#include <numeric>
#include "lightextinction.h"
#include "biophysicsutils.h"
#include "forestutils.h"
#include "hydraulics.h"
#include "soil.h"
#include "swb.h"
#include <Rcpp.h>
#include <meteoland.h>
using namespace Rcpp;

const double leafCperDry = 0.3; //g C · g dry-1
const double rootCperDry = 0.4959; //g C · g dry-1

const double leaf_RR = 0.005;
const double stem_RR = 0.0001;
const double root_RR = 0.003;
const double Q10_resp = 2.0;

const double dailySAturnoverProportion = 0.0001261398; //Equivalent to annual 4.5% 1-(1-0.045)^(1.0/365)

NumericVector carbonCompartments(double SA, double LAI, double H, double Z, double N, double SLA, double WoodDens, double WoodC) {
  double B_leaf = leafCperDry*1000.0*(LAI/(N/10000.0))/SLA; //Biomass in g C · ind-1
  double B_stem = WoodC*SA*(H+Z)*WoodDens;
  double B_fineroot = B_leaf/2.5;
  return(NumericVector::create(B_leaf, B_stem, B_fineroot)); 
}
double dailyRespiration(double B_leaf,double B_stem,double B_root, double Tmean) {
  double q = pow(Q10_resp,(Tmean-20.0)/10.0);
  return((B_leaf*leaf_RR + B_stem*stem_RR + B_root*root_RR)*q);
}

double storageTransferMaximumRate(double fastCstorage, double fastCstoragemax) {
  double f = ((2.0/(1.0+exp(-5.0*((fastCstorage/fastCstoragemax)-0.5)/0.5)))-1.0);
  return((fastCstoragemax/30.0)*f);
}
double temperatureGrowthFactor(double Tmean) {
  double Tlow = 5.0;
  double Thigh = 40.0;
  double Topt = 25.0;
  double f = ((Tmean-Tlow)*(Thigh-Tmean))/((Topt-Tlow)*(Thigh-Topt));
  return(std::min(std::max(f,0.0),1.0));
}
double turgorGrowthFactor(double psi, double psi_tlp) {
  return(std::max(0.0,1.0 - pow(exp((psi/psi_tlp)-1),5.0)));
}

void checkgrowthInput(List x, List soil, String transpirationMode) {
  if(!x.containsElementNamed("above")) stop("above missing in growthInput");
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  if(!above.containsElementNamed("LAI_live")) stop("LAI_live missing in growthInput$above");
  if(!above.containsElementNamed("LAI_expanded")) stop("LAI_expanded missing in growthInput$above");
  if(!above.containsElementNamed("LAI_dead")) stop("LAI_dead missing in growthInput$above");
  if(!above.containsElementNamed("LAI_predrought")) stop("LAI_predrought missing in growthInput$above");
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
  if(!paramsBase.containsElementNamed("k")) stop("k missing in growthInput$paramsBase");
  if(!paramsBase.containsElementNamed("g")) stop("g missing in growthInput$paramsBase");
  
  if(!x.containsElementNamed("paramsGrowth")) stop("paramsGrowth missing in growthInput");
  DataFrame paramsGrowth = Rcpp::as<Rcpp::DataFrame>(x["paramsGrowth"]);
  if(!paramsGrowth.containsElementNamed("SLA")) stop("SLA missing in growthInput$paramsGrowth");
  if(!paramsGrowth.containsElementNamed("Al2As")) stop("Al2As missing in growthInput$paramsGrowth");
  if(!paramsGrowth.containsElementNamed("WoodC")) stop("WoodC missing in growthInput$paramsGrowth");
  if(!paramsGrowth.containsElementNamed("WoodDens")) stop("WoodDens missing in growthInput$paramsGrowth");
  // if(!paramsGrowth.containsElementNamed("Cstoragepmax")) stop("Cstoragepmax missing in growthInput$paramsGrowth");
  if(!paramsGrowth.containsElementNamed("RGRmax")) stop("RGRmax missing in growthInput$paramsGrowth");
  
  if(!x.containsElementNamed("paramsTransp")) stop("paramsTransp missing in growthInput");
  DataFrame paramsTransp = Rcpp::as<Rcpp::DataFrame>(x["paramsTransp"]);
  if(transpirationMode=="Simple") {
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
  if(!soil.containsElementNamed("psi")) stop("psi missing in soil");
  if(!soil.containsElementNamed("dVec")) stop("dVec missing in soil");
  if(!soil.containsElementNamed("Theta_FC")) stop("Theta_FC missing in soil");
  if(!soil.containsElementNamed("Water_FC")) stop("Water_FC missing in soil");
  if(!soil.containsElementNamed("macro")) stop("macro missing in soil");
  if(!soil.containsElementNamed("clay")) stop("clay missing in soil");
  if(!soil.containsElementNamed("sand")) stop("sand missing in soil");
}

// [[Rcpp::export("growth")]]
List growth(List x, List soil, DataFrame meteo, double latitude = NA_REAL, double elevation = NA_REAL, double slope = NA_REAL, double aspect = NA_REAL) {
  String transpirationMode = x["TranspirationMode"];
  String storagePool = x["storagePool"];
  bool verbose = x["verbose"];
  bool allowEmbolism = x["allowEmbolism"];
  checkgrowthInput(x, soil, transpirationMode);
  
  NumericVector Precipitation = meteo["Precipitation"];
  NumericVector MeanTemperature = meteo["MeanTemperature"];
  IntegerVector DOY = meteo["DOY"];
  NumericVector MinTemperature, MaxTemperature, MinRelativeHumidity, MaxRelativeHumidity, Radiation, WindSpeed, PET;
  int numDays = Precipitation.size();
  if(transpirationMode=="Sperry") {
    if(NumericVector::is_na(latitude)) stop("Value for 'latitude' should not be missing.");
    if(NumericVector::is_na(elevation)) stop("Value for 'elevation' should not be missing.");
    MinTemperature = meteo["MinTemperature"];
    MaxTemperature = meteo["MaxTemperature"];
    MinRelativeHumidity = meteo["MinRelativeHumidity"];
    MaxRelativeHumidity = meteo["MaxRelativeHumidity"];
    Radiation = meteo["Radiation"];
    WindSpeed = meteo["WindSpeed"];
    PET = NumericVector(numDays);
  } else {
    PET = meteo["PET"];
    WindSpeed = meteo["WindSpeed"];
  }
  
  CharacterVector dateStrings = meteo.attr("row.names");
  
  NumericVector GDD = gdd(DOY, MeanTemperature, 5.0);
  NumericVector ER = er(DOY);

  //Aboveground parameters  
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  NumericVector SP = above["SP"];
  NumericVector DBH = above["DBH"];
  NumericVector Cover = above["Cover"];
  NumericVector H = above["H"];
  NumericVector N = above["N"];
  NumericVector CR = above["CR"];
  NumericVector LAI_live = above["LAI_live"];
  NumericVector LAI_expanded = above["LAI_expanded"];
  NumericVector LAI_dead = above["LAI_dead"];
  NumericVector SA = above["SA"];
  NumericVector Psi_leafmin = above["Psi_leafmin"];
  NumericVector LAI_predrought = above["LAI_predrought"];
  NumericVector fastCstorage, slowCstorage;
  int numCohorts = SP.size();

  //Belowground parameters  
  List below = Rcpp::as<Rcpp::List>(x["below"]);
  NumericVector Z = Rcpp::as<Rcpp::NumericVector>(below["Z"]);

  //Base parameters
  DataFrame paramsBase = Rcpp::as<Rcpp::DataFrame>(x["paramsBase"]);
  NumericVector Sgdd = paramsBase["Sgdd"];
  NumericVector k = paramsBase["k"];
  

  //Transpiration parameters
  DataFrame paramsTransp = Rcpp::as<Rcpp::DataFrame>(x["paramsTransp"]);
  NumericVector xylem_kmax, VCstem_kmax, Psi_Extract, VCstem_c, VCstem_d;
  if(transpirationMode=="Sperry") {
    xylem_kmax = paramsTransp["xylem_kmax"];
    VCstem_kmax = paramsTransp["VCstem_kmax"];
    VCstem_c = paramsTransp["VCstem_c"];
    VCstem_d = paramsTransp["VCstem_d"];
  } else if(transpirationMode == "Simple"){
    Psi_Extract = paramsTransp["Psi_Extract"];
  }

  NumericVector Water_FC = soil["Water_FC"];
  NumericVector W = soil["W"];
  
  int nlayers = W.size();
  
  //Growth parameters
  DataFrame paramsGrowth = Rcpp::as<Rcpp::DataFrame>(x["paramsGrowth"]);
  NumericVector SLA = Rcpp::as<Rcpp::NumericVector>(paramsGrowth["SLA"]);
  NumericVector Al2As = Rcpp::as<Rcpp::NumericVector>(paramsGrowth["Al2As"]);
  NumericVector WoodC = Rcpp::as<Rcpp::NumericVector>(paramsGrowth["WoodC"]);
  NumericVector WoodDens = Rcpp::as<Rcpp::NumericVector>(paramsGrowth["WoodDens"]);
  NumericVector RGRmax = Rcpp::as<Rcpp::NumericVector>(paramsGrowth["RGRmax"]);
  NumericVector Cstoragepmax, slowCstorage_max(numCohorts), fastCstorage_max(numCohorts);
  if(storagePool=="one") {
    fastCstorage = above["fastCstorage"];
    Cstoragepmax = Rcpp::as<Rcpp::NumericVector>(paramsGrowth["Cstoragepmax"]);
  } else if(storagePool=="two") {
    fastCstorage = above["fastCstorage"];
    slowCstorage = above["slowCstorage"];
    Cstoragepmax = Rcpp::as<Rcpp::NumericVector>(paramsGrowth["Cstoragepmax"]);
  }
  
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
  
  //Water balance output variables
  NumericMatrix PlantPsi(numDays, numCohorts);
  NumericMatrix PlantStress(numDays, numCohorts);
  NumericMatrix PlantTranspiration(numDays, numCohorts);
  NumericMatrix PlantRespiration(numDays, numCohorts);
  NumericMatrix PlantCstorageFast(numDays, numCohorts);
  NumericMatrix PlantCstorageSlow(numDays, numCohorts);
  NumericMatrix PlantSA(numDays, numCohorts);
  NumericMatrix PlantSAgrowth(numDays, numCohorts);
  NumericMatrix PlantPhotosynthesis(numDays, numCohorts);
  NumericMatrix PlantLAIdead(numDays, numCohorts);
  NumericMatrix PlantLAIlive(numDays, numCohorts);
  NumericVector SAgrowthcum(numCohorts, 0.0);
  
  //Water balance output variables
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
  
  Wdays(0,_) = W;
  
  if(verbose) {
    for(int l=0;l<nlayers;l++) Rcout << "W"<<(l+1)<<"i:"<< round(100*W[l])/100<<" ";
    Rcout<<"\n";
  }
  
  if(verbose) Rcout << "Daily growth:";
  List s;
  for(int i=0;i<numDays;i++) {
    if(verbose) Rcout<<".";
    double ws = WindSpeed[i];
    if(NumericVector::is_na(ws)) ws = 1.0; //Default 1 m/s -> 10% of fall every day
    
    //1. Phenology and leaf fall
    NumericVector phe = leafDevelopmentStatus(Sgdd, GDD[i]);
    for(int j=0;j<numCohorts;j++) {
      LAI_dead[j] *= exp(-1.0*(ws/10.0)); //Decrease dead leaf area according to wind speed
      double LAI_exp_prev=0.0;
      if(i>0) LAI_exp_prev= LAI_expanded[j]; //Store previous value
      LAI_expanded[j] = LAI_live[j]*phe[j]; //Update expanded leaf area (will decrease if LAI_live decreases)
      LAI_dead[j] += std::max(0.0, LAI_exp_prev-LAI_expanded[j]);//Check increase dead leaf area if expanded leaf area has decreased
    }
    
    //2. Water balance and photosynthesis
    if(transpirationMode=="Simple") {
      s = swbDay1(x, soil, MeanTemperature[i], PET[i], Precipitation[i], ER[i], 0.0, false); //No Runon in simulations for a single cell
    } else if(transpirationMode=="Sperry") {
      std::string c = as<std::string>(dateStrings[i]);
      int J = meteoland::radiation_julianDay(std::atoi(c.substr(0, 4).c_str()),std::atoi(c.substr(5,2).c_str()),std::atoi(c.substr(8,2).c_str()));
      double delta = meteoland::radiation_solarDeclination(J);
      s = swbDay2(x, soil, MinTemperature[i], MaxTemperature[i], 
                       MinRelativeHumidity[i], MaxRelativeHumidity[i], Radiation[i], WindSpeed[i], 
                       latitude, elevation, slope, aspect, delta, Precipitation[i], ER[i], 0.0, false);
      
      PET[i] = s["PET"];
    }    
    
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
    NumericVector An =  Rcpp::as<Rcpp::NumericVector>(x["Photosynthesis"]);
    PlantPhotosynthesis(i,_) = An;
    PlantTranspiration(i,_) = EplantCoh;
    PlantStress(i,_) = Rcpp::as<Rcpp::NumericVector>(s["DDS"]);
    NumericVector psiCoh =  Rcpp::as<Rcpp::NumericVector>(s["psiCoh"]);;
    PlantPsi(i,_) = psiCoh;
      
    //3. Carbon balance and growth
    double B_leaf_expanded, B_stem, B_fineroot;
    for(int j=0;j<numCohorts;j++){
      //3.1 Live biomass and maximum C pool capacity
      NumericVector compartments = carbonCompartments(SA[j], LAI_expanded[j], H[j], Z[j], N[j], SLA[j], WoodDens[j], WoodC[j]);
      B_leaf_expanded = compartments[0];
      B_stem = compartments[1];
      B_fineroot = compartments[2];
      if(storagePool == "one") {
        fastCstorage_max[j] = Cstoragepmax[j]*(B_leaf_expanded+B_stem+B_fineroot);
      } else if(storagePool == "two") {
        fastCstorage_max[j] = 0.05*(B_leaf_expanded+B_stem+B_fineroot);
        slowCstorage_max[j] = std::max(slowCstorage_max[j],(Cstoragepmax[j]-0.05)*(B_leaf_expanded+B_stem+B_fineroot)); //Slow pool capacity cannot decrease
      }
      
      //3.2 Respiration and photosynthesis 
      double Anj = An[j]/(N[j]/10000.0); //Translate g C · m-2 to g C · ind-1
      double Rj = dailyRespiration(B_leaf_expanded, B_stem, B_fineroot, MeanTemperature[i]);
      //3.3. Carbon balance, update of fast C pool and C available for growth
      double growthAvailableC = 0.0;
      if(storagePool=="none") {
        growthAvailableC = std::max(0.0,(Anj-Rj));
      } else  {
        fastCstorage[j] = std::max(0.0,fastCstorage[j]+(Anj-Rj));
        growthAvailableC = fastCstorage[j];
      }
      //3.4 Growth in LAI_live and SA
      double deltaSAturnover = (dailySAturnoverProportion/(1.0+15*exp(-0.01*H[j])))*SA[j];
      double f_turgor = turgorGrowthFactor(psiCoh[j],-1.1);
      double deltaSAgrowth = 0.0;
      if(f_turgor>0.0) { //Growth is possible
        double costLA = 0.1*leafCperDry*(Al2As[j]/SLA[j]); //Construction cost in g C·cm-2 of sapwood
        double costSA = WoodC[j]*(H[j]+Z[j])*WoodDens[j];  //Construction cost in g C·cm-2 of sapwood
        double costFR = costLA/2.5;
        double cost = 1.3*(costLA+costSA+costFR);  //Construction cost in g C·cm-2 of sapwood (including 30% growth respiration)
        double deltaSAsource = growthAvailableC/cost;
        double f_temp = temperatureGrowthFactor(MeanTemperature[i]);
        double deltaSAsink = RGRmax[j]*SA[j]*f_temp*f_turgor;
        if(storagePool!="none") {
          double f_sugars = (1.0/(1.0+exp(-5.0*((fastCstorage[j]/fastCstorage_max[j])-0.5)/0.5)));
          deltaSAsink *=f_sugars;
        } 
        deltaSAgrowth = std::min(deltaSAsink, deltaSAsource);
        
        //update pools
        if(storagePool != "none") {
          fastCstorage[j] = fastCstorage[j]-deltaSAgrowth*cost; //Remove construction costs from (fast) C pool
        }
        SA[j] = SA[j] + deltaSAgrowth - deltaSAturnover; //Update sapwood area
        LAI_dead[j] += (N[j]/10000.0)*(deltaSAturnover/10000.0)*Al2As[j]; //Update dead LAI
        LAI_live[j] += (N[j]/10000.0)*((deltaSAgrowth-deltaSAturnover)/10000.0)*Al2As[j]; //Update live LAI
        Psi_leafmin[j] = 0.0;
        LAI_predrought[j] = LAI_live[j];
        SAgrowthcum[j] += deltaSAgrowth; //Store cumulative SA growth (for structural variable update)
      } else if(allowEmbolism) { //Growth is not possible, evaluate embolism
        Psi_leafmin[j] = std::min(Psi_leafmin[j], psiCoh[j]);
        double propEmb = 0.0;
        if(transpirationMode == "Simple") {
          propEmb = 1.0-Psi2K(Psi_leafmin[j], Psi_Extract[j]);
        } else if(transpirationMode == "Sperry") {
          propEmb = 1.0-exp(-pow(Psi_leafmin[j]/VCstem_d[j],VCstem_c[j]));
        }
        SA[j] = SA[j] - deltaSAturnover; //Update sapwood area (only turnover)
        double LAIturnover = LAI_live[j] - (N[j]/10000.0)*(deltaSAturnover/10000.0)*Al2As[j];
        double LAIembolism = LAI_predrought[j]*(1.0-propEmb);
        double prevLive = LAI_live[j];
        LAI_live[j] = std::min(LAIturnover, LAIembolism);
        LAI_dead[j] += prevLive - LAI_live[j];
      }
      LAI_expanded[j] = LAI_live[j]*phe[j]; //Update expanded leaf area
      
      //update C pools
      if(storagePool == "one") {
        fastCstorage[j] = std::max(0.0,std::min(fastCstorage[j], fastCstorage_max[j]));
      } else if(storagePool == "two") {
        double maxtransferRate = storageTransferMaximumRate(fastCstorage[j], fastCstorage_max[j]);//Maximum transfer rate (maximum daily 5% of fast storage capacity)
        if(maxtransferRate>0.0) { //Transfer from fast to slow
          double transfer = std::min(maxtransferRate, fastCstorage[j]);
          fastCstorage[j] -= transfer;
          slowCstorage[j] += transfer*0.9; //10% cost in respiration (not added to slow pool)
        } else { //Transfer from slow to fast
          double transfer = std::min(-maxtransferRate, slowCstorage[j]);
          fastCstorage[j] += transfer*0.9; //10% cost in respiration (removed from slow pool)
          slowCstorage[j] -= transfer;
        }
        //Trim pools to maximum capacity
        fastCstorage[j] = std::max(0.0,std::min(fastCstorage[j], fastCstorage_max[j]));
        slowCstorage[j] = std::max(0.0,std::min(slowCstorage[j], slowCstorage_max[j]));
      }
      
      //3.5 Update stem conductance (Sperry mode)
      if(transpirationMode=="Sperry") {
        double hubberValue = (LAI_expanded[j]/(N[j]/10000.0))/(SA[j]/10000.0);
        VCstem_kmax[j]=maximumStemHydraulicConductance(xylem_kmax[j], hubberValue,H[j]);
        // Rcout<<Al2As[j]<<" "<< hubberValue<<" "<<VCstem_kmax[j]<<"\n";
      }
      
      //Output variables
      PlantRespiration(i,j) = Rj*(N[j]/10000.0); //Scaled to cohort level: Translate g C · ind-1 to g C · m-2
      if(storagePool=="one") {
        PlantCstorageFast(i,j) = fastCstorage[j];
      } else if(storagePool == "two") {
        PlantCstorageFast(i,j) = fastCstorage[j];
        PlantCstorageSlow(i,j) = slowCstorage[j];
      }
      PlantSA(i,j) = SA[j];
      PlantLAIlive(i,j) = LAI_live[j];
      PlantLAIdead(i,j) = LAI_dead[j];
      PlantSAgrowth(i,j) = deltaSAgrowth;
    }
    //4 Update structural variables
    if(((DOY[i]==1) & (i>0)) | (i==(numDays-1))) { 
      if(verbose) Rcout<<" [update structural variables] ";
      NumericVector deltaDBH(numCohorts, 0.0);
      for(int j=0;j<numCohorts; j++) {
        if(!NumericVector::is_na(DBH[j])) {
          deltaDBH[j] = 2.0*sqrt(pow(DBH[j]/2.0,2.0)+(SAgrowthcum[j]/PI)) - DBH[j];
          DBH[j] = DBH[j] + deltaDBH[j];
        } 
        SAgrowthcum[j] = 0.0; //Reset cumulative growth
      }

      NumericVector L = parcohortC(H, LAI_live, LAI_dead, k, CR);
      for(int j=0;j<numCohorts; j++) {
        if(!NumericVector::is_na(DBH[j])) {
          double fHmod = std::max(0.0,std::min(1.0,(1.0-((H[j]-137.0)/(Hmax[j]-137.0)))));
          double fHD = (fHDmin[j]*(L[j]/100.0) + fHDmax[j]*(1.0-(L[j]/100.0)))*fHmod;
          // Rcout << fHmod<<" "<< fHD<<" "<< L[j]<<"\n";
          H[j] = H[j] + fHD*deltaDBH[j];
        }
      }
      NumericVector crNew = treeCrownRatio(N, DBH, H, Acw, Bcw, Acr, B1cr, B2cr, B3cr, C1cr, C2cr);
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
          H[j] = pow(pow(10,6.0)*PV/(Aash[j]*CR[j]), 1.0/3.0); //Updates shrub height
          if(H[j]> Hmax[j]) { //Limit height (and update the former variables)
            H[j] = Hmax[j];
            PV = (Aash[j]*pow(H[j],2.0)/10000.0)*(H[j]/100.0)*CR[j];
            Wleaves = Absh[j]*pow(PV, Bbsh[j])/r635[j];
            double prevLive = LAI_live[j];
            LAI_live[j] = Wleaves*((N[j]/10000)*SLA[j]); //Update LAI_live to the maximum
            LAI_dead[j] += prevLive - LAI_live[j]; //Increment dead LAI with the difference
            LAI_expanded[j] = LAI_live[j]*phe[j]; //Update expanded leaf area
          }
          Cover[j] = (N[j]*Aash[j]*pow(H[j],2.0)/pow(10, 6.0)); //Updates shrub cover
        }
      }
    }
    if(i<(numDays-1)) Wdays(i+1,_) = W;
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
    double Eplantsum = sum(Eplanttot);
    
    Rcout<<"Total Precipitation (mm) "  <<round(Precipitationsum) <<"\n";
    for(int l=0;l<nlayers;l++) Rcout << "W"<<(l+1)<<"f:"<< round(100*W[l])/100<<" ";
    Rcout<<"\n";

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
                        Named("PlantRespiration") = PlantRespiration,
                        Named("PlantCstorageFast") = PlantCstorageFast,
                        Named("PlantCstorageSlow") = PlantCstorageSlow,
                        Named("PlantSAgrowth") = PlantSAgrowth,
                        Named("PlantSA")=PlantSA,
                        Named("PlantPsi") = PlantPsi, 
                        Named("PlantStress") = PlantStress,
                        Named("PlantLAIdead") = PlantLAIdead,
                        Named("PlantLAIlive") = PlantLAIlive);
  l.attr("class") = CharacterVector::create("growth","list");
  if(verbose) Rcout<<"done.\n";
  return(l);
}