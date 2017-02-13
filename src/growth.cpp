#include <numeric>
#include "lightextinction.h"
#include "biophysicsutils.h"
#include "forestutils.h"
#include "soil.h"
#include "swb.h"
#include <Rcpp.h>
#include <meteoland.h>
using namespace Rcpp;

const double leafCperDry = 0.3; //g C · g dry-1
const double rootCperDry = 0.4959; //g C · g dry-1

const double leaf_RR = 0.001;
const double stem_RR = 0.00001;
const double root_RR = 0.001;
const double Q10_resp = 2.0;

const double dailySAturnoverProportion = 0.0001261398; //Equals to annual 4.5% 1-(1-0.045)^(1.0/365)

double dailyRespiration(double B_leaf,double B_stem,double B_root, double Tmean) {
  double q = pow(Q10_resp,(Tmean-20.0)/10.0);
  return((B_leaf*leaf_RR + B_stem*stem_RR + B_root*root_RR)*q);
}

void checkgrowthInput(List x, List soil, String transpirationMode) {
  if(!x.containsElementNamed("above")) stop("above missing in growthInput");
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  if(!above.containsElementNamed("LAI_live")) stop("LAI_live missing in growthInput$above");
  if(!above.containsElementNamed("LAI_expanded")) stop("LAI_expanded missing in growthInput$above");
  if(!above.containsElementNamed("LAI_dead")) stop("LAI_dead missing in growthInput$above");
  if(!above.containsElementNamed("LAI_predrought")) stop("LAI_predrought missing in growthInput$above");
  if(!above.containsElementNamed("SA")) stop("SA missing in growthInput$above");
  if(!above.containsElementNamed("Cstorage")) stop("Cstorage missing in growthInput$above");
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
  if(!paramsGrowth.containsElementNamed("Cstoragepmax")) stop("Cstoragepmax missing in growthInput$paramsGrowth");
  
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
  bool verbose = x["verbose"];
  
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
  NumericVector H = above["H"];
  NumericVector N = above["N"];
  NumericVector LAI_live = above["LAI_live"];
  NumericVector LAI_expanded = above["LAI_expanded"];
  NumericVector LAI_dead = above["LAI_dead"];
  NumericVector SA = above["SA"];
  NumericVector Cstorage = above["Cstorage"];
  int numCohorts = SP.size();

  //Aboveground parameters  
  List below = Rcpp::as<Rcpp::List>(x["below"]);
  NumericVector Z = Rcpp::as<Rcpp::NumericVector>(below["Z"]);
  
  //Base parameters
  DataFrame paramsBase = Rcpp::as<Rcpp::DataFrame>(x["paramsBase"]);
  NumericVector Sgdd = paramsBase["Sgdd"];
  
  NumericVector Water_FC = soil["Water_FC"];
  NumericVector W = soil["W"];
  
  int nlayers = W.size();
  
  //Growth parameters
  DataFrame paramsGrowth = Rcpp::as<Rcpp::DataFrame>(x["paramsGrowth"]);
  NumericVector SLA = Rcpp::as<Rcpp::NumericVector>(paramsGrowth["SLA"]);
  NumericVector Al2As = Rcpp::as<Rcpp::NumericVector>(paramsGrowth["Al2As"]);
  NumericVector WoodC = Rcpp::as<Rcpp::NumericVector>(paramsGrowth["WoodC"]);
  NumericVector WoodDens = Rcpp::as<Rcpp::NumericVector>(paramsGrowth["WoodDens"]);
  NumericVector Cstoragepmax = Rcpp::as<Rcpp::NumericVector>(paramsGrowth["Cstoragepmax"]);

  //Water balance output variables
  NumericMatrix PlantPsi(numDays, numCohorts);
  NumericMatrix PlantStress(numDays, numCohorts);
  NumericMatrix PlantTranspiration(numDays, numCohorts);
  NumericMatrix PlantRespiration(numDays, numCohorts);
  NumericMatrix PlantCstorage(numDays, numCohorts);
  NumericMatrix PlantSA(numDays, numCohorts);
  NumericMatrix PlantPhotosynthesis(numDays, numCohorts);
  NumericMatrix PlantLAIdead(numDays, numCohorts);
  NumericMatrix PlantLAIlive(numDays, numCohorts);
  
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
    NumericVector EplantCoh = s["EplantCoh"];
    NumericVector EplantVec = s["EplantVec"];
    NumericVector An =  Rcpp::as<Rcpp::NumericVector>(x["Photosynthesis"]);
    PlantPhotosynthesis(i,_) = An;
    PlantTranspiration(i,_) = EplantCoh;
    PlantStress(i,_) = Rcpp::as<Rcpp::NumericVector>(s["DDS"]);
    PlantPsi(i,_) = Rcpp::as<Rcpp::NumericVector>(s["psiCoh"]);
    PlantLAIdead(i,_) = LAI_dead;
      
    //3. Carbon balance and growth
    double B_leaf, B_stem, B_root;
    for(int j=0;j<numCohorts;j++){
      //3.1 Live biomass
      B_leaf = leafCperDry*1000.0*(LAI_expanded[j]/(N[j]/10000.0))/SLA[j]; //Biomass in g C · ind-1
      B_stem = WoodC[j]*SA[j]*(H[j]+Z[j])*WoodDens[j];
      B_root = B_leaf/2.5;
      //3.2 Respiration and photosynthesis 
      double Anj = An[j]/(N[j]/10000.0); //Translate g C · m-2 to g C · ind-1
      double Rj = dailyRespiration(B_leaf, B_stem, B_root, MeanTemperature[i]);
      PlantRespiration(i,j) = Rj*(N[j]/10000.0); //Scaled to cohort level
      //3.3. Carbon balance and C storage update
      double Cstorage_max = Cstoragepmax[j]*(B_leaf+B_stem+B_root);
      Cstorage[j] = std::min(Cstorage[j]+(Anj-Rj), Cstorage_max);
      PlantCstorage(i,j) = Cstorage[j];
      //3.4 Growth in LAI_live and SA
      double deltaSA = (dailySAturnoverProportion/(1.0+15*exp(-0.01*H[j])))*SA[j];
      SA[j] = SA[j] - deltaSA;
      LAI_live[j] = LAI_live[j] - (N[j]/10000.0)*(deltaSA/10000.0)*Al2As[j];
      PlantSA(i,j) = SA[j];
      PlantLAIlive(i,_) = LAI_live;
    }
    
  }
  if(verbose) Rcout << "done\n";
  
  NumericVector Eplanttot(numDays,0.0);


  if(verbose) {
    double Precipitationsum = sum(Precipitation);
    double Eplantsum = sum(Eplanttot);
    
    Rcout<<"Total Precipitation (mm) "  <<round(Precipitationsum) <<"\n";
    for(int l=0;l<nlayers;l++) Rcout << "W"<<(l+1)<<"f:"<< round(100*W[l])/100<<" ";
    Rcout<<"\n";

  }

  if(verbose) Rcout<<"plant output ...";
  
  PlantTranspiration.attr("dimnames") = List::create(meteo.attr("row.names"), SP.attr("names"));
  PlantStress.attr("dimnames") = List::create(meteo.attr("row.names"), SP.attr("names")) ;
  
  if(verbose) Rcout<<"list ...";
  List l = List::create(Named("TranspirationMode") = transpirationMode,
                        Named("NumSoilLayers") = nlayers,
                        Named("PlantTranspiration") = PlantTranspiration,
                        Named("PlantPhotosynthesis") = PlantPhotosynthesis,
                        Named("PlantRespiration") = PlantRespiration,
                        Named("PlantCstorage") = PlantCstorage,
                        Named("PlantSA")=PlantSA,
                        Named("PlantPsi") = PlantPsi, 
                        Named("PlantStress") = PlantStress,
                        Named("PlantLAIdead") = PlantLAIdead,
                        Named("PlantLAIlive") = PlantLAIlive);
  l.attr("class") = CharacterVector::create("growth","list");
  if(verbose) Rcout<<"done.\n";
  return(l);
}