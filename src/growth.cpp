#include <numeric>
#include "lightextinction.h"
#include "windextinction.h"
#include "hydraulics.h"
#include "biophysicsutils.h"
#include "photosynthesis.h"
#include "soil.h"
#include "swb.h"
#include <Rcpp.h>
#include <meteoland.h>
using namespace Rcpp;

void checkgrowthInput(List x, List soil, String transpirationMode) {
  if(!x.containsElementNamed("plants")) stop("plants missing in growthInput");
  DataFrame plants = Rcpp::as<Rcpp::DataFrame>(x["plants"]);
  if(!plants.containsElementNamed("LAI_live")) stop("LAI_live missing in growthInput$plants");
  if(!plants.containsElementNamed("LAI_dead")) stop("LAI_dead missing in growthInput$plants");
  if(!plants.containsElementNamed("LA_live")) stop("LA_live missing in growthInput$plants");
  if(!plants.containsElementNamed("LA_dead")) stop("LA_dead missing in growthInput$plants");
  if(!plants.containsElementNamed("LA_predrought")) stop("LA_predrought missing in growthInput$plants");
  if(!plants.containsElementNamed("SA")) stop("SA missing in growthInput$plants");
  if(!plants.containsElementNamed("Cstorage")) stop("Cstorage missing in growthInput$plants");
  if(!plants.containsElementNamed("CR")) stop("CR missing in growthInput$plants");
  if(!plants.containsElementNamed("H")) stop("H missing in growthInput$plants");
  if(!plants.containsElementNamed("Z")) stop("Z missing in growthInput$plants");
  if(!plants.containsElementNamed("N")) stop("N missing in growthInput$plants");
  if(!plants.containsElementNamed("DBH")) stop("DBH missing in growthInput$plants");
  
  if(!x.containsElementNamed("V")) stop("V missing in growthInput");

  if(!x.containsElementNamed("params")) stop("params missing in growthInput");
  DataFrame params = Rcpp::as<Rcpp::DataFrame>(x["params"]);
  if(!params.containsElementNamed("Sgdd")) stop("Sgdd missing in growthInput$params");
  if(!params.containsElementNamed("k")) stop("k missing in growthInput$params");
  if(!params.containsElementNamed("g")) stop("g missing in growthInput$params");
  if(!params.containsElementNamed("SLA")) stop("SLA missing in growthInput$params");
  
  if(!x.containsElementNamed("paramsTransp")) stop("paramsTransp missing in growthInput");
  DataFrame paramsTransp = Rcpp::as<Rcpp::DataFrame>(x["paramsTransp"]);
  if(transpirationMode=="Simple") {
    if(!paramsTransp.containsElementNamed("Psi_Extract")) stop("Psi_Extract missing in growthInput$paramsTransp");
    if(!paramsTransp.containsElementNamed("WUE")) stop("WUE missing in growthInput$paramsTransp");
  } else if(transpirationMode=="Sperry") {
    if(!soil.containsElementNamed("VG_n")) stop("VG_n missing in soil");
    if(!soil.containsElementNamed("VG_alpha")) stop("VG_alpha missing in soil");
    if(!x.containsElementNamed("VGrhizo_kmax")) stop("VCstem_kmax missing in growthInput");
    if(!x.containsElementNamed("VCroot_kmax")) stop("VCroot_kmax missing in growthInput");
    
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
  }
  
  CharacterVector dateStrings = meteo.attr("row.names");
  
  NumericVector GDD = gdd(DOY, MeanTemperature, 5.0);
  NumericVector ER = er(DOY);
  
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  NumericVector SP = above["SP"];
  int numCohorts = SP.size();
  
  NumericVector Water_FC = soil["Water_FC"];
  NumericVector W = soil["W"];
  
  int nlayers = W.size();
  
  //Growth parameters
  DataFrame paramsGrowth = Rcpp::as<Rcpp::DataFrame>(x["paramsGrowth"]);
  NumericVector SLA = Rcpp::as<Rcpp::NumericVector>(paramsGrowth["SLA"]);
  NumericVector Al2As = Rcpp::as<Rcpp::NumericVector>(paramsGrowth["Al2As"]);
  
  //Water balance output variables
  NumericMatrix PlantPsi(numDays, numCohorts);
  NumericMatrix PlantStress(numDays, numCohorts);
  NumericMatrix PlantTranspiration(numDays, numCohorts);
  NumericMatrix PlantPhotosynthesis(numDays, numCohorts);

  if(verbose) {
    for(int l=0;l<nlayers;l++) Rcout << "W"<<(l+1)<<"i:"<< round(100*W[l])/100<<" ";
    Rcout<<"\n";
  }
  
  if(verbose) Rcout << "Daily growth:";
  List s;
  for(int i=0;i<numDays;i++) {
    if(verbose) Rcout<<".";
    //Phenology
    
    //Water balance and photosynthesis
    if(transpirationMode=="Simple") {
      s = swbDay1(x, soil, GDD[i], MeanTemperature[i], PET[i], Precipitation[i], ER[i], 0.0, false); //No Runon in simulations for a single cell
    } else if(transpirationMode=="Sperry") {
      std::string c = as<std::string>(dateStrings[i]);
      int J = meteoland::radiation_julianDay(std::atoi(c.substr(0, 4).c_str()),std::atoi(c.substr(5,2).c_str()),std::atoi(c.substr(8,2).c_str()));
      double delta = meteoland::radiation_solarDeclination(J);
      s = swbDay2(x, soil, GDD[i], MinTemperature[i], MaxTemperature[i], 
                       MinRelativeHumidity[i], MaxRelativeHumidity[i], Radiation[i], WindSpeed[i], 
                       latitude, elevation, slope, aspect, delta, Precipitation[i], ER[i], 0.0, false);
      
      PET[i] = s["PET"];
    }    
    NumericVector EplantCoh = s["EplantCoh"];
    NumericVector EplantVec = s["EplantVec"];
    PlantPhotosynthesis(i,_) = Rcpp::as<Rcpp::NumericVector>(x["Photosynthesis"]);
    PlantTranspiration(i,_) = EplantCoh;
    PlantStress(i,_) = Rcpp::as<Rcpp::NumericVector>(s["DDS"]);
    PlantPsi(i,_) = Rcpp::as<Rcpp::NumericVector>(s["psiCoh"]);
    
    //3. Carbon balance and growth
    double B_leaf, B_stem, B_root;
    for(int j=0;j<numCohorts;j++){
      //3.1 Live biomass and maximum C storage
      
      //3.2 Respiration
      
      //3.3. Carbon balance and C storage update
      
      //3.4 Growth in LAI_live and SA
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
                        Named("PlantPsi") = PlantPsi, 
                        Named("PlantStress") = PlantStress);
  l.attr("class") = CharacterVector::create("growth","list");
  if(verbose) Rcout<<"done.\n";
  return(l);
}