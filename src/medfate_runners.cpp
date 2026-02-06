#include <RcppArmadillo.h>
#include "medfate_runners.h"  
#include <string>
#include "spwb_day_c.h"
#include "lowlevel_structures_c.h"
using namespace Rcpp;

SPWBrunner::SPWBrunner(Rcpp::List x_list, 
                             double latitude, double elevation, double slope, double aspect) : 
  x(x_list),
  latitude(latitude), 
  elevation(elevation),
  slope(slope),
  aspect(aspect),
  SPWBcomm(x.cohorts.SpeciesIndex.size(), x.soil.getNlayers(), x.canopy.zlow.size(), x.control.advancedWB.ndailysteps,
            x.control.soilDomains) {
}

void SPWBrunner::run_day(Rcpp::CharacterVector date, Rcpp::NumericVector meteovec, 
                            double runon, Rcpp::Nullable<Rcpp::NumericVector> lateralFlows, double waterTableDepth) {
  int nlayers = x.soil.getNlayers();
  std::vector<double> lateralFlows_c(nlayers, 0.0);
  NumericVector lateralFlows_mm;
  if(lateralFlows.isNotNull()) {
    lateralFlows_mm = NumericVector(lateralFlows);
    for(int l=0;l<nlayers;l++) {
      lateralFlows_c[l] = lateralFlows_mm[l];
    }
  }
  spwbDay_inner_c(SPWBcomm, x, 
                  Rcpp::as<std::string>(date[0]),
                  WeatherInputVector(meteovec), 
                  latitude, elevation, slope, aspect,
                  runon, 
                  lateralFlows_c, waterTableDepth);
}

Rcpp::List SPWBrunner::get_output() {
  List modelOutput;
  if(x.control.transpirationMode=="Granier") {
    modelOutput = copyBasicSPWBResult_c(SPWBcomm.BSPWBres, x);
  } else {
    modelOutput = copyAdvancedSPWBResult_c(SPWBcomm.ASPWBres, x);
  }
  return(modelOutput);
}

RCPP_MODULE(mod_spwb) {
  class_<SPWBrunner>( "SPWBrunner" )
  .constructor<Rcpp::List, double, double, double, double>()
  .method( "run_day", &SPWBrunner::run_day )
  .method( "get_output", &SPWBrunner::get_output)
  ;
}