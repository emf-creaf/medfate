#include <RcppArmadillo.h>
#include "medfate_runners.h"  
#include <string>
#include "spwb_day_c.h"
#include "lowlevel_structures_c.h"
using namespace Rcpp;

//  SPWB_runner

SPWB_runner::SPWB_runner(Rcpp::List x_list, 
                             double latitude, double elevation, double slope, double aspect) : 
  x(x_list),
  latitude(latitude), 
  elevation(elevation),
  slope(slope),
  aspect(aspect),
  SPWBcomm(x.cohorts.SpeciesIndex.size(), x.soil.getNlayers(), x.canopy.zlow.size(), x.control.advancedWB.ndailysteps) {
}

void SPWB_runner::run_day(Rcpp::CharacterVector date, Rcpp::NumericVector meteovec, 
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

Rcpp::List SPWB_runner::get_output() {
  List modelOutput;
  if(x.control.transpirationMode=="Granier") {
    modelOutput = copyBasicSPWBResult_c(SPWBcomm.BSPWBres, x);
  } else {
    modelOutput = copyAdvancedSPWBResult_c(SPWBcomm.ASPWBres, x);
  }
  return(modelOutput);
}

RCPP_MODULE(mod_spwb) {
  class_<SPWB_runner>( "SPWB_runner" )
  .constructor<Rcpp::List, double, double, double, double>()
  .method( "run_day", &SPWB_runner::run_day )
  .method( "get_output", &SPWB_runner::get_output)
  ;
}

//  SPWB_multiple_runner
SPWB_multiple_runner::SPWB_multiple_runner(List spwbInput_vec, 
                                           NumericVector latitude, 
                                           NumericVector elevation, 
                                           NumericVector slope, 
                                           NumericVector aspect) :
  x_vec(spwbInput_vec.size()), 
  latitude_vec(spwbInput_vec.size()), 
  topo_vec(spwbInput_vec.size()) {
    for(int i=0; i<spwbInput_vec.size();i++) {
      x_vec[i] = ModelInput(spwbInput_vec[i]);
      latitude_vec[i] = latitude[i];
      topo_vec[i] = Topography();
      topo_vec[i].elevation = elevation[i];
      topo_vec[i].slope = slope[i];
      topo_vec[i].aspect = aspect[i];
    }
}

void SPWB_multiple_runner::run_day(Rcpp::CharacterVector date, Rcpp::List meteovec_list) {
  std::vector<WeatherInputVector> weather_vec(x_vec.size());
  for(int i=0; i<x_vec.size();i++) {
    NumericVector meteovec_i = meteovec_list[i];
    weather_vec[i] = WeatherInputVector(meteovec_i);
  }
  int ncohorts_max = 0;
  int nlayers_max = 0;
  int ncanlayers_max = 0;
  int ntimesteps_max = 0;
  SPWBCommunicationStructures SPWBcomm(ncohorts_max, nlayers_max, ncanlayers_max, ntimesteps_max);
  std::string date_str = Rcpp::as<std::string>(date[0]);
  SPWB_worker worker(SPWBcomm,
                     date_str, 
                     x_vec,
                     latitude_vec,
                     topo_vec,
                     weather_vec,
                     SPWB_res_vec);
  
  // call it with parallelFor
  parallelFor(0, x_vec.size(), worker);
  
  // do I have to copy back the result
}

Rcpp::List SPWB_multiple_runner::get_output_at(int i) {
  return(copySPWBResult_c(SPWB_res_vec[i], x_vec[i]));
}

RCPP_MODULE(mod_multiple_spwb) {
  class_<SPWB_multiple_runner>( "SPWB_multiple_runner" )
  .constructor<Rcpp::List, NumericVector, NumericVector, NumericVector, NumericVector>()
  .method( "run_day", &SPWB_multiple_runner::run_day )
  .method( "get_output_at", &SPWB_multiple_runner::get_output_at)
  ;
}