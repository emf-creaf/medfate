#include <RcppArmadillo.h>
#include "medfate_runners.h"  
#include <string>
#include "spwb_day_c.h"
#include "lowlevel_structures_c.h"
using namespace Rcpp;

//  Constructor for SPWB_runner
SPWB_runner::SPWB_runner(Rcpp::List x_list, 
                             double latitude, double elevation, double slope, double aspect) : 
  x(x_list),
  latitude(latitude), 
  elevation(elevation),
  slope(slope),
  aspect(aspect),
  SPWBcomm(x.cohorts.SpeciesIndex.size(), x.soil.getNlayers(), x.canopy.zlow.size(), x.control.advancedWB.ndailysteps) {
  int numCohorts = x.cohorts.SpeciesIndex.size();
  int nlayers = x.soil.getNlayers();
  int ncanlayers = x.canopy.zlow.size();
  int ntimesteps = x.control.advancedWB.ndailysteps;
  if(x.control.transpirationMode=="Granier") {
    BasicTranspiration_RESULT BTres(numCohorts, nlayers);
    BasicSPWB_RESULT BSPWBres(BTres);
    SPWBres = std::make_unique<BasicSPWB_RESULT>(BSPWBres);
  } else {
    AdvancedTranspiration_RESULT ATres(numCohorts, nlayers, ncanlayers, ntimesteps);
    AdvancedSPWB_RESULT ASPWBres(ATres);
    SPWBres = std::make_unique<AdvancedSPWB_RESULT>(ASPWBres);
  }
}
SPWB_runner::~SPWB_runner(){
  // delete SPWBres;
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
  spwbDay_inner_c(*SPWBres, SPWBcomm, x,
                  Rcpp::as<std::string>(date[0]),
                  WeatherInputVector(meteovec),
                  latitude, elevation, slope, aspect,
                  runon,
                  lateralFlows_c, waterTableDepth);}

Rcpp::List SPWB_runner::get_output() {
  return(copySPWBResult_c(*SPWBres, x));
}

RCPP_MODULE(mod_spwb) {
  class_<SPWB_runner>( "SPWB_runner" )
  .constructor<Rcpp::List, double, double, double, double>()
  .method( "run_day", &SPWB_runner::run_day )
  .method( "get_output", &SPWB_runner::get_output)
  ;
}

//  Constructor for SPWB_multiple_runner
SPWB_multiple_runner::SPWB_multiple_runner(List spwbInput_vec, 
                                           NumericVector latitude, 
                                           NumericVector elevation, 
                                           NumericVector slope, 
                                           NumericVector aspect) :
  n(spwbInput_vec.size()), latitude_vec(n), topo_vec(n), x_vec(n), SPWB_res_vec(n), SPWBcomm(0,0,0,0) {
  size_t numCohorts_max = 0;
  size_t nlayers_max = 0;
  size_t ncanlayers_max = 0;
  size_t ntimesteps_max = 0;
  for(int i=0; i<n;i++) {
    ModelInput x_i = ModelInput(spwbInput_vec[i]);
    x_vec[i] = std::make_unique<ModelInput>(x_i);
    int numCohorts_i = x_i.cohorts.SpeciesIndex.size();
    int nlayers_i = x_i.soil.getNlayers();
    int ncanlayers_i = x_i.canopy.zlow.size();
    int ntimesteps_i = x_i.control.advancedWB.ndailysteps;
    if(x_i.control.transpirationMode=="Granier") {
      BasicTranspiration_RESULT BTres(numCohorts_i, nlayers_i);
      BasicSPWB_RESULT BSPWBres(BTres);
      SPWB_res_vec[i] = std::make_unique<BasicSPWB_RESULT>(BSPWBres);
    } else {
      AdvancedTranspiration_RESULT ATres(numCohorts_i, nlayers_i, ncanlayers_i, ntimesteps_i);
      AdvancedSPWB_RESULT ASPWBres(ATres);
      SPWB_res_vec[i] = std::make_unique<AdvancedSPWB_RESULT>(ASPWBres);
    }
    numCohorts_max = std::max(numCohorts_max, x_i.cohorts.SpeciesIndex.size());
    nlayers_max = std::max(nlayers_max, (size_t) x_i.soil.getNlayers());
    ncanlayers_max = std::max(ncanlayers_max, x_i.canopy.zlow.size());
    ntimesteps_max = std::max(ntimesteps_max, (size_t) x_i.control.advancedWB.ndailysteps);
    latitude_vec[i] = latitude[i];
    Topography topo_i = Topography();
    topo_i.elevation = elevation[i];
    topo_i.slope = slope[i];
    topo_i.aspect = aspect[i];
    topo_vec[i] = std::make_unique<Topography>(topo_i);
  }
  SPWBcomm = SPWBCommunicationStructures(numCohorts_max, nlayers_max, ncanlayers_max, ntimesteps_max);
}
SPWB_multiple_runner::~SPWB_multiple_runner() {
}

void SPWB_multiple_runner::run_day(Rcpp::CharacterVector date, Rcpp::List meteovec_list, bool parallelize) {

  std::vector<WeatherInputVector> weather_vec(x_vec.size());
  for(size_t i=0; i<x_vec.size();i++) {
    NumericVector meteovec_i = meteovec_list[i];
    weather_vec[i] = WeatherInputVector(meteovec_i);
  }
  std::string date_str = Rcpp::as<std::string>(date[0]);
  if(parallelize) {
    //build worker
    SPWB_worker worker(SPWBcomm,
                       date_str, 
                       x_vec,
                       latitude_vec,
                       topo_vec,
                       weather_vec,
                       SPWB_res_vec);
    // call it with parallelFor
    parallelFor(0, x_vec.size(), worker, std::min(100, (int) x_vec.size()));
  } else {
    for(size_t i=0;i<x_vec.size();i++) {
      std::vector<double> lateralFlows_c(x_vec[i]->soil.getNlayers(), 0.0);
      double runon = 0.0;
      double waterTableDepth = medfate::NA_DOUBLE;
      spwbDay_inner_c(*SPWB_res_vec[i], SPWBcomm, *x_vec[i],
                      date_str,
                      weather_vec[i],
                      latitude_vec[i], topo_vec[i]->elevation, topo_vec[i]->slope, topo_vec[i]->aspect,
                      runon,
                      lateralFlows_c, waterTableDepth);
    }
  }
}

//Returs output (decreases index)
Rcpp::List SPWB_multiple_runner::get_output_at(int i) {
  return(copySPWBResult_c(*SPWB_res_vec[i-1], *x_vec[i-1]));
}

RCPP_MODULE(mod_multiple_spwb) {
  class_<SPWB_multiple_runner>( "SPWB_multiple_runner" )
  .constructor<Rcpp::List, NumericVector, NumericVector, NumericVector, NumericVector>()
  .method( "run_day", &SPWB_multiple_runner::run_day )
  .method( "get_output_at", &SPWB_multiple_runner::get_output_at)
  ;
}