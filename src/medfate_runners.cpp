#include <RcppArmadillo.h>
#include "medfate_runners.h"  
#include <string>
#include "spwb_day_c.h"
#include "growth_day_c.h"
#include "lowlevel_structures_c.h"
using namespace Rcpp;

//  Constructor for WB_runner
WB_runner::WB_runner(Rcpp::List x_list, 
                     double latitude, double elevation, double slope, double aspect) : 
  latitude(latitude),
  elevation(elevation),
  slope(slope),
  aspect(aspect),
  WBcomm(0,0,0,0) {
  Rcpp::CharacterVector classVector = x_list.attr("class");
  Rcpp::String s = classVector[0];
  std::string input_classIn = s.get_cstring();
  if(input_classIn=="spwbInput") {
    ModelInput x_m = ModelInput(x_list);
    x = std::make_unique<ModelInput>(x_m);
    int numCohorts = x_m.cohorts.SpeciesIndex.size();
    int nlayers = x_m.soil.getNlayers();
    int ncanlayers = x_m.canopy.zlow.size();
    int ntimesteps = x_m.control.advancedWB.ndailysteps;
    if(x_m.control.transpirationMode=="Granier") {
      BasicTranspiration_RESULT BTres(numCohorts, nlayers);
      BasicSPWB_RESULT BSPWBres(BTres);
      WBres = std::make_unique<BasicSPWB_RESULT>(BSPWBres);
    } else {
      AdvancedTranspiration_RESULT ATres(numCohorts, nlayers, ncanlayers, ntimesteps);
      AdvancedSPWB_RESULT ASPWBres(ATres);
      WBres = std::make_unique<AdvancedSPWB_RESULT>(ASPWBres);
    }
    WBcomm = WBCommunicationStructures(x_m.cohorts.SpeciesIndex.size(), x_m.soil.getNlayers(), x_m.canopy.zlow.size(), x_m.control.advancedWB.ndailysteps);
  } else if (input_classIn == "aspwbInput") {
    AgricultureModelInput x_m = AgricultureModelInput(x_list);
    x = std::make_unique<AgricultureModelInput>(x_m);
    int nlayers = x_m.soil.getNlayers();
    AgricultureWB_RESULT AgrWBres(nlayers);
    WBres = std::make_unique<AgricultureWB_RESULT>(AgrWBres);
    WBcomm = WBCommunicationStructures(0, nlayers,0, 0);
  } else {
    throw medfate::MedfateInternalError("Wrong model input class (should be spwbInput or aspwbInput)");
  }
}
WB_runner::~WB_runner(){
  // delete SPWBres;
}

void WB_runner::run_day(Rcpp::CharacterVector date, Rcpp::NumericVector meteovec, 
                        double runon, Rcpp::Nullable<Rcpp::NumericVector> lateralFlows, double waterTableDepth) {
  int nlayers = x->soil.getNlayers();
  std::vector<double> lateralFlows_c(nlayers, 0.0);
  NumericVector lateralFlows_mm;
  if(lateralFlows.isNotNull()) {
    lateralFlows_mm = NumericVector(lateralFlows);
    for(int l=0;l<nlayers;l++) {
      lateralFlows_c[l] = lateralFlows_mm[l];
    }
  }
  wb_day_inner_c(*WBres, WBcomm, *x,
                  Rcpp::as<std::string>(date[0]),
                  WeatherInputVector(meteovec),
                  latitude, elevation, slope, aspect,
                  runon,
                  lateralFlows_c, waterTableDepth);}

Rcpp::List WB_runner::get_output() {
  return(copyWBResult_c(*WBres, *x));
}
void WB_runner::update_input(List x_list) {
  x->copyStateToList(x_list);
}

//  Constructor for SPWB_multiple_runner
WB_multiple_runner::WB_multiple_runner(List wbInput_vec, 
                                       NumericVector latitude, 
                                       NumericVector elevation, 
                                       NumericVector slope, 
                                       NumericVector aspect) :
  n(wbInput_vec.size()), latitude_vec(n), topo_vec(n), x_vec(n), WBres_vec(n), WBcomm(0,0,0,0) {
  size_t numCohorts_max = 0;
  size_t nlayers_max = 0;
  size_t ncanlayers_max = 0;
  size_t ntimesteps_max = 0;
  for(int i=0; i<n;i++) {
    Rcpp::List wbInput_i = wbInput_vec[i];
    Rcpp::CharacterVector classVector_i = wbInput_i.attr("class");
    Rcpp::String s_i = classVector_i[0];
    std::string input_class_i = s_i.get_cstring();
    if(input_class_i=="spwbInput") {
      ModelInput x_i = ModelInput(wbInput_i);
      x_vec[i] = std::make_unique<ModelInput>(x_i);
      int numCohorts_i = x_i.cohorts.SpeciesIndex.size();
      int nlayers_i = x_i.soil.getNlayers();
      int ncanlayers_i = x_i.canopy.zlow.size();
      int ntimesteps_i = x_i.control.advancedWB.ndailysteps;
      if(x_i.control.transpirationMode=="Granier") {
        BasicTranspiration_RESULT BTres(numCohorts_i, nlayers_i);
        BasicSPWB_RESULT BSPWBres(BTres);
        WBres_vec[i] = std::make_unique<BasicSPWB_RESULT>(BSPWBres);
      } else {
        AdvancedTranspiration_RESULT ATres(numCohorts_i, nlayers_i, ncanlayers_i, ntimesteps_i);
        AdvancedSPWB_RESULT ASPWBres(ATres);
        WBres_vec[i] = std::make_unique<AdvancedSPWB_RESULT>(ASPWBres);
      }
      numCohorts_max = std::max(numCohorts_max, x_i.cohorts.SpeciesIndex.size());
      nlayers_max = std::max(nlayers_max, (size_t) x_i.soil.getNlayers());
      ncanlayers_max = std::max(ncanlayers_max, x_i.canopy.zlow.size());
      ntimesteps_max = std::max(ntimesteps_max, (size_t) x_i.control.advancedWB.ndailysteps);
    } else if (input_class_i == "aspwbInput") {
      AgricultureModelInput x_i = AgricultureModelInput(wbInput_i);
      x_vec[i] = std::make_unique<AgricultureModelInput>(x_i);
      int nlayers_i = x_i.soil.getNlayers();
      nlayers_max = std::max(nlayers_max, (size_t) nlayers_i);
      AgricultureWB_RESULT AgrWBres(nlayers_i);
      WBres_vec[i] = std::make_unique<AgricultureWB_RESULT>(AgrWBres);
    } else {
      throw medfate::MedfateInternalError("Wrong model input class (should be spwbInput or aspwbInput)");
    }
    latitude_vec[i] = latitude[i];
    Topography topo_i = Topography();
    topo_i.elevation = elevation[i];
    topo_i.slope = slope[i];
    topo_i.aspect = aspect[i];
    topo_vec[i] = std::make_unique<Topography>(topo_i);
  } 
  WBcomm = WBCommunicationStructures(numCohorts_max, nlayers_max, ncanlayers_max, ntimesteps_max);
}
WB_multiple_runner::~WB_multiple_runner() {
}

void WB_multiple_runner::run_day(Rcpp::CharacterVector date, Rcpp::List meteovec_list, bool parallelize) {

  std::vector<WeatherInputVector> weather_vec(x_vec.size());
  for(size_t i=0; i<x_vec.size();i++) {
    NumericVector meteovec_i = meteovec_list[i];
    weather_vec[i] = WeatherInputVector(meteovec_i);
  }
  std::string date_str = Rcpp::as<std::string>(date[0]);
  if(parallelize) {
    //build worker
    WATERBALANCE_worker worker(WBcomm,
                               date_str, 
                               x_vec,
                               latitude_vec,
                               topo_vec,
                               weather_vec,
                               WBres_vec);
    // call it with parallelFor
    parallelFor(0, x_vec.size(), worker, std::min(100, (int) x_vec.size()));
  } else {
    for(size_t i=0;i<x_vec.size();i++) {
      std::vector<double> lateralFlows_c(x_vec[i]->soil.getNlayers(), 0.0);
      double runon = 0.0;
      double waterTableDepth = medfate::NA_DOUBLE;
      wb_day_inner_c(*WBres_vec[i], WBcomm, *x_vec[i],
                      date_str,
                      weather_vec[i],
                      latitude_vec[i], topo_vec[i]->elevation, topo_vec[i]->slope, topo_vec[i]->aspect,
                      runon,
                      lateralFlows_c, waterTableDepth);
    }
  }
}

//Returs output (decreases index)
Rcpp::List WB_multiple_runner::get_output_at(int i) {
  return(copyWBResult_c(*WBres_vec[i-1], *x_vec[i-1]));
}
void WB_multiple_runner::update_input_at(int i, List x_list) {
  x_vec[i-1]->copyStateToList(x_list);
}



//  Constructor for WB_runner
GROWTH_runner::GROWTH_runner(Rcpp::List x_list, 
                             double latitude, double elevation, double slope, double aspect) : 
  latitude(latitude),
  elevation(elevation),
  slope(slope),
  aspect(aspect),
  GROWTHcomm(0,0,0,0) {
  Topography top;
  Rcpp::CharacterVector classVector = x_list.attr("class");
  Rcpp::String s = classVector[0];
  std::string input_classIn = s.get_cstring();
  if(input_classIn=="spwbInput") {
    ModelInput x_m = ModelInput(x_list);
    x = std::make_unique<ModelInput>(x_m);
    int numCohorts = x_m.cohorts.SpeciesIndex.size();
    int nlayers = x_m.soil.getNlayers();
    int ncanlayers = x_m.canopy.zlow.size();
    int ntimesteps = x_m.control.advancedWB.ndailysteps;
    if(x_m.control.transpirationMode=="Granier") {
      BasicTranspiration_RESULT BTres(numCohorts, nlayers);
      BasicSPWB_RESULT BSPWBres(BTres);
      BasicGROWTH_RESULT BGROWTHres(BSPWBres, numCohorts);
      GROWTHres = std::make_unique<BasicGROWTH_RESULT>(BGROWTHres);
    } else {
      AdvancedTranspiration_RESULT ATres(numCohorts, nlayers, ncanlayers, ntimesteps);
      AdvancedSPWB_RESULT ASPWBres(ATres);
      AdvancedGROWTH_RESULT AGROWTHres(ASPWBres, numCohorts, ntimesteps);
      GROWTHres = std::make_unique<AdvancedGROWTH_RESULT>(AGROWTHres);
    }
    GROWTHcomm = GROWTHCommunicationStructures(x_m.cohorts.SpeciesIndex.size(), x_m.soil.getNlayers(), x_m.canopy.zlow.size(), x_m.control.advancedWB.ndailysteps);
  }
}
GROWTH_runner::~GROWTH_runner(){}

void GROWTH_runner::run_day(Rcpp::CharacterVector date, Rcpp::NumericVector meteovec, 
                            double runon, Rcpp::Nullable<Rcpp::NumericVector> lateralFlows, double waterTableDepth) {
  int nlayers = x->soil.getNlayers();
  std::vector<double> lateralFlows_c(nlayers, 0.0);
  NumericVector lateralFlows_mm;
  if(lateralFlows.isNotNull()) {
    lateralFlows_mm = NumericVector(lateralFlows);
    for(int l=0;l<nlayers;l++) {
      lateralFlows_c[l] = lateralFlows_mm[l];
    }
  }
  growthDay_inner_c(*GROWTHres, GROWTHcomm, *x,
                 Rcpp::as<std::string>(date[0]),
                 WeatherInputVector(meteovec),
                 latitude, elevation, slope, aspect,
                 runon,
                 lateralFlows_c, waterTableDepth);}

Rcpp::List GROWTH_runner::get_output() {
  return(copyGROWTHResult_c(*GROWTHres, *x));
}
void GROWTH_runner::update_input(List x_list) {
  x->copyStateToList(x_list);
}



//  Constructor for GROWTH_multiple_runner
GROWTH_multiple_runner::GROWTH_multiple_runner(List wbInput_vec, 
                                       NumericVector latitude, 
                                       NumericVector elevation, 
                                       NumericVector slope, 
                                       NumericVector aspect) :
  n(wbInput_vec.size()), latitude_vec(n), topo_vec(n), x_vec(n), GROWTHres_vec(n), GROWTHcomm(0,0,0,0) {
  size_t numCohorts_max = 0;
  size_t nlayers_max = 0;
  size_t ncanlayers_max = 0;
  size_t ntimesteps_max = 0;
  for(int i=0; i<n;i++) {
    ModelInput x_i = ModelInput(wbInput_vec[i]);
    x_vec[i] = std::make_unique<ModelInput>(x_i);
    int numCohorts_i = x_i.cohorts.SpeciesIndex.size();
    int nlayers_i = x_i.soil.getNlayers();
    int ncanlayers_i = x_i.canopy.zlow.size();
    int ntimesteps_i = x_i.control.advancedWB.ndailysteps;
    if(x_i.control.transpirationMode=="Granier") {
      BasicTranspiration_RESULT BTres(numCohorts_i, nlayers_i);
      BasicSPWB_RESULT BSPWBres(BTres);
      BasicGROWTH_RESULT BGROWTHres(BSPWBres, numCohorts_i);
      GROWTHres_vec[i] = std::make_unique<BasicGROWTH_RESULT>(BGROWTHres);
    } else {
      AdvancedTranspiration_RESULT ATres(numCohorts_i, nlayers_i, ncanlayers_i, ntimesteps_i);
      AdvancedSPWB_RESULT ASPWBres(ATres);
      AdvancedGROWTH_RESULT AGROWTHres(ASPWBres, numCohorts_i, ntimesteps_i);
      GROWTHres_vec[i] = std::make_unique<AdvancedGROWTH_RESULT>(AGROWTHres);
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
  GROWTHcomm = GROWTHCommunicationStructures(numCohorts_max, nlayers_max, ncanlayers_max, ntimesteps_max);
}
GROWTH_multiple_runner::~GROWTH_multiple_runner() {
}

void GROWTH_multiple_runner::run_day(Rcpp::CharacterVector date, Rcpp::List meteovec_list, bool parallelize) {
  
  std::vector<WeatherInputVector> weather_vec(x_vec.size());
  for(size_t i=0; i<x_vec.size();i++) {
    NumericVector meteovec_i = meteovec_list[i];
    weather_vec[i] = WeatherInputVector(meteovec_i);
  }
  std::string date_str = Rcpp::as<std::string>(date[0]);
  if(parallelize) {
    //build worker
    GROWTH_worker worker(GROWTHcomm,
                               date_str, 
                               x_vec,
                               latitude_vec,
                               topo_vec,
                               weather_vec,
                               GROWTHres_vec);
    // call it with parallelFor
    parallelFor(0, x_vec.size(), worker, std::min(100, (int) x_vec.size()));
  } else {
    for(size_t i=0;i<x_vec.size();i++) {
      std::vector<double> lateralFlows_c(x_vec[i]->soil.getNlayers(), 0.0);
      double runon = 0.0;
      double waterTableDepth = medfate::NA_DOUBLE;
      growthDay_inner_c(*GROWTHres_vec[i], GROWTHcomm, *x_vec[i],
                        date_str,
                        weather_vec[i],
                        latitude_vec[i], topo_vec[i]->elevation, topo_vec[i]->slope, topo_vec[i]->aspect,
                        runon,
                        lateralFlows_c, waterTableDepth);
    }
  }
}

//Returns output (decreases index)
Rcpp::List GROWTH_multiple_runner::get_output_at(int i) {
  return(copyGROWTHResult_c(*GROWTHres_vec[i-1], *x_vec[i-1]));
}
void GROWTH_multiple_runner::update_input_at(int i, List x_list) {
  x_vec[i-1]->copyStateToList(x_list);
}

/* 
 * CREATE RCPP MODULES 
 */

RCPP_MODULE(mod_wb) {
  class_<WB_runner>( "WB_runner" )
  .constructor<Rcpp::List, double, double, double, double>()
  .method( "run_day", &WB_runner::run_day )
  .method( "get_output", &WB_runner::get_output)
  .method( "update_input", &WB_runner::update_input)
  ;
}
RCPP_MODULE(mod_multiple_wb) {
  class_<WB_multiple_runner>( "WB_multiple_runner" )
  .constructor<Rcpp::List, NumericVector, NumericVector, NumericVector, NumericVector>()
  .method( "run_day", &WB_multiple_runner::run_day )
  .method( "get_output_at", &WB_multiple_runner::get_output_at)
  .method( "update_input_at", &WB_multiple_runner::update_input_at)
  ;
}
RCPP_MODULE(mod_growth) {
  class_<GROWTH_runner>( "GROWTH_runner" )
  .constructor<Rcpp::List, double, double, double, double>()
  .method( "run_day", &GROWTH_runner::run_day )
  .method( "get_output", &GROWTH_runner::get_output)
  .method( "update_input", &GROWTH_runner::update_input)
  ;
}
RCPP_MODULE(mod_multiple_growth) {
  class_<GROWTH_multiple_runner>( "GROWTH_multiple_runner" )
  .constructor<Rcpp::List, NumericVector, NumericVector, NumericVector, NumericVector>()
  .method( "run_day", &GROWTH_multiple_runner::run_day )
  .method( "get_output_at", &GROWTH_multiple_runner::get_output_at)
  .method( "update_input_at", &GROWTH_multiple_runner::update_input_at)
  ;
}
