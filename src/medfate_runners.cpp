#include <RcppArmadillo.h>
#include "medfate_runners.h"  
#include <string>
#include "spwb_day_c.h"
#include "growth_day_c.h"
#include "lowlevel_structures_c.h"
using namespace Rcpp;

//  Constructor for single_runner
single_runner::single_runner(Rcpp::List x_list, 
                             double latitude, double elevation, double slope, double aspect) : 
  latitude(latitude),
  elevation(elevation),
  slope(slope),
  aspect(aspect) {
  Rcpp::CharacterVector classVector = x_list.attr("class");
  Rcpp::String s = classVector[0];
  std::string input_classIn = s.get_cstring();
  // Rcpp::Rcout << input_classIn <<"\n";
  if(input_classIn=="spwbInput") {
    ModelInput x = ModelInput(x_list);
    p_x = std::make_unique<ModelInput>(x);
    int numCohorts = x.cohorts.SpeciesIndex.size();
    int nlayers = x.soil.getNlayers();
    int ncanlayers = x.canopy.zlow.size();
    int ntimesteps = x.control.advancedWB.ndailysteps;
    if(x.control.transpirationMode=="Granier") {
      BasicTranspiration_RESULT BTres(numCohorts, nlayers);
      BasicSPWB_RESULT BSPWBres(BTres);
      p_result = std::make_unique<BasicSPWB_RESULT>(BSPWBres);
    } else {
      AdvancedTranspiration_RESULT ATres(numCohorts, nlayers, ncanlayers, ntimesteps);
      AdvancedSPWB_RESULT ASPWBres(ATres);
      p_result = std::make_unique<AdvancedSPWB_RESULT>(ASPWBres);
    }
    WBCommunicationStructures WBcomm = WBCommunicationStructures(numCohorts, nlayers, ncanlayers, ntimesteps);
    p_WBcomm = std::make_unique<WBCommunicationStructures>(WBcomm);
  } else if (input_classIn == "growthInput") {
    ModelInput x = ModelInput(x_list);
    p_x = std::make_unique<ModelInput>(x);
    int numCohorts = x.cohorts.SpeciesIndex.size();
    int nlayers = x.soil.getNlayers();
    int ncanlayers = x.canopy.zlow.size();
    int ntimesteps = x.control.advancedWB.ndailysteps;
    if(x.control.transpirationMode=="Granier") {
      BasicTranspiration_RESULT BTres(numCohorts, nlayers);
      BasicSPWB_RESULT BSPWBres(BTres);
      BasicGROWTH_RESULT BGROWTHres(BSPWBres, numCohorts);
      p_result = std::make_unique<BasicGROWTH_RESULT>(BGROWTHres);
    } else if(x.control.transpirationMode=="Sperry" || x.control.transpirationMode=="Sureau")  {
      AdvancedTranspiration_RESULT ATres(numCohorts, nlayers, ncanlayers, ntimesteps);
      AdvancedSPWB_RESULT ASPWBres(ATres);
      AdvancedGROWTH_RESULT AGROWTHres(ASPWBres, numCohorts, ntimesteps);
      p_result = std::make_unique<AdvancedGROWTH_RESULT>(AGROWTHres);
    }
    GROWTHCommunicationStructures GROWTHcomm = GROWTHCommunicationStructures(numCohorts, nlayers, ncanlayers, ntimesteps);
    p_GROWTHcomm = std::make_unique<GROWTHCommunicationStructures>(GROWTHcomm);
  } else if (input_classIn == "aspwbInput") {
    AgricultureModelInput x = AgricultureModelInput(x_list);
    p_x = std::make_unique<AgricultureModelInput>(x);
    int nlayers = x.soil.getNlayers();
    AgricultureWB_RESULT AgrWBres(nlayers);
    p_result = std::make_unique<AgricultureWB_RESULT>(AgrWBres);
    WBCommunicationStructures WBcomm = WBCommunicationStructures(0, nlayers,0, 0);
    p_WBcomm = std::make_unique<WBCommunicationStructures>(WBcomm);
  } else {
    throw medfate::MedfateInternalError("Wrong model input class (should be spwbInput or aspwbInput)");
  }
}
single_runner::~single_runner(){
  // delete SPWBres;
}

// Performs simulation
void single_runner::run_day(Rcpp::CharacterVector date, Rcpp::NumericVector meteovec, 
                            double runon, Rcpp::Nullable<Rcpp::NumericVector> lateralFlows, double waterTableDepth) {
  if((p_x->getInputClass() == "spwbInput") || (p_x->getInputClass() == "aspwbInput")) {
    WaterBalanceModelInput& x = dynamic_cast<WaterBalanceModelInput&>(*p_x);
    WB_RESULT& WBres = dynamic_cast<WB_RESULT&>(*p_result);
    int nlayers = x.soil.getNlayers();
    std::vector<double> lateralFlows_c(nlayers, 0.0);
    NumericVector lateralFlows_mm;
    if(lateralFlows.isNotNull()) {
      lateralFlows_mm = NumericVector(lateralFlows);
      for(int l=0;l<nlayers;l++) {
        lateralFlows_c[l] = lateralFlows_mm[l];
      }
    }
    wb_day_inner_c(WBres, *p_WBcomm, x,
                   Rcpp::as<std::string>(date[0]),
                   WeatherInputVector(meteovec),
                   latitude, elevation, slope, aspect,
                   runon,
                   lateralFlows_c, waterTableDepth);
  } else if(p_x->getInputClass() == "growthInput") {
    ModelInput& x = dynamic_cast<ModelInput&>(*p_x);
    GROWTH_RESULT& GROWTHres = dynamic_cast<GROWTH_RESULT&>(*p_result);
    int nlayers = x.soil.getNlayers();
    std::vector<double> lateralFlows_c(nlayers, 0.0);
    NumericVector lateralFlows_mm;
    if(lateralFlows.isNotNull()) {
      lateralFlows_mm = NumericVector(lateralFlows);
      for(int l=0;l<nlayers;l++) {
        lateralFlows_c[l] = lateralFlows_mm[l];
      }
    }
    growthDay_inner_c(GROWTHres, *p_GROWTHcomm, x,
                      Rcpp::as<std::string>(date[0]),
                      WeatherInputVector(meteovec),
                      latitude, elevation, slope, aspect,
                      runon,
                      lateralFlows_c, waterTableDepth);
  } else {
    throw medfate::MedfateInternalError("Wrong input class for simulation launching");
  }
}

// Retrieves output
Rcpp::List single_runner::get_output() {
  List l;
  if(p_x->getInputClass() == "spwbInput") {
    ModelInput& x_m = dynamic_cast<ModelInput&>(*p_x);
    WB_RESULT& WBres = dynamic_cast<WB_RESULT&>(*p_result);
    l = copyWBResult_c(WBres, x_m);
  } else if(p_x->getInputClass() == "growthInput") {
    ModelInput& x_m = dynamic_cast<ModelInput&>(*p_x);
    GROWTH_RESULT& GROWTHres = dynamic_cast<GROWTH_RESULT&>(*p_result);
    l = copyGROWTHResult_c(GROWTHres, x_m);
  } else if(p_x->getInputClass() == "aspwbInput") {
    AgricultureModelInput& x_m = dynamic_cast<AgricultureModelInput&>(*p_x);
    AgricultureWB_RESULT& WBres = dynamic_cast<AgricultureWB_RESULT&>(*p_result);
    l = copyWBResult_c(WBres, x_m);
  } else {
    throw medfate::MedfateInternalError("Wrong input class for result copying");
  }
  return(l);
}

// Updates (external) input object
void single_runner::update_input(List x_list) {
  if(p_x->getInputClass()=="spwbInput" || p_x->getInputClass()=="growthInput") {
    ModelInput& x = dynamic_cast<ModelInput&>(*p_x);
    x.copyStateToList(x_list);
  } else if(p_x->getInputClass()=="aspwbInput") {
    AgricultureModelInput& x = dynamic_cast<AgricultureModelInput&>(*p_x);
    x.copyStateToList(x_list);
  }
}

//  Constructor for multiple_runner
multiple_runner::multiple_runner(List input_vec, 
                                 NumericVector latitude, 
                                 NumericVector elevation, 
                                 NumericVector slope, 
                                 NumericVector aspect) :
  n(input_vec.size()), latitude_vec(n), p_topo_vec(n), p_x_vec(n), p_result_vec(n) {
  size_t numCohorts_max = 0;
  size_t nlayers_max = 0;
  size_t ncanlayers_max = 0;
  size_t ntimesteps_max = 0;
  for(int i=0; i< n;i++) {
    Rcpp::List input_i = input_vec[i];
    Rcpp::CharacterVector classVector_i = input_i.attr("class");
    Rcpp::String s_i = classVector_i[0];
    std::string input_class_i = s_i.get_cstring();
    if(input_class_i=="spwbInput") {
      ModelInput x_i = ModelInput(input_i);
      p_x_vec[i] = std::make_unique<ModelInput>(x_i);
      int numCohorts_i = x_i.cohorts.SpeciesIndex.size();
      int nlayers_i = x_i.soil.getNlayers();
      int ncanlayers_i = x_i.canopy.zlow.size();
      int ntimesteps_i = x_i.control.advancedWB.ndailysteps;
      if(x_i.control.transpirationMode=="Granier") {
        BasicTranspiration_RESULT BTres(numCohorts_i, nlayers_i);
        BasicSPWB_RESULT BSPWBres(BTres);
        p_result_vec[i] = std::make_unique<BasicSPWB_RESULT>(BSPWBres);
      } else {
        AdvancedTranspiration_RESULT ATres(numCohorts_i, nlayers_i, ncanlayers_i, ntimesteps_i);
        AdvancedSPWB_RESULT ASPWBres(ATres);
        p_result_vec[i] = std::make_unique<AdvancedSPWB_RESULT>(ASPWBres);
      }
      numCohorts_max = std::max(numCohorts_max, x_i.cohorts.SpeciesIndex.size());
      nlayers_max = std::max(nlayers_max, (size_t) x_i.soil.getNlayers());
      ncanlayers_max = std::max(ncanlayers_max, x_i.canopy.zlow.size());
      ntimesteps_max = std::max(ntimesteps_max, (size_t) x_i.control.advancedWB.ndailysteps);
    } else if (input_class_i == "growthInput") {
      ModelInput x_i = ModelInput(input_i);
      p_x_vec[i] = std::make_unique<ModelInput>(x_i);
      int numCohorts_i = x_i.cohorts.SpeciesIndex.size();
      int nlayers_i = x_i.soil.getNlayers();
      int ncanlayers_i = x_i.canopy.zlow.size();
      int ntimesteps_i = x_i.control.advancedWB.ndailysteps;
      if(x_i.control.transpirationMode=="Granier") {
        BasicTranspiration_RESULT BTres(numCohorts_i, nlayers_i);
        BasicSPWB_RESULT BSPWBres(BTres);
        BasicGROWTH_RESULT BGROWTHres(BSPWBres, numCohorts_i);
        p_result_vec[i] = std::make_unique<BasicGROWTH_RESULT>(BGROWTHres);
      } else {
        AdvancedTranspiration_RESULT ATres(numCohorts_i, nlayers_i, ncanlayers_i, ntimesteps_i);
        AdvancedSPWB_RESULT ASPWBres(ATres);
        AdvancedGROWTH_RESULT AGROWTHres(ASPWBres, numCohorts_i, ntimesteps_i);
        p_result_vec[i] = std::make_unique<AdvancedGROWTH_RESULT>(AGROWTHres);
      }
      numCohorts_max = std::max(numCohorts_max, x_i.cohorts.SpeciesIndex.size());
      nlayers_max = std::max(nlayers_max, (size_t) x_i.soil.getNlayers());
      ncanlayers_max = std::max(ncanlayers_max, x_i.canopy.zlow.size());
      ntimesteps_max = std::max(ntimesteps_max, (size_t) x_i.control.advancedWB.ndailysteps);
    } else if (input_class_i == "aspwbInput") {
      AgricultureModelInput x_i = AgricultureModelInput(input_i);
      p_x_vec[i] = std::make_unique<AgricultureModelInput>(x_i);
      int nlayers_i = x_i.soil.getNlayers();
      nlayers_max = std::max(nlayers_max, (size_t) nlayers_i);
      AgricultureWB_RESULT AgrWBres(nlayers_i);
      p_result_vec[i] = std::make_unique<AgricultureWB_RESULT>(AgrWBres);
    } else {
      throw medfate::MedfateInternalError("Wrong model input class (should be spwbInput or aspwbInput)");
    }
    latitude_vec[i] = latitude[i];
    Topography topo_i = Topography();
    topo_i.elevation = elevation[i];
    topo_i.slope = slope[i];
    topo_i.aspect = aspect[i];
    p_topo_vec[i] = std::make_unique<Topography>(topo_i);
  } 
  WBCommunicationStructures WBcomm = WBCommunicationStructures(numCohorts_max, nlayers_max, ncanlayers_max, ntimesteps_max);
  p_WBcomm = std::make_unique<WBCommunicationStructures>(WBcomm);
  GROWTHCommunicationStructures GROWTHcomm = GROWTHCommunicationStructures(numCohorts_max, nlayers_max, ncanlayers_max, ntimesteps_max);
  p_GROWTHcomm = std::make_unique<GROWTHCommunicationStructures>(GROWTHcomm);
}
multiple_runner::~multiple_runner() {
}

void multiple_runner::run_day(Rcpp::CharacterVector date, Rcpp::List meteovec_list, bool parallelize) {

  if(meteovec_list.size() != n) throw medfate::MedfateInternalError("Wrong weather list size");
  
  std::vector<WeatherInputVector> weather_vec(n);
  for(int i=0; i<n;i++) {
    NumericVector meteovec_i = meteovec_list[i];
    weather_vec[i] = WeatherInputVector(meteovec_i);
  }
  std::string date_str = Rcpp::as<std::string>(date[0]);
  if(parallelize) {
    //build worker
    DAY_worker worker(*p_WBcomm,
                      *p_GROWTHcomm,
                       date_str,
                       p_x_vec,
                       latitude_vec,
                       p_topo_vec,
                       weather_vec,
                       p_result_vec);
    // call it with parallelFor
    parallelFor(0, n, worker);
  } else {
    double runon = 0.0;
    double waterTableDepth = medfate::NA_DOUBLE;
    for(int i=0;i<n;i++) {
      if((p_x_vec[i]->getInputClass() == "spwbInput") || (p_x_vec[i]->getInputClass() == "aspwbInput")) {
        WaterBalanceModelInput& x_i = dynamic_cast<WaterBalanceModelInput&>(*p_x_vec[i]);
        WB_RESULT& WBres_i = dynamic_cast<WB_RESULT&>(*p_result_vec[i]);
        int nlayers_i = x_i.soil.getNlayers();
        std::vector<double> lateralFlows_c(nlayers_i, 0.0);
        wb_day_inner_c(WBres_i, *p_WBcomm, x_i,
                       Rcpp::as<std::string>(date[0]),
                       WeatherInputVector(weather_vec[i]),
                       latitude_vec[i], p_topo_vec[i]->elevation, p_topo_vec[i]->slope, p_topo_vec[i]->aspect,
                       runon,
                       lateralFlows_c, waterTableDepth);
      } else if(p_x_vec[i]->getInputClass() == "growthInput") {
        ModelInput& x_i = dynamic_cast<ModelInput&>(*p_x_vec[i]);
        GROWTH_RESULT& GROWTHres_i = dynamic_cast<GROWTH_RESULT&>(*p_result_vec[i]);
        int nlayers_i = x_i.soil.getNlayers();
        std::vector<double> lateralFlows_c(nlayers_i, 0.0);
        growthDay_inner_c(GROWTHres_i, *p_GROWTHcomm, x_i,
                          Rcpp::as<std::string>(date[0]),
                          WeatherInputVector(weather_vec[i]),
                          latitude_vec[i], p_topo_vec[i]->elevation, p_topo_vec[i]->slope, p_topo_vec[i]->aspect,
                          runon,
                          lateralFlows_c, waterTableDepth);
      } else {
        throw medfate::MedfateInternalError("Wrong input class for simulation launching");
      }
    }
  }
}

//Returs output (decreases index)
Rcpp::List multiple_runner::get_output_at(int i) {
  List l;
  if(p_x_vec[i-1]->getInputClass() == "spwbInput") {
    ModelInput& x_i = dynamic_cast<ModelInput&>(*p_x_vec[i-1]);
    WB_RESULT& WBres = dynamic_cast<WB_RESULT&>(*p_result_vec[i-1]);
    l = copyWBResult_c(WBres, x_i);
  } else if(p_x_vec[i-1]->getInputClass() == "growthInput") {
    ModelInput& x_i = dynamic_cast<ModelInput&>(*p_x_vec[i-1]);
    GROWTH_RESULT& GROWTHres = dynamic_cast<GROWTH_RESULT&>(*p_result_vec[i-1]);
    l = copyGROWTHResult_c(GROWTHres, x_i);
  } else if(p_x_vec[i-1]->getInputClass() == "aspwbInput") {
    AgricultureModelInput& x_i = dynamic_cast<AgricultureModelInput&>(*p_x_vec[i-1]);
    WB_RESULT& WBres = dynamic_cast<WB_RESULT&>(*p_result_vec[i-1]);
    l = copyWBResult_c(WBres, x_i);
  } else {
    throw medfate::MedfateInternalError("Wrong input class for result copying");
  }
  return(l);
}

void multiple_runner::update_input_at(int i, List x_list) {
  if(p_x_vec[i-1]->getInputClass()=="spwbInput" || p_x_vec[i-1]->getInputClass()=="growthInput") {
    ModelInput& x_i = dynamic_cast<ModelInput&>(*p_x_vec[i-1]);
    x_i.copyStateToList(x_list);
  } else if(p_x_vec[i-1]->getInputClass()=="aspwbInput") {
    AgricultureModelInput& x_i = dynamic_cast<AgricultureModelInput&>(*p_x_vec[i-1]);
    x_i.copyStateToList(x_list);
  }
}


/* 
 * CREATE RCPP MODULES 
 */

RCPP_MODULE(mod_single) {
  class_<single_runner>( "single_runner" )
  .constructor<Rcpp::List, double, double, double, double>()
  .method( "run_day", &single_runner::run_day )
  .method( "get_output", &single_runner::get_output)
  .method( "update_input", &single_runner::update_input)
  ;
}
RCPP_MODULE(mod_multiple) {
  class_<multiple_runner>( "multiple_runner" )
  .constructor<Rcpp::List, NumericVector, NumericVector, NumericVector, NumericVector>()
  .method( "run_day", &multiple_runner::run_day )
  .method( "get_output_at", &multiple_runner::get_output_at)
  .method( "update_input_at", &multiple_runner::update_input_at)
  ;
}
