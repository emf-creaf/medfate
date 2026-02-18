#include <RcppArmadillo.h>
#include "medfate_runners.h"  
#include "spwbland_const.h"
#include <string>
#include "spwb_day_c.h"
#include "growth_day_c.h"
#include "lowlevel_structures_c.h"
#include "numerical_solving_c.h"
using namespace Rcpp;

////////////////////////////////////
//    SINGLE-SITE RUNNER
////////////////////////////////////
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
  if (input_classIn == "aspwbInput") {
    AgricultureModelInput x = AgricultureModelInput(x_list);
    p_x = std::make_unique<AgricultureModelInput>(x);
    int nlayers = x.soil.getNlayers();
    AgricultureWB_RESULT AgrWBres(nlayers);
    p_result = std::make_unique<AgricultureWB_RESULT>(AgrWBres);
    WBCommunicationStructures WBcomm = WBCommunicationStructures(0, nlayers,0, 0);
    p_WBcomm = std::make_unique<WBCommunicationStructures>(WBcomm);
  } else if(input_classIn=="spwbInput") {
    ModelInput x = ModelInput(x_list);
    p_x = std::make_unique<ModelInput>(x);
    int numCohorts = x.cohorts.SpeciesIndex.size();
    int nlayers = x.soil.getNlayers();
    int ncanlayers = x.canopy.zlow.size();
    int ntimesteps = x.control.advancedWB.ndailysteps;
    if(x.control.transpirationMode=="Granier") {
      BasicTranspiration_RESULT BTres(numCohorts, nlayers);
      BasicSPWB_RESULT BSPWBres(BTres, nlayers);
      p_result = std::make_unique<BasicSPWB_RESULT>(BSPWBres);
    } else {
      AdvancedTranspiration_RESULT ATres(numCohorts, nlayers, ncanlayers, ntimesteps);
      AdvancedSPWB_RESULT ASPWBres(ATres, nlayers);
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
      BasicSPWB_RESULT BSPWBres(BTres, nlayers);
      BasicGROWTH_RESULT BGROWTHres(BSPWBres, numCohorts);
      p_result = std::make_unique<BasicGROWTH_RESULT>(BGROWTHres);
    } else if(x.control.transpirationMode=="Sperry" || x.control.transpirationMode=="Sureau")  {
      AdvancedTranspiration_RESULT ATres(numCohorts, nlayers, ncanlayers, ntimesteps);
      AdvancedSPWB_RESULT ASPWBres(ATres, nlayers);
      AdvancedGROWTH_RESULT AGROWTHres(ASPWBres, numCohorts, ntimesteps);
      p_result = std::make_unique<AdvancedGROWTH_RESULT>(AGROWTHres);
    }
    GROWTHCommunicationStructures GROWTHcomm = GROWTHCommunicationStructures(numCohorts, nlayers, ncanlayers, ntimesteps);
    p_GROWTHcomm = std::make_unique<GROWTHCommunicationStructures>(GROWTHcomm);
  } else {
    throw medfate::MedfateInternalError("Wrong model input class (should be spwbInput or aspwbInput)");
  }
}
single_runner::~single_runner(){}
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
    ModelInput& x = dynamic_cast<ModelInput&>(*p_x);
    WB_RESULT& WBres = dynamic_cast<WB_RESULT&>(*p_result);
    l = copyWBResult_c(WBres, x);
  } else if(p_x->getInputClass() == "growthInput") {
    ModelInput& x = dynamic_cast<ModelInput&>(*p_x);
    GROWTH_RESULT& GROWTHres = dynamic_cast<GROWTH_RESULT&>(*p_result);
    l = copyGROWTHResult_c(GROWTHres, x);
  } else if(p_x->getInputClass() == "aspwbInput") {
    AgricultureModelInput& x = dynamic_cast<AgricultureModelInput&>(*p_x);
    AgricultureWB_RESULT& WBres = dynamic_cast<AgricultureWB_RESULT&>(*p_result);
    l = copyWBResult_c(WBres, x);
  } else {
    throw medfate::MedfateInternalError("Wrong input class for result copying");
  }
  return(l);
}
// Updates state of (external) input object
void single_runner::update_input(List x_list) {
  if(p_x->getInputClass()=="spwbInput" || p_x->getInputClass()=="growthInput") {
    ModelInput& x = dynamic_cast<ModelInput&>(*p_x);
    x.copyStateToList(x_list);
  } else if(p_x->getInputClass()=="aspwbInput") {
    AgricultureModelInput& x = dynamic_cast<AgricultureModelInput&>(*p_x);
    x.copyStateToList(x_list);
  } else {
    throw medfate::MedfateInternalError("Wrong input class for input update");
  }
}

////////////////////////////////////
//    MULTIPLE-SITE RUNNER
////////////////////////////////////
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
        BasicSPWB_RESULT BSPWBres(BTres, nlayers_i);
        p_result_vec[i] = std::make_unique<BasicSPWB_RESULT>(BSPWBres);
      } else {
        AdvancedTranspiration_RESULT ATres(numCohorts_i, nlayers_i, ncanlayers_i, ntimesteps_i);
        AdvancedSPWB_RESULT ASPWBres(ATres, nlayers_i);
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
        BasicSPWB_RESULT BSPWBres(BTres, nlayers_i);
        BasicGROWTH_RESULT BGROWTHres(BSPWBres, numCohorts_i);
        p_result_vec[i] = std::make_unique<BasicGROWTH_RESULT>(BGROWTHres);
      } else {
        AdvancedTranspiration_RESULT ATres(numCohorts_i, nlayers_i, ncanlayers_i, ntimesteps_i);
        AdvancedSPWB_RESULT ASPWBres(ATres, nlayers_i);
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
//  Destructor for multiple_runner
multiple_runner::~multiple_runner() {}
//  Run one day for multiple cells
void multiple_runner::run_day(Rcpp::CharacterVector date, Rcpp::List meteovec_list, bool parallelize) {

  if(meteovec_list.size() != n) throw medfate::MedfateInternalError("Wrong weather list size");
  
  std::vector<std::unique_ptr<WeatherInputVector>> p_weather_vec(n);
  for(int i=0; i<n;i++) {
    NumericVector meteovec_i = meteovec_list[i];
    p_weather_vec[i] = std::make_unique<WeatherInputVector>(WeatherInputVector(meteovec_i));
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
                       p_weather_vec,
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
                       *p_weather_vec[i],
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
                          *p_weather_vec[i],
                          latitude_vec[i], p_topo_vec[i]->elevation, p_topo_vec[i]->slope, p_topo_vec[i]->aspect,
                          runon,
                          lateralFlows_c, waterTableDepth);
      } else {
        throw medfate::MedfateInternalError("Wrong input class for simulation launching");
      }
    }
  }
}
//Returns output (decreases index)
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
//Updates state of (external) input object
void multiple_runner::update_input_at(int i, List x_list) {
  if(p_x_vec[i-1]->getInputClass()=="spwbInput" || p_x_vec[i-1]->getInputClass()=="growthInput") {
    ModelInput& x_i = dynamic_cast<ModelInput&>(*p_x_vec[i-1]);
    x_i.copyStateToList(x_list);
  } else if(p_x_vec[i-1]->getInputClass()=="aspwbInput") {
    AgricultureModelInput& x_i = dynamic_cast<AgricultureModelInput&>(*p_x_vec[i-1]);
    x_i.copyStateToList(x_list);
  } else {
    throw medfate::MedfateInternalError("Wrong input class for input update");
  }
}

////////////////////////////////////
//    WATERSHED RUNNER
////////////////////////////////////
//  Constructor for watershed_runner
watershed_runner::watershed_runner(List input_vec, 
                                   NumericVector latitude, 
                                   NumericVector elevation, 
                                   NumericVector slope, 
                                   NumericVector aspect,
                                   List sf_routingIn) :
  n(input_vec.size()), latitude_vec(n), p_topo_vec(n), p_x_vec(n), p_result_vec(n), sf_routing(sf_routingIn) {
  size_t numCohorts_max = 0;
  size_t nlayers_max = 0;
  size_t ncanlayers_max = 0;
  size_t ntimesteps_max = 0;
  for(int i=0; i< n;i++) {
    Rcpp::List input_i = input_vec[i];
    Rcpp::CharacterVector classVector_i = input_i.attr("class");
    Rcpp::String s_i = classVector_i[0];
    std::string input_class_i = s_i.get_cstring();
    if(input_class_i=="nswbInput") {
      NonSoilWaterBalanceModelInput x_i = NonSoilWaterBalanceModelInput(input_i);
      p_x_vec[i] = std::make_unique<NonSoilWaterBalanceModelInput>(x_i);
      NSWB_RESULT NSWBres;
      p_result_vec[i] = std::make_unique<NSWB_RESULT>(NSWBres);
    } else if(input_class_i=="spwbInput") {
      ModelInput x_i = ModelInput(input_i);
      p_x_vec[i] = std::make_unique<ModelInput>(x_i);
      int numCohorts_i = x_i.cohorts.SpeciesIndex.size();
      int nlayers_i = x_i.soil.getNlayers();
      int ncanlayers_i = x_i.canopy.zlow.size();
      int ntimesteps_i = x_i.control.advancedWB.ndailysteps;
      if(x_i.control.transpirationMode=="Granier") {
        BasicTranspiration_RESULT BTres(numCohorts_i, nlayers_i);
        BasicSPWB_RESULT BSPWBres(BTres, nlayers_i);
        p_result_vec[i] = std::make_unique<BasicSPWB_RESULT>(BSPWBres);
      } else {
        AdvancedTranspiration_RESULT ATres(numCohorts_i, nlayers_i, ncanlayers_i, ntimesteps_i);
        AdvancedSPWB_RESULT ASPWBres(ATres, nlayers_i);
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
        BasicSPWB_RESULT BSPWBres(BTres, nlayers_i);
        BasicGROWTH_RESULT BGROWTHres(BSPWBres, numCohorts_i);
        p_result_vec[i] = std::make_unique<BasicGROWTH_RESULT>(BGROWTHres);
      } else {
        AdvancedTranspiration_RESULT ATres(numCohorts_i, nlayers_i, ncanlayers_i, ntimesteps_i);
        AdvancedSPWB_RESULT ASPWBres(ATres, nlayers_i);
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
//  Destructor for watershed_runner
watershed_runner::~watershed_runner() {}
//  Run one day for watershed
void watershed_runner::run_day(Rcpp::CharacterVector date, DataFrame gridMeteo, 
                               NumericVector waterTableDepth, 
                               List output,
                               double rock_max_infiltration,
                               bool free_drainage_outlets,
                               bool standSummary, bool fireHazardSummary, bool carbonBalanceSummary, bool biomassBalanceSummary) {
  
  if(waterTableDepth.size() != n) throw medfate::MedfateInternalError("Wrong waterTableDepth size");
  if(gridMeteo.nrow() != n) throw medfate::MedfateInternalError("Wrong grid input size");
  
  DataFrame outWB = Rcpp::as<Rcpp::DataFrame>(output["WatershedWaterBalance"]);
  List localResults = output["LocalResults"];
  
  //WATERSHED-LEVEL OUTPUT
  NumericVector Runoff=  outWB[WBCOM_Runoff];
  NumericVector Runon=  outWB[WBCOM_Runon];
  NumericVector MinTemperature = outWB[WBCOM_MinTemperature];
  NumericVector MaxTemperature = outWB[WBCOM_MaxTemperature];
  NumericVector PET = outWB[WBCOM_PET];
  NumericVector Rain = outWB[WBCOM_Rain];
  NumericVector Interception = outWB[WBCOM_Interception];
  NumericVector NetRain = outWB[WBCOM_NetRain];
  NumericVector Infiltration = outWB[WBCOM_Infiltration];
  NumericVector InfiltrationExcess = outWB[WBCOM_InfiltrationExcess];
  NumericVector InterflowBalance = outWB[WBCOM_InterflowBalance];
  NumericVector DeepDrainage = outWB[WBCOM_DeepDrainage];
  NumericVector Snow = outWB[WBCOM_Snow];
  NumericVector Snowmelt = outWB[WBCOM_Snowmelt];
  NumericVector ChannelExport = outWB[WBCOM_ChannelExport];
  NumericVector WatershedExport = outWB[WBCOM_WatershedExport];
  NumericVector CapillarityRise = outWB[WBCOM_CapillarityRise];
  NumericVector SoilEvaporation = outWB[WBCOM_SoilEvaporation];
  NumericVector SaturationExcess = outWB[WBCOM_SaturationExcess];
  NumericVector Transpiration = outWB[WBCOM_Transpiration];
  NumericVector HerbTranspiration = outWB[WBCOM_HerbTranspiration];
  NumericVector LAI, LAIherb, LAIlive, LAIexpanded, LAIdead, Cm, LgroundPAR, LgroundSWR;
  NumericVector StructuralBalance, LabileBalance, PlantBalance, MortalityLoss, CohortBalance;
  NumericVector GrossPrimaryProduction, MaintenanceRespiration, SynthesisRespiration, NetPrimaryProduction;
  NumericVector Loading_understory, Loading_overstory, CFMC_understory, CFMC_overstory, DFMC, ROS_surface, I_b_surface, t_r_surface, FL_surface, Ic_ratio, ROS_crown, I_b_crown, t_r_crown, FL_crown, SFP, CFP;
  if(standSummary) {
    DataFrame outStand = as<DataFrame>(output["WatershedStand"]);
    LAI = outStand[STCOM_LAI];
    LAIherb = outStand[STCOM_LAIherb];
    LAIlive = outStand[STCOM_LAIlive];
    LAIexpanded = outStand[STCOM_LAIexpanded];
    LAIdead = outStand[STCOM_LAIdead];
    Cm = outStand[STCOM_Cm];
    LgroundPAR = outStand[STCOM_LgroundPAR];
    LgroundSWR = outStand[STCOM_LgroundSWR];
  }
  if(fireHazardSummary) {
    DataFrame fireStand = as<DataFrame>(output["WatershedFireHazard"]);
    Loading_overstory = fireStand[FHCOM_Loading_overstory];
    Loading_understory = fireStand[FHCOM_Loading_understory];
    CFMC_overstory = fireStand[FHCOM_CFMC_overstory];
    CFMC_understory = fireStand[FHCOM_CFMC_understory];
    DFMC = fireStand[FHCOM_DFMC];
    ROS_surface = fireStand[FHCOM_ROS_surface];
    I_b_surface = fireStand[FHCOM_I_b_surface];
    t_r_surface = fireStand[FHCOM_t_r_surface];
    FL_surface = fireStand[FHCOM_FL_surface];
    Ic_ratio = fireStand[FHCOM_Ic_ratio];
    ROS_crown = fireStand[FHCOM_ROS_crown];
    I_b_crown = fireStand[FHCOM_I_b_crown];
    t_r_crown = fireStand[FHCOM_t_r_crown];
    FL_crown = fireStand[FHCOM_FL_crown];
    SFP = fireStand[FHCOM_SFP];
    CFP = fireStand[FHCOM_CFP];
  }
  if(carbonBalanceSummary) {
    DataFrame outCB = as<DataFrame>(output["WatershedCarbonBalance"]);
    GrossPrimaryProduction = outCB[CBCOM_GrossPrimaryProduction];
    MaintenanceRespiration = outCB[CBCOM_MaintenanceRespiration];
    SynthesisRespiration = outCB[CBCOM_SynthesisRespiration];
    NetPrimaryProduction = outCB[CBCOM_NetPrimaryProduction];
  }
  if(biomassBalanceSummary) {
    DataFrame outBB = as<DataFrame>(output["WatershedBiomassBalance"]);
    StructuralBalance = outBB[BBCOM_StructuralBalance];
    LabileBalance = outBB[BBCOM_LabileBalance];
    PlantBalance = outBB[BBCOM_PlantBalance];
    MortalityLoss = outBB[BBCOM_MortalityLoss];
    CohortBalance = outBB[BBCOM_CohortBalance];
  }
  
  // WEATHER INPUT VECTORS
  NumericVector tminVec = Rcpp::as<Rcpp::NumericVector>(gridMeteo["MinTemperature"]);
  NumericVector tmaxVec = Rcpp::as<Rcpp::NumericVector>(gridMeteo["MaxTemperature"]);
  NumericVector rhminVec = Rcpp::as<Rcpp::NumericVector>(gridMeteo["MinRelativeHumidity"]);
  NumericVector rhmaxVec = Rcpp::as<Rcpp::NumericVector>(gridMeteo["MaxRelativeHumidity"]);
  NumericVector precVec = Rcpp::as<Rcpp::NumericVector>(gridMeteo["Precipitation"]);
  NumericVector radVec = Rcpp::as<Rcpp::NumericVector>(gridMeteo["Radiation"]);
  NumericVector wsVec = Rcpp::as<Rcpp::NumericVector>(gridMeteo["WindSpeed"]);
  NumericVector C02Vec = Rcpp::as<Rcpp::NumericVector>(gridMeteo["CO2"]);
  
  // ROUTING INFORMATION
  IntegerVector waterOrder = sf_routing["waterOrder"];
  List queenNeigh = sf_routing["queenNeigh"];
  List waterQ = sf_routing["waterQ"];
  LogicalVector isChannel = sf_routing["channel"];
  LogicalVector isOutlet = sf_routing["outlet"];

  // WEATHER VECTOR TO BE RECYCLED
  WeatherInputVector meteovec;
  
  for(int i=0;i<n;i++) {
    //get next cell in order
    int iCell = waterOrder[i]-1; //Decrease index!!!!
    const std::string& input_class_iCell = p_x_vec[iCell]->getInputClass();
    if((input_class_iCell == "spwbInput") || (input_class_iCell == "growthInput") || (input_class_iCell == "aspwbInput")) {
      
      WaterBalanceModelInput& x_iCell = dynamic_cast<WaterBalanceModelInput&>(*p_x_vec[iCell]);
      
      //Soil cell: Prepare input
      meteovec.tmin = tminVec[iCell];
      meteovec.tmax = tmaxVec[iCell];
      meteovec.rhmin = rhminVec[iCell];
      meteovec.rhmax = rhmaxVec[iCell];
      meteovec.prec = precVec[iCell];
      meteovec.rad = radVec[iCell];
      meteovec.wind = wsVec[iCell];
      meteovec.Catm = C02Vec[iCell];
      
      
      
      double wtd_c = waterTableDepth[iCell];
      // This effectively uncouples capillarity rise from aquifer elevation but helps avoiding instabilities 
      // in the local balance of outlet cells 
      if(free_drainage_outlets && (isOutlet[iCell] || isChannel[iCell])) wtd_c = medfate::NA_DOUBLE;

      //Assume water will exit/enter faster in soil layers that have more water
      //This is better than assuming weights equal to layer widths, which can result
      //In water leaving bottom layers when rain fills the top layers
      int nlayers = x_iCell.soil.getNlayers();
      std::vector<double> wl(nlayers);
      std::vector<double> lateralFlows_c(nlayers);
      double wsum = 0.0;
      for(int l=0;l<nlayers;l++) {
        wl[l] = x_iCell.soil.getWater(l);
        wsum += wl[l];
      }
      for(int l=0;l<nlayers;l++) {
        lateralFlows_c[l] = InterflowBalance[iCell]*(wl[l]/wsum);
      }
      if((input_class_iCell == "spwbInput") || (input_class_iCell == "aspwbInput")) {
        WaterBalanceModelInput& x_iCell = dynamic_cast<WaterBalanceModelInput&>(*p_x_vec[iCell]);
        WB_RESULT& WBres_iCell = dynamic_cast<WB_RESULT&>(*p_result_vec[iCell]);
        wb_day_inner_c(WBres_iCell, *p_WBcomm, x_iCell,
                       Rcpp::as<std::string>(date[0]),
                       meteovec,
                       latitude_vec[iCell], p_topo_vec[iCell]->elevation, p_topo_vec[iCell]->slope, p_topo_vec[iCell]->aspect,
                       Runon[iCell],
                       lateralFlows_c, wtd_c);
        if(input_class_iCell == "spwbInput") {
          SPWB_RESULT& WBres_iCell = dynamic_cast<SPWB_RESULT&>(*p_result_vec[iCell]);
          Snow[iCell] = WBres_iCell.WaterBalance.Snow;
          Snowmelt[iCell] = WBres_iCell.WaterBalance.Snowmelt;
          PET[iCell] = WBres_iCell.WaterBalance.PET;
          Rain[iCell] = WBres_iCell.WaterBalance.Rain;
          SoilEvaporation[iCell] = WBres_iCell.WaterBalance.SoilEvaporation;
          NetRain[iCell] = WBres_iCell.WaterBalance.NetRain;
          Infiltration[iCell] = WBres_iCell.WaterBalance.Infiltration;
          Runoff[iCell] = WBres_iCell.WaterBalance.Runoff;
          InfiltrationExcess[iCell] = WBres_iCell.WaterBalance.InfiltrationExcess;
          SaturationExcess[iCell] = WBres_iCell.WaterBalance.SaturationExcess;
          DeepDrainage[iCell] = WBres_iCell.WaterBalance.DeepDrainage;
          CapillarityRise[iCell] = WBres_iCell.WaterBalance.CapillarityRise;
          Transpiration[iCell] = WBres_iCell.WaterBalance.Transpiration;
          HerbTranspiration[iCell] = WBres_iCell.WaterBalance.HerbTranspiration; 
          if(standSummary) {
            LAI[iCell] = WBres_iCell.Stand.LAI;
            LAIherb[iCell] = WBres_iCell.Stand.LAIherb;
            LAIlive[iCell] = WBres_iCell.Stand.LAIlive;
            LAIexpanded[iCell] = WBres_iCell.Stand.LAIexpanded;
            LAIdead[iCell] = WBres_iCell.Stand.LAIdead;
            Cm[iCell] = WBres_iCell.Stand.Cm;
            LgroundPAR[iCell] = WBres_iCell.Stand.LgroundPAR;
            LgroundSWR[iCell] = WBres_iCell.Stand.LgroundSWR;
          }
          if(fireHazardSummary) {
            Loading_overstory[iCell] = WBres_iCell.fccs.Loading_overstory;
            Loading_understory[iCell] = WBres_iCell.fccs.Loading_understory;
            CFMC_understory[iCell] = WBres_iCell.fccs.CFMC_understory;
            CFMC_overstory[iCell] = WBres_iCell.fccs.CFMC_overstory;
            DFMC[iCell] = WBres_iCell.fccs.DFMC;
            ROS_surface[iCell] = WBres_iCell.fccs.ROS_surface;
            I_b_surface[iCell] = WBres_iCell.fccs.I_b_surface;
            t_r_surface[iCell] = WBres_iCell.fccs.t_r_surface;
            FL_surface[iCell] = WBres_iCell.fccs.FL_surface;
            Ic_ratio[iCell] = WBres_iCell.fccs.Ic_ratio;
            ROS_crown[iCell] = WBres_iCell.fccs.ROS_crown;
            I_b_crown[iCell] = WBres_iCell.fccs.I_b_crown;
            t_r_crown[iCell] = WBres_iCell.fccs.t_r_crown;
            FL_crown[iCell] = WBres_iCell.fccs.FL_crown;
            SFP[iCell] = WBres_iCell.fccs.SFP;
            CFP[iCell] = WBres_iCell.fccs.CFP;
          }
        } else if(input_class_iCell == "aspwbInput") {
          AgricultureWB_RESULT& AgrWBres_iCell = dynamic_cast<AgricultureWB_RESULT&>(*p_result_vec[iCell]);
          Snow[iCell] = AgrWBres_iCell.WaterBalance.Snow;
          Snowmelt[iCell] = AgrWBres_iCell.WaterBalance.Snowmelt;
          PET[iCell] = AgrWBres_iCell.WaterBalance.PET;
          Rain[iCell] = AgrWBres_iCell.WaterBalance.Rain;
          SoilEvaporation[iCell] = AgrWBres_iCell.WaterBalance.SoilEvaporation;
          NetRain[iCell] = AgrWBres_iCell.WaterBalance.NetRain;
          Infiltration[iCell] = AgrWBres_iCell.WaterBalance.Infiltration;
          Runoff[iCell] = AgrWBres_iCell.WaterBalance.Runoff;
          InfiltrationExcess[iCell] = AgrWBres_iCell.WaterBalance.InfiltrationExcess;
          SaturationExcess[iCell] = AgrWBres_iCell.WaterBalance.SaturationExcess;
          DeepDrainage[iCell] = AgrWBres_iCell.WaterBalance.DeepDrainage;
          CapillarityRise[iCell] = AgrWBres_iCell.WaterBalance.CapillarityRise;   
        }
      } else if(input_class_iCell == "growthInput") {
        ModelInput& x_iCell = dynamic_cast<ModelInput&>(*p_x_vec[iCell]);
        GROWTH_RESULT& GROWTHres_iCell = dynamic_cast<GROWTH_RESULT&>(*p_result_vec[iCell]);
        growthDay_inner_c(GROWTHres_iCell, *p_GROWTHcomm, x_iCell,
                          Rcpp::as<std::string>(date[0]),
                          meteovec,
                          latitude_vec[iCell], p_topo_vec[iCell]->elevation, p_topo_vec[iCell]->slope, p_topo_vec[iCell]->aspect,
                          Runon[iCell],
                          lateralFlows_c, wtd_c);
        if(x_iCell.control.transpirationMode=="Granier") {
          BasicGROWTH_RESULT& BGROWTHres_iCell = dynamic_cast<BasicGROWTH_RESULT&>(*p_result_vec[iCell]);
          Snow[iCell] = BGROWTHres_iCell.BSPWBres.WaterBalance.Snow;
          Snowmelt[iCell] = BGROWTHres_iCell.BSPWBres.WaterBalance.Snowmelt;
          PET[iCell] = BGROWTHres_iCell.BSPWBres.WaterBalance.PET;
          Rain[iCell] = BGROWTHres_iCell.BSPWBres.WaterBalance.Rain;
          SoilEvaporation[iCell] = BGROWTHres_iCell.BSPWBres.WaterBalance.SoilEvaporation;
          NetRain[iCell] = BGROWTHres_iCell.BSPWBres.WaterBalance.NetRain;
          Infiltration[iCell] = BGROWTHres_iCell.BSPWBres.WaterBalance.Infiltration;
          Runoff[iCell] = BGROWTHres_iCell.BSPWBres.WaterBalance.Runoff;
          InfiltrationExcess[iCell] = BGROWTHres_iCell.BSPWBres.WaterBalance.InfiltrationExcess;
          SaturationExcess[iCell] = BGROWTHres_iCell.BSPWBres.WaterBalance.SaturationExcess;
          DeepDrainage[iCell] = BGROWTHres_iCell.BSPWBres.WaterBalance.DeepDrainage;
          CapillarityRise[iCell] = BGROWTHres_iCell.BSPWBres.WaterBalance.CapillarityRise;   
          HerbTranspiration[iCell] = BGROWTHres_iCell.BSPWBres.WaterBalance.HerbTranspiration; 
          if(standSummary) {
            LAI[iCell] =  BGROWTHres_iCell.BSPWBres.Stand.LAI;
            LAIherb[iCell] =  BGROWTHres_iCell.BSPWBres.Stand.LAIherb;
            LAIlive[iCell] =  BGROWTHres_iCell.BSPWBres.Stand.LAIlive;
            LAIexpanded[iCell] =  BGROWTHres_iCell.BSPWBres.Stand.LAIexpanded;
            LAIdead[iCell] =  BGROWTHres_iCell.BSPWBres.Stand.LAIdead;
            Cm[iCell] =  BGROWTHres_iCell.BSPWBres.Stand.Cm;
            LgroundPAR[iCell] =  BGROWTHres_iCell.BSPWBres.Stand.LgroundPAR;
            LgroundSWR[iCell] =  BGROWTHres_iCell.BSPWBres.Stand.LgroundSWR;
          }
          if(fireHazardSummary) {
            Loading_overstory[iCell] = BGROWTHres_iCell.BSPWBres.fccs.Loading_overstory;
            Loading_understory[iCell] = BGROWTHres_iCell.BSPWBres.fccs.Loading_understory;
            CFMC_understory[iCell] = BGROWTHres_iCell.BSPWBres.fccs.CFMC_understory;
            CFMC_overstory[iCell] = BGROWTHres_iCell.BSPWBres.fccs.CFMC_overstory;
            DFMC[iCell] = BGROWTHres_iCell.BSPWBres.fccs.DFMC;
            ROS_surface[iCell] = BGROWTHres_iCell.BSPWBres.fccs.ROS_surface;
            I_b_surface[iCell] = BGROWTHres_iCell.BSPWBres.fccs.I_b_surface;
            t_r_surface[iCell] = BGROWTHres_iCell.BSPWBres.fccs.t_r_surface;
            FL_surface[iCell] = BGROWTHres_iCell.BSPWBres.fccs.FL_surface;
            Ic_ratio[iCell] = BGROWTHres_iCell.BSPWBres.fccs.Ic_ratio;
            ROS_crown[iCell] = BGROWTHres_iCell.BSPWBres.fccs.ROS_crown;
            I_b_crown[iCell] = BGROWTHres_iCell.BSPWBres.fccs.I_b_crown;
            t_r_crown[iCell] = BGROWTHres_iCell.BSPWBres.fccs.t_r_crown;
            FL_crown[iCell] = BGROWTHres_iCell.BSPWBres.fccs.FL_crown;
            SFP[iCell] = BGROWTHres_iCell.BSPWBres.fccs.SFP;
            CFP[iCell] = BGROWTHres_iCell.BSPWBres.fccs.CFP;
          }
        } else {
          AdvancedGROWTH_RESULT& AGROWTHres_iCell = dynamic_cast<AdvancedGROWTH_RESULT&>(*p_result_vec[iCell]);
          Snow[iCell] = AGROWTHres_iCell.ASPWBres.WaterBalance.Snow;
          Snowmelt[iCell] = AGROWTHres_iCell.ASPWBres.WaterBalance.Snowmelt;
          PET[iCell] = AGROWTHres_iCell.ASPWBres.WaterBalance.PET;
          Rain[iCell] = AGROWTHres_iCell.ASPWBres.WaterBalance.Rain;
          SoilEvaporation[iCell] = AGROWTHres_iCell.ASPWBres.WaterBalance.SoilEvaporation;
          NetRain[iCell] = AGROWTHres_iCell.ASPWBres.WaterBalance.NetRain;
          Infiltration[iCell] = AGROWTHres_iCell.ASPWBres.WaterBalance.Infiltration;
          Runoff[iCell] = AGROWTHres_iCell.ASPWBres.WaterBalance.Runoff;
          InfiltrationExcess[iCell] = AGROWTHres_iCell.ASPWBres.WaterBalance.InfiltrationExcess;
          SaturationExcess[iCell] = AGROWTHres_iCell.ASPWBres.WaterBalance.SaturationExcess;
          DeepDrainage[iCell] = AGROWTHres_iCell.ASPWBres.WaterBalance.DeepDrainage;
          CapillarityRise[iCell] = AGROWTHres_iCell.ASPWBres.WaterBalance.CapillarityRise;   
          HerbTranspiration[iCell] = AGROWTHres_iCell.ASPWBres.WaterBalance.HerbTranspiration;   
          if(standSummary) {
            LAI[iCell] =  AGROWTHres_iCell.ASPWBres.Stand.LAI;
            LAIherb[iCell] =  AGROWTHres_iCell.ASPWBres.Stand.LAIherb;
            LAIlive[iCell] =  AGROWTHres_iCell.ASPWBres.Stand.LAIlive;
            LAIexpanded[iCell] =  AGROWTHres_iCell.ASPWBres.Stand.LAIexpanded;
            LAIdead[iCell] = AGROWTHres_iCell.ASPWBres.Stand.LAIdead;
            Cm[iCell] =  AGROWTHres_iCell.ASPWBres.Stand.Cm;
            LgroundPAR[iCell] =  AGROWTHres_iCell.ASPWBres.Stand.LgroundPAR;
            LgroundSWR[iCell] =  AGROWTHres_iCell.ASPWBres.Stand.LgroundSWR;
          }
          if(fireHazardSummary) {
            Loading_overstory[iCell] = AGROWTHres_iCell.ASPWBres.fccs.Loading_overstory;
            Loading_understory[iCell] = AGROWTHres_iCell.ASPWBres.fccs.Loading_understory;
            CFMC_understory[iCell] = AGROWTHres_iCell.ASPWBres.fccs.CFMC_understory;
            CFMC_overstory[iCell] = AGROWTHres_iCell.ASPWBres.fccs.CFMC_overstory;
            DFMC[iCell] = AGROWTHres_iCell.ASPWBres.fccs.DFMC;
            ROS_surface[iCell] = AGROWTHres_iCell.ASPWBres.fccs.ROS_surface;
            I_b_surface[iCell] = AGROWTHres_iCell.ASPWBres.fccs.I_b_surface;
            t_r_surface[iCell] = AGROWTHres_iCell.ASPWBres.fccs.t_r_surface;
            FL_surface[iCell] = AGROWTHres_iCell.ASPWBres.fccs.FL_surface;
            Ic_ratio[iCell] = AGROWTHres_iCell.ASPWBres.fccs.Ic_ratio;
            ROS_crown[iCell] = AGROWTHres_iCell.ASPWBres.fccs.ROS_crown;
            I_b_crown[iCell] = AGROWTHres_iCell.ASPWBres.fccs.I_b_crown;
            t_r_crown[iCell] = AGROWTHres_iCell.ASPWBres.fccs.t_r_crown;
            FL_crown[iCell] = AGROWTHres_iCell.ASPWBres.fccs.FL_crown;
            SFP[iCell] = AGROWTHres_iCell.ASPWBres.fccs.SFP;
            CFP[iCell] = AGROWTHres_iCell.ASPWBres.fccs.CFP;
          }
        }
        if(carbonBalanceSummary){
          GrossPrimaryProduction[iCell] = GROWTHres_iCell.standCB.GrossPrimaryProduction;
          MaintenanceRespiration[iCell] = GROWTHres_iCell.standCB.MaintenanceRespiration;
          SynthesisRespiration[iCell] = GROWTHres_iCell.standCB.SynthesisRespiration;
          NetPrimaryProduction[iCell] = GROWTHres_iCell.standCB.NetPrimaryProduction;
        }
        if(biomassBalanceSummary) {
          StructuralBalance[iCell] = vecsum(GROWTHres_iCell.PBBres.StructuralBiomassBalance);
          LabileBalance[iCell] = vecsum(GROWTHres_iCell.PBBres.LabileBiomassBalance);
          PlantBalance[iCell] = vecsum(GROWTHres_iCell.PBBres.PlantBiomassBalance);
          MortalityLoss[iCell] = vecsum(GROWTHres_iCell.PBBres.MortalityBiomassLoss);
          CohortBalance[iCell] = vecsum(GROWTHres_iCell.PBBres.CohortBiomassBalance);
        }
      } else {
        throw medfate::MedfateInternalError("Wrong input class for simulation launching");
      }
    } else if(input_class_iCell == "nswbInput") {
      NonSoilWaterBalanceModelInput& x_iCell = dynamic_cast<NonSoilWaterBalanceModelInput&>(*p_x_vec[iCell]);
      NSWB_RESULT& NSWBres_iCell = dynamic_cast<NSWB_RESULT&>(*p_result_vec[iCell]);
      nswbDay_c(NSWBres_iCell, x_iCell,
                meteovec,
                p_topo_vec[iCell]->elevation, p_topo_vec[iCell]->slope, p_topo_vec[iCell]->aspect,
                Runon[iCell],
                rock_max_infiltration);
      PET[iCell] = medfate::NA_DOUBLE;
      Snow[iCell] = NSWBres_iCell.WaterBalance.Snow;
      Snowmelt[iCell] = NSWBres_iCell.WaterBalance.Snowmelt;
      Rain[iCell] =  NSWBres_iCell.WaterBalance.Rain;
      Infiltration[iCell] = NSWBres_iCell.WaterBalance.Infiltration;
      InfiltrationExcess[iCell] = NSWBres_iCell.WaterBalance.InfiltrationExcess;
      DeepDrainage[iCell] = NSWBres_iCell.WaterBalance.DeepDrainage;
      Runoff[iCell] = NSWBres_iCell.WaterBalance.Runoff;
      NetRain[iCell] = NSWBres_iCell.WaterBalance.NetRain;
    }
    //Common cell output
    MinTemperature[iCell] = tminVec[iCell];
    MaxTemperature[iCell] = tmaxVec[iCell];
    Interception[iCell] = Rain[iCell] - NetRain[iCell];
    
    // OVERLAND RUNOFF
    // Assign runoff to runon of downhill neighbours
    double ri_tot =  Runoff[iCell];
    NumericVector qi = Rcpp::as<Rcpp::NumericVector>(waterQ[iCell]);
    IntegerVector ni = Rcpp::as<Rcpp::IntegerVector>(queenNeigh[iCell]);
    
    // Add aquifer exfiltration to the water to be distributed if not an outlet or channel
    if(ri_tot>0.0) {
      double ri = ri_tot;
      if(isChannel[iCell]) { // If is channel then export
        ChannelExport[iCell] += ri;
        ri = 0.0;
      } else if(sum(qi)==0.0) {// If is outlet then export
        WatershedExport[iCell] += ri;
        ri = 0.0;
      } else { // Otherwise, distribute among waterQ neighbours
        if(ni.size()>0) {
          for(int j=0;j<ni.size();j++)  {
            Runon[ni[j]-1] += (qi[j]*ri_tot); //decrease index
            ri -= (qi[j]*ri_tot);
          }
        }
        if(ri > 0.000001) {
          Rcout<< i <<ni.size()<< " "<<qi.size()<<" "<<iCell<< " "<< sum(qi)<< " "<< ri<<"\n";
          stop("Non-outlet or channel cell with runoff export");
        }
      }
    }
  }
}


//Returns output (decreases index)
Rcpp::List watershed_runner::get_output_at(int i) {
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
  } else if(p_x_vec[i-1]->getInputClass() == "nswbInput") {
    NonSoilWaterBalanceModelInput& x_i = dynamic_cast<NonSoilWaterBalanceModelInput&>(*p_x_vec[i-1]);
    NSWB_RESULT& NSWBres = dynamic_cast<NSWB_RESULT&>(*p_result_vec[i-1]);
    l = copyNSWBResult_c(NSWBres, x_i);
  } else {
    throw medfate::MedfateInternalError("Wrong input class for result copying");
  }
  return(l);
}
//Updates state of (external) input object
void watershed_runner::update_input_at(int i, List x_list) {
  if(p_x_vec[i-1]->getInputClass()=="spwbInput" || p_x_vec[i-1]->getInputClass()=="growthInput") {
    ModelInput& x_i = dynamic_cast<ModelInput&>(*p_x_vec[i-1]);
    x_i.copyStateToList(x_list);
  } else if(p_x_vec[i-1]->getInputClass()=="aspwbInput") {
    AgricultureModelInput& x_i = dynamic_cast<AgricultureModelInput&>(*p_x_vec[i-1]);
    x_i.copyStateToList(x_list);
  } else if(p_x_vec[i-1]->getInputClass()=="nswbInput") {
    NonSoilWaterBalanceModelInput& x_i = dynamic_cast<NonSoilWaterBalanceModelInput&>(*p_x_vec[i-1]);
    x_i.copyStateToList(x_list);
  } else {
    throw medfate::MedfateInternalError("Wrong input class for input update");
  }
}


/////////////////////////
// CREATE RCPP MODULE 
/////////////////////////
RCPP_MODULE(runners) {
  class_<single_runner>( "single_runner" )
  .constructor<Rcpp::List, double, double, double, double>()
  .method( "run_day", &single_runner::run_day )
  .method( "get_output", &single_runner::get_output)
  .method( "update_input", &single_runner::update_input)
  ;
  class_<multiple_runner>( "multiple_runner" )
    .constructor<Rcpp::List, NumericVector, NumericVector, NumericVector, NumericVector>()
    .method( "run_day", &multiple_runner::run_day )
    .method( "get_output_at", &multiple_runner::get_output_at)
    .method( "update_input_at", &multiple_runner::update_input_at)
  ;
  class_<watershed_runner>( "watershed_runner" )
    .constructor<Rcpp::List, NumericVector, NumericVector, NumericVector, NumericVector, Rcpp::List>()
    .method( "run_day", &watershed_runner::run_day )
    .method( "get_output_at", &watershed_runner::get_output_at)
    .method( "update_input_at", &watershed_runner::update_input_at)
  ;
}
