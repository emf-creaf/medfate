#include "RcppArmadillo.h"
#include "medfate.h"
#include "medfate_workers.h"
#include "modelInput_c.h"
#include "spwb_day_c.h"
#include "growth_day_c.h"
using namespace Rcpp;

#ifndef MEDFATE_RUNNERS_H
#define MEDFATE_RUNNERS_H

class single_runner {
private:
  std::unique_ptr<AbstractModelInput> p_x;
  double latitude;
  double elevation, slope, aspect;
  std::unique_ptr<WBCommunicationStructures> p_WBcomm;
  std::unique_ptr<GROWTHCommunicationStructures> p_GROWTHcomm;
  std::unique_ptr<ABSTRACTMODEL_RESULT> p_result;
public:
  single_runner(List x_list, 
                double latitude, double elevation, double slope, double aspect);
  void run_day(CharacterVector date, NumericVector meteovec, 
               double runon = 0, Nullable<NumericVector> lateralFlows = R_NilValue, double waterTableDepth = NA_REAL);
  ~single_runner();
  List get_output();
  void update_input(List x_list);
};

class multiple_runner {
private:
  int n;
  size_t numCohorts_max;
  size_t nlayers_max;
  size_t ncanlayers_max;
  size_t ntimesteps_max;
  std::vector<double> latitude_vec;
  std::vector<std::unique_ptr<Topography>> p_topo_vec;
  std::vector<std::unique_ptr<AbstractModelInput>> p_x_vec;
  std::vector<std::unique_ptr<ABSTRACTMODEL_RESULT>> p_result_vec;
public:
  multiple_runner(List x_vec, 
                  NumericVector latitude_vec, NumericVector elevation_vec, NumericVector slope_vec, NumericVector aspect_vec);
  ~multiple_runner();
  void run_day(CharacterVector date, List meteovec_list, bool parallelize = false);
  List get_output_at(int i);
  void store_output_at(int i, List l_vec);
  void update_input_at(int i, List x_list);
};

class watershed_runner {
private:
  int n;
  std::vector<double> latitude_vec;
  std::vector<std::unique_ptr<Topography>> p_topo_vec;
  std::vector<std::unique_ptr<AbstractModelInput>> p_x_vec;
  std::vector<std::unique_ptr<ABSTRACTMODEL_RESULT>> p_result_vec;
  std::unique_ptr<WBCommunicationStructures> p_WBcomm;
  std::unique_ptr<GROWTHCommunicationStructures> p_GROWTHcomm;
  List sf_routing;

public:
  watershed_runner(List x_vec, 
                   NumericVector latitude_vec, NumericVector elevation_vec, NumericVector slope_vec, NumericVector aspect_vec,
                   NumericVector snowpack_vec,
                   List sf_routing,
                   double KsatMultiplier = 1.0);
  ~watershed_runner();
  void run_day(CharacterVector date, DataFrame gridMeteo, 
               NumericVector waterTableDepth, 
               List output,
               double rock_max_infiltration,
               bool free_drainage_outlets,
               bool standSummary, bool fireHazardSummary, bool carbonBalanceSummary, bool biomassBalanceSummary);
  List get_output_at(int i);
  void store_output_at(int i, List l_vec);
  void update_input_at(int i, List x_list);
};
#endif
