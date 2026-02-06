#include "RcppArmadillo.h"
#include "medfate.h"
#include "medfate_workers.h"
#include "modelInput_c.h"
#include "spwb_day_c.h"
using namespace Rcpp;

#ifndef MEDFATE_RUNNERS_H
#define MEDFATE_RUNNERS_H

class SPWB_runner {
private:
  ModelInput x;
  double latitude, elevation, slope, aspect;
  SPWBCommunicationStructures SPWBcomm;
public:
  SPWB_runner(List x_list, 
             double latitude, double elevation, double slope, double aspect);
  void run_day(CharacterVector date, NumericVector meteovec, 
               double runon = 0, Nullable<NumericVector> lateralFlows = R_NilValue, double waterTableDepth = NA_REAL);
  List get_output();
};

class SPWB_multiple_runner {
private:
  std::vector<double> latitude_vec;
  std::vector<Topography> topo_vec;
  std::vector<ModelInput> x_vec;
  std::vector<SPWB_RESULT> SPWB_res_vec;
public:
  SPWB_multiple_runner(List x_vec, 
                       NumericVector latitude_vec, NumericVector elevation_vec, NumericVector slope_vec, NumericVector aspect_vec);
  void run_day(CharacterVector date, List meteovec_list);
  List get_output_at(int i);
};

#endif