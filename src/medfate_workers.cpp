// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
#include "medfate_workers.h"
using namespace RcppParallel;

// Implements simulation for the range of elements 
void DAY_worker::operator()(std::size_t begin, std::size_t end) {
  for(size_t i = begin; i< end; i++) {
    double runon = 0.0;
    double waterTableDepth = medfate::NA_DOUBLE;
    if((p_x_vec[i]->getInputClass() == "spwbInput") || (p_x_vec[i]->getInputClass() == "aspwbInput")) {
      WaterBalanceModelInput& x_i = dynamic_cast<WaterBalanceModelInput&>(*p_x_vec[i]);
      WB_RESULT& WBres_i = dynamic_cast<WB_RESULT&>(*p_result_vec[i]);
      int nlayers_i = x_i.soil.getNlayers();
      std::vector<double> lateralFlows_c(nlayers_i, 0.0);
      wb_day_inner_c(WBres_i, WBcomm, x_i,
                     date,
                     *p_weather_vec[i],
                     latitude_vec[i], p_topo_vec[i]->elevation, p_topo_vec[i]->slope, p_topo_vec[i]->aspect,
                     runon,
                     lateralFlows_c, waterTableDepth);
    } else if(p_x_vec[i]->getInputClass() == "growthInput") {
      ModelInput& x_i = dynamic_cast<ModelInput&>(*p_x_vec[i]);
      GROWTH_RESULT& GROWTHres_i = dynamic_cast<GROWTH_RESULT&>(*p_result_vec[i]);
      int nlayers_i = x_i.soil.getNlayers();
      std::vector<double> lateralFlows_c(nlayers_i, 0.0);
      growthDay_inner_c(GROWTHres_i, GROWTHcomm, x_i,
                        date,
                        *p_weather_vec[i],
                        latitude_vec[i], p_topo_vec[i]->elevation, p_topo_vec[i]->slope, p_topo_vec[i]->aspect,
                        runon,
                        lateralFlows_c, waterTableDepth);
    } else {
      throw medfate::MedfateInternalError("Wrong input class for simulation launching");
    }
  }
}