// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
#include "medfate_workers.h"
using namespace RcppParallel;

// Implements simulation for the range of elements 
void WATERBALANCE_worker::operator()(std::size_t begin, std::size_t end) {
  double runon = 0.0;
  double waterTableDepth = 0.0;
  for(size_t i = begin; i< end; i++) {
    std::vector<double> lateralFlows_c(input_vec[i]->soil.getNlayers(), 0.0);

    //Call simulation
    wb_day_inner_c(*output_vec[i], WBcomm, *input_vec[i],
                   date,
                   weather_vec[i],
                   latitude_vec[i], topo_vec[i]->elevation, topo_vec[i]->slope, topo_vec[i]->aspect,
                   runon,
                   lateralFlows_c, waterTableDepth);
  }
}
