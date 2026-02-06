// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
#include "medfate_workers.h"
using namespace RcppParallel;

// Implements simulation for the range of elements 
void SPWB_worker::operator()(std::size_t begin, std::size_t end) {
  double runon = 0.0;
  double waterTableDepth = 0.0;
  for(size_t i = begin; i< end; i++) {
    std::vector<double> lateralFlows_c(input_vec[i].soil.getNlayers(), 0.0);
    spwbDay_inner_c(SPWBcomm, input_vec[i], 
                    date[0].c_string(),
                    weather_vec[i],
                    latitude_vec[i], topo_vec[i].elevation, topo_vec[i].slope, topo_vec[i].aspect,
                    runon, 
                    lateralFlows_c, waterTableDepth);
    if(input_vec[i].control.transpirationMode=="Granier") {
      output_vec[i] = SPWBcomm.BSPWBres;
    } else {
      output_vec[i] = SPWBcomm.ASPWBres;
    }
  }
}
