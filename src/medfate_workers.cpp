// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
#include "medfate_workers.h"
using namespace RcppParallel;

// Implements simulation for the range of elements 
void SPWB_worker::operator()(std::size_t begin, std::size_t end) {
  double runon = 0.0;
  double waterTableDepth = 0.0;
  for(size_t i = begin; i< end; i++) {
    int numCohorts_i = input_vec[i]->cohorts.CohortCode.size();
    int nlayers_i = input_vec[i]->soil.getNlayers();
    int ncanlayers_i = input_vec[i]->canopy.zlow.size();
    int ntimesteps_i = input_vec[i]->control.advancedWB.ndailysteps;
    
    std::vector<double> lateralFlows_c(input_vec[i]->soil.getNlayers(), 0.0);
    //Initialises a result
    if(input_vec[i]->control.transpirationMode=="Granier") {
      BasicTranspiration_RESULT BTres_i(numCohorts_i, nlayers_i);
      BasicSPWB_RESULT BSPWBres(BTres_i);
      *output_vec[i] = BSPWBres;
    } else {
      AdvancedTranspiration_RESULT ATres_i(numCohorts_i, nlayers_i, ncanlayers_i, ntimesteps_i);
      AdvancedSPWB_RESULT ASPWBres(ATres_i);
      *output_vec[i] = ASPWBres;
    }

    //Call simulation
    spwbDay_inner_c(*output_vec[i], SPWBcomm, *input_vec[i], 
                    date,
                    weather_vec[i],
                    latitude_vec[i], topo_vec[i]->elevation, topo_vec[i]->slope, topo_vec[i]->aspect,
                    runon, 
                    lateralFlows_c, waterTableDepth);
  }
}
