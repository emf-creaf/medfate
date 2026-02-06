// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
#include "modelInput_c.h"
#include "spwb_day_c.h"
#include "lowlevel_structures_c.h"
using namespace RcppParallel;


#ifndef MEDFATE_WORKERS_H
#define MEDFATE_WORKERS_H

struct SPWB_worker : public Worker
{
  // source vectors and date
  SPWBCommunicationStructures SPWBcomm;
  std::string date;
  std::vector<ModelInput> input_vec;
  std::vector<double> latitude_vec;
  std::vector<Topography> topo_vec;
  std::vector<WeatherInputVector> weather_vec;
  
  // destination vector
  std::vector<SPWB_RESULT> output_vec;
  
  // initialize with source and destination
  SPWB_worker(SPWBCommunicationStructures& SPWBcomm,
             std::string& date, 
             std::vector<ModelInput>& input_vec,
             std::vector<double>& latitude_vec,
             std::vector<Topography>& topo_vec,
             std::vector<WeatherInputVector>& weather_vec,
             std::vector<SPWB_RESULT>& output_vec) : 
    SPWBcomm(SPWBcomm), 
    date(date), input_vec(input_vec), 
    latitude_vec(latitude_vec), 
    topo_vec(topo_vec), 
    weather_vec(weather_vec), 
    output_vec(output_vec) {} 
  
  // take the square root of the range of elements requested
  void operator()(std::size_t begin, std::size_t end);
};

#endif