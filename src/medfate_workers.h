// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
#include "modelInput_c.h"
#include "spwb_day_c.h"
#include "lowlevel_structures_c.h"
using namespace RcppParallel;


#ifndef MEDFATE_WORKERS_H
#define MEDFATE_WORKERS_H

struct WATERBALANCE_worker : public Worker
{
  // source vectors and date
  WBCommunicationStructures& WBcomm;
  const std::string& date;
  std::vector<std::unique_ptr<WaterBalanceModelInput>>& input_vec;
  std::vector<double>& latitude_vec;
  std::vector<std::unique_ptr<Topography>>& topo_vec;
  std::vector<WeatherInputVector>& weather_vec;
  
  // destination vector
  std::vector<std::unique_ptr<WB_RESULT>>& output_vec;
  
  // initialize with source and destination
  WATERBALANCE_worker(WBCommunicationStructures& WBcomm,
                      std::string& date, 
                      std::vector<std::unique_ptr<WaterBalanceModelInput>>& input_vec,
                      std::vector<double>& latitude_vec,
                      std::vector<std::unique_ptr<Topography>>& topo_vec,
                      std::vector<WeatherInputVector>& weather_vec,
                      std::vector<std::unique_ptr<WB_RESULT>>& output_vec) : 
    WBcomm(WBcomm), 
    date(date), input_vec(input_vec), 
    latitude_vec(latitude_vec), 
    topo_vec(topo_vec), 
    weather_vec(weather_vec), 
    output_vec(output_vec) {} 
  
  // take the square root of the range of elements requested
  void operator()(std::size_t begin, std::size_t end);
};

#endif