// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
#include "medfate.h"
#include "modelInput_c.h"
#include "spwb_day_c.h"
#include "growth_day_c.h"
#include "lowlevel_structures_c.h"
using namespace RcppParallel;


#ifndef MEDFATE_WORKERS_H
#define MEDFATE_WORKERS_H

struct DAY_worker : public Worker
{
  // source vectors and date
  size_t numCohorts_max;
  size_t nlayers_max;
  size_t ncanlayers_max;
  size_t ntimesteps_max;
  const std::string& date;
  std::vector<std::unique_ptr<AbstractModelInput>>& p_x_vec;
  std::vector<double>& latitude_vec;
  std::vector<std::unique_ptr<Topography>>& p_topo_vec;
  std::vector<std::unique_ptr<WeatherInputVector>>& p_weather_vec;
  
  // destination vector
  std::vector<std::unique_ptr<ABSTRACTMODEL_RESULT>>& p_result_vec;
  
  // initialize with source and destination
  DAY_worker(size_t numCohorts_maxIn,
             size_t nlayers_maxIn,
             size_t ncanlayers_maxIn,
             size_t ntimesteps_maxIn,
             std::string& dateIn, 
             std::vector<std::unique_ptr<AbstractModelInput>>& p_x_vecIn,
             std::vector<double>& latitude_vecIn,
             std::vector<std::unique_ptr<Topography>>& p_topo_vecIn,
             std::vector<std::unique_ptr<WeatherInputVector>>& p_weather_vecIn,
             std::vector<std::unique_ptr<ABSTRACTMODEL_RESULT>>& p_result_vecIn) : 
    numCohorts_max(numCohorts_maxIn),
    nlayers_max(nlayers_maxIn),
    ncanlayers_max(ncanlayers_maxIn),
    ntimesteps_max(ntimesteps_maxIn),
    date(dateIn), p_x_vec(p_x_vecIn), 
    latitude_vec(latitude_vecIn), 
    p_topo_vec(p_topo_vecIn), 
    p_weather_vec(p_weather_vecIn), 
    p_result_vec(p_result_vecIn) {} 
  
  // take the square root of the range of elements requested
  void operator()(std::size_t begin, std::size_t end);
};

#endif