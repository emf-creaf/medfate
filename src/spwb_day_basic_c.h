#include "Rcpp.h"
#include "communication_structures_c.h"
#include "lightextinction_basic_c.h"
#include "transpiration_basic_c.h"
#include "modelInput_c.h"
#include "firebehaviour_c.h"
#include "hydrology_c.h"

#ifndef SPWB_DAY_BASIC_C_H
#define SPWB_DAY_BASIC_C_H

struct BasicSPWB_RESULT {
  Topography topo;
  WeatherInputVector meteovec;
  StandWB_RESULT WaterBalance;
  Soil_RESULT Soil;
  Stand_RESULT Stand;
  PlantsBasicTranspiration_RESULT Plants;
  FCCS_RESULT fccs;
  
  BasicSPWB_RESULT(PlantsBasicTranspiration_RESULT& PlantsIN, size_t nlayers) : Soil(nlayers) {
    Plants = PlantsIN;
  }
};
Rcpp::List copyBasicSPWBResult_c(const BasicSPWB_RESULT& BSPWBres, ModelInput& x);

struct BasicSPWB_COMM {
  WaterInputs_COMM waterInputs;
  SoilWaterBalance_COMM SWBcomm;
  BasicTranspiration_COMM BTcomm;
  BasicSPWB_COMM(size_t numCohorts = 0, size_t ncanlayers = 0, size_t nlayers= 0,
                 std::string soilDomains = "buckets") : 
    BTcomm(numCohorts, ncanlayers, nlayers),
    SWBcomm(nlayers, soilDomains){}
};

void spwbDay_basic_c(BasicSPWB_RESULT& BTres, BasicSPWB_COMM& BT_comm, ModelInput& x, 
                     const std::vector<double>& meteovec, const double elevation, const double slope, const double aspect,
                     const double runon, 
                     const std::vector<double>& lateralFlows, const double waterTableDepth);

#endif