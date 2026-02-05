#include "Rcpp.h"
#include "communication_structures_c.h"
#include "lightextinction_basic_c.h"
#include "transpiration_basic_c.h"
#include "modelInput_c.h"
#include "firebehaviour_c.h"
#include "hydrology_c.h"

#ifndef SPWB_DAY_C_H
#define SPWB_DAY_C_H

struct BasicSPWB_RESULT {
  Topography topo;
  WeatherInputVector meteovec;
  StandWB_RESULT WaterBalance;
  Soil_RESULT Soil;
  Stand_RESULT Stand;
  SoilWaterBalance_RESULT SWBres;
  BasicTranspiration_RESULT BTres;
  FCCSBehaviour_RESULT fccsbeh;
  FCCS_RESULT fccs;
  
  BasicSPWB_RESULT(BasicTranspiration_RESULT& BTresIN, size_t nlayers) : Soil(nlayers) {
    BTres = BTresIN;
  }
};
Rcpp::List copyBasicSPWBResult_c(const BasicSPWB_RESULT& BSPWBres, ModelInput& x);

struct BasicSPWB_COMM {
  WaterInputs_COMM waterInputs;
  BasicTranspiration_COMM BTcomm;
  SoilWaterBalance_COMM SWBcomm;
  BasicSPWB_COMM(size_t numCohorts = 0, size_t ncanlayers = 0, size_t nlayers= 0,
                 std::string soilDomains = "buckets") : 
    BTcomm(numCohorts, ncanlayers, nlayers),
    SWBcomm(nlayers, soilDomains){}
};

void spwbDay_basic_c(BasicSPWB_RESULT& BTres, BasicSPWB_COMM& BT_comm, ModelInput& x, 
                     const WeatherInputVector& meteovec, const double elevation, const double slope, const double aspect,
                     const double runon, 
                     const std::vector<double>& lateralFlows, const double waterTableDepth);

#endif