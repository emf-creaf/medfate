#include "Rcpp.h"
#include "communication_structures_c.h"
#include "lightextinction_basic_c.h"
#include "transpiration_basic_c.h"
#include "modelInput_c.h"

#ifndef SPWB_DAY_BASIC_C_H
#define SPWB_DAY_BASIC_C_H

struct BasicSPWB_RESULT {
  Topography topo;
  WeatherInputVector meteovec;
  StandWB_RESULT WaterBalance;
  Soil_RESULT Soil;
  Stand_RESULT Stand;
  PlantsBasicTranspiration_RESULT Plants;
  
  BasicSPWB_RESULT(PlantsBasicTranspiration_RESULT& PlantsIN, size_t nlayers) : Soil(nlayers) {
    Plants = PlantsIN;
  }
};
Rcpp::List copyBasicSPWBResult_c(const BasicSPWB_RESULT& BSPWBres, ModelInput& x);

struct BasicSPWB_COMM {
  BasicTranspiration_COMM BTcomm;
};


#endif