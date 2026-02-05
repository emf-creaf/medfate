#include "Rcpp.h"
#include "lowlevel_structures_c.h"
#include "lightextinction_basic_c.h"
#include "transpiration_basic_c.h"
#include "transpiration_advanced_c.h"
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
  BasicSPWB_COMM(BasicTranspiration_COMM& BTcommIn, size_t nlayers= 0, std::string soilDomains = "buckets") : 
    SWBcomm(nlayers, soilDomains) {
    BTcomm = BTcommIn;
  }
};

void spwbDay_basic_c(BasicSPWB_RESULT& BTres, BasicSPWB_COMM& BT_comm, ModelInput& x, 
                     const WeatherInputVector& meteovec, const double elevation, const double slope, const double aspect,
                     const double runon, 
                     const std::vector<double>& lateralFlows, const double waterTableDepth);


struct AdvancedSPWB_RESULT {
  Topography topo;
  WeatherInputVector meteovec;
  StandWB_RESULT WaterBalance;
  Soil_RESULT Soil;
  Stand_RESULT Stand;
  SoilWaterBalance_RESULT SWBres;
  AdvancedTranspiration_RESULT ATres;
  FCCSBehaviour_RESULT fccsbeh;
  FCCS_RESULT fccs;
  
  AdvancedSPWB_RESULT(AdvancedTranspiration_RESULT& ATresIN, size_t nlayers) : Soil(nlayers) {
    ATres = ATresIN;
  }
};
Rcpp::List copyAdvancedSPWBResult_c(const AdvancedSPWB_RESULT& ASPWBres, ModelInput& x);

struct AdvancedSPWB_COMM {
  WaterInputs_COMM waterInputs;
  AdvancedTranspiration_COMM ATcomm;
  SoilWaterBalance_COMM SWBcomm;
  AdvancedSPWB_COMM(AdvancedTranspiration_COMM& ATcommIn, size_t nlayers= 0, size_t ntimesteps = 0,
                 std::string soilDomains = "buckets") : 
    SWBcomm(nlayers, soilDomains){
    ATcomm = ATcommIn;
  }
};

void spwbDay_advanced_c(AdvancedSPWB_RESULT& ASPWBres, AdvancedSPWB_COMM& ASPWB_comm, ModelInput& x, 
                        const WeatherInputVector& meteovec, 
                        const double latitude, const double elevation, const double slope, const double aspect,
                        const double solarConstant, const double delta, 
                        const double runon, 
                        const std::vector<double>& lateralFlows, const double waterTableDepth);


struct SPWBCommunicationStructures {
  
  BasicTranspiration_RESULT BTres;
  BasicTranspiration_COMM BTcomm;
  AdvancedTranspiration_RESULT ATres;
  AdvancedTranspiration_COMM ATcomm;
  BasicSPWB_RESULT BSPWBres;
  AdvancedSPWB_RESULT ASPWBres;
  BasicSPWB_COMM BSPWBcomm;
  AdvancedSPWB_COMM ASPWBcomm;
  
  SPWBCommunicationStructures(size_t numCohorts, size_t nlayers, size_t ncanlayers, size_t ntimesteps,
                              std::string soilDomains) :
    BTres(numCohorts, nlayers),
    BTcomm(numCohorts, ncanlayers, nlayers),
    ATres(numCohorts, nlayers, ncanlayers, ntimesteps),
    ATcomm(numCohorts, nlayers, ncanlayers, ntimesteps),
    BSPWBres(BTres, nlayers),
    ASPWBres(ATres, nlayers),
    BSPWBcomm(BTcomm, nlayers, soilDomains),
    ASPWBcomm(ATcomm, nlayers, ntimesteps, soilDomains) {
  } 
};

void spwbDay_inner_c(SPWBCommunicationStructures& SPWBcomm, ModelInput& x, 
                     std::string date,
                     WeatherInputVector meteovec, 
                     double latitude, double elevation, double slope, double aspect,
                     const double runon, 
                     const std::vector<double>& lateralFlows, const double waterTableDepth);
#endif
