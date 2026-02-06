#include "Rcpp.h"
#include "lowlevel_structures_c.h"
#include "lightextinction_basic_c.h"
#include "transpiration_basic_c.h"
#include "hydrology_c.h"
#include "transpiration_advanced_c.h"
#include "modelInput_c.h"
#include "firebehaviour_c.h"
#include "hydrology_c.h"

#ifndef SPWB_DAY_C_H
#define SPWB_DAY_C_H

struct SPWB_RESULT {
  Topography topo;
  WeatherInputVector meteovec;
  Soil_RESULT Soil;
  Stand_RESULT Stand;
  StandWB_RESULT WaterBalance;
  SoilWaterBalance_RESULT SWBres;
  FCCSBehaviour_RESULT fccsbeh;
  FCCS_RESULT fccs;
  
  SPWB_RESULT(size_t nlayers) : Soil(nlayers) {}
  
  virtual ~SPWB_RESULT() = default;
};
Rcpp::List copySPWBResult_c(const SPWB_RESULT& SPWBres, ModelInput& x);

struct BasicSPWB_RESULT : SPWB_RESULT {
  BasicTranspiration_RESULT BTres;
  
  BasicSPWB_RESULT(BasicTranspiration_RESULT& BTres, size_t nlayers) : SPWB_RESULT(nlayers), BTres(BTres) {}
};
Rcpp::List copyBasicSPWBResult_c(const BasicSPWB_RESULT& BSPWBres, ModelInput& x);

struct AdvancedSPWB_RESULT : SPWB_RESULT {
  AdvancedTranspiration_RESULT ATres;
  
  AdvancedSPWB_RESULT(AdvancedTranspiration_RESULT& ATresIN, size_t nlayers) : SPWB_RESULT(nlayers), ATres(ATresIN) {}
};
Rcpp::List copyAdvancedSPWBResult_c(const AdvancedSPWB_RESULT& ASPWBres, ModelInput& x);

struct BasicSPWB_COMM {
  WaterInputs_COMM waterInputs;
  SoilWaterBalance_COMM SWBcomm;
  BasicTranspiration_COMM BTcomm;
  BasicSPWB_COMM(SoilWaterBalance_COMM& SWBcommIn, BasicTranspiration_COMM& BTcommIn) : 
    SWBcomm(SWBcommIn), BTcomm(BTcommIn) {
  }
};

void spwbDay_basic_c(BasicSPWB_RESULT& BTres, BasicSPWB_COMM& BT_comm, ModelInput& x, 
                     const WeatherInputVector& meteovec, const double elevation, const double slope, const double aspect,
                     const double runon, 
                     const std::vector<double>& lateralFlows, const double waterTableDepth);



struct AdvancedSPWB_COMM {
  WaterInputs_COMM waterInputs;
  SoilWaterBalance_COMM SWBcomm;
  AdvancedTranspiration_COMM ATcomm;
  AdvancedSPWB_COMM(SoilWaterBalance_COMM& SWBcommIn, AdvancedTranspiration_COMM& ATcommIn) : 
    SWBcomm(SWBcommIn), ATcomm(ATcommIn) {}
};

void spwbDay_advanced_c(AdvancedSPWB_RESULT& ASPWBres, AdvancedSPWB_COMM& ASPWB_comm, ModelInput& x, 
                        const WeatherInputVector& meteovec, 
                        const double latitude, const double elevation, const double slope, const double aspect,
                        const double solarConstant, const double delta, 
                        const double runon, 
                        const std::vector<double>& lateralFlows, const double waterTableDepth);


struct SPWBCommunicationStructures {
  
  SoilWaterBalance_COMM SWBcomm;
  BasicTranspiration_RESULT BTres;
  BasicTranspiration_COMM BTcomm;
  AdvancedTranspiration_RESULT ATres;
  AdvancedTranspiration_COMM ATcomm;
  BasicSPWB_RESULT BSPWBres;
  AdvancedSPWB_RESULT ASPWBres;
  BasicSPWB_COMM BSPWBcomm;
  AdvancedSPWB_COMM ASPWBcomm;
  
  SPWBCommunicationStructures(size_t numCohorts, size_t nlayers, size_t ncanlayers, size_t ntimesteps) :
    SWBcomm(nlayers),
    BTres(numCohorts, nlayers),
    BTcomm(numCohorts, ncanlayers, nlayers),
    ATres(numCohorts, nlayers, ncanlayers, ntimesteps),
    ATcomm(numCohorts, nlayers, ncanlayers, ntimesteps),
    BSPWBres(BTres, nlayers),
    ASPWBres(ATres, nlayers),
    BSPWBcomm(SWBcomm, BTcomm),
    ASPWBcomm(SWBcomm, ATcomm) {} 
};

void spwbDay_inner_c(SPWBCommunicationStructures& SPWBcomm, ModelInput& x, 
                     std::string date,
                     WeatherInputVector meteovec, 
                     double latitude, double elevation, double slope, double aspect,
                     const double runon, 
                     const std::vector<double>& lateralFlows, const double waterTableDepth);
#endif
