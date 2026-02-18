#include "Rcpp.h"
#include "medfate.h"
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

struct NSWB_RESULT : public ABSTRACTMODEL_RESULT {
  StandWB_RESULT WaterBalance;
};
Rcpp::List copyNSWBResult_c(NSWB_RESULT& WBres, NonSoilWaterBalanceModelInput& x);


void nswbDay_c(NSWB_RESULT& NSWBres, NonSoilWaterBalanceModelInput& x,
               const WeatherInputVector& meteovec, 
               const double elevation, const double slope, const double aspect,
               const double runon, 
               const double rock_max_infiltration);
  
  
struct WB_RESULT : public ABSTRACTMODEL_RESULT {
  Topography topo;
  WeatherInputVector meteovec;
  Soil_RESULT Soil;
  WB_RESULT() : Soil(0) {}
  WB_RESULT(size_t nlayers) : Soil(nlayers) {}
  virtual ~WB_RESULT() = default;
};
Rcpp::List copyWBResult_c(WB_RESULT& WBres, WaterBalanceModelInput& x);

struct AgricultureWB_RESULT: public WB_RESULT {
  SoilWaterBalance_RESULT SWBres;
  StandWB_RESULT WaterBalance;
  AgricultureWB_RESULT() : WB_RESULT(0) {}
  AgricultureWB_RESULT(size_t nlayers) : WB_RESULT(nlayers) {}
};
Rcpp::List copyAgricultureWBResult_c(AgricultureWB_RESULT& AgrWBres, AgricultureModelInput& x);

struct SPWB_RESULT : public WB_RESULT {
  Stand_RESULT Stand;
  StandWB_RESULT WaterBalance;
  SoilWaterBalance_RESULT SWBres;
  FCCSBehaviour_RESULT fccsbeh;
  FCCS_RESULT fccs;
  
  SPWB_RESULT() : WB_RESULT(0) {}
  SPWB_RESULT(size_t nlayers) : WB_RESULT(nlayers) {}
  virtual ~SPWB_RESULT() = default;
};
Rcpp::List copySPWBResult_c(SPWB_RESULT& SPWBres, ModelInput& x);

struct BasicSPWB_RESULT : public SPWB_RESULT {
  BasicTranspiration_RESULT BTres;
  
  BasicSPWB_RESULT(BasicTranspiration_RESULT& BTres, size_t nlayers) : SPWB_RESULT(nlayers), BTres(BTres) {}
};
Rcpp::List copyBasicSPWBResult_c(const BasicSPWB_RESULT& BSPWBres, ModelInput& x);

struct AdvancedSPWB_RESULT : public SPWB_RESULT {
  AdvancedTranspiration_RESULT ATres;
  
  AdvancedSPWB_RESULT(AdvancedTranspiration_RESULT& ATresIN, size_t nlayers) : SPWB_RESULT(nlayers), ATres(ATresIN) {}
};
Rcpp::List copyAdvancedSPWBResult_c(const AdvancedSPWB_RESULT& ASPWBres, ModelInput& x);

struct AgricultureWB_COMM : public AbstractCommunicationStructures {
  WaterInputs_COMM waterInputs;
  SoilWaterBalance_COMM SWBcomm;
  AgricultureWB_COMM(SoilWaterBalance_COMM& SWBcommIn) : SWBcomm(SWBcommIn) {}
};

void aspwbDay_c(AgricultureWB_RESULT& AgrWBres, AgricultureWB_COMM& AgrWBcomm, AgricultureModelInput& x, 
                const WeatherInputVector& meteovec, 
                const double elevation, const double slope, const double aspect,
                const double runon, 
                const std::vector<double>& lateralFlows, const double waterTableDepth);

struct BasicSPWB_COMM : public AbstractCommunicationStructures {
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


struct WBCommunicationStructures : public AbstractCommunicationStructures {
  
  SoilWaterBalance_COMM SWBcomm;
  BasicTranspiration_COMM BTcomm;
  AdvancedTranspiration_COMM ATcomm;
  AgricultureWB_COMM AgrWBcomm;
  BasicSPWB_COMM BSPWBcomm;
  AdvancedSPWB_COMM ASPWBcomm;

  WBCommunicationStructures(size_t numCohorts, size_t nlayers, size_t ncanlayers, size_t ntimesteps) :
    SWBcomm(nlayers),
    BTcomm(numCohorts, ncanlayers, nlayers),
    ATcomm(numCohorts, nlayers, ncanlayers, ntimesteps),
    AgrWBcomm(SWBcomm),
    BSPWBcomm(SWBcomm, BTcomm),
    ASPWBcomm(SWBcomm, ATcomm)
     {} 
};

void wb_day_inner_c(WB_RESULT& SPWBres, WBCommunicationStructures& WBcomm, WaterBalanceModelInput& x, 
                    std::string date,
                    WeatherInputVector meteovec, 
                    double latitude, double elevation, double slope, double aspect,
                    const double runon, 
                    const std::vector<double>& lateralFlows, const double waterTableDepth);
#endif
