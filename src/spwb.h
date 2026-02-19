#include <RcppArmadillo.h>
#include "spwb_day_c.h"
#include "soil_c.h"

#ifndef SPWB_H
#define SPWB_H
#endif
using namespace Rcpp;

void checkspwbInput(List x, String transpirationMode);

CharacterVector getWeatherDates(DataFrame meteo);

DataFrame defineStandDailyOutput(CharacterVector dateStrings);
DataFrame defineWaterBalanceDailyOutput(CharacterVector dateStrings, String transpirationMode);
List defineSoilDailyOutput(CharacterVector dateStrings, DataFrame soil, bool includePlants = true);
DataFrame defineSnowDailyOutput(CharacterVector dateStrings);
DataFrame defineEnergyBalanceDailyOutput(CharacterVector dateStrings);
DataFrame defineTemperatureDailyOutput(CharacterVector dateStrings);
DataFrame defineFireHazardOutput(CharacterVector dateStrings);
NumericMatrix defineTemperatureLayersDailyOutput(CharacterVector dateStrings, DataFrame canopy);
List defineSunlitShadeLeavesDailyOutput(CharacterVector dateStrings, DataFrame above);
List definePlantWaterDailyOutput(CharacterVector dateStrings, DataFrame above, DataFrame soil, List control);

void fillStandDailyOutput(DataFrame Stand, List sDay, int iday);
void fillStandDailyOutput_c(DataFrame Stand, Stand_RESULT& stand, int iday);
void fillWaterBalanceDailyOutput(DataFrame DWB, List sDay, int iday, String transpirationMode);
void fillWaterBalanceDailyOutput_c(DataFrame DWB, StandWB_RESULT& db, int iday);
void fillSoilDailyOutput(List Soil, DataFrame soil, List sDay, 
                         int iday, int numDays, String soilFunctions,
                         bool includePlants = true);
void fillSoilDailyOutput_c(List SWB, Soil& soil, Soil_RESULT& sb, 
                           int iday, int numDays,
                           bool includePlants = true);
void fillSoilPoolDailyOutput_c(List soilPools, Soil& soil, arma::mat& Wpool, int iday);
void fillSnowDailyOutput(DataFrame Snow, List x, int iday);
void fillSnowDailyOutput_c(DataFrame Snow, WaterBalanceModelInput& x, int iday);
void fillEnergyBalanceDailyOutput(DataFrame DEB, List sDay, int iday, int ntimesteps);
void fillEnergyBalanceDailyOutput_c(DataFrame DEB, EnergyBalance_RESULT& EB, int iday, int ntimesteps);
void fillTemperatureDailyOutput(DataFrame DT, List sDay, int iday, int ntimesteps);
void fillTemperatureDailyOutput_c(DataFrame DT, EnergyBalance_RESULT& EB, int iday, int ntimesteps);
void fillTemperatureLayersDailyOutput(NumericMatrix DLT, List sDay, int iday, int ncanlayers, int ntimesteps);
void fillTemperatureLayersDailyOutput_c(NumericMatrix DLT, EnergyBalance_RESULT& EB, int iday, int ncanlayers, int ntimesteps);
void fillPlantWaterDailyOutput(List x, List sDay, int iday, String transpirationMode);
void fillPlantWaterDailyOutput_c(List x, SPWB_RESULT& sDay, int iday, String transpirationMode, int ntimesteps);
void fillSunlitShadeLeavesDailyOutput(List sunlit, List shade, List sDay, int day, int numCohorts);
void fillSunlitShadeLeavesDailyOutput_c(List sunlit, List shade, AdvancedTranspiration_RESULT& sDay, int iday, int numCohorts);
void fillFireHazardOutput(DataFrame fireHazard, List sDay, int iday);
void fillFireHazardOutput_c(DataFrame fireHazard, FCCS_RESULT& fhd, int iday);
void printWaterBalanceResult(List outputList, List x,
                             NumericVector initialPlantContent, NumericVector initialSoilContent, double initialSnowContent,
                             String transpirationMode);
