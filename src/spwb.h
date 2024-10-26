#include <Rcpp.h>

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
void fillWaterBalanceDailyOutput(DataFrame DWB, List sDay, int iday, String transpirationMode);
void fillSoilDailyOutput(List Soil, DataFrame soil, List sDay, 
                         int iday, int numDays, String soilFunctions,
                         bool includePlants = true);
void fillSnowDailyOutput(DataFrame Snow, List x, int iday);
void fillEnergyBalanceDailyOutput(DataFrame DEB, List sDay, int iday, int ntimesteps);
void fillTemperatureDailyOutput(DataFrame DT, List sDay, int iday, int ntimesteps);
void fillTemperatureLayersDailyOutput(NumericMatrix DLT, List sDay, int iday, int ncanlayers, int ntimesteps);
void fillPlantWaterDailyOutput(List x, List sDay, int iday, String transpirationMode);
void fillSunlitShadeLeavesDailyOutput(List sunlit, List shade, List sDay, int day, int numCohorts);
void fillFireHazardOutput(DataFrame fireHazard, List sDay, int iday);

void printWaterBalanceResult(List outputList, List x,
                             NumericVector initialPlantContent, NumericVector initialSoilContent, double initialSnowContent,
                             String transpirationMode);
