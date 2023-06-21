#include <Rcpp.h>

#ifndef SPWB_H
#define SPWB_H
#endif
using namespace Rcpp;

void checkspwbInput(List x, String transpirationMode);

CharacterVector getWeatherDates(DataFrame meteo);

DataFrame defineWaterBalanceDailyOutput(CharacterVector dateStrings, NumericVector PET, String transpirationMode);
DataFrame defineSoilWaterBalanceDailyOutput(CharacterVector dateStrings, List soil, String transpirationMode);
DataFrame defineEnergyBalanceDailyOutput(CharacterVector dateStrings);
DataFrame defineTemperatureDailyOutput(CharacterVector dateStrings);
DataFrame defineFireHazardOutput(CharacterVector dateStrings);
NumericMatrix defineTemperatureLayersDailyOutput(CharacterVector dateStrings, DataFrame canopy);
List defineSunlitShadeLeavesDailyOutput(CharacterVector dateStrings, DataFrame above);
List definePlantWaterDailyOutput(CharacterVector dateStrings, DataFrame above, List soil, List control);

void fillWaterBalanceDailyOutput(DataFrame DWB, List sDay, int iday, String transpirationMode);
void fillSoilWaterBalanceDailyOutput(DataFrame SWB, List soil, List sDay, 
                                     int iday, int numDays, String transpirationMode,
                                     String soilFunctions);
void fillEnergyBalanceTemperatureDailyOutput(DataFrame DEB, DataFrame DT, NumericMatrix DLT, List sDay, 
                                             int iday, bool multiLayerBalance);
void fillPlantWaterDailyOutput(List x, List sunlit, List shade, List sDay, int day, String transpirationMode);
void fillFireHazardOutput(DataFrame fireHazard, List sDay, int iday);

void printWaterBalanceResult(DataFrame DWB, List plantDWOL, List x,
                             NumericVector initialPlantContent, NumericVector initialSoilContent, double initialSnowContent,
                             String transpirationMode);

NumericVector fccsHazard(List x, NumericVector meteovec, List transp, double slope);

List spwbDay(List x, CharacterVector date, NumericVector meteovec,
             double latitude, double elevation, double slope, double aspect,  
             double runon=0.0, bool modifyInput = true);
List spwbDay_basic(List x, NumericVector meteovec, 
                   double elevation, double slope, double aspect, 
                   double runon=0.0, bool verbose=false);
List spwbDay_advanced(List x, NumericVector meteovec, 
                      double latitude, double elevation, double slope, double aspect,
                      double solarConstant, double delta, 
                      double runon=0.0, bool verbose = false);