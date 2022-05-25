#include <Rcpp.h>

#ifndef SPWB_H
#define SPWB_H
#endif
using namespace Rcpp;

void checkspwbInput(List x, String transpirationMode);

DataFrame defineWaterBalanceDailyOutput(DataFrame meteo, NumericVector PET, String transpirationMode);
DataFrame defineSoilWaterBalanceDailyOutput(DataFrame meteo, List soil, String transpirationMode);
DataFrame defineEnergyBalanceDailyOutput(DataFrame meteo);
DataFrame defineTemperatureDailyOutput(DataFrame meteo);
NumericMatrix defineTemperatureLayersDailyOutput(DataFrame meteo, DataFrame canopy);
List defineSunlitShadeLeavesDailyOutput(DataFrame meteo, DataFrame above);
List definePlantWaterDailyOutput(DataFrame meteo, DataFrame above, List soil, List control);

void fillWaterBalanceDailyOutput(DataFrame DWB, List sDay, int iday, String transpirationMode);
void fillSoilWaterBalanceDailyOutput(DataFrame SWB, List soil, List sDay, 
                                     int iday, int numDays, String transpirationMode,
                                     String soilFunctions);
void fillEnergyBalanceTemperatureDailyOutput(DataFrame DEB, DataFrame DT, NumericMatrix DLT, List sDay, 
                                             int iday, bool multiLayerBalance);
void fillPlantWaterDailyOutput(List x, List sunlit, List shade, List sDay, int day, String transpirationMode);

void printWaterBalanceResult(DataFrame DWB, List plantDWOL, 
                             List soil, String soilFunctions,
                             NumericVector initialContent, double initialSnowContent,
                             String transpirationMode);

List spwbDay(List x, CharacterVector date, double tmin, double tmax, 
             double rhmin, double rhmax, double rad, double wind, 
             double latitude, double elevation, double slope, double aspect,  
             double prec, double runon=0.0, bool modifyInput = true);
List spwbDay1(List x, NumericVector meteovec, 
              double elevation, double slope, double aspect, 
              double runon=0.0, bool verbose=false);
List spwbDay2(List x, NumericVector meteovec, 
              double latitude, double elevation, double slope, double aspect,
              double solarConstant, double delta, 
              double runon=0.0, bool verbose = false);