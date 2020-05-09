#include <Rcpp.h>

#ifndef SPWB_H
#define SPWB_H
#endif
using namespace Rcpp;

void checkspwbInput(List x, List soil, String transpirationMode);

DataFrame defineWaterBalanceDailyOutput(DataFrame meteo, NumericVector PET, String transpirationMode);
DataFrame defineSoilWaterBalanceDailyOutput(DataFrame meteo, List soil, String transpirationMode);
DataFrame defineEnergyBalanceDailyOutput(DataFrame meteo);
DataFrame defineTemperatureDailyOutput(DataFrame meteo);
List defineSunlitShadeLeavesDailyOutput(DataFrame meteo, DataFrame above);
List definePlantWaterDailyOutput(DataFrame meteo, DataFrame above, List soil, List control);

void fillWaterBalanceDailyOutput(DataFrame DWB, List sDay, int iday, String transpirationMode);
void fillSoilWaterBalanceDailyOutput(DataFrame SWB, List soil, List sDay, 
                                     int iday, int numDays, String transpirationMode,
                                     String soilFunctions);
void fillEnergyBalanceTemperatureDailyOutput(DataFrame DEB, DataFrame DT, List sDay, int iday);
void fillPlantWaterDailyOutput(List x, List sunlit, List shade, List sDay, int day, String transpirationMode);

void printWaterBalanceResult(DataFrame DWB, List plantDWOL, 
                             List soil, String soilFunctions,
                             NumericVector initialContent, double initialSnowContent,
                             String transpirationMode);

List spwbDay(List x, List soil, CharacterVector date, double tmin, double tmax, 
             double rhmin, double rhmax, double rad, double wind, 
             double latitude, double elevation, double slope, double aspect,  
             double prec, double runon=0.0);
List spwbDay1(List x, List soil, double tday, double pet, double prec, double er, double runon=0.0, 
              double rad = NA_REAL, double elevation = NA_REAL, bool verbose=false);
List spwbDay2(List x, List soil, double tmin, double tmax, double tminPrev, double tmaxPrev, double tminNext, 
              double rhmin, double rhmax, double rad, double wind, 
              double latitude, double elevation, double slope, double aspect,
              double solarConstant, double delta, 
              double prec, double pet, double er, double runon=0.0, bool verbose = false);