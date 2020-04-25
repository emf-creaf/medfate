#include <Rcpp.h>

#ifndef SPWB_H
#define SPWB_H
#endif
using namespace Rcpp;

void checkspwbInput(List x, List soil, String transpirationMode);

List definePlantWaterDailyOutputList(DataFrame meteo, DataFrame above, List soil, List control);
void fillPlantWaterDailyOutputList(List x, List sDay, int day, String transpirationMode);
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