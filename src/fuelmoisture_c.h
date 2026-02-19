#include "RcppArmadillo.h"

#ifndef FUELMOISTURE_C_H
#define FUELMOISTURE_C_H

struct FuelConditions {
  double temperature, humidity;
};

void fuelConditions_c(FuelConditions& fc, double airTemp, double airHumidity, double fuelRadiation, double fuelWindSpeed);
double EMCSimard_c(double fuelTemperature, double fuelHumidity);
double fine1hday_c(double m0, double fuelTemp, double fuelHumidity, double fuelWindSpeed, double effRain);
double coarse10hday_c(double m0, 
                      double prevFuelTempMax, double prevFuelHumidityMin,
                      double currFuelTempMin, double currFuelHumidityMax,
                      double rainDuration);
double coarse100hday_c(double m0, double fuelTempMin, double fuelHumidityMax, 
                       double fuelTempMax, double fuelHumidityMin, 
                       double numSunHours, double rainDuration);
#endif