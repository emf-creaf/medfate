#include <RcppArmadillo.h>
#include "medfate.h"
#include "modelInput_c.h"
#include <vector>
#include <string>
#include <limits>
#include <cmath>

#ifndef LOWLEVEL_STRUCTURES_C_H
#define LOWLEVEL_STRUCTURES_C_H

struct WeatherInputVector {
  double pet;
  double rhmax;
  double rhmin;
  double tmax;
  double tmin;
  double tmaxPrev;
  double tminPrev;
  double tminNext;
  double tday;
  double Catm;
  double Patm;
  double rint;
  double prec;
  double wind;
  double rad;
  double pfire;
  WeatherInputVector() {
    pet = medfate::NA_DOUBLE;
    rhmax = medfate::NA_DOUBLE;
    rhmin = medfate::NA_DOUBLE;
    tmax = medfate::NA_DOUBLE;
    tmin = medfate::NA_DOUBLE;
    tmaxPrev = medfate::NA_DOUBLE;
    tminPrev = medfate::NA_DOUBLE;
    tminNext = medfate::NA_DOUBLE;
    tday = medfate::NA_DOUBLE;
    Catm = medfate::NA_DOUBLE;
    Patm = medfate::NA_DOUBLE;
    rint = medfate::NA_DOUBLE;
    prec = medfate::NA_DOUBLE;
    wind = medfate::NA_DOUBLE;
    rad = medfate::NA_DOUBLE;
    pfire = medfate::NA_DOUBLE;
  }
  WeatherInputVector(Rcpp::NumericVector meteovec) {
    if(meteovec.containsElementNamed("MinTemperature")) tmin = meteovec["MinTemperature"];
    else tmin = medfate::NA_DOUBLE;
  
    if(meteovec.containsElementNamed("MaxTemperature")) tmax = meteovec["MaxTemperature"];
    else tmax = medfate::NA_DOUBLE;
  
    tminNext = medfate::NA_DOUBLE;
    tminPrev = medfate::NA_DOUBLE;
    tmaxPrev = medfate::NA_DOUBLE;
    tday = medfate::NA_DOUBLE;
    
    if(meteovec.containsElementNamed("Precipitation")) prec = meteovec["Precipitation"];
    else prec = medfate::NA_DOUBLE;
  
    if(meteovec.containsElementNamed("MaxRelativeHumidity")) rhmax = meteovec["MaxRelativeHumidity"];
    else rhmax = medfate::NA_DOUBLE;
  
    if(meteovec.containsElementNamed("MinRelativeHumidity")) rhmin = meteovec["MinRelativeHumidity"];
    else rhmin = medfate::NA_DOUBLE;
  
    if(meteovec.containsElementNamed("Radiation")) rad = meteovec["Radiation"];
    else rad = medfate::NA_DOUBLE;
  
    if(meteovec.containsElementNamed("WindSpeed")) wind = meteovec["WindSpeed"];
    else wind = medfate::NA_DOUBLE;
  
    if(meteovec.containsElementNamed("CO2")) Catm = meteovec["CO2"];
    else Catm = medfate::NA_DOUBLE;
    
    if(meteovec.containsElementNamed("Patm")) Patm = meteovec["Patm"];
    else Patm = medfate::NA_DOUBLE;
    
    if(meteovec.containsElementNamed("RainfallIntensity")) rint = meteovec["RainfallIntensity"];
    else rint = medfate::NA_DOUBLE;
    
    if(meteovec.containsElementNamed("PET")) pet = meteovec["PET"];
    else pet = medfate::NA_DOUBLE;
    
  }
};
Rcpp::NumericVector copyWeather_c(const WeatherInputVector& meteo, const std::string& transpirationMode);

struct Topography {
  double elevation;
  double slope;
  double aspect;
  Topography() {
    elevation = medfate::NA_DOUBLE;
    slope = medfate::NA_DOUBLE;
    aspect = medfate::NA_DOUBLE;
  }
};
Rcpp::NumericVector copyTopo_c(const Topography& topo);

struct StandWB_RESULT {
  double PET, Rain, Snow;
  double NetRain, Snowmelt, Runon;
  double Infiltration, InfiltrationExcess, SaturationExcess, Runoff;
  double DeepDrainage, CapillarityRise, SoilEvaporation, HerbTranspiration;
  double PlantExtraction, Transpiration, HydraulicRedistribution;
  
  StandWB_RESULT() {
    PET = Rain = Snow = medfate::NA_DOUBLE;
    NetRain = Snowmelt = Runon = medfate::NA_DOUBLE;
    Infiltration = InfiltrationExcess = SaturationExcess = Runoff = medfate::NA_DOUBLE;
    DeepDrainage = CapillarityRise = SoilEvaporation = HerbTranspiration = medfate::NA_DOUBLE;
    PlantExtraction = Transpiration = HydraulicRedistribution = medfate::NA_DOUBLE;
  }
  
};
Rcpp::NumericVector copyWaterBalanceResult_c(const StandWB_RESULT& SWBres);

struct Soil_RESULT {
  std::vector<double> Psi, HerbTranspiration, HydraulicInput, HydraulicOutput, PlantExtraction;
  Soil_RESULT(size_t nlayers) {
    Psi = std::vector<double>(nlayers, medfate::NA_DOUBLE);
    HerbTranspiration = std::vector<double>(nlayers, medfate::NA_DOUBLE);
    HydraulicInput = std::vector<double>(nlayers, medfate::NA_DOUBLE);
    HydraulicOutput = std::vector<double>(nlayers, medfate::NA_DOUBLE);
    PlantExtraction = std::vector<double>(nlayers, medfate::NA_DOUBLE);
  }
};
Rcpp::DataFrame copySoilResult_c(const Soil_RESULT& Soil);

struct Stand_RESULT {
  double LAI, LAIherb, LAIlive, LAIexpanded, LAIdead;
  double Cm, LgroundPAR, LgroundSWR;
  Stand_RESULT() {
    LAI = LAIherb = LAIlive = LAIexpanded = LAIdead = medfate::NA_DOUBLE;
    Cm = LgroundPAR = LgroundSWR = medfate::NA_DOUBLE;
  }
};
Rcpp::NumericVector copyStandResult_c(const Stand_RESULT& Sres);

struct StandCB_RESULT {
  double GrossPrimaryProduction;
  double MaintenanceRespiration;
  double SynthesisRespiration;
  double NetPrimaryProduction; 
  double HeterotrophicRespiration;
  double FireCombustion;
  double NetEcosystemProduction;
  StandCB_RESULT() {
    GrossPrimaryProduction = MaintenanceRespiration = SynthesisRespiration = medfate::NA_DOUBLE;
    NetPrimaryProduction = HeterotrophicRespiration = FireCombustion = medfate::NA_DOUBLE;
    NetEcosystemProduction = medfate::NA_DOUBLE;
  }
};
Rcpp::NumericVector copyCarbonBalanceResult_c(const StandCB_RESULT& SWBres);

Rcpp::NumericMatrix copyNumericMatrix_c(arma::mat comm, int rows, int cols);


#endif