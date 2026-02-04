#include <RcppArmadillo.h>
#include "medfate.h"
#include "modelInput_c.h"
#include <vector>
#include <string>
#include <limits>
#include <cmath>

#ifndef COMMUNICATION_STRUCTURES_C_H
#define COMMUNICATION_STRUCTURES_C_H

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
  }
};
Rcpp::NumericVector copyWeather_c(const WeatherInputVector& meteo);

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

Rcpp::NumericMatrix copyNumericMatrix_c(arma::mat comm, int rows, int cols);


#endif