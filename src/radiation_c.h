#include "RcppArmadillo.h"
#include "medfate.h"

#ifndef RADIATION_C_H
#define RADIATION_C_H

struct SunriseSet {
  double sunrise;
  double set;
};

struct DirectDiffuseInstant_RESULT {
  double SolarElevation;
  double Rpot;
  double Rpot_flat;
  double Rg;
  double SWR_direct;
  double SWR_diffuse;
  double PAR_direct;
  double PAR_diffuse;
};

struct DirectDiffuseDay_RESULT {
  std::vector<double> SolarHour;
  std::vector<double> SolarElevation;
  std::vector<double> Rpot;
  std::vector<double> Rpot_flat;
  std::vector<double> Rg;
  std::vector<double> SWR_direct;
  std::vector<double> SWR_diffuse;
  std::vector<double> PAR_direct;
  std::vector<double> PAR_diffuse;
  
  DirectDiffuseDay_RESULT(size_t ntimesteps) {
    SolarHour = std::vector<double>(ntimesteps, medfate::NA_DOUBLE);
    SolarElevation = std::vector<double>(ntimesteps, medfate::NA_DOUBLE);
    Rpot = std::vector<double>(ntimesteps, medfate::NA_DOUBLE);
    Rpot_flat = std::vector<double>(ntimesteps, medfate::NA_DOUBLE);
    Rg = std::vector<double>(ntimesteps, medfate::NA_DOUBLE);
    SWR_direct = std::vector<double>(ntimesteps, medfate::NA_DOUBLE);
    SWR_diffuse = std::vector<double>(ntimesteps, medfate::NA_DOUBLE);
    PAR_direct = std::vector<double>(ntimesteps, medfate::NA_DOUBLE);
    PAR_diffuse = std::vector<double>(ntimesteps, medfate::NA_DOUBLE);
  }
};

int julianDay_c(int year, int month, int day);
double solarDeclination_c(int J);
double solarElevation_c(double latrad, double delta, double hrad);
double solarConstant_c(int J);
double daylength_c(double latrad, double slorad, double asprad, double delta);
double daylengthseconds_c(double latrad, double slorad, double asprad, double delta);
double RpotDay_c(double solarConstant, double latrad,  double slorad, double asprad, double delta);
double RDay_c(double solarConstant, double latrad, double elevation, double slorad, double asprad, double delta,
              double diffTemp, double diffTempMonth, double vpa, double precipitation);
double skyLongwaveRadiation_c(double Tair, double vpa, double c = 0);

void directDiffuseDay_c(DirectDiffuseDay_RESULT& res, 
                        double solarConstant, double latrad, double slorad, double asprad, double delta,
                        double R_s, bool clearday);

double outgoingLongwaveRadiation_c(double solarConstant, double latrad, double elevation,  double slorad,  double asprad, double delta,
                                   double vpa, double tmin, double tmax, double R_s);

double netRadiation_c(double solarConstant, double latrad,  double elevation, double slorad, double asprad, double delta,
                      double vpa, double tmin, double tmax, double R_s,
                      double alpha = 0.08);

Rcpp::DataFrame copyDirectDiffuseDayResult_c(const DirectDiffuseDay_RESULT& ddd);

#endif
