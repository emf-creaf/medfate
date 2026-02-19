#include <RcppArmadillo.h>
#include "medfate.h"

#ifndef WINDKATUL_C_H
#define WINDKATUL_C_H


struct CanopyTurbulenceModel_RESULT {
  std::vector<double> z1;
  std::vector<double> U1;
  std::vector<double> dU1;
  std::vector<double> epsilon1;
  std::vector<double> k1;
  std::vector<double> uw1;
  std::vector<double> Lmix1;
  CanopyTurbulenceModel_RESULT(size_t N) {
    z1 = std::vector<double>(N, medfate::NA_DOUBLE);
    U1 = std::vector<double>(N, medfate::NA_DOUBLE);
    dU1 = std::vector<double>(N, medfate::NA_DOUBLE);
    epsilon1 = std::vector<double>(N, medfate::NA_DOUBLE);
    k1 = std::vector<double>(N, medfate::NA_DOUBLE);
    uw1 = std::vector<double>(N, medfate::NA_DOUBLE);
    Lmix1 = std::vector<double>(N, medfate::NA_DOUBLE);
  }
};
Rcpp::DataFrame copyCanopyTurbulenceModelResult_c(const CanopyTurbulenceModel_RESULT& canopyTurbulenceModel);

struct CanopyTurbulence_RESULT {
  std::vector<double> zmid;
  std::vector<double> u;
  std::vector<double> du;
  std::vector<double> epsilon;
  std::vector<double> k;
  std::vector<double> uw;
  CanopyTurbulence_RESULT(size_t N) {
    zmid = std::vector<double>(N, medfate::NA_DOUBLE);
    u = std::vector<double>(N, medfate::NA_DOUBLE);
    du = std::vector<double>(N, medfate::NA_DOUBLE);
    epsilon = std::vector<double>(N, medfate::NA_DOUBLE);
    k = std::vector<double>(N, medfate::NA_DOUBLE);
    uw = std::vector<double>(N, medfate::NA_DOUBLE);
  }
};
Rcpp::DataFrame copyCanopyTurbulenceResult_c(const CanopyTurbulence_RESULT& canopyTurbulence);

void windCanopyTurbulenceModel_inner_c(CanopyTurbulenceModel_RESULT& comm, 
                                       const std::vector<double>& zm, 
                                       const std::vector<double>& Cx, 
                                       double hm, double d0, double z0,
                                       std::string model);
void windCanopyTurbulence_inner_c(CanopyTurbulence_RESULT& canopyTurbulence, CanopyTurbulenceModel_RESULT& canopyTurbulenceModel, 
                                  const std::vector<double>& zmid, 
                                  const std::vector<double>& LAD, 
                                  double canopyHeight,
                                  double u, double windMeasurementHeight, 
                                  std::string model);

#endif
