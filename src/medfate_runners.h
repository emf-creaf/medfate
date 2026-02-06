#include "RcppArmadillo.h"
#include "medfate.h"
#include "modelInput_c.h"
#include "spwb_day_c.h"


#ifndef MEDFATE_RUNNERS_H
#define MEDFATE_RUNNERS_H

class SPWBrunner {
private:
  ModelInput x;
  double latitude, elevation, slope, aspect;
  SPWBCommunicationStructures SPWBcomm;
public:
  SPWBrunner(Rcpp::List x_list, 
                double latitude, double elevation, double slope, double aspect);
  void run_day(Rcpp::CharacterVector date, Rcpp::NumericVector meteovec, 
               double runon = 0, Rcpp::Nullable<Rcpp::NumericVector> lateralFlows = R_NilValue, double waterTableDepth = NA_REAL);
  Rcpp::List get_output();
};

#endif