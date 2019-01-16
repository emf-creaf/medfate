#include <Rcpp.h>

#ifndef SPWB_H
#define SPWB_H
#endif
using namespace Rcpp;

void checkspwbInput(List x, List soil, String transpirationMode);
void resetInputs(List x, List soil, List from = R_NilValue, int day = NA_INTEGER);

List spwbDay1(List x, List soil, double tday, double pet, double prec, double er, double runon=0.0, 
              double rad = NA_REAL, double elevation = NA_REAL, bool verbose=false);
List spwbDay2(List x, List soil, double tmin, double tmax, double rhmin, double rhmax, double rad, double wind, 
             double latitude, double elevation, double solarConstant, double delta, 
             double prec, double er, double runon=0.0, bool verbose = false);