#include <Rcpp.h>

#ifndef INNER_SPERRY_H
#define INNER_SPERRY_H
#endif
using namespace Rcpp;

List profitMaximization2(List supplyFunction, int initialPos,
                         double Catm, double Patm, double Tair, double vpa, double u, 
                         double SWRabs, double LWRnet, double Q, double Vmax298, double Jmax298, 
                         double leafWidth, double refLeafArea,
                         double Gswmin, double Gswmax);

void innerSperry(List x, List input, List output, int n, double tstep, 
                 bool verbose = false, int stepFunctions = NA_INTEGER, bool modifyInput = true);