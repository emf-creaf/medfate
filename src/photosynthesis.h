#include <Rcpp.h>

#ifndef PHOTOSYNTHESIS_H
#define PHOTOSYNTHESIS_H
#endif
using namespace Rcpp;

List profitMaximization(List supplyFunction, List photosynthesisFunction, double Gwmin, double Gwmax);
List profitMaximization2(List supplyFunction, List photosynthesisFunction, double Gwmin, double Gwmax, double kstemmax);
List profitMaximization3(List supplyFunction, List photosynthesisFunction, double Gwmin, double Gwmax, double kstemmax);
List photosynthesisFunction(List supplyFunction, double Catm, double Patm, double Tair, double vpa, double u, 
                            double absRad, double Q, double Vmax298, double Jmax298, bool verbose = false);