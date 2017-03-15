#include <Rcpp.h>

#ifndef PHOTOSYNTHESIS_H
#define PHOTOSYNTHESIS_H
#endif
using namespace Rcpp;

List profitMaximization(List supplyFunction, List photosynthesisFunction, double Gwmin, double Gwmax);
List profitMaximization2(List supplyFunction, List photosynthesisFunction, double Gwmin, double Gwmax, double kstemmax);
List profitMaximization3(List supplyFunction, List photosynthesisFunction, double Gwmin, double Gwmax, double kstemmax);

List leafPhotosynthesisFunction(List supplyFunction, double Catm, double Patm, double Tair, double vpa, double u, 
                            double absRad, double Q, double Vmax298, double Jmax298, bool verbose = false);

List canopyPhotosynthesisFunction(List supplyFunction, double Catm, double Patm, double Tair, double vpa, 
                                  NumericVector SLarea, NumericVector SHarea,
                                  NumericVector u, NumericVector absRadSL, NumericVector absRadSH,
                                  NumericVector QSL, NumericVector QSH, 
                                  NumericVector Vmax298, NumericVector Jmax298, 
                                  bool verbose = false);