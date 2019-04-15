#include <Rcpp.h>

#ifndef PHOTOSYNTHESIS_H
#define PHOTOSYNTHESIS_H
#endif
using namespace Rcpp;

NumericVector leafphotosynthesis(double Q, double Catm, double Gc, double leaf_temp, double Vmax298, double Jmax298, bool verbose=false);

double VmaxTemp(double Vmax298, double leaf_temp);
double JmaxTemp(double Jmax298, double leaf_temp);

DataFrame leafPhotosynthesisFunction(NumericVector E, double Catm, double Patm, double Tair, double vpa, double u, 
                            double absRad, double Q, double Vmax298, double Jmax298, double Gwmin, double Gwmax, 
                            double leafWidth = 1.0, double refLeafArea = 1.0, bool verbose = false);

DataFrame sunshadePhotosynthesisFunction(NumericVector E, double Catm, double Patm, double Tair, double vpa, 
                                    double SLarea, double SHarea,
                                    double u, double absRadSL, double absRadSH,
                                    double QSL, double QSH, 
                                    double Vmax298SL, double Vmax298SH, 
                                    double Jmax298SL, double Jmax298SH, 
                                    double Gwmin, double Gwmax, double leafWidth = 1.0, bool verbose = false);

DataFrame multilayerPhotosynthesisFunction(NumericVector E, double Catm, double Patm, double Tair, double vpa, 
                                  NumericVector SLarea, NumericVector SHarea,
                                  NumericVector u, NumericVector absRadSL, NumericVector absRadSH,
                                  NumericVector QSL, NumericVector QSH, 
                                  NumericVector Vmax298, NumericVector Jmax298, 
                                  double Gwmin, double Gwmax, double leafWidth = 1.0, bool verbose = false);

