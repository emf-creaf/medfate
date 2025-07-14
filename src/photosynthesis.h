#include <Rcpp.h>
#include "struct_photosynthesis.h"

#ifndef PHOTOSYNTHESIS_H
#define PHOTOSYNTHESIS_H
#endif
using namespace Rcpp;

double gLeafBoundary(double u, double leafWidth, double gBound0 = 0.397);
double gCrown(double u, double gCrown0 = 0.150);

NumericVector leafphotosynthesis(double Q, double Catm, double Gc, double Tleaf, double Vmax298, double Jmax298, bool verbose=false);

void photosynthesisBaldocchi_inner(BaldocchiPhoto &photoOut,
                                   double Q, 
                                   double Catm, 
                                   double Tleaf, 
                                   double u,
                                   double Vmax298, 
                                   double Jmax298, 
                                   double leafWidth,
                                   double Gsw_AC_slope,
                                   double Gsw_AC_intercept);

double VmaxTemp(double Vmax298, double Tleaf);
double JmaxTemp(double Jmax298, double Tleaf);

DataFrame leafPhotosynthesisFunction(NumericVector E, NumericVector psiLeaf, double Catm, double Patm, double Tair, double vpa, double u, 
                                     double absRad, double Q, double Vmax298, double Jmax298, 
                                     double leafWidth = 1.0, double refLeafArea = 1.0, bool verbose = false);

NumericVector leafPhotosynthesisOneFunction2(double E, double psiLeaf, double Catm, double Patm, double Tair, double vpa, double u, 
                                             double SWRabs, double LWRnet, double Q, double Vmax298, double Jmax298, 
                                             double leafWidth = 1.0, double refLeafArea = 1.0, bool verbose = false);

DataFrame leafPhotosynthesisFunction2(NumericVector E, NumericVector psiLeaf, double Catm, double Patm, double Tair, double vpa, double u, 
                                     double SWRabs, double LWRnet, double Q, double Vmax298, double Jmax298, 
                                     double leafWidth = 1.0, double refLeafArea = 1.0, bool verbose = false);

DataFrame sunshadePhotosynthesisFunction(NumericVector E, NumericVector psiLeaf, double Catm, double Patm, double Tair, double vpa, 
                                    double SLarea, double SHarea,
                                    double u, double absRadSL, double absRadSH,
                                    double QSL, double QSH, 
                                    double Vmax298SL, double Vmax298SH, 
                                    double Jmax298SL, double Jmax298SH, 
                                    double leafWidth = 1.0, bool verbose = false);

DataFrame multilayerPhotosynthesisFunction(NumericVector E, NumericVector psiLeaf, 
                                           double Catm, double Patm, double Tair, double vpa, 
                                  NumericVector SLarea, NumericVector SHarea,
                                  NumericVector u, NumericVector absRadSL, NumericVector absRadSH,
                                  NumericVector QSL, NumericVector QSH, 
                                  NumericVector Vmax298, NumericVector Jmax298, 
                                  double leafWidth = 1.0, bool verbose = false);

