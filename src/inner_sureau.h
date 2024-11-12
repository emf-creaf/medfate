#include <Rcpp.h>

#ifndef INNER_COCHARD_H
#define INNER_COCHARD_H
#endif
using namespace Rcpp;

List initSureauNetwork(int c, NumericVector LAIphe,
                       DataFrame internalWater, 
                       DataFrame paramsAnatomy, DataFrame paramsTranspiration, DataFrame paramsWaterStorage,
                       NumericVector VCroot_kmax, NumericVector VGrhizo_kmax,
                       NumericVector PsiSoil, NumericVector VG_n, NumericVector VG_alpha,
                       List control, double sapFluidityDay = 1.0);

void innerSureau(List x, List input, List output, int n, double tstep, 
                 bool verbose = false);
