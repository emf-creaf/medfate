#include <Rcpp.h>

#ifndef INNER_COCHARD_H
#define INNER_COCHARD_H
#endif
using namespace Rcpp;

List initCochardNetwork(int c,
                       DataFrame internalWater, DataFrame paramsTranspiration, DataFrame paramsWaterStorage,
                       NumericVector VCroot_kmax, NumericVector VGrhizo_kmax,
                       double sapFluidityDay = 1.0);

void innerCochard(List x, List input, List output, int n, double tstep, 
                 bool verbose = false, bool modifyInput = true);
