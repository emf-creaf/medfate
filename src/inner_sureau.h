#include <Rcpp.h>
#include "struct_sureau.h"

#ifndef INNER_COCHARD_H
#define INNER_COCHARD_H
#endif
using namespace Rcpp;

void initSureauNetwork_inner(SureauNetwork &network, int c, NumericVector LAIphe,
                             DataFrame internalWater, 
                             DataFrame paramsAnatomy, DataFrame paramsTranspiration, DataFrame paramsWaterStorage,
                             NumericVector VCroot_kmax, NumericVector VGrhizo_kmax,
                             NumericVector PsiSoil, NumericVector VG_n, NumericVector VG_alpha, 
                             List control, double sapFluidityDay = 1.0);
void deleteSureauNetworkPointers(SureauNetwork &network);


void innerSureau(List x, SureauNetwork* networks, List input, List output, int n, double tstep, 
                 bool verbose = false);
