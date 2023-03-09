#include <Rcpp.h>

#ifndef MODELINPUT_H
#define MODELINPUT_H
#endif
using namespace Rcpp;

List forest2spwbInput(List x, List soil, DataFrame SpParams, List control);
List forest2growthInput(List x, List soil, DataFrame SpParams, List control);

void resetInputs(List x);
