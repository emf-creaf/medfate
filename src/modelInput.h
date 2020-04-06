#include <Rcpp.h>

#ifndef MODELINPUT_H
#define MODELINPUT_H
#endif
using namespace Rcpp;

void checkSpeciesParameters(DataFrame SpParams, CharacterVector params);

List forest2spwbInput(List x, List soil, DataFrame SpParams, List control);
List forest2growthInput(List x, List soil, DataFrame SpParams, List control);

void resetInputs(List x, List soil, List from = R_NilValue, int day = NA_INTEGER);
