#include <Rcpp.h>

#ifndef PARAMUTILS_H
#define PARAMUTILS_H
#endif
using namespace Rcpp;

int findRowIndex(int sp, DataFrame SpParams);
void checkSpeciesParameters(DataFrame SpParams, CharacterVector params);

NumericVector speciesNumericParameter(IntegerVector SP, DataFrame SpParams, String parName);
CharacterVector speciesCharacterParameter(IntegerVector SP, DataFrame SpParams, String parName);
NumericVector cohortNumericParameter(List x, DataFrame SpParams, String parName);
CharacterVector cohortCharacterParameter(List x, DataFrame SpParams, String parName);

NumericVector speciesNumericParameterWithImputation(IntegerVector SP, DataFrame SpParams, String parName, bool fillMissing = true);
NumericVector cohortNumericParameterWithImputation(List x, DataFrame SpParams, String parName, bool fillMissing = true);
