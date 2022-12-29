#include <Rcpp.h>

#ifndef PARAMUTILS_H
#define PARAMUTILS_H
#endif
using namespace Rcpp;

int findSpParamsRowByName(String spname, DataFrame SpParams);
int findSpParamsRowBySpIndex(int sp, DataFrame SpParams);

IntegerVector speciesIndex(CharacterVector species, DataFrame SpParams);

void checkSpeciesParameters(DataFrame SpParams, CharacterVector params);

NumericVector speciesNumericParameterFromIndex(IntegerVector SP, DataFrame SpParams, String parName);
NumericVector speciesNumericParameter(CharacterVector species, DataFrame SpParams, String parName);

CharacterVector speciesCharacterParameterFromIndex(IntegerVector SP, DataFrame SpParams, String parName);
CharacterVector speciesCharacterParameter(CharacterVector species, DataFrame SpParams, String parName);

NumericVector cohortNumericParameter(List x, DataFrame SpParams, String parName);
CharacterVector cohortCharacterParameter(List x, DataFrame SpParams, String parName);

NumericVector speciesNumericParameterWithImputation(IntegerVector SP, DataFrame SpParams, String parName, bool fillMissing = true);
NumericVector speciesNumericParameterWithImputation(CharacterVector species, DataFrame SpParams, String parName, bool fillMissing = true);

NumericVector cohortNumericParameterWithImputation(List x, DataFrame SpParams, String parName, bool fillMissing = true);
