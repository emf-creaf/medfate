#include <Rcpp.h>

#ifndef PARAMUTILS_H
#define PARAMUTILS_H
#endif
using namespace Rcpp;

int findRowIndex(int sp, DataFrame SpParams);
void checkSpeciesParameters(DataFrame SpParams, CharacterVector params);

NumericVector speciesNumericParameter(IntegerVector SP, DataFrame SpParams, String parName);
CharacterVector speciesCharacterParameter(IntegerVector SP, DataFrame SpParams, String parName);
NumericVector fineFoliarRatioWithImputation(IntegerVector SP, DataFrame SpParams);
NumericVector specificLeafAreaWithImputation(IntegerVector SP, DataFrame SpParams);
NumericVector treeAllometricCoefficientWithImputation(IntegerVector SP, DataFrame SpParams, String parName);
NumericVector shrubAllometricCoefficientWithImputation(IntegerVector SP, DataFrame SpParams, String parName);

NumericVector kPARWithImputation(IntegerVector SP, DataFrame SpParams);
NumericVector gammaSWRWithImputation(IntegerVector SP, DataFrame SpParams);
NumericVector alphaSWRWithImputation(IntegerVector SP, DataFrame SpParams);
NumericVector gWithImputation(IntegerVector SP, DataFrame SpParams);