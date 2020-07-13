#include <Rcpp.h>

#ifndef FORESTUTILS_H
#define FORESTUTILS_H
#endif
using namespace Rcpp;

int findRowIndex(int sp, DataFrame SpParams);

DataFrame forest2aboveground(List x, DataFrame SpParams, double gdd = NA_REAL, String mode = "MED");
NumericMatrix forest2belowground(List x, List soil);

NumericVector cohortNumericParameter(List x, DataFrame SpParams, String parName);
NumericVector cohortNumericParameter(IntegerVector SP, DataFrame SpParams, String parName);
CharacterVector cohortCharacterParameter(List x, DataFrame SpParams, String parName);
CharacterVector cohortCharacterParameter(IntegerVector SP, DataFrame SpParams, String parName);

double leafAreaProportion(double z1, double z2, double zmin, double zmax);

CharacterVector cohortIDs(List x);

NumericVector cohortHeight(List x);

NumericVector cohortDensity(List x, DataFrame SpParams, String mode = "MED");
NumericVector speciesDensity(List x, DataFrame SpParams, String mode = "MED");

NumericVector treeBasalArea(NumericVector N, NumericVector dbh);
NumericVector treeCohortBasalArea(List x);
NumericVector cohortBasalArea(List x);
NumericVector dbhClassBasalArea(List x, NumericVector DBHbreaks);
double forestBasalArea(List x);
double forestBasalAreaForMinDBH(List x, double minDBH);

double treeDensity(List x);
double minDBHDensity(List x, double minDBH);
NumericVector dbhClassDensity(List x, NumericVector DBHbreaks);

NumericVector treeCrownRatioMED(NumericVector N, NumericVector dbh, NumericVector H, 
                                NumericVector Acw, NumericVector Bcw,
                                NumericVector Acr, NumericVector B1cr, NumericVector B2cr, NumericVector B3cr,
                                NumericVector C1cr, NumericVector C2cr);

NumericVector treeFuel(IntegerVector SP, NumericVector N, NumericVector dbh, DataFrame SpParams, double gdd = NA_REAL, bool includeDead = true);
NumericVector shrubFuelMED(IntegerVector SP, NumericVector Cover, NumericVector H, DataFrame SpParams, double gdd = NA_REAL, bool includeDead = true);
NumericVector shrubFuelUS(IntegerVector SP, NumericVector Cover, NumericVector H, DataFrame SpParams, double gdd = NA_REAL, bool includeDead = true);

NumericVector cohortFuel(List x, DataFrame SpParams, double gdd = NA_REAL, bool includeDead = true, String mode = "MED");

NumericVector cohortCrownRatio(List x, DataFrame SpParams, String mode = "MED");
NumericVector cohortCrownBaseHeight(List x, DataFrame SpParams, String mode = "MED");
NumericVector cohortCrownLength(List x, DataFrame SpParams, String mode = "MED");

NumericVector cohortFoliarBiomass(List x, DataFrame SpParams, double gdd = NA_REAL, String mode = "MED");

NumericVector cohortEquilibriumLeafLitter(List x, DataFrame SpParams, double AET = 800);
NumericVector cohortEquilibriumSmallBranchLitter(List x, DataFrame SpParams, double smallBranchDecompositionRate = 0.81);

NumericVector treeLAI(IntegerVector SP, NumericVector N, NumericVector dbh, DataFrame SpParams, NumericVector pEmb=NumericVector(0), double gdd = NA_REAL);
NumericVector shrubLAI(IntegerVector SP, NumericVector Cover, NumericVector H, DataFrame SpParams, double gdd = NA_REAL);
NumericVector cohortLAI(List x, DataFrame SpParams, double gdd = NA_REAL, String mode = "MED");
NumericMatrix LAIdistributionVectors(NumericVector z, NumericVector LAI, NumericVector H, NumericVector CR);
NumericMatrix LAIdistribution(NumericVector z, List x, DataFrame SpParams, double gdd = NA_REAL, String mode = "MED");

double shrubCover(List x, double excludeMinHeight = 0.0);

void deleteTreeCohort(List x, int treeCohort);
void deleteShrubCohort(List x, int shrubCohort);

int minDBHTreeCohort(List x, double excludeMin = 0.0);
