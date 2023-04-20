#include <Rcpp.h>

#ifndef FORESTUTILS_H
#define FORESTUTILS_H
#endif
using namespace Rcpp;

DataFrame forest2aboveground(List x, DataFrame SpParams, double gdd = NA_REAL, bool loading = false);
NumericMatrix forest2belowground(List x, List soil, DataFrame SpParams);

CharacterVector cohortIDs(List x, DataFrame SpParams, int treeOffset = 0, int shrubOffset = 0);


double leafAreaProportion(double z1, double z2, double zmin, double zmax);

NumericVector cohortHeight(List x, DataFrame SpParams);

NumericVector cohortDensity(List x, DataFrame SpParams);
NumericVector speciesDensity(List x, DataFrame SpParams);

NumericVector treeBasalArea(NumericVector N, NumericVector dbh);
NumericVector largerTreeBasalArea(NumericVector N, NumericVector dbh, double self_include_prop = 0.5);

NumericVector treeCohortBasalArea(List x);
NumericVector cohortBasalArea(List x, DataFrame SpParams);
NumericVector dbhClassBasalArea(List x, NumericVector DBHbreaks);
double forestBasalArea(List x);
double forestBasalAreaForMinDBH(List x, double minDBH);

double treeDensity(List x);
double minDBHDensity(List x, double minDBH);
NumericVector dbhClassDensity(List x, NumericVector DBHbreaks);

NumericVector treeCrownRatio(NumericVector N, NumericVector dbh, NumericVector H, 
                                NumericVector Acw, NumericVector Bcw,
                                NumericVector Acr, NumericVector B1cr, NumericVector B2cr, NumericVector B3cr,
                                NumericVector C1cr, NumericVector C2cr);

NumericVector cohortFuel(List x, DataFrame SpParams, double gdd = NA_REAL, bool includeDead = true);

NumericVector cohortCrownRatio(List x, DataFrame SpParams);
NumericVector cohortCrownBaseHeight(List x, DataFrame SpParams);
NumericVector cohortCrownLength(List x, DataFrame SpParams);

double herbFoliarBiomass(double herbCover, double herbHeight);
NumericVector cohortFoliarBiomass(List x, DataFrame SpParams, double gdd = NA_REAL);

NumericVector cohortEquilibriumLeafLitter(List x, DataFrame SpParams, double AET = 800);
NumericVector cohortEquilibriumSmallBranchLitter(List x, DataFrame SpParams, double smallBranchDecompositionRate = 0.81);

NumericVector cohortCover(List x, DataFrame SpParams);

double herbLAI(double herbCover, double herbHeight);
NumericVector cohortLAI(List x, DataFrame SpParams, double gdd = NA_REAL, bool bounded = true);

NumericMatrix LAIdistributionVectors(NumericVector z, NumericVector LAI, NumericVector H, NumericVector CR);
NumericVector LAIprofileVectors(NumericVector z, NumericVector LAI, NumericVector H, NumericVector CR);

NumericMatrix LAIdistribution(NumericVector z, List x, DataFrame SpParams, double gdd = NA_REAL, bool bounded = true);

IntegerVector uniqueSpp(IntegerVector sp);

double shrubCover(List x, double excludeMinHeight = 0.0);

