#include <RcppArmadillo.h>

#ifndef CARBON_H
#define CARBON_H
using namespace Rcpp;

double sapwoodStructuralBiomass(double SA, double H, NumericVector L, NumericVector V, 
                                double woodDensity);
double sapwoodStorageVolume(double SA, double H, NumericVector L, NumericVector V, 
                            double woodDensity, double conduit2sapwood);
double sapwoodStarchCapacity(double SA, double H, NumericVector L, NumericVector V, 
                             double woodDensity, double conduit2sapwood);

DataFrame carbonCompartments(List x, String biomassUnits);
void fillCarbonCompartments(DataFrame cc, List x, String biomassUnits);

#endif
