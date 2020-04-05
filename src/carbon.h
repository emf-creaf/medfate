#include <Rcpp.h>

#ifndef CARBON_H
#define CARBON_H
#endif
using namespace Rcpp;

const double glucoseMolarWeight = 180.156; //g*mol-1
const double starchMolarWeight = 162.1406; //g*mol-1
const double starchDensity = 1.5; //g·cm-3
const double leafCperDry = 0.3; //g C · g dry-1
const double rootCperDry = 0.4959; //g C · g dry-1

double leafStarchCapacity(double LAI, double N, double SLA, double leafDensity);
double leafStorageVolume(double LAI, double N, double SLA, double leafDensity);
double leafCstructural(double LAI, double N, double SLA);

double sapwoodStorageVolume(double SA, double H, double Z, double woodDensity, double vessel2sapwood);
double sapwoodCstructural(double SA, double H, double Z, double woodDensity, double woodCperDry);
double sapwoodStarchCapacity(double SA, double H, double Z, double woodDensity, double vessel2sapwood);

NumericVector carbonCompartments(double SA, double LAI, double H, double Z, double N, double SLA, double WoodDensity, double WoodC);