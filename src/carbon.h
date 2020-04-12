#include <Rcpp.h>

#ifndef CARBON_H
#define CARBON_H
#endif
using namespace Rcpp;

const double carbonMolarMass = 12.0107; //g*mol-1
const double glucoseMolarMass = 180.156; //g*mol-1
const double starchMolarMass = 162.1406; //g*mol-1
const double starchDensity = 1.5; //g·cm-3
const double leafCperDry = 0.3; //g C · g dry-1
const double rootCperDry = 0.4959; //g C · g dry-1


double osmoticWaterPotential(double conc, double temp, double nonSugarConc = 0.4);
double sugarConcentration(double osmoticWP, double temp, double nonSugarConc = 0.4);
double turgor(double psi, double conc, double temp);
double relativeSapViscosity(double conc, double temp);

double leafArea(double LAI, double N);

double leafStorageVolume(double LAI, double N, double SLA, double leafDensity);
double leafStructuralBiomass(double LAI, double N, double SLA);
double leafStarchCapacity(double LAI, double N, double SLA, double leafDensity);

double sapwoodStorageVolume(double SA, double H, double Z, double woodDensity, double vessel2sapwood);
double sapwoodStructuralBiomass(double SA, double H, double Z, double woodDensity);
double sapwoodStructuralLivingBiomass(double SA, double H, double Z, double woodDensity, double vessel2sapwood);
double sapwoodStarchCapacity(double SA, double H, double Z, double woodDensity, double vessel2sapwood);

double sugarStarchDynamicsLeaf(double sugarConc, double starchConc, double eqSugarConc);
double sugarStarchDynamicsStem(double sugarConc, double starchConc, double eqSugarConc);
double sugarStarchDynamicsRoot(double sugarConc, double starchConc, double eqSugarConc);

// NumericVector carbonCompartments(double SA, double LAI, double H, double Z, double N, double SLA, double WoodDensity, double WoodC);