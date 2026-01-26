#include <RcppArmadillo.h>
#include "medfate.h"

#ifndef CARBON_H
#define CARBON_H
using namespace Rcpp;


struct CarbonCompartments {
  std::vector<double> LeafStorageVolume;
  std::vector<double> SapwoodStorageVolume;
  std::vector<double> LeafStarchMaximumConcentration;
  std::vector<double> SapwoodStarchMaximumConcentration;
  std::vector<double> LeafStarchCapacity;
  std::vector<double> SapwoodStarchCapacity;
  std::vector<double> LeafStructuralBiomass;
  std::vector<double> TwigStructuralBiomass;
  std::vector<double> SapwoodStructuralBiomass;
  std::vector<double> SapwoodLivingStructuralBiomass;
  std::vector<double> HeartwoodStructuralBiomass;
  std::vector<double> AbovegroundWoodBiomass;
  std::vector<double> BelowgroundWoodBiomass;
  std::vector<double> FineRootBiomass;
  std::vector<double> StructuralBiomass;
  std::vector<double> LabileBiomass;
  std::vector<double> TotalLivingBiomass;
  std::vector<double> TotalBiomass;
  
  CarbonCompartments(size_t numCohorts = 0) {
    auto na = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    LeafStorageVolume = SapwoodStorageVolume = na;
    LeafStarchMaximumConcentration = SapwoodStarchMaximumConcentration = na;
    LeafStarchCapacity = SapwoodStarchCapacity = na;
    LeafStructuralBiomass = TwigStructuralBiomass = na;
    SapwoodStructuralBiomass = SapwoodLivingStructuralBiomass = na;
    HeartwoodStructuralBiomass = AbovegroundWoodBiomass = na;
    BelowgroundWoodBiomass = FineRootBiomass = na;
    StructuralBiomass = LabileBiomass = na;
    TotalLivingBiomass = TotalBiomass = na;
  }
};

const double carbonMolarMass = 12.0107; //g*mol-1
const double glucoseMolarMass = 180.156; //g*mol-1
const double starchMolarMass = 162.1406; //g*mol-1
const double starchDensity = 1.5; //g·cm-3
const double leafCperDry = 0.3; //g C · g dry-1
const double rootCperDry = 0.4959; //g C · g dry-1


double osmoticWaterPotential(double sugarConc, double temp, double nonSugarConc);
double sugarConcentration(double osmoticWP, double temp, double nonSugarConc);
double turgor(double psi, double sugarConc, double temp, double nonSugarConc);
double relativeSapViscosity(double sugarConc, double temp);

double leafArea(double LAI, double N);

double leafStorageVolume(double LAI, double N, double SLA, double leafDensity);
double leafStructuralBiomass(double LAI, double N, double SLA);
double leafStarchCapacity(double LAI, double N, double SLA, double leafDensity);

double sapwoodStorageVolume(double SA, double H, NumericVector L, NumericVector V, 
                            double woodDensity, double conduit2sapwood);
double sapwoodStructuralBiomass(double SA, double H, NumericVector L, NumericVector V, 
                                double woodDensity);
double sapwoodStarchCapacity(double SA, double H, NumericVector L, NumericVector V, 
                             double woodDensity, double conduit2sapwood);

double sugarStarchDynamicsLeaf(double sugarConc, double starchConc, double eqSugarConc);
double sugarStarchDynamicsStem(double sugarConc, double starchConc, double eqSugarConc);
double sugarStarchDynamicsRoot(double sugarConc, double starchConc, double eqSugarConc);

DataFrame carbonCompartments(List x, String biomassUnits = "g_m2");
void fillCarbonCompartments(DataFrame cc, List x, String biomassUnits);

#endif
