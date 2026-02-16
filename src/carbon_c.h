#include <RcppArmadillo.h>
#include "medfate.h"
#include "modelInput_c.h"

#ifndef CARBON_C_H
#define CARBON_C_H


const double carbonMolarMass = 12.0107; //g*mol-1
const double glucoseMolarMass = 180.156; //g*mol-1
const double starchMolarMass = 162.1406; //g*mol-1
const double starchDensity = 1.5; //g·cm-3
const double leafCperDry = 0.3; //g C · g dry-1
const double rootCperDry = 0.4959; //g C · g dry-1


struct CarbonCompartments {
  std::vector<double> LeafStorageVolume;
  std::vector<double> SapwoodStorageVolume;
  std::vector<double> LeafStarchMaximumConcentration;
  std::vector<double> SapwoodStarchMaximumConcentration;
  std::vector<double> LeafStarchCapacity;
  std::vector<double> SapwoodStarchCapacity;
  std::vector<double> LeafStructuralBiomass;
  std::vector<double> TwigStructuralBiomass;
  std::vector<double> TwigLivingStructuralBiomass;
  std::vector<double> DeadLeafStructuralBiomass;
  std::vector<double> DeadTwigStructuralBiomass;
  std::vector<double> SapwoodStructuralBiomass;
  std::vector<double> SapwoodLivingStructuralBiomass;
  std::vector<double> HeartwoodStructuralBiomass;
  std::vector<double> AbovegroundWoodBiomass;
  std::vector<double> BelowgroundWoodBiomass;
  std::vector<double> FineRootBiomass;
  std::vector<double> StructuralBiomass;
  std::vector<double> DeadBiomass;
  std::vector<double> LabileBiomass;
  std::vector<double> TotalLivingBiomass;
  std::vector<double> TotalBiomass;
  
  CarbonCompartments(size_t numCohorts = 0) {
    LeafStorageVolume = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    SapwoodStorageVolume = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    LeafStarchMaximumConcentration = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    SapwoodStarchMaximumConcentration = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    LeafStarchCapacity = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    SapwoodStarchCapacity = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    LeafStructuralBiomass = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    TwigStructuralBiomass = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    TwigLivingStructuralBiomass = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    DeadLeafStructuralBiomass = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    DeadTwigStructuralBiomass = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    SapwoodStructuralBiomass = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    SapwoodLivingStructuralBiomass = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    HeartwoodStructuralBiomass = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    AbovegroundWoodBiomass = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    BelowgroundWoodBiomass = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    FineRootBiomass = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    StructuralBiomass = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    DeadBiomass = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    LabileBiomass = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    TotalLivingBiomass = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    TotalBiomass = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
  }
};


double osmoticWaterPotential_c(double sugarConc, double temp, double nonSugarConc);
double sugarConcentration_c(double osmoticWP, double temp, double nonSugarConc);
double turgor_c(double psi, double sugarConc, double temp, double nonSugarConc);
double relativeSapViscosity_c(double sugarConc, double temp);

double leafArea_c(double LAI, double N);

double leafStorageVolume_c(double LAI, double N, double SLA, double leafDensity);
double leafStructuralBiomass_c(double LAI, double N, double SLA);
double leafStarchCapacity_c(double LAI, double N, double SLA, double leafDensity);

double twigStructuralBiomass_c(double LAI, double N, double SLA, double r635);

double abovegroundSapwoodStructuralBiomass_c(double SA, double H, double woodDensity);
double belowgroundSapwoodStructuralBiomass_c(double SA, const std::vector<double>& L, const std::vector<double>& V, 
                                             double woodDensity);
double sapwoodStructuralBiomass_c(double SA, double H, const std::vector<double>& L, const std::vector<double>& V, 
                                  double woodDensity);
double abovegroundHeartwoodStructuralBiomass_c(double DBH, double SA, double H, double woodDensity);
double belowgroundHeartwoodStructuralBiomass_c(double DBH, double SA, const std::vector<double>& L, const std::vector<double>& V, 
                                               double woodDensity);
double heartwoodStructuralBiomass_c(double DBH, double SA, double H, const std::vector<double>& L, const std::vector<double>& V, 
                                    double woodDensity);

double sapwoodStorageVolume_c(double SA, double H, const std::vector<double>& L, const std::vector<double>& V, 
                            double woodDensity, double conduit2sapwood);
double sapwoodStructuralBiomass_c(double SA, double H, const std::vector<double>& L, const std::vector<double>& V, 
                                double woodDensity);
double sapwoodStarchCapacity_c(double SA, double H, const std::vector<double>& L, const std::vector<double>& V, 
                             double woodDensity, double conduit2sapwood);

double sugarStarchDynamicsLeaf_c(double sugarConc, double starchConc, double eqSugarConc);
double sugarStarchDynamicsStem_c(double sugarConc, double starchConc, double eqSugarConc);
double sugarStarchDynamicsRoot_c(double sugarConc, double starchConc, double eqSugarConc);


Rcpp::DataFrame copyCarbonCompartments_c(const CarbonCompartments& cc, ModelInput& x);
void fillCarbonCompartments_c(CarbonCompartments& cc, ModelInput& x, const std::string& biomassUnits);

#endif
