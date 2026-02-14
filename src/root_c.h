#include "medfate.h"
#include <vector>

#ifndef ROOT_C_H
#define ROOT_C_H

struct RadialAxialLengths {
  std::vector<double> radial;
  std::vector<double> axial;
  
  RadialAxialLengths(size_t nlayers) {
    radial = std::vector<double>(nlayers, medfate::NA_DOUBLE);
    axial = std::vector<double>(nlayers, medfate::NA_DOUBLE);
  }
};

void ldrRS_one_c(std::vector<double>& ldr, 
                 double Z50, double Z95, double Z100, const std::vector<double>&  d);

double fineRootRadius_c(double specificRootLength, double rootTissueDensity);

double specificRootSurfaceArea_c(double specificRootLength, double rootTissueDensity);

double fineRootSoilVolume_c(double fineRootBiomass, double specificRootLength, double rootLengthDensity );

double coarseRootSoilVolumeFromConductance_c(double Kmax_rootxylem, double VCroot_kmax, double Al2As,
                                             const std::vector<double>& v, const std::vector<double>& d, const std::vector<double>& rfc);

std::vector<double> coarseRootLengthsFromVolume_c(double VolInd, const std::vector<double>& v, const std::vector<double>& d, const std::vector<double>& rfc);

RadialAxialLengths coarseRootRadialAxialLengths_c(const std::vector<double>& v, const std::vector<double>& d, double depthWidthRatio);

std::vector<double> coarseRootLengths_c(const std::vector<double>& v, const std::vector<double>& d, double depthWidthRatio);

double coarseRootSoilVolume_c(const std::vector<double>& v, const std::vector<double>& d, double depthWidthRatio);

double fineRootLengthPerArea_c(double Ksoil, double krhizo, double lai,
                               double radius, double rootLengthDensity);
  
double fineRootBiomassPerIndividual_c(const std::vector<double>& Ksoil, const std::vector<double>& krhizo, double lai, double N,
                                      double specificRootLength, double rootTissueDensity,  
                                      double rootLengthDensity);

std::vector<double> rhizosphereMaximumConductance_c(const std::vector<double>& Ksoil, const std::vector<double>& fineRootBiomass, double lai, double N,
                                                    double specificRootLength, double rootTissueDensity,  
                                                    double rootLengthDensity);
#endif
