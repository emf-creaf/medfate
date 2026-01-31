#include "medfate.h"
#include <vector>

#ifndef ROOT_C_H
#define ROOT_C_H

void ldrRS_one_c(std::vector<double>& ldr, 
                 double Z50, double Z95, double Z100, const std::vector<double>&  d);

double fineRootRadius_c(double specificRootLength, double rootTissueDensity);

double specificRootSurfaceArea_c(double specificRootLength, double rootTissueDensity);

double fineRootSoilVolume_c(double fineRootBiomass, double specificRootLength, double rootLengthDensity );

#endif
