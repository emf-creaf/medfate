#include <RcppArmadillo.h>
#include "modelInput_c.h"

#ifndef PHENOLOGY_C_H
#define PHENOLOGY_C_H

double leafDevelopmentStatus_c(double gdd, double Sgdd, double Ugdd);

void updatePhenology_c(ModelInput& x, int doy, double photoperiod, double tmean);
void updateLeaves_c(ModelInput& x, double wind, bool fromGrowthModel);
#endif