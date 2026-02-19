#include <RcppArmadillo.h>
#include "modelInput_c.h"

#ifndef PHENOLOGY_C_H
#define PHENOLOGY_C_H

double leafDevelopmentStatus_c(double Sgdd, double gdd, double unfoldingDD = 300.0);
bool leafSenescenceStatus_c(double Ssen, double sen);

void updatePhenology_c(ModelInput& x, int doy, double photoperiod, double tmean);
void updateLeaves_c(ModelInput& x, double wind, bool fromGrowthModel);
#endif