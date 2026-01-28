#include "communication_structures_c.h"
#include "modelInput_c.h"

#ifndef TRANSPIRATION_BASIC_C_H
#define TRANSPIRATION_BASIC_C_H

struct ParamsVolume{
  double leafpi0;
  double leafeps;
  double leafaf;
  double stempi0;
  double stemeps;
  double stem_c;
  double stem_d;
  double stemaf;
  double Vsapwood;
  double Vleaf;
  double LAI;
  double LAIlive;
}; 

double plantVol_c(double plantPsi, ParamsVolume pars);
double findNewPlantPsiConnected_c(double flowFromRoots, double plantPsi, double rootCrownPsi,
                                  ParamsVolume parsVol);
void transpirationBasic_c(BasicTranspirationOutput& transpOutput, ModelInput& x, 
                          const WeatherInputVector& meteovec,  const double elevation);
#endif
