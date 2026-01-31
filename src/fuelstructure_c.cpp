#include "medfate.h"
#include "forestutils_c.h"
#include "fuelstructure_c.h"

/**
 * Returns the proportion of the crown (Hbc - H) that lies within given interval (zLow-zHigh)
 */
double crownProportionInLayer_c(double zLow, double zHigh, double H, double Hbc) {
  return(leafAreaProportion_c(zLow, zHigh, Hbc, H));
}


/**
 * Returns the fuel loading (kg/m2) of the crown (Hbc - H) that lies within given interval (zLow-zHigh)
 */
double crownFuelInLayer_c(double zLow, double zHigh, double fb, double H, double Hbc) {
  return(fb*crownProportionInLayer_c(zLow, zHigh, H, Hbc));
}

/**
 * Returns the crown length (cm) of the crown (Hbc - H) that lies within given interval (zLow-zHigh)
 */
double crownLengthInLayer_c(double zLow, double zHigh, double cl, double H, double Hbc) {
  return(cl*crownProportionInLayer_c(zLow, zHigh, H, Hbc));
}

double layerFuelAverageParameter_c(double minHeight, double maxHeight, 
                                   const std::vector<double>& cohortParameter, 
                                   const std::vector<double>& cohortLoading, 
                                   const std::vector<double>& H, 
                                   const std::vector<double>& CR) {
  double num = 0.0, den = 0.0, cfl = 0.0;
  int nCoh = cohortLoading.size();
  for(int i=0;i<nCoh; i++) {
    cfl = crownFuelInLayer_c(minHeight, maxHeight, cohortLoading[i], H[i], H[i]*(1.0 - CR[i]));
    if(!std::isnan(cohortParameter[i])) {
      num +=(cohortParameter[i]*cfl);
      den += cfl;
    }
  }
  if(den>0) return(num/den);
  return(medfate::NA_DOUBLE);
}