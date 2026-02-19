#include <RcppArmadillo.h>
#include "hydraulics.h"
#include "forestutils.h"
#include "tissuemoisture.h"
#include <meteoland.h>
using namespace Rcpp;








double layerLiveFuelMoisture(double minHeight, double maxHeight, NumericVector cohortFMC, NumericVector cohortLoading, NumericVector H, NumericVector CR) {
  double num = 0.0, den = 0.0, pin;
  int nCoh = cohortLoading.size();
  for(int i=0;i<nCoh; i++) {
    pin = (std::min(H[i], maxHeight)-std::max((1.0-CR[i])*H[i], minHeight))/(CR[i]*H[i]);
    if(pin<0.0) pin = 0.0;
    num +=(cohortFMC[i]*cohortLoading[i]*pin);
    den += (cohortLoading[i]*pin);
    // Rcout<<cohortFMC[i]<< " "<<H[i]<<" "<<maxHeight<<" "<< pBole[i]*H[i]<< " "<<minHeight<< ": "<<pin<<"\n";
  }
  if(den>0) return(num/den);
  return(NA_REAL);
}

double canopyLiveFuelMoisture(double canopyBaseHeight, double canopyTopHeight, NumericVector cohortFMC, NumericVector cohortLoading, NumericVector H, NumericVector CR) {
  return(layerLiveFuelMoisture(canopyBaseHeight, canopyTopHeight, cohortFMC, cohortLoading, H, CR));
}

double fuelbedLiveFuelMoisture(double fuelbedHeight, NumericVector cohortFMC, NumericVector cohortLoading, NumericVector H, NumericVector CR) {
  return(layerLiveFuelMoisture(0, fuelbedHeight, cohortFMC, cohortLoading, H, CR));
}
