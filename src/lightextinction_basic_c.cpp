#include "lightextinction_basic_c.h"
#include "forestutils_c.h"
#include <math.h>
#include <meteoland.h>

/***
 * FUNCTIONS FOR LIGHT EXTINCTION (BASIC)
 */
double availableLight_c(double h, 
                        const std::vector<double>& H, 
                        const std::vector<double>& LAI_expanded, 
                        const std::vector<double>& LAI_dead, 
                        const std::vector<double>& k, 
                        const std::vector<double>& CR) {
  double s= 0.0, p=0.0;
  int n = H.size();
  for(int j=0; j< n; j++) {
    double cbh = H[j]*(1.0-CR[j]);
    // p = (H[j]-h)/(H[j]*CR[j]);
    p = leafAreaProportion_c(h, H[j], cbh, H[j]);
    if(p<0.0) p = 0.0;
    else if(p>1.0) p=1.0;
    s = s + k[j]*p*(LAI_expanded[j]+LAI_dead[j]);
  }
  return(100*exp((-1)*s));
}

void parcohortC_c(std::vector<double>& PARcohort,
                  const std::vector<double>& H, 
                  const std::vector<double>& LAI_expanded, 
                  const std::vector<double>& LAI_dead, 
                  const std::vector<double>& k, 
                  const std::vector<double>& CR){
  int n = H.size();
  for(int i=0; i<n;i++) PARcohort[i] = availableLight_c( H[i]*(1.0-(1.0-CR[i])/2.0), H, LAI_expanded, LAI_dead, k, CR);
}


/**
 * Fraction of the incident SWR radiation in a layer that is absorbed
 */
void layerAbsorbedSWRFractionIncident_c(std::vector<double>& f,
                                        const arma::mat& LAIme, 
                                        const arma::mat& LAImd, 
                                        const std::vector<double>& kSWR) {
  int nlayer = LAIme.n_rows;
  int ncoh = LAIme.n_cols;
  double s;
  for(int l = 0;l<nlayer;l++) {
    s = 0.0;
    for(int c=0;c<ncoh; c++) s+=kSWR[c]*(LAIme(l,c)+LAImd(l,c));
    f[l] = 1.0 - exp(-1.0*s);
  }
}

/**
 * Fraction of the incident SWR radiation in a layer that is absorbed by live leaves of each cohort
 */
void cohortLayerAbsorbedSWRFractionIncident_c(arma::mat& fij,
                                              const std::vector<double>& fi, 
                                              const arma::mat& LAIme, 
                                              const arma::mat& LAImd, 
                                              const std::vector<double>& kSWR) {
  int nlayer = LAIme.n_rows;
  int ncoh = LAIme.n_cols;
  double s = 0.0;
  for(int l = 0;l<nlayer;l++) {
    s = 0.0;
    for(int c=0;c<ncoh; c++) s+=kSWR[c]*(LAIme(l,c)+LAImd(l,c));
    if(s>0.0) {
      for(int c=0;c<ncoh; c++) fij(l,c) = fi[l]*kSWR[c]*LAIme(l,c)/s; 
    }
  }
}

/**
 * Fraction of the SWR radiation that is absorbed by each cohort
 */
void cohortAbsorbedSWRFraction_c(std::vector<double>& SWRfraction, 
                                 const arma::mat& LAIme, 
                                 const arma::mat& LAImd, 
                                 const std::vector<double>& kSWR) {
  int ncoh = LAIme.n_cols;
  int nlayer = LAIme.n_rows;
  
  // TO DO: INCLUDE IN COMMUNICATION STRUCTURES
  std::vector<double> fi(nlayer,0.0); 
  arma::mat fij(nlayer, ncoh); 
  std::vector<double> rem(nlayer);
  
  layerAbsorbedSWRFractionIncident_c(fi, LAIme, LAImd, kSWR);
  cohortLayerAbsorbedSWRFractionIncident_c(fij, fi, LAIme, LAImd, kSWR);
  for(int i = 0;i<nlayer;i++) {
    rem[i] = 1.0;
    for(int h = (nlayer-1);h>i;h--) {
      rem[i] = rem[i]*(1.0-fi[h]);
    }
  }
  double s;
  for(int j=0;j<ncoh;j++) {
    s = 0.0;
    for(int i = 0;i<nlayer;i++) {
      s = s + fij(i,j)*rem[i];
    }
    SWRfraction[j] = s;
  }
}

