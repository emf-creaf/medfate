#include "RcppArmadillo.h"
#include "forestutils_c.h"
#include "incgamma_c.h"

double leafAreaProportion_c(double z1, double z2, double zmin, double zmax) {
  double mu = (zmax+zmin)/2.0;
  double sd15 = (zmax-zmin)/2.0;
  double sd = sd15/1.5;
  z1 = std::max(z1, zmin);
  z2 = std::max(z2, zmin);
  z1 = std::min(z1, zmax);
  z2 = std::min(z2, zmax);
  double x1 = (z1-mu)/sd;
  double x2 = (z2-mu)/sd;
  double p1 = 0.5*(1.0+errorfunction_c(x1/sqrt(2.0), false, false));
  double p2 = 0.5*(1.0+errorfunction_c(x2/sqrt(2.0), false, false));
  double v = (p2-p1)/0.8663856; //truncated to -1.5 to 1.5
  return(v);
}

void updateLAIdistributionVectors_c(arma::mat& LAIdist, 
                                    const std::vector<double>& z, 
                                    const std::vector<double>& LAI, 
                                    const std::vector<double>& H, 
                                    const std::vector<double>& CR) {
  int ncanlayers = LAIdist.n_rows;
  int numCohorts = LAIdist.n_cols;
  for(int ci=0;ci<numCohorts;ci++) {
    double cbh = H[ci]*(1.0-CR[ci]);
    for(int hi=0;hi<ncanlayers;hi++) {
      if(z[hi]<= H[ci]) {
        LAIdist(hi,ci) = LAI[ci]*leafAreaProportion_c(z[hi],z[hi+1], cbh,H[ci]);
      } else {
        LAIdist(hi,ci) = 0.0;
      }
    }
  }
}