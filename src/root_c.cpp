#include <RcppArmadillo.h>
#include <vector>
#include "root_c.h"

void ldrRS_one_c(std::vector<double>& ldr, 
                 double Z50, double Z95, double Z100, const std::vector<double>&  d){
  int nlayer = d.size();
  double* dCum = new double[nlayer];
  double c = 2.94/log(Z50/Z95);
  ldr[0] = 1.0/(1.0+pow(d[0]/Z50,c));
  //Cumulate d
  dCum[0] = d[0];
  for(int i=1;i<nlayer;i++) dCum[i] = d[i]+dCum[i-1];
  //LDR equation
  for(int i=1;i<nlayer;i++){
    ldr[i] = 1.0/(1.0+pow(dCum[i]/Z50,c)) -1.0/(1.0+pow(dCum[i-1]/Z50,c));
  }
  //Truncate distribution if cumulative depth of the previous layer is larger than maximum rooting depth
  if(!std::isnan(Z100)) {
    for(int i=1;i<nlayer;i++){
      if(dCum[i-1]>Z100) ldr[i] = 0.0;
    }
  }
  //Rescale proportions so that they sum 1
  double Vtot = std::accumulate(ldr.begin(), ldr.end(), 0.0);
  for(int i=0;i<nlayer; i++) {
    ldr[i] = ldr[i]/Vtot;
  }
  delete[] dCum;
}


/**
 * Fine root radius in cm
 */
//' @rdname root
//' @keywords internal
// [[Rcpp::export("root_fineRootRadius")]]
double fineRootRadius_c(double specificRootLength, double rootTissueDensity) {
  return(sqrt(1.0/(M_PI*specificRootLength*rootTissueDensity)));
}

/**
 *  specificRootSurfaceArea (SRSA; cm2/g) as function of: 
 *    . specific root length (SRL; cm/g) e.g. 3870 cm/g
 *    . root tissue density (RTD; g/cm3) e.g. 0.165 g/cm3
 */
//' @rdname root
//' @keywords internal
// [[Rcpp::export("root_specificRootSurfaceArea")]]
double specificRootSurfaceArea_c(double specificRootLength, double rootTissueDensity) {
  return(2.0*sqrt(M_PI*specificRootLength/rootTissueDensity));
}

/**
 *   Estimates soil volume (m3) occupied with fine roots
 *    . fine root biomass (g dry)
 *    . specific root length (SRL; cm/g) e.g. 3870 cm/g
 *    . root length density (RLD; cm/cm3) e.g. 10 cm/cm3 = 0.1 mm/mm3
 */
//' @rdname root
//' @keywords internal
// [[Rcpp::export("root_fineRootSoilVolume")]]
double fineRootSoilVolume_c(double fineRootBiomass, double specificRootLength, double rootLengthDensity) {
  return(fineRootBiomass*(specificRootLength/rootLengthDensity)*1e-6);
}

