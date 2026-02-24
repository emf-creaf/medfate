#include <RcppArmadillo.h>
#include "carbon.h"
#include "carbon_c.h"
#include "phenology_c.h"
using namespace Rcpp;

// [[Rcpp::export(".gdd")]]
NumericVector gdd(IntegerVector DOY, NumericVector Temp, double Tbase = 5.0, double cum = 0.0){
  int nDays = Temp.size();
  NumericVector GDD(nDays);
  for(int i=0;i<nDays;i++){
    if((Temp[i]-Tbase < 0.0) && (DOY[i]>200)) {
      cum = 0.0;
    } else if (DOY[i]<200){ //Only increase in the first part of the year
      if(Temp[i]-Tbase>0.0) cum = cum + (Temp[i]-Tbase);
    }
    GDD[i] = cum;
    if(DOY[i] >= 365) cum = 0.0;
  }
  return(GDD);
}

//' Leaf phenology
//'
//' Function \code{pheno_leafDevelopmentStatus} returns the expanded status (0 to 1) of leaves according to the growth degree days required to start bud burst and leaf unfolding, as dictated by a simple ecodormancy (one-phase) model (Chuine et al. 2013). 
//' Function \code{pheno_leafSenescenceStatus} returns the 0/1 senescence status of leaves according to the one-phase senescence model of Delpierre et al. (2009) on the basis of photoperiod and temperature.
//' Function \code{pheno_updateLeaves} updates the status of expanded leaves and dead leaves of object \code{x} given the photoperiod, temperature and wind of a given day. It applies the development model for 1 < doy < 180 and the senescence model for 181 > doy > 365.
//' 
//' @param Sgdd Degree days required for leaf budburst (in Celsius).
//' @param gdd Cumulative degree days (in Celsius)
//' @param unfoldingDD Degree-days for complete leaf unfolding after budburst has occurred.
//' 
//' @return Function \code{pheno_leafDevelopmentStatus} returns a vector of values between 0 and 1, 
//' whereas function \code{pheno_leafSenescenceStatus} returns a vector of 0 (senescent) and 1 (expanded) values. 
//' The other two functions do not return any value (see note).
//' 
//' @note Functions \code{pheno_updatePhenology} and \code{pheno_updateLeaves} modify the input object \code{x}. The first modifies phenological state and the second modifies the leaf area accordingly.
//' 
//' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
//' 
//' @references
//' Chuine, I., De Cortazar-Atauri, I.G., Kramer, K., \enc{Hänninen}{Hanninen}, H., 2013. Plant development models. Phenology: An Integrative Environmental Science. Springer, pp. 275–293.
//' 
//' Delpierre N, Dufrêne E, Soudani K et al (2009) Modelling interannual and spatial variability of leaf senescence for three deciduous tree species in France. Agric For Meteorol 149:938–948. doi:10.1016/j.agrformet.2008.11.014
//' 
//' @seealso \code{\link{spwb}}, \code{\link{spwbInput}}
//' 
//' @name pheno_updateLeaves
//' @keywords internal
// [[Rcpp::export("pheno_leafDevelopmentStatus")]]
NumericVector leafDevelopmentStatus(NumericVector Sgdd, NumericVector gdd, double unfoldingDD = 300.0) {
  NumericVector phe(Sgdd.size());
  for(int i=0;i<Sgdd.size();i++) phe[i] = leafDevelopmentStatus_c(Sgdd[i], gdd[i], unfoldingDD);
  return(phe);
}

//' @param Ssen Threshold to start leaf senescence.
//' @param sen Cumulative senescence variable.
//' @rdname pheno_updateLeaves
//' @keywords internal
// [[Rcpp::export("pheno_leafSenescenceStatus")]]
LogicalVector leafSenescenceStatus(NumericVector Ssen, NumericVector sen) {
  LogicalVector phe(Ssen.size());
  for(int i=0;i<Ssen.size();i++) phe[i] = leafSenescenceStatus_c(Ssen[i], sen[i]);
  return(phe);
}

//' @param x An object of class \code{\link{spwbInput}}.
//' @param doy Day of the year.
//' @param photoperiod Day length (in hours).
//' @param tmean Average day temperature (in Celsius).
//' @rdname pheno_updateLeaves
//' @keywords internal
// [[Rcpp::export("pheno_updatePhenology")]]
void updatePhenology(List x, int doy, double photoperiod, double tmean) {
  ModelInput x_c(x);
  updatePhenology_c(x_c, doy, photoperiod, tmean);
  x_c.copyStateToList(x);
}

//' @param wind Average day wind speed (in m/s).
//' @param fromGrowthModel Boolean flag to indicate that routine is called from \code{\link{growth}} simulation function.
//' @rdname pheno_updateLeaves
//' 
//' @keywords internal
// [[Rcpp::export("pheno_updateLeaves")]]
void updateLeaves(List x, double wind, bool fromGrowthModel) {
  ModelInput x_c(x);
  updateLeaves_c(x_c, wind, fromGrowthModel);
  x_c.copyStateToList(x);
}
