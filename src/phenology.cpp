#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export(".gdd")]]
NumericVector gdd(IntegerVector DOY, NumericVector Temp, double Tbase = 5.0, double cum = 0.0){
  int nDays = Temp.size();
  NumericVector GDD(nDays);
  for(int i=0;i<nDays;i++){
    if((Temp[i]-Tbase < 0.0) & (DOY[i]>180)) {
      cum = 0.0;
    } else if (DOY[i]<180){ //Only increase in the first part of the year
      if(Temp[i]-Tbase>0.0) cum = cum + (Temp[i]-Tbase);
    }
    GDD[i] = cum;
    if(DOY[i] >= 365) cum = 0.0;
  }
  return(GDD);
}

double leafDevelopmentStatus(double Sgdd, double gdd) {
  if(Sgdd>0.0) return(std::min(std::max(gdd/Sgdd,0.0),1.0));
  return 1.0;
}
// [[Rcpp::export("pheno_leafDevelopmentStatus")]]
NumericVector leafDevelopmentStatus(NumericVector Sgdd, double gdd) {
  NumericVector phe(Sgdd.size());
  for(int i=0;i<Sgdd.size();i++) phe[i] = leafDevelopmentStatus(Sgdd[i], gdd);
  return(phe);
}

// [[Rcpp::export("pheno_updateLeaves")]]
void updateLeaves(List x, double doy, double tmean, double wind, double Tbase = 5.0) {
  DataFrame paramsBase = Rcpp::as<Rcpp::DataFrame>(x["paramsBase"]);
  NumericVector Sgdd = paramsBase["Sgdd"];
  
  //Plant input
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  DataFrame cohorts = Rcpp::as<Rcpp::DataFrame>(x["cohorts"]);
  NumericVector SP = cohorts["SP"];
  NumericVector LAI_live = above["LAI_live"];
  NumericVector LAI_dead = above["LAI_dead"];
  NumericVector LAI_expanded = above["LAI_expanded"];
  int numCohorts = SP.size();
  
  List canopyParams = x["canopy"];
  double gddday = canopyParams["gdd"];
  if((tmean-Tbase < 0.0) & (doy>180)) {
    gddday = 0.0;
  } else if (doy<180){ //Only increase in the first part of the year
    if(tmean-Tbase>0.0) gddday = gddday + (tmean - Tbase);
  }
  canopyParams["gdd"] = gddday;
  
  //Update phenological status
  NumericVector phe = leafDevelopmentStatus(Sgdd, gddday);
  for(int j=0;j<numCohorts;j++) {
    LAI_dead[j] *= exp(-1.0*(wind/10.0)); //Decrease dead leaf area according to wind speed
    double LAI_exp_prev= LAI_expanded[j]; //Store previous value
    LAI_expanded[j] = LAI_live[j]*phe[j]; //Update expanded leaf area (will decrease if LAI_live decreases)
    LAI_dead[j] += std::max(0.0, LAI_exp_prev-LAI_expanded[j]);//Check increase dead leaf area if expanded leaf area has decreased
  }
}
