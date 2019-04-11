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
NumericVector leafDevelopmentStatus(NumericVector Sgdd, double gdd) {
  NumericVector phe(Sgdd.size());
  for(int i=0;i<Sgdd.size();i++) phe[i] = leafDevelopmentStatus(Sgdd[i], gdd);
  return(phe);
}

