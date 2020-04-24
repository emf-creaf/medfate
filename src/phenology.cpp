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

double leafDevelopmentStatus(double Sgdd, double gdd, double unfoldingDD = 300.0) {
  double ds = 0.0;
  if(Sgdd>0.0) {
    if(gdd>Sgdd) ds = std::min(1.0, (gdd - Sgdd)/unfoldingDD);
  } else {
    ds = 1.0;
  }
  return(ds);
}
double leafSenescenceStatus(double Ssen, double sen) {
  if(sen>Ssen) return(0.0);
  return 1.0;
}

// [[Rcpp::export("pheno_leafDevelopmentStatus")]]
NumericVector leafDevelopmentStatus(NumericVector Sgdd, NumericVector gdd, double unfoldingDD = 300.0) {
  NumericVector phe(Sgdd.size());
  for(int i=0;i<Sgdd.size();i++) phe[i] = leafDevelopmentStatus(Sgdd[i], gdd[i], unfoldingDD);
  return(phe);
}

// [[Rcpp::export("pheno_leafSenescenceStatus")]]
NumericVector leafSenescenceStatus(NumericVector Ssen, NumericVector sen) {
  NumericVector phe(Ssen.size());
  for(int i=0;i<Ssen.size();i++) phe[i] = leafSenescenceStatus(Ssen[i], sen[i]);
  return(phe);
}

// [[Rcpp::export("pheno_updateLeaves")]]
void updateLeaves(List x, int doy, double photoperiod, double tmean, double wind) {
  
  List control = x["control"];
  double unfoldingDD = control["unfoldingDD"];
  
  DataFrame paramsPhenology = Rcpp::as<Rcpp::DataFrame>(x["paramsPhenology"]);
  CharacterVector phenoType = paramsPhenology["type"];
  NumericVector Sgdd = paramsPhenology["Sgdd"];
  NumericVector Tbgdd = paramsPhenology["Tbgdd"];
  NumericVector Ssen = paramsPhenology["Ssen"];
  NumericVector Psen = paramsPhenology["Psen"];
  NumericVector Tbsen = paramsPhenology["Tbsen"];
  
  //Plant input
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  DataFrame cohorts = Rcpp::as<Rcpp::DataFrame>(x["cohorts"]);
  NumericVector SP = cohorts["SP"];
  NumericVector LAI_live = above["LAI_live"];
  NumericVector LAI_dead = above["LAI_dead"];
  NumericVector LAI_expanded = above["LAI_expanded"];
  int numCohorts = SP.size();
  
  DataFrame internalPhenology =  Rcpp::as<Rcpp::DataFrame>(x["internalPhenology"]);
  NumericVector gdd = internalPhenology["gdd"];
  NumericVector sen = internalPhenology["sen"];
  for(int j=0;j<numCohorts;j++) {
    if(Sgdd[j]>0.0) {
      if(doy>180) {
        gdd[j] = 0.0;
      } else { //Only increase in the first part of the year
        if(tmean-Tbgdd[j]>0.0) gdd[j] = gdd[j] + (tmean - Tbgdd[j]);
      }
    }
    if(Ssen[j]>0.0) {
      // Rcout<< doy <<" "<< photoperiod<< " "<<Psen[j]<<"\n";
      if(photoperiod>Psen[j]) {
        sen[j] = 0.0;
      } else {
        double rsen = 0.0;
        if(tmean-Tbsen[j]<0.0) {
          rsen = pow(Tbsen[j]-tmean,2.0) * pow(photoperiod/Psen[j],2.0);
        }
        sen[j] = sen[j] + rsen;
      }
    }
  }

  //Update phenological status
  NumericVector phe;
  if(doy<180) phe = leafDevelopmentStatus(Sgdd, gdd, unfoldingDD);
  else phe = leafSenescenceStatus(Ssen, sen);
  for(int j=0;j<numCohorts;j++) {
    LAI_dead[j] *= exp(-1.0*(wind/10.0)); //Decrease dead leaf area according to wind speed
    if(phenoType[j] == "winter-deciduous") {
      double LAI_exp_prev= LAI_expanded[j]; //Store previous value
      LAI_expanded[j] = LAI_live[j]*phe[j]; //Update expanded leaf area (will decrease if LAI_live decreases)
      LAI_dead[j] += std::max(0.0, LAI_exp_prev-LAI_expanded[j]);//Check increase dead leaf area if expanded leaf area has decreased
    }
  }
}
