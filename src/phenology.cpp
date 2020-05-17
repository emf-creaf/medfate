#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export(".gdd")]]
NumericVector gdd(IntegerVector DOY, NumericVector Temp, double Tbase = 5.0, double cum = 0.0){
  int nDays = Temp.size();
  NumericVector GDD(nDays);
  for(int i=0;i<nDays;i++){
    if((Temp[i]-Tbase < 0.0) & (DOY[i]>200)) {
      cum = 0.0;
    } else if (DOY[i]<200){ //Only increase in the first part of the year
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
bool leafSenescenceStatus(double Ssen, double sen) {
  if(sen>Ssen) return(true);
  return false;
}

// [[Rcpp::export("pheno_leafDevelopmentStatus")]]
NumericVector leafDevelopmentStatus(NumericVector Sgdd, NumericVector gdd, double unfoldingDD = 300.0) {
  NumericVector phe(Sgdd.size());
  for(int i=0;i<Sgdd.size();i++) phe[i] = leafDevelopmentStatus(Sgdd[i], gdd[i], unfoldingDD);
  return(phe);
}

// [[Rcpp::export("pheno_leafSenescenceStatus")]]
LogicalVector leafSenescenceStatus(NumericVector Ssen, NumericVector sen) {
  LogicalVector phe(Ssen.size());
  for(int i=0;i<Ssen.size();i++) phe[i] = leafSenescenceStatus(Ssen[i], sen[i]);
  return(phe);
}

// [[Rcpp::export("pheno_updatePhenology")]]
void updatePhenology(List x, int doy, double photoperiod, double tmean) {
  List control = x["control"];
  double unfoldingDD = control["unfoldingDD"];
  
  DataFrame paramsPhenology = Rcpp::as<Rcpp::DataFrame>(x["paramsPhenology"]);
  CharacterVector phenoType = paramsPhenology["PhenologyType"];
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
  NumericVector phi = internalPhenology["phi"];
  LogicalVector budFormation = internalPhenology["budFormation"];
  LogicalVector leafUnfolding = internalPhenology["leafUnfolding"];
  LogicalVector leafSenescence = internalPhenology["leafSenescence"];
  LogicalVector leafDormancy = internalPhenology["leafDormancy"];
  
  for(int j=0;j<numCohorts;j++) {
    if(phenoType[j] == "winter-deciduous" || phenoType[j] == "winter-semideciduous") {
      if(doy>200) {
        double rsen = 0.0;
        phi[j] = 0.0;
        gdd[j] = 0.0;
        leafUnfolding[j] = false;
        if(photoperiod>Psen[j]) {
          sen[j] = 0.0;
          leafSenescence[j] = false;
          budFormation[j] = false;
          leafDormancy[j] = false;
        } else if(!leafDormancy[j]) {
          if(tmean-Tbsen[j]<0.0) {
            rsen = pow(Tbsen[j]-tmean,2.0);
            // rsen = pow(Tbsen[j]-tmean,2.0) * pow(photoperiod/Psen[j],2.0);
          }
          sen[j] = sen[j] + rsen;
          leafSenescence[j] = leafSenescenceStatus(Ssen[j],sen[j]);
          budFormation[j] = leafSenescence[j];
          leafDormancy[j] = leafSenescence[j];
        }
        // Rcout << doy<< " "<< photoperiod<<" "<< rsen <<" "<< sen[j]<<" "<<  leafSenescence[j] << "\n";
      } else if (doy<=200) { //Only increase in the first part of the year
        sen[j] = 0.0;
        budFormation[j] = false;
        leafSenescence[j] = false;
        if(tmean-Tbgdd[j]>0.0) gdd[j] = gdd[j] + (tmean - Tbgdd[j]);
        phi[j] = leafDevelopmentStatus(Sgdd[j], gdd[j], unfoldingDD);
        leafUnfolding[j] = (phi[j]>0.0);
        leafDormancy[j] = (phi[j]==0.0);
        // Rcout << doy<< " "<< photoperiod<<" "<< gdd[j]<<" "<<  leafUnfolding[j] << "\n";
      }
    }
    else if(phenoType[j] == "oneflush-evergreen") {
      if(doy>200) {
        gdd[j] = 0.0;
        leafUnfolding[j] = false;
        leafSenescence[j] = false;
        if(photoperiod>Psen[j]) {
          sen[j] = 0.0;
          budFormation[j] = true;
          leafDormancy[j] = false;
        } else if (!leafDormancy[j]){
          double rsen = 0.0;
          if(tmean-Tbsen[j]<0.0) {
            rsen = pow(Tbsen[j]-tmean,2.0);
            // rsen = pow(Tbsen[j]-tmean,2.0) * pow(photoperiod/Psen[j],2.0);
          }
          sen[j] = sen[j] + rsen;
          leafDormancy[j] = leafSenescenceStatus(Ssen[j],sen[j]);
          budFormation[j] = !leafDormancy[j];
        }
      } else if (doy<=200) { //Only increase in the first part of the year
        sen[j] = 0.0;
        budFormation[j] = false;
        if(!leafUnfolding[j]) {
          if(tmean-Tbgdd[j]>0.0) gdd[j] = gdd[j] + (tmean - Tbgdd[j]);
          double ph = leafDevelopmentStatus(Sgdd[j], gdd[j],unfoldingDD);
          leafSenescence[j] = (ph>0.0);
          leafUnfolding[j] = (ph>0.0);
          leafDormancy[j] = (ph==0.0);
        }
        // Rcout<<j<< " phi: "<< ph<<"\n";
      }
    }
    else if(phenoType[j] == "progressive-evergreen") {
      leafSenescence[j] = true;
      leafUnfolding[j] = true;
      budFormation[j] = true;
      leafDormancy[j] = false;
    }
  }
}

// [[Rcpp::export("pheno_updateLeaves")]]
void updateLeaves(List x, double wind, bool fromGrowthModel) {
  List control = x["control"];
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  NumericVector LAI_live = above["LAI_live"];
  NumericVector LAI_dead = above["LAI_dead"];
  NumericVector LAI_expanded = above["LAI_expanded"];
  int numCohorts = LAI_live.size();

  DataFrame paramsPhenology = Rcpp::as<Rcpp::DataFrame>(x["paramsPhenology"]);
  CharacterVector phenoType = paramsPhenology["PhenologyType"];
  NumericVector leafDuration = paramsPhenology["LeafDuration"];
  
  DataFrame internalPhenology =  Rcpp::as<Rcpp::DataFrame>(x["internalPhenology"]);
  NumericVector phi = internalPhenology["phi"];
  LogicalVector budFormation = internalPhenology["budFormation"];
  LogicalVector leafUnfolding = internalPhenology["leafUnfolding"];
  LogicalVector leafSenescence = internalPhenology["leafSenescence"];
  LogicalVector leafDormancy = internalPhenology["leafDormancy"];
  
  for(int j=0;j<numCohorts;j++) {
    bool leafFall = true;
    if(phenoType[j] == "winter-semideciduous") leafFall = leafUnfolding[j];
    if(leafFall) LAI_dead[j] *= exp(-1.0*(wind/10.0)); //Decrease dead leaf area according to wind speed
    if(phenoType[j] == "winter-deciduous" || phenoType[j] == "winter-semideciduous") {
      if(leafSenescence[j]) {
        double LAI_exp_prev= LAI_expanded[j]; //Store previous value
        LAI_expanded[j] = 0.0; //Update expanded leaf area (will decrease if LAI_live decreases)
        if(fromGrowthModel) LAI_live[j] = LAI_expanded[j];
        LAI_dead[j] += LAI_exp_prev;//Check increase dead leaf area if expanded leaf area has decreased
        leafSenescence[j] = false;
        budFormation[j] = false;
        leafDormancy[j] = true;
      } 
      else if(leafDormancy[j]) {
        if(phenoType[j] == "winter-semideciduous") LAI_dead[j] += LAI_expanded[j];
        LAI_expanded[j] = 0.0;
        if(fromGrowthModel) LAI_live[j] = LAI_expanded[j];
      }
      else if(!fromGrowthModel && leafUnfolding[j]) {
        LAI_expanded[j] = LAI_live[j]*phi[j]; //Update expanded leaf area (will decrease if LAI_live decreases)
      }
    } 
  }    
}
