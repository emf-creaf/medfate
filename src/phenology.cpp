#include <Rcpp.h>
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
// [[Rcpp::export("pheno_leafDevelopmentStatus")]]
NumericVector leafDevelopmentStatus(NumericVector Sgdd, NumericVector gdd, double unfoldingDD = 300.0) {
  NumericVector phe(Sgdd.size());
  for(int i=0;i<Sgdd.size();i++) phe[i] = leafDevelopmentStatus(Sgdd[i], gdd[i], unfoldingDD);
  return(phe);
}

//' @param Ssen Threshold to start leaf senescence.
//' @param sen Cumulative senescence variable.
//' @rdname pheno_updateLeaves
// [[Rcpp::export("pheno_leafSenescenceStatus")]]
LogicalVector leafSenescenceStatus(NumericVector Ssen, NumericVector sen) {
  LogicalVector phe(Ssen.size());
  for(int i=0;i<Ssen.size();i++) phe[i] = leafSenescenceStatus(Ssen[i], sen[i]);
  return(phe);
}

//' @param x An object of class \code{\link{spwbInput}}.
//' @param doy Day of the year.
//' @param photoperiod Day length (in hours).
//' @param tmean Average day temperature (in Celsius).
//' @rdname pheno_updateLeaves
// [[Rcpp::export("pheno_updatePhenology")]]
void updatePhenology(List x, int doy, double photoperiod, double tmean) {
  List control = x["control"];
  double unfoldingDD = control["unfoldingDD"];
  
  DataFrame paramsPhenology = Rcpp::as<Rcpp::DataFrame>(x["paramsPhenology"]);
  CharacterVector phenoType = paramsPhenology["PhenologyType"];
  NumericVector t0gdd = paramsPhenology["t0gdd"];
  NumericVector Sgdd = paramsPhenology["Sgdd"];
  NumericVector Tbgdd = paramsPhenology["Tbgdd"];
  NumericVector Ssen = paramsPhenology["Ssen"];
  NumericVector Phsen = paramsPhenology["Phsen"];
  NumericVector Tbsen = paramsPhenology["Tbsen"];
  NumericVector xsen = paramsPhenology["xsen"];
  NumericVector ysen = paramsPhenology["ysen"];
  
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
        gdd[j] = 0.0;
        if(photoperiod>Phsen[j]) { //Primary growth still possible until decrease of photoperiod
          sen[j] = 0.0;
          leafUnfolding[j] = true;
          leafSenescence[j] = false;
          budFormation[j] = false;
          leafDormancy[j] = false;
        } else { //Start (or continue) accumulation
          leafUnfolding[j] = false; //Primary growth has arrested
          if(!leafDormancy[j]) { // Check temperature accumulation until dormancy occurs
            double rsen = 0.0;
            if(tmean-Tbsen[j]<0.0) {
              rsen = pow(Tbsen[j]-tmean, xsen[j])*pow(photoperiod/Phsen[j], ysen[j]);
            }
            sen[j] = sen[j] + rsen;
            leafSenescence[j] = leafSenescenceStatus(Ssen[j],sen[j]);
            leafDormancy[j] = leafSenescence[j];
          }
          if(leafDormancy[j]) phi[j] = 0.0;
          budFormation[j] = !leafDormancy[j]; //Buds can be formed (i.e target leaf area) until dormancy occurs
        }
        // Rcout << doy<< " "<< photoperiod<<" "<< rsen <<" "<< sen[j]<<" "<<  leafSenescence[j] << "\n";
      } else if (doy<=200) { //Only increase in the first part of the year and if doy > t0gdd
        sen[j] = 0.0;
        budFormation[j] = false;
        leafSenescence[j] = false;
        if((tmean-Tbgdd[j]>0.0) && (doy >= ((int) t0gdd[j]))) gdd[j] = gdd[j] + (tmean - Tbgdd[j]);
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
        if(photoperiod>Phsen[j]) {
          sen[j] = 0.0;
          budFormation[j] = true;
          leafDormancy[j] = false;
        } else if (!leafDormancy[j]){
          double rsen = 0.0;
          if(tmean-Tbsen[j]<0.0) {
            rsen = pow(Tbsen[j]-tmean,2.0);
            // rsen = pow(Tbsen[j]-tmean,2.0) * pow(photoperiod/Phsen[j],2.0);
          }
          sen[j] = sen[j] + rsen;
          leafDormancy[j] = leafSenescenceStatus(Ssen[j],sen[j]);
          budFormation[j] = !leafDormancy[j];
        }
      } else if (doy<=200) { //Only increase in the first part of the year
        sen[j] = 0.0;
        budFormation[j] = false;
        if(!leafUnfolding[j]) { //Check until unfolding starts
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
    // Rcout<< j << " phi "<< phi[j] <<" ";
  }
  // Rcout<<"\n";
}

//' @param wind Average day wind speed (in m/s).
//' @param fromGrowthModel Boolean flag to indicate that routine is called from \code{\link{growth}} simulation function.
//' @rdname pheno_updateLeaves
// [[Rcpp::export("pheno_updateLeaves")]]
void updateLeaves(List x, double wind, bool fromGrowthModel) {
  List control = x["control"];
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  NumericVector LAI_live = above["LAI_live"];
  NumericVector LAI_dead = above["LAI_dead"];
  NumericVector LAI = above["LAI_expanded"];
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
    //Leaf unfolding and senescence only dealt with if called from spwb
    if(!fromGrowthModel) {
      if(phenoType[j] == "winter-deciduous" || phenoType[j] == "winter-semideciduous") {
        if((leafSenescence[j]) && (LAI[j]>0.0)) {
          double LAI_exp_prev= LAI[j]; //Store previous value
          LAI[j] = 0.0; //Update expanded leaf area (will decrease if LAI_live decreases)
          LAI_dead[j] += LAI_exp_prev;//Check increase dead leaf area if expanded leaf area has decreased
          leafSenescence[j] = false;
        } 
        else if(leafUnfolding[j]) {
          LAI[j] = LAI_live[j]*phi[j]; //Update expanded leaf area (will decrease if LAI_live decreases)
        }
      } 
    }
  }    
}
