#include <numeric>
#include <Rcpp.h>
#include <string.h>
#include <stdio.h>
#include "phenology.h"
#include "root.h"
#include "incgamma.h"
using namespace Rcpp;


/**
 *  Utils
 */
int findRowIndex(int sp, DataFrame SpParams) {
  IntegerVector spIndexSP = SpParams["SpIndex"];
  for(int i=0;i<spIndexSP.length();i++) if(spIndexSP[i]==sp) return(i);
  stop("Species code not found in SpParams");
  return(NA_INTEGER);
}

double leafAreaProportion(double z1, double z2, double zmin, double zmax) {
  double mu = (zmax+zmin)/2.0;
  double sd15 = (zmax-zmin)/2.0;
  double sd = sd15/1.5;
  z1 = std::max(z1, zmin);
  z2 = std::max(z2, zmin);
  z1 = std::min(z1, zmax);
  z2 = std::min(z2, zmax);
  double x1 = (z1-mu)/sd;
  double x2 = (z2-mu)/sd;
  double p1 = 0.5*(1.0+errorfunction(x1/sqrt(2.0), false, false));
  double p2 = 0.5*(1.0+errorfunction(x2/sqrt(2.0), false, false));
  double v = (p2-p1)/0.8663856; //truncated to -1.5 to 1.5
  return(v);
}



// [[Rcpp::export("plant_ID")]]
CharacterVector cohortIDs(List x) {
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  DataFrame shrubData = Rcpp::as<Rcpp::DataFrame>(x["shrubData"]);
  int ntree = treeData.nrows();
  int nshrub = shrubData.nrows();
  int numCohorts  = ntree+nshrub;
  IntegerVector treeSP = treeData["Species"];
  IntegerVector shrubSP = shrubData["Species"];  
  
  CharacterVector IDs(numCohorts);
  for(int i=0;i<ntree;i++) {
    char Result[16];
    sprintf(Result, "T%d_%d", i+1, treeSP[i]);
    IDs[i] = Result;
  }
  for(int i=0;i<nshrub;i++) {
    char Result[16];
    sprintf(Result, "S%d_%d", i+1, shrubSP[i]);
    IDs[ntree+i] =Result;
  }
  return(IDs);
}

// [[Rcpp::export("plant_parameter")]]
NumericVector cohortNumericParameter(List x, DataFrame SpParams, String parName){
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  DataFrame shrubData = Rcpp::as<Rcpp::DataFrame>(x["shrubData"]);
  NumericVector tSP = treeData["Species"];
  NumericVector shSP = shrubData["Species"];
  NumericVector par(tSP.size()+shSP.size());
  NumericVector parSP = SpParams[parName];
  for(int i=0;i<tSP.size();i++) {
    int iSP = findRowIndex(tSP[i], SpParams);
    par[i] = parSP[iSP];
  }
  for(int i=0;i<shSP.size();i++) {
    int iSP = findRowIndex(shSP[i], SpParams);
    par[i+tSP.size()] = parSP[iSP];
  }
  par.attr("names") = cohortIDs(x);
  return(par);
}

NumericVector cohortNumericParameter(IntegerVector SP, DataFrame SpParams, String parName){
  NumericVector par(SP.size());
  NumericVector parSP = SpParams[parName];
  for(int i=0;i<SP.size();i++) {
    int iSP = findRowIndex(SP[i], SpParams);
    par[i] = parSP[iSP];
  }
  return(par);
}

// [[Rcpp::export("plant_characterParameter")]]
CharacterVector cohortCharacterParameter(List x, DataFrame SpParams, String parName){
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  DataFrame shrubData = Rcpp::as<Rcpp::DataFrame>(x["shrubData"]);
  NumericVector tSP = treeData["Species"];
  NumericVector shSP = shrubData["Species"];
  CharacterVector par(tSP.size()+shSP.size());
  CharacterVector parSP = SpParams[parName];
  for(int i=0;i<tSP.size();i++) {
    int iSP = findRowIndex(tSP[i], SpParams);
    par[i] = parSP[iSP];
  }
  for(int i=0;i<shSP.size();i++) {
    int iSP = findRowIndex(shSP[i], SpParams);
    par[i+tSP.size()] = parSP[iSP];
  }
  par.attr("names") = cohortIDs(x);
  return(par);
}

CharacterVector cohortCharacterParameter(IntegerVector SP, DataFrame SpParams, String parName){
  CharacterVector par(SP.size());
  CharacterVector parSP = SpParams[parName];
  for(int i=0;i<SP.size();i++) {
    int iSP = findRowIndex(SP[i], SpParams);
    par[i] = parSP[iSP];
  }
  return(par);
}

IntegerVector uniqueSpp(IntegerVector sp) {
  int nunique = 0;
  IntegerVector spUnique(sp.size());
  if(sp.size()>0) {
    nunique = 1;
    spUnique[0] = sp[0];
    for(int i=1;i<sp.size();i++) {
      bool f = false;
      for(int j=0; j<i;j++) if(sp[i]==sp[j]) f = true;
      if(!f) {
        spUnique[nunique] = sp[i];
        nunique +=1;
      }
    }
  }
  if(nunique>0) {
    IntegerVector spUniqueDef(nunique);
    for(int i=0;i<nunique;i++) spUniqueDef[i] = spUnique[i];
    spUnique = spUniqueDef;
  }
  return(spUnique);
}

NumericVector sumBySpecies(NumericVector x, IntegerVector sp, DataFrame SpParams) {
  IntegerVector uniqueSp = uniqueSpp(sp);
  NumericVector sba(uniqueSp.size(),0.0);
  for(int i=0;i<sp.size();i++) {
    for(int j=0;j<uniqueSp.size();j++) {
      if(sp[i]==uniqueSp[j]) {
        sba[j] +=x[i];
      }
    }
  }
  sba.attr("names") = cohortCharacterParameter(uniqueSp, SpParams, "Name");
  return(sba);
}


/** 
 * Species
 */
// [[Rcpp::export("plant_species")]]
IntegerVector cohortSpecies(List x) {
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  DataFrame shrubData = Rcpp::as<Rcpp::DataFrame>(x["shrubData"]);
  int ntree = treeData.nrows();
  int nshrub = shrubData.nrows();
  IntegerVector treeSP = treeData["Species"];
  IntegerVector shrubSP = shrubData["Species"];  
  int numCohorts  = ntree+nshrub;
  IntegerVector SP(numCohorts);
  for(int i=0;i<ntree;i++) {
    SP[i] = treeSP[i];
  }
  for(int i=0;i<nshrub;i++) {
    SP[ntree+i] = shrubSP[i];
  }
  SP.attr("names") = cohortIDs(x);
  return(SP);
}

// [[Rcpp::export("plant_speciesName")]]
CharacterVector cohortSpeciesName(List x, DataFrame SpParams) {
  CharacterVector sn = cohortCharacterParameter(x,SpParams, "Name");
  sn.attr("names") = cohortIDs(x);
  return(sn);
}

/**
 *  Basal area (BA, in m2/ha)
 */
// [[Rcpp::export(".treeBasalArea")]]
NumericVector treeBasalArea(NumericVector N, NumericVector dbh) {
  int ncoh = N.size(); //N is density of individuals (ind/ha) in the cell
  NumericVector BA(ncoh); 
  for(int i=0;i<ncoh;i++) {
    BA[i] = N[i]*3.141593*pow(dbh[i]/200,2.0); //Basal area in m2/ha
  }
  return(BA);
}
NumericVector largerTreeBasalArea(NumericVector N, NumericVector dbh) {
  int ncoh = N.size();
  NumericVector BA = treeBasalArea(N, dbh); 
  NumericVector ltBA(ncoh);
  for(int i=0;i<ncoh;i++) {
    ltBA[i] = 0.0;
    for(int j=0;j<ncoh;j++) {
      if(i==j) ltBA[i] += BA[j]/2.0; //add half of its own basal area
      else if(dbh[j]>dbh[i]) ltBA[i] += BA[j];
    }
  }
  return(ltBA);
}
NumericVector treeCohortBasalArea(List x) {
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  NumericVector tba = treeBasalArea(treeData["N"], treeData["DBH"]);
  return(tba);  
}

// Basal area of plant cohorts. NA values are returned for shrubs
// [[Rcpp::export("plant_basalArea")]]
NumericVector cohortBasalArea(List x) {
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  DataFrame shrubData = Rcpp::as<Rcpp::DataFrame>(x["shrubData"]);
  NumericVector tba = treeBasalArea(treeData["N"], treeData["DBH"]);
  NumericVector ba(tba.size()+shrubData.nrows(), NA_REAL);
  for(int i=0;i<tba.size();i++) {
    ba[i] = tba[i];
  }
  ba.attr("names") = cohortIDs(x);
  return(ba);
}


// [[Rcpp::export("species_basalArea")]]
NumericVector speciesBasalArea(List x, DataFrame SpParams) {
  NumericVector cBA = cohortBasalArea(x);
  return(sumBySpecies(cBA, cohortSpecies(x), SpParams));
}

// [[Rcpp::export("plant_largerTreeBasalArea")]]
NumericVector cohortLargerTreeBasalArea(List x) {
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  DataFrame shrubData = Rcpp::as<Rcpp::DataFrame>(x["shrubData"]);
  NumericVector tba = largerTreeBasalArea(treeData["N"], treeData["DBH"]);
  NumericVector ba(tba.size()+shrubData.nrows(), NA_REAL);
  for(int i=0;i<tba.size();i++) {
    ba[i] = tba[i];
  }
  ba.attr("names") = cohortIDs(x);
  return(ba);
}


NumericVector dbhClassBasalArea(List x, NumericVector DBHbreaks) {
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  NumericVector treeDBH = treeData["DBH"];
  NumericVector tba = treeBasalArea(treeData["N"], treeDBH);
  int nclasses = DBHbreaks.size()-1;
  int ntree = treeData.nrows();
  NumericVector dcBA(nclasses);
  for(int i=0;i<ntree;i++) {
    for(int j=0;j<nclasses;j++) if((treeDBH[i]>=DBHbreaks[j]) & (treeDBH[i]<DBHbreaks[j+1])) dcBA[j] += tba[i]; 
  }
  return(dcBA);
}
// [[Rcpp::export("stand_basalArea")]]
double standBasalArea(List x) {
  NumericVector ba = cohortBasalArea(x);
  double tba = 0.0;
  for(int i=0;i<ba.size();i++){if(!NumericVector::is_na(ba[i])) tba+=ba[i];}
  return(tba);
}
double standBasalAreaForMinDBH(List x, double minDBH) {
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  NumericVector tba = treeBasalArea(treeData["N"], treeData["DBH"]);
  NumericVector treeDBH = treeData["DBH"];
  double ba = 0.0;
  for(int i=0;i<tba.size();i++) {
    if(treeDBH[i]>=minDBH) ba += tba[i];
  }
  return(ba);  
}

/**
 *  Density (in number of individuals/ha)
 */
double treeDensity(List x) {
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  NumericVector N = treeData["N"];
  double density = sum(N); //Density in indivudals per hectare (10000 m2)
  return(density);  
}
double minDBHDensity(List x, double minDBH) {
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  NumericVector treeDBH = treeData["DBH"];
  NumericVector treeN = treeData["N"];
  NumericVector treeDens = treeN;
  int ntree = treeData.nrows();
  double dens =0.0;
  for(int i=0;i<ntree;i++) {
    if(treeDBH[i]>=minDBH) dens += treeDens[i]; 
  }
  return(dens);
}
NumericVector dbhClassDensity(List x, NumericVector DBHbreaks) {
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  NumericVector treeDBH = treeData["DBH"];
  NumericVector treeN = treeData["N"];
  NumericVector treeDens = treeN;
  int nclasses = DBHbreaks.size()-1;
  int ntree = treeData.nrows();
  NumericVector dcDens(nclasses);
  for(int i=0;i<ntree;i++) {
    for(int j=0;j<nclasses;j++) if((treeDBH[i]>=DBHbreaks[j]) & (treeDBH[i]<DBHbreaks[j+1])) dcDens[j] += treeDens[i]; 
  }
  return(dcDens);
}

//area of an individual (in m2)
NumericVector shrubIndividualArea(IntegerVector SP, NumericVector Cover, NumericVector H, DataFrame SpParams){
  NumericVector aShrubArea = cohortNumericParameter(SP,SpParams, "a_ash");
  int ncoh = SP.size();
  NumericVector areaind(ncoh);
  for(int i=0;i<ncoh;i++) {
    if((!NumericVector::is_na(Cover[i])) & (!NumericVector::is_na(H[i]))) {
      areaind[i] = aShrubArea[i]*pow(H[i],2.0)/10000.0; 
    }
  }
  return(areaind);
}

/*
 * Cohort density in ind/ha
 */
// [[Rcpp::export("plant_density")]]
NumericVector cohortDensity(List x, DataFrame SpParams) {
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  DataFrame shrubData = Rcpp::as<Rcpp::DataFrame>(x["shrubData"]);
  int ntree = treeData.nrows();
  int nshrub = shrubData.nrows();
  NumericVector treeN = treeData["N"];
  IntegerVector shrubSP = shrubData["Species"];
  NumericVector shrubCover = shrubData["Cover"];
  NumericVector shrubHeight = shrubData["Height"];
  int numCohorts  = ntree+nshrub;
  NumericVector N(numCohorts);
  for(int i=0;i<ntree;i++) {
    N[i] = treeN[i];
  }
  NumericVector shrubArea = shrubIndividualArea(shrubSP,shrubCover, shrubHeight, SpParams);
  for(int i=0;i<nshrub;i++) {
    N[ntree+i] = 10000.0*(shrubCover[i]/(100.0*shrubArea[i]));
  }
  N.attr("names") = cohortIDs(x);
  return(N);
}

// [[Rcpp::export("species_density")]]
NumericVector speciesDensity(List x, DataFrame SpParams) {
  NumericVector d = cohortDensity(x, SpParams);
  return(sumBySpecies(d, cohortSpecies(x), SpParams));
}
/** 
 * Height
 */
double maxCohortHeight(List x) {
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  DataFrame shrubData = Rcpp::as<Rcpp::DataFrame>(x["shrubData"]);
  NumericVector treeH = treeData["Height"];
  NumericVector shrubH = shrubData["Height"];  

  int ntree = treeData.nrows();
  int nshrub = shrubData.nrows();
  double height =0.0;
  for(int i=0;i<ntree;i++) if(treeH[i]>=height) height += treeH[i]; 
  for(int i=0;i<nshrub;i++) if(shrubH[i]>=height) height += shrubH[i]; 
  return(height);
}
// [[Rcpp::export("plant_height")]]
NumericVector cohortHeight(List x) {
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  DataFrame shrubData = Rcpp::as<Rcpp::DataFrame>(x["shrubData"]);
  int ntree = treeData.nrows();
  int nshrub = shrubData.nrows();
  NumericVector treeH = treeData["Height"];
  NumericVector shrubH = shrubData["Height"];  
  int numCohorts  = ntree+nshrub;
  NumericVector H(numCohorts);
  for(int i=0;i<ntree;i++) {
    H[i] = treeH[i];
  }
  for(int i=0;i<nshrub;i++) {
    H[ntree+i] = shrubH[i];
  }
  H.attr("names") = cohortIDs(x);
  return(H);
}



/** 
 * Crown ratio, crown base height, crown length
 */
// [[Rcpp::export(".shrubCrownRatio")]]
NumericVector shrubCrownRatio(IntegerVector SP, DataFrame SpParams) {
  return(cohortNumericParameter(SP, SpParams, "cr"));
}

double crownCompetitionFactor(NumericVector N, NumericVector dbh, NumericVector Acw, NumericVector Bcw) {
  int ntree = N.size();
  double ccf = 0.0;
  for(int i=0;i<ntree;i++) {
    if(!NumericVector::is_na(dbh[i])) {
      double cw = Acw[i]*pow(dbh[i], Bcw[i]);
      ccf = ccf + (N[i]*PI*pow(cw/2.0,2.0)/100.0);
    }
  }
  return(ccf);
}

double crownCompetitionFactor(IntegerVector SP, NumericVector N, NumericVector dbh, DataFrame SpParams) {
  NumericVector Acw = cohortNumericParameter(SP, SpParams, "a_cw");
  NumericVector Bcw = cohortNumericParameter(SP, SpParams, "b_cw");
  int ntree = SP.size();
  double ccf = 0.0;
  for(int i=0;i<ntree;i++) {
    double cw = Acw[i]*pow(dbh[i], Bcw[i]);
    ccf = ccf + (N[i]*PI*pow(cw/2.0,2.0)/100.0);
  }
  return(ccf);
}

NumericVector treeCrownRatio(NumericVector N, NumericVector dbh, NumericVector H, 
                             NumericVector Acw, NumericVector Bcw,
                             NumericVector Acr, NumericVector B1cr, NumericVector B2cr, NumericVector B3cr,
                             NumericVector C1cr, NumericVector C2cr) {
  NumericVector BAL = largerTreeBasalArea(N, dbh);
  double ccf = crownCompetitionFactor(N, dbh, Acw, Bcw);
  // Rcout<<ccf<<"\n";
  int ntree = N.size();
  NumericVector treeCR(ntree, NA_REAL);
  for(int i=0;i<ntree;i++) {
    if(!NumericVector::is_na(dbh[i])) {
      double lm = Acr[i]+ B1cr[i]*(H[i]/(100.0*dbh[i]))+B2cr[i]*(H[i]/100.0)+B3cr[i]*pow(dbh[i],2.0)+C1cr[i]*BAL[i]+C2cr[i]*log(ccf);
      treeCR[i] = 1.0/(1.0+ exp(-1.0*lm));
    }
  }
  return(treeCR);
}

NumericVector treeCrownRatio(IntegerVector SP, NumericVector N, NumericVector dbh, NumericVector H, DataFrame SpParams) {
  NumericVector BAL = largerTreeBasalArea(N, dbh);
  double ccf = crownCompetitionFactor(SP, N, dbh, SpParams);
  // Rcout<<ccf<<"\n";
  NumericVector Acr = cohortNumericParameter(SP, SpParams, "a_cr");
  NumericVector B1cr = cohortNumericParameter(SP, SpParams, "b_1cr");
  NumericVector B2cr = cohortNumericParameter(SP, SpParams, "b_2cr");
  NumericVector B3cr = cohortNumericParameter(SP, SpParams, "b_3cr");
  NumericVector C1cr = cohortNumericParameter(SP, SpParams, "c_1cr");
  NumericVector C2cr = cohortNumericParameter(SP, SpParams, "c_2cr");
  int ntree = SP.size();
  NumericVector treeCR(ntree);
  for(int i=0;i<ntree;i++) {
    double lm = Acr[i]+ B1cr[i]*(H[i]/(100.0*dbh[i]))+B2cr[i]*(H[i]/100.0)+B3cr[i]*pow(dbh[i],2.0)+C1cr[i]*BAL[i]+C2cr[i]*log(ccf);
    treeCR[i] = 1.0/(1.0+ exp(-1.0*lm));
  }
  return(treeCR);
}

// [[Rcpp::export("plant_crownRatio")]]
NumericVector cohortCrownRatio(List x, DataFrame SpParams) {
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  DataFrame shrubData = Rcpp::as<Rcpp::DataFrame>(x["shrubData"]);
  int ntree = treeData.nrows();
  int nshrub = shrubData.nrows();
  IntegerVector shrubSP = shrubData["Species"];  
  NumericVector crSh = cohortNumericParameter(shrubSP, SpParams, "cr");
  int numCohorts  = ntree+nshrub;
  NumericVector treeCR = treeCrownRatio(treeData["Species"],treeData["N"], treeData["DBH"], treeData["Height"], SpParams);
  NumericVector CR(numCohorts);
  for(int i=0;i<ntree;i++) {
    CR[i] = treeCR[i];
  }
  for(int i=0;i<nshrub;i++) {
    CR[ntree+i] = crSh[i];
  }
  CR.attr("names") = cohortIDs(x);
  return(CR);
}
// [[Rcpp::export("plant_crownBaseHeight")]]
NumericVector cohortCrownBaseHeight(List x, DataFrame SpParams) {
  NumericVector CR = cohortCrownRatio(x, SpParams);
  NumericVector H = cohortHeight(x);
  int numCohorts = H.size();
  NumericVector CBH(numCohorts);
  for(int i=0;i<numCohorts;i++) {
    CBH[i] = H[i]*(1.0-CR[i]);
  }
  CBH.attr("names") = cohortIDs(x);
  return(CBH);
}
// [[Rcpp::export("plant_crownLength")]]
NumericVector cohortCrownLength(List x, DataFrame SpParams) {
  NumericVector CR = cohortCrownRatio(x, SpParams);
  NumericVector H = cohortHeight(x);
  int numCohorts = H.size();
  NumericVector CL(numCohorts);
  for(int i=0;i<numCohorts;i++) {
    CL[i] = H[i]*CR[i];
  }
  CL.attr("names") = cohortIDs(x);
  return(CL);
}


/**
 * Foliar biomass (in kg/m2)
 */
// [[Rcpp::export(".treeFoliarBiomass")]]
NumericVector treeFoliarBiomass(IntegerVector SP, NumericVector N, NumericVector dbh, DataFrame SpParams, double gdd = NA_REAL){
  NumericVector afbt = cohortNumericParameter(SP, SpParams, "a_fbt");
  NumericVector bfbt = cohortNumericParameter(SP, SpParams, "b_fbt");
  NumericVector cfbt = cohortNumericParameter(SP, SpParams, "c_fbt");
  NumericVector dfbt = cohortNumericParameter(SP, SpParams, "d_fbt");
  NumericVector ltba = largerTreeBasalArea(N,dbh);
  int ncoh = N.size();
  NumericVector lb(ncoh);
  for(int i=0;i<ncoh;i++) {
    lb[i] = ((N[i]/10000)*afbt[i]*pow(dbh[i], bfbt[i])*exp(cfbt[i]*ltba[i])*pow(dbh[i], dfbt[i]*ltba[i]));
  }
  if(!NumericVector::is_na(gdd)) {
    NumericVector Sgdd = cohortNumericParameter(SP, SpParams, "Sgdd");
    for(int i=0;i<ncoh;i++) {
      if(!NumericVector::is_na(SP[i])) lb[i] = lb[i]*leafDevelopmentStatus(Sgdd[i], gdd);
    }
  }
  
  return(lb);
}

// [[Rcpp::export(".shrubFoliarBiomass")]]
NumericVector shrubFoliarBiomass(IntegerVector SP, NumericVector Cover, NumericVector H, NumericVector CR, DataFrame SpParams, double gdd = NA_REAL){
  NumericVector aShrubArea = cohortNumericParameter(SP, SpParams, "a_ash");
  NumericVector aShrubFuel = cohortNumericParameter(SP, SpParams, "a_bsh");
  NumericVector bShrubFuel = cohortNumericParameter(SP, SpParams, "b_bsh");
  NumericVector Sgdd = cohortNumericParameter(SP, SpParams, "Sgdd");
  NumericVector pDead = cohortNumericParameter(SP, SpParams, "pDead");
  NumericVector fTreeFuel = cohortNumericParameter(SP, SpParams, "r635");
  int ncoh = SP.size();
  double W = 0.0; //Fine fuel
  NumericVector fb(ncoh);
  NumericVector areaind = shrubIndividualArea(SP,Cover,H,SpParams); 
  double volind = NA_REAL,weightkgind = NA_REAL;
  for(int i=0;i<ncoh;i++) {
    if((!NumericVector::is_na(Cover[i])) & (!NumericVector::is_na(H[i]))) {
      volind = areaind[i]*((H[i]-(H[i]*(1.0-CR[i])))/100.0); //Crown phytovolume of an individual (in m3)
      weightkgind = aShrubFuel[i]*pow(volind,bShrubFuel[i]); //Fuel (in kg) of an individual (includes dead fuels)
      weightkgind = weightkgind - (weightkgind*pDead[i]); //Removes dead fuels
      if(areaind[i]>0.0) {
        // multiply by 'number of individuals' per m2 
        W = weightkgind*(Cover[i]/(100*areaind[i]));  //Fine fuel (kg/m2)
        fb[i] = W/fTreeFuel[i]; //Foliar biomass (kg/m2)
        // Rcout<<Cover[i]<<" "<<(Cover[i]/(100*areaind))<<" "<< W<< " "<< fb[i]<<"\n";
        if(!NumericVector::is_na(gdd)) { //Apply phenology correction to foliar fuels
          fb[i] = fb[i]*leafDevelopmentStatus(Sgdd[i], gdd); 
        } 
      }
    }
    else fb[i] = NA_REAL;
  }
  return(fb);
}
// [[Rcpp::export("plant_foliarBiomass")]]
NumericVector cohortFoliarBiomass(List x, DataFrame SpParams, double gdd = NA_REAL) {
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  DataFrame shrubData = Rcpp::as<Rcpp::DataFrame>(x["shrubData"]);
  NumericVector tFB = treeFoliarBiomass(treeData["Species"], treeData["N"], treeData["DBH"], SpParams, gdd);
  IntegerVector shSP = shrubData["Species"];
  NumericVector shCR = shrubCrownRatio(shSP, SpParams);
  NumericVector shFB = shrubFoliarBiomass(shSP, shrubData["Cover"], shrubData["Height"], shCR, SpParams, gdd);
  NumericVector FB(tFB.size()+shFB.size());
  for(int i=0;i<tFB.size();i++) {
    FB[i] = tFB[i];
  }
  for(int i=0;i<shFB.size();i++) {
    FB[i+tFB.size()] = shFB[i];
  }
  FB.attr("names") = cohortIDs(x);
  return(FB);
}

// [[Rcpp::export("species_foliarBiomass")]]
NumericVector speciesFoliarBiomass(List x, DataFrame SpParams, double gdd = NA_REAL) {
  NumericVector fb = cohortFoliarBiomass(x, SpParams, gdd);
  return(sumBySpecies(fb, cohortSpecies(x), SpParams));
}

// [[Rcpp::export("stand_foliarBiomass")]]
double standFoliarBiomass(List x, DataFrame SpParams, double gdd = NA_REAL) {
  NumericVector fb = cohortFoliarBiomass(x, SpParams, gdd);
  double tfb= 0.0;
  for(int i=0;i<fb.size();i++){if(!NumericVector::is_na(fb[i])) tfb+=fb[i];}
  return(tfb);
}

/**
 *  Cover (percent)
 */
// [[Rcpp::export(".shrubCover")]]
double shrubCover(List x, double excludeMinHeight = 0.0) {
  DataFrame shrubData = Rcpp::as<Rcpp::DataFrame>(x["shrubData"]);
  NumericVector height = shrubData["Height"];
  NumericVector cover = shrubData["Cover"];
  double cov = 0.0;
  int n = height.size();
  for(int i=0;i<n;i++) {
    if(height[i]>excludeMinHeight) cov += cover[i];
  }
  return(cov);
}
// [[Rcpp::export("plant_cover")]]
NumericVector cohortCover(List x) {
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  DataFrame shrubData = Rcpp::as<Rcpp::DataFrame>(x["shrubData"]);
  NumericVector cover = shrubData["Cover"];
  NumericVector vol(treeData.nrows()+shrubData.nrows(), NA_REAL);
  for(int i=0;i<cover.size();i++) {
    vol[i+treeData.nrows()] = cover[i];
  }
  vol.attr("names") = cohortIDs(x);
  return(vol);
}

// [[Rcpp::export("species_cover")]]
NumericVector speciesCover(List x, DataFrame SpParams) {
  NumericVector cc = cohortCover(x);
  return(sumBySpecies(cc, cohortSpecies(x), SpParams));
}

/**
 *  Shrub phytovolume (in m3/m2)
 */
// [[Rcpp::export(".shrubCrownPhytovolume")]]
NumericVector shrubCrownPhytovolume(IntegerVector SP, NumericVector Cover, NumericVector H, NumericVector CR, DataFrame SpParams){
  NumericVector aShrubArea = cohortNumericParameter(SP, SpParams, "a_ash");
  int ncoh = Cover.size();
  NumericVector vol(ncoh);
  for(int i=0;i<ncoh;i++) {
    if((!NumericVector::is_na(Cover[i]))& (!NumericVector::is_na(H[i]))) {
      double areaind = aShrubArea[i]*pow(H[i],2.0)/10000.0; //area of an individual (in m2)
      double volind = areaind*((H[i]-(H[i]*(1.0-CR[i])))/100.0); //Crown phytovolume of an individual (in m3)
      // Rcout <<areaind<<" "<< volind<<"\n";
      vol[i] = volind * (Cover[i]/(100*areaind));
    } else {
      vol[i] = NA_REAL;
    }
  }
  return(vol);
}
// [[Rcpp::export("plant_phytovolume")]]
NumericVector cohortPhytovolume(List x, DataFrame SpParams) {
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  DataFrame shrubData = Rcpp::as<Rcpp::DataFrame>(x["shrubData"]);
  IntegerVector SP = shrubData["Species"];
  NumericVector CR = shrubCrownRatio(SP, SpParams);
  NumericVector shvol = shrubCrownPhytovolume(SP, shrubData["Cover"], shrubData["Height"], CR, SpParams);
  NumericVector vol(treeData.nrows()+shrubData.nrows(), NA_REAL);
  for(int i=0;i<shvol.size();i++) {
    vol[i+treeData.nrows()] = shvol[i];
  }
  vol.attr("names") = cohortIDs(x);
  return(vol);
}

// [[Rcpp::export("species_phytovolume")]]
NumericVector speciesPhytovolume(List x, DataFrame SpParams) {
  NumericVector cp = cohortPhytovolume(x, SpParams);
  return(sumBySpecies(cp, cohortSpecies(x), SpParams));
}

// [[Rcpp::export("stand_phytovolume")]]
double standPhytovolume(List x, DataFrame SpParams, double gdd = NA_REAL) {
  NumericVector cp = cohortPhytovolume(x, SpParams);
  double tp= 0.0;
  for(int i=0;i<cp.size();i++){if(!NumericVector::is_na(cp[i])) tp+=cp[i];}
  return(tp);
}
/**
 * Fine fuel loading (in kg/m2)
 */
// [[Rcpp::export(".treeFuel")]]
NumericVector treeFuel(IntegerVector SP, NumericVector N, NumericVector dbh, DataFrame SpParams, double gdd = NA_REAL, bool includeDead = true){
  NumericVector fb = treeFoliarBiomass(SP, N, dbh, SpParams, NA_REAL); //Do not include phenology (to have correct estimates of branch biomass)
  NumericVector Sgdd = cohortNumericParameter(SP, SpParams, "Sgdd");
  NumericVector fTreeFuel = cohortNumericParameter(SP, SpParams, "r635");
  NumericVector pDead = cohortNumericParameter(SP, SpParams, "pDead");
  int ncoh = N.size();
  double ftf = 0.0, btf = 0.0;
  NumericVector fuel(ncoh);
  for(int i=0;i<ncoh;i++) {
    if((!NumericVector::is_na(dbh[i]))& (!NumericVector::is_na(N[i]))) {
      ftf = fb[i]; //Foliar biomass (kg per m2)
      btf = ftf*(fTreeFuel[i]-1.0); // Small branch fuels (proportion of foliar fuels)
      if(!NumericVector::is_na(gdd)) { //Apply phenology correction to foliar fuels
        ftf = ftf*leafDevelopmentStatus(Sgdd[i], gdd); 
      } 
      fuel[i] =  ftf + btf; //Tree fuel (kg per m2) is sum of both fuels
      if(includeDead) fuel[i] = fuel[i] + fuel[i]*pDead[i]; //If required add fine dead fuels (proportion of live fuels)
    }
    else fuel[i] = NA_REAL;
  }
  return(fuel);
}
// [[Rcpp::export(".shrubFuel")]]
NumericVector shrubFuel(IntegerVector SP, NumericVector Cover, NumericVector H, NumericVector CR, DataFrame SpParams, double gdd = NA_REAL, bool includeDead = true){
  NumericVector aShrubArea = cohortNumericParameter(SP, SpParams, "a_ash");
  NumericVector aShrubFuel = cohortNumericParameter(SP, SpParams, "a_bsh");
  NumericVector bShrubFuel = cohortNumericParameter(SP, SpParams, "b_bsh");
  NumericVector pDead = cohortNumericParameter(SP, SpParams, "pDead");
  NumericVector fTreeFuel = cohortNumericParameter(SP, SpParams, "r635");

  int ncoh = SP.size();
  double areaind = NA_REAL, volind = NA_REAL, weightkgind = NA_REAL;
  //W in kg/m2. Fine fuel, does not include phenology 
  NumericVector W(ncoh);
  for(int i=0;i<ncoh;i++) {
    if((!NumericVector::is_na(Cover[i])) & (!NumericVector::is_na(H[i]))) {
      areaind = aShrubArea[i]*pow(H[i],2.0)/10000.0; //area of an individual (in m2)
      volind = areaind*((H[i]-(H[i]*(1.0-CR[i])))/100.0); //Crown phytovolume of an individual (in m3)
      weightkgind = aShrubFuel[i]*pow(volind,bShrubFuel[i]); //Dry weight (in kg) of an individual
      if(!includeDead) weightkgind = weightkgind - (weightkgind*pDead[i]); //Remove dead fuels if asked
      if(areaind>0.0) {
        // multiply by 'number of individuals' per m2 
        W[i] = weightkgind*(Cover[i]/(100*areaind)); 
      }
    }
    else W[i] = NA_REAL;
  }
  //Remove (if necessary), the weight due to leaves that are not there
   if(!NumericVector::is_na(gdd)) {
    double fsf = 0.0, bsf = 0.0;
    NumericVector Sgdd = cohortNumericParameter(SP, SpParams, "Sgdd");
    for(int i=0;i<ncoh;i++) {
      fsf = W[i]/fTreeFuel[i]; //foliar biomass
      bsf = W[i] - fsf; //branch biomass
      W[i] = bsf + fsf*leafDevelopmentStatus(Sgdd[i], gdd); 
    }
  }
  return(W);
}
// [[Rcpp::export("plant_fuel")]]
NumericVector cohortFuel(List x, DataFrame SpParams, double gdd = NA_REAL, bool includeDead = true){
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  DataFrame shrubData = Rcpp::as<Rcpp::DataFrame>(x["shrubData"]);
  NumericVector tFuel = treeFuel(treeData["Species"], treeData["N"], treeData["DBH"], SpParams, gdd, includeDead);
  IntegerVector shSP = shrubData["Species"];
  NumericVector shCR = shrubCrownRatio(shSP, SpParams);
  NumericVector shFuel = shrubFuel(shSP, shrubData["Cover"], shrubData["Height"], shCR, SpParams, gdd, includeDead);
  NumericVector fuel(tFuel.size()+shFuel.size());
  for(int i=0;i<tFuel.size();i++) {
    fuel[i] = tFuel[i];
  }
  for(int i=0;i<shFuel.size();i++) {
    fuel[i+tFuel.size()] = shFuel[i];
  }
  fuel.attr("names") = cohortIDs(x);
  return(fuel);
}


// [[Rcpp::export("species_fuel")]]
NumericVector speciesFuel(List x, DataFrame SpParams, double gdd = NA_REAL, bool includeDead = true) {
  NumericVector cf = cohortFuel(x, SpParams, gdd, includeDead);
  return(sumBySpecies(cf, cohortSpecies(x), SpParams));
}
// [[Rcpp::export("stand_fuel")]]
double standFuel(List x, DataFrame SpParams, double gdd = NA_REAL, bool includeDead = true) {
  NumericVector cf = cohortFuel(x, SpParams, gdd, includeDead);
  double tf= 0.0;
  for(int i=0;i<cf.size();i++){if(!NumericVector::is_na(cf[i])) tf+=cf[i];}
  return(tf);
}

/**
 * Cohort equilibrium leaf litter (in kg/m2)
 */
// [[Rcpp::export("plant_equilibriumLeafLitter")]]
NumericVector cohortEquilibriumLeafLitter(List x, DataFrame SpParams, double AET = 800) {
  NumericVector fb = cohortFoliarBiomass(x, SpParams);
  NumericVector ld = cohortNumericParameter(x, SpParams, "LeafDuration");
  NumericVector lignin = cohortNumericParameter(x, SpParams, "LigninPercent");
  int ncoh = fb.size();
  NumericVector eqli(ncoh);
  double ki = 0.0;
  for(int i=0;i<ncoh;i++) {
    ki = (-0.5365+0.00241*AET) - (-0.01586+0.000056*AET)*lignin[i];//Meentemeyer (1978)
    // Rcout<<ki<<"\n";
    eqli[i] = fb[i]/(ld[i]*ki);
  }
  eqli.attr("names") = cohortIDs(x);
  return(eqli);
}

/**
 * Cohort equilibrium small branch (6.35mm) litter (in kg/m2)
 */
// [[Rcpp::export("plant_equilibriumSmallBranchLitter")]]
NumericVector cohortEquilibriumSmallBranchLitter(List x, DataFrame SpParams, double smallBranchDecompositionRate = 0.81) {
  NumericVector fu = cohortFuel(x, SpParams);
  NumericVector fb = cohortFoliarBiomass(x, SpParams);
  NumericVector ld = cohortNumericParameter(x, SpParams, "LeafDuration");
  int ncoh = fb.size();
  NumericVector eqli(ncoh);
  for(int i=0;i<ncoh;i++) {
    eqli[i] = (fu[i]-fb[i])/((ld[i]*2.0)*smallBranchDecompositionRate);
  }
  eqli.attr("names") = cohortIDs(x);
  return(eqli);
}


/**
 *  Leaf Area Index (LAI)
 */

// [[Rcpp::export(".treeLAI")]]
NumericVector treeLAI(IntegerVector SP, NumericVector N, NumericVector dbh, DataFrame SpParams, double gdd = NA_REAL){
  NumericVector SLA = cohortNumericParameter(SP, SpParams, "SLA"); // m2/kg (=mg/mm2)
  NumericVector lb = treeFoliarBiomass(SP, N, dbh, SpParams, gdd); //kg per m2
  int ncoh = N.size();
  NumericVector lai(ncoh);
  for(int i=0;i<ncoh;i++) {
     lai[i] = (SLA[i]*lb[i]); 
  }
  return(lai);
}
// [[Rcpp::export(".shrubLAI")]]
NumericVector shrubLAI(IntegerVector SP, NumericVector Cover, NumericVector H, DataFrame SpParams, double gdd = NA_REAL){
  NumericVector SLA = cohortNumericParameter(SP, SpParams, "SLA"); // m2/kg (=mg/mm2)
  NumericVector CR = shrubCrownRatio(SP, SpParams);
  NumericVector lb = shrubFoliarBiomass(SP, Cover, H, CR, SpParams, gdd); //kg per m2
  int ncoh = SP.size();
  NumericVector lai(ncoh);
  for(int i=0;i<ncoh;i++) {
    lai[i] = (SLA[i]*lb[i]); 
  }
  return(lai);
}
// [[Rcpp::export("plant_LAI")]]
NumericVector cohortLAI(List x, DataFrame SpParams, double gdd = NA_REAL){
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  DataFrame shrubData = Rcpp::as<Rcpp::DataFrame>(x["shrubData"]);
  NumericVector tLAI = treeLAI(treeData["Species"], treeData["N"], treeData["DBH"], SpParams, gdd);
  NumericVector shLAI = shrubLAI(shrubData["Species"], shrubData["Cover"], shrubData["Height"], SpParams, gdd);
  NumericVector lai(tLAI.size()+shLAI.size());
  for(int i=0;i<tLAI.size();i++) {
    lai[i] = tLAI[i];
  }
  for(int i=0;i<shLAI.size();i++) {
    lai[i+tLAI.size()] = shLAI[i];
  }
  lai.attr("names") = cohortIDs(x);
  return(lai);
}


// [[Rcpp::export("species_LAI")]]
NumericVector speciesLAI(List x, DataFrame SpParams, double gdd = NA_REAL) {
  NumericVector cl = cohortLAI(x, SpParams, gdd);
  return(sumBySpecies(cl, cohortSpecies(x), SpParams));
}

// [[Rcpp::export("stand_LAI")]]
double standLAI(List x, DataFrame SpParams, double gdd = NA_REAL) {
  NumericVector cl = cohortLAI(x, SpParams, gdd);
  double tl= 0.0;
  for(int i=0;i<cl.size();i++){if(!NumericVector::is_na(cl[i])) tl+=cl[i];}
  return(tl);
}

// [[Rcpp::export(".LAIdistributionVectors")]]
NumericMatrix LAIdistributionVectors(NumericVector z, NumericVector LAI, NumericVector H, NumericVector CR) {
  int nh = z.size();
  int ncoh = LAI.size();
  // double h1, h2;
  NumericMatrix LAIdist(nh-1, ncoh);
  for(int ci=0;ci<ncoh;ci++) {
    double cbh = H[ci]*(1.0-CR[ci]);
    // double laih = LAI[ci]/(H[ci]-cbh);
    for(int hi=0;hi<(nh-1);hi++) {
      // h1 = std::max(z[hi],cbh);
      // h2 = std::min(z[hi+1],H[ci]);
      LAIdist(hi,ci) = LAI[ci]*leafAreaProportion(z[hi],z[hi+1], cbh,H[ci]);
      // LAIdist(hi,ci) =laih*std::max(0.0, h2 - h1);
    }
  }
  return(LAIdist);
}
// [[Rcpp::export(".LAIdistribution")]]
NumericMatrix LAIdistribution(NumericVector z, List x, DataFrame SpParams, double gdd = NA_REAL) {
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  DataFrame shrubData = Rcpp::as<Rcpp::DataFrame>(x["shrubData"]);
  int ntree = treeData.nrows();
  int nshrub = shrubData.nrows();

  IntegerVector treeSP = Rcpp::as<Rcpp::IntegerVector>(treeData["Species"]);
  NumericVector treeH = treeData["Height"];
  IntegerVector shrubSP = Rcpp::as<Rcpp::IntegerVector>(shrubData["Species"]);
  NumericVector shrubH = shrubData["Height"];  
  
  NumericVector LAI = cohortLAI(x, SpParams, gdd);
  int numCohorts  = ntree+nshrub;
  NumericVector H(numCohorts);
  IntegerVector SP(numCohorts);
  for(int i=0;i<ntree;i++) {
    SP[i] = (int) treeSP[i];
    H[i] = treeH[i];
  }
  for(int i=0;i<nshrub;i++) {
    SP[ntree+i] = (int) shrubSP[i];
    H[ntree+i] = shrubH[i];
  }
  NumericVector CR = cohortCrownRatio(x, SpParams);

  return(LAIdistributionVectors(z, LAI, H, CR));
}
// [[Rcpp::export(".LAIprofileVectors")]]
NumericVector LAIprofileVectors(NumericVector z, NumericVector LAI, NumericVector H, NumericVector CR) {
  int nh = z.size();
  int ncoh = LAI.size();
  // double h1, h2;
  NumericVector LAIprof(nh-1);
  for(int ci=0;ci<ncoh;ci++) {
    double cbh = H[ci]*(1.0-CR[ci]);
    // double laih = LAI[ci]/(H[ci]-cbh);
    for(int hi=0;hi<(nh-1);hi++) {
      // h1 = std::max(z[hi],cbh);
      // h2 = std::min(z[hi+1],H[ci]);
      LAIprof[hi] +=LAI[ci]*leafAreaProportion(z[hi],z[hi+1], cbh,H[ci]);
    }
  }
  return(LAIprof);
}
// [[Rcpp::export(".LAIprofile")]]
NumericVector LAIprofile(NumericVector z, List x, DataFrame SpParams, double gdd = NA_REAL) {
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  DataFrame shrubData = Rcpp::as<Rcpp::DataFrame>(x["shrubData"]);
  int ntree = treeData.nrows();
  int nshrub = shrubData.nrows();

  IntegerVector treeSP = Rcpp::as<Rcpp::IntegerVector>(treeData["Species"]);
  NumericVector treeH = treeData["Height"];
  IntegerVector shrubSP = Rcpp::as<Rcpp::IntegerVector>(shrubData["Species"]);
  NumericVector shrubH = shrubData["Height"];  

  NumericVector CR = cohortCrownRatio(x, SpParams);
  NumericVector LAI = cohortLAI(x, SpParams, gdd);
  int numCohorts  = ntree+nshrub;
  NumericVector H(numCohorts);
  IntegerVector SP(numCohorts);
  for(int i=0;i<ntree;i++) {
    SP[i] = (int) treeSP[i];
    H[i] = treeH[i];
  }
  for(int i=0;i<nshrub;i++) {
    SP[ntree+i] = (int) shrubSP[i];
    H[ntree+i] = shrubH[i];
  }
  return(LAIprofileVectors(z, LAI, H, CR));
}




/**
 *  Delete plant cohorts
 */
void deleteTreeCohort(List x, int treeCohort) {
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  IntegerVector SP = treeData["Species"];
  NumericVector N  = treeData["N"];
  NumericVector DBH = treeData["DBH"];
  NumericVector Height = treeData["Height"];
  NumericVector Z50 = treeData["Z50"];
  NumericVector Z95 = treeData["Z95"];
  NumericVector NSC = treeData["NSC"];
  NumericVector pEmb = treeData["pEmb"];
  int n = SP.size();
  IntegerVector newSP = IntegerVector(n-1);
  NumericVector newN  = NumericVector(n-1);
  NumericVector newDBH = NumericVector(n-1);
  NumericVector newHeight = NumericVector(n-1);
  NumericVector newZ50 = NumericVector(n-1);
  NumericVector newZ95 = NumericVector(n-1);
  NumericVector newNSC = NumericVector(n-1);
  NumericVector newpEmb = NumericVector(n-1);
  for(int i=0;i<treeCohort;i++) {
    newSP[i] = SP[i];
    newN[i] = N[i];
    newDBH[i] = DBH[i];
    newHeight[i] = Height[i];
    newZ50[i] = Z50[i];
    newZ95[i] = Z95[i];
    newNSC[i] = NSC[i];
    newpEmb[i] = pEmb[i];
  }
  for(int i=(treeCohort+1);i<n;i++) {
    newSP[i-1] = SP[i];
    newN[i-1] = N[i];
    newDBH[i-1] = DBH[i];
    newHeight[i-1] = Height[i];
    newZ50[i-1] = Z50[i];
    newZ95[i-1] = Z95[i];
    newNSC[i-1] = NSC[i];
    newpEmb[i-1] = pEmb[i];
  }
  x["treeData"] = DataFrame::create(_["Species"] = newSP, _["N"] = newN,
                                    _["DBH"] = newDBH, _["Height"] = newHeight, 
                                    _["Z50"] = newZ50, _["Z95"] = newZ95,
                                    _["NSC"] = newNSC, _["pEmb"] = newpEmb);
}

void deleteShrubCohort(List x, int shrubCohort) {
  DataFrame shrubData = Rcpp::as<Rcpp::DataFrame>(x["shrubData"]);
  IntegerVector SP = shrubData["Species"];
  NumericVector Cover  = shrubData["Cover"];
  NumericVector Height = shrubData["Height"];
  NumericVector Z = shrubData["Z"];
  NumericVector NSC = shrubData["NSC"];
  int n = SP.size();
  IntegerVector newSP = IntegerVector(n-1);
  NumericVector newCover  = NumericVector(n-1);
  NumericVector newHeight = NumericVector(n-1);
  NumericVector newZ = NumericVector(n-1);
  NumericVector newNSC = NumericVector(n-1);
  for(int i=0;i<shrubCohort;i++) {
    newSP[i] = SP[i];
    newCover[i] = Cover[i];
    newHeight[i] = Height[i];
    newZ[i] = Z[i];
    newNSC[i] = NSC[i];
  }
  for(int i=(shrubCohort+1);i<n;i++) {
    newSP[i-1] = SP[i];
    newCover[i-1] = Cover[i];
    newHeight[i-1] = Height[i];
    newZ[i-1] = Z[i];
    newNSC[i-1] = NSC[i];
  }
  x["shrubData"] = DataFrame::create(_["Species"] = newSP, _["Cover"] = newCover,
                                    _["Height"] = newHeight, _["Z"] = newZ,
                                    _["NSC"] = newNSC);
}

/**
 *  Find plant cohorts
 */
// Finds tree cohort with minimum DBH (can exclude very small trees)
int minDBHTreeCohort(List x, double excludeMin = 0.0) {
  int treeCohort = NA_INTEGER;
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  NumericVector DBH = treeData["DBH"];
  int n = DBH.size();
  double min = 999999.0;
  for(int i=0;i<n;i++) {
    if((DBH[i]<min) & (DBH[i]>excludeMin)) {
      min = DBH[i];
      treeCohort = i;
    }
  }
  return(treeCohort);
}





// [[Rcpp::export("forest2aboveground")]]
DataFrame forest2aboveground(List x, DataFrame SpParams, double gdd = NA_REAL) {
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  DataFrame shrubData = Rcpp::as<Rcpp::DataFrame>(x["shrubData"]);
  int ntree = treeData.nrows();
  int nshrub = shrubData.nrows();
  
  NumericVector CR = cohortCrownRatio(x, SpParams);
  NumericVector LAI_live = cohortLAI(x, SpParams, NA_REAL);
  NumericVector LAI_expanded = cohortLAI(x, SpParams, gdd);
  IntegerVector SP(ntree+nshrub);
  NumericVector H(ntree+nshrub);
  
  IntegerVector treeSP = Rcpp::as<Rcpp::IntegerVector>(treeData["Species"]);
  NumericVector treeH = treeData["Height"];
  NumericVector treeDBH = treeData["DBH"];
  IntegerVector shrubSP = Rcpp::as<Rcpp::IntegerVector>(shrubData["Species"]);
  NumericVector shrubH = shrubData["Height"];  
  NumericVector shrubCover = shrubData["Cover"];  
  
  NumericVector N = cohortDensity(x, SpParams);
    
  NumericVector LAI_dead(ntree+nshrub);
  NumericVector DBH(ntree+nshrub);
  NumericVector Cover(ntree+nshrub);
  
  for(int i=0;i<ntree;i++) {
    SP[i] = treeSP[i];
    H[i] = treeH[i];
    DBH[i] = treeDBH[i];
    LAI_dead[i] = 0.0;
    Cover[i] = NA_REAL;
  }
  for(int i=0;i<nshrub;i++) {
    SP[ntree+i] = shrubSP[i];
    H[ntree+i] = shrubH[i];
    DBH[ntree+i] = NA_REAL;
    Cover[ntree+i] = shrubCover[i];
  }
  DataFrame above = DataFrame::create(_["SP"]=SP, _["N"] = N,  _["DBH"] = DBH,_["Cover"] = Cover, _["H"]=H, _["CR"] = CR, 
                    _["LAI_live"]=LAI_live, _["LAI_expanded"] = LAI_expanded, _["LAI_dead"] = LAI_dead);
  above.attr("row.names") = cohortIDs(x); //Assign cohort IDs to row.names
    
  return(above);
  
}


// [[Rcpp::export("forest2belowground")]]
NumericMatrix forest2belowground(List x, List soil, DataFrame SpParams) {
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  DataFrame shrubData = Rcpp::as<Rcpp::DataFrame>(x["shrubData"]);
  NumericVector d = soil["dVec"];
  int ntree = treeData.nrows();
  int nshrub = shrubData.nrows();
  
  int nlayers = d.size();
  NumericMatrix V(ntree+nshrub,nlayers);
  NumericVector treeZ50 = treeData["Z50"];
  NumericVector treeZ95 = treeData["Z95"];
  NumericVector shrubZ50 = shrubData["Z50"];  
  NumericVector shrubZ95 = shrubData["Z95"];  
  NumericVector Vi;
  CharacterVector ln(nlayers);
  for(int l=0;l<nlayers;l++){
    char Result[16]; 
    sprintf(Result, "%d", l+1);
    ln[l] = Result;
  }
  
  for(int i=0;i<ntree;i++) {
    Vi = ldrRS_one(treeZ50[i], treeZ95[i],d);
    V(i,_) = Vi;
  }
  for(int i=0;i<nshrub;i++) {
    Vi = ldrRS_one(shrubZ50[i],shrubZ95[i],d);
    V(ntree+i,_) = Vi;
  }
  V.attr("dimnames") = List::create(cohortIDs(x),ln);
  return(V);
}