#include <numeric>
#include <Rcpp.h>
using namespace Rcpp;




/**
 * Phenology
 */
double leafDevelopmentStatus(double Sgdd, double gdd) {
  if(Sgdd>0.0) return(std::min(std::max(gdd/Sgdd,0.0),1.0));
  return 1.0;
}
NumericVector leafDevelopmentStatus(NumericVector Sgdd, double gdd) {
  NumericVector phe(Sgdd.size());
  for(int i=0;i<Sgdd.size();i++) phe[i] = leafDevelopmentStatus(Sgdd[i], gdd);
  return(phe);
}


// [[Rcpp::export("plant.Parameter")]]
NumericVector cohortParameter(List x, DataFrame SpParams, String parName){
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  DataFrame shrubData = Rcpp::as<Rcpp::DataFrame>(x["shrubData"]);
  NumericVector tSP = treeData["Species"];
  NumericVector shSP = shrubData["Species"];
  NumericVector par(tSP.size()+shSP.size());
  NumericVector parSP = SpParams[parName];
  for(int i=0;i<tSP.size();i++) {
    par[i] = parSP[tSP[i]];
  }
  for(int i=0;i<shSP.size();i++) {
    par[i+tSP.size()] = parSP[shSP[i]];
  }
  return(par);
}
// [[Rcpp::export("plant.CharacterParameter")]]
CharacterVector cohortCharacterParameter(List x, DataFrame SpParams, String parName){
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  DataFrame shrubData = Rcpp::as<Rcpp::DataFrame>(x["shrubData"]);
  NumericVector tSP = treeData["Species"];
  NumericVector shSP = shrubData["Species"];
  CharacterVector par(tSP.size()+shSP.size());
  CharacterVector parSP = SpParams[parName];
  for(int i=0;i<tSP.size();i++) {
    par[i] = parSP[tSP[i]];
  }
  for(int i=0;i<shSP.size();i++) {
    par[i+tSP.size()] = parSP[shSP[i]];
  }
  return(par);
}

/** 
 * Species
 */
// [[Rcpp::export("plant.Species")]]
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
  return(SP);
}
// [[Rcpp::export("plant.SpeciesName")]]
CharacterVector cohortSpeciesName(List x, DataFrame SpParams) {
  CharacterVector nameSP = SpParams["Name"];
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  DataFrame shrubData = Rcpp::as<Rcpp::DataFrame>(x["shrubData"]);
  int ntree = treeData.nrows();
  int nshrub = shrubData.nrows();
  IntegerVector treeSP = treeData["Species"];
  IntegerVector shrubSP = shrubData["Species"];  
  int numCohorts  = ntree+nshrub;
  CharacterVector names(numCohorts);
  for(int i=0;i<ntree;i++) {
    names[i] = nameSP[treeSP[i]];
  }
  for(int i=0;i<nshrub;i++) {
    names[ntree+i] = nameSP[shrubSP[i]];
  }
  return(names);
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
// [[Rcpp::export("plant.BasalArea")]]
NumericVector cohortBasalArea(List x) {
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  DataFrame shrubData = Rcpp::as<Rcpp::DataFrame>(x["shrubData"]);
  NumericVector tba = treeBasalArea(treeData["N"], treeData["DBH"]);
  NumericVector ba(tba.size()+shrubData.nrows(), NA_REAL);
  for(int i=0;i<tba.size();i++) {
    ba[i] = tba[i];
  }
  return(ba);
}
// [[Rcpp::export("species.BasalArea")]]
NumericVector speciesBasalArea(List x, DataFrame SpParams) {
  NumericVector cBA = cohortBasalArea(x);
  IntegerVector sp = cohortSpecies(x);
  
  int nspp= SpParams.nrows();
  NumericVector sba(nspp);
  for(int i=0;i<sp.size();i++) sba[sp[i]] += cBA[i];
  return(sba);
}

// [[Rcpp::export("plant.LargerTreeBasalArea")]]
NumericVector cohortLargerTreeBasalArea(List x) {
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  DataFrame shrubData = Rcpp::as<Rcpp::DataFrame>(x["shrubData"]);
  NumericVector tba = largerTreeBasalArea(treeData["N"], treeData["DBH"]);
  NumericVector ba(tba.size()+shrubData.nrows(), NA_REAL);
  for(int i=0;i<tba.size();i++) {
    ba[i] = tba[i];
  }
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
// [[Rcpp::export("forest.BasalArea")]]
double forestBasalArea(List x) {
  NumericVector ba = cohortBasalArea(x);
  double tba = 0.0;
  for(int i=0;i<ba.size();i++){if(!NumericVector::is_na(ba[i])) tba+=ba[i];}
  return(tba);
}
double forestBasalAreaForMinDBH(List x, double minDBH) {
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
  NumericVector aShrubAreaSP = SpParams["a_ash"];
  int ncoh = SP.size();
  NumericVector areaind(ncoh);
  for(int i=0;i<ncoh;i++) {
    if((!NumericVector::is_na(SP[i])) & (!NumericVector::is_na(Cover[i])) & (!NumericVector::is_na(H[i]))) {
      areaind[i] = aShrubAreaSP[SP[i]]*pow(H[i],2.0)/10000.0; 
    }
  }
  return(areaind);
}

/*
 * Cohort density in ind/ha
 */
// [[Rcpp::export("plant.Density")]]
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
  return(N);
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
// [[Rcpp::export("plant.Height")]]
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
  return(H);
}



/** 
 * Crown ratio, crown base height, crown length
 */
// [[Rcpp::export(".shrubCrownRatio")]]
NumericVector shrubCrownRatio(IntegerVector SP, DataFrame SpParams) {
  NumericVector crSP = SpParams["cr"];
  int nshrub = SP.size();
  NumericVector CR(nshrub);
  for(int i=0;i<nshrub;i++) {
    CR[i] = crSP[SP[i]];
  }
  return(CR);
}
double crownCompetitionFactor(IntegerVector SP, NumericVector N, NumericVector dbh, DataFrame SpParams) {
  NumericVector spAcw = SpParams["a_cw"];
  NumericVector spBcw = SpParams["b_cw"];
  int ntree = SP.size();
  double ccf = 0.0;
  for(int i=0;i<ntree;i++) {
    double cw = spAcw[SP[i]]*pow(dbh[i], spBcw[SP[i]]);
    ccf = ccf + (N[i]*PI*pow(cw/2.0,2.0)/100.0);
  }
  return(ccf);
}
NumericVector treeCrownRatio(IntegerVector SP, NumericVector N, NumericVector dbh, NumericVector H, DataFrame SpParams) {
  NumericVector BAL = largerTreeBasalArea(N, dbh);
  double ccf = crownCompetitionFactor(SP, N, dbh, SpParams);
  // Rcout<<ccf<<"\n";
  NumericVector spAcr = SpParams["a_cr"];
  NumericVector spB1cr = SpParams["b_1cr"];
  NumericVector spB2cr = SpParams["b_2cr"];
  NumericVector spB3cr = SpParams["b_3cr"];
  NumericVector spC1cr = SpParams["c_1cr"];
  NumericVector spC2cr = SpParams["c_2cr"];
  int ntree = SP.size();
  NumericVector treeCR(ntree);
  for(int i=0;i<ntree;i++) {
    double lm = spAcr[SP[i]]+ spB1cr[SP[i]]*(H[i]/(100.0*dbh[i]))+spB2cr[SP[i]]*(H[i]/100.0)+spB3cr[SP[i]]*pow(dbh[i],2.0)+spC1cr[SP[i]]*BAL[i]+spC2cr[SP[i]]*log(ccf);
    treeCR[i] = 1.0/(1.0+ exp(-1.0*lm));
  }
  return(treeCR);
}
// [[Rcpp::export("plant.CrownRatio")]]
NumericVector cohortCrownRatio(List x, DataFrame SpParams) {
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  DataFrame shrubData = Rcpp::as<Rcpp::DataFrame>(x["shrubData"]);
  NumericVector crSP = SpParams["cr"];
  int ntree = treeData.nrows();
  int nshrub = shrubData.nrows();
  NumericVector shrubSP = shrubData["Species"];  
  int numCohorts  = ntree+nshrub;
  NumericVector treeCR = treeCrownRatio(treeData["Species"],treeData["N"], treeData["DBH"], treeData["Height"], SpParams);
  NumericVector CR(numCohorts);
  for(int i=0;i<ntree;i++) {
    CR[i] = treeCR[i];
  }
  for(int i=0;i<nshrub;i++) {
    CR[ntree+i] = crSP[shrubSP[i]];
  }
  return(CR);
}
// [[Rcpp::export("plant.CrownBaseHeight")]]
NumericVector cohortCrownBaseHeight(List x, DataFrame SpParams) {
  NumericVector CR = cohortCrownRatio(x, SpParams);
  NumericVector H = cohortHeight(x);
  int numCohorts = H.size();
  NumericVector CBH(numCohorts);
  for(int i=0;i<numCohorts;i++) {
    CBH[i] = H[i]*(1.0-CR[i]);
  }
  return(CBH);
}
// [[Rcpp::export("plant.CrownLength")]]
NumericVector cohortCrownLength(List x, DataFrame SpParams) {
  NumericVector CR = cohortCrownRatio(x, SpParams);
  NumericVector H = cohortHeight(x);
  int numCohorts = H.size();
  NumericVector CL(numCohorts);
  for(int i=0;i<numCohorts;i++) {
    CL[i] = H[i]*CR[i];
  }
  return(CL);
}


/**
 * Foliar biomass (in kg/m2)
 */
// [[Rcpp::export(".treeFoliarBiomass")]]
NumericVector treeFoliarBiomass(IntegerVector SP, NumericVector N, NumericVector dbh, DataFrame SpParams, double gdd = NA_REAL){
  NumericVector afbtSP = SpParams["a_fbt"];
  NumericVector bfbtSP = SpParams["b_fbt"];
  NumericVector cfbtSP = SpParams["c_fbt"];
  NumericVector dfbtSP = SpParams["d_fbt"];
  NumericVector ltba = largerTreeBasalArea(N,dbh);
  int ncoh = N.size();
  NumericVector lb(ncoh);
  for(int i=0;i<ncoh;i++) {
    lb[i] = ((N[i]/10000)*afbtSP[SP[i]]*pow(dbh[i], bfbtSP[SP[i]])*exp(cfbtSP[SP[i]]*ltba[i])*pow(dbh[i], dfbtSP[SP[i]]*ltba[i]));
  }
  if(!NumericVector::is_na(gdd)) {
    NumericVector SgddSP = SpParams["Sgdd"];
    for(int i=0;i<ncoh;i++) {
      if(!NumericVector::is_na(SP[i])) lb[i] = lb[i]*leafDevelopmentStatus(SgddSP[SP[i]], gdd);
    }
  }
  
  return(lb);
}

// [[Rcpp::export(".shrubFoliarBiomass")]]
NumericVector shrubFoliarBiomass(IntegerVector SP, NumericVector Cover, NumericVector H, NumericVector CR, DataFrame SpParams, double gdd = NA_REAL){
  NumericVector aShrubAreaSP = SpParams["a_ash"];
  NumericVector aShrubFuelSP = SpParams["a_bsh"];
  NumericVector bShrubFuelSP = SpParams["b_bsh"];
  NumericVector SgddSP = SpParams["Sgdd"];
  NumericVector pDeadSP = SpParams["pDead"];
  NumericVector fTreeFuelSP = SpParams["r635"];
  int ncoh = SP.size();
  double W = 0.0; //Fine fuel
  NumericVector fb(ncoh);
  NumericVector areaind = shrubIndividualArea(SP,Cover,H,SpParams); 
  double volind = NA_REAL,weightkgind = NA_REAL;
  for(int i=0;i<ncoh;i++) {
    if((!NumericVector::is_na(SP[i])) & (!NumericVector::is_na(Cover[i])) & (!NumericVector::is_na(H[i]))) {
      volind = areaind[i]*((H[i]-(H[i]*(1.0-CR[i])))/100.0); //Crown phytovolume of an individual (in m3)
      weightkgind = aShrubFuelSP[SP[i]]*pow(volind,bShrubFuelSP[SP[i]]); //Fuel (in kg) of an individual (includes dead fuels)
      weightkgind = weightkgind - (weightkgind*pDeadSP[SP[i]]); //Removes dead fuels
      if(areaind[i]>0.0) {
        // multiply by 'number of individuals' per m2 
        W = weightkgind*(Cover[i]/(100*areaind[i]));  //Fine fuel (kg/m2)
        fb[i] = W/fTreeFuelSP[SP[i]]; //Foliar biomass (kg/m2)
        // Rcout<<Cover[i]<<" "<<(Cover[i]/(100*areaind))<<" "<< W<< " "<< fb[i]<<"\n";
        if(!NumericVector::is_na(gdd)) { //Apply phenology correction to foliar fuels
          fb[i] = fb[i]*leafDevelopmentStatus(SgddSP[SP[i]], gdd); 
        } 
      }
    }
    else fb[i] = NA_REAL;
  }
  return(fb);
}
// [[Rcpp::export("plant.FoliarBiomass")]]
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
  return(FB);
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
// [[Rcpp::export("plant.Cover")]]
NumericVector cohortCover(List x) {
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  DataFrame shrubData = Rcpp::as<Rcpp::DataFrame>(x["shrubData"]);
  NumericVector cover = shrubData["Cover"];
  NumericVector vol(treeData.nrows()+shrubData.nrows(), NA_REAL);
  for(int i=0;i<cover.size();i++) {
    vol[i+treeData.nrows()] = cover[i];
  }
  return(vol);
}


/**
 *  Shrub phytovolume (in m3/m2)
 */
// [[Rcpp::export(".shrubCrownPhytovolume")]]
NumericVector shrubCrownPhytovolume(IntegerVector SP, NumericVector Cover, NumericVector H, NumericVector CR, DataFrame SpParams){
  NumericVector aShrubAreaSP = SpParams["a_ash"];
  int ncoh = Cover.size();
  NumericVector vol(ncoh);
  for(int i=0;i<ncoh;i++) {
    if((!NumericVector::is_na(Cover[i]))& (!NumericVector::is_na(H[i]))) {
      double areaind = aShrubAreaSP[SP[i]]*pow(H[i],2.0)/10000.0; //area of an individual (in m2)
      double volind = areaind*((H[i]-(H[i]*(1.0-CR[i])))/100.0); //Crown phytovolume of an individual (in m3)
      // Rcout <<areaind<<" "<< volind<<"\n";
      vol[i] = volind * (Cover[i]/(100*areaind));
    } else {
      vol[i] = NA_REAL;
    }
  }
  return(vol);
}
// [[Rcpp::export("plant.Phytovolume")]]
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
  return(vol);
}


/**
 * Fine fuel loading (in kg/m2)
 */
// [[Rcpp::export(".treeFuel")]]
NumericVector treeFuel(IntegerVector SP, NumericVector N, NumericVector dbh, DataFrame SpParams, double gdd = NA_REAL, bool includeDead = true){
  NumericVector fb = treeFoliarBiomass(SP, N, dbh, SpParams, NA_REAL); //Do not include phenology (to have correct estimates of branch biomass)
  NumericVector SgddSP = SpParams["Sgdd"];
  NumericVector fTreeFuelSP = SpParams["r635"];
  NumericVector pDeadSP = SpParams["pDead"];
  int ncoh = N.size();
  double ftf = 0.0, btf = 0.0;
  NumericVector fuel(ncoh);
  for(int i=0;i<ncoh;i++) {
    if((!NumericVector::is_na(SP[i])) & (!NumericVector::is_na(dbh[i]))& (!NumericVector::is_na(N[i]))) {
      ftf = fb[i]; //Foliar biomass (kg per m2)
      btf = ftf*(fTreeFuelSP[SP[i]]-1.0); // Small branch fuels (proportion of foliar fuels)
      if(!NumericVector::is_na(gdd)) { //Apply phenology correction to foliar fuels
        ftf = ftf*leafDevelopmentStatus(SgddSP[SP[i]], gdd); 
      } 
      fuel[i] =  ftf + btf; //Tree fuel (kg per m2) is sum of both fuels
      if(includeDead) fuel[i] = fuel[i] + fuel[i]*pDeadSP[SP[i]]; //If required add fine dead fuels (proportion of live fuels)
    }
    else fuel[i] = NA_REAL;
  }
  return(fuel);
}
// [[Rcpp::export(".shrubFuel")]]
NumericVector shrubFuel(IntegerVector SP, NumericVector Cover, NumericVector H, NumericVector CR, DataFrame SpParams, double gdd = NA_REAL, bool includeDead = true){
  NumericVector aShrubAreaSP = SpParams["a_ash"];
  NumericVector aShrubFuelSP = SpParams["a_bsh"];
  NumericVector bShrubFuelSP = SpParams["b_bsh"];
  NumericVector fTreeFuelSP = SpParams["r635"];
  NumericVector pDeadSP = SpParams["pDead"];
  int ncoh = SP.size();
  double areaind = NA_REAL, volind = NA_REAL, weightkgind = NA_REAL;
  //W in kg/m2. Fine fuel, does not include phenology 
  NumericVector W(ncoh);
  for(int i=0;i<ncoh;i++) {
    if((!NumericVector::is_na(SP[i])) & (!NumericVector::is_na(Cover[i])) & (!NumericVector::is_na(H[i]))) {
      areaind = aShrubAreaSP[SP[i]]*pow(H[i],2.0)/10000.0; //area of an individual (in m2)
      volind = areaind*((H[i]-(H[i]*(1.0-CR[i])))/100.0); //Crown phytovolume of an individual (in m3)
      weightkgind = aShrubFuelSP[SP[i]]*pow(volind,bShrubFuelSP[SP[i]]); //Dry weight (in kg) of an individual
      if(!includeDead) weightkgind = weightkgind - (weightkgind*pDeadSP[SP[i]]); //Remove dead fuels if asked
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
    NumericVector SgddSP = SpParams["Sgdd"];
    for(int i=0;i<ncoh;i++) {
      if(!NumericVector::is_na(SP[i])) {
        fsf = W[i]/fTreeFuelSP[SP[i]]; //foliar biomass
        bsf = W[i] - fsf; //branch biomass
        W[i] = bsf + fsf*leafDevelopmentStatus(SgddSP[SP[i]], gdd); 
      }
    }
  }
  return(W);
}
// [[Rcpp::export("plant.Fuel")]]
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
  return(fuel);
}



/**
 * Cohort equilibrium leaf litter (in kg/m2)
 */
// [[Rcpp::export("plant.EquilibriumLeafLitter")]]
NumericVector cohortEquilibriumLeafLitter(List x, DataFrame SpParams, double AET = 800) {
  NumericVector fb = cohortFoliarBiomass(x, SpParams);
  NumericVector ld = cohortParameter(x, SpParams, "LeafDuration");
  NumericVector lignin = cohortParameter(x, SpParams, "LigninPercent");
  int ncoh = fb.size();
  NumericVector eqli(ncoh);
  double ki = 0.0;
  for(int i=0;i<ncoh;i++) {
    ki = (-0.5365+0.00241*AET) - (-0.01586+0.000056*AET)*lignin[i];//Meentemeyer (1978)
    // Rcout<<ki<<"\n";
    eqli[i] = fb[i]/(ld[i]*ki);
  }
  return(eqli);
}

/**
 * Cohort equilibrium small branch (6.35mm) litter (in kg/m2)
 */
// [[Rcpp::export("plant.EquilibriumSmallBranchLitter")]]
NumericVector cohortEquilibriumSmallBranchLitter(List x, DataFrame SpParams, double smallBranchDecompositionRate = 0.81) {
  NumericVector fu = cohortFuel(x, SpParams);
  NumericVector fb = cohortFoliarBiomass(x, SpParams);
  NumericVector ld = cohortParameter(x, SpParams, "LeafDuration");
  int ncoh = fb.size();
  NumericVector eqli(ncoh);
  for(int i=0;i<ncoh;i++) {
    eqli[i] = (fu[i]-fb[i])/((ld[i]*2.0)*smallBranchDecompositionRate);
  }
  return(eqli);
}


/**
 *  Leaf Area Index (LAI)
 */

// [[Rcpp::export(".treeLAI")]]
NumericVector treeLAI(IntegerVector SP, NumericVector N, NumericVector dbh, DataFrame SpParams, double gdd = NA_REAL){
  NumericVector SLASP = SpParams["SLA"]; // m2/kg (=mg/mm2)
  NumericVector lb = treeFoliarBiomass(SP, N, dbh, SpParams, gdd); //kg per m2
  int ncoh = N.size();
  NumericVector lai(ncoh);
  for(int i=0;i<ncoh;i++) {
     lai[i] = (SLASP[SP[i]]*lb[i]); 
  }
  return(lai);
}
// [[Rcpp::export(".shrubLAI")]]
NumericVector shrubLAI(IntegerVector SP, NumericVector Cover, NumericVector H, DataFrame SpParams, double gdd = NA_REAL){
  NumericVector SLASP = SpParams["SLA"]; // m2/kg (=mg/mm2)
  NumericVector CR = shrubCrownRatio(SP, SpParams);
  NumericVector lb = shrubFoliarBiomass(SP, Cover, H, CR, SpParams, gdd); //kg per m2
  int ncoh = SP.size();
  NumericVector lai(ncoh);
  for(int i=0;i<ncoh;i++) {
    lai[i] = (SLASP[SP[i]]*lb[i]); 
  }
  return(lai);
}
// [[Rcpp::export("plant.LAI")]]
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
  return(lai);
}

NumericMatrix LAIdistribution(NumericVector z, NumericVector LAI, NumericVector H, NumericVector CR) {
  int nh = z.size();
  int ncoh = LAI.size();
  double h1, h2;
  NumericMatrix LAIdist(nh-1, ncoh);
  for(int ci=0;ci<ncoh;ci++) {
    double cbh = H[ci]*(1.0-CR[ci]);
    double laih = LAI[ci]/(H[ci]-cbh);
    for(int hi=0;hi<(nh-1);hi++) {
      h1 = std::max(z[hi],cbh);
      h2 = std::min(z[hi+1],H[ci]);
      LAIdist(hi,ci) =laih*std::max(0.0, h2 - h1);
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

  return(LAIdistribution(z, LAI, H, CR));
}
NumericVector LAIprofile(NumericVector z, NumericVector LAI, NumericVector H, NumericVector CR) {
  int nh = z.size();
  int ncoh = LAI.size();
  double h1, h2;
  NumericVector LAIprof(nh-1);
  for(int ci=0;ci<ncoh;ci++) {
    double cbh = H[ci]*(1.0-CR[ci]);
    double laih = LAI[ci]/(H[ci]-cbh);
    for(int hi=0;hi<(nh-1);hi++) {
      h1 = std::max(z[hi],cbh);
      h2 = std::min(z[hi+1],H[ci]);
      LAIprof[hi] +=laih*std::max(0.0, h2 - h1);
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
  return(LAIprofile(z, LAI, H, CR));
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
  NumericVector LAI_live = cohortLAI(x, SpParams, gdd);
  IntegerVector SP(ntree+nshrub);
  NumericVector H(ntree+nshrub);
  
  IntegerVector treeSP = Rcpp::as<Rcpp::IntegerVector>(treeData["Species"]);
  NumericVector treeH = treeData["Height"];
  NumericVector treeDBH = treeData["DBH"];
  IntegerVector shrubSP = Rcpp::as<Rcpp::IntegerVector>(shrubData["Species"]);
  NumericVector shrubH = shrubData["Height"];  
  
  NumericVector N = cohortDensity(x, SpParams);
    
  NumericVector LAI_dead(ntree+nshrub);
  NumericVector DBH(ntree+nshrub);
  
  for(int i=0;i<ntree;i++) {
    SP[i] = treeSP[i];
    H[i] = treeH[i];
    DBH[i] = treeDBH[i];
    LAI_dead[i] = 0.0;
  }
  for(int i=0;i<nshrub;i++) {
    SP[ntree+i] = shrubSP[i];
    H[ntree+i] = shrubH[i];
    DBH[ntree+i] = NA_REAL;
  }
  return(DataFrame::create(_["SP"]=SP, _["N"] = N,  _["DBH"] = DBH, _["H"]=H, _["CR"] = CR, 
                           _["LAI_live"]=LAI_live, _["LAI_dead"] = LAI_dead));
  
}