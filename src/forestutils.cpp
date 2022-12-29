#define STRICT_R_HEADERS
#include <numeric>
#include <Rcpp.h>
#include <string.h>
#include <stdio.h>
#include "paramutils.h"
#include "phenology.h"
#include "root.h"
#include "soil.h"
#include "incgamma.h"
using namespace Rcpp;


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
        if(!NumericVector::is_na(x[i])) sba[j] +=x[i];
      }
    }
  }
  sba.attr("names") = speciesCharacterParameter(uniqueSp, SpParams, "Name");
  return(sba);
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
NumericVector largerTreeBasalArea(NumericVector N, NumericVector dbh, double self_include_prop = 0.5) {
  int ncoh = N.size();
  NumericVector BA = treeBasalArea(N, dbh); 
  NumericVector ltBA(ncoh);
  for(int i=0;i<ncoh;i++) {
    ltBA[i] = 0.0;
    for(int j=0;j<ncoh;j++) {
      if(i==j) ltBA[i] += (BA[j]*self_include_prop); //add half of its own basal area
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




NumericVector dbhClassBasalArea(List x, NumericVector DBHbreaks) {
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  NumericVector treeDBH = treeData["DBH"];
  NumericVector tba = treeBasalArea(treeData["N"], treeDBH);
  int nclasses = DBHbreaks.size()-1;
  int ntree = treeData.nrows();
  NumericVector dcBA(nclasses);
  for(int i=0;i<ntree;i++) {
    for(int j=0;j<nclasses;j++) if((treeDBH[i]>=DBHbreaks[j]) && (treeDBH[i]<DBHbreaks[j+1])) dcBA[j] += tba[i]; 
  }
  return(dcBA);
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
    for(int j=0;j<nclasses;j++) if((treeDBH[i]>=DBHbreaks[j]) && (treeDBH[i]<DBHbreaks[j+1])) dcDens[j] += treeDens[i]; 
  }
  return(dcDens);
}

//area of an individual (in m2)
NumericVector shrubIndividualAreaMED(IntegerVector SP, NumericVector Cover, NumericVector H, DataFrame SpParams){
  NumericVector aShrubArea = speciesNumericParameterWithImputation(SP,SpParams, "a_ash",true);
  NumericVector bShrubArea = speciesNumericParameterWithImputation(SP,SpParams, "b_ash",true);
  int ncoh = SP.size();
  NumericVector areaind(ncoh);
  for(int i=0;i<ncoh;i++) {
    if((!NumericVector::is_na(Cover[i])) && (!NumericVector::is_na(H[i]))) {
      areaind[i] = aShrubArea[i]*pow(H[i],bShrubArea[i])/10000.0; 
    }
  }
  return(areaind);
}

NumericVector shrubIndividualAreaUS(NumericVector H, NumericVector SingleShrubCrownArea){
  int ncoh = SingleShrubCrownArea.size();
  NumericVector areaind(ncoh);
  for(int i=0;i<ncoh;i++) {
    if((!NumericVector::is_na(SingleShrubCrownArea[i])) && (!NumericVector::is_na(H[i]))) {
      areaind[i] = SingleShrubCrownArea[i];   //This is the new equation with single shrub crown area modelled from RVS directly. Note RVS PCH (Projected area, crown, horizontal surface) is in cm^2, but in medfate input file, I have converted it to m2
    }
  }
  return(areaind);
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


/** 
 * Crown ratio, crown base height, crown length
 */
// [[Rcpp::export(".shrubCrownRatio")]]
NumericVector shrubCrownRatio(IntegerVector SP, DataFrame SpParams) {
  return(speciesNumericParameterWithImputation(SP, SpParams, "cr", true));
}

double crownCompetitionFactorMED(NumericVector N, NumericVector dbh, NumericVector Acw, NumericVector Bcw) {
  int ntree = N.size();
  double ccf = 0.0;
  for(int i=0;i<ntree;i++) {
    if(!NumericVector::is_na(dbh[i])) {
      double cw = Acw[i]*pow(dbh[i], Bcw[i]);
      ccf = ccf + (N[i]*M_PI*pow(cw/2.0,2.0)/100.0);
    }
  }
  return(ccf);
}
double crownCompetitionFactorMED(IntegerVector SP, NumericVector N, NumericVector dbh, DataFrame SpParams) {
  NumericVector Acw = speciesNumericParameterWithImputation(SP, SpParams, "a_cw",true);
  NumericVector Bcw = speciesNumericParameterWithImputation(SP, SpParams, "b_cw",true);
  return(crownCompetitionFactorMED(N,dbh,Acw,Bcw));
}

double crownCompetitionFactorUS(NumericVector N, NumericVector dbh, NumericVector CrownWidth) {
  int ntree = N.size();
  double ccf = 0.0;
  for(int i=0;i<ntree;i++) {
    if(!NumericVector::is_na(dbh[i])) {
      double cw = CrownWidth[i]; //Shengli: This is tree crown width. We will get the info from treedata, so I modify the sentence here. Note Crown width (unit in meter) that a tree of cohort i would have in open-ground conditions
      ccf = ccf + (N[i]*M_PI*pow(cw/2.0,2.0)/100.0);  //I do not understand why there is 100.0 in this sentence. I believe the cw unit in medfate in meter (see I need to ask what is the cw unit in medfate.)
    }
  }
  return(ccf);
}


NumericVector treeCrownRatioMED(NumericVector N, NumericVector dbh, NumericVector H, 
                                NumericVector Acw, NumericVector Bcw,
                                NumericVector Acr, NumericVector B1cr, NumericVector B2cr, NumericVector B3cr,
                                NumericVector C1cr, NumericVector C2cr) {
  NumericVector BAL = largerTreeBasalArea(N, dbh);
  double ccf = crownCompetitionFactorMED(N, dbh, Acw, Bcw);
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

NumericVector treeCrownRatioMED(IntegerVector SP, NumericVector N, NumericVector dbh, NumericVector H, DataFrame SpParams) {
  NumericVector BAL = largerTreeBasalArea(N, dbh);
  double ccf = crownCompetitionFactorMED(SP, N, dbh, SpParams);
  // Rcout<<ccf<<"\n";
  NumericVector Acr = speciesNumericParameterWithImputation(SP, SpParams, "a_cr",true);
  NumericVector B1cr = speciesNumericParameterWithImputation(SP, SpParams, "b_1cr",true);
  NumericVector B2cr = speciesNumericParameterWithImputation(SP, SpParams, "b_2cr",true);
  NumericVector B3cr = speciesNumericParameterWithImputation(SP, SpParams, "b_3cr",true);
  NumericVector C1cr = speciesNumericParameterWithImputation(SP, SpParams, "c_1cr",true);
  NumericVector C2cr = speciesNumericParameterWithImputation(SP, SpParams, "c_2cr",true);
  int ntree = SP.size();
  NumericVector treeCR(ntree);
  for(int i=0;i<ntree;i++) {
    double lm = Acr[i]+ B1cr[i]*(H[i]/(100.0*dbh[i]))+B2cr[i]*(H[i]/100.0)+B3cr[i]*pow(dbh[i],2.0)+C1cr[i]*BAL[i]+C2cr[i]*log(ccf);
    treeCR[i] = 1.0/(1.0+ exp(-1.0*lm));
  }
  return(treeCR);
}



/**
 * Foliar biomass (in kg/m2)
 */
NumericVector treeFoliarBiomassMED(IntegerVector SP, NumericVector N, NumericVector dbh, DataFrame SpParams, double gdd = NA_REAL){
  NumericVector afbt = speciesNumericParameterWithImputation(SP, SpParams, "a_fbt",true);
  NumericVector bfbt = speciesNumericParameterWithImputation(SP, SpParams, "b_fbt",true);
  NumericVector cfbt = speciesNumericParameterWithImputation(SP, SpParams, "c_fbt",true);
  NumericVector dfbt = speciesNumericParameterWithImputation(SP, SpParams, "d_fbt",true);
  NumericVector ltba = largerTreeBasalArea(N,dbh);
  int ncoh = N.size();
  NumericVector lb(ncoh);
  for(int i=0;i<ncoh;i++) {
    lb[i] = ((N[i]/10000)*afbt[i]*pow(dbh[i], bfbt[i])*exp(cfbt[i]*ltba[i])*pow(dbh[i], dfbt[i]*ltba[i]));
  }
  if(!NumericVector::is_na(gdd)) {
    NumericVector Sgdd = speciesNumericParameterWithImputation(SP, SpParams, "Sgdd");
    for(int i=0;i<ncoh;i++) {
      if(!NumericVector::is_na(SP[i])) lb[i] = lb[i]*leafDevelopmentStatus(Sgdd[i], gdd);
    }
  }
  
  return(lb);
}

NumericVector treeFoliarBiomassUS(IntegerVector SP, NumericVector N, NumericVector FoliageBiomass, DataFrame SpParams, double gdd = NA_REAL){
  int ncoh = N.size();
  NumericVector lb(ncoh);
  for(int i=0;i<ncoh;i++) {
    //note FoliageBiomass is live tree foliage biomass in KG of the SINGLE tree. N[i] is TreesPerHa, and 10000 is 1 ha (10000 m2). Therefore, lb[i] is kg/m2
    lb[i] = N[i] * FoliageBiomass[i] / 10000.0;
  }
  if(!NumericVector::is_na(gdd)) {
    NumericVector Sgdd = speciesNumericParameterWithImputation(SP, SpParams, "Sgdd");
    for(int i=0;i<ncoh;i++) {
      if(!NumericVector::is_na(SP[i])) lb[i] = lb[i]*leafDevelopmentStatus(Sgdd[i], gdd);
    }
  }
  
  return(lb);
}

NumericVector shrubFoliarBiomassMED(IntegerVector SP, NumericVector Cover, NumericVector H, NumericVector CR, 
                                    DataFrame SpParams, double gdd = NA_REAL){
  NumericVector aShrubFuel = speciesNumericParameterWithImputation(SP, SpParams, "a_bsh",true);
  NumericVector bShrubFuel = speciesNumericParameterWithImputation(SP, SpParams, "b_bsh",true);
  NumericVector Sgdd = speciesNumericParameterWithImputation(SP, SpParams, "Sgdd");
  NumericVector pDead = speciesNumericParameterWithImputation(SP, SpParams, "pDead");
  NumericVector fTreeFuel = speciesNumericParameterWithImputation(SP, SpParams, "r635");
  int ncoh = SP.size();
  double W = 0.0; //Fine fuel
  NumericVector fb(ncoh);
  NumericVector areaind = shrubIndividualAreaMED(SP,Cover,H,SpParams); 
  double volind = NA_REAL,weightkgind = NA_REAL;
  for(int i=0;i<ncoh;i++) {
    // Rcout<<i<<": "<< H[i]<<" "<<CR[i]<<" "<<aShrubFuel[i]<<" "<<bShrubFuel[i]<< " "<< pDead[i]<<" "<<fTreeFuel[i]<<" "<< areaind[i]<<".\n";
    if((!NumericVector::is_na(Cover[i])) && (!NumericVector::is_na(H[i]))) {
      volind = areaind[i]*((H[i]-(H[i]*(1.0-CR[i])))/100.0); //Crown phytovolume of an individual (in m3)
      weightkgind = aShrubFuel[i]*pow(volind,bShrubFuel[i]); //Fuel (in kg) of an individual (includes dead fuels)
      weightkgind = weightkgind - (weightkgind*pDead[i]); //Removes dead fuels
      // Rcout<<volind<< " "<<weightkgind;
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
NumericVector shrubFoliarBiomassUS(IntegerVector SP, NumericVector H, 
                                   NumericVector SingleShrubCrownArea, NumericVector FoliageBiomassPerUnitArea, 
                                   DataFrame SpParams, double gdd = NA_REAL){
  NumericVector Sgdd = speciesNumericParameterWithImputation(SP, SpParams, "Sgdd");
  int ncoh = SP.size();
  NumericVector fb(ncoh);
  NumericVector areaind = shrubIndividualAreaUS(H,SingleShrubCrownArea); //SingleShrubCrownArea added here
  for(int i=0;i<ncoh;i++) {
    if((!NumericVector::is_na(H[i]))) {
      if(areaind[i]>0.0) {
        fb[i] = FoliageBiomassPerUnitArea[i];  // This is from medfate input file, defined as "FoliageBiomassPerUnitArea" in kg/m2. Note "the number of shrub stems" represented by this record is already taken into account.
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

NumericVector treeCoverMED(IntegerVector SP, NumericVector N, NumericVector dbh, DataFrame SpParams){
  NumericVector Acw = speciesNumericParameterWithImputation(SP, SpParams, "a_cw",true);
  NumericVector Bcw = speciesNumericParameterWithImputation(SP, SpParams, "b_cw",true);
  int ncoh = N.size();
  NumericVector cov(ncoh);
  for(int i=0;i<ncoh;i++) {
    if(!NumericVector::is_na(dbh[i])) {
      double cw = Acw[i]*pow(dbh[i], Bcw[i]);
      cov[i] = std::min(100.0,(N[i]*M_PI*pow(cw/2.0,2.0)/100.0));
    }
  }
  return(cov);
}

NumericVector treeCoverUS(IntegerVector SP, NumericVector N, NumericVector CrownWidth, DataFrame SpParams){
  int ncoh = N.size();
  NumericVector cov(ncoh);
  for(int i=0;i<ncoh;i++) {
    if(!NumericVector::is_na(CrownWidth[i])) {
      double cw = CrownWidth[i]; //Shengli: This is tree crown width. We will get the info from treedata, so I modify the sentence here. Note Crown width (unit in meter) that a tree of cohort i would have in open-ground conditions
      cov[i] = std::min(100.0,(N[i]*M_PI*pow(cw/2.0,2.0)/100.0));  //I do not understand why there is 100.0 in this sentence. I believe the cw unit in medfate in meter (see I need to ask what is the cw unit in medfate.)
    }
  }
  return(cov);
}

/**
 *  Shrub phytovolume (in m3/m2)
 */
// [[Rcpp::export(".shrubPhytovolume")]]
NumericVector shrubPhytovolume(IntegerVector SP, NumericVector Cover, NumericVector H, NumericVector CR, DataFrame SpParams){
  int ncoh = Cover.size();
  NumericVector vol(ncoh);
  NumericVector areaind = shrubIndividualAreaMED(SP,Cover,H,SpParams); //area of an individual (in m2)
  for(int i=0;i<ncoh;i++) {
    if((!NumericVector::is_na(Cover[i])) && (!NumericVector::is_na(H[i]))) {
      double volind = areaind[i]*(H[i]/100.0); //Phytovolume of an individual (in m3)
      // Rcout <<areaind<<" "<< volind<<"\n";
      vol[i] = volind * (Cover[i]/(100.0*areaind[i]));
    } else {
      vol[i] = NA_REAL;
    }
  }
  return(vol);
}


/**
 * Fine fuel loading (in kg/m2)
 */
NumericVector treeFuelMED(IntegerVector SP, NumericVector N, NumericVector dbh, DataFrame SpParams, double gdd = NA_REAL, bool includeDead = true){
  NumericVector fb = treeFoliarBiomassMED(SP, N, dbh, SpParams, NA_REAL); //Do not include phenology (to have correct estimates of branch biomass)
  NumericVector Sgdd = speciesNumericParameterWithImputation(SP, SpParams, "Sgdd");
  NumericVector fTreeFuel = speciesNumericParameterWithImputation(SP, SpParams, "r635", true);
  NumericVector pDead = speciesNumericParameterWithImputation(SP, SpParams, "pDead");
  int ncoh = N.size();
  double ftf = 0.0, btf = 0.0;
  NumericVector fuel(ncoh);
  for(int i=0;i<ncoh;i++) {
    if((!NumericVector::is_na(dbh[i])) && (!NumericVector::is_na(N[i]))) {
      ftf = fb[i]; //Foliar biomass (kg per m2)
      btf = ftf*(fTreeFuel[i]-1.0); // Small branch fuels (proportion of foliar fuels)
      if(!NumericVector::is_na(gdd)) { //Apply phenology correction to foliar fuels
        ftf = ftf*leafDevelopmentStatus(Sgdd[i], gdd); 
      } 
      fuel[i] =  ftf + btf; //Tree fuel (kg per m2) is sum of both fuels
      if(includeDead) {
        fuel[i] = fuel[i] + fuel[i]*pDead[i]; //If required add fine dead fuels (proportion of live fuels)
      }
    }
    else fuel[i] = NA_REAL;
  }
  return(fuel);
}

NumericVector treeFuelUS(IntegerVector SP, NumericVector N, NumericVector FoliageBiomass, DataFrame SpParams, double gdd = NA_REAL, bool includeDead = true){
  NumericVector fb = treeFoliarBiomassUS(SP, N, FoliageBiomass, SpParams, NA_REAL); //Do not include phenology (to have correct estimates of branch biomass)
  NumericVector Sgdd = speciesNumericParameterWithImputation(SP, SpParams, "Sgdd");
  NumericVector fTreeFuel = speciesNumericParameterWithImputation(SP, SpParams, "r635", true);
  NumericVector pDead = speciesNumericParameterWithImputation(SP, SpParams, "pDead");
  int ncoh = N.size();
  double ftf = 0.0, btf = 0.0;
  NumericVector fuel(ncoh);
  for(int i=0;i<ncoh;i++) {
    if((!NumericVector::is_na(fb[i]))) {
      ftf = fb[i]; //Foliar biomass (kg per m2)
      btf = ftf*(fTreeFuel[i]-1.0); // Small branch fuels (proportion of foliar fuels)
      if(!NumericVector::is_na(gdd)) { //Apply phenology correction to foliar fuels
        ftf = ftf*leafDevelopmentStatus(Sgdd[i], gdd); 
      } 
      fuel[i] =  ftf + btf; //Tree fuel (kg per m2) is sum of both fuels
      if(includeDead) {
        fuel[i] = fuel[i] + fuel[i]*pDead[i]; //If required add fine dead fuels (proportion of live fuels) 
      }
    }
    else fuel[i] = NA_REAL;
  }
  return(fuel);
}

NumericVector shrubFuelMED(IntegerVector SP, NumericVector Cover, NumericVector H, NumericVector CR, DataFrame SpParams, double gdd = NA_REAL, bool includeDead = true){
  NumericVector aShrubFuel = speciesNumericParameterWithImputation(SP, SpParams, "a_bsh",true);
  NumericVector bShrubFuel = speciesNumericParameterWithImputation(SP, SpParams, "b_bsh",true);
  NumericVector pDead = speciesNumericParameterWithImputation(SP, SpParams, "pDead");
  NumericVector fTreeFuel = speciesNumericParameterWithImputation(SP, SpParams, "r635", true);

  int ncoh = SP.size();
  double volind = NA_REAL, weightkgind = NA_REAL;
  //W in kg/m2. Fine fuel, does not include phenology 
  NumericVector W(ncoh);
  NumericVector areaind = shrubIndividualAreaMED(SP,Cover,H,SpParams); //area of an individual (in m2)
  for(int i=0;i<ncoh;i++) {
    if((!NumericVector::is_na(Cover[i])) && (!NumericVector::is_na(H[i]))) {
      volind = areaind[i]*(H[i]/100.0); //Phytovolume of an individual (in m3)
      weightkgind = aShrubFuel[i]*pow(volind,bShrubFuel[i]); //Dry weight (in kg) of an individual
      if(!includeDead) {
        weightkgind = weightkgind - (weightkgind*pDead[i]); //Remove dead fuels if asked 
      }
      if(areaind[i]>0.0) {
        // multiply by 'number of individuals' per m2 
        W[i] = weightkgind*(Cover[i]/(100.0*areaind[i])); 
      }
    }
    else W[i] = NA_REAL;
  }
  //Remove (if necessary), the weight due to leaves that are not there
   if(!NumericVector::is_na(gdd)) {
    double fsf = 0.0, bsf = 0.0;
    NumericVector Sgdd = speciesNumericParameterWithImputation(SP, SpParams, "Sgdd");
    for(int i=0;i<ncoh;i++) {
      fsf = W[i]/fTreeFuel[i]; //foliar biomass
      bsf = W[i] - fsf; //branch biomass
      W[i] = bsf + fsf*leafDevelopmentStatus(Sgdd[i], gdd); 
    }
  }
  return(W);
}

NumericVector shrubFuelUS(IntegerVector SP, NumericVector H, 
                          NumericVector SingleShrubCrownArea, NumericVector FoliageBiomassPerUnitArea,
                          DataFrame SpParams, double gdd = NA_REAL, bool includeDead = true){
  NumericVector fTreeFuel = speciesNumericParameterWithImputation(SP, SpParams, "r635", true);
  
  int ncoh = SP.size();
  //W in kg/m2. Fine fuel, does not include phenology 
  NumericVector W(ncoh);
  for(int i=0;i<ncoh;i++) {
    if((!NumericVector::is_na(H[i]))) {
      W[i] = FoliageBiomassPerUnitArea[i]*fTreeFuel[i];
    }
    else W[i] = NA_REAL;
  }
  //Remove (if necessary), the weight due to leaves that are not there
  if(!NumericVector::is_na(gdd)) {
    double fsf = 0.0, bsf = 0.0;
    NumericVector Sgdd = speciesNumericParameterWithImputation(SP, SpParams, "Sgdd");
    for(int i=0;i<ncoh;i++) {
      fsf = W[i]/fTreeFuel[i]; //foliar biomass
      bsf = W[i] - fsf; //branch biomass
      W[i] = bsf + fsf*leafDevelopmentStatus(Sgdd[i], gdd); 
    }
  }
  return(W);
}

/**
 *  Leaf Area Index (LAI)
 */
NumericVector treeLAIMED(IntegerVector SP, NumericVector N, NumericVector dbh, DataFrame SpParams, double gdd = NA_REAL){
  NumericVector SLA = speciesNumericParameterWithImputation(SP, SpParams, "SLA", true); // m2/kg (=mg/mm2)
  NumericVector lb = treeFoliarBiomassMED(SP, N, dbh, SpParams, gdd); //kg per m2
  int ncoh = N.size();
  NumericVector lai(ncoh);
  for(int i=0;i<ncoh;i++) {
     lai[i] = (SLA[i]*lb[i]); 
  }
  return(lai);
}
NumericVector treeLAIUS(IntegerVector SP, NumericVector N, NumericVector FoliageBiomass, DataFrame SpParams, double gdd = NA_REAL){
  NumericVector SLA = speciesNumericParameterWithImputation(SP, SpParams, "SLA", true); // m2/kg (=mg/mm2)
  NumericVector lb = treeFoliarBiomassUS(SP, N, FoliageBiomass, SpParams, gdd); //kg per m2
  int ncoh = N.size();
  NumericVector lai(ncoh);
  for(int i=0;i<ncoh;i++) {
    lai[i] = (SLA[i]*lb[i]); 
  }
  return(lai);
}
NumericVector shrubLAIMED(IntegerVector SP, NumericVector Cover, NumericVector H, DataFrame SpParams, double gdd = NA_REAL){
  NumericVector SLA = speciesNumericParameterWithImputation(SP, SpParams, "SLA", true); // m2/kg (=mg/mm2)
  NumericVector CR = shrubCrownRatio(SP, SpParams);
  NumericVector lb = shrubFoliarBiomassMED(SP, Cover, H, CR, SpParams, gdd); //kg per m2
  int ncoh = SP.size();
  NumericVector lai(ncoh);
  for(int i=0;i<ncoh;i++) {
    lai[i] = (SLA[i]*lb[i]); 
  }
  return(lai);
}

NumericVector shrubLAIUS(IntegerVector SP, NumericVector H, 
                         NumericVector SingleShrubCrownArea, NumericVector FoliageBiomassPerUnitArea, 
                         DataFrame SpParams, double gdd = NA_REAL){
  NumericVector SLA = speciesNumericParameterWithImputation(SP, SpParams, "SLA", true); // m2/kg (=mg/mm2)
  NumericVector lb = shrubFoliarBiomassUS(SP, H, 
                                          SingleShrubCrownArea, FoliageBiomassPerUnitArea, 
                                          SpParams, gdd); //kg per m2   ////Shengli: THis need to be revised according to new function
  int ncoh = SP.size();
  NumericVector lai(ncoh);
  for(int i=0;i<ncoh;i++) {
    lai[i] = (SLA[i]*lb[i]); 
  }
  return(lai);
}


//' Plant description functions
//'
//' Functions to calculate attributes of plants in a \code{\link{forest}} object.
//' 
//' @param x An object of class \code{\link{forest}}.
//' @param SpParams A data frame with species parameters (see \code{\link{SpParamsMED}}).
//' @param parName A string with a parameter name.
//' @param mode Calculation mode, either "MED" or "US".
//' @param gdd Growth degree days (to account for leaf phenology effects).
//' @param AET Actual annual evapotranspiration (in mm).
//' @param smallBranchDecompositionRate Decomposition rate of small branches.
//' @param includeDead A flag to indicate that standing dead fuels (dead branches) are included.
//' @param treeOffset,shrubOffset Integers to offset cohort IDs.
//' @param fillMissing A boolean flag to try imputation on missing values.
//' @param self_proportion Proportion of the target cohort included in the assessment
//' 
//' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
//' 
//' @return
//' A vector with values for each plant of the input \code{\link{forest}} object:
//' \itemize{
//'   \item{\code{plant_basalArea}: Tree basal area (m2/ha).}
//'   \item{\code{plant_largerTreeBasalArea}: Basal area (m2/ha) of trees larger (in diameter) than the tree. Half of the trees of the same record are included.}
//'   \item{\code{plant_characterParameter}: The parameter values of each plant, as strings.}
//'   \item{\code{plant_cover}: Shrub cover (in percent).}
//'   \item{\code{plant_crownBaseHeight}: The height corresponding to the start of the crown (in cm).}
//'   \item{\code{plant_crownLength}: The difference between crown base height and total height (in cm).}
//'   \item{\code{plant_crownRatio}: The ratio between crown length and total height (between 0 and 1).}
//'   \item{\code{plant_density}: Plant density (ind/ha). Tree density is directly taken from the forest object, while the shrub density is estimated from cover and height by calculating the area of a single individual.}
//'   \item{\code{plant_equilibriumLeafLitter}: Litter biomass of leaves at equilibrium (in kg/m2).}
//'   \item{\code{plant_equilibriumSmallBranchLitter}: Litter biomass of small branches (< 6.35 mm diameter) at equilibrium (in kg/m2).}
//'   \item{\code{plant_foliarBiomass}: Standing biomass of leaves (in kg/m2).}
//'   \item{\code{plant_fuel}: Fine fuel load (in kg/m2).}
//'   \item{\code{plant_height}: Total height (in cm).}
//'   \item{\code{plant_ID}: Cohort coding for simulation functions (concatenation of 'T' (Trees) or 'S' (Shrub), cohort index and species index).}
//'   \item{\code{plant_LAI}: Leaf area index (m2/m2).}
//'   \item{\code{plant_individualArea}: Area (m2) occupied by a shrub individual.}
//'   \item{\code{plant_parameter}: The parameter values of each plant, as numeric.}
//'   \item{\code{plant_phytovolume}: Shrub phytovolume (m3/m2).}
//'   \item{\code{plant_species}: Species identity integer (indices start with 0).}
//'   \item{\code{plant_speciesName}: String with species taxonomic name (or a functional group).}
//' }
//' 
//' @seealso  \code{\link{spwb}}, \code{\link{forest}}, \code{\link{summary.forest}}
//'
//' @examples
//' #Default species parameterization
//' data(SpParamsMED)
//' 
//' #Load example plot
//' data(exampleforestMED)
//' 
//' #A short way to obtain total basal area
//' sum(plant_basalArea(exampleforestMED), na.rm=TRUE)
//' 
//' #The same forest level function for LAI
//' sum(plant_LAI(exampleforestMED, SpParamsMED))
//'   
//' #The same forest level function for fuel loading
//' sum(plant_fuel(exampleforestMED, SpParamsMED))
//'       
//' #Summary function for 'forest' objects can be also used
//' summary(exampleforestMED, SpParamsMED)
//' 
//' plant_speciesName(exampleforestMED, SpParamsMED)
//' 
//' plant_ID(exampleforestMED)
//'       
//' @name plant_values
// [[Rcpp::export("plant_ID")]]
CharacterVector cohortIDs(List x, int treeOffset = 0, int shrubOffset = 0) {
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  DataFrame shrubData = Rcpp::as<Rcpp::DataFrame>(x["shrubData"]);
  int ntree = treeData.nrows();
  int nshrub = shrubData.nrows();
  int numCohorts  = ntree+nshrub;
  IntegerVector treeSP = treeData["Species"];
  IntegerVector shrubSP = shrubData["Species"];  
  
  CharacterVector IDs(numCohorts);
  for(int i=0;i<ntree;i++) {
    String s("T");
    s += (i+treeOffset+1);
    s += "_";
    s += treeSP[i];
    IDs[i] = s;
  }
  for(int i=0;i<nshrub;i++) {
    String s("S");
    s += (i+shrubOffset+1);
    s += "_";
    s += shrubSP[i];
    IDs[ntree+i] =s;
  }
  return(IDs);
}

//' @rdname plant_values
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

//' @rdname plant_values
// [[Rcpp::export("plant_largerTreeBasalArea")]]
NumericVector cohortLargerTreeBasalArea(List x, double self_proportion = 0.5) {
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  DataFrame shrubData = Rcpp::as<Rcpp::DataFrame>(x["shrubData"]);
  NumericVector tba = largerTreeBasalArea(treeData["N"], treeData["DBH"], self_proportion);
  NumericVector ba(tba.size()+shrubData.nrows(), NA_REAL);
  for(int i=0;i<tba.size();i++) {
    ba[i] = tba[i];
  }
  ba.attr("names") = cohortIDs(x);
  return(ba);
}

//' @rdname plant_values
// [[Rcpp::export("plant_cover")]]
NumericVector cohortCover(List x, DataFrame SpParams, String mode = "MED") {
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  DataFrame shrubData = Rcpp::as<Rcpp::DataFrame>(x["shrubData"]);
  NumericVector cov(treeData.nrows()+shrubData.nrows(), NA_REAL);
  NumericVector tcover;
  if(mode=="MED") {
    tcover = treeCoverMED(treeData["Species"], treeData["N"], treeData["DBH"],
                          SpParams);
  } else {
    tcover = treeCoverUS(treeData["Species"], treeData["N"], treeData["CrownWidth"], 
                         SpParams);
  }
  for(int i=0;i<tcover.size();i++) {
    cov[i] = tcover[i];
  }
  NumericVector shcover = shrubData["Cover"];
  for(int i=0;i<shcover.size();i++) {
    cov[i+treeData.nrows()] = shcover[i];
  }
  cov.attr("names") = cohortIDs(x);
  return(cov);
}

//' @rdname plant_values
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

//' @rdname plant_values
// [[Rcpp::export("plant_speciesName")]]
CharacterVector cohortSpeciesName(List x, DataFrame SpParams) {
  CharacterVector sn = cohortCharacterParameter(x,SpParams, "Name");
  sn.attr("names") = cohortIDs(x);
  return(sn);
}

//' @rdname plant_values
// [[Rcpp::export("plant_density")]]
NumericVector cohortDensity(List x, DataFrame SpParams, String mode = "MED") {
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  DataFrame shrubData = Rcpp::as<Rcpp::DataFrame>(x["shrubData"]);
  int ntree = treeData.nrows();
  int nshrub = shrubData.nrows();
  NumericVector treeN = treeData["N"];
  IntegerVector shrubSP = shrubData["Species"];
  NumericVector shrubHeight = shrubData["Height"];
  int numCohorts  = ntree+nshrub;
  NumericVector N(numCohorts);
  for(int i=0;i<ntree;i++) {
    N[i] = treeN[i];
  }
  NumericVector shrubArea;
  if(mode=="MED") {
    NumericVector shrubCover = shrubData["Cover"];
    shrubArea = shrubIndividualAreaMED(shrubSP,shrubCover, shrubHeight, SpParams); 
    for(int i=0;i<nshrub;i++) {
      N[ntree+i] = 10000.0*(shrubCover[i]/(100.0*shrubArea[i]));
    }
  }
  else if(mode=="US") {
    NumericVector SingleShrubCrownArea = shrubData["SingleShrubCrownArea"];  //This is a newly added sentence to read "SingleShrubCrownArea"
    NumericVector NumberOfShrub = shrubData["NumberOfShrub"];  //This is a newly added sentence to read "NumberOfShrub" 
    shrubArea = shrubIndividualAreaUS(shrubHeight, SingleShrubCrownArea);
    for(int i=0;i<nshrub;i++) {
      N[ntree+i] = NumberOfShrub[i];  //RVS output unit is stems number per Acre (note an acre=4000 m2. In medfate input file, I have converted to stems number per ha. note 1 ha=10000 m2
    }
  }
  else stop("Wrong mode.");
  N.attr("names") = cohortIDs(x);
  return(N);
}

//' @rdname plant_values
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

//' @rdname plant_values
// [[Rcpp::export("plant_individualArea")]]
NumericVector cohortIndividualArea(List x, DataFrame SpParams, String mode = "MED"){
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  DataFrame shrubData = Rcpp::as<Rcpp::DataFrame>(x["shrubData"]);
  int ntree = treeData.nrows();
  int nshrub = shrubData.nrows();
  NumericVector treeN = treeData["N"];
  IntegerVector shrubSP = shrubData["Species"];
  NumericVector shrubHeight = shrubData["Height"];
  int numCohorts  = ntree+nshrub;
  NumericVector indArea(numCohorts, NA_REAL);
  NumericVector shrubArea;
  if(mode=="MED") {
    NumericVector shrubCover = shrubData["Cover"];
    shrubArea = shrubIndividualAreaMED(shrubSP,shrubCover, shrubHeight, SpParams); 
    for(int i=0;i<nshrub;i++) {
      indArea[ntree+i] = shrubArea[i];
    }
  }
  else if(mode=="US") {
    NumericVector SingleShrubCrownArea = shrubData["SingleShrubCrownArea"];  //This is a newly added sentence to read "SingleShrubCrownArea"
    shrubArea = shrubIndividualAreaUS(shrubHeight, SingleShrubCrownArea);
    for(int i=0;i<nshrub;i++) {
      indArea[ntree+i] = shrubArea[i];  //RVS output unit is stems number per Acre (note an acre=4000 m2. In medfate input file, I have converted to stems number per ha. note 1 ha=10000 m2
    }
  }
  indArea.attr("names") = cohortIDs(x);
  return(indArea);
}

//' @rdname plant_values
// [[Rcpp::export("plant_crownRatio")]]
NumericVector cohortCrownRatio(List x, DataFrame SpParams, String mode = "MED") {
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  DataFrame shrubData = Rcpp::as<Rcpp::DataFrame>(x["shrubData"]);
  int ntree = treeData.nrows();
  int nshrub = shrubData.nrows();
  IntegerVector shrubSP = shrubData["Species"];  
  NumericVector crSh = speciesNumericParameterWithImputation(shrubSP, SpParams, "cr",true);
  int numCohorts  = ntree+nshrub;
  NumericVector treeCR;
  if(mode=="MED") {
    treeCR = treeCrownRatioMED(treeData["Species"],treeData["N"], treeData["DBH"], treeData["Height"], SpParams); 
  } else if(mode=="US") {
    treeCR = treeData["CrownRatio"];
  }
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

//' @rdname plant_values
// [[Rcpp::export("plant_crownBaseHeight")]]
NumericVector cohortCrownBaseHeight(List x, DataFrame SpParams, String mode = "MED") {
  NumericVector CR = cohortCrownRatio(x, SpParams, mode);
  NumericVector H = cohortHeight(x);
  int numCohorts = H.size();
  NumericVector CBH(numCohorts);
  for(int i=0;i<numCohorts;i++) {
    CBH[i] = H[i]*(1.0-CR[i]);
  }
  CBH.attr("names") = cohortIDs(x);
  return(CBH);
}

//' @rdname plant_values
// [[Rcpp::export("plant_crownLength")]]
NumericVector cohortCrownLength(List x, DataFrame SpParams, String mode = "MED") {
  NumericVector CR = cohortCrownRatio(x, SpParams, mode);
  NumericVector H = cohortHeight(x);
  int numCohorts = H.size();
  NumericVector CL(numCohorts);
  for(int i=0;i<numCohorts;i++) {
    CL[i] = H[i]*CR[i];
  }
  CL.attr("names") = cohortIDs(x);
  return(CL);
}

//' @rdname plant_values
// [[Rcpp::export("plant_foliarBiomass")]]
NumericVector cohortFoliarBiomass(List x, DataFrame SpParams, double gdd = NA_REAL, String mode = "MED") {
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  DataFrame shrubData = Rcpp::as<Rcpp::DataFrame>(x["shrubData"]);
  NumericVector tFB;
  IntegerVector shSP = shrubData["Species"];
  NumericVector shCR = shrubCrownRatio(shSP, SpParams);
  NumericVector shFB;
  if(mode=="MED") {
    tFB = treeFoliarBiomassMED(treeData["Species"], treeData["N"], treeData["DBH"], SpParams, gdd);
    shFB= shrubFoliarBiomassMED(shSP, shrubData["Cover"], shrubData["Height"], shCR, 
                                SpParams, gdd);
  } else if(mode=="US") {
    tFB = treeFoliarBiomassUS(treeData["Species"], treeData["N"], treeData["FoliageBiomass"], SpParams, gdd);
    shFB= shrubFoliarBiomassUS(shSP, shrubData["Height"],
                               shrubData["SingleShrubCrownArea"], shrubData["FoliageBiomassPerUnitArea"],
                                                                           SpParams, gdd);
  }
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

//' @rdname plant_values
// [[Rcpp::export("plant_fuel")]]
NumericVector cohortFuel(List x, DataFrame SpParams, double gdd = NA_REAL, bool includeDead = true, String mode = "MED"){
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  DataFrame shrubData = Rcpp::as<Rcpp::DataFrame>(x["shrubData"]);
  NumericVector tFuel;
  IntegerVector shSP = shrubData["Species"];
  NumericVector shCR = shrubCrownRatio(shSP, SpParams);
  NumericVector shFuel;
  if(mode=="MED") {
    tFuel = treeFuelMED(treeData["Species"], treeData["N"], treeData["DBH"], SpParams, gdd, includeDead);
    shFuel = shrubFuelMED(shSP, shrubData["Cover"], shrubData["Height"], shCR, 
                          SpParams, gdd, includeDead);
  }
  else if(mode == "US") {
    tFuel = treeFuelUS(treeData["Species"], treeData["N"], treeData["FoliageBiomass"], SpParams, gdd, includeDead);
    shFuel = shrubFuelUS(shSP, shrubData["Height"], 
                         shrubData["SingleShrubCrownArea"], shrubData["FoliageBiomassPerUnitArea"],
                                                                     SpParams, gdd, includeDead);
    
  }
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

//' @rdname plant_values
// [[Rcpp::export("plant_equilibriumLeafLitter")]]
NumericVector cohortEquilibriumLeafLitter(List x, DataFrame SpParams, double AET = 800, String mode = "MED") {
  NumericVector fb = cohortFoliarBiomass(x, SpParams, NA_REAL, mode);
  NumericVector ld = cohortNumericParameter(x, SpParams, "LeafDuration");
  NumericVector lignin = cohortNumericParameterWithImputation(x, SpParams, "LigninPercent", true);
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

//' @rdname plant_values
// [[Rcpp::export("plant_equilibriumSmallBranchLitter")]]
NumericVector cohortEquilibriumSmallBranchLitter(List x, DataFrame SpParams, double smallBranchDecompositionRate = 0.81, String mode = "MED") {
  NumericVector fu = cohortFuel(x, SpParams,NA_REAL,true, mode);
  NumericVector fb = cohortFoliarBiomass(x, SpParams, NA_REAL, mode);
  NumericVector ld = cohortNumericParameter(x, SpParams, "LeafDuration");
  int ncoh = fb.size();
  NumericVector eqli(ncoh);
  for(int i=0;i<ncoh;i++) {
    eqli[i] = (fu[i]-fb[i])/((ld[i]*2.0)*smallBranchDecompositionRate);
  }
  eqli.attr("names") = cohortIDs(x);
  return(eqli);
}

//' @rdname plant_values
// [[Rcpp::export("plant_phytovolume")]]
NumericVector cohortPhytovolume(List x, DataFrame SpParams) {
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  DataFrame shrubData = Rcpp::as<Rcpp::DataFrame>(x["shrubData"]);
  IntegerVector SP = shrubData["Species"];
  NumericVector CR = shrubCrownRatio(SP, SpParams);
  NumericVector shvol = shrubPhytovolume(SP, shrubData["Cover"], shrubData["Height"], CR, SpParams);
  NumericVector vol(treeData.nrows()+shrubData.nrows(), NA_REAL);
  for(int i=0;i<shvol.size();i++) {
    vol[i+treeData.nrows()] = shvol[i];
  }
  vol.attr("names") = cohortIDs(x);
  return(vol);
}

//' @rdname plant_values
// [[Rcpp::export("plant_LAI")]]
NumericVector cohortLAI(List x, DataFrame SpParams, double gdd = NA_REAL, String mode = "MED"){
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  DataFrame shrubData = Rcpp::as<Rcpp::DataFrame>(x["shrubData"]);
  NumericVector tLAI;
  NumericVector shLAI;
  if(mode=="MED") {
    tLAI = treeLAIMED(treeData["Species"], treeData["N"], treeData["DBH"], SpParams, gdd);
    shLAI = shrubLAIMED(shrubData["Species"], shrubData["Cover"], shrubData["Height"], 
                        SpParams, gdd);
  } 
  else if(mode=="US") {
    tLAI = treeLAIUS(treeData["Species"], treeData["N"], treeData["FoliageBiomass"], SpParams, gdd);
    shLAI = shrubLAIUS(shrubData["Species"], shrubData["Height"], 
                       shrubData["SingleShrubCrownArea"], shrubData["FoliageBiomassPerUnitArea"],
                                                                   SpParams, gdd);
  } 
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





//' Species description functions
//'
//' Functions to calculate attributes of a \code{\link{forest}} object by species or to extract species parameters from a species parameter table (\code{\link{SpParamsMED}}).
//' 
//' @param x An object of class \code{\link{forest}}.
//' @param SpParams A data frame with species parameters (see \code{\link{SpParamsMED}}).
//' @param gdd Growth degree days (to account for leaf phenology effects).
//' @param includeDead A flag to indicate that standing dead fuels (dead branches) are included.
//' @param mode Calculation mode, either "MED" or "US".
//' @param SP An integer vector of species codes.
//' @param parName A string with a parameter name.
//' @param fillMissing A boolean flag to try imputation on missing values.
//' 
//' @return
//' A vector with values for each species in \code{SpParams}:
//' \itemize{
//'   \item{\code{species_basalArea}: Species basal area (m2/ha).}
//'   \item{\code{species_cover}: Shrub cover (in percent).}
//'   \item{\code{species_density}: Plant density (ind/ha). Tree density is directly taken from the forest object, while the shrub density is estimated from cover and height by calculating the area of a single individual.}
//'   \item{\code{species_foliarBiomass}: Standing biomass of leaves (in kg/m2).}
//'   \item{\code{species_fuel}: Fine fuel load (in kg/m2).}
//'   \item{\code{species_LAI}: Leaf area index (m2/m2).}
//'   \item{\code{species_phytovolume}: Shrub phytovolume (m3/m2).}
//'   \item{\code{species_parameter}: A numeric vector with the parameter values of each input species.}
//'   \item{\code{species_characterParameter}: A character vector with the parameter values of each input species.}
//' }
//' 
//' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
//' 
//' @seealso \code{\link{spwb}}, \code{\link{forest}}, \code{\link{plant_basalArea}}, \code{\link{summary.forest}}
//' 
//' @examples
//' # Default species parameterization
//' data(SpParamsMED)
//' 
//' # Load example plot
//' data(exampleforestMED)
//' 
//' # Species basal area in the forest plot
//' species_basalArea(exampleforestMED, SpParamsMED)
//'   
//' # Value of parameter "Psi_Extract" for species 157 (Pinus halepensis)
//' # and 176 (Quercus ilex)
//' species_parameter(c(157,176), SpParamsMED, "Psi_Extract")
//'     
//' @name species_values
// [[Rcpp::export("species_basalArea")]]
NumericVector speciesBasalArea(List x, DataFrame SpParams) {
  NumericVector cBA = cohortBasalArea(x);
  return(sumBySpecies(cBA, cohortSpecies(x), SpParams));
}

//' @rdname species_values
// [[Rcpp::export("species_cover")]]
NumericVector speciesCover(List x, DataFrame SpParams, String mode = "MED") {
  NumericVector cc = cohortCover(x, SpParams, mode);
  NumericVector sc = sumBySpecies(cc, cohortSpecies(x), SpParams);
  for(int i=0;i<sc.length();i++) sc[i] = std::min(100.0, sc[i]);
  return(sc);
}

//' @rdname species_values
// [[Rcpp::export("species_density")]]
NumericVector speciesDensity(List x, DataFrame SpParams, String mode = "MED") {
  NumericVector d = cohortDensity(x, SpParams, mode);
  return(sumBySpecies(d, cohortSpecies(x), SpParams));
}

//' @rdname species_values
// [[Rcpp::export("species_foliarBiomass")]]
NumericVector speciesFoliarBiomass(List x, DataFrame SpParams, double gdd = NA_REAL, String mode = "MED") {
  NumericVector fb = cohortFoliarBiomass(x, SpParams, gdd, mode);
  return(sumBySpecies(fb, cohortSpecies(x), SpParams));
}

//' @rdname species_values
// [[Rcpp::export("species_fuel")]]
NumericVector speciesFuel(List x, DataFrame SpParams, double gdd = NA_REAL, bool includeDead = true, String mode = "MED") {
  NumericVector cf = cohortFuel(x, SpParams, gdd, includeDead, mode);
  return(sumBySpecies(cf, cohortSpecies(x), SpParams));
}

//' @rdname species_values
// [[Rcpp::export("species_phytovolume")]]
NumericVector speciesPhytovolume(List x, DataFrame SpParams) {
  NumericVector cp = cohortPhytovolume(x, SpParams);
  return(sumBySpecies(cp, cohortSpecies(x), SpParams));
}

//' @rdname species_values
// [[Rcpp::export("species_LAI")]]
NumericVector speciesLAI(List x, DataFrame SpParams, double gdd = NA_REAL, String mode = "MED") {
  NumericVector cl = cohortLAI(x, SpParams, gdd, mode);
  return(sumBySpecies(cl, cohortSpecies(x), SpParams));
}

//' @rdname stand_values
// [[Rcpp::export("stand_basalArea")]]
double standBasalArea(List x, double minDBH = 7.5) {
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  NumericVector tba = treeBasalArea(treeData["N"], treeData["DBH"]);
  NumericVector treeDBH = treeData["DBH"];
  double ba = 0.0;
  for(int i=0;i<tba.size();i++) {
    if(treeDBH[i]>=minDBH) ba += tba[i];
  }
  return(ba);  
}

//' @rdname stand_values
// [[Rcpp::export("stand_foliarBiomass")]]
double standFoliarBiomass(List x, DataFrame SpParams, double gdd = NA_REAL, String mode = "MED") {
  NumericVector fb = cohortFoliarBiomass(x, SpParams, gdd, mode);
  double tfb= 0.0;
  for(int i=0;i<fb.size();i++){if(!NumericVector::is_na(fb[i])) tfb+=fb[i];}
  return(tfb);
}

//' @rdname stand_values
// [[Rcpp::export("stand_phytovolume")]]
double standPhytovolume(List x, DataFrame SpParams) {
  NumericVector cp = cohortPhytovolume(x, SpParams);
  double tp= 0.0;
  for(int i=0;i<cp.size();i++){if(!NumericVector::is_na(cp[i])) tp+=cp[i];}
  return(tp);
}

//' @rdname stand_values
// [[Rcpp::export("stand_fuel")]]
double standFuel(List x, DataFrame SpParams, double gdd = NA_REAL, bool includeDead = true, String mode = "MED") {
  NumericVector cf = cohortFuel(x, SpParams, gdd, includeDead, mode);
  double tf= 0.0;
  for(int i=0;i<cf.size();i++){if(!NumericVector::is_na(cf[i])) tf+=cf[i];}
  return(tf);
}



//' @rdname stand_values
// [[Rcpp::export("stand_LAI")]]
double standLAI(List x, DataFrame SpParams, double gdd = NA_REAL, String mode = "MED") {
  NumericVector cl = cohortLAI(x, SpParams, gdd, mode);
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
NumericMatrix LAIdistribution(NumericVector z, List x, DataFrame SpParams, double gdd = NA_REAL, String mode = "MED") {
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  DataFrame shrubData = Rcpp::as<Rcpp::DataFrame>(x["shrubData"]);
  int ntree = treeData.nrows();
  int nshrub = shrubData.nrows();
  
  IntegerVector treeSP = Rcpp::as<Rcpp::IntegerVector>(treeData["Species"]);
  NumericVector treeH = treeData["Height"];
  IntegerVector shrubSP = Rcpp::as<Rcpp::IntegerVector>(shrubData["Species"]);
  NumericVector shrubH = shrubData["Height"];  
  
  NumericVector LAI = cohortLAI(x, SpParams, gdd, mode);
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
NumericVector LAIprofile(NumericVector z, List x, DataFrame SpParams, double gdd = NA_REAL, String mode = "MED") {
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  DataFrame shrubData = Rcpp::as<Rcpp::DataFrame>(x["shrubData"]);
  int ntree = treeData.nrows();
  int nshrub = shrubData.nrows();
  
  IntegerVector treeSP = Rcpp::as<Rcpp::IntegerVector>(treeData["Species"]);
  NumericVector treeH = treeData["Height"];
  IntegerVector shrubSP = Rcpp::as<Rcpp::IntegerVector>(shrubData["Species"]);
  NumericVector shrubH = shrubData["Height"];  
  
  NumericVector CR = cohortCrownRatio(x, SpParams, mode);
  NumericVector LAI = cohortLAI(x, SpParams, gdd, mode);
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



//' @rdname modelInput
// [[Rcpp::export("forest2aboveground")]]
DataFrame forest2aboveground(List x, DataFrame SpParams, double gdd = NA_REAL, String mode = "MED") {
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  DataFrame shrubData = Rcpp::as<Rcpp::DataFrame>(x["shrubData"]);
  int ntree = treeData.nrows();
  int nshrub = shrubData.nrows();
  
  NumericVector CR = cohortCrownRatio(x, SpParams, mode);
  NumericVector LAI_live = cohortLAI(x, SpParams, NA_REAL, mode);
  NumericVector LAI_expanded = cohortLAI(x, SpParams, gdd, mode);
  IntegerVector SP(ntree+nshrub);
  NumericVector H(ntree+nshrub);
  
  IntegerVector treeSP = Rcpp::as<Rcpp::IntegerVector>(treeData["Species"]);
  NumericVector treeH = treeData["Height"];
  NumericVector treeDBH = treeData["DBH"];
  IntegerVector shrubSP = Rcpp::as<Rcpp::IntegerVector>(shrubData["Species"]);
  NumericVector shrubH = shrubData["Height"];  
  NumericVector shrubCover = shrubData["Cover"];  
  
  NumericVector N = cohortDensity(x, SpParams, mode);
    
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


//' @rdname modelInput
// [[Rcpp::export("forest2belowground")]]
NumericMatrix forest2belowground(List x, List soil) {
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  DataFrame shrubData = Rcpp::as<Rcpp::DataFrame>(x["shrubData"]);
  NumericVector d = soil["dVec"];
  int nlayers = d.size();
  NumericVector treeZ50 = treeData["Z50"];
  NumericVector treeZ95 = treeData["Z95"];
  NumericVector shrubZ50 = shrubData["Z50"];  
  NumericVector shrubZ95 = shrubData["Z95"];  
  NumericMatrix V = ldrDistribution(treeZ50, shrubZ50, treeZ95, shrubZ95, d);
  V.attr("dimnames") = List::create(cohortIDs(x), layerNames(nlayers));
  return(V);
}