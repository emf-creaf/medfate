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
  sba.attr("names") = speciesCharacterParameterFromIndex(uniqueSp, SpParams, "Name");
  return(sba);
}


/**
 *  Basal area (BA, in m2/ha)
 */
// [[Rcpp::export(".treeBasalArea")]]
NumericVector treeBasalArea(NumericVector N, NumericVector dbh) {
  int ncoh = N.size(); //N is density of individuals (ind/ha) in the cell
  NumericVector BA(ncoh, NA_REAL); 
  for(int i=0;i<ncoh;i++) {
    if(!NumericVector::is_na(dbh[i])) BA[i] = N[i]*3.141593*pow(dbh[i]/200,2.0); //Basal area in m2/ha
  }
  return(BA);
}
NumericVector largerTreeBasalArea(NumericVector N, NumericVector dbh, double self_include_prop = 0.5) {
  int ncoh = N.size();
  NumericVector BA = treeBasalArea(N, dbh); 
  NumericVector ltBA(ncoh, NA_REAL);
  for(int i=0;i<ncoh;i++) {
    if(!NumericVector::is_na(BA[i])) {
      ltBA[i] = 0.0;
      for(int j=0;j<ncoh;j++) {
        if(i==j) ltBA[i] += (BA[j]*self_include_prop); //add half of its own basal area
        else if(dbh[j]>dbh[i]) ltBA[i] += BA[j];
      }
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
NumericVector shrubIndividualAreaAllometric(IntegerVector SP, NumericVector Cover, NumericVector H, DataFrame SpParams){
  NumericVector aShrubArea = speciesNumericParameterWithImputation(SP,SpParams, "a_ash",true, true);
  NumericVector bShrubArea = speciesNumericParameterWithImputation(SP,SpParams, "b_ash",true, true);
  int ncoh = SP.size();
  NumericVector areaind(ncoh);
  for(int i=0;i<ncoh;i++) {
    if((!NumericVector::is_na(Cover[i])) && (!NumericVector::is_na(H[i]))) {
      areaind[i] = aShrubArea[i]*pow(H[i],bShrubArea[i])/10000.0; 
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
NumericVector shrubCrownRatioAllometric(IntegerVector SP, DataFrame SpParams) {
  return(speciesNumericParameterWithImputation(SP, SpParams, "cr", true, true));
}

double crownCompetitionFactorAllometric(NumericVector N, NumericVector dbh, NumericVector Acw, NumericVector Bcw) {
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
double crownCompetitionFactorAllometric(IntegerVector SP, NumericVector N, NumericVector dbh, DataFrame SpParams) {
  NumericVector Acw = speciesNumericParameterWithImputation(SP, SpParams, "a_cw",true, true);
  NumericVector Bcw = speciesNumericParameterWithImputation(SP, SpParams, "b_cw",true, true);
  return(crownCompetitionFactorAllometric(N,dbh,Acw,Bcw));
}


NumericVector treeCrownRatioAllometric(NumericVector N, NumericVector dbh, NumericVector H, 
                                NumericVector Acw, NumericVector Bcw,
                                NumericVector Acr, NumericVector B1cr, NumericVector B2cr, NumericVector B3cr,
                                NumericVector C1cr, NumericVector C2cr) {
  NumericVector BAL = largerTreeBasalArea(N, dbh);
  double ccf = crownCompetitionFactorAllometric(N, dbh, Acw, Bcw);
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

NumericVector treeCrownRatioAllometric(IntegerVector SP, NumericVector N, NumericVector dbh, NumericVector H, DataFrame SpParams) {
  NumericVector BAL = largerTreeBasalArea(N, dbh);
  double ccf = crownCompetitionFactorAllometric(SP, N, dbh, SpParams);
  // Rcout<<ccf<<"\n";
  NumericVector Acr = speciesNumericParameterWithImputation(SP, SpParams, "a_cr",true, true);
  NumericVector B1cr = speciesNumericParameterWithImputation(SP, SpParams, "b_1cr",true, true);
  NumericVector B2cr = speciesNumericParameterWithImputation(SP, SpParams, "b_2cr",true, true);
  NumericVector B3cr = speciesNumericParameterWithImputation(SP, SpParams, "b_3cr",true, true);
  NumericVector C1cr = speciesNumericParameterWithImputation(SP, SpParams, "c_1cr",true, true);
  NumericVector C2cr = speciesNumericParameterWithImputation(SP, SpParams, "c_2cr",true, true);
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
NumericVector treeFoliarBiomassAllometric(IntegerVector SP, NumericVector N, NumericVector dbh, DataFrame SpParams, double gdd = NA_REAL,
                                          bool competitionEffect = true){
  NumericVector afbt = speciesNumericParameterWithImputation(SP, SpParams, "a_fbt",true, true);
  NumericVector bfbt = speciesNumericParameterWithImputation(SP, SpParams, "b_fbt",true, true);
  NumericVector cfbt = speciesNumericParameterWithImputation(SP, SpParams, "c_fbt",true, true);
  int ncoh = N.size();
  NumericVector lb(ncoh);
  for(int i=0;i<ncoh;i++) {
    lb[i] = ((N[i]/10000.0)*afbt[i]*pow(std::min(100.0,dbh[i]), bfbt[i]));
    lb[i] = lb[i] * exp(-0.0001*N[i]);//Correct for high density packing
  }
  if(competitionEffect) {
    NumericVector ltba = largerTreeBasalArea(N,dbh, 1.0); //Allometries were calibrated including the target cohort
    for(int i=0;i<ncoh;i++) {
      lb[i] = lb[i]*exp(cfbt[i]*std::min(ltba[i], 100.0));
    }
  }
  if(!NumericVector::is_na(gdd)) {
    NumericVector Sgdd = speciesNumericParameterWithImputation(SP, SpParams, "Sgdd", true, true);
    for(int i=0;i<ncoh;i++) {
      if(!NumericVector::is_na(SP[i])) lb[i] = lb[i]*leafDevelopmentStatus(Sgdd[i], gdd);
    }
  }
  
  return(lb);
}

NumericVector shrubFoliarBiomassAllometric(IntegerVector SP, NumericVector Cover, NumericVector H, 
                                 DataFrame SpParams, double gdd = NA_REAL, double treeLAI = 0.0,
                                 bool competitionEffect = true){
  NumericVector aShrubFuel = speciesNumericParameterWithImputation(SP, SpParams, "a_bsh",true, true);
  NumericVector bShrubFuel = speciesNumericParameterWithImputation(SP, SpParams, "b_bsh",true, true);
  NumericVector Sgdd = speciesNumericParameterWithImputation(SP, SpParams, "Sgdd", true, true);
  NumericVector pDead = speciesNumericParameterWithImputation(SP, SpParams, "pDead", true, true);
  NumericVector fTreeFuel = speciesNumericParameterWithImputation(SP, SpParams, "r635", true, true);
  int ncoh = SP.size();
  double W = 0.0; //Fine fuel
  NumericVector fb(ncoh);
  NumericVector areaind = shrubIndividualAreaAllometric(SP,Cover,H,SpParams); 
  double volind = NA_REAL,weightkgind = NA_REAL;
  for(int i=0;i<ncoh;i++) {
    // Rcout<<i<<": "<< H[i]<<" "<<CR[i]<<" "<<aShrubFuel[i]<<" "<<bShrubFuel[i]<< " "<< pDead[i]<<" "<<fTreeFuel[i]<<" "<< areaind[i]<<".\n";
    if((!NumericVector::is_na(Cover[i])) && (!NumericVector::is_na(H[i]))) {
      volind = areaind[i]*(H[i]/100.0); // Phytovolume of an individual (in m3)
      weightkgind = aShrubFuel[i]*pow(volind,bShrubFuel[i]); //Fuel (in kg) of an individual (includes dead fuels)
      weightkgind = weightkgind - (weightkgind*pDead[i]); //Removes dead fuels
      // Rcout<<volind<< " "<<weightkgind;
      if(areaind[i]>0.0) {
        // multiply by 'number of individuals' per m2 
        W = weightkgind*(Cover[i]/(100*areaind[i]));  //Fine fuel (kg/m2)
        fb[i] = W/fTreeFuel[i]; //Foliar biomass (kg/m2)
        if(competitionEffect) fb[i] = fb[i]*exp(-0.235*treeLAI); //Correct depending on tree leaf area
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

NumericVector treeCoverAllometric(IntegerVector SP, NumericVector N, NumericVector dbh, DataFrame SpParams){
  NumericVector Acw = speciesNumericParameterWithImputation(SP, SpParams, "a_cw",true, true);
  NumericVector Bcw = speciesNumericParameterWithImputation(SP, SpParams, "b_cw",true, true);
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

/**
 *  Shrub phytovolume (in m3/m2)
 */
// [[Rcpp::export(".shrubPhytovolume")]]
NumericVector shrubPhytovolumeAllometric(IntegerVector SP, NumericVector Cover, NumericVector H, DataFrame SpParams){
  int ncoh = Cover.size();
  NumericVector vol(ncoh);
  NumericVector areaind = shrubIndividualAreaAllometric(SP,Cover,H,SpParams); //area of an individual (in m2)
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
NumericVector treeFuelAllometric(IntegerVector SP, NumericVector FB,
                                 DataFrame SpParams, double gdd = NA_REAL, bool includeDead = true){
  NumericVector Sgdd = speciesNumericParameterWithImputation(SP, SpParams, "Sgdd", true, true);
  NumericVector fTreeFuel = speciesNumericParameterWithImputation(SP, SpParams, "r635", true, true);
  NumericVector pDead = speciesNumericParameterWithImputation(SP, SpParams, "pDead", true, true);
  int ncoh = SP.size();
  double ftf = 0.0, btf = 0.0;
  NumericVector fuel(ncoh, NA_REAL);
  for(int i=0;i<ncoh;i++) {
    if(!NumericVector::is_na(FB[i])){
      ftf = FB[i]; //Foliar biomass (kg per m2)
      btf = ftf*(fTreeFuel[i]-1.0); // Small branch fuels (proportion of foliar fuels)
      if(!NumericVector::is_na(gdd)) { //Apply phenology correction to foliar fuels
        ftf = ftf*leafDevelopmentStatus(Sgdd[i], gdd); 
      } 
      fuel[i] =  ftf + btf; //Tree fuel (kg per m2) is sum of both fuels
      if(includeDead) {
        fuel[i] = fuel[i] + fuel[i]*pDead[i]; //If required add fine dead fuels (proportion of live fuels)
      }
    }
  }
  return(fuel);
}


NumericVector shrubFuelAllometric(IntegerVector SP, NumericVector FB, DataFrame SpParams, 
                                  double gdd = NA_REAL, bool includeDead = true){

  int ncoh = SP.size();
  NumericVector fShrubFuelRatio = speciesNumericParameterWithImputation(SP, SpParams, "r635", true, true);
  NumericVector Sgdd = speciesNumericParameterWithImputation(SP, SpParams, "Sgdd", true, true);
  NumericVector W(ncoh,NA_REAL);
    
  for(int i=0;i<ncoh;i++) {
    W[i] = FB[i]*fShrubFuelRatio[i]; //fine fuel biomass
    //Remove (if necessary), the weight due to leaves that are not there
    if(!NumericVector::is_na(gdd)) {
      double bsf = W[i] - FB[i]; //branch biomass
      W[i] = bsf + FB[i]*leafDevelopmentStatus(Sgdd[i], gdd); 
    }
  }
  return(W);
}


/**
 *  Leaf Area Index (LAI)
 */
NumericVector treeLAIAllometric(IntegerVector SP, NumericVector N, NumericVector dbh, DataFrame SpParams, double gdd = NA_REAL,
                                bool competitionEffect = true){
  NumericVector SLA = speciesNumericParameterWithImputation(SP, SpParams, "SLA", true, true); // m2/kg (=mg/mm2)
  NumericVector lb = treeFoliarBiomassAllometric(SP, N, dbh, SpParams, gdd, competitionEffect); //kg per m2
  int ncoh = N.size();
  NumericVector lai(ncoh);
  for(int i=0;i<ncoh;i++) {
     lai[i] = (SLA[i]*lb[i]); 
  }
  return(lai);
}

NumericVector shrubLAIAllometric(IntegerVector SP, NumericVector Cover, NumericVector H, DataFrame SpParams, 
                       double gdd = NA_REAL, double treeLAI = 0.0, bool competitionEffect = true){
  NumericVector SLA = speciesNumericParameterWithImputation(SP, SpParams, "SLA", true, true); // m2/kg (=mg/mm2)
  NumericVector lb = shrubFoliarBiomassAllometric(SP, Cover, H, SpParams, gdd, treeLAI, competitionEffect); //kg per m2
  int ncoh = SP.size();
  NumericVector lai(ncoh);
  for(int i=0;i<ncoh;i++) {
    lai[i] = (SLA[i]*lb[i]); 
  }
  return(lai);
}







//' Woody plant cohort description functions
//'
//' Functions to calculate attributes of woody plants in a \code{\link{forest}} object.
//' 
//' @param x An object of class \code{\link{forest}}.
//' @param SpParams A data frame with species parameters (see \code{\link{SpParamsMED}}).
//' @param parName A string with a parameter name.
//' @param gdd Growth degree days (to account for leaf phenology effects).
//' @param AET Actual annual evapotranspiration (in mm).
//' @param smallBranchDecompositionRate Decomposition rate of small branches.
//' @param includeDead A flag to indicate that standing dead fuels (dead branches) are included.
//' @param treeOffset,shrubOffset Integers to offset cohort IDs.
//' @param fillMissing A boolean flag to try imputation on missing values.
//' @param fillWithGenus A boolean flag to try imputation of missing values using genus values.
//' @param self_proportion Proportion of the target cohort included in the assessment
//' @param bounded A boolean flag to indicate that extreme values should be prevented (maximum tree LAI = 7 and maximum shrub LAI = 3)
//' 
//' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
//' 
//' @return
//' A vector with values for each woody plant cohort of the input \code{\link{forest}} object:
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
//'   \item{\code{plant_fuelLoading}: Fine fuel load (in kg/m2).}
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
//' data(exampleforest)
//' 
//' #A plant-level way to obtain stand basal area
//' sum(plant_basalArea(exampleforest, SpParamsMED), na.rm=TRUE)
//' 
//' #The analogous plant-level function for LAI
//' sum(plant_LAI(exampleforest, SpParamsMED))
//'   
//' #The analogous plant-level function for fuel loading
//' sum(plant_fuelLoading(exampleforest, SpParamsMED))
//'       
//' #Summary function for 'forest' objects can be also used
//' summary(exampleforest, SpParamsMED)
//' 
//' #Cohort IDs in the models
//' plant_ID(exampleforest, SpParamsMED)
//'       
//' @name plant_values
//' @keywords internal
// [[Rcpp::export("plant_ID")]]
CharacterVector cohortIDs(List x, DataFrame SpParams, int treeOffset = 0, int shrubOffset = 0) {
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  DataFrame shrubData = Rcpp::as<Rcpp::DataFrame>(x["shrubData"]);
  int ntree = treeData.nrows();
  int nshrub = shrubData.nrows();
  int numCohorts  = ntree+nshrub;
  IntegerVector treeSP, shrubSP;
  if((TYPEOF(treeData["Species"]) == INTSXP) || (TYPEOF(treeData["Species"]) == REALSXP)) {
    treeSP = Rcpp::as<Rcpp::IntegerVector>(treeData["Species"]);
  } else {
    CharacterVector tspecies = Rcpp::as<Rcpp::CharacterVector>(treeData["Species"]);
    treeSP = speciesIndex(tspecies, SpParams);
  }
  if((TYPEOF(shrubData["Species"]) == INTSXP) || (TYPEOF(shrubData["Species"]) == REALSXP)) {
    shrubSP = Rcpp::as<Rcpp::IntegerVector>(shrubData["Species"]);  
  } else {
    CharacterVector sspecies = Rcpp::as<Rcpp::CharacterVector>(shrubData["Species"]);
    shrubSP = speciesIndex(sspecies, SpParams);
  }
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
//' @keywords internal
// [[Rcpp::export("plant_basalArea")]]
NumericVector cohortBasalArea(List x, DataFrame SpParams) {
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  DataFrame shrubData = Rcpp::as<Rcpp::DataFrame>(x["shrubData"]);
  NumericVector tba = treeBasalArea(treeData["N"], treeData["DBH"]);
  NumericVector ba(tba.size()+shrubData.nrows(), NA_REAL);
  for(int i=0;i<tba.size();i++) {
    ba[i] = tba[i];
  }
  ba.attr("names") = cohortIDs(x, SpParams);
  return(ba);
}

//' @rdname plant_values
//' @keywords internal
// [[Rcpp::export("plant_largerTreeBasalArea")]]
NumericVector cohortLargerTreeBasalArea(List x, DataFrame SpParams, double self_proportion = 0.5) {
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  DataFrame shrubData = Rcpp::as<Rcpp::DataFrame>(x["shrubData"]);
  NumericVector tba = largerTreeBasalArea(treeData["N"], treeData["DBH"], self_proportion);
  NumericVector ba(tba.size()+shrubData.nrows(), NA_REAL);
  for(int i=0;i<tba.size();i++) {
    ba[i] = tba[i];
  }
  ba.attr("names") = cohortIDs(x, SpParams);
  return(ba);
}

//' @rdname plant_values
//' @keywords internal
// [[Rcpp::export("plant_cover")]]
NumericVector cohortCover(List x, DataFrame SpParams) {
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  DataFrame shrubData = Rcpp::as<Rcpp::DataFrame>(x["shrubData"]);
  NumericVector cov(treeData.nrows()+shrubData.nrows(), NA_REAL);
  NumericVector tcover;
  IntegerVector treeSP;
  if((TYPEOF(treeData["Species"]) == INTSXP) || (TYPEOF(treeData["Species"]) == REALSXP)) {
    treeSP = Rcpp::as<Rcpp::IntegerVector>(treeData["Species"]);
  } else {
    CharacterVector tspecies = Rcpp::as<Rcpp::CharacterVector>(treeData["Species"]);
    treeSP = speciesIndex(tspecies, SpParams);
  }
  tcover = treeCoverAllometric(treeSP, treeData["N"], treeData["DBH"],
                               SpParams);
  for(int i=0;i<tcover.size();i++) {
    cov[i] = tcover[i];
  }
  NumericVector shcover = shrubData["Cover"];
  for(int i=0;i<shcover.size();i++) {
    cov[i+treeData.nrows()] = shcover[i];
  }
  cov.attr("names") = cohortIDs(x, SpParams);
  return(cov);
}

//' @rdname plant_values
//' @keywords internal
// [[Rcpp::export("plant_species")]]
IntegerVector cohortSpecies(List x, DataFrame SpParams) {
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  DataFrame shrubData = Rcpp::as<Rcpp::DataFrame>(x["shrubData"]);
  int ntree = treeData.nrows();
  int nshrub = shrubData.nrows();
  IntegerVector treeSP, shrubSP;
  if((TYPEOF(treeData["Species"]) == INTSXP) || (TYPEOF(treeData["Species"]) == REALSXP)) {
    treeSP = Rcpp::as<Rcpp::IntegerVector>(treeData["Species"]);
  } else {
    CharacterVector tspecies = Rcpp::as<Rcpp::CharacterVector>(treeData["Species"]);
    treeSP = speciesIndex(tspecies, SpParams);
  }
  if((TYPEOF(shrubData["Species"]) == INTSXP) || (TYPEOF(shrubData["Species"]) == REALSXP)) {
    shrubSP = Rcpp::as<Rcpp::IntegerVector>(shrubData["Species"]);  
  } else {
    CharacterVector sspecies = Rcpp::as<Rcpp::CharacterVector>(shrubData["Species"]);
    shrubSP = speciesIndex(sspecies, SpParams);
  }
  int numCohorts  = ntree+nshrub;
  IntegerVector SP(numCohorts);
  for(int i=0;i<ntree;i++) {
    SP[i] = treeSP[i];
  }
  for(int i=0;i<nshrub;i++) {
    SP[ntree+i] = shrubSP[i];
  }
  SP.attr("names") = cohortIDs(x, SpParams);
  return(SP);
}

//' @rdname plant_values
//' @keywords internal
// [[Rcpp::export("plant_speciesName")]]
CharacterVector cohortSpeciesName(List x, DataFrame SpParams) {
  CharacterVector sn = cohortCharacterParameter(x, SpParams, "Name");
  sn.attr("names") = cohortIDs(x, SpParams);
  return(sn);
}

//' @rdname plant_values
//' @keywords internal
// [[Rcpp::export("plant_density")]]
NumericVector cohortDensity(List x, DataFrame SpParams) {
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  DataFrame shrubData = Rcpp::as<Rcpp::DataFrame>(x["shrubData"]);
  int ntree = treeData.nrows();
  int nshrub = shrubData.nrows();
  NumericVector treeN = treeData["N"];
  IntegerVector shrubSP;
  if((TYPEOF(shrubData["Species"]) == INTSXP) || (TYPEOF(shrubData["Species"]) == REALSXP)) {
    shrubSP = Rcpp::as<Rcpp::IntegerVector>(shrubData["Species"]);  
  } else {
    CharacterVector sspecies = Rcpp::as<Rcpp::CharacterVector>(shrubData["Species"]);
    shrubSP = speciesIndex(sspecies, SpParams);
  }

  NumericVector shrubHeight = shrubData["Height"];
  int numCohorts  = ntree+nshrub;
  NumericVector N(numCohorts);
  for(int i=0;i<ntree;i++) {
    N[i] = treeN[i];
  }
  NumericVector shrubCover = shrubData["Cover"];
  NumericVector shrubArea = shrubIndividualAreaAllometric(shrubSP,shrubCover, shrubHeight, SpParams); 
  for(int i=0;i<nshrub;i++) {
    N[ntree+i] = 10000.0*(shrubCover[i]/(100.0*shrubArea[i]));
  }
  N.attr("names") = cohortIDs(x, SpParams);
  return(N);
}

//' @rdname plant_values
//' @keywords internal
// [[Rcpp::export("plant_height")]]
NumericVector cohortHeight(List x, DataFrame SpParams) {
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
  H.attr("names") = cohortIDs(x, SpParams);
  return(H);
}

//' @rdname plant_values
//' @keywords internal
// [[Rcpp::export("plant_individualArea")]]
NumericVector cohortIndividualArea(List x, DataFrame SpParams){
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  DataFrame shrubData = Rcpp::as<Rcpp::DataFrame>(x["shrubData"]);
  int ntree = treeData.nrows();
  int nshrub = shrubData.nrows();
  NumericVector treeN = treeData["N"];
  IntegerVector shrubSP;
  if((TYPEOF(shrubData["Species"]) == INTSXP) || (TYPEOF(shrubData["Species"]) == REALSXP)) {
    shrubSP = Rcpp::as<Rcpp::IntegerVector>(shrubData["Species"]);  
  } else {
    CharacterVector sspecies = Rcpp::as<Rcpp::CharacterVector>(shrubData["Species"]);
    shrubSP = speciesIndex(sspecies, SpParams);
  }
  NumericVector shrubHeight = shrubData["Height"];
  int numCohorts  = ntree+nshrub;
  NumericVector indArea(numCohorts, NA_REAL);
  NumericVector shrubCover = shrubData["Cover"];
  NumericVector shrubArea = shrubIndividualAreaAllometric(shrubSP,shrubCover, shrubHeight, SpParams); 
  for(int i=0;i<nshrub;i++) {
    indArea[ntree+i] = shrubArea[i];
  }
  indArea.attr("names") = cohortIDs(x, SpParams);
  return(indArea);
}

//' @rdname plant_values
//' @keywords internal
// [[Rcpp::export("plant_crownRatio")]]
NumericVector cohortCrownRatio(List x, DataFrame SpParams) {
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  DataFrame shrubData = Rcpp::as<Rcpp::DataFrame>(x["shrubData"]);
  int ntree = treeData.nrows();
  int nshrub = shrubData.nrows();
  NumericVector treeCR(ntree, NA_REAL), shrubCR(nshrub, NA_REAL);
  
  IntegerVector treeSP, shrubSP;
  if((TYPEOF(treeData["Species"]) == INTSXP) || (TYPEOF(treeData["Species"]) == REALSXP)) {
    treeSP = Rcpp::as<Rcpp::IntegerVector>(treeData["Species"]);
  } else {
    CharacterVector tspecies = Rcpp::as<Rcpp::CharacterVector>(treeData["Species"]);
    treeSP = speciesIndex(tspecies, SpParams);
  }
  if((TYPEOF(shrubData["Species"]) == INTSXP) || (TYPEOF(shrubData["Species"]) == REALSXP)) {
    shrubSP = Rcpp::as<Rcpp::IntegerVector>(shrubData["Species"]);  
  } else {
    CharacterVector sspecies = Rcpp::as<Rcpp::CharacterVector>(shrubData["Species"]);
    shrubSP = speciesIndex(sspecies, SpParams);
  }
  
  if(treeData.containsElementNamed("CrownRatio")) {
    treeCR = treeData["CrownRatio"];
  } 
  NumericVector treeCRAllom = treeCrownRatioAllometric(treeSP,treeData["N"], treeData["DBH"], treeData["Height"], SpParams); 
  for(int i=0;i<ntree;i++) {
    if(NumericVector::is_na(treeCR[i])) treeCR[i] = treeCRAllom[i];
  }
  if(shrubData.containsElementNamed("CrownRatio")) {
    shrubCR = shrubData["CrownRatio"];
  } 
  NumericVector shrubCRAllom = shrubCrownRatioAllometric(shrubSP, SpParams);
  for(int i=0;i<nshrub;i++) {
    if(NumericVector::is_na(shrubCR[i])) shrubCR[i] = shrubCRAllom[i];
  }
  
  int numCohorts  = ntree+nshrub;
  NumericVector CR(numCohorts);
  for(int i=0;i<ntree;i++) {
    CR[i] = treeCR[i];
  }
  for(int i=0;i<nshrub;i++) {
    CR[ntree+i] = shrubCR[i];
  }
  CR.attr("names") = cohortIDs(x, SpParams);
  return(CR);
}

//' @rdname plant_values
//' @keywords internal
// [[Rcpp::export("plant_crownBaseHeight")]]
NumericVector cohortCrownBaseHeight(List x, DataFrame SpParams) {
  NumericVector CR = cohortCrownRatio(x, SpParams);
  NumericVector H = cohortHeight(x, SpParams);
  int numCohorts = H.size();
  NumericVector CBH(numCohorts);
  for(int i=0;i<numCohorts;i++) {
    CBH[i] = H[i]*(1.0-CR[i]);
  }
  CBH.attr("names") = cohortIDs(x, SpParams);
  return(CBH);
}

//' @rdname plant_values
//' @keywords internal
// [[Rcpp::export("plant_crownLength")]]
NumericVector cohortCrownLength(List x, DataFrame SpParams) {
  NumericVector CR = cohortCrownRatio(x, SpParams);
  NumericVector H = cohortHeight(x, SpParams);
  int numCohorts = H.size();
  NumericVector CL(numCohorts);
  for(int i=0;i<numCohorts;i++) {
    CL[i] = H[i]*CR[i];
  }
  CL.attr("names") = cohortIDs(x, SpParams);
  return(CL);
}

//' @rdname plant_values
//' @keywords internal
// [[Rcpp::export("plant_foliarBiomass")]]
NumericVector cohortFoliarBiomass(List x, DataFrame SpParams, double gdd = NA_REAL, bool competitionEffect = true) {
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  DataFrame shrubData = Rcpp::as<Rcpp::DataFrame>(x["shrubData"]);
  NumericVector tFB(treeData.nrows(), NA_REAL), shFB(shrubData.nrows(), NA_REAL);
  NumericVector tLAI(treeData.nrows(), NA_REAL);
  
  IntegerVector treeSP, shrubSP;
  if((TYPEOF(treeData["Species"]) == INTSXP) || (TYPEOF(treeData["Species"]) == REALSXP)) {
    treeSP = Rcpp::as<Rcpp::IntegerVector>(treeData["Species"]);
  } else {
    CharacterVector tspecies = Rcpp::as<Rcpp::CharacterVector>(treeData["Species"]);
    treeSP = speciesIndex(tspecies, SpParams);
  }
  if((TYPEOF(shrubData["Species"]) == INTSXP) || (TYPEOF(shrubData["Species"]) == REALSXP)) {
    shrubSP = Rcpp::as<Rcpp::IntegerVector>(shrubData["Species"]);  
  } else {
    CharacterVector sspecies = Rcpp::as<Rcpp::CharacterVector>(shrubData["Species"]);
    shrubSP = speciesIndex(sspecies, SpParams);
  }
  
  NumericVector tSgdd = speciesNumericParameterWithImputation(treeSP, SpParams, "Sgdd", true, true);
  NumericVector tSLA = speciesNumericParameterWithImputation(treeSP, SpParams, "SLA", true, true); // m2/kg (=mg/mm2)
  if(treeData.containsElementNamed("FoliarBiomass")) {
    // If FoliarBiomass is already in treeData, calculate LAI and only apply phenology
    tFB = treeData["FoliarBiomass"];
    for(int i=0;i<tFB.size();i++) {
      tLAI[i] = tFB[i]*tSLA[i];
    }
    if(!NumericVector::is_na(gdd)) {
      for(int i=0;i<tFB.size();i++) {
        tFB[i] = tFB[i]*leafDevelopmentStatus(tSgdd[i], gdd);
      }
    }
  } else if(treeData.containsElementNamed("LAI")) {
    //If LAI is in treeData, apply SLA and phenology
    tLAI = treeData["LAI"];
    for(int i=0;i<tLAI.size();i++) {
      tFB[i] = tLAI[i]/tSLA[i];
    }
    if(!NumericVector::is_na(gdd)) {
      for(int i=0;i<tLAI.size();i++) {
        tFB[i] = tFB[i]*leafDevelopmentStatus(tSgdd[i], gdd);
      }
    }
  } 
  // Apply foliar biomass allometries to fill gaps
  NumericVector tFBAllom = treeFoliarBiomassAllometric(treeSP, treeData["N"], treeData["DBH"], SpParams, gdd, competitionEffect);
  for(int i=0;i<treeData.nrows();i++) {
    if(NumericVector::is_na(tFB[i])) {
      tFB[i] = tFBAllom[i];
      tLAI[i] = tFB[i]*tSLA[i];
      if(!NumericVector::is_na(gdd)) {
        tFB[i] = tFB[i]*leafDevelopmentStatus(tSgdd[i], gdd);
      }
    }
  }
  
  NumericVector shSgdd = speciesNumericParameterWithImputation(shrubSP, SpParams, "Sgdd", true, true);
  NumericVector shSLA = speciesNumericParameterWithImputation(shrubSP, SpParams, "SLA", true, true); // m2/kg (=mg/mm2)
  if(shrubData.containsElementNamed("FoliarBiomass")) {
    // If FoliarBiomass is already in shrubData, only apply phenology
    shFB = shrubData["FoliarBiomass"];
    if(!NumericVector::is_na(gdd)) {
      for(int i=0;i<shFB.size();i++) {
        shFB[i] = shFB[i]*leafDevelopmentStatus(shSgdd[i], gdd);
      }
    }
  } else if(shrubData.containsElementNamed("LAI")) {
    //If LAI is in shrubData, apply SLA and phenology
    NumericVector shLAI = shrubData["LAI"];
    for(int i=0;i<shLAI.size();i++) {
      shFB[i] = shLAI[i]/shSLA[i];
    }
    if(!NumericVector::is_na(gdd)) {
      for(int i=0;i<shFB.size();i++) {
        shFB[i] = shFB[i]*leafDevelopmentStatus(shSgdd[i], gdd);
      }
    }
  } 
  
  // Apply foliar biomass allometries to fill gaps
  NumericVector tba = treeBasalArea(treeData["N"], treeData["DBH"]);
  //Sum tree LAI for shrub correction
  double treeLAI = sum(tLAI);
  NumericVector shFBAllom= shrubFoliarBiomassAllometric(shrubSP, shrubData["Cover"], shrubData["Height"], 
                                                        SpParams, gdd, treeLAI, competitionEffect);
  for(int i=0;i<shrubData.nrows();i++) {
    if(NumericVector::is_na(shFB[i])) {
      shFB[i] = shFBAllom[i];
      if(!NumericVector::is_na(gdd)) shFB[i] = shFB[i]*leafDevelopmentStatus(shSgdd[i], gdd);
    }
  }
  
  //Copy values to single vector
  NumericVector FB(tFB.size()+shFB.size());
  for(int i=0;i<tFB.size();i++) {
    FB[i] = tFB[i];
  }
  for(int i=0;i<shFB.size();i++) {
    FB[i+tFB.size()] = shFB[i];
  }
  FB.attr("names") = cohortIDs(x, SpParams);
  return(FB);
}

//' @rdname plant_values
//' @keywords internal
// [[Rcpp::export("plant_fuelLoading")]]
NumericVector cohortFuelLoading(List x, DataFrame SpParams, double gdd = NA_REAL, bool includeDead = true){
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  DataFrame shrubData = Rcpp::as<Rcpp::DataFrame>(x["shrubData"]);
  NumericVector tFuel(treeData.nrows(), NA_REAL), shFuel(shrubData.nrows(), NA_REAL);
  NumericVector tFB(treeData.nrows(), NA_REAL), shFB(shrubData.nrows(), NA_REAL);
  NumericVector tLAI(treeData.nrows(), NA_REAL);
  
  IntegerVector treeSP, shrubSP;
  if((TYPEOF(treeData["Species"]) == INTSXP) || (TYPEOF(treeData["Species"]) == REALSXP)) {
    treeSP = Rcpp::as<Rcpp::IntegerVector>(treeData["Species"]);
  } else {
    CharacterVector tspecies = Rcpp::as<Rcpp::CharacterVector>(treeData["Species"]);
    treeSP = speciesIndex(tspecies, SpParams);
  }
  if((TYPEOF(shrubData["Species"]) == INTSXP) || (TYPEOF(shrubData["Species"]) == REALSXP)) {
    shrubSP = Rcpp::as<Rcpp::IntegerVector>(shrubData["Species"]);  
  } else {
    CharacterVector sspecies = Rcpp::as<Rcpp::CharacterVector>(shrubData["Species"]);
    shrubSP = speciesIndex(sspecies, SpParams);
  }
  
  NumericVector tSLA = speciesNumericParameterWithImputation(treeSP, SpParams, "SLA", true, true); // m2/kg (=mg/mm2)
  if(treeData.containsElementNamed("FuelLoading")) {
    tFuel = treeData["FuelLoading"];
  } else if(treeData.containsElementNamed("FoliarBiomass")) {
    tFB = treeData["FoliarBiomass"];
    for(int i=0;i<tFB.size();i++) {
      tLAI[i] = tFB[i]*tSLA[i];
    }
  } else if(treeData.containsElementNamed("LAI")) {
    //If LAI is in treeData, apply SLA
    tLAI = treeData["LAI"];
    for(int i=0;i<tLAI.size();i++) {
      tFB[i] = tLAI[i]/tSLA[i];
    }
  }
  //Use foliar allometries to fill foliar biomass gaps
  NumericVector tFBAllom = treeFoliarBiomassAllometric(treeSP,treeData["N"], treeData["DBH"], SpParams, NA_REAL, true); //Do not include phenology (to have correct estimates of branch biomass)
  for(int i=0;i<tFB.size();i++) {
    if(NumericVector::is_na(tFB[i])) {
      tFB[i] = tFBAllom[i];
      tLAI[i] = tFB[i]*tSLA[i];
    }
  }
  
  //Estimate fuel derived from foliar allometries or measured foliar biomass
  NumericVector tFuelAllom = treeFuelAllometric(treeSP, tFB, SpParams, gdd, includeDead);
  
  //Replace missing fuel values with allometric estimates
  for(int i=0;i<tFuel.size();i++) {
    if(NumericVector::is_na(tFuel[i])) tFuel[i] = tFuelAllom[i];
  } 

    
  if(shrubData.containsElementNamed("FuelLoading")) {
    shFuel = shrubData["FuelLoading"];
  } else if(shrubData.containsElementNamed("FoliarBiomass")) {
    shFB = shrubData["FoliarBiomass"];
  } else if(shrubData.containsElementNamed("LAI")){
    //If LAI is in shrubData, apply SLA and phenology
    NumericVector shLAI = shrubData["LAI"];
    NumericVector shSLA = speciesNumericParameterWithImputation(shrubSP, SpParams, "SLA", true, true); // m2/kg (=mg/mm2)
    for(int i=0;i<shLAI.size();i++) {
      shFB[i] = shLAI[i]/shSLA[i];
    }
  }
  
  //Use foliar allometries to fill foliar biomass gaps
  double treeLAI = sum(tLAI);
  NumericVector shFBAllom= shrubFoliarBiomassAllometric(shrubSP, shrubData["Cover"], shrubData["Height"], 
                                                        SpParams, gdd, treeLAI, true);
  for(int i=0;i<shrubData.nrows();i++) {
    if(NumericVector::is_na(shFB[i])) shFB[i] = shFBAllom[i];
  }
  
  //Estimate fuel derived from foliar allometries or measured foliar biomass
  NumericVector shFuelAllom = shrubFuelAllometric(shrubSP, shFB, 
                                                  SpParams, gdd, includeDead);
  
  //Replace missing fuel values with allometric estimates
  for(int i=0;i<shFuel.size();i++) {
    if(NumericVector::is_na(shFuel[i])) shFuel[i] = shFuelAllom[i];
  }
  
  
  NumericVector fuel(tFuel.size()+shFuel.size());
  for(int i=0;i<tFuel.size();i++) {
    fuel[i] = tFuel[i];
  }
  for(int i=0;i<shFuel.size();i++) {
    fuel[i+tFuel.size()] = shFuel[i];
  }
  fuel.attr("names") = cohortIDs(x, SpParams);
  return(fuel);
}

//' @rdname plant_values
//' @keywords internal
// [[Rcpp::export("plant_equilibriumLeafLitter")]]
NumericVector cohortEquilibriumLeafLitter(List x, DataFrame SpParams, double AET = 800) {
  NumericVector fb = cohortFoliarBiomass(x, SpParams, NA_REAL);
  NumericVector ld = cohortNumericParameterWithImputation(x, SpParams, "LeafDuration");
  NumericVector lignin = cohortNumericParameterWithImputation(x, SpParams, "LigninPercent", true);
  int ncoh = fb.size();
  NumericVector eqli(ncoh);
  double ki = 0.0;
  for(int i=0;i<ncoh;i++) {
    ki = (-0.5365+0.00241*AET) - (-0.01586+0.000056*AET)*lignin[i];//Meentemeyer (1978)
    // Rcout<<ki<<"\n";
    eqli[i] = fb[i]/(ld[i]*ki);
  }
  eqli.attr("names") = cohortIDs(x, SpParams);
  return(eqli);
}

//' @rdname plant_values
//' @keywords internal
// [[Rcpp::export("plant_equilibriumSmallBranchLitter")]]
NumericVector cohortEquilibriumSmallBranchLitter(List x, DataFrame SpParams, double smallBranchDecompositionRate = 0.81) {
  NumericVector fu = cohortFuelLoading(x, SpParams,NA_REAL,true);
  NumericVector fb = cohortFoliarBiomass(x, SpParams, NA_REAL);
  NumericVector ld = cohortNumericParameterWithImputation(x, SpParams, "LeafDuration");
  int ncoh = fb.size();
  NumericVector eqli(ncoh);
  for(int i=0;i<ncoh;i++) {
    eqli[i] = (fu[i]-fb[i])/((ld[i]*2.0)*smallBranchDecompositionRate);
  }
  eqli.attr("names") = cohortIDs(x, SpParams);
  return(eqli);
}

//' @rdname plant_values
//' @keywords internal
// [[Rcpp::export("plant_phytovolume")]]
NumericVector cohortPhytovolume(List x, DataFrame SpParams) {
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  DataFrame shrubData = Rcpp::as<Rcpp::DataFrame>(x["shrubData"]);
  IntegerVector shrubSP;
  if((TYPEOF(shrubData["Species"]) == INTSXP) || (TYPEOF(shrubData["Species"]) == REALSXP)) {
    shrubSP = Rcpp::as<Rcpp::IntegerVector>(shrubData["Species"]);  
  } else {
    CharacterVector sspecies = Rcpp::as<Rcpp::CharacterVector>(shrubData["Species"]);
    shrubSP = speciesIndex(sspecies, SpParams);
  }
  
  NumericVector shvol = shrubPhytovolumeAllometric(shrubSP, shrubData["Cover"], shrubData["Height"], SpParams);
  NumericVector vol(treeData.nrows()+shrubData.nrows(), NA_REAL);
  for(int i=0;i<shvol.size();i++) {
    vol[i+treeData.nrows()] = shvol[i];
  }
  vol.attr("names") = cohortIDs(x, SpParams);
  return(vol);
}

//' @rdname plant_values
//' @param competitionEffect Logical flag to indicate the inclusion of competition effect on LAI estimates.
//' @keywords internal
// [[Rcpp::export("plant_LAI")]]
NumericVector cohortLAI(List x, DataFrame SpParams, double gdd = NA_REAL, 
                        bool bounded = true, bool competitionEffect = true){
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  DataFrame shrubData = Rcpp::as<Rcpp::DataFrame>(x["shrubData"]);
  NumericVector tLAI(treeData.nrows(), NA_REAL);
  NumericVector shLAI(shrubData.nrows(), NA_REAL);
  IntegerVector treeSP, shrubSP;
  
  if((TYPEOF(treeData["Species"]) == INTSXP) || (TYPEOF(treeData["Species"]) == REALSXP)) {
    treeSP = Rcpp::as<Rcpp::IntegerVector>(treeData["Species"]);
  } else {
    CharacterVector tspecies = Rcpp::as<Rcpp::CharacterVector>(treeData["Species"]);
    treeSP = speciesIndex(tspecies, SpParams);
  }
  if((TYPEOF(shrubData["Species"]) == INTSXP) || (TYPEOF(shrubData["Species"]) == REALSXP)) {
    shrubSP = Rcpp::as<Rcpp::IntegerVector>(shrubData["Species"]);  
  } else {
    CharacterVector sspecies = Rcpp::as<Rcpp::CharacterVector>(shrubData["Species"]);
    shrubSP = speciesIndex(sspecies, SpParams);
  }
  
  if(treeData.containsElementNamed("LAI")) {
    // If LAI is already in treeData, only apply phenology
    NumericVector inputTreeLAI = treeData["LAI"];
    tLAI = clone(inputTreeLAI);
  } else if(treeData.containsElementNamed("FoliarBiomass")) {
    //If FoliarBiomass is already in treeData, applySLA
    NumericVector tFB = treeData["FoliarBiomass"];
    NumericVector tSLA = speciesNumericParameterWithImputation(treeSP, SpParams, "SLA", true, true); // m2/kg (=mg/mm2)
    for(int i=0;i<tLAI.size();i++) tLAI[i] = tFB[i]*tSLA[i];
  }
  // Apply allometries to fill LAI gaps
  NumericVector tLAIAllom = treeLAIAllometric(treeSP, treeData["N"], treeData["DBH"], SpParams, NA_REAL, competitionEffect); //NO phenology yet
  for(int i=0;i<tLAI.size();i++) {
    if(NumericVector::is_na(tLAI[i])) tLAI[i] = tLAIAllom[i];
  }
  //Correct for maximum tree LAI boundaries
  double max_tree_lai = 7.0;
  double cum_tree = sum(tLAI);
  for(int i=0;i<tLAI.size();i++) {
    if(bounded && (cum_tree > max_tree_lai)) tLAI[i] = tLAI[i]*(max_tree_lai/cum_tree);
  }

  //Shrub LAI
  if(shrubData.containsElementNamed("LAI")) {
    // If LAI is already in treeData, only apply phenology
    NumericVector inputShrubLAI = shrubData["LAI"];
    shLAI = clone(inputShrubLAI);
  } else if(treeData.containsElementNamed("FoliarBiomass")) {
    //If FoliarBiomass is already in shrubData, apply phenology and SLA
    NumericVector shFB = shrubData["FoliarBiomass"];
    NumericVector shSLA = speciesNumericParameterWithImputation(shrubSP, SpParams, "SLA", true, true); // m2/kg (=mg/mm2)
    for(int i=0;i<shLAI.size();i++) {
      shLAI[i] = shFB[i]*shSLA[i];
    }
  } 
  // Apply allometries for missing values
  double treeLAI = sum(tLAI);
  NumericVector shLAIAllom = shrubLAIAllometric(shrubSP, shrubData["Cover"], shrubData["Height"], 
                                                SpParams, NA_REAL, treeLAI, competitionEffect); //NO phenology yet
  for(int i=0;i<shLAI.size();i++) {
    if(NumericVector::is_na(shLAI[i])) shLAI[i] = shLAIAllom[i];
  }
  
  //Correct for maximum shrub LAI boundaries
  double max_shrub_lai = 3.0;
  double cum_shrub = sum(shLAI);
  for(int i=0;i<shLAI.size();i++) {
    if(bounded && (cum_shrub > max_shrub_lai)) shLAI[i] = shLAI[i]*(max_shrub_lai/cum_shrub);
  }
  

  //Apply phenology to tree and shrub LAI
  if(!NumericVector::is_na(gdd)) {
    NumericVector tSgdd = speciesNumericParameterWithImputation(treeSP, SpParams, "Sgdd", true, true);
    for(int i=0;i<tLAI.size();i++) {
      tLAI[i] = tLAI[i]*leafDevelopmentStatus(tSgdd[i], gdd);
    }
    NumericVector shSgdd = speciesNumericParameterWithImputation(shrubSP, SpParams, "Sgdd", true, true);
    for(int i=0;i<shLAI.size();i++) {
      shLAI[i] = shLAI[i]*leafDevelopmentStatus(shSgdd[i], gdd);
    }
  }
  
  //Copy values
  NumericVector lai(tLAI.size()+shLAI.size());
  for(int i=0;i<tLAI.size();i++) {
    lai[i] = tLAI[i];
  }
  for(int i=0;i<shLAI.size();i++) {
    lai[i+tLAI.size()] = shLAI[i];
  }
  lai.attr("names") = cohortIDs(x, SpParams);
  return(lai);
}

double herbFoliarBiomassAllometric(double herbCover, double herbHeight, double woodyLAI){
  double herbFB = 0.014*herbCover*(herbHeight/100.0); // From piropinus
  herbFB = herbFB*exp(-0.235*woodyLAI); //Apply correction for LAI of woody elements
  if(NumericVector::is_na(herbFB)) herbFB = 0.0;
  return(herbFB);
}


double herbLAIAllometric(double herbCover, double herbHeight, double woodyLAI, double sla_herb = 9.0){
  return(std::min(2.0,herbFoliarBiomassAllometric(herbCover, herbHeight, woodyLAI)*sla_herb)); // SLA = 9 m2/kg from Brachypodium retusum in BROT2
}


//' Herbaceous description functions
//'
//' Functions to calculate attributes of the herbaceous component of a \code{\link{forest}} object 
//' @param SpParams A data frame with species parameters (see \code{\link{SpParamsMED}}).
//'  
//' @param x An object of class \code{\link{forest}}.
//' 
//' @return
//' A single scalar:
//' \itemize{
//'   \item{\code{herb_foliarBiomass}: Herbaceous biomass of leaves (in kg/m2).}
//'   \item{\code{herb_fuelLoading}: Herbaceous fine fuel loading (in kg/m2).}
//'   \item{\code{herb_LAI}: Herbaceous leaf area index (m2/m2).}
//' }
//' 
//' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
//' 
//' @seealso \code{\link{spwb}}, \code{\link{forest}}, \code{\link{plant_basalArea}}, \code{\link{summary.forest}}
//' 
//' @name herb_values
//' @keywords internal
// [[Rcpp::export("herb_foliarBiomass")]]
double herbFoliarBiomass(List x, DataFrame SpParams){
   double herbFB = NA_REAL;
   if(x.containsElementNamed("herbFoliarBiomass")) {
     herbFB = x["herbFoliarBiomass"];
   } else if(x.containsElementNamed("herbFuelLoading")) {
     herbFB = x["herbFuelLoading"];
   } else if(x.containsElementNamed("herbLAI")) {
     double herbLAI = x["herbLAI"];
     herbFB = herbLAI/9.0; //assume SLA = 9
   }
   cohortLAI(x, SpParams);
   if(NumericVector::is_na(herbFB)) {
     NumericVector LAIlive = cohortLAI(x, SpParams);
     double woodyLAI = sum(LAIlive);
     herbFB = herbFoliarBiomassAllometric(x["herbCover"], x["herbHeight"], woodyLAI);
   }
   return(herbFB);
 }

//' @rdname herb_values
//' @keywords internal
// [[Rcpp::export("herb_fuelLoading")]]
double herbFuelLoading(List x, DataFrame SpParams){
   return(herbFoliarBiomass(x, SpParams));
}

//' @rdname herb_values
//' @keywords internal
// [[Rcpp::export("herb_LAI")]]
double herbLAI(List x, DataFrame SpParams){
  double herbLAI = NA_REAL;
  if(x.containsElementNamed("herbLAI")) {
    herbLAI = x["herbLAI"];
  } 
  if(NumericVector::is_na(herbLAI)) {
    NumericVector LAIlive = cohortLAI(x, SpParams);
    double woodyLAI = sum(LAIlive);
    herbLAI = herbLAIAllometric(x["herbCover"], x["herbHeight"], woodyLAI);
  }
  return(herbLAI);
}


//' Species description functions
//'
//' Functions to calculate attributes of a \code{\link{forest}} object by species or to extract species parameters from a species parameter table (\code{\link{SpParamsMED}}).
//' 
//' @param x An object of class \code{\link{forest}}.
//' @param SpParams A data frame with species parameters (see \code{\link{SpParamsMED}}).
//' @param gdd Growth degree days (to account for leaf phenology effects).
//' @param includeDead A flag to indicate that standing dead fuels (dead branches) are included.
//' @param species A character vector of species names.
//' @param parName A string with a parameter name.
//' @param fillMissing A boolean flag to try imputation on missing values.
//' @param fillWithGenus A boolean flag to try imputation of missing values using genus values.
//' @param bounded A boolean flag to indicate that extreme values should be prevented (maximum tree LAI = 7 and maximum shrub LAI = 3)
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
//' data(exampleforest)
//' 
//' # Species basal area in the forest plot
//' species_basalArea(exampleforest, SpParamsMED)
//'   
//' # Value of parameter "Psi_Extract" for two species
//' species_parameter(c("Pinus halepensis", "Quercus ilex"), SpParamsMED, "Psi_Extract")
//'     
//' @name species_values
//' @keywords internal
// [[Rcpp::export("species_basalArea")]]
NumericVector speciesBasalArea(List x, DataFrame SpParams) {
  NumericVector cBA = cohortBasalArea(x, SpParams);
  return(sumBySpecies(cBA, cohortSpecies(x, SpParams), SpParams));
}

//' @rdname species_values
//' @keywords internal
// [[Rcpp::export("species_cover")]]
NumericVector speciesCover(List x, DataFrame SpParams) {
  NumericVector cc = cohortCover(x, SpParams);
  NumericVector sc = sumBySpecies(cc, cohortSpecies(x, SpParams), SpParams);
  for(int i=0;i<sc.length();i++) sc[i] = std::min(100.0, sc[i]);
  return(sc);
}

//' @rdname species_values
//' @keywords internal
// [[Rcpp::export("species_density")]]
NumericVector speciesDensity(List x, DataFrame SpParams) {
  NumericVector d = cohortDensity(x, SpParams);
  return(sumBySpecies(d, cohortSpecies(x, SpParams), SpParams));
}

//' @rdname species_values
//' @keywords internal
// [[Rcpp::export("species_foliarBiomass")]]
NumericVector speciesFoliarBiomass(List x, DataFrame SpParams, double gdd = NA_REAL) {
  NumericVector fb = cohortFoliarBiomass(x, SpParams, gdd);
  return(sumBySpecies(fb, cohortSpecies(x, SpParams), SpParams));
}

//' @rdname species_values
// [[Rcpp::export("species_fuelLoading")]]
NumericVector speciesFuelLoading(List x, DataFrame SpParams, double gdd = NA_REAL, bool includeDead = true) {
  NumericVector cf = cohortFuelLoading(x, SpParams, gdd, includeDead);
  return(sumBySpecies(cf, cohortSpecies(x, SpParams), SpParams));
}

//' @rdname species_values
// [[Rcpp::export("species_LAI")]]
NumericVector speciesLAI(List x, DataFrame SpParams, double gdd = NA_REAL, 
                         bool bounded = true) {
  NumericVector cl = cohortLAI(x, SpParams, gdd, bounded);
  return(sumBySpecies(cl, cohortSpecies(x, SpParams), SpParams));
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
double standFoliarBiomass(List x, DataFrame SpParams, double gdd = NA_REAL) {
  NumericVector fb = cohortFoliarBiomass(x, SpParams, gdd);
  double tfb= 0.0;
  for(int i=0;i<fb.size();i++){if(!NumericVector::is_na(fb[i])) tfb+=fb[i];}
  double hFB = herbFoliarBiomass(x, SpParams);
  tfb += hFB;
  return(tfb);
}


//' @rdname stand_values
//' @keywords internal
// [[Rcpp::export("stand_fuelLoading")]]
double standFuelLoading(List x, DataFrame SpParams, double gdd = NA_REAL, bool includeDead = true) {
  NumericVector cf = cohortFuelLoading(x, SpParams, gdd, includeDead);
  double tf= 0.0;
  for(int i=0;i<cf.size();i++){if(!NumericVector::is_na(cf[i])) tf+=cf[i];}
  double hFuel = herbFuelLoading(x, SpParams);
  tf += hFuel;
  return(tf);
}

//' @rdname stand_values
//' @keywords internal
// [[Rcpp::export("stand_shrubVolume")]]
 double standShrubVolume(List x, DataFrame SpParams) {
   NumericVector ph = cohortPhytovolume(x, SpParams);
   double sv= 0.0;
   for(int i=0;i<ph.size();i++){if(!NumericVector::is_na(ph[i])) sv+=ph[i];}
   return(sv);
 }



//' @rdname stand_values
//' @keywords internal
// [[Rcpp::export("stand_LAI")]]
double standLAI(List x, DataFrame SpParams, double gdd = NA_REAL, 
                bool bounded = true) {
  NumericVector cl = cohortLAI(x, SpParams, gdd, bounded);
  double tl= 0.0;
  for(int i=0;i<cl.size();i++){if(!NumericVector::is_na(cl[i])) tl+=cl[i];}
  double hLAI = herbLAI(x, SpParams);
  tl += hLAI;
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
NumericMatrix LAIdistribution(NumericVector z, List x, DataFrame SpParams, double gdd = NA_REAL, 
                              bool bounded = true) {
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  DataFrame shrubData = Rcpp::as<Rcpp::DataFrame>(x["shrubData"]);
  int ntree = treeData.nrows();
  int nshrub = shrubData.nrows();
  
  IntegerVector treeSP, shrubSP;
  if((TYPEOF(treeData["Species"]) == INTSXP) || (TYPEOF(treeData["Species"]) == REALSXP)) {
    treeSP = Rcpp::as<Rcpp::IntegerVector>(treeData["Species"]);
  } else {
    CharacterVector tspecies = Rcpp::as<Rcpp::CharacterVector>(treeData["Species"]);
    treeSP = speciesIndex(tspecies, SpParams);
  }
  if((TYPEOF(shrubData["Species"]) == INTSXP) || (TYPEOF(shrubData["Species"]) == REALSXP)) {
    shrubSP = Rcpp::as<Rcpp::IntegerVector>(shrubData["Species"]);  
  } else {
    CharacterVector sspecies = Rcpp::as<Rcpp::CharacterVector>(shrubData["Species"]);
    shrubSP = speciesIndex(sspecies, SpParams);
  }
  
  NumericVector treeH = treeData["Height"];
  NumericVector shrubH = shrubData["Height"];  
  
  NumericVector LAI = cohortLAI(x, SpParams, gdd, bounded);
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
NumericVector LAIprofile(NumericVector z, List x, DataFrame SpParams, double gdd = NA_REAL, 
                         bool bounded = true) {
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  DataFrame shrubData = Rcpp::as<Rcpp::DataFrame>(x["shrubData"]);
  int ntree = treeData.nrows();
  int nshrub = shrubData.nrows();
  
  IntegerVector treeSP = Rcpp::as<Rcpp::IntegerVector>(treeData["Species"]);
  NumericVector treeH = treeData["Height"];
  IntegerVector shrubSP = Rcpp::as<Rcpp::IntegerVector>(shrubData["Species"]);
  NumericVector shrubH = shrubData["Height"];  
  
  NumericVector CR = cohortCrownRatio(x, SpParams);
  NumericVector LAI = cohortLAI(x, SpParams, gdd, bounded);
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



//' Input for simulation models (deprecated)
//'
//' Functions \code{forest2spwbInput()} and \code{forest2growthInput()} take an object of class \code{\link{forest}} 
//' and a soil data input to create input objects for simulation functions \code{\link{spwb}} (or \code{\link{pwb}}) and \code{\link{growth}}, respectively. 
//' Function \code{forest2aboveground()} calculates aboveground variables such as leaf area index. 
//' Function \code{forest2belowground()} calculates belowground variables such as fine root distribution.
//' 
//' @param x An object of class \code{\link{forest}}.
//' @param SpParams A data frame with species parameters (see \code{\link{SpParamsDefinition}} and \code{\link{SpParamsMED}}).
//' @param gdd Growth degree days to account for leaf phenology effects (in Celsius). This should be left \code{NA} in most applications.
//' @param loading A logical flag to indicate that fuel loading should be included (for fire hazard calculations). 
//' @param soil An object of class \code{\link{data.frame}} or \code{\link{soil}}, containing soil parameters per soil layer.
//' @param control A list with default control parameters (see \code{\link{defaultControl}}).
//' 
//' @details
//' Function \code{forest2aboveground()} extracts height and species identity from plant cohorts of \code{x}, 
//' and calculate leaf area index and crown ratio. 
//' 
//' \emph{IMPORTANT NOTE}: Function names \code{forest2spwbInput()} and \code{forest2growthInput()} are now internal and deprecated, but 
//' they can still be used for back-compatibility. They correspond to functions \code{\link{spwbInput}} and \code{\link{growthInput}} 
//' 
//' @return 
//' Function \code{forest2aboveground()} returns a data frame with the following columns (rows are identified as specified by function \code{\link{plant_ID}}):
//' \itemize{
//'   \item{\code{SP}: Species identity (an integer) (first species is 0).}
//'   \item{\code{N}: Cohort density (ind/ha) (see function \code{\link{plant_density}}).}
//'   \item{\code{DBH}: Tree diameter at breast height (cm).}
//'   \item{\code{H}: Plant total height (cm).}
//'   \item{\code{CR}: Crown ratio (crown length to total height) (between 0 and 1).}
//'   \item{\code{LAI_live}: Live leaf area index (m2/m2) (one-side leaf area relative to plot area), includes leaves in winter dormant buds.}
//'   \item{\code{LAI_expanded}: Leaf area index of expanded leaves (m2/m2) (one-side leaf area relative to plot area).}
//'   \item{\code{LAI_dead}: Dead leaf area index (m2/m2) (one-side leaf area relative to plot area).}
//'   \item{\code{Loading}: Fine fuel loading (kg/m2), only if \code{loading = TRUE}.}
//' }
//' 
//' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
//' 
//' @seealso \code{\link{spwbInput}}, \code{\link{soil}},  
//' \code{\link{forest}}, \code{\link{SpParamsMED}}, \code{\link{defaultSoilParams}}, \code{\link{plant_ID}}
//' 
//' @examples
//' #Load example plot plant data
//' data(exampleforest)
//' 
//' #Default species parameterization
//' data(SpParamsMED)
//' 
//' # Aboveground parameters
//' forest2aboveground(exampleforest, SpParamsMED)
//' 
//' # Example of aboveground parameters taken from a forest
//' # described using LAI and crown ratio
//' data(exampleforest2)
//' forest2aboveground(exampleforest2, SpParamsMED)
//' 
//' # Define soil with default soil params (4 layers)
//' examplesoil <- defaultSoilParams(4)
//'
//' # Bewowground parameters (distribution of fine roots)
//' forest2belowground(exampleforest, examplesoil, SpParamsMED)
//' 
//' 
//' @name forest2aboveground
//' @keywords internal
// [[Rcpp::export("forest2aboveground")]]
DataFrame forest2aboveground(List x, DataFrame SpParams, double gdd = NA_REAL, bool loading = false) {
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  DataFrame shrubData = Rcpp::as<Rcpp::DataFrame>(x["shrubData"]);
  int ntree = treeData.nrows();
  int nshrub = shrubData.nrows();
  
  NumericVector CR = cohortCrownRatio(x, SpParams);
  // Maximum LAI
  NumericVector LAI_live = cohortLAI(x, SpParams, NA_REAL, true, true);
  // Actual LAI accounting for phenology
  NumericVector LAI_expanded = cohortLAI(x, SpParams, gdd, true, true);
  // Maximum LAI assuming no competition effect
  NumericVector LAI_nocomp = cohortLAI(x, SpParams, NA_REAL, true, false);
  IntegerVector SP(ntree+nshrub);
  NumericVector H(ntree+nshrub);
  
  IntegerVector treeSP, shrubSP;
  if((TYPEOF(treeData["Species"]) == INTSXP) || (TYPEOF(treeData["Species"]) == REALSXP)) {
    treeSP = Rcpp::as<Rcpp::IntegerVector>(treeData["Species"]);
  } else {
    CharacterVector tspecies = Rcpp::as<Rcpp::CharacterVector>(treeData["Species"]);
    treeSP = speciesIndex(tspecies, SpParams);
  }
  
  if((TYPEOF(shrubData["Species"]) == INTSXP) || (TYPEOF(shrubData["Species"]) == REALSXP)) {
    shrubSP = Rcpp::as<Rcpp::IntegerVector>(shrubData["Species"]);  
  } else {
    CharacterVector sspecies = Rcpp::as<Rcpp::CharacterVector>(shrubData["Species"]);
    shrubSP = speciesIndex(sspecies, SpParams);
  }
  
  NumericVector treeH = treeData["Height"];
  NumericVector treeDBH = treeData["DBH"];
  NumericVector shrubH = shrubData["Height"];  
  NumericVector shrubCover = shrubData["Cover"];  
  CharacterVector treeObsID(ntree, NA_STRING);
  CharacterVector shrubObsID(nshrub, NA_STRING);
  if(treeData.containsElementNamed("ObsID")) treeObsID = treeData["ObsID"];
  if(shrubData.containsElementNamed("ObsID")) shrubObsID = shrubData["ObsID"];
  NumericVector N = cohortDensity(x, SpParams);
    
  NumericVector LAI_dead(ntree+nshrub, 0.0);
  NumericVector DBH(ntree+nshrub, NA_REAL);
  NumericVector Cover(ntree+nshrub, NA_REAL);
  CharacterVector ObsID(ntree+nshrub, NA_STRING);
  
  for(int i=0;i<ntree;i++) {
    SP[i] = treeSP[i];
    H[i] = treeH[i];
    DBH[i] = treeDBH[i];
    LAI_dead[i] = 0.0;
    Cover[i] = NA_REAL;
    ObsID[i] = treeObsID[i];
  }
  for(int i=0;i<nshrub;i++) {
    SP[ntree+i] = shrubSP[i];
    H[ntree+i] = shrubH[i];
    DBH[ntree+i] = NA_REAL;
    LAI_dead[i] = 0.0;
    ObsID[ntree+i] = shrubObsID[i];
    Cover[ntree+i] = shrubCover[i];
  }
  DataFrame above = DataFrame::create(_["SP"]=SP, _["N"] = N,  _["DBH"] = DBH,_["Cover"] = Cover, _["H"]=H, _["CR"] = CR, 
                    _["LAI_live"]=LAI_live, _["LAI_expanded"] = LAI_expanded, _["LAI_dead"] = LAI_dead, 
                    _["LAI_nocomp"] = LAI_nocomp, _["ObsID"] = ObsID);
  if(loading) {
    NumericVector cohLoading = cohortFuelLoading(x, SpParams, gdd, true);
    above.push_back(cohLoading, "Loading");
  }
  above.attr("row.names") = cohortIDs(x, SpParams); //Assign cohort IDs to row.names
    
  return(above);
  
}


//' @rdname forest2aboveground
//' @keywords internal
// [[Rcpp::export("forest2belowground")]]
NumericMatrix forest2belowground(List x, DataFrame soil, DataFrame SpParams) {
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  DataFrame shrubData = Rcpp::as<Rcpp::DataFrame>(x["shrubData"]);
  NumericVector widths = soil["widths"];
  int nlayers = widths.size();
  NumericVector treeZ50 = treeData["Z50"];
  NumericVector treeZ95 = treeData["Z95"];
  NumericVector treeZ100(treeZ50.size(), NA_REAL);
  if(treeData.containsElementNamed("Z100")) treeZ100 = Rcpp::as<Rcpp::NumericVector>(treeData["Z100"]);
  
  NumericVector shrubZ50 = shrubData["Z50"];  
  NumericVector shrubZ95 = shrubData["Z95"];  
  NumericVector shrubZ100(shrubData.size(), NA_REAL);
  if(shrubData.containsElementNamed("Z100")) shrubZ100 = Rcpp::as<Rcpp::NumericVector>(shrubData["Z100"]);
  
  NumericMatrix V = ldrDistribution(treeZ50, shrubZ50, 
                                    treeZ95, shrubZ95, 
                                    treeZ100, shrubZ100, 
                                    widths);
  V.attr("dimnames") = List::create(cohortIDs(x, SpParams), layerNames(nlayers));
  return(V);
}
