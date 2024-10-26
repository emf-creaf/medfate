#include <numeric>
#include <Rcpp.h>
#include "forestutils.h"
#include "paramutils.h"
#include "fuelmoisture.h"
using namespace Rcpp;

const double woodyBulkDensity = 26.43; //kg/m3
const double shortLinearBulkDensity = 26.43; //kg/m3
const double longLinearBulkDensity = 26.43; //kg/m3
const double scaleBulkDensity = 26.43; //kg/m3
const double broadleavedBulkDensity = 13.30; //kg/m3

const double defaultParticleDensity = 400.0; //kg/m3
const double defaultLowHeatContent = 18608.0; //kJ/kg

const double AET = 700; //mm (leads to k = 0.6837 for Lignin = 20)
const double smallBranchDecompositionRate = 0.3336; //year^-1 (AET = 700, Lignin = 35)

const double herbSAV = 11483.0; //m2/m3
const double woodySAV = 1601.05; //m2/m3
const double shortLinearSAV = 6562.0; //m2/m3
const double longLinearSAV = 4921.0; //m2/m3
const double scaleSAV = 2000.0; //m2/m3
const double broadleavedSAV = 8202.0; //m2/m3

const double shortLinearWmax = 0.3248; //kg/m2
const double longLinearWmax = 0.6496; //kg/m2
const double scaleWmax = 0.3248; //kg/m2
const double broadleavedWmax = 0.3472; //kg/m2

const double shortLinearReactionEfficiency = 0.18; //unitless
const double longLinearReactionEfficiency = 0.27; //unitless
const double scaleReactionEfficiency = 0.18; //unitless
const double broadleavedReactionEfficiency = 0.11; //unitless

const double shortLinearRPR = 8.03; //unitless
const double longLinearRPR = 6.35; //unitless
const double scaleRPR = 8.03; //unitless
const double broadleavedRPR = 10.31; //unitless


/**
 * Returns the proportion of the crown (Hbc - H) that lies within given interval (zLow-zHigh)
 */
double crownProportionInLayer(double zLow, double zHigh, double H, double Hbc) {
  return(leafAreaProportion(zLow, zHigh, Hbc, H));
}
/**
 * Returns the fuel loading (kg/m2) of the crown (Hbc - H) that lies within given interval (zLow-zHigh)
 */
double crownFuelInLayer(double zLow, double zHigh, double fb, double H, double Hbc) {
  return(fb*crownProportionInLayer(zLow, zHigh, H, Hbc));
}
/**
 * Returns the crown length (cm) of the crown (Hbc - H) that lies within given interval (zLow-zHigh)
 */
double crownLengthInLayer(double zLow, double zHigh, double cl, double H, double Hbc) {
  return(cl*crownProportionInLayer(zLow, zHigh, H, Hbc));
}

/**
 * Calculates fuel bulk density profile (kg/m3)
 */
NumericVector woodyFuelProfile(NumericVector z, NumericVector fuelBiomass, NumericVector H, NumericVector CR) {
  int nh = z.size();
  int ncoh = fuelBiomass.size(); //input in kg/m2
  NumericVector wfp(nh-1, 0.0);
  for(int ci=0;ci<ncoh;ci++) {
    double cbh = H[ci]*(1.0-CR[ci]);
    for(int hi=0;hi<(nh-1);hi++) {
      wfp[hi] += crownFuelInLayer(z[hi], z[hi+1], fuelBiomass[ci], H[ci], cbh)/(z[hi+1]-z[hi]);
    }
  }
  //Change units from kg/(m3*cm) to kg/m3
  wfp = wfp*100.0; // cm/m = 100; 
  return(wfp);
}
// [[Rcpp::export(".woodyFuelProfile")]]
NumericVector woodyFuelProfile(NumericVector z, List x, DataFrame SpParams, 
                               double gdd = NA_REAL) {
  NumericVector Fuel = cohortFuelLoading(x, SpParams, gdd, true); //in kg/m2
  NumericVector H = cohortHeight(x, SpParams);
  NumericVector CR = cohortCrownRatio(x, SpParams);
  return(woodyFuelProfile(z, Fuel, H, CR));
}


// [[Rcpp::export(".layerCohortFuelLoading")]]
NumericVector layerCohortFuelLoading(double minHeight, double maxHeight, NumericVector cohortLoading, NumericVector H, NumericVector CR) {
  int nCoh = cohortLoading.size();
  NumericVector l(nCoh);
  for(int i=0;i<nCoh; i++) {
    l[i] =crownFuelInLayer(minHeight, maxHeight, cohortLoading[i], H[i], H[i]*(1.0- CR[i]));
  }
  return(l);
}
// [[Rcpp::export(".layerFuelLoading")]]
double layerFuelLoading(double minHeight, double maxHeight, NumericVector cohortLoading, NumericVector H, NumericVector CR) {
  double sum = 0.0;
  int nCoh = cohortLoading.size();
  for(int i=0;i<nCoh; i++) {
    sum +=crownFuelInLayer(minHeight, maxHeight, cohortLoading[i], H[i], H[i]*(1.0-CR[i]));
  }
  return(sum);
}
// [[Rcpp::export(".layerLAI")]]
double layerLAI(double minHeight, double maxHeight, NumericVector cohortLAI, NumericVector H, NumericVector CR) {
  double sum = 0.0;
  int nCoh = cohortLAI.size();
  for(int i=0;i<nCoh; i++) {
    sum += cohortLAI[i]*crownProportionInLayer(minHeight, maxHeight, H[i], H[i]*(1.0-CR[i]));
  }
  return(sum);
}
// [[Rcpp::export(".layerFuelAverageSpeciesParameter")]]
double layerFuelAverageSpeciesParameter(String spParName, double minHeight, double maxHeight, List x, DataFrame SpParams, double gdd = NA_REAL) {
  NumericVector cohortLoading = cohortFuelLoading(x, SpParams, gdd); //in kg/m2
  NumericVector parValues = cohortNumericParameterWithImputation(x, SpParams, spParName);
  NumericVector CR = cohortCrownRatio(x, SpParams);
  NumericVector H = cohortHeight(x, SpParams);
  double num = 0.0, den = 0.0, cfl = 0.0;
  int nCoh = cohortLoading.size();
  for(int i=0;i<nCoh; i++) {
    cfl = crownFuelInLayer(minHeight, maxHeight, cohortLoading[i], H[i], H[i]*(1.0-CR[i]));
    num +=(parValues[i]*cfl);
    den += cfl;
  }
  if(den>0) return(num/den);
  return(NA_REAL);
}

// [[Rcpp::export(".layerFuelAverageParameter")]]
double layerFuelAverageParameter(double minHeight, double maxHeight, NumericVector cohortParameter, NumericVector cohortLoading, NumericVector H, NumericVector CR) {
  double num = 0.0, den = 0.0, cfl = 0.0;
  int nCoh = cohortLoading.size();
  for(int i=0;i<nCoh; i++) {
    cfl = crownFuelInLayer(minHeight, maxHeight, cohortLoading[i], H[i], H[i]*(1.0 - CR[i]));
    if(!NumericVector::is_na(cohortParameter[i])) {
      num +=(cohortParameter[i]*cfl);
      den += cfl;
    }
  }
  if(den>0) return(num/den);
  return(NA_REAL);
}
// [[Rcpp::export(".layerFuelAverageCrownLength")]]
double layerFuelAverageCrownLength(double minHeight, double maxHeight, NumericVector cohortCrownLength, NumericVector cohortLoading, NumericVector H, NumericVector CR) {
  double num = 0.0, den = 0.0, cf = 0.0, cl = 0.0;
  int nCoh = cohortLoading.size();
  for(int i=0;i<nCoh; i++) {
    cf = crownFuelInLayer(minHeight, maxHeight, cohortLoading[i], H[i], H[i]*(1.0 - CR[i]));
    cl = crownLengthInLayer(minHeight, maxHeight, cohortCrownLength[i], H[i], H[i]*(1.0 - CR[i]));
    num +=(cl*cf);
    den += cf;
  }
  if(den>0.0) return(num/den);
  return(0.0);
}

//' Fuel stratification and fuel characteristics
//' 
//' Function \code{fuel_stratification} provides a stratification of the stand into understory and canopy strata. 
//' Function \code{fuel_FCCS} calculates fuel characteristics from a \code{forest} object 
//' following an adaptation of the protocols described for the Fuel Characteristics Classification System (Prichard et al. 2013). 
//' 
//' @param object An object of class \code{\link{forest}}
//' @param SpParams A data frame with species parameters (see \code{\link{SpParamsMED}}).
//' @param cohortFMC A numeric vector of (actual) fuel moisture content by cohort.
//' @param loadingOffset A vector of length five with fine fuel loading values (canopy, shrub, herb, woody and litter) to be added to loading estimations from \code{forest}.
//' @param gdd Growth degree-days.
//' @param heightProfileStep Precision for the fuel bulk density profile.
//' @param maxHeightProfile Maximum height for the fuel bulk density profile.
//' 
//' @return 
//' Function \code{fuel_FCCS} returns a data frame with five rows corresponding to  fuel layers: 
//' \code{canopy}, \code{shrub}, \code{herb}, \code{woody} and \code{litter}. Columns correspond fuel properties:
//' \itemize{
//'   \item{\code{w}: Fine fuel loading (in kg/m2).}
//'   \item{\code{cover}: Percent cover.}
//'   \item{\code{hbc}: Height to base of crowns (in m).}
//'   \item{\code{htc}: Height to top of crowns (in m).}
//'   \item{\code{delta}: Fuel depth (in m).}
//'   \item{\code{rhob}: Fuel bulk density (in kg/m3).}
//'   \item{\code{rhop}: Fuel particle density (in kg/m3).}
//'   \item{\code{PV}: Particle volume (in m3/m2).}
//'   \item{\code{beta}: Packing ratio (unitless).}
//'   \item{\code{betarel}: Relative packing ratio (unitless).}
//'   \item{\code{etabetarel}: Reaction efficiency (unitless).}
//'   \item{\code{sigma}: Surface area-to-volume ratio (m2/m3).}
//'   \item{\code{pDead}: Proportion of dead fuels.}
//'   \item{\code{FAI}: Fuel area index (unitless).}
//'   \item{\code{h}: High heat content (in kJ/kg).}
//'   \item{\code{RV}: Reactive volume (in m3/m2).}
//'   \item{\code{MinFMC}: Minimum fuel moisture content (as percent over dry weight).}
//'   \item{\code{MaxFMC}: Maximum fuel moisture content (as percent over dry weight).}
//'   \item{\code{ActFMC}: Actual fuel moisture content (as percent over dry weight). These are set to \code{NA} if parameter \code{cohortFMC} is empty.}
//' }
//' 
//' Function \code{fuel_stratification} returns a list with the following items:
//'   \itemize{
//'     \item{\code{surfaceLayerBaseHeight}: Base height of crowns of shrubs in the surface layer (in cm).}
//'     \item{\code{surfaceLayerTopHeight}: Top height of crowns of shrubs in the surface layer (in cm).}
//'     \item{\code{understoryLAI}: Cumulated LAI of the understory layer (i.e. leaf area comprised between surface layer base and top heights).}
//'     \item{\code{canopyBaseHeight}: Base height of tree crowns in the canopy (in cm).}
//'     \item{\code{canopyTopHeight}: Top height of tree crowns in the canopy (in cm).}
//'     \item{\code{canopyLAI}: Cumulated LAI of the canopy (i.e. leaf area comprised between canopy base and top heights).}
//'   }
//'   
//' 
//' 
//' @references
//' Prichard, S. J., D. V Sandberg, R. D. Ottmar, E. Eberhardt, A. Andreu, P. Eagle, and K. Swedin. 2013. Classification System Version 3.0: Technical Documentation.
//' 
//' Reinhardt, E., D. Lutes, and J. Scott. 2006. FuelCalc: A method for estimating fuel characteristics. Pages 273–282.
//' 
//' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
//' 
//' @seealso \code{\link{fire_FCCS}}, \code{\link{spwb}}
//' 
//' @examples
//' #Load example plot plant data
//' data(exampleforest)
//' 
//' #Default species parameterization
//' data(SpParamsMED)
//' 
//' #Show stratification of fuels
//' fuel_stratification(exampleforest, SpParamsMED)
//'   
//' #Calculate fuel properties according to FCCS
//' fccs <- fuel_FCCS(exampleforest, SpParamsMED)
//' fccs
//' 
//' 
//' @name fuel_properties
// [[Rcpp::export("fuel_stratification")]]
List fuelLiveStratification(List object, DataFrame SpParams, double gdd = NA_REAL,  
                            double heightProfileStep = 10.0, double maxHeightProfile = 5000.0,double bulkDensityThreshold = 0.05) {
  int numSteps = (int) (maxHeightProfile/heightProfileStep);
  NumericVector z(numSteps);
  for(int i=0;i<numSteps; i++)  {
   if(i==0) z[i] = 0.0;
   else z[i] = z[i-1] + heightProfileStep;
  }
  //Profile has length numSteps-1
  NumericVector wfp = woodyFuelProfile(z, object, SpParams, gdd); //kg/m3
  // Rcout<<wfp.size()<< " "<< numSteps-1<<"\n";
  // for(int i=0;i<wfp.size();i++) Rcout<< i << " : "<<wfp[i]<<" - ";
  // Rcout<<"\n";
    
  int index0 = 0;
  int index0abs = 0;
  
  //Minimum height between 0 and 200 cm where BD > BDT (BD > 0 in abs)
  while((wfp[index0]<bulkDensityThreshold) && (z[index0+1]<=200.0)){index0++;}
  while((wfp[index0abs] < 0.000001) && (z[index0abs+1]<=200.0)){index0abs++;}
  
  int index1 = 0;
  int index1abs = 0;
  //Maximum height between 0 and 200 cm where BD > BDT (BD > 0 in abs)
  for(int i=0;z[i+1]<=200.0;i++) {
    if((wfp[i]>bulkDensityThreshold) && (i>index1)) {index1 = i;}
    if((wfp[i]>0.000001) && (i>index1abs)) {index1abs = i;}
  }
  //Ensure minimum height <= maximum height
  if(index0>index1) index0 = index1;
  if(index0abs>index1abs) index0abs = index1abs;
  // Rcout<<" 0: " << index0<< " "<< index0abs<<"\n";
  // Rcout<<" 1: " << index1<< " "<< index1abs<<"\n";
  
  //Minimum height above index1 where BD > BDT (BD > 0 in abs)
  int index2 = index1+1;
  int index2abs = index1abs+1;
  while((wfp[index2] < bulkDensityThreshold) && (index2<(numSteps-2))) {index2 = index2 +1;}
  while((wfp[index2abs] < 0.000001) && (index2abs<(numSteps-2))) {index2abs = index2abs +1;}
  // Rcout<< " 2: " <<index2<< " "<< index2abs<<"\n";
  
  //Maximum height where BD > BDT (BD > 0 in abs)
  int index3 = numSteps-2;
  int index3abs = numSteps-2;
  while((wfp[index3]<bulkDensityThreshold) && (index3>index2)) {index3 = index3 - 1;}  
  while((wfp[index3abs] < 0.000001) && (index3abs>index2abs)) {index3abs = index3abs - 1;}  
  //Ensure minimum height <= maximum height
  if(index2>index3) index0 = index1;
  if(index0abs>index1abs) index0abs = index1abs;
  // Rcout<<" 3: " << index3<< " "<< index3abs<<"\n";
  
  //Fuelbed height
  double fbbh = z[index0];
  double fbh = z[index1+1];
  double fbbhabs = z[index0abs];
  double fbhabs = z[index1abs+1];
  
  //Crown base and top heights (cm)
  double cbh = NA_REAL, cth = NA_REAL, cbhabs = NA_REAL, cthabs = NA_REAL;
  if(index2<(numSteps-2)) cbh = z[index2];
  if(index3<(numSteps-2)) cth = z[index3];
  // Rcout<<" cbh: " << cbh<< " "<<" cth: " << cth<< "\n";
  if(index2abs<(numSteps-2)) cbhabs = z[index2abs];
  if(index3abs<(numSteps-2)) cthabs = z[index3abs];
  // Rcout<<" cbhabs: " << cbhabs<< " "<<" cthabs: " << cthabs<< "\n";
  
  NumericVector cLAI = cohortLAI(object,SpParams, NA_REAL, true, true);
  NumericVector cH = cohortHeight(object, SpParams);
  NumericVector cCR = cohortCrownRatio(object,SpParams);
  double understoryLAI = layerLAI(fbbh, fbh, cLAI, cH, cCR);
  double canopyLAI = layerLAI(cbh, cth, cLAI, cH, cCR);
  return(List::create(_["surfaceLayerBaseHeight"] = fbbh,
                      _["surfaceLayerTopHeight"] = fbh,
                      _["surfaceLayerAbsoluteBaseHeight"] = fbbhabs,
                      _["surfaceLayerAbsoluteTopHeight"] = fbhabs,
                      _["understoryLAI"] = understoryLAI,
                      _["canopyBaseHeight"] = cbh,
                      _["canopyTopHeight"] = cth,
                      _["canopyAbsoluteBaseHeight"] = cbhabs,
                      _["canopyAbsoluteTopHeight"] = cthabs,
                      _["canopyLAI"] = canopyLAI));
}  


CharacterVector leafLitterFuelType(List object, DataFrame SpParams) {
  CharacterVector leafShape = cohortCharacterParameter(object, SpParams, "LeafShape");
  CharacterVector leafSize = cohortCharacterParameter(object, SpParams, "LeafSize");
  CharacterVector leafLitterType(leafShape.size(), NA_STRING);
  for(int i=0;i<leafShape.size();i++) {
    if(leafShape[i]=="Linear" || leafShape[i] =="Needle") {
      if(leafSize[i]=="Small") {
        leafLitterType[i] = "ShortLinear";
      } else {
        leafLitterType[i] = "LongLinear";
      }
    } else if(leafShape[i] == "Scale") {
      leafLitterType[i] = "Scale";
    } else {
      leafLitterType[i] = "Broadleaved";
    }
  }
  return(leafLitterType);
}


//' @rdname fuel_properties
//' 
//' @param bulkDensityThreshold Minimum fuel bulk density to delimit fuel strata.
//' @param depthMode Specifies how fuel depth (and therefore canopy and understory bulk density) should be estimated: 
//'   \itemize{
//'     \item{\code{"crownaverage"}: As weighed average of crown lengths using loadings as weights.}
//'     \item{\code{"profile"}: As the difference of base and top heights in bulk density profiles.}  
//'     \item{\code{"absoluteprofile"}: As the difference of absolute base and absolute top heights in bulk density profiles.}  
//'   }
// [[Rcpp::export("fuel_FCCS")]]
DataFrame FCCSproperties(List object, DataFrame SpParams, NumericVector cohortFMC = NumericVector::create(), 
                         NumericVector loadingOffset = NumericVector::create(0.0, 0.0, 0.0, 0.0, 0.0),
                         double gdd = NA_REAL,  
                         double heightProfileStep = 10.0, double maxHeightProfile = 5000, double bulkDensityThreshold = 0.05,
                         String depthMode = "crownaverage") {
  List liveStrat = fuelLiveStratification(object, SpParams, gdd, 
                                          heightProfileStep, maxHeightProfile, bulkDensityThreshold);
  
  NumericVector cc = cohortCover(object, SpParams);
  
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(object["treeData"]);
  DataFrame shrubData = Rcpp::as<Rcpp::DataFrame>(object["shrubData"]);
  int ntree = treeData.nrow();
  int nshrub = shrubData.nrow();
  double CanopyCover = 0.0;
  double ShrubCover = 0.0;
  for(int i=0;i<ntree;i++) CanopyCover = CanopyCover + cc[i];
  for(int i=0;i<nshrub;i++) ShrubCover = ShrubCover + cc[i+ntree];
  CanopyCover = std::min(100.0, CanopyCover);
  ShrubCover = std::min(100.0, ShrubCover);
  
  NumericVector cohLoading = cohortFuelLoading(object, SpParams, gdd, true);
  NumericVector cohLeafLitter = cohortEquilibriumLeafLitter(object, SpParams, AET);
  NumericVector cohSmallBranchLitter = cohortEquilibriumSmallBranchLitter(object, SpParams, smallBranchDecompositionRate);
  NumericVector cohHeight = cohortHeight(object, SpParams);
  NumericVector cohCL = cohortCrownLength(object, SpParams);
  for(int i=0;i<cohLoading.size();i++) {
    if(NumericVector::is_na(cohLoading[i])) cohLoading[i] = 0.0;
    if(NumericVector::is_na(cohLeafLitter[i])) cohLeafLitter[i] = 0.0;
    if(NumericVector::is_na(cohSmallBranchLitter[i])) cohSmallBranchLitter[i] = 0.0;
    if(NumericVector::is_na(cohHeight[i])) cohHeight[i] = 0.0;
    if(NumericVector::is_na(cohCL[i])) cohCL[i] = 0.0;
  }
  
  NumericVector cohSAV = cohortNumericParameterWithImputation(object, SpParams, "SAV", true);
  NumericVector cohHeatContent = cohortNumericParameterWithImputation(object, SpParams, "HeatContent", true);
  NumericVector cohpDead = cohortNumericParameterWithImputation(object, SpParams, "pDead");
  
  NumericVector r635 = cohortNumericParameterWithImputation(object, SpParams, "r635", true);
  NumericVector cohWoodDensity = cohortNumericParameterWithImputation(object, SpParams, "WoodDensity", true);
  NumericVector cohLeafDensity = cohortNumericParameterWithImputation(object, SpParams, "LeafDensity", true);
  NumericVector cohParticleDensity(cohWoodDensity.size(), NA_REAL);
  for(int i=0;i<cohLoading.size();i++) {
    double leaves_branches_ratio_weight = r635[i] - 1.0;
    double leaves_branches_ratio_volume = leaves_branches_ratio_weight*(cohWoodDensity[i]/cohLeafDensity[i]);
    double f_leaves_volume = 1.0/(leaves_branches_ratio_volume+1.0);
    cohParticleDensity[i] = 1000.0*(cohLeafDensity[i]*(f_leaves_volume) + cohWoodDensity[i]*(1.0 - f_leaves_volume));
  }
  NumericVector cohCR = cohortCrownRatio(object, SpParams);
  NumericVector cohMinFMC = cohortNumericParameterWithImputation(object, SpParams, "minFMC");
  NumericVector cohMaxFMC = cohortNumericParameterWithImputation(object, SpParams, "maxFMC");
  CharacterVector leafLitterType = leafLitterFuelType(object, SpParams);
  //Canopy limits and loading  
  double canopyBaseHeight = liveStrat["canopyBaseHeight"];
  double canopyTopHeight = liveStrat["canopyTopHeight"];
  double canopyAbsoluteBaseHeight = liveStrat["canopyAbsoluteBaseHeight"];
  double canopyAbsoluteTopHeight = liveStrat["canopyAbsoluteTopHeight"];
  double canopyDepth = NA_REAL;
  if(depthMode=="profile") {
    canopyDepth = (canopyTopHeight - canopyBaseHeight)/100.0;
  } else if(depthMode=="absoluteprofile") {
    canopyDepth = (canopyAbsoluteTopHeight - canopyAbsoluteBaseHeight)/100.0;
  } else if(depthMode == "crownaverage"){
    canopyDepth = layerFuelAverageCrownLength(200.0,10000.0,  cohCL, cohLoading, cohHeight, cohCR)/100.0; 
  }
  NumericVector cohCanopyLoading = layerCohortFuelLoading(200.0, 10000.0, cohLoading, cohHeight, cohCR);
  double canopyLoading = loadingOffset[0] + std::accumulate(cohCanopyLoading.begin(),cohCanopyLoading.end(),0.0);
  
  //Shrub limits and loading  
  double shrubBaseHeight = liveStrat["surfaceLayerBaseHeight"];
  double shrubTopHeight = liveStrat["surfaceLayerTopHeight"];
  double shrubAbsoluteBaseHeight = liveStrat["surfaceLayerAbsoluteBaseHeight"];
  double shrubAbsoluteTopHeight = liveStrat["surfaceLayerAbsoluteTopHeight"];
  double shrubDepth = NA_REAL;
  if(depthMode=="profile") {
    shrubDepth = (shrubTopHeight - shrubBaseHeight)/100.0;
  } else if(depthMode=="absoluteprofile") {
    shrubDepth = (shrubAbsoluteTopHeight - shrubAbsoluteBaseHeight)/100.0;
  } else if(depthMode == "crownaverage"){
    shrubDepth = layerFuelAverageCrownLength(0.0,200.0,  cohCL, cohLoading, cohHeight, cohCR)/100.0; //in m
  }

  NumericVector cohShrubLoading = layerCohortFuelLoading(0.0, 200.0, cohLoading, cohHeight, cohCR);
  double shrubLoading = loadingOffset[1] + std::accumulate(cohShrubLoading.begin(),cohShrubLoading.end(),0.0);
  //Herb limits and loading  
  double herbCover = object["herbCover"];
  if(NumericVector::is_na(herbCover)) herbCover = 0.0;
  double herbHeight = object["herbHeight"];
  if(NumericVector::is_na(herbHeight)) herbHeight = 0.0;
  double herbDepth = herbHeight/100.0; //in cm
  NumericVector LAIlive = cohortLAI(object, SpParams, NA_REAL, true, true);//Without correction
  double woodyLAI = sum(LAIlive);
  double herbLoading = loadingOffset[2] + herbFoliarBiomassAllometric(herbCover, herbHeight, woodyLAI); // From piropinus
  
  //Woody loading  
  double woodyLoading = loadingOffset[3] + std::accumulate(cohSmallBranchLitter.begin(),cohSmallBranchLitter.end(),0.0);
  //Litter loading  
  double litterLoading = loadingOffset[4] + std::accumulate(cohLeafLitter.begin(),cohLeafLitter.end(),0.0);
  
  //Properties
  NumericVector cover(5,NA_REAL); //Percent cover
  cover[0] = CanopyCover;
  cover[1] = ShrubCover;
  cover[2] = herbCover;
    
  NumericVector hbc(5,NA_REAL); //Crown base height (m)
  hbc[0] = canopyBaseHeight/100.0;
  hbc[1] = shrubBaseHeight/100.0; //in m
  hbc[2] = hbc[3] = hbc[4] = 0.0;
  NumericVector habc(5,NA_REAL); //Crown absolute base height (m)
  habc[0] = canopyAbsoluteBaseHeight/100.0;
  habc[1] = shrubAbsoluteBaseHeight/100.0; //in m
  habc[2] = habc[3] = habc[4] = 0.0;
  
  NumericVector htc(5,NA_REAL); //Crown top height (m)
  htc[0] = canopyTopHeight/100.0;//in m
  htc[1] = shrubTopHeight/100.0; //in m
  NumericVector hatc(5,NA_REAL); //Crown absolute top height (m)
  hatc[0] = canopyAbsoluteTopHeight/100.0;//in m
  hatc[1] = shrubAbsoluteTopHeight/100.0; //in m

  NumericVector rhob(5,0.0); //Bulk density (kg/m3)
  if(canopyDepth>0.0) rhob[0] = canopyLoading/canopyDepth;
  if(shrubDepth>0.0) rhob[1] = shrubLoading/shrubDepth;
  if(herbDepth>0.0) rhob[2] = herbLoading/herbDepth;
  rhob[3] = woodyBulkDensity;
  double woodyDepth = woodyLoading/woodyBulkDensity;
  int ncoh = cohLoading.size();
  rhob[4] = 0.0;
  NumericVector rhop(5,defaultParticleDensity); //Particle density (kg/m3)
  if(canopyLoading>0.0) rhop[0] = layerFuelAverageParameter(200.0, 10000.0, cohParticleDensity, cohLoading, cohHeight, cohCR);
  if(shrubLoading>0.0) rhop[1] = layerFuelAverageParameter(0.0, 200.0, cohParticleDensity, cohLoading, cohHeight, cohCR);
  if(woodyLoading>0.0) rhop[3] = 1000.0*layerFuelAverageParameter(0.0, 200.0, cohWoodDensity, cohLoading, cohHeight, cohCR);
  if(litterLoading>0.0) rhop[4] = 1000.0*layerFuelAverageParameter(0.0, 200.0, cohLeafDensity, cohLoading, cohHeight, cohCR);
    
  NumericVector pDead(5,1.0); //Proportion of dead fuels
  if(canopyLoading>0.0) pDead[0] = layerFuelAverageParameter(200.0, 10000.0, cohpDead, cohLoading, cohHeight, cohCR);
  else pDead[0] = NA_REAL;
  if(shrubLoading>0.0) pDead[1] = layerFuelAverageParameter(0.0, 200.0, cohpDead, cohLoading, cohHeight, cohCR);
  else pDead[1] = NA_REAL;
  pDead[2] = 0.0;
  
  NumericVector beta(5,0.0); //Packing ratio (unitless)
  NumericVector betarel(5,0.0); //Relative packing ratio (unitless)
  betarel[4] = 0.0;
  NumericVector etabetarel(5,0.0); //Reaction efficiency
  etabetarel[4] = 0.0;
  NumericVector PV(5,0.0); //Particle volume (m3/m2)
  PV[0] = PV[1] = PV[4] = 0.0;
  PV[2] = herbLoading/rhop[2];
  PV[3] = woodyLoading/rhop[3];
  NumericVector sigma(5,0.0); //Surface-to-area-ratio (m2/m3)
  if(canopyLoading>0.0) sigma[0] = layerFuelAverageParameter(200.0, 10000.0, cohSAV, cohLoading, cohHeight, cohCR);
  else sigma[0] = NA_REAL;
  if(shrubLoading>0.0) sigma[1] = layerFuelAverageParameter(0.0, 200.0, cohSAV, cohLoading, cohHeight, cohCR);
  else sigma[1] = NA_REAL;
  sigma[4] = 0.0;
  sigma[2] = herbSAV;
  sigma[3] = woodySAV;
  NumericVector FAI(5,0.0); //Fuel area index (m2/m2)
  FAI[0] = FAI[1] = FAI[4] = 0.0;
  NumericVector h(5,0.0); //Default low heat content (kJ/kg)
  if(canopyLoading>0.0) h[0] = layerFuelAverageParameter(200.0, 10000.0, cohHeatContent, cohLoading, cohHeight, cohCR);
  else h[0] = NA_REAL;
  if(shrubLoading>0.0) h[1] = layerFuelAverageParameter(0.0, 200.0, cohHeatContent, cohLoading, cohHeight, cohCR);
  else h[1] = NA_REAL;
  h[2] = h[3] = h[4] = defaultLowHeatContent;

  NumericVector minFMC(5,NA_REAL); //Minimum FMC
  if(canopyLoading>0.0) minFMC[0] = layerFuelAverageParameter(200.0, 10000.0, cohMinFMC, cohLoading, cohHeight, cohCR);
  else minFMC[0] = NA_REAL;
  if(shrubLoading>0.0) minFMC[1] = layerFuelAverageParameter(0.0, 200.0, cohMinFMC, cohLoading, cohHeight, cohCR);
  else minFMC[1] = NA_REAL;
  NumericVector maxFMC(5,NA_REAL); //Maximum FMC
  if(canopyLoading>0.0) maxFMC[0] = layerFuelAverageParameter(200.0, 10000.0, cohMaxFMC, cohLoading, cohHeight, cohCR);
  else maxFMC[0] = NA_REAL;
  if(shrubLoading>0.0) maxFMC[1] = layerFuelAverageParameter(0.0, 200.0, cohMaxFMC, cohLoading, cohHeight, cohCR);
  else maxFMC[1] = NA_REAL;

  NumericVector actFMC(5,NA_REAL); //Actual FMC
  if(cohortFMC.size()==cohLoading.size()) {
    if(canopyLoading>0.0) actFMC[0] = layerFuelAverageParameter(200.0, 10000.0, cohortFMC, cohLoading, cohHeight, cohCR);
    else actFMC[0] = NA_REAL;
    if(shrubLoading>0.0) actFMC[1] = layerFuelAverageParameter(0.0, 200.0, cohortFMC, cohLoading, cohHeight, cohCR);
    else actFMC[1] = NA_REAL;
  }
  NumericVector RV(5,0.0); //Reactive volume (m3/m2)
  double wmaxli = 0.0;
  double rhoblitter = 0.0, SAVlitter = 0.0, wmaxlitter = 0.0, betarellitter = 0.0, etalitter = 0.0;
  for(int i =0;i<ncoh;i++) {
    if(leafLitterType[i]=="Broadleaved")  {
      rhoblitter = broadleavedBulkDensity;
      SAVlitter = broadleavedSAV;
      wmaxlitter = broadleavedWmax;
      etalitter = broadleavedReactionEfficiency;
      betarellitter= broadleavedRPR;
    }
    else if(leafLitterType[i]=="ShortLinear")  {
      rhoblitter = shortLinearBulkDensity;
      SAVlitter = shortLinearSAV;
      wmaxlitter = shortLinearWmax;
      etalitter = shortLinearReactionEfficiency;
      betarellitter=shortLinearRPR;
    }
    else if(leafLitterType[i]=="LongLinear") {
      rhoblitter = longLinearBulkDensity; 
      SAVlitter = longLinearSAV;
      wmaxlitter = longLinearWmax;
      etalitter = longLinearReactionEfficiency;
      betarellitter=longLinearRPR;
    }
    else {
      rhoblitter = scaleBulkDensity; 
      SAVlitter = scaleSAV;
      wmaxlitter = scaleWmax;
      etalitter = scaleReactionEfficiency;
      betarellitter=scaleRPR;
    }
    etabetarel[4]+=etalitter*cohLeafLitter[i];
    betarel[4]+=betarellitter*cohLeafLitter[i];
    rhob[4] += rhoblitter *cohLeafLitter[i];
    PV[0] += cohCanopyLoading[i]/cohParticleDensity[i];
    PV[1] += cohShrubLoading[i]/cohParticleDensity[i];
    PV[4] += cohLeafLitter[i]/rhop[4];
    sigma[4] += cohLeafLitter[i]*SAVlitter;
    FAI[0] += (cohCanopyLoading[i]/cohParticleDensity[i])*cohSAV[i];
    FAI[1] += (cohShrubLoading[i]/cohParticleDensity[i])*cohSAV[i];
    FAI[4] += (cohLeafLitter[i]/rhop[4])*SAVlitter;
    wmaxli +=wmaxlitter *cohLeafLitter[i];
  }
  if(herbLoading>0.0) {
    if(herbDepth>0.0) beta[2] = PV[2]/herbDepth; 
    FAI[2] = PV[2]*sigma[2];
    if(rhop[2]>0.0) RV[2] = herbLoading/rhop[2];
  }
  if(woodyLoading>0.0) {
    if(woodyDepth>0.0) beta[3] = PV[3]/woodyDepth; 
    FAI[3] = PV[3]*sigma[3];
    if(rhop[3]) RV[3] = woodyLoading/rhop[3];
  }
  double litterDepth = 0.0;
  if(litterLoading>0.0) {
    rhob[4] = rhob[4]/litterLoading;
    if(rhob[4]>0.0) litterDepth = litterLoading/rhob[4];
    if(litterDepth>0.0) beta[4] = PV[4]/litterDepth; 
    sigma[4] = sigma[4]/litterLoading;
    etabetarel[4] = etabetarel[4]/litterLoading;
    betarel[4] = betarel[4]/litterLoading;
    wmaxli = wmaxli/litterLoading;
    if(rhob[4]>0.0) RV[4] = std::min(litterLoading, wmaxli)/rhop[4];
  }
  if(canopyLoading>0.0) {
    if(canopyDepth>0.0) beta[0] = PV[0]/canopyDepth;
    if(rhop[0]>0.0) RV[0] = canopyLoading/rhop[0];
  }
  if(shrubLoading>0.0) {
    if(shrubDepth>0.0) beta[1] = PV[1]/shrubDepth;
    if(rhop[1]>0.0) RV[1] = shrubLoading/rhop[1];
  }
  
  //Depth (m)
  NumericVector delta = NumericVector::create(canopyDepth, shrubDepth, herbDepth, woodyDepth, litterDepth);
  NumericVector w = NumericVector::create(canopyLoading, shrubLoading, herbLoading, woodyLoading, litterLoading);
  
  //Optimum depth (m). Calculations in british units and back-transformed to metric
  double delta_opt_allsurf = 0.3048*(((PV[1]+PV[2]+PV[3])* 3.2808399) +45.0*3.2808399*(RV[1]+RV[2]+RV[3]));
  double delta_opt_lowsurf = 0.3048*(((PV[2]+PV[3])* 3.2808399) +45.0*3.2808399*(RV[2]+RV[3]));
  double delta_opt_canopy = 0.3048*(0.4*FAI[0] + beta[0]*((delta[0]*3.2808399)*(cover[0]/100)));
  
  //Effective depth (m)
  double delta_eff_allsurf = ((RV[1]*delta[1])+(RV[2]*delta[2])+(RV[3]*delta[3]))/(RV[1]+RV[2]+RV[3]);
  double delta_eff_lowsurf = ((RV[2]*delta[2])+(RV[3]*delta[3]))/(RV[2]+RV[3]);
  
  //Relative packing ratio (unitless)
  if(delta[0]>0.0) betarel[0] = delta_opt_canopy/delta[0];
  if(delta_eff_allsurf>0.0) betarel[1] = delta_opt_allsurf/delta_eff_allsurf; 
  if(delta_eff_lowsurf>0.0) {
    betarel[2] = delta_opt_lowsurf/delta_eff_lowsurf; 
    betarel[3] = delta_opt_lowsurf/delta_eff_lowsurf; 
  }  
  //Reaction efficiency
  etabetarel[0] = betarel[0]*exp(1.0 - betarel[0]);
  etabetarel[1] = betarel[1]*exp(1.0 - betarel[1]);
  etabetarel[2] = betarel[2]*exp(1.0 - betarel[2]);
  etabetarel[3] = betarel[3]*exp(1.0 - betarel[3]);
  
  
  Rcpp::DataFrame FCCSProperties = DataFrame::create(_["w"] = w,
                                                     _["cover"] = cover,
                                                     _["hbc"] = hbc,
                                                     _["htc"] = htc,
                                                     _["habc"] = habc,
                                                     _["hatc"] = hatc,
                                                     _["delta"] = delta,
                                                     _["rhob"]=rhob, 
                                                     _["rhop"]=rhop,
                                                     _["PV"] = PV,
                                                     _["beta"] = beta,
                                                     _["betarel"] = betarel,
                                                     _["etabetarel"] = etabetarel,
                                                     _["sigma"]=sigma,
                                                     _["pDead"]=pDead,
                                                     _["FAI"] = FAI,
                                                     _["h"] = h,
                                                     _["RV"] = RV);
  FCCSProperties.push_back(minFMC, "MinFMC");
  FCCSProperties.push_back(maxFMC, "MaxFMC");
  FCCSProperties.push_back(actFMC,"ActFMC");
  FCCSProperties.attr("row.names") = CharacterVector::create("canopy","shrub", "herb", "woody","litter");
  return(FCCSProperties);
}

/** ROTHERMEL FUEL COMPLEX
*  Recodified from function 'ros' in package 'Rothermel' (Vacchiano & Ascoli)
*  Fuel classes are: 1-hour, 10-hour, 100-hour, live herbs and live woody
* 
*  wSI: vector of fuel load (t/ha) for five fuel classes
*  sSI: vector of surface-to-volume ratio (m2/m3) for five fuel classes
*  delta: a value of fuel bed depth (cm)
*  hSI: a vector of heat content (kJ/kg) for five fuel classes
*/
// List rothermelFuelComplex(NumericVector wSI, NumericVector sSI, double delta, NumericVector hSI) {
//   //Rescale variables to units of rothermel model
//   NumericVector w = wSI*0.02048161; //from t/ha to lbs/ft2
//   NumericVector s = sSI/3.281; //from m-1 to ft-1
//   delta = delta*0.0328084; //from cm to feet
//   NumericVector h = hSI*0.429922614; //from kJ/kg to ?
//   //Constants 
//   double rhop= 32.0; // = 513*0.0624279606 Scott and Burgan (2005)
//   
//   //Area fractions
//   NumericVector a = s*w;
//   a = a/rhop;
//   double a_dead=a[0]+a[1]+a[2];
//   double a_live=a[3]+a[4];
//   double a_tot=(a_dead+a_live);
//   NumericVector f = NumericVector::create(a[0]/a_dead,a[1]/a_dead,a[2]/a_dead,a[3]/a_live,a[4]/a_live);
//   if(a_live==0.0) {
//     f[3] = 0.0;
//     f[4] = 0.0;
//   }
//   if(a_dead==0.0) {
//     f[0] = 0.0;
//     f[1] = 0.0;
//     f[2] = 0.0;
//   }
//   double f_dead=a_dead/a_tot;
//   double f_live=a_live/a_tot;
//   if(a_tot==0.0) {
//     f_dead = 0.0;
//     f_live = 0.0;
//   }
//   
//   //weighted SAV ratio
//   double sigma_dead=f[0]*s[0]+f[1]*s[1]+f[2]*s[2];
//   double sigma_live=f[3]*s[3]+f[4]*s[4];
//   double sigma_tot=(f_dead*sigma_dead+f_live*sigma_live); //characteristic SAV
//   //weighted heat content for fuel complex
//   double h_dead=f[0]*h[0]+f[1]*h[1]+f[2]*h[2];
//   double h_live=f[3]*h[3]+f[4]*h[4];
//   //mean packing ratio for fuel complex
//   double beta=(1.0/delta)*(w[0]/rhop+w[1]/rhop+w[2]/rhop+w[3]/rhop+w[4]/rhop);
//   if(delta==0.0) beta = 0.0;
//   
//   // Andrews 2013 pag.E
//   //optimum packing ratio
//   double beta_op=3.348*pow(sigma_tot,-0.8189);
//   if(sigma_tot==0.0) beta_op = 0.0;
//   //  Rcout<<beta_op<<"\n";
//   
//   //relative packing ratio
//   double rpr=beta/beta_op; 
//   if(beta_op==0.0) rpr = 0.0;
//   
//   //for heat sink (denominator ROS)
//   double rhob=(w[0]+w[1]+w[2]+w[3]+w[4])/delta; 
//   if(delta==0.0) rhob = 0.0;
//   
//   //return values for rothermel function
//   List output=List::create(Named("Characteristic SAV [m2/m3]")=sigma_tot*3.281,
//                            Named("Bulk density [kg/m3]") = rhob*16.0184634,
//                            Named("Packing ratio [dimensionless]")=beta,
//                            Named("Relative packing ratio [dimensionless]")=rpr,
//                            Named("Dead heat content [kJ/kg]")=h_dead,
//                            Named("Live heat content [kJ/kg]")=h_live);
//   return(output);
// }

/**
 * Calculates fuel biomass (in tons/hectare) for five fuel classes:
 * Fuel classes are: 1-hour, 10-hour, 100-hour, live herbs and live woody
 */
// List fuelStructure(List object, DataFrame SpParams, DataFrame FuelModelParams, double gdd = NA_REAL, 
//                    double heightProfileStep = 10.0, double maxHeightProfile = 5000, double bulkDensityThreshold = 0.05,
//                    bool useModelForLive = false) {
//     
//   List fs = fuelStructure(object, SpParams, gdd, heightProfileStep, maxHeightProfile, bulkDensityThreshold);
// 
//     //Fuel biomass tons/ha of five classes
//   CharacterVector FMcode = object["FuelModelCode"];
//   CharacterVector models = FuelModelParams.attr("row.names");
//   NumericVector load1h = FuelModelParams["Load_1h"];
//   NumericVector load10h = FuelModelParams["Load_10h"];
//   NumericVector load100h = FuelModelParams["Load_100h"];
//   NumericVector loadLiveHerb = FuelModelParams["Load_Live_Herb"];
//   NumericVector loadLiveWoody = FuelModelParams["Load_Live_Woody"];
//   NumericVector fuelbedDepth = FuelModelParams["Fuel_Bed_Depth"];
//   NumericVector fbLoading(5), fbSurfaceToVolumeRatio(5),fbHeatContent(5); 
//   double patchsize = object["patchsize"]; //in square meters
//   double fbHeight = NA_REAL;
//   for(int i=0;i<models.length();i++) {
//     if(models[i]==FMcode[0]) {
//       fbLoading[0] = (load1h[i]*(10000.0/patchsize));
//       fbLoading[1] = (load10h[i]*(10000.0/patchsize));
//       fbLoading[2] = (load100h[i]*(10000.0/patchsize));
//       fbLoading[3] = (loadLiveHerb[i]*(10000.0/patchsize));
//       fbLoading[4] = (loadLiveWoody[i]*(10000.0/patchsize));
//       fbSurfaceToVolumeRatio[0] = Rcpp::as<Rcpp::NumericVector>(FuelModelParams["SA/V_1h"])[i];
//       fbSurfaceToVolumeRatio[1] = Rcpp::as<Rcpp::NumericVector>(FuelModelParams["SA/V_10h"])[i];
//       fbSurfaceToVolumeRatio[2] = Rcpp::as<Rcpp::NumericVector>(FuelModelParams["SA/V_100h"])[i];
//       fbSurfaceToVolumeRatio[3] = Rcpp::as<Rcpp::NumericVector>(FuelModelParams["SA/V_Live_Herb"])[i];
//       fbSurfaceToVolumeRatio[4] = Rcpp::as<Rcpp::NumericVector>(FuelModelParams["SA/V_Live_Woody"])[i];
//       fbHeatContent[0] = Rcpp::as<Rcpp::NumericVector>(FuelModelParams["Heat_1h"])[i];
//       fbHeatContent[1] = Rcpp::as<Rcpp::NumericVector>(FuelModelParams["Heat_10h"])[i];
//       fbHeatContent[2] = Rcpp::as<Rcpp::NumericVector>(FuelModelParams["Heat_100h"])[i];
//       fbHeatContent[3] = Rcpp::as<Rcpp::NumericVector>(FuelModelParams["Heat_Live_Herb"])[i];
//       fbHeatContent[4] = Rcpp::as<Rcpp::NumericVector>(FuelModelParams["Heat_Live_Woody"])[i];
//       fbHeight = fuelbedDepth[i];
//     }
//   }
//   
//   NumericVector cohortLoading = cohortFuelLoading(object, SpParams, gdd);
//   NumericVector minFMC = cohortNumericParameter(object, SpParams, "minFMC");
//   NumericVector maxFMC = cohortNumericParameter(object, SpParams, "maxFMC");
//   NumericVector CR = cohortCrownRatio(object, SpParams);
//   NumericVector H = cohortHeight(object);
//   
//   double canopyBaseHeight = fs["canopyBaseHeight"];
//   double canopyTopHeight =  fs["canopyTopHeight"];
//   
//   if(!useModelForLive) {
//     double herbCover = object["herbCover"];
//     double herbHeight = object["herbHeight"];
//     if((!NumericVector::is_na(herbCover)) && (!NumericVector::is_na(herbHeight))) { //If plot data is available
//       fbLoading[3] = 0.14*herbCover*(herbHeight/100.0); // From piropinus
//     }
//     fbLoading[4] = fs["fuelBedWoodyBiomass"];
//     fbHeight = fs["fuelbedHeight"];
//     NumericVector cohortSAV = cohortNumericParameter(object, SpParams, "SAV");
//     NumericVector cohortHeatContent = cohortNumericParameter(object, SpParams, "HeatContent");
//     //Herbaceous fuel
//     fbSurfaceToVolumeRatio[3] = 4921.0;
//     fbHeatContent[3] = 18622;
//     //Woody fuel
//     fbSurfaceToVolumeRatio[4] = layerFuelAverageParameter( 0.0, fbHeight, cohortSAV, cohortLoading, H, CR);
//     fbHeatContent[4] = layerFuelAverageParameter( 0.0, fbHeight, cohortHeatContent, cohortLoading, H, CR);
//   }
//   
//   
//   //Rothermel fuel complex
//   List fc = rothermelFuelComplex(fbLoading, fbSurfaceToVolumeRatio, fbHeight, fbHeatContent);
//   
//   //Fuel minimum and maximum moisture
//   double minCanopyFMC = canopyLiveFuelMoisture( canopyBaseHeight, canopyTopHeight, minFMC, 
//                                            cohortLoading, H, CR);
//   double maxCanopyFMC = canopyLiveFuelMoisture( canopyBaseHeight, canopyTopHeight, maxFMC, 
//                                            cohortLoading, H, CR);
//   double minFuelBedFMC = fuelbedLiveFuelMoisture( fbHeight, minFMC, 
//                                            cohortLoading, H, CR);
//   double maxFuelBedFMC = fuelbedLiveFuelMoisture( fbHeight, maxFMC, 
//                                            cohortLoading, H, CR);
//   
//   return(List::create(_["fuelbedHeight [cm]"] = fbHeight,
//                       _["fuelbedLoading [kg/m2]"] = fbLoading,
//                       _["fuelbedSAV [m2/m3]"] = fbSurfaceToVolumeRatio,
//                       _["fuelbedHeatContent [kJ/kg]"] = fbHeatContent,
//                       _["fuelcomplexSAV [m2/m3]"] = fc["Characteristic SAV [m2/m3]"],
//                       _["fuelcomplexBulkDensity [kg/m3]"] = fc["Bulk density [kg/m3]"],
//                       _["fuelcomplexPackingRatio [dimensionless]"] = fc["Packing ratio [dimensionless]"],
//                       _["fuelcomplexRelativePackingRatio [dimensionless]"] = fc["Relative packing ratio [dimensionless]"],
//                       _["fuelcomplexDeadHeatContent [kJ/kg]"] = fc["Dead heat content [kJ/kg]"],
//                       _["fuelcomplexLiveHeatContent [kJ/kg]"] = fc["Live heat content [kJ/kg]"],
//                       _["canopyBaseHeight [cm]"] = canopyBaseHeight,
//                       _["canopyTopHeight [cm]"] = canopyTopHeight,
//                       _["canopyLength [cm]"] = fs["canopyLength"],
//                       _["canopyBulkDensity [kg/m3]"] = fs["canopyBulkDensity"],
//                       _["canopyLAI [dimensionless]"] = fs["canopyLAI"],
//                       _["canopyMinFMC [%]"] = minCanopyFMC,
//                       _["canopyMaxFMC [%]"] = maxCanopyFMC,
//                       _["fuelbedMinFMC [%]"] = minFuelBedFMC,
//                       _["fuelbedMaxFMC [%]"] = maxFuelBedFMC));
// }
