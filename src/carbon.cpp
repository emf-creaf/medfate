#include <RcppArmadillo.h>
#include "communication_structures.h"
#include "carbon.h"
#include "forestutils.h"
#include "carbon_c.h"
#include "biophysicsutils.h"
#include "biophysicsutils_c.h"
#include "meteoland/utils_c.hpp"
using namespace Rcpp;

double sapwoodStorageVolume(double SA, double H, NumericVector L, NumericVector V, 
                              double woodDensity, double conduit2sapwood) { 
  return(sapwoodStorageVolume_c(SA, H,
                                Rcpp::as<std::vector<double>>(L),
                                Rcpp::as<std::vector<double>>(V),
                                woodDensity, conduit2sapwood));
}

//' @rdname carbon
//' @keywords internal
// [[Rcpp::export("carbon_belowgroundSapwoodStructuralBiomass")]]
double belowgroundSapwoodStructuralBiomass(double SA, NumericVector L, NumericVector V, 
                                             double woodDensity) {
  return(belowgroundSapwoodStructuralBiomass_c(SA, 
                                               Rcpp::as<std::vector<double>>(L),
                                               Rcpp::as<std::vector<double>>(V),
                                               woodDensity));
}

//' @rdname carbon
//' @keywords internal
// [[Rcpp::export("carbon_sapwoodStructuralBiomass")]]
double sapwoodStructuralBiomass(double SA, double H, NumericVector L, NumericVector V, 
                                double woodDensity) {
  return(sapwoodStructuralBiomass_c(SA, H,
                                    Rcpp::as<std::vector<double>>(L),
                                    Rcpp::as<std::vector<double>>(V),
                                    woodDensity));
}

//' @rdname carbon
//' @keywords internal
// [[Rcpp::export("carbon_belowgroundHeartwoodStructuralBiomass")]]
double belowgroundHeartwoodStructuralBiomass(double DBH, double SA, NumericVector L, NumericVector V, 
                                             double woodDensity) {
  return(belowgroundHeartwoodStructuralBiomass_c(DBH, SA,
                                                 Rcpp::as<std::vector<double>>(L),
                                                 Rcpp::as<std::vector<double>>(V),
                                                 woodDensity));
}

//' @rdname carbon
//' @keywords internal
// [[Rcpp::export("carbon_heartwoodStructuralBiomass")]]
double heartwoodStructuralBiomass(double DBH, double SA, double H, NumericVector L, NumericVector V, 
                                  double woodDensity) {
  return(heartwoodStructuralBiomass_c(DBH, SA, H,
                                      Rcpp::as<std::vector<double>>(L),
                                      Rcpp::as<std::vector<double>>(V),
                                      woodDensity));
}

//' @rdname carbon
//' @keywords internal
// [[Rcpp::export("carbon_sapwoodStarchCapacity")]]
double sapwoodStarchCapacity(double SA, double H, NumericVector L, NumericVector V, 
                             double woodDensity, double conduit2sapwood) {
  return(sapwoodStarchCapacity_c(SA, H,
                                 Rcpp::as<std::vector<double>>(L),
                                 Rcpp::as<std::vector<double>>(V),
                                 woodDensity, conduit2sapwood));
}

void fillCarbonCompartments(DataFrame cc, List x, String biomassUnits) {
  
  if((biomassUnits!="g_m2") && (biomassUnits !="g_ind") && (biomassUnits !="gC_m2")) stop("Wrong biomass units");
  //Aboveground parameters  
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  NumericVector H = above["H"];
  NumericVector N = above["N"];
  NumericVector LAI_expanded = above["LAI_expanded"];
  NumericVector LAI_dead = above["LAI_dead"];
  NumericVector SA = above["SA"];
  NumericVector DBH = above["DBH"];
  int numCohorts = above.nrow();
  CharacterVector ctype = cohortType(above.attr("row.names"));
  
  DataFrame belowdf = Rcpp::as<Rcpp::DataFrame>(x["below"]);
  NumericVector fineRootBiomassIn = Rcpp::as<Rcpp::NumericVector>(belowdf["fineRootBiomass"]);
  List belowLayers = Rcpp::as<Rcpp::List>(x["belowLayers"]);
  NumericMatrix V = Rcpp::as<Rcpp::NumericMatrix>(belowLayers["V"]);
  NumericMatrix L = Rcpp::as<Rcpp::NumericMatrix>(belowLayers["L"]);
  
  DataFrame internalCarbon = Rcpp::as<Rcpp::DataFrame>(x["internalCarbon"]);
  NumericVector sugarLeaf = internalCarbon["sugarLeaf"]; //Concentrations assuming RWC = 1
  NumericVector starchLeaf = internalCarbon["starchLeaf"];
  NumericVector sugarSapwood = internalCarbon["sugarSapwood"];
  NumericVector starchSapwood = internalCarbon["starchSapwood"];
  
  //Internal state variables
  DataFrame internalWater = Rcpp::as<Rcpp::DataFrame>(x["internalWater"]);
  //Values at the end of the day (after calling spwb)
  NumericVector StemPLC = Rcpp::as<Rcpp::NumericVector>(internalWater["StemPLC"]);
  
  //Anatomy parameters
  DataFrame paramsAnatomy = Rcpp::as<Rcpp::DataFrame>(x["paramsAnatomy"]);
  NumericVector SLA = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["SLA"]);
  NumericVector Al2As = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["Al2As"]);
  NumericVector r635 = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["r635"]);
  NumericVector WoodDensity = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["WoodDensity"]);
  NumericVector LeafDensity = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["LeafDensity"]);
  NumericVector conduit2sapwood = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["conduit2sapwood"]);
  
  //Anatomy parameters
  DataFrame paramsGrowth = Rcpp::as<Rcpp::DataFrame>(x["paramsGrowth"]);
  NumericVector WoodC = Rcpp::as<Rcpp::NumericVector>(paramsGrowth["WoodC"]);
  
  //output vectors
  NumericVector Volume_leaves = cc["LeafStorageVolume"];
  NumericVector Volume_sapwood = cc["SapwoodStorageVolume"];
  NumericVector Starch_max_leaves = cc["LeafStarchMaximumConcentration"];
  NumericVector Starch_max_sapwood = cc["SapwoodStarchMaximumConcentration"];
  NumericVector LeafStarchCapacity = cc["LeafStarchCapacity"];
  NumericVector SapwoodStarchCapacity = cc["SapwoodStarchCapacity"];
  
  NumericVector leafStructBiomass = cc["LeafStructuralBiomass"];
  NumericVector twigStructBiomass = cc["TwigStructuralBiomass"];
  NumericVector twigStructLivingBiomass = cc["TwigLivingStructuralBiomass"];
  NumericVector deadLeafStructBiomass = cc["DeadLeafStructuralBiomass"];
  NumericVector deadTwigStructBiomass = cc["DeadTwigStructuralBiomass"];
  NumericVector sapwoodStructBiomass = cc["SapwoodStructuralBiomass"];
  NumericVector sapwoodStructLivingBiomass = cc["SapwoodLivingStructuralBiomass"];
  NumericVector heartwoodStructBiomass = cc["HeartwoodStructuralBiomass"];
  NumericVector abovegroundWoodBiomass = cc["AbovegroundWoodBiomass"];
  NumericVector belowgroundWoodBiomass = cc["BelowgroundWoodBiomass"];
  NumericVector fineRootBiomass = cc["FineRootBiomass"];
  NumericVector labileBiomass = cc["LabileBiomass"];
  
  NumericVector structuralBiomass = cc["StructuralBiomass"];
  NumericVector deadBiomass = cc["DeadBiomass"];
  NumericVector totalLivingBiomass = cc["TotalLivingBiomass"]; 
  NumericVector totalBiomass = cc["TotalBiomass"];
  
  for(int j=0;j<numCohorts;j++){
    fineRootBiomass[j] = fineRootBiomassIn[j];
    leafStructBiomass[j] = leafStructuralBiomass_c(LAI_expanded[j],N[j],SLA[j]);
    twigStructBiomass[j] = twigStructuralBiomass_c(LAI_expanded[j],N[j],SLA[j], r635[j]);
    deadLeafStructBiomass[j] = leafStructuralBiomass_c(LAI_dead[j],N[j],SLA[j]);
    deadTwigStructBiomass[j] = twigStructuralBiomass_c(LAI_dead[j],N[j],SLA[j], r635[j]);
    twigStructLivingBiomass[j] = twigStructBiomass[j]*(1.0 - StemPLC[j])*(1.0-conduit2sapwood[j]);
    sapwoodStructBiomass[j] = sapwoodStructuralBiomass(SA[j], H[j], L(j,_),V(j,_), WoodDensity[j]);
    sapwoodStructLivingBiomass[j] = sapwoodStructBiomass[j]*(1.0 - StemPLC[j])*(1.0-conduit2sapwood[j]);
    abovegroundWoodBiomass[j] = abovegroundSapwoodStructuralBiomass_c(SA[j], H[j], WoodDensity[j]);
    belowgroundWoodBiomass[j] = belowgroundSapwoodStructuralBiomass(SA[j], L(j,_), V(j,_), WoodDensity[j]);
    if(ctype[j] =="tree") {
      heartwoodStructBiomass[j] = heartwoodStructuralBiomass(DBH[j], SA[j], H[j], L(j,_),V(j,_), WoodDensity[j]);
      abovegroundWoodBiomass[j] += abovegroundHeartwoodStructuralBiomass_c(DBH[j], SA[j], H[j], WoodDensity[j]);
      belowgroundWoodBiomass[j] += belowgroundHeartwoodStructuralBiomass(DBH[j], SA[j], L(j,_),V(j,_), WoodDensity[j]);
    } else {
      heartwoodStructBiomass[j] = 0.0;
    }
    Volume_leaves[j] = leafStorageVolume_c(LAI_expanded[j],  N[j], SLA[j], LeafDensity[j]);
    Volume_sapwood[j] = sapwoodStorageVolume(SA[j], H[j], L(j,_),V(j,_),WoodDensity[j], conduit2sapwood[j]);
    LeafStarchCapacity[j] = leafStarchCapacity_c(LAI_expanded[j],  N[j], SLA[j], LeafDensity[j]);
    SapwoodStarchCapacity[j] =  sapwoodStarchCapacity(SA[j], H[j], L(j,_), V(j,_),WoodDensity[j], conduit2sapwood[j]);
    Starch_max_leaves[j] = LeafStarchCapacity[j]/Volume_leaves[j];
    Starch_max_sapwood[j] = SapwoodStarchCapacity[j]/Volume_sapwood[j];
    if(Volume_leaves[j]==0.0) Starch_max_leaves[j] = 0.0;
    if(Volume_sapwood[j]==0.0) Starch_max_sapwood[j] = 0.0;
    
    double labileMassLeaf = (sugarLeaf[j]+starchLeaf[j])*(glucoseMolarMass*Volume_leaves[j]);
    double labileMassSapwood = (sugarSapwood[j]+starchSapwood[j])*(glucoseMolarMass*Volume_sapwood[j]);
    labileBiomass[j] = labileMassSapwood+labileMassLeaf;
    
    if(biomassUnits=="g_m2") {
      double f = N[j]/(10000.0);
      leafStructBiomass[j] = leafStructBiomass[j]*f;
      twigStructBiomass[j] = twigStructBiomass[j]*f;
      twigStructLivingBiomass[j] = twigStructLivingBiomass[j]*f;
      deadLeafStructBiomass[j] = deadLeafStructBiomass[j]*f;
      deadTwigStructBiomass[j] = deadTwigStructBiomass[j]*f;
      sapwoodStructBiomass[j] = sapwoodStructBiomass[j]*f;
      sapwoodStructLivingBiomass[j] = sapwoodStructLivingBiomass[j]*f;
      heartwoodStructBiomass[j] = heartwoodStructBiomass[j]*f;
      abovegroundWoodBiomass[j] = abovegroundWoodBiomass[j]*f;
      belowgroundWoodBiomass[j] = belowgroundWoodBiomass[j]*f;
      fineRootBiomass[j] = fineRootBiomass[j]*f;
      labileBiomass[j] = labileBiomass[j]*f;
    } else if(biomassUnits=="gC_m2") {
      double f = N[j]/(10000.0);
      leafStructBiomass[j] = leafStructBiomass[j]*f*leafCperDry;
      twigStructBiomass[j] = twigStructBiomass[j]*f*WoodC[j];
      twigStructLivingBiomass[j] = twigStructLivingBiomass[j]*f*WoodC[j];
      deadLeafStructBiomass[j] = deadLeafStructBiomass[j]*f*leafCperDry;
      deadTwigStructBiomass[j] = deadTwigStructBiomass[j]*f*WoodC[j];
      sapwoodStructBiomass[j] = sapwoodStructBiomass[j]*f*WoodC[j];
      sapwoodStructLivingBiomass[j] = sapwoodStructLivingBiomass[j]*f*WoodC[j];
      heartwoodStructBiomass[j] = heartwoodStructBiomass[j]*f*WoodC[j];
      abovegroundWoodBiomass[j] = abovegroundWoodBiomass[j]*f*WoodC[j];
      belowgroundWoodBiomass[j] = belowgroundWoodBiomass[j]*f*WoodC[j];
      fineRootBiomass[j] = fineRootBiomass[j]*f*rootCperDry;
      labileBiomass[j] = labileBiomass[j]*f*(6.0*carbonMolarMass/glucoseMolarMass);
    }
    
    structuralBiomass[j] = leafStructBiomass[j]  + twigStructBiomass[j]   + sapwoodStructBiomass[j]  + fineRootBiomass[j];
    deadBiomass[j] = deadLeafStructBiomass[j] + deadTwigStructBiomass[j] + heartwoodStructBiomass[j];
    totalLivingBiomass[j] = leafStructBiomass[j] + twigStructLivingBiomass[j] + sapwoodStructLivingBiomass[j] + fineRootBiomass[j] + labileBiomass[j];
    totalBiomass[j] = structuralBiomass[j] + deadBiomass[j] + labileBiomass[j];
  }
  // if(biomassUnits=="gC_m2"){
  //   Rcout << " B leaves " << sum(leafStructBiomass) << " fineroots "<< sum(fineRootBiomass) <<" sapwood " << sum(sapwoodStructBiomass) << " heartwood "<< sum(heartwoodStructBiomass)<<" struct " << sum(structuralBiomass)<< " dead " << sum(deadBiomass) << " labile " << sum(labileBiomass) << " total " << sum(totalBiomass) <<"\n";
  // }
}

//' @rdname carbon
//' @param x An object of class \code{\link{growthInput}}.
//' @param biomassUnits A string for output biomass units, either "g_ind" (g per individual) or "g_m2" (g per square meter).
//' @keywords internal
// [[Rcpp::export("carbon_carbonCompartments")]]
DataFrame carbonCompartments(List x, String biomassUnits = "g_m2") {
  DataFrame above = as<DataFrame>(x["above"]);
  int numCohorts = above.nrow();
  DataFrame cc = communicationCarbonCompartments(numCohorts);
  fillCarbonCompartments(cc, x, biomassUnits);
  return(cc);
}
