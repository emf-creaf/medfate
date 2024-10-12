#include <Rcpp.h>
#include "communication_structures.h"
#include "tissuemoisture.h"
#include "biophysicsutils.h"
using namespace Rcpp;

const double carbonMolarMass = 12.0107; //g*mol-1
const double glucoseMolarMass = 180.156; //g*mol-1
const double starchMolarMass = 162.1406; //g*mol-1
const double starchDensity = 1.5; //g·cm-3

const double leafCperDry = 0.3; //g C · g dry-1
const double rootCperDry = 0.4959; //g C · g dry-1

const double Rn = 0.008314472; // The perfect gas constant MPa·l/K·mol = kJ/K·mol

//Sains Malaysiana 44(7)(2015): 973–977
//Modeling of Sago Starch Hydrolysis Using Glucoamylase
// double starchHydrolysis(double starchConc) {
//   double ks_hydrolysis_starch_vol = 2.324; // = g·L-1
//   double vmax_hydrolysis_starch_vol = 0.424; // = g·L-1·min-1
//   starchConc = starchConc*starchMolarMass; // from mol·L-1 to g·L-1
//   double v = (vmax_hydrolysis_starch_vol*starchConc)/(ks_hydrolysis_starch_vol + starchConc);
//   v = v/(starchMolarMass*60.0); // from g·L-1·min to mol·L-1·s-1
//   return(v); //return mol·L-1·s-1
// }

double sugarStarchDynamics(double sugarConc, double starchConc,
                           double kmsyn, double vmaxsyn, double khyd, double eqSugarConc) {
  double sugarM = sugarConc*(18.0/1000.0); //from mol·L-1 to mol·mol-1
  double starchM = starchConc*(18.0/1000.0); //from mol·L-1 to mol·mol-1
  double STsyn = vmaxsyn*(sugarM/(kmsyn + sugarM));
  double SThyd = khyd*starchM;
  double dSdt;
  if(sugarConc > eqSugarConc) {
    dSdt = STsyn*(1000.0/18.0);
  } else { //Downregulate starch synthesis
    dSdt = -SThyd*(1000.0/18.0);
  }
  return(dSdt/(3600.0*24.0)); //return mol·l-1·s-1
}

//' Carbon-related functions
//' 
//' Low-level functions used in the calculation of carbon balance.
//' 
//' @param LAI Leaf area index.
//' @param N Density (ind·ha-1).
//' @param SLA Specific leaf area (mm2/mg = m2/kg).
//' @param leafDensity  Density of leaf tissue (dry weight over volume).
//' @param SA Sapwood area (cm2).
//' @param H Plant height (cm).
//' @param L Coarse root length (mm) for each soil layer.
//' @param V Proportion of fine roots in each soil layer.
//' @param woodDensity Wood density (dry weight over volume).
//' @param conduit2sapwood Proportion of sapwood corresponding to conducive elements (vessels or tracheids) as opposed to parenchymatic tissue.
//' @param osmoticWP Osmotic water potential (MPa).
//' @param temp Temperature (degrees Celsius).
//' @param nonSugarConc Concentration of inorganic solutes (mol/l).
//' @param sugarConc Concentration of soluble sugars (mol/l).
//' @param starchConc Concentration of starch (mol/l)
//' @param eqSugarConc Equilibrium concentration of soluble sugars (mol/l).
//' 
//' @return Values returned for each function are:
//' \itemize{
//'   \item{\code{carbon_leafStarchCapacity}: Capacity of storing starch in the leaf compartment (mol gluc/ind.).}
//'   \item{\code{carbon_leafStructuralBiomass}: Leaf structural biomass (g dry/ind.)}
//'   \item{\code{carbon_sapwoodStarchCapacity}: Capacity of storing starch in the sapwood compartment (mol gluc/ind.).}
//'   \item{\code{carbon_sapwoodStructuralBiomass}: Sapwood structural biomass (g dry/ind.)}
//'   \item{\code{carbon_sapwoodStructuralLivingBiomass}: Living sapwood (parenchyma) structural biomass (g dry/ind.)}
//'   \item{\code{carbon_sugarConcentration}: Sugar concentration (mol gluc/l)}
//'   \item{\code{carbon_osmoticWaterPotential}: Osmotic component of water potential (MPa)}
//'   \item{\code{carbon_relativeSapViscosity}: Relative viscosity of sapwood with respect to pure water (according to Forst et al. (2002)).}
//'   \item{\code{carbon_sugarStarchDynamicsLeaf}: Rate of conversion from sugar to starch in leaf (mol gluc/l/s).}
//'   \item{\code{carbon_sugarStarchDynamicsStem}: Rate of conversion from sugar to starch in leaf (mol gluc/l/s).}
//'   \item{\code{carbon_carbonCompartments}: A data frame with the size of compartments for each plant cohort, in the specified units.}
//' }
//' 
//' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
//' 
//' @references
//' Forst P, Wermer F, Delgado A (2002). On the pressure dependence of the viscosity of aqueous sugar solutions. Rheol Acta 41: 369–374 DOI 10.1007/s00397-002-0238-y
//' 
//' @seealso \code{\link{growth}}
//' 
//' @name carbon
//' @keywords internal
// [[Rcpp::export("carbon_sugarStarchDynamicsLeaf")]]
double sugarStarchDynamicsLeaf(double sugarConc, double starchConc, double eqSugarConc) {
  return(sugarStarchDynamics(sugarConc, starchConc, 0.1, 0.3, 1, eqSugarConc));
}

//' @rdname carbon
//' @keywords internal
// [[Rcpp::export("carbon_sugarStarchDynamicsStem")]]
double sugarStarchDynamicsStem(double sugarConc, double starchConc, double eqSugarConc) {
  return(sugarStarchDynamics(sugarConc, starchConc, 0.1, 0.15, 0.4, eqSugarConc));
}
double sugarStarchDynamicsRoot(double sugarConc, double starchConc, double eqSugarConc) {
  return(sugarStarchDynamics(sugarConc, starchConc, 0.1, 0.6, 0.4, eqSugarConc));
}

/**
 * Van 't Hoff equation
 *  conc - mol/l 
 *  temp - deg C
 *  wp - MPa
 */
//' @rdname carbon
//' @keywords internal
// [[Rcpp::export("carbon_osmoticWaterPotential")]]
double osmoticWaterPotential(double sugarConc, double temp, double nonSugarConc) {
  sugarConc = std::max(0.0, sugarConc); // To avoid positive osmotic water potential
  return(- (sugarConc + nonSugarConc)*Rn*(temp + 273.15));
}
//' @rdname carbon
//' @keywords internal
// [[Rcpp::export("carbon_sugarConcentration")]]
double sugarConcentration(double osmoticWP, double temp, double nonSugarConc) {
  return(- osmoticWP/(Rn*(temp + 273.15)) - nonSugarConc);
}


/**
* 
*  sugarConc - sugar concentration (mol/l)
*  temp - temperature (degrees C)
*/
//' @rdname carbon
//' @keywords internal
// [[Rcpp::export("carbon_relativeSapViscosity")]]
double relativeSapViscosity(double sugarConc, double temp) {
  double x = sugarConc*glucoseMolarMass/1e3; //from mol/l to g*cm-3
  double Tkelvin = temp + 273.15;
  double q0a = 1.12; //g*cm-3
  double q1 = -0.248;
  double Ea = 2.61; //kJ*mol-1 energy of activation
  double va = x/(q0a*exp(-1.0*Ea/(Rn*Tkelvin)));
  double relVisc = exp(va/(1.0 + q1*va)); // relative viscosity
  double relWat = waterDynamicViscosity(temp); 
  return(relWat*relVisc);
}

/**
 *  Turgor (MPa)
 *  conc - mol/l 
 *  temp - deg C
 *  psi - water potential (MPa)
 */
double turgor(double psi, double sugarConc, double temp, double nonSugarConc) {
  return(std::max(0.0, psi-osmoticWaterPotential(sugarConc,temp, nonSugarConc)));
}

/**
 * Leaf area in m2 · ind-1
 */
double leafArea(double LAI, double N) {
  if(N==0.0) return(0.0);
  return(10000.0*LAI/N);
}
/**
 * leaf volume in l
 */
double leafStorageVolume(double LAI, double N, double SLA, double leafDensity) {
  return(leafArea(LAI,N)*leafWaterCapacity(SLA, leafDensity)); 
}

/*
 * Leaf structural biomass in g dw · ind-1
 */
//' @rdname carbon
//' @keywords internal
// [[Rcpp::export("carbon_leafStructuralBiomass")]]
double leafStructuralBiomass(double LAI, double N, double SLA) {
  return(1000.0*leafArea(LAI,N)/SLA);  
}

/*
 * Leaf starch storage capacity in mol · ind-1
 * Up to 10% of leaf cell volume
 */
//' @rdname carbon
//' @keywords internal
// [[Rcpp::export("carbon_leafStarchCapacity")]]
double leafStarchCapacity(double LAI, double N, double SLA, double leafDensity) {
  return(0.1*1000.0*leafStorageVolume(LAI,N,SLA,leafDensity)*starchDensity/starchMolarMass);
}


/**
 * sapwood volume in l = dm3
 * 
 * SA - cm2
 * H - cm
 * Z - mm
 */
double sapwoodVolume(double SA, double H, NumericVector L, NumericVector V) {
  int nlayers = V.size();
  double vAbove = 0.001*SA*H;
  double vBelow = 0.0;
  for(int i=0;i<nlayers;i++) {
    vBelow += 0.001*SA*V[i]*(L[i]/10.0);
  }
  return(vAbove+vBelow);
}
/**
 * sapwood storage volume in l
 *  
 *  SA - cm2
 *  H - cm
 *  Z - mm
 *  woodDensity - g/cm3
 */
double sapwoodStorageVolume(double SA, double H, NumericVector L, NumericVector V, 
                            double woodDensity, double conduit2sapwood) { 
  double woodPorosity = (1.0- (woodDensity/1.54));
  return((1.0 - conduit2sapwood)*sapwoodVolume(SA,H,L,V)*woodPorosity);
}

/**
 * sapwood structural biomass in g dw / ind
 * 
 *  SA - cm2·ind-1
 *  H - cm
 *  Z - mm
 *  woodDensity - g/cm3
 */
//' @rdname carbon
// [[Rcpp::export("carbon_sapwoodStructuralBiomass")]]
double sapwoodStructuralBiomass(double SA, double H, NumericVector L, NumericVector V, 
                                double woodDensity) {
  return(1000.0*sapwoodVolume(SA,H,L,V)*woodDensity);
}

//' @rdname carbon
//' @keywords internal
// [[Rcpp::export("carbon_sapwoodStructuralLivingBiomass")]]
double sapwoodStructuralLivingBiomass(double SA, double H, NumericVector L, NumericVector V, 
                                      double woodDensity, double conduit2sapwood) {
  return(sapwoodStructuralBiomass(SA,H,L,V,woodDensity)*(1.0-conduit2sapwood));
}

/*
 *  Sapwood starch storage capacity in mol · ind-1
 *  Up to 50% of volume of non-conductive cells
 *  
 *  SA - cm2
 *  H - cm
 *  Z - mm
 *  woodDensity - g/cm3
 */
//' @rdname carbon
//' @keywords internal
// [[Rcpp::export("carbon_sapwoodStarchCapacity")]]
double sapwoodStarchCapacity(double SA, double H, NumericVector L, NumericVector V, 
                             double woodDensity, double conduit2sapwood) {
  return(0.5*1000.0*sapwoodStorageVolume(SA,H,L,V,woodDensity,conduit2sapwood)*starchDensity/starchMolarMass);
}

void fillCarbonCompartments(DataFrame cc, List x, String biomassUnits) {
  
  if((biomassUnits!="g_m2") && (biomassUnits !="g_ind")) stop("Wrong biomass units");
  //Cohort info
  DataFrame cohorts = Rcpp::as<Rcpp::DataFrame>(x["cohorts"]);
  IntegerVector SP = Rcpp::as<Rcpp::IntegerVector>(cohorts["SP"]);
  int numCohorts = SP.size();
  
  //Aboveground parameters  
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  NumericVector H = above["H"];
  NumericVector N = above["N"];
  NumericVector LAI_expanded = above["LAI_expanded"];
  NumericVector SA = above["SA"];
  
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
  NumericVector WoodDensity = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["WoodDensity"]);
  NumericVector LeafDensity = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["LeafDensity"]);
  NumericVector conduit2sapwood = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["conduit2sapwood"]);
  
  //output vectors
  NumericVector Volume_leaves = cc["LeafStorageVolume"];
  NumericVector Volume_sapwood = cc["SapwoodStorageVolume"];
  NumericVector Starch_max_leaves = cc["LeafStarchMaximumConcentration"];
  NumericVector Starch_max_sapwood = cc["SapwoodStarchMaximumConcentration"];
  NumericVector LeafStarchCapacity = cc["LeafStarchCapacity"];
  NumericVector SapwoodStarchCapacity = cc["SapwoodStarchCapacity"];
  NumericVector leafStructBiomass = cc["LeafStructuralBiomass"];
  NumericVector sapwoodStructBiomass = cc["SapwoodStructuralBiomass"];
  NumericVector sapwoodStructLivingBiomass = cc["SapwoodLivingStructuralBiomass"];
  NumericVector fineRootBiomass = cc["FineRootBiomass"];
  NumericVector structuralBiomass = cc["StructuralBiomass"];
  NumericVector labileBiomass = cc["LabileBiomass"];
  NumericVector totalLivingBiomass = cc["TotalLivingBiomass"]; 
  NumericVector totalBiomass = cc["TotalBiomass"];
  
  for(int j=0;j<numCohorts;j++){
    fineRootBiomass[j] = fineRootBiomassIn[j];
    leafStructBiomass[j] = leafStructuralBiomass(LAI_expanded[j],N[j],SLA[j]);
    sapwoodStructBiomass[j] = sapwoodStructuralBiomass(SA[j], H[j], L(j,_),V(j,_), WoodDensity[j]);
    sapwoodStructLivingBiomass[j] = sapwoodStructuralLivingBiomass((1.0 - StemPLC[j])*SA[j], H[j], L(j,_),V(j,_), WoodDensity[j], conduit2sapwood[j]);
      
    Volume_leaves[j] = leafStorageVolume(LAI_expanded[j],  N[j], SLA[j], LeafDensity[j]);
    Volume_sapwood[j] = sapwoodStorageVolume(SA[j], H[j], L(j,_),V(j,_),WoodDensity[j], conduit2sapwood[j]);
    LeafStarchCapacity[j] = leafStarchCapacity(LAI_expanded[j],  N[j], SLA[j], LeafDensity[j]);
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
      sapwoodStructBiomass[j] = sapwoodStructBiomass[j]*f;
      sapwoodStructLivingBiomass[j] = sapwoodStructLivingBiomass[j]*f;
      fineRootBiomass[j] = fineRootBiomass[j]*f;
      labileBiomass[j] = labileBiomass[j]*f;
    }
    
    structuralBiomass[j] = leafStructBiomass[j] + sapwoodStructBiomass[j] + fineRootBiomass[j];
    totalLivingBiomass[j] = leafStructBiomass[j] + sapwoodStructLivingBiomass[j] + fineRootBiomass[j] + labileBiomass[j];
    totalBiomass[j] = leafStructBiomass[j] + sapwoodStructBiomass[j] + fineRootBiomass[j] + labileBiomass[j];
  }
}
//' @rdname carbon
//' @param x An object of class \code{\link{growthInput}}.
//' @param biomassUnits A string for output biomass units, either "g_ind" (g per individual) or "g_m2" (g per square meter).
//' @keywords internal
// [[Rcpp::export("carbon_carbonCompartments")]]
DataFrame carbonCompartments(List x, String biomassUnits = "g_m2") {
  DataFrame above = as<DataFrame>(x["above"]);
  DataFrame cc = internalCarbonCompartments(above);
  fillCarbonCompartments(cc, x, biomassUnits);
  return(cc);
}
