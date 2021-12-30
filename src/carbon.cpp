#include <Rcpp.h>
#include "hydraulics.h"
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

// [[Rcpp::export("carbon_sugarStarchDynamicsLeaf")]]
double sugarStarchDynamicsLeaf(double sugarConc, double starchConc, double eqSugarConc) {
  return(sugarStarchDynamics(sugarConc, starchConc, 0.1, 0.3, 1, eqSugarConc));
}
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
// [[Rcpp::export("carbon_osmoticWaterPotential")]]
double osmoticWaterPotential(double sugarConc, double temp, double nonSugarConc) {
  return(- (sugarConc + nonSugarConc)*Rn*(temp + 273.15));
}
// [[Rcpp::export("carbon_sugarConcentration")]]
double sugarConcentration(double osmoticWP, double temp, double nonSugarConc) {
  return(- osmoticWP/(Rn*(temp + 273.15)) - nonSugarConc);
}


/**
* On the pressure dependence of the viscosity of aqueous sugar solutions
* Rheol Acta (2002) 41: 369–374 DOI 10.1007/s00397-002-0238-y
* 
*  sugarConc - sugar concentration (mol/l)
*  temp - temperature (degrees C)
*/
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
// [[Rcpp::export("carbon_leafStructuralBiomass")]]
double leafStructuralBiomass(double LAI, double N, double SLA) {
  return(1000.0*leafArea(LAI,N)/SLA);  
}

/*
 * Leaf starch storage capacity in mol · ind-1
 * Up to 10% of leaf cell volume
 */
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
// [[Rcpp::export("carbon_sapwoodStructuralBiomass")]]
double sapwoodStructuralBiomass(double SA, double H, NumericVector L, NumericVector V, 
                                double woodDensity) {
  return(1000.0*sapwoodVolume(SA,H,L,V)*woodDensity);
}

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
// [[Rcpp::export("carbon_sapwoodStarchCapacity")]]
double sapwoodStarchCapacity(double SA, double H, NumericVector L, NumericVector V, 
                             double woodDensity, double conduit2sapwood) {
  return(0.5*1000.0*sapwoodStorageVolume(SA,H,L,V,woodDensity,conduit2sapwood)*starchDensity/starchMolarMass);
}

// [[Rcpp::export("carbon_carbonCompartments")]]
DataFrame carbonCompartments(List x, String units = "kg_m2") {
  
  if((units!="kg_m2") & (units !="g_ind")) stop("Wrong units");
  //Cohort info
  DataFrame cohorts = Rcpp::as<Rcpp::DataFrame>(x["cohorts"]);
  IntegerVector SP = Rcpp::as<Rcpp::IntegerVector>(cohorts["SP"]);
  int numCohorts = SP.size();
  
  //Aboveground parameters  
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  NumericVector DBH = above["DBH"];
  NumericVector Cover = above["Cover"];
  NumericVector H = above["H"];
  NumericVector N = above["N"];
  NumericVector CR = above["CR"];
  NumericVector LAI_live = above["LAI_live"];
  NumericVector LAI_expanded = above["LAI_expanded"];
  NumericVector LAI_dead = above["LAI_dead"];
  NumericVector SA = above["SA"];
  
  DataFrame belowdf = Rcpp::as<Rcpp::DataFrame>(x["below"]);
  NumericVector fineRootBiomass(numCohorts,NA_REAL);
  if(belowdf.containsElementNamed("fineRootBiomass")) fineRootBiomass = clone(Rcpp::as<Rcpp::NumericVector>(belowdf["fineRootBiomass"]));
  List belowLayers = Rcpp::as<Rcpp::List>(x["belowLayers"]);
  NumericMatrix V = Rcpp::as<Rcpp::NumericMatrix>(belowLayers["V"]);
  NumericMatrix L = Rcpp::as<Rcpp::NumericMatrix>(belowLayers["L"]);
  
  DataFrame internalCarbon = Rcpp::as<Rcpp::DataFrame>(x["internalCarbon"]);
  NumericVector sugarLeaf = internalCarbon["sugarLeaf"]; //Concentrations assuming RWC = 1
  NumericVector starchLeaf = internalCarbon["starchLeaf"];
  NumericVector sugarSapwood = internalCarbon["sugarSapwood"];
  NumericVector starchSapwood = internalCarbon["starchSapwood"];
  
  //Anatomy parameters
  DataFrame paramsAnatomy = Rcpp::as<Rcpp::DataFrame>(x["paramsAnatomy"]);
  NumericVector SLA = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["SLA"]);
  NumericVector Al2As = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["Al2As"]);
  NumericVector WoodDensity = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["WoodDensity"]);
  NumericVector LeafDensity = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["LeafDensity"]);
  NumericVector conduit2sapwood = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["conduit2sapwood"]);
  
  //Storage volume and maximum starch capacity for leaves and sapwood  
  NumericVector Volume_leaves(numCohorts,0.0);
  NumericVector Volume_sapwood(numCohorts,0.0);
  
  NumericVector LeafStructBiomass(numCohorts,0.0);
  NumericVector SapwoodStructBiomass(numCohorts,0.0);
  NumericVector LabileBiomass(numCohorts, 0.0);
  NumericVector TotalBiomass(numCohorts, 0.0);
  NumericVector TotalBiomass_kgm2(numCohorts, 0.0);
  
  for(int j=0;j<numCohorts;j++){
    LeafStructBiomass[j] = leafStructuralBiomass(LAI_expanded[j],N[j],SLA[j]);
    SapwoodStructBiomass[j] = sapwoodStructuralBiomass(SA[j], H[j], L(j,_),V(j,_), WoodDensity[j]);
    Volume_leaves[j] = leafStorageVolume(LAI_expanded[j],  N[j], SLA[j], LeafDensity[j]);
    Volume_sapwood[j] = sapwoodStorageVolume(SA[j], H[j], L(j,_),V(j,_),WoodDensity[j], conduit2sapwood[j]);
    double labileMassLeaf = (sugarLeaf[j]+starchLeaf[j])*(glucoseMolarMass*Volume_leaves[j]);
    double labileMassSapwood = (sugarSapwood[j]+starchSapwood[j])*(glucoseMolarMass*Volume_sapwood[j]);
    if(NumericVector::is_na(fineRootBiomass[j]))  fineRootBiomass[j] = LeafStructBiomass[j]/2.0; //If missing
    LabileBiomass[j] = labileMassSapwood+labileMassLeaf;
    TotalBiomass[j] = LeafStructBiomass[j] + SapwoodStructBiomass[j] + fineRootBiomass[j]+ LabileBiomass[j];
  }
  if(units=="kg_m2") {
    for(int j=0;j<numCohorts;j++){
      double f = N[j]/(10000.0*1000.0);
      LeafStructBiomass[j] = LeafStructBiomass[j]*f;
      SapwoodStructBiomass[j] = SapwoodStructBiomass[j]*f;
      fineRootBiomass[j] = fineRootBiomass[j]*f;
      LabileBiomass[j] = LabileBiomass[j]*f;
      TotalBiomass[j] = TotalBiomass[j]*f;
    }
  }
  
  DataFrame df = DataFrame::create(
    _["LeafStructBiomass"] = LeafStructBiomass,
    _["SapwoodStructBiomass"] = SapwoodStructBiomass,
    _["FineRootBiomass"] = fineRootBiomass,
    _["LabileBiomass"] = LabileBiomass,
    _["TotalBiomass"] = TotalBiomass
  );
  df.attr("row.names") = above.attr("row.names");
  return(df);
}