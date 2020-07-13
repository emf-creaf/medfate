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
                            double woodDensity, double vessel2sapwood) { 
  double woodPorosity = (1.0- (woodDensity/1.54));
  return((1.0 - vessel2sapwood)*sapwoodVolume(SA,H,L,V)*woodPorosity);
}

/**
 * sapwood structural biomass in g dw
 * 
 *  SA - cm2
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
                                      double woodDensity, double vessel2sapwood) {
  return(sapwoodStructuralBiomass(SA,H,L,V,woodDensity)*(1.0-vessel2sapwood));
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
                             double woodDensity, double vessel2sapwood) {
  return(0.5*1000.0*sapwoodStorageVolume(SA,H,L,V,woodDensity,vessel2sapwood)*starchDensity/starchMolarMass);
}

// NumericVector carbonCompartments(double SA, double LAI, double H, double Z, double N, double SLA, double WoodDensity, double WoodC) {
//   double B_leaf = leafStructural(LAI,N,SLA);
//   double B_stem = sapwoodCstructural(SA,H,Z,WoodDensity, WoodC);
//   double B_fineroot = B_leaf/2.5;
//   return(NumericVector::create(B_leaf, B_stem, B_fineroot)); 
// }