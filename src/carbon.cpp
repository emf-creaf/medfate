#include <Rcpp.h>
#include "hydraulics.h"
using namespace Rcpp;


const double glucoseMolarWeight = 180.156; //g*mol-1
const double starchMolarWeight = 162.1406; //g*mol-1
const double starchDensity = 1.5; //g·cm-3

const double leafCperDry = 0.3; //g C · g dry-1
const double rootCperDry = 0.4959; //g C · g dry-1

double leafArea(double LAI, double N) {
  return(10000.0*LAI/N); //Leaf area in m2 · ind-1
}

// [[Rcpp::export("carbon_leafStorageVolume")]]
double leafStorageVolume(double LAI, double N, double SLA, double leafDensity) {
  return(1000.0*leafArea(LAI,N)*leafWaterCapacity(SLA, leafDensity)); 
}

// [[Rcpp::export("carbon_leafCstructural")]]
double leafCstructural(double LAI, double N, double SLA) {
  return(leafCperDry*1000.0*leafArea(LAI,N)/SLA);  //Leaf structural biomass in g C · ind-1
}
/**
 * sapwood volume in cm3
 */
double sapwoodVolume(double SA, double H, double Z) { //SA in cm2, H and Z in cm
  return(SA*(H+Z));
}
/**
 * sapwood storage volume in cm3
 */
// [[Rcpp::export("carbon_sapwoodStorageVolume")]]
double sapwoodStorageVolume(double SA, double H, double Z, double woodDensity, double vessel2sapwood) { //SA in cm2, H and Z in cm
  double woodPorosity = (1.0- (woodDensity/1.54));
  return((1.0 - vessel2sapwood)*sapwoodVolume(SA,H,Z)*woodPorosity);
}
// [[Rcpp::export("carbon_sapwoodCstructural")]]
double sapwoodCstructural(double SA, double H, double Z, double woodDensity, double woodCperDry) {//SA in cm2, H and Z in cm
  return(woodCperDry*sapwoodVolume(SA,H,Z)*woodDensity);
}
/*
 * Leaf starch storage capacity in g C · ind-1
 * Up to 30% of leaf cell volume
 */
// [[Rcpp::export("carbon_leafStarchCapacity")]]
double leafStarchCapacity(double LAI, double N, double SLA, double leafDensity) {
  return(0.3*leafStorageVolume(LAI,N,SLA,leafDensity)*starchDensity);
}

/*
 *  Sapwood starch storage capacity in g C · ind-1
 *  Up to 50% of volume of non-conductive cells
 */
// [[Rcpp::export("carbon_sapwoodStarchCapacity")]]
double sapwoodStarchCapacity(double SA, double H, double Z, double woodDensity, double vessel2sapwood) {
  return(0.5*sapwoodStorageVolume(SA,H,Z,woodDensity,vessel2sapwood)*starchDensity);
}

NumericVector carbonCompartments(double SA, double LAI, double H, double Z, double N, double SLA, double WoodDensity, double WoodC) {
  double B_leaf = leafCstructural(LAI,N,SLA);
  double B_stem = sapwoodCstructural(SA,H,Z,WoodDensity, WoodC);
  double B_fineroot = B_leaf/2.5;
  return(NumericVector::create(B_leaf, B_stem, B_fineroot)); 
}