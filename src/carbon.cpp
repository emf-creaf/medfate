#include <Rcpp.h>
using namespace Rcpp;


const double sucroseMolarWeight = 342.3; //g*mol-1
const double starchMolarWeight = 162.1406; //g*mol-1
const double leafCperDry = 0.3; //g C · g dry-1
const double rootCperDry = 0.4959; //g C · g dry-1



NumericVector carbonCompartments(double SA, double LAI, double H, double Z, double N, double SLA, double WoodDensity, double WoodC) {
  double B_leaf = leafCperDry*1000.0*(LAI/(N/10000.0))/SLA; //Biomass in g C · ind-1
  double B_stem = WoodC*SA*(H+Z)*WoodDensity;
  double B_fineroot = B_leaf/2.5;
  return(NumericVector::create(B_leaf, B_stem, B_fineroot)); 
}