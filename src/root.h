#include <Rcpp.h>

#ifndef ROOT_H
#define ROOT_H
#endif
using namespace Rcpp;


NumericVector ldrRS_one(double Z50, double Z95, NumericVector d);
NumericVector conicRS_one(double Z, NumericVector d);
NumericMatrix conicDistribution(NumericVector Z, NumericVector d);
NumericMatrix ldrDistribution(NumericVector Z50, NumericVector Z95, NumericVector d);

double fineRootRadius(double specificRootLength, double rootTissueDensity);
double specificRootSurfaceArea(double specificRootLength, double rootTissueDensity);
double fineRootAreaIndex(NumericVector Ksoil, NumericVector krhizo, double lai,
                         double specificRootLength, double rootTissueDensity,  
                         double rootLengthDensity = 10.0);
double fineRootBiomassPerIndividual(NumericVector Ksoil, NumericVector krhizo,  double lai, double N,
                                    double specificRootLength, double rootTissueDensity,  
                                    double rootLengthDensity = 10.0);
double fineRootSoilVolume(double fineRootBiomass, double specificRootLength, double rootLengthDensity = 10.0);

double coarseRootSoilVolume(double Kmax_rootxylem, double VCroot_kmax, double Al2As,
                            NumericVector V, NumericVector d, NumericVector rfc);

NumericVector coarseRootLengths(double VolInd, NumericVector V, NumericVector d, NumericVector rfc);

List horizontalProportionsBasic(NumericVector poolProportions, NumericMatrix V, 
                                double LAIcell, double poolOverlapFactor);
List horizontalProportionsAdvanced(NumericVector poolProportions, NumericVector VolInd, 
                                   NumericVector N, NumericMatrix V, 
                                   NumericVector d, NumericVector rfc);