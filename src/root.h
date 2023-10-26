#include <Rcpp.h>

#ifndef ROOT_H
#define ROOT_H
#endif
using namespace Rcpp;


NumericVector ldrRS_one(double Z50, double Z95, NumericVector d);
NumericVector conicRS_one(double Z, NumericVector d);
NumericMatrix conicDistribution(NumericVector Z, NumericVector d);
NumericMatrix ldrDistribution(NumericVector Z50, NumericVector Z95, NumericVector d);
NumericMatrix ldrDistribution(NumericVector treeZ50, NumericVector shrubZ50, 
                              NumericVector treeZ95, NumericVector shrubZ95, NumericVector d);

double fineRootRadius(double specificRootLength, double rootTissueDensity);
double specificRootSurfaceArea(double specificRootLength, double rootTissueDensity);
double fineRootAreaIndex(NumericVector Ksoil, NumericVector krhizo, double lai,
                         double specificRootLength, double rootTissueDensity,  
                         double rootLengthDensity );
double fineRootBiomassPerIndividual(NumericVector Ksoil, NumericVector krhizo,  double lai, double N,
                                    double specificRootLength, double rootTissueDensity,  
                                    double rootLengthDensity );
NumericVector rhizosphereMaximumConductance(NumericVector Ksoil, NumericVector fineRootBiomass, double lai, double N,
                                         double specificRootLength, double rootTissueDensity,  
                                         double rootLengthDensity);

double fineRootSoilVolume(double fineRootBiomass, double specificRootLength, double rootLengthDensity );

double coarseRootSoilVolumeFromConductance(double Kmax_rootxylem, double VCroot_kmax, double Al2As,
                                    NumericVector v, NumericVector d, NumericVector rfc);
double coarseRootSoilVolume(NumericVector v, NumericVector d, double depthWidthRatio = 1.0);

NumericVector coarseRootLengthsFromVolume(double VolInd, NumericVector v, NumericVector d, NumericVector rfc);
NumericVector coarseRootLengths(NumericVector v, NumericVector d, double depthWidthRatio = 1.0);

List nonoverlapHorizontalProportions(NumericMatrix V);
List equaloverlapHorizontalProportions(NumericVector poolProportions, NumericMatrix V);
List horizontalProportions(NumericVector poolProportions, NumericVector VolInd, 
                           NumericVector N, NumericMatrix V, 
                           NumericVector d, NumericVector rfc);