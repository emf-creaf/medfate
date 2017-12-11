#include <Rcpp.h>

#ifndef HYDRAULICS_H
#define HYDRAULICS_H
#endif
using namespace Rcpp;

double K2Psi(double K, double Psi_extract, double ws= 3.0);
NumericVector K2Psi(NumericVector K, NumericVector Psi_extract, double ws= 3.0);
double Psi2K(double psi, double Psi_extract, double ws= 3.0);
NumericVector Psi2K(double psi, NumericVector Psi_extract, double ws= 3.0);
double averagePsi(NumericVector psi, NumericVector v, double c, double d);

double xylemConductance(double psi, double kxylemmax, double c, double d);
double xylemPsi(double kxylem, double kxylemmax, double c, double d);

double taperFactor(double height);
double maximumStemHydraulicConductance(double xylemConductivity, double Al2As, double height, bool taper = false);
double maximumRootHydraulicConductance(double xylemConductivity, double Al2As, NumericVector v, NumericVector d, double depthWidthRatio = 1.0);

NumericVector regulatedPsiTwoElements(double Emax, double psiSoil, double krhizomax, double kxylemmax, double n, double alpha, double c, double d, double dE = 0.1, double psiMax = -10.0);

double findRhizosphereMaximumConductance(double averageResistancePercent, double n, double alpha,
                                         double krootmax, double rootc, double rootd,
                                         double kstemmax, double stemc, double stemd);
List supplyFunctionNetwork(NumericVector psiSoil, 
                           NumericVector krhizomax, NumericVector nsoil, NumericVector alphasoil,
                           NumericVector krootmax, double rootc, double rootd, 
                           double kstemmax, double stemc, double stemd,
                           double psiCav = 0.0,
                           double minFlow = 0.0, int maxNsteps=200, double psiStep = -0.001, double psiMax = -10.0, int ntrial = 10, double psiTol = 0.0001, double ETol = 0.001);
List E2psiNetwork(double E, NumericVector psiSoil, 
                  NumericVector krhizomax, NumericVector nsoil, NumericVector alphasoil,
                  NumericVector krootmax, double rootc, double rootd, 
                  double kstemmax, double stemc, double stemd,
                  NumericVector psiIni = NumericVector::create(0),
                  double psiCav = 0.0,
                  double psiStep = -0.001, double psiMax = -10.0, int ntrial = 10, double psiTol = 0.0001, double ETol = 0.001);