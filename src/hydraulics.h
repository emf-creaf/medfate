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

double vanGenuchtenConductance(double psi, double krhizomax, double n, double alpha);
double xylemConductance(double psi, double kxylemmax, double c, double d);
double xylemPsi(double kxylem, double kxylemmax, double c, double d);

double taperFactorSavage(double height);
double terminalConduitRadius(double height);
double referenceConductivityHeightFactor(double refheight, double height);
double maximumStemHydraulicConductance(double xylemConductivity, double refheight, double Al2As,  double height, bool angiosperm = true, bool taper = false);
double maximumRootHydraulicConductance(double xylemConductivity, double Al2As, NumericVector v, NumericVector d, double depthWidthRatio = 1.0);

double stemWaterCapacity(double Al2As, double height, double wd);
double leafWaterCapacity(double SLA, double ld);

NumericVector regulatedPsiTwoElements(double Emax, double psiSoil, double krhizomax, double kxylemmax, double n, double alpha, double c, double d, double dE = 0.1, double psiMax = -10.0);

double findRhizosphereMaximumConductance(double averageResistancePercent, double n, double alpha,
                                         double krootmax, double rootc, double rootd,
                                         double kstemmax, double stemc, double stemd,
                                         double kleafmax, double leafc, double leafd);
List supplyFunctionNetwork(NumericVector psiSoil, 
                           NumericVector krhizomax, NumericVector nsoil, NumericVector alphasoil,
                           NumericVector krootmax, double rootc, double rootd, 
                           double kstemmax, double stemc, double stemd,
                           double kleafmax, double leafc, double leafd,                           
                           double psiCav = 0.0,
                           double minFlow = 0.0, int maxNsteps=200, double psiStep = -0.001, double psiMax = -10.0, int ntrial = 10, double psiTol = 0.0001, double ETol = 0.001);
List supplyFunctionAboveground(NumericVector Erootcrown, NumericVector psiRootCrown, 
                               double EPrev, double psiRootCrownPrev,
                               NumericVector psiStemPrev, NumericVector PLCstemPrev, NumericVector RWCsympstemPrev, 
                               double psiLeafPrev, double RWCsympleafPrev,
                               double kstemmax, double stemc, double stemd,
                               double kleafmax, double leafc, double leafd,
                               double Vsapwood, double stemfapo, double stempi0, double stemeps,
                               double Vleaf, double leaffapo, double leafpi0, double leafeps,
                               double klat, double ksto,
                               double tstep = 3600.0, int nSubSteps = 1000,
                               double psiStep = -0.0001, double psiMax = -10.0);
List supplyFunctionBelowground(NumericVector psiSoil, 
                              NumericVector krhizomax, NumericVector nsoil, NumericVector alphasoil,
                              NumericVector krootmax, double rootc, double rootd, 
                              double minFlow = 0.0, int maxNsteps=400, double psiStep = -0.0001, double psiMax = -10.0, 
                              int ntrial = 10, double psiTol = 0.0001, double ETol = 0.0001);

double E2psiXylem(double E, double psiUpstream, double kxylemmax, double c, double d, double psiCav = 0.0, 
                  double psiStep = -0.01, double psiMax = -10.0);


List E2psiNetwork(double E, NumericVector psiSoil, 
                  NumericVector krhizomax, NumericVector nsoil, NumericVector alphasoil,
                  NumericVector krootmax, double rootc, double rootd, 
                  double kstemmax, double stemc, double stemd,
                  double kleafmax, double leafc, double leafd,
                  NumericVector psiIni = NumericVector::create(0),
                  double psiCav = 0.0,
                  double psiStep = -0.001, double psiMax = -10.0, int ntrial = 10, double psiTol = 0.0001, double ETol = 0.001);

List E2psiAboveGround(double E, double psiRootCrown,      
                      double EPrev, double psiRootCrownPrev,
                      NumericVector psiStemPrev, NumericVector PLCstemPrev, NumericVector RWCsympstemPrev, 
                      double psiLeafPrev, double RWCsympleafPrev,
                      double kstemmax, double stemc, double stemd,
                      double kleafmax, double leafc, double leafd,
                      double Vsapwood, double stemfapo, double stempi0, double stemeps,
                      double Vleaf, double leaffapo, double leafpi0, double leafeps,
                      double klat, double ksto,
                      double tstep = 3600.0, int nSubSteps = 1000,
                      double psiStep = -0.0001, double psiMax = -10.0);

List E2psiAboveGroundDisconnected(double E,                           
                                  NumericVector PLCstem, NumericVector RWCsympstem, 
                                  double psiLeaf, double RWCsympleaf,
                                  double kstemmax, double stemc, double stemd,
                                  double kleafmax, double leafc, double leafd,
                                  double Vsapwood, double stemfapo, double stempi0, double stemeps,
                                  double Vleaf, double leaffapo, double leafpi0, double leafeps,
                                  double klat, double ksto,
                                  double tstep = 3600.0, int nSubSteps = 1000);