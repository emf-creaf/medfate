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
double averagePsiPool(NumericMatrix Psi, NumericMatrix RHOPcoh, double c, double d);

double correctConductanceForViscosity(double kxylem, double temp);

double vanGenuchtenConductance(double psi, double krhizomax, double n, double alpha);
double xylemConductance(double psi, double kxylemmax, double c, double d);
double xylemPsi(double kxylem, double kxylemmax, double c, double d);

double taperFactorSavage(double height);
double terminalConduitRadius(double height);
double referenceConductivityHeightFactor(double refheight, double height);
double maximumStemHydraulicConductance(double xylemConductivity, double refheight, double Al2As,  double height, bool taper = false);
NumericVector rootxylemConductanceProportions(NumericVector V, NumericVector L);

NumericVector psi2Weibull(double psi50, double psi88 = NA_REAL, double psi12 = NA_REAL);


NumericVector regulatedPsiTwoElements(double Emax, double psiSoil, double krhizomax, double kxylemmax, double n, double alpha, double c, double d, double dE = 0.1, double psiMax = -10.0);

double maximumSoilPlantConductance(NumericVector krhizomax, NumericVector krootmax, 
                                   double kstemmax, double kleafmax);
double findRhizosphereMaximumConductance(double averageResistancePercent, double n, double alpha,
                                         double krootmax, double rootc, double rootd,
                                         double kstemmax, double stemc, double stemd,
                                         double kleafmax, double leafc, double leafd,
                                         double initialValue = 0.0);


double EXylem(double psiPlant, double psiUpstream, 
              double kxylemmax, double c, double d, 
              bool allowNegativeFlux = true, double psiCav = 0.0);

double E2psiXylem(double E, double psiUpstream, double kxylemmax, double c, double d, double psiCav = 0.0);
double E2psiXylemUp(double E, double psiDownstream, double kxylemmax, double c, double d, double psiCav = 0.0);



List E2psiBelowground(double E, List hydraulicNetwork,
                      NumericVector psiIni = NumericVector::create(0),
                      int ntrial = 10,
                      double psiTol = 0.0001, double ETol = 0.0001);

List E2psiAboveground(double E, double psiRootCrown, List hydraulicNetwork);

List E2psiAbovegroundCapacitance(double E, double psiRootCrown, 
                                 NumericVector psiStemPrev, NumericVector PLCstem,
                                 double psiLeafPrev, 
                                 double kstemmax, double stemc, double stemd,
                                 double kleafmax, double leafc, double leafd,
                                 double Vsapwood, double stemfapo, double stempi0, double stemeps,
                                 double Vleaf, double leaffapo, double leafpi0, double leafeps,
                                 double tstep = 3600.0);

List E2psiAbovegroundCapacitanceDisconnected(double E,                           
                                             NumericVector psiStemPrev, NumericVector PLCstem, NumericVector RWCsympstemPrev, 
                                             double psiLeafPrev, double RWCsympleafPrev,
                                             double kstemmax, double stemc, double stemd,
                                             double kleafmax, double leafc, double leafd,
                                             double Vsapwood, double stemfapo, double stempi0, double stemeps,
                                             double Vleaf, double leaffapo, double leafpi0, double leafeps,
                                             double klat,
                                             double tstep = 3600.0);

List E2psiFineRootLeaf(double E, double psiFineRoot, 
                       double krootmax, double rootc, double rootd,
                       double kstemmax, double stemc, double stemd,
                       double kleafmax, double leafc, double leafd,
                       double PLCstem);

List E2psiNetworkStem1(double E, NumericVector psiSoil, 
                       NumericVector krhizomax, NumericVector nsoil, NumericVector alphasoil,
                       NumericVector krootmax, double rootc, double rootd, 
                       double kstemmax, double stemc, double stemd,
                       double PLCstem,
                       NumericVector psiIni = NumericVector::create(0),
                       int ntrial = 10, 
                       double psiTol = 0.0001, double ETol = 0.0001);

List E2psiNetwork(double E, List hydraulicNetwork,
                  NumericVector psiIni = NumericVector::create(0),
                  int ntrial = 10,
                  double psiTol = 0.0001, double ETol = 0.0001);

List E2psiNetworkCapacitance(double E, NumericVector psiSoil, 
                             NumericVector psiStemPrev, NumericVector PLCstem,
                             double psiLeafPrev, 
                             NumericVector krhizomax, NumericVector nsoil, NumericVector alphasoil,
                             NumericVector krootmax, double rootc, double rootd, 
                             double kstemmax, double stemc, double stemd,
                             double kleafmax, double leafc, double leafd,
                             double Vsapwood, double stemfapo, double stempi0, double stemeps,
                             double Vleaf, double leaffapo, double leafpi0, double leafeps,
                             double tstep = 3600.0,
                             NumericVector psiIni = NumericVector::create(0),
                             int ntrial = 10, 
                             double psiTol = 0.0001, double ETol = 0.0001);



List supplyFunctionAboveground(NumericVector Erootcrown, NumericVector psiRootCrown, 
                               double kstemmax, double stemc, double stemd,
                               double kleafmax, double leafc, double leafd,
                               NumericVector PLCstem);

List supplyFunctionAbovegroundCapacitance(NumericVector Erootcrown, NumericVector psiRootCrown,
                                          NumericVector psiStemPrev, NumericVector PLCstemPrev,
                                          double psiLeafPrev, 
                                          double kstemmax, double stemc, double stemd,
                                          double kleafmax, double leafc, double leafd,
                                          double Vsapwood, double stemfapo, double stempi0, double stemeps,
                                          double Vleaf, double leaffapo, double leafpi0, double leafeps,
                                          double tstep = 3600.0);


List supplyFunctionBelowground(NumericVector psiSoil,
                              NumericVector krhizomax, NumericVector nsoil, NumericVector alphasoil,
                              NumericVector krootmax, double rootc, double rootd,
                              double minFlow = 0.0, int maxNsteps=400,
                              int ntrial = 10, double psiTol = 0.0001, double ETol = 0.0001,
                              double pCrit = 0.001);


List supplyFunctionFineRootLeaf(double psiFineRoot,
                                double krootmax, double rootc, double rootd,
                                double kstemmax, double stemc, double stemd,
                                double kleafmax, double leafc, double leafd,
                                double PLCstem,
                                double minFlow = 0.0, int maxNsteps=400, 
                                double ETol = 0.0001, double pCrit = 0.001);
List supplyFunctionNetworkStem1(NumericVector psiSoil, 
                                NumericVector krhizomax, NumericVector nsoil, NumericVector alphasoil,
                                NumericVector krootmax, double rootc, double rootd, 
                                double kstemmax, double stemc, double stemd,
                                double PLCstem,
                                double minFlow = 0.0, int maxNsteps=400, 
                                int ntrial = 200, double psiTol = 0.0001, double ETol = 0.0001,
                                double pCrit = 0.001);

List supplyFunctionNetwork(NumericVector psiSoil, 
                           NumericVector krhizomax, NumericVector nsoil, NumericVector alphasoil,
                           NumericVector krootmax, double rootc, double rootd, 
                           double kstemmax, double stemc, double stemd,
                           double kleafmax, double leafc, double leafd,
                           NumericVector PLCstem,
                           double minFlow = 0.0, int maxNsteps=400, 
                           int ntrial = 200, double psiTol = 0.0001, double ETol = 0.0001,
                           double pCrit = 0.001);

List supplyFunctionNetworkCapacitance(NumericVector psiSoil, 
                                      NumericVector psiStemPrev, NumericVector PLCstemPrev,
                                      double psiLeafPrev, 
                                      NumericVector krhizomax, NumericVector nsoil, NumericVector alphasoil,
                                      NumericVector krootmax, double rootc, double rootd, 
                                      double kstemmax, double stemc, double stemd,
                                      double kleafmax, double leafc, double leafd,
                                      double Vsapwood, double stemfapo, double stempi0, double stemeps,
                                      double Vleaf, double leaffapo, double leafpi0, double leafeps,
                                      double tstep = 3600.0,
                                      double minFlow = 0.0, int maxNsteps=400, 
                                      int ntrial = 200, double psiTol = 0.0001, double ETol = 0.0001,
                                      double pCrit = 0.001);