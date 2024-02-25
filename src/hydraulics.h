#include <Rcpp.h>

#ifndef HYDRAULICS_H
#define HYDRAULICS_H
#endif
using namespace Rcpp;

double K2Psi(double K, double psi_extract, double exp_extract= 3.0);
NumericVector K2Psi(NumericVector K, NumericVector psi_extract, double exp_extract= 3.0);
double Psi2K(double psi, double psi_extract, double exp_extract= 3.0);
NumericVector Psi2K(double psi, NumericVector psi_extract, double exp_extract= 3.0);

double gmin(double leafTemperature, double gmin_20, 
            double TPhase, double Q10_1, double Q10_2);

double averagePsi(NumericVector psi, NumericVector v, double exp_extract, double psi_extract);
double averagePsiPool(NumericMatrix Psi, NumericMatrix RHOPcoh, double exp_extract, double psi_extract);

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

double proportionDefoliationSigmoid(double psiLeaf, double P50, double slope, 
                                    double PLC_crit = 0.88, double P50_cv = 10.0);

double proportionDefoliationWeibull(double psiLeaf, double c, double d, 
                                    double PLC_crit = 0.88, double P50_cv = 10.0);

double EXylem(double psiPlant, double psiUpstream, 
              double kxylemmax, double c, double d, 
              bool allowNegativeFlux = true, double psiCav = 0.0);

double E2psiXylem(double E, double psiUpstream, double kxylemmax, double c, double d, double psiCav = 0.0);
double E2psiXylemUp(double E, double psiDownstream, double kxylemmax, double c, double d, double psiCav = 0.0);



List E2psiBelowground(double E, List hydraulicNetwork,
                      NumericVector psiIni = NumericVector::create(0));

List E2psiAboveground(double E, double psiRootCrown, List hydraulicNetwork);


List E2psiNetwork(double E, List hydraulicNetwork,
                  NumericVector psiIni = NumericVector::create(0));


List supplyFunctionAboveground(NumericVector Erootcrown, NumericVector psiRootCrown, 
                               List hydraulicNetwork);


List supplyFunctionBelowground(List hydraulicNetwork,
                              double minFlow = 0.0,
                              double pCrit = 0.001);

List supplyFunctionNetwork(List hydraulicNetwork,
                           double minFlow = 0.0,
                           double pCrit = 0.001);
