#include "RcppArmadillo.h"

#ifndef HYDRAULICS_C_H
#define HYDRAULICS_C_H

double const maxPsi = -0.000001; //Maximum Psi value to avoid numerical problems
double const cmhead2MPa = 0.00009804139; //Constant to transform cm head to MPa

double K2Psi_c(double K, double psi_extract, double exp_extract);
double Psi2K_c(double psi, double psi_extract, double exp_extract);

double gmin_c(double leafTemperature, double gmin_20, 
              double TPhase, double Q10_1, double Q10_2);

double xylemConductance_c(double psi, double kxylemmax, double c, double d);
double xylemPsi_c(double kxylem, double kxylemmax, double c, double d);
double psiCrit_c(double c, double d, double pCrit);
double vanGenuchtenConductance_c(double psi, double krhizomax, double n, double alpha);
double xylemConductanceSigmoid_c(double psi, double kxylemmax, double P50, double slope);

double correctConductanceForViscosity_c(double kxylem, double temp);

double averagePsi_c(const std::vector<double>& psi, const std::vector<double>& v, double exp_extract, double psi_extract);
double averagePsiPool_c(const arma::mat& Psi, const arma::mat& RHOPcohV, double exp_extract, double psi_extract);
#endif