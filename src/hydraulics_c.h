#ifndef HYDRAULICS_C_H
#define HYDRAULICS_C_H

double const maxPsi = -0.000001; //Maximum Psi value to avoid numerical problems
double const cmhead2MPa = 0.00009804139; //Constant to transform cm head to MPa

double K2Psi(double K, double psi_extract, double exp_extract);
double Psi2K(double psi, double psi_extract, double exp_extract);

double gmin(double leafTemperature, double gmin_20, 
            double TPhase, double Q10_1, double Q10_2);

double xylemConductance(double psi, double kxylemmax, double c, double d);
double xylemPsi(double kxylem, double kxylemmax, double c, double d);
double psiCrit(double c, double d, double pCrit);
double vanGenuchtenConductance(double psi, double krhizomax, double n, double alpha);
double xylemConductanceSigmoid(double psi, double kxylemmax, double P50, double slope);

#endif