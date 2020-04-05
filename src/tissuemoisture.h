#include <Rcpp.h>

#ifndef TISSUEMOISTURE_H
#define TISSUEMOISTURE_H
#endif
using namespace Rcpp;

double sugarConcentration(double osmoticWP, double temp);
double osmoticWaterPotential(double conc, double temp);
double turgor(double psi, double conc, double temp);
double relativeSapViscosity(double conc, double temp);

double symplasticRelativeWaterContent(double psiSym, double pi0, double epsilon);
double symplasticWaterPotential(double RWC, double pi0, double epsilon);
double apoplasticRelativeWaterContent(double psiApo, double c, double d);
double apoplasticWaterPotential(double RWC, double c, double d);
double segmentRelativeWaterContent(double psiSym, double pi0, double epsilon, 
                                  double psiApo, double c, double d, 
                                  double af);
double tissueFMC(double RWC, double density, double d0 = 1.54);
double tissueRelativeWaterContent(double psiSym, double pi0, double epsilon, 
                                  double psiApo, double c, double d, 
                                  double af, double femb = 0.0);