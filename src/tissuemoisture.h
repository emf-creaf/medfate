#include <Rcpp.h>

#ifndef TISSUEMOISTURE_H
#define TISSUEMOISTURE_H
#endif
using namespace Rcpp;


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

double turgorLossPoint(double pi0, double epsilon);

List cohortFMC(List spwb);
