#include <Rcpp.h>

#ifndef TISSUEMOISTURE_H
#define TISSUEMOISTURE_H
#endif
using namespace Rcpp;

double sapwoodWaterCapacity(double Al2As, double height, NumericVector V, NumericVector L, double wd);
double leafWaterCapacity(double SLA, double ld);

double symplasticRelativeWaterContent(double psiSym, double pi0, double epsilon);
double symplasticWaterPotential(double RWC, double pi0, double epsilon);
double apoplasticRelativeWaterContent(double psiApo, double c, double d);
double apoplasticWaterPotential(double RWC, double c, double d);
double segmentRelativeWaterContent(double psiSym, double pi0, double epsilon, 
                                  double psiApo, double c, double d, 
                                  double af);
double tissueRelativeWaterContent(double psiSym, double pi0, double epsilon, 
                                  double psiApo, double c, double d, 
                                  double af, double femb = 0.0);

double turgorLossPoint(double pi0, double epsilon);

List cohortFMC(List spwb, DataFrame SpParams);
