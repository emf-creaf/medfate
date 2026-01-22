#ifndef TISSUEMOISTURE_C_H
#define TISSUEMOISTURE_C_H

double leafWaterCapacity_c(double SLA, double ld);

double symplasticRelativeWaterContent_c(double psiSym, double pi0, double epsilon);
double symplasticWaterPotential_c(double RWC, double pi0, double epsilon);
double apoplasticRelativeWaterContent_c(double psiApo, double c, double d);
double apoplasticWaterPotential_c(double RWC, double c, double d);
double tissueRelativeWaterContent_c(double psiSym, double pi0, double epsilon, 
                                  double psiApo, double c, double d, 
                                  double af);

double turgorLossPoint_c(double pi0, double epsilon);

#endif
