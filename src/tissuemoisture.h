#include <Rcpp.h>

#ifndef TISSUEMOISTURE_H
#define TISSUEMOISTURE_H
#endif
using namespace Rcpp;

double symplasticRelativeWaterContent(double psi, double pi0, double epsilon);
double symplasticWaterPotential(double RWC, double pi0, double epsilon);
double apoplasticRelativeWaterContent(double psi, double c, double d, double cellWallFraction = 0.07);
double leafRelativeWaterContent(double psi, double pi0, double epsilon, double rwc_res);
double branchRelativeWaterContent(double psi, double wd, double c, double d, double af = 0.80);
double stemRelativeWaterContent(double psi, double pi0, double epsilon, double c, double d, double af);
double fineFuelRelativeWaterContent(double psi, double leaf_pi0, double leaf_eps, double leaf_af, 
                                    double wd, double c, double d, double r635);