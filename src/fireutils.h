#include <Rcpp.h>

#ifndef FIREUTILS_H
#define FIREUTILS_H
#endif
using namespace Rcpp;


NumericVector vectorAddition(NumericVector v1, NumericVector v2);
NumericVector ellipseROS(NumericVector phi, double theta, double vws, double ros);
NumericVector doubleEllipseROS(NumericVector phi, double theta, double vws, double ros);
int getEllipseIndex(int angle);
List rothermel(String modeltype, NumericVector wSI, NumericVector sSI, double delta, double mx_dead,
                  NumericVector hSI, NumericVector mSI, double u, double windDir, double slope, double slopeDir);

double findFireBrandLoftedHeight(double t0, double z0, double zF,double Dp);
double fireBrandFallingHeight(double initialHeight, double timeFalling, double Dp);
double fireBrandFlameHeightFromCanopyStructure(double crownLength, double LAIc);
double fireBrandBurningTimeFromCanopyStructure(double LAIc);
bool willBurnWhenHitFloor(double zIni, double Dp);

double criticalFirelineIntensity(double CBH, double M);
