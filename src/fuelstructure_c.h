#include <RcppArmadillo.h>

#ifndef FUELSTRUCTURE_C_H
#define FUELSTRUCTURE_C_H


const double woodyBulkDensity = 26.43; //kg/m3
const double shortLinearBulkDensity = 26.43; //kg/m3
const double longLinearBulkDensity = 26.43; //kg/m3
const double scaleBulkDensity = 26.43; //kg/m3
const double broadleavedBulkDensity = 13.30; //kg/m3

const double defaultParticleDensity = 400.0; //kg/m3
const double defaultLowHeatContent = 18608.0; //kJ/kg

const double AET = 700; //mm (leads to k = 0.6837 for Lignin = 20)
const double smallBranchDecompositionRate = 0.3336; //year^-1 (AET = 700, Lignin = 35)

const double herbSAV = 11483.0; //m2/m3
const double woodySAV = 1601.05; //m2/m3
const double shortLinearSAV = 6562.0; //m2/m3
const double longLinearSAV = 4921.0; //m2/m3
const double scaleSAV = 2000.0; //m2/m3
const double broadleavedSAV = 8202.0; //m2/m3

const double shortLinearWmax = 0.3248; //kg/m2
const double longLinearWmax = 0.6496; //kg/m2
const double scaleWmax = 0.3248; //kg/m2
const double broadleavedWmax = 0.3472; //kg/m2

const double shortLinearReactionEfficiency = 0.18; //unitless
const double longLinearReactionEfficiency = 0.27; //unitless
const double scaleReactionEfficiency = 0.18; //unitless
const double broadleavedReactionEfficiency = 0.11; //unitless

const double shortLinearRPR = 8.03; //unitless
const double longLinearRPR = 6.35; //unitless
const double scaleRPR = 8.03; //unitless
const double broadleavedRPR = 10.31; //unitless

double crownProportionInLayer_c(double zLow, double zHigh, double H, double Hbc);
double crownFuelInLayer_c(double zLow, double zHigh, double fb, double H, double Hbc);
double crownLengthInLayer_c(double zLow, double zHigh, double cl, double H, double Hbc);
double layerFuelAverageParameter_c(double minHeight, double maxHeight, 
                                   const std::vector<double>& cohortParameter, 
                                   const std::vector<double>& cohortLoading, 
                                   const std::vector<double>& H, 
                                   const std::vector<double>& CR);

#endif
