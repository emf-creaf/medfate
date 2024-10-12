#include <Rcpp.h>

#ifndef COMMUNICATION_STRUCTURES_H
#define COMMUNICATION_STRUCTURES_H
#endif
using namespace Rcpp;

void addCommunicationStructures(List x);
void clearCommunicationStructures(List x);

List internalLongWaveRadiation(int ncanlayers);
DataFrame internalCarbonCompartments(DataFrame above);