#include <Rcpp.h>

#ifndef COMMUNICATION_STRUCTURES_H
#define COMMUNICATION_STRUCTURES_H
#endif
using namespace Rcpp;

List instanceCommunicationStructures(List x);
List generalCommunicationStructures();

List basicTranspirationCommunicationOutput(int numCohorts, int nlayers);
List advancedTranspirationCommunicationOutput(int numCohorts, int nlayers, int ncanlayers, int ntimesteps);
List copyBasicTranspirationOutput(List btc, List x);
List copyAdvancedTranspirationOutput(List atc, List x);
List copyBasicSPWBOutput(List boc, List x);
List copyAdvancedSPWBOutput(List aoc, List x);
List copyBasicGROWTHOutput(List boc, List x);
List copyAdvancedGROWTHOutput(List aoc, List x);