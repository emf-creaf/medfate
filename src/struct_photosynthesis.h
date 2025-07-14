#include <Rcpp.h>

#ifndef STRUCT_PHOTOSYNTHESIS_H
#define STRUCT_PHOTOSYNTHESIS_H
#endif
using namespace Rcpp;

struct BaldocchiPhoto{
  double Gsw, Cs, Ci, An, Ag;
};

