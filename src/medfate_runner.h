#include "Rcpp.h"
#include "medfate.h"
#include "soil_c.h"
class MedfateRunner {
private:
  Soil soil;
public:
  MedfateRunner(Rcpp::DataFrame soilParams, Rcpp::String soilFunctions);
  double run_soil_evaporation(double snowpack, double pet, double LgroundSWR, bool modifySoil);
};