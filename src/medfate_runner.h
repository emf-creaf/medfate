#include "RcppArmadillo.h"
#include "medfate.h"
#include "modelInput_c.h"
#include "soil_c.h"
class MedfateRunner {
private:
  ModelInput input;
public:
  MedfateRunner(Rcpp::List x);
  double run_soil_evaporation(double snowpack, double pet, double LgroundSWR, bool modifySoil);
};