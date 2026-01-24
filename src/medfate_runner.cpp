#include <Rcpp.h>
#include "medfate_runner.h"  
#include <string>
#include "soil_c.h"
#include "soil.h"
using namespace Rcpp;

MedfateRunner::MedfateRunner(Rcpp::DataFrame soilParams, String soilFunctions) : soil(soilParams, soilFunctions) {
}
