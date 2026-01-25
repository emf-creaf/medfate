#include <Rcpp.h>
#include "medfate_runner.h"  
#include <string>
#include "soil_c.h"
#include "hydrology_c.h"
#include "soil.h"
using namespace Rcpp;

MedfateRunner::MedfateRunner(Rcpp::DataFrame soilParams, String soilFunctions) : soil(soilParams, soilFunctions) {
}
double MedfateRunner::run_soil_evaporation(double snowpack, double pet, double LgroundSWR, bool modifySoil) {
  double Esoil = soilEvaporation_c(soil, 
                                   snowpack, pet, LgroundSWR, modifySoil);
  Rcout<< snowpack << " " << pet << " " << LgroundSWR << " " << Esoil<< "\n";
  return(Esoil);
}

RCPP_MODULE(mod_medfate) {
  class_<MedfateRunner>( "MedfateRunner" )
  .constructor<Rcpp::DataFrame, Rcpp::String>()
  .method( "run_soil_evaporation", &MedfateRunner::run_soil_evaporation )
  ;
}