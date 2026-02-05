// [[Rcpp::interfaces(r,cpp)]]
#include <RcppArmadillo.h>
#include "modelInput_c.h"
#include "lowlevel_structures_c.h"
using namespace Rcpp;

Rcpp::NumericVector copyWeather_c(const WeatherInputVector& meteo) {
  NumericVector meteovec = NumericVector::create(
    Named("tday") = meteo.tday, 
    Named("prec") = meteo.prec,
    Named("tmin") = meteo.tmin, 
    Named("tmax") = meteo.tmax,
    Named("rhmin") = meteo.rhmin, 
    Named("rhmax") = meteo.rhmax, 
    Named("rad") = meteo.rad, 
    Named("wind") = meteo.wind, 
    Named("Catm") = meteo.Catm,
    Named("Patm") = meteo.Patm,
    Named("pet") = meteo.pet,
    Named("rint") = meteo.rint);
  return(meteovec);
}

Rcpp::NumericVector copyTopo_c(const Topography& topo) {
  NumericVector topovec = NumericVector::create(topo.elevation, topo.slope, topo.aspect);
  topovec.attr("names") = CharacterVector::create("elevation", "slope", "aspect");
  return(topovec);
}

Rcpp::NumericVector copyStandResult_c(const Stand_RESULT& Sres) {
  NumericVector Stand = NumericVector::create(_["LAI"] = Sres.LAI, 
                                              _["LAIherb"] = Sres.LAIherb, 
                                              _["LAIlive"] = Sres.LAIlive,  
                                              _["LAIexpanded"] = Sres.LAIexpanded,
                                              _["LAIdead"] = Sres.LAIdead,
                                              _["Cm"] = Sres.Cm, 
                                              _["LgroundPAR"] = Sres.LgroundPAR, 
                                              _["LgroundSWR"] = Sres.LgroundSWR);
  return(Stand);  
}

Rcpp::DataFrame copySoilResult_c(const Soil_RESULT& Soilres) {
  DataFrame soilDF = DataFrame::create(_["Psi"] = Rcpp::wrap(Soilres.Psi),
                                       _["HerbTranspiration"] = Rcpp::wrap(Soilres.HerbTranspiration),
                                       _["HydraulicInput"] = Rcpp::wrap(Soilres.HydraulicInput),
                                       _["HydraulicOutput"] = Rcpp::wrap(Soilres.HydraulicOutput), 
                                       _["PlantExtraction"] = Rcpp::wrap(Soilres.PlantExtraction));
  return(soilDF);  
}

Rcpp::NumericVector copyWaterBalanceResult_c(const StandWB_RESULT& SWBres) {
  NumericVector WaterBalance = NumericVector::create(_["PET"] = SWBres.PET, 
                                                     _["Rain"] = SWBres.Rain, 
                                                     _["Snow"] = SWBres.Snow, 
                                                     _["NetRain"] = SWBres.NetRain, 
                                                     _["Snowmelt"] = SWBres.Snowmelt,
                                                     _["Runon"] = SWBres.Runon, 
                                                     _["Infiltration"] = SWBres.Infiltration, 
                                                     _["InfiltrationExcess"] = SWBres.InfiltrationExcess, 
                                                     _["SaturationExcess"] = SWBres.SaturationExcess, 
                                                     _["Runoff"] = SWBres.Runoff, 
                                                     _["DeepDrainage"] = SWBres.DeepDrainage, 
                                                     _["CapillarityRise"] = SWBres.CapillarityRise,
                                                     _["SoilEvaporation"] = SWBres.SoilEvaporation, 
                                                     _["HerbTranspiration"] = SWBres.HerbTranspiration,
                                                     _["PlantExtraction"] = SWBres.PlantExtraction, 
                                                     _["Transpiration"] = SWBres.Transpiration,
                                                     _["HydraulicRedistribution"] = SWBres.HydraulicRedistribution);
  return(WaterBalance);
}

Rcpp::NumericMatrix copyNumericMatrix_c(arma::mat comm, int rows, int cols) {
  NumericMatrix out(rows, cols);
  for(int r=0;r<rows;r++) {
    for(int c=0;c<cols; c++) {
      out(r, c) = comm(r, c);
    }
  }
  return(out);
}
