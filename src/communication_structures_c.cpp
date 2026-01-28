// [[Rcpp::interfaces(r,cpp)]]
#include <RcppArmadillo.h>
#include "modelInput_c.h"
#include "communication_structures_c.h"
using namespace Rcpp;

Rcpp::NumericMatrix copyNumericMatrix_c(arma::mat comm, int rows, int cols) {
  NumericMatrix out(rows, cols);
  for(int r=0;r<rows;r++) {
    for(int c=0;c<cols; c++) {
      out(r, c) = comm(r, c);
    }
  }
  return(out);
}

Rcpp::List copyBasicTranspirationOutput_c(const BasicTranspirationOutput& btc, ModelInput& x) {
  const std::string& rhizosphereOverlap = x.control.rhizosphereOverlap;
  bool plantWaterPools = (rhizosphereOverlap!="total");
  int nlayers = x.soil.getNlayers();
  int numCohorts = x.cohorts.CohortCode.size();
  
  const arma::mat& extractionComm = btc.extraction;
  NumericMatrix Extraction = copyNumericMatrix_c(extractionComm, numCohorts, nlayers); // this is final extraction of each cohort from each layer
  Extraction.attr("dimnames") = List::create(x.cohorts.CohortCode, seq(1,nlayers));
  
  List ExtractionPools(numCohorts);
  const std::vector< arma::mat>& ExtractionPoolsComm = btc.extractionPools;
  if(plantWaterPools) {
    for(int c=0;c<numCohorts;c++) {
      const arma::mat& extractionPoolsCohComm = ExtractionPoolsComm[c];
      NumericMatrix ExtractionPoolsCohComm_c = copyNumericMatrix_c(extractionPoolsCohComm, numCohorts, nlayers); // this is final extraction of each cohort from each layer
      ExtractionPoolsCohComm_c.attr("dimnames") = List::create(x.cohorts.CohortCode, seq(1,nlayers));
      ExtractionPools[c] = ExtractionPoolsCohComm_c;
    }
    ExtractionPools.attr("names") = x.cohorts.CohortCode;
  }
  
  NumericVector standVEC = NumericVector::create(_["LAI"] = btc.stand.LAI,
                                              _["LAIlive"] = btc.stand.LAIlive, 
                                              _["LAIexpanded"] = btc.stand.LAIexpanded, 
                                              _["LAIdead"] = btc.stand.LAIdead);
  
  DataFrame cohortsDF = DataFrame::create(
    _["SP"] = Rcpp::wrap(x.cohorts.SpeciesIndex),
    _["Name"]= Rcpp::wrap(x.cohorts.SpeciesName)
  );
  cohortsDF.attr("row.names") = x.cohorts.CohortCode;
  
  DataFrame plantsDF = DataFrame::create(
    _["LAI"] = Rcpp::wrap(btc.plants.LAI),
    _["LAIlive"] = Rcpp::wrap(btc.plants.LAIlive),
    _["FPAR"] = Rcpp::wrap(btc.plants.FPAR),
    _["AbsorbedSWRFraction"] = Rcpp::wrap(btc.plants.AbsorbedSWRFraction),
    _["Extraction"] = Rcpp::wrap(btc.plants.Extraction),
    _["Transpiration"] = Rcpp::wrap(btc.plants.Transpiration),
    _["GrossPhotosynthesis"] = Rcpp::wrap(btc.plants.GrossPhotosynthesis),
    _["PlantPsi"] = Rcpp::wrap(btc.plants.PlantPsi),
    _["DDS"] = Rcpp::wrap(btc.plants.DDS),
    _["StemRWC"] = Rcpp::wrap(btc.plants.StemRWC),
    _["LeafRWC"] = Rcpp::wrap(btc.plants.LeafRWC),
    _["LFMC"] = Rcpp::wrap(btc.plants.LFMC),
    _["StemPLC"] = Rcpp::wrap(btc.plants.StemPLC),
    _["LeafPLC"] = Rcpp::wrap(btc.plants.LeafPLC),
    _["WaterBalance"] = Rcpp::wrap(btc.plants.WaterBalance)
  );
  plantsDF.attr("row.names") = x.cohorts.CohortCode;
  
  List l = List::create(_["cohorts"] = cohortsDF,
                        _["Stand"] = standVEC,
                        _["Plants"] = plantsDF,
                        _["Extraction"] = Extraction,
                        _["ExtractionPools"] = ExtractionPools);
  return(l);
}
