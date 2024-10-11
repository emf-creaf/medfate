#include <Rcpp.h>
using namespace Rcpp;



List internalLongWaveRadiation(int ncanlayers) {
  NumericVector Lup(ncanlayers, NA_REAL), Ldown(ncanlayers, NA_REAL), Lnet(ncanlayers, NA_REAL);
  NumericVector tau(ncanlayers, NA_REAL), sumTauComp(ncanlayers, NA_REAL);
  DataFrame LWR_layer = DataFrame::create(_["tau"] = tau,
                                          _["sumTauComp"] = sumTauComp,
                                          _["Ldown"] = Ldown, 
                                          _["Lup"] = Lup,
                                          _["Lnet"] = Lnet);
  List lwr_struct = List::create(_["LWR_layer"] = LWR_layer,
                                 _["Ldown_ground"] = NA_REAL,
                                 _["Lup_ground"] = NA_REAL,
                                 _["Lnet_ground"] = NA_REAL,
                                 _["Ldown_canopy"] = NA_REAL,
                                 _["Lup_canopy"] = NA_REAL,
                                 _["Lnet_canopy"] = NA_REAL,
                                 _["Lnet_cohort_layer"] = NA_REAL);
  return(lwr_struct);
}

// [[Rcpp::export(".addSPWBCommunicationStructures")]]
void addSPWBCommunicationStructures(List x) {
  List control = x["control"];
  String transpirationMode = control["transpirationMode"];
  List ic = as<List>(x["internalCommunication"]);
  if(transpirationMode!="Granier") {
    DataFrame paramsCanopydf = as<DataFrame>(x["canopy"]);
    if(!ic.containsElementNamed("internalLWR")) ic.push_back(internalLongWaveRadiation(paramsCanopydf.nrow()), "internalLWR"); 
  }
  x["internalCommunication"] = ic;
}

// [[Rcpp::export(".clearCommunicationStructures")]]
void clearCommunicationStructures(List x) {
  List control = x["control"];
  x["internalCommunication"] = List::create();
}