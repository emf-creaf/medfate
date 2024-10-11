#include <Rcpp.h>
using namespace Rcpp;

List basicTranspirationOutput(List x) {
  List control = x["control"];
  String rhizosphereOverlap = control["rhizosphereOverlap"];
  bool plantWaterPools = (rhizosphereOverlap!="total");
  DataFrame cohorts = Rcpp::as<Rcpp::DataFrame>(x["cohorts"]);
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  DataFrame soil = Rcpp::as<Rcpp::DataFrame>(x["soil"]);
  int nlayers = Rcpp::as<Rcpp::NumericVector>(soil["widths"]).size();
  int numCohorts = above.nrow();
  
  NumericMatrix Extraction(numCohorts, nlayers); // this is final extraction of each cohort from each layer
  Extraction.attr("dimnames") = List::create(above.attr("row.names"), seq(1,nlayers));
  
  List ExtractionPools(numCohorts);
  if(plantWaterPools) {
    for(int c=0;c<numCohorts;c++) {
      NumericMatrix ExtractionPoolsCoh(numCohorts, nlayers);
      ExtractionPools[c] = ExtractionPoolsCoh;
    }
  }
  
  NumericVector Stand = NumericVector::create(_["LAI"] = NA_REAL,
                                              _["LAIlive"] = NA_REAL, 
                                              _["LAIexpanded"] = NA_REAL, 
                                              _["LAIdead"] = NA_REAL);
  
  NumericVector LAI(numCohorts, NA_REAL);
  NumericVector LAIlive(numCohorts, NA_REAL);
  NumericVector FPAR(numCohorts, NA_REAL);
  NumericVector AbsorbedSWRFraction(numCohorts, NA_REAL);
  NumericVector ExtractionByPlant(numCohorts, NA_REAL);
  NumericVector Transpiration(numCohorts, NA_REAL);
  NumericVector GrossPhotosynthesis(numCohorts, NA_REAL);
  NumericVector PlantPsi(numCohorts, NA_REAL);
  NumericVector DDS(numCohorts, NA_REAL);
  NumericVector StemRWC(numCohorts, NA_REAL);
  NumericVector LeafRWC(numCohorts, NA_REAL);
  NumericVector LFMC(numCohorts, NA_REAL);
  NumericVector StemPLC(numCohorts, NA_REAL);
  NumericVector LeafPLC(numCohorts, NA_REAL);
  NumericVector PWB(numCohorts, NA_REAL);
  
  DataFrame Plants = DataFrame::create(_["LAI"] = LAI,
                                       _["LAIlive"] = LAIlive,
                                       _["FPAR"] = FPAR,
                                       _["AbsorbedSWRFraction"] = AbsorbedSWRFraction, 
                                       _["Extraction"] = ExtractionByPlant,
                                       _["Transpiration"] = Transpiration, 
                                       _["GrossPhotosynthesis"] = GrossPhotosynthesis,
                                       _["PlantPsi"] = PlantPsi, 
                                       _["DDS"] = DDS,
                                       _["StemRWC"] = StemRWC,
                                       _["LeafRWC"] = LeafRWC,
                                       _["LFMC"] = LFMC,
                                       _["StemPLC"] = StemPLC,
                                       _["LeafPLC"] = LeafPLC,
                                       _["WaterBalance"] = PWB);
  Plants.attr("row.names") = above.attr("row.names");
  
  List l = List::create(_["cohorts"] = clone(cohorts),
                        _["Stand"] = Stand,
                        _["Plants"] = Plants,
                        _["Extraction"] = Extraction,
                        _["ExtractionPools"] = ExtractionPools);
  return(l);
}

List basicSPWBOutput(List x) {
  List control = x["control"];
  DataFrame cohorts = Rcpp::as<Rcpp::DataFrame>(x["cohorts"]);
  DataFrame soil = Rcpp::as<Rcpp::DataFrame>(x["soil"]);
  int nlayers = Rcpp::as<Rcpp::NumericVector>(soil["widths"]).size();
  
  
  NumericVector topo = NumericVector::create(NA_REAL, NA_REAL, NA_REAL);
  topo.attr("names") = CharacterVector::create("elevation", "slope", "aspect");
  NumericVector meteovec_bas = NumericVector::create(
    Named("tday") = NA_REAL, 
    Named("prec") = NA_REAL,
    Named("tmin") = NA_REAL, 
    Named("tmax") = NA_REAL,
    Named("rhmin") = NA_REAL, 
    Named("rhmax") = NA_REAL, 
    Named("rad") = NA_REAL, 
    Named("wind") = NA_REAL, 
    Named("Catm") = NA_REAL,
    Named("Patm") = NA_REAL,
    Named("pet") = NA_REAL,
    Named("rint") = NA_REAL);
  NumericVector DB = NumericVector::create(_["PET"] = NA_REAL, 
                                           _["Rain"] = NA_REAL, 
                                           _["Snow"] = NA_REAL, 
                                           _["NetRain"] = NA_REAL, _["Snowmelt"] = NA_REAL,
                                           _["Runon"] = NA_REAL, 
                                           _["Infiltration"] = NA_REAL, _["InfiltrationExcess"] = NA_REAL, _["SaturationExcess"] = NA_REAL, _["Runoff"] = NA_REAL, 
                                           _["DeepDrainage"] = NA_REAL, _["CapillarityRise"] = NA_REAL,
                                           _["SoilEvaporation"] = NA_REAL, _["HerbTranspiration"] = NA_REAL,
                                           _["PlantExtraction"] = NA_REAL, _["Transpiration"] = NA_REAL,
                                           _["HydraulicRedistribution"] = NA_REAL);
  
  
  NumericVector Stand = NumericVector::create(_["LAI"] = NA_REAL, _["LAIherb"] = NA_REAL, 
                                              _["LAIlive"] = NA_REAL,  _["LAIexpanded"] = NA_REAL, _["LAIdead"] = NA_REAL,
                                              _["Cm"] = NA_REAL, _["LgroundPAR"] = NA_REAL, _["LgroundSWR"] = NA_REAL);
  
  NumericVector psiVec(nlayers, 0.0);
  NumericVector EherbVec(nlayers, 0.0);
  NumericVector ExtractionVec(nlayers, 0.0);
  NumericVector soilHydraulicInput(nlayers, 0.0); //Water that entered into the layer across all time steps
  NumericVector soilHydraulicOutput(nlayers, 0.0);  //Water that left the layer across all time steps
  
  DataFrame SB = DataFrame::create(_["Psi"] = psiVec,
                                   _["HerbTranspiration"] = EherbVec,
                                   _["HydraulicInput"] = soilHydraulicInput, 
                                   _["HydraulicOutput"] = soilHydraulicOutput, 
                                   _["PlantExtraction"] = ExtractionVec);
  List l = List::create(_["cohorts"] = clone(cohorts),
                        _["topography"] = topo,
                        _["weather"] = meteovec_bas,
                        _["WaterBalance"] = DB, 
                        _["Soil"] = SB,
                        _["Stand"] = Stand,
                        _["Plants"] = List::create());
  if(control["fireHazardResults"]) l.push_back(List::create(), "FireHazard");
  l.attr("class") = CharacterVector::create("spwb_day","list");
  return(l);
}

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
  if(transpirationMode=="Granier") {
    if(!ic.containsElementNamed("basicTranspirationOutput")) ic.push_back(basicTranspirationOutput(x), "basicTranspirationOutput"); 
    if(!ic.containsElementNamed("basicSPWBOutput")) ic.push_back(basicSPWBOutput(x), "basicSPWBOutput"); 
    List basicTranspirationOutput = ic["basicTranspirationOutput"];
    List basicSPWBOutput = ic["basicSPWBOutput"];
    basicSPWBOutput["Plants"] = basicTranspirationOutput["Plants"];
  } else {
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