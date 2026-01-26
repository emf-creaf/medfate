#include "modelInput_c.h"
#include "control.h"
#include <Rcpp.h>

/**
 * Implementation of ModelInput class
 *
 * Initialises parameters and internal variables from an R list.
 */
ModelInput::ModelInput(Rcpp::List x) {

  //Control parameters
  Rcpp::List controlList = x["control"];
  control = ControlParameters(controlList);
  
  //Soil
  Rcpp::String soilFunctions = controlList["soilFunctions"];
  soil = Soil(Rcpp::as<Rcpp::DataFrame>(x["soil"]), soilFunctions);

  snowpack = Rcpp::as<double>(x["snowpack"]);
  herbLAI = Rcpp::as<double>(x["herbLAI"]);
  herbLAImax = Rcpp::as<double>(x["herbLAImax"]);

  //Phenology parameters
  Rcpp::DataFrame phenoDF = Rcpp::as<Rcpp::DataFrame>(x["paramsPhenology"]);
  paramsPhenology.phenoType = Rcpp::as< std::vector<std::string> >(phenoDF["PhenologyType"]);
  paramsPhenology.leafDuration = Rcpp::as< std::vector<double> >(phenoDF["LeafDuration"]);
  paramsPhenology.t0gdd = Rcpp::as< std::vector<double> >(phenoDF["t0gdd"]);
  paramsPhenology.Sgdd = Rcpp::as< std::vector<double> >(phenoDF["Sgdd"]);
  paramsPhenology.Tbgdd = Rcpp::as< std::vector<double> >(phenoDF["Tbgdd"]);
  paramsPhenology.Ssen = Rcpp::as< std::vector<double> >(phenoDF["Ssen"]);
  paramsPhenology.Phsen = Rcpp::as< std::vector<double> >(phenoDF["Phsen"]);
  paramsPhenology.Tbsen = Rcpp::as< std::vector<double> >(phenoDF["Tbsen"]);
  paramsPhenology.xsen = Rcpp::as< std::vector<double> >(phenoDF["xsen"]);
  paramsPhenology.ysen = Rcpp::as< std::vector<double> >(phenoDF["ysen"]);

  //Anatomy parameters
  Rcpp::DataFrame anatomyDF = Rcpp::as<Rcpp::DataFrame>(x["paramsAnatomy"]);
  if(anatomyDF.containsElementNamed("Hmax")) paramsAnatomy.Hmax = Rcpp::as< std::vector<double> >(anatomyDF["Hmax"]);
  if(anatomyDF.containsElementNamed("Hmed")) paramsAnatomy.Hmed = Rcpp::as< std::vector<double> >(anatomyDF["Hmed"]);
  paramsAnatomy.Al2As = Rcpp::as< std::vector<double> >(anatomyDF["Al2As"]);
  if(anatomyDF.containsElementNamed("Ar2Al")) paramsAnatomy.Ar2Al = Rcpp::as< std::vector<double> >(anatomyDF["Ar2Al"]);
  paramsAnatomy.SLA = Rcpp::as< std::vector<double> >(anatomyDF["SLA"]);
  if(anatomyDF.containsElementNamed("LeafWidth")) paramsAnatomy.LeafWidth = Rcpp::as< std::vector<double> >(anatomyDF["LeafWidth"]);
  paramsAnatomy.LeafDensity = Rcpp::as< std::vector<double> >(anatomyDF["LeafDensity"]);
  paramsAnatomy.WoodDensity = Rcpp::as< std::vector<double> >(anatomyDF["WoodDensity"]);
  paramsAnatomy.FineRootDensity = Rcpp::as< std::vector<double> >(anatomyDF["FineRootDensity"]);
  if(anatomyDF.containsElementNamed("conduit2sapwood")) paramsAnatomy.conduit2sapwood = Rcpp::as< std::vector<double> >(anatomyDF["conduit2sapwood"]);
  paramsAnatomy.SRL = Rcpp::as< std::vector<double> >(anatomyDF["SRL"]);
  paramsAnatomy.RLD = Rcpp::as< std::vector<double> >(anatomyDF["RLD"]);
  paramsAnatomy.r635 = Rcpp::as< std::vector<double> >(anatomyDF["r635"]);

  // //Water storage parameters
  // Rcpp::DataFrame waterStorageDF = Rcpp::as<Rcpp::DataFrame>(x["paramsWaterStorage"]);
  // paramsWaterStorage.maxFMC = Rcpp::as< std::vector<double> >(waterStorageDF["maxFMC"]);
  // paramsWaterStorage.maxMCleaf = Rcpp::as< std::vector<double> >(waterStorageDF["maxMCleaf"]);
  // paramsWaterStorage.maxMCstem = Rcpp::as< std::vector<double> >(waterStorageDF["maxMCstem"]);
  // paramsWaterStorage.LeafPI0 = Rcpp::as< std::vector<double> >(waterStorageDF["LeafPI0"]);
  // paramsWaterStorage.LeafEPS = Rcpp::as< std::vector<double> >(waterStorageDF["LeafEPS"]);
  // paramsWaterStorage.LeafAF = Rcpp::as< std::vector<double> >(waterStorageDF["LeafAF"]);
  // paramsWaterStorage.Vleaf = Rcpp::as< std::vector<double> >(waterStorageDF["Vleaf"]);
  // paramsWaterStorage.StemPI0 = Rcpp::as< std::vector<double> >(waterStorageDF["StemPI0"]);
  // paramsWaterStorage.StemEPS = Rcpp::as< std::vector<double> >(waterStorageDF["StemEPS"]);
  // paramsWaterStorage.StemAF = Rcpp::as< std::vector<double> >(waterStorageDF["StemAF"]);
  // paramsWaterStorage.Vsapwood = Rcpp::as< std::vector<double> >(waterStorageDF["Vsapwood"]);
  // 
}