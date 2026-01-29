#include "Rcpp.h"
#include "communication_structures_c.h"
#include "lightextinction_basic_c.h"
#include "modelInput_c.h"

#ifndef TRANSPIRATION_BASIC_C_H
#define TRANSPIRATION_BASIC_C_H


// ----------------------------------------------------------------------------
// Basic Transpiration Output Structure
// ----------------------------------------------------------------------------
struct PlantsBasicTranspiration_RESULT {
  std::vector<double> LAI;
  std::vector<double> LAIlive;
  std::vector<double> FPAR;
  std::vector<double> AbsorbedSWRFraction;
  std::vector<double> Extraction;
  std::vector<double> Transpiration;
  std::vector<double> GrossPhotosynthesis;
  std::vector<double> PlantPsi;
  std::vector<double> DDS;
  std::vector<double> StemRWC;
  std::vector<double> LeafRWC;
  std::vector<double> LFMC;
  std::vector<double> StemPLC;
  std::vector<double> LeafPLC;
  std::vector<double> WaterBalance;
  
  PlantsBasicTranspiration_RESULT(size_t numCohorts = 0) {
    LAI = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    LAIlive = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    FPAR = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    AbsorbedSWRFraction = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    Extraction = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    Transpiration = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    GrossPhotosynthesis = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    PlantPsi = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    DDS = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    StemRWC = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    LeafRWC = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    LFMC = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    StemPLC = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    LeafPLC = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
    WaterBalance = std::vector<double>(numCohorts, medfate::NA_DOUBLE);
  }
};
struct StandBasicTranspiration_RESULT {
  double LAI;
  double LAIlive;
  double LAIexpanded;
  double LAIdead;
};
struct BasicTranspiration_RESULT {
  // Stand-level (4 fields)
  StandBasicTranspiration_RESULT stand;
  
  // Plants data frame
  PlantsBasicTranspiration_RESULT plants;
  
  // Extraction matrices
  arma::mat extraction;
  std::vector<arma::mat> extractionPools;
  
  BasicTranspiration_RESULT(size_t numCohorts = 0, size_t nlayers = 0) : 
    plants(numCohorts), 
    extraction(numCohorts, nlayers, arma::fill::zeros),
    extractionPools(numCohorts) {
    for(size_t c = 0; c < numCohorts; c++) {
      extractionPools[c] = arma::mat(numCohorts, nlayers, arma::fill::zeros);
    }
  }
};

// ----------------------------------------------------------------------------
// Basic Transpiration Communication Structure
// ----------------------------------------------------------------------------
struct BasicTranspiration_COMM {
  AbsorbedSWR_COMM AbSWRcomm;
  std::vector<double> CohASWRF;
  std::vector<double> Tmax;
  std::vector<double> TmaxCoh;
  BasicTranspiration_COMM(size_t numCohorts = 0, size_t ncanlayers = 0) : 
    AbSWRcomm(numCohorts, ncanlayers), 
    CohASWRF(numCohorts),
    Tmax(numCohorts),
    TmaxCoh(numCohorts){}
};



struct ParamsVolume{
  double leafpi0;
  double leafeps;
  double leafaf;
  double stempi0;
  double stemeps;
  double stem_c;
  double stem_d;
  double stemaf;
  double Vsapwood;
  double Vleaf;
  double LAI;
  double LAIlive;
}; 

double plantVol_c(double plantPsi, ParamsVolume pars);
double findNewPlantPsiConnected_c(double flowFromRoots, double plantPsi, double rootCrownPsi,
                                  ParamsVolume parsVol);

void transpirationBasic_c(BasicTranspiration_RESULT& BTres, BasicTranspiration_COMM& BT_comm, ModelInput& x, 
                          const WeatherInputVector& meteovec,  const double elevation);

Rcpp::List copyBasicTranspirationOutput_c(const BasicTranspiration_RESULT& BTres, ModelInput& x);

#endif
