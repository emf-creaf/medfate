#include <RcppArmadillo.h>
#include "communication_structures_c.h"
#include "transpiration_advanced_c.h"
using namespace Rcpp;


Rcpp::DataFrame copyPlantAdvancedTranspirationResult_c(const PlantsAdvancedTranspiration_RESULT& plants, ModelInput& x) {
  DataFrame plantsDF = DataFrame::create(
    _["LAI"] = Rcpp::wrap(plants.LAI),
    _["LAIlive"] = Rcpp::wrap(plants.LAIlive),
    _["FPAR"] = Rcpp::wrap(plants.FPAR),
    _["Extraction"] = Rcpp::wrap(plants.Extraction),
    _["Transpiration"] = Rcpp::wrap(plants.Transpiration),
    _["GrossPhotosynthesis"] = Rcpp::wrap(plants.GrossPhotosynthesis),
    _["NetPhotosynthesis"] = Rcpp::wrap(plants.NetPhotosynthesis),
    _["RootPsi"] = Rcpp::wrap(plants.RootPsi),
    _["StemPsi"] = Rcpp::wrap(plants.StemPsi),
    _["StemPLC"] = Rcpp::wrap(plants.StemPLC),
    _["LeafPLC"] = Rcpp::wrap(plants.LeafPLC),
    _["LeafPsiMin"] = Rcpp::wrap(plants.LeafPsiMin),
    _["LeafPsiMax"] = Rcpp::wrap(plants.LeafPsiMax),
    _["dEdP"] = Rcpp::wrap(plants.dEdP),
    _["DDS"] = Rcpp::wrap(plants.DDS),
    _["StemRWC"] = Rcpp::wrap(plants.StemRWC),
    _["LeafRWC"] = Rcpp::wrap(plants.LeafRWC),
    _["LFMC"] = Rcpp::wrap(plants.LFMC),
    _["WaterBalance"] = Rcpp::wrap(plants.WaterBalance)
  );
  plantsDF.attr("row.names") = x.cohorts.CohortCode;
  return(plantsDF);
}

Rcpp::DataFrame copyLeafAdvancedTranspirationResult_c(const LeafAdvancedTranspiration_RESULT& leaf, ModelInput& x) {
  DataFrame leafDF = DataFrame::create(
    _["LeafPsiMin"] = Rcpp::wrap(leaf.LeafPsiMin),
    _["LeafPsiMax"] = Rcpp::wrap(leaf.LeafPsiMax),
    _["GSWMin"] = Rcpp::wrap(leaf.GSWMin),
    _["GSWMax"] = Rcpp::wrap(leaf.GSWMax),
    _["TempMin"] = Rcpp::wrap(leaf.TempMin),
    _["TempMax"] = Rcpp::wrap(leaf.TempMax)
  );
  return(leafDF);
}

Rcpp::List copyEnergyBalanceResult_c(const EnergyBalance_RESULT& EBres, ModelInput& x) {
  int ncanlayers = x.canopy.zlow.size(); //Number of canopy layers
  int ntimesteps = EBres.Ebalcan.size();
  int nlayers = x.soil.getNlayers();
  DataFrame Tinst = DataFrame::create(
    _["SolarHour"] = Rcpp::wrap(EBres.SolarHour),
    _["Tatm"] = Rcpp::wrap(EBres.Tatm),
    _["Tcan"] = Rcpp::wrap(EBres.Tcan)
  );
  DataFrame CEBinst = DataFrame::create(
    _["SolarHour"] = Rcpp::wrap(EBres.SolarHour),
    _["SWRcan"] = Rcpp::wrap(EBres.SWRcan),
    _["LWRcan"] = Rcpp::wrap(EBres.LWRcan),
    _["LEVcan"] = Rcpp::wrap(EBres.LEVcan),
    _["LEFsnow"] = Rcpp::wrap(EBres.LEFsnow),
    _["Hcan"] = Rcpp::wrap(EBres.Hcan),
    _["Ebalcan"] = Rcpp::wrap(EBres.Ebalcan)
  );
  DataFrame SEBinst = DataFrame::create(
    _["SolarHour"] = Rcpp::wrap(EBres.SolarHour),
    _["Hcansoil"] = Rcpp::wrap(EBres.Hcansoil),
    _["LEVsoil"] = Rcpp::wrap(EBres.LEVsoil),
    _["SWRsoil"] = Rcpp::wrap(EBres.SWRsoil),
    _["LEFsnow"] = Rcpp::wrap(EBres.LEFsnow),
    _["LWRsoil"] = Rcpp::wrap(EBres.LWRsoil),
    _["Ebalsoil"] = Rcpp::wrap(EBres.Ebalsoil)
  );
  NumericMatrix Tcan_mat= copyNumericMatrix_c(EBres.TemperatureLayers, ntimesteps, ncanlayers);
  Tcan_mat.attr("dimnames") = List::create(seq(1,ntimesteps), seq(1,ncanlayers));
  NumericMatrix VPcan_mat= copyNumericMatrix_c(EBres.VaporPressureLayers, ntimesteps, ncanlayers);
  VPcan_mat.attr("dimnames") = List::create(seq(1,ntimesteps), seq(1,ncanlayers));
  NumericMatrix Tsoil_mat= copyNumericMatrix_c(EBres.SoilTemperature, ntimesteps, nlayers);
  Tsoil_mat.attr("dimnames") = List::create(seq(1,ntimesteps), seq(1,nlayers));
  List EnergyBalance = List::create(_["Temperature"]=Tinst, 
                                    _["SoilTemperature"] = Tsoil_mat,
                                    _["CanopyEnergyBalance"] = CEBinst, 
                                    _["SoilEnergyBalance"] = SEBinst,
                                    _["TemperatureLayers"] = Tcan_mat, 
                                    _["VaporPressureLayers"] = VPcan_mat);
  return(EnergyBalance);
}

Rcpp::List copyAdvancedTranspirationResult_c(const AdvancedTranspiration_RESULT& BTres, ModelInput& x) {
  const std::string& rhizosphereOverlap = x.control.rhizosphereOverlap;
  bool plantWaterPools = (rhizosphereOverlap!="total");
  int nlayers = x.soil.getNlayers();
  int numCohorts = x.cohorts.CohortCode.size();
  
  const arma::mat& extractionComm = BTres.extraction;
  Rcpp::NumericMatrix Extraction = copyNumericMatrix_c(extractionComm, numCohorts, nlayers); // this is final extraction of each cohort from each layer
  Extraction.attr("dimnames") = Rcpp::List::create(x.cohorts.CohortCode, Rcpp::seq(1,nlayers));
  
  Rcpp::List ExtractionPools(numCohorts);
  const std::vector< arma::mat>& ExtractionPoolsComm = BTres.extractionPools;
  if(plantWaterPools) {
    for(int c=0;c<numCohorts;c++) {
      const arma::mat& extractionPoolsCohComm = ExtractionPoolsComm[c];
      Rcpp::NumericMatrix ExtractionPoolsCohComm_c = copyNumericMatrix_c(extractionPoolsCohComm, numCohorts, nlayers); // this is final extraction of each cohort from each layer
      ExtractionPoolsCohComm_c.attr("dimnames") = Rcpp::List::create(x.cohorts.CohortCode, Rcpp::seq(1,nlayers));
      ExtractionPools[c] = ExtractionPoolsCohComm_c;
    }
    ExtractionPools.attr("names") = x.cohorts.CohortCode;
  }
  
  NumericVector standVEC = copyStandBasicTranspirationResult_c(BTres.stand);
  
  List EnergyBalance = copyEnergyBalanceResult_c(BTres.energy, x);
  
  List l = List::create(_["cohorts"] = copyCohorts_c(x.cohorts),
                        _["EnergyBalance"] = EnergyBalance,
                        _["Stand"] = standVEC,
                        _["Plants"] = copyPlantAdvancedTranspirationResult_c(BTres.plants, x),
                        _["SunlitLeaves"] = copyLeafAdvancedTranspirationResult_c(BTres.sunlit, x),
                        _["ShadeLeaves"] = copyLeafAdvancedTranspirationResult_c(BTres.shade, x),
                        _["Extraction"] = Extraction,
                        _["ExtractionPools"] = ExtractionPools);
  return(l);
}


void transpirationAdvanced_c(AdvancedTranspiration_RESULT& BTres, AdvancedTranspiration_COMM& BT_comm, ModelInput& x, 
                             const WeatherInputVector& meteovec,  const double elevation) {
  
}