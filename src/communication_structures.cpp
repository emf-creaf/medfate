// [[Rcpp::interfaces(r,cpp)]]
#include <Rcpp.h>
using namespace Rcpp;

List communicationLongWaveRadiation(int ncanlayers) {
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

DataFrame communicationCanopyTurbulence(int ncanlayers) {
  DataFrame output = DataFrame::create(Named("zmid") = NumericVector(ncanlayers, NA_REAL),
                                       Named("u") = NumericVector(ncanlayers, NA_REAL),
                                       Named("du") = NumericVector(ncanlayers, NA_REAL),
                                       Named("epsilon") = NumericVector(ncanlayers, NA_REAL),
                                       Named("k") = NumericVector(ncanlayers, NA_REAL),
                                       Named("uw") = NumericVector(ncanlayers, NA_REAL));
  return(output);
}

const int SOILWBCOM_dZ_m = 0;
const int SOILWBCOM_dZUp = 1;
const int SOILWBCOM_dZDown = 2;
const int SOILWBCOM_lambda = 3;
const int SOILWBCOM_theta_micro = 4;
const int SOILWBCOM_theta_b = 5;
const int SOILWBCOM_theta_macro = 6;
const int SOILWBCOM_theta_sat_fict = 7;
const int SOILWBCOM_Ksat_b = 8;
const int SOILWBCOM_Ksat_b_ms = 9;
const int SOILWBCOM_Ksat = 10;
const int SOILWBCOM_Ksat_ms = 11;
const int SOILWBCOM_Psi = 12;
const int SOILWBCOM_K = 13;
const int SOILWBCOM_C = 14;
const int SOILWBCOM_Psi_m = 15;
const int SOILWBCOM_K_ms = 16;
const int SOILWBCOM_Kbc = 17;
const int SOILWBCOM_Kbc_ms = 18;
const int SOILWBCOM_C_m = 19;
const int SOILWBCOM_S_macro = 20;
const int SOILWBCOM_e_macro = 21;
const int SOILWBCOM_Kmacro_ms = 22;
const int SOILWBCOM_waterFluidity = 23;
const int SOILWBCOM_a = 24;
const int SOILWBCOM_b = 25;
const int SOILWBCOM_c = 26;
const int SOILWBCOM_d = 27;
const int SOILWBCOM_e = 28;
const int SOILWBCOM_f = 29;
const int SOILWBCOM_K_step_ms05 = 30;
const int SOILWBCOM_C_step_m05 = 31;
const int SOILWBCOM_C_step = 32;
const int SOILWBCOM_C_step_m = 33;
const int SOILWBCOM_K_step_ms = 34;
const int SOILWBCOM_K_step = 35;
const int SOILWBCOM_Psi_step = 36;
const int SOILWBCOM_Psi_step_m = 37;
const int SOILWBCOM_S_macro_step = 38;
const int SOILWBCOM_Kmacro_step_ms = 39;
const int SOILWBCOM_theta_macro_step = 40;
const int SOILWBCOM_theta_micro_step = 41;
const int SOILWBCOM_finalSourceSinks_m3s = 42;
const int SOILWBCOM_capill_below = 43;
const int SOILWBCOM_drain_above = 44;
const int SOILWBCOM_drain_below = 45;
const int SOILWBCOM_lateral_flows_step_mm = 46;

List communicationSoilWaterBalance(int nlayers) {
  int ncol = 50;
  List out(ncol);
  CharacterVector colnames(ncol);
  for(int i = 0; i<ncol; i++) out[i] = NumericVector(nlayers, NA_REAL);
  colnames[SOILWBCOM_dZ_m] = "dZ_m";
  colnames[SOILWBCOM_dZUp] = "dZUp";
  colnames[SOILWBCOM_dZDown] = "dZDown";
  colnames[SOILWBCOM_lambda] = "lambda";
  colnames[SOILWBCOM_theta_micro] = "theta_micro";
  colnames[SOILWBCOM_theta_b] = "theta_b";
  colnames[SOILWBCOM_theta_macro] = "theta_macro";
  colnames[SOILWBCOM_theta_sat_fict] = "theta_sat_fict";
  colnames[SOILWBCOM_Ksat_b] = "Ksat_b";
  colnames[SOILWBCOM_Ksat_b_ms] = "Ksat_b_ms";
  colnames[SOILWBCOM_Ksat] = "Ksat";
  colnames[SOILWBCOM_Ksat_ms] = "Ksat_ms";
  colnames[SOILWBCOM_Psi] = "Psi";
  colnames[SOILWBCOM_K] = "K";
  colnames[SOILWBCOM_C] = "C";
  colnames[SOILWBCOM_Psi_m] = "Psi_m";
  colnames[SOILWBCOM_K_ms] = "K_ms";
  colnames[SOILWBCOM_Kbc] = "Kbc";
  colnames[SOILWBCOM_Kbc_ms] = "Kbc_ms";
  colnames[SOILWBCOM_C_m] = "C_m";
  colnames[SOILWBCOM_S_macro] = "S_macro";
  colnames[SOILWBCOM_e_macro] = "e_macro";
  colnames[SOILWBCOM_Kmacro_ms] = "Kmacro_ms";
  colnames[SOILWBCOM_waterFluidity] = "waterFluidity";
  colnames[SOILWBCOM_a] = "a";
  colnames[SOILWBCOM_b] = "b";
  colnames[SOILWBCOM_c] = "c";
  colnames[SOILWBCOM_d] = "d";
  colnames[SOILWBCOM_e] = "e";
  colnames[SOILWBCOM_f] = "f";
  colnames[SOILWBCOM_K_step_ms05] = "K_step_ms05";
  colnames[SOILWBCOM_C_step_m05] = "C_step_m05";
  colnames[SOILWBCOM_C_step] = "C_step";
  colnames[SOILWBCOM_C_step_m] = "C_step_m";
  colnames[SOILWBCOM_K_step_ms] = "K_step_ms";
  colnames[SOILWBCOM_K_step] = "K_step";
  colnames[SOILWBCOM_Psi_step] = "Psi_step";
  colnames[SOILWBCOM_Psi_step_m] = "Psi_step_m";
  colnames[SOILWBCOM_S_macro_step] = "S_macro_step";
  colnames[SOILWBCOM_Kmacro_step_ms] = "Kmacro_step_ms";
  colnames[SOILWBCOM_theta_macro_step] = "theta_macro_step";
  colnames[SOILWBCOM_theta_micro_step] = "theta_micro_step";
  colnames[SOILWBCOM_finalSourceSinks_m3s] = "finalSourceSinks_m3s";
  colnames[SOILWBCOM_capill_below] = "capill_below";
  colnames[SOILWBCOM_drain_above] = "drain_above";
  colnames[SOILWBCOM_drain_below] = "drain_below";
  colnames[SOILWBCOM_lateral_flows_step_mm] = "lateral_flows_step_mm";
  out.attr("names") = colnames;
  return(out);
}

const int SOILEBCOM_dZ_m = 0;
const int SOILEBCOM_dZUp = 1;
const int SOILEBCOM_dZDown = 2;
const int SOILEBCOM_Zup = 3;
const int SOILEBCOM_Zdown = 4;
const int SOILEBCOM_Zcent = 5;
const int SOILEBCOM_a = 6;
const int SOILEBCOM_b = 7;
const int SOILEBCOM_c = 8;
const int SOILEBCOM_d = 9;
const int SOILEBCOM_e = 10;
const int SOILEBCOM_f = 11;
const int SOILEBCOM_k_up = 12;
const int SOILEBCOM_k_down = 13;
const int SOILEBCOM_tempch = 14;
List communicationSoilEnergyBalance(int nlayers) {
  int ncol = 15;
  List out(ncol);
  CharacterVector colnames(ncol);
  for(int i = 0; i<ncol; i++) out[i] = NumericVector(nlayers, NA_REAL);
  colnames[SOILEBCOM_dZ_m] = "dZ_m";
  colnames[SOILEBCOM_dZUp] = "dZUp";
  colnames[SOILEBCOM_dZDown] = "dZDown";
  colnames[SOILEBCOM_Zup] = "Zup";
  colnames[SOILEBCOM_Zdown] = "Zdown";
  colnames[SOILEBCOM_Zcent] = "Zcent";
  colnames[SOILEBCOM_a] = "a";
  colnames[SOILEBCOM_b] = "b";
  colnames[SOILEBCOM_c] = "c";
  colnames[SOILEBCOM_d] = "d";
  colnames[SOILEBCOM_e] = "e";
  colnames[SOILEBCOM_f] = "f";
  colnames[SOILEBCOM_k_up] = "k_up";
  colnames[SOILEBCOM_k_down] = "k_down";
  colnames[SOILEBCOM_tempch] = "tempch";
  out.attr("names") = colnames;
  return(out);
}

List basicTranspirationCommunicationOutput(int numCohorts, int nlayers) {

  NumericMatrix Extraction(numCohorts, nlayers); // this is final extraction of each cohort from each layer

  List ExtractionPools(numCohorts);
  for(int c=0;c<numCohorts;c++) {
    NumericMatrix ExtractionPoolsCoh(numCohorts, nlayers);
    std::fill(ExtractionPoolsCoh.begin(), ExtractionPoolsCoh.end(), 0.0);
    ExtractionPools[c] = ExtractionPoolsCoh;
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

  List l = List::create(_["Stand"] = Stand,
                        _["Plants"] = Plants,
                        _["Extraction"] = Extraction,
                        _["ExtractionPools"] = ExtractionPools);
  return(l);
}
List basicSPWBCommunicationOutput(List outputTransp, int nlayers) {
  
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
  NumericVector WaterBalance = NumericVector::create(_["PET"] = NA_REAL, 
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
  
  DataFrame Soil = DataFrame::create(_["Psi"] = NumericVector(nlayers, 0.0),
                                     _["HerbTranspiration"] = NumericVector(nlayers, 0.0),
                                     _["HydraulicInput"] = NumericVector(nlayers, 0.0), //Water that entered into the layer across all time steps
                                     _["HydraulicOutput"] = NumericVector(nlayers, 0.0), //Water that left the layer across all time steps
                                     _["PlantExtraction"] = NumericVector(nlayers, 0.0));
  
  List l = List::create(_["topography"] = topo,
                        _["weather"] = meteovec_bas,
                        _["WaterBalance"] = WaterBalance, 
                        _["Soil"] = Soil,
                        _["Stand"] = Stand,
                        _["Plants"] = outputTransp["Plants"]);
  l.push_back(List::create(), "FireHazard");
  l.attr("class") = CharacterVector::create("spwb_day","list");
  return(l);
}
List basicGROWTHCommunicationOutput(List spwbOut, int numCohorts, int nlayers) {
  
  DataFrame labileCarbonBalance = DataFrame::create(_["GrossPhotosynthesis"] = NumericVector(numCohorts, NA_REAL),
                                                    _["MaintenanceRespiration"] = NumericVector(numCohorts, NA_REAL),
                                                    _["GrowthCosts"] = NumericVector(numCohorts, NA_REAL),
                                                    _["RootExudation"] = NumericVector(numCohorts, NA_REAL),
                                                    _["LabileCarbonBalance"] = NumericVector(numCohorts, NA_REAL),
                                                    _["SugarLeaf"] = NumericVector(numCohorts, NA_REAL),
                                                    _["StarchLeaf"] = NumericVector(numCohorts, NA_REAL),
                                                    _["SugarSapwood"] = NumericVector(numCohorts, NA_REAL),
                                                    _["StarchSapwood"] = NumericVector(numCohorts, NA_REAL),
                                                    _["SugarTransport"] = NumericVector(numCohorts, NA_REAL));
  
  NumericVector standCB = {NA_REAL, NA_REAL, NA_REAL, NA_REAL};
  standCB.attr("names") = CharacterVector({"GrossPrimaryProduction", "MaintenanceRespiration", "SynthesisRespiration", "NetPrimaryProduction"});
  
  //Final Biomass compartments
  DataFrame plantStructure = DataFrame::create(
    _["LeafBiomass"] = NumericVector(numCohorts, NA_REAL),
    _["SapwoodBiomass"] = NumericVector(numCohorts, NA_REAL),
    _["FineRootBiomass"] = NumericVector(numCohorts, NA_REAL),
    _["LeafArea"] = NumericVector(numCohorts, NA_REAL),
    _["SapwoodArea"] = NumericVector(numCohorts, NA_REAL),
    _["FineRootArea"] = NumericVector(numCohorts, NA_REAL),
    _["HuberValue"] = NumericVector(numCohorts, NA_REAL),
    _["RootAreaLeafArea"] = NumericVector(numCohorts, NA_REAL),
    _["DBH"] = NumericVector(numCohorts, NA_REAL),
    _["Height"] = NumericVector(numCohorts, NA_REAL)
  );
  
  DataFrame growthMortality = DataFrame::create(
    _["SAgrowth"] = NumericVector(numCohorts, NA_REAL),
    _["LAgrowth"] = NumericVector(numCohorts, NA_REAL),
    _["FRAgrowth"] = NumericVector(numCohorts, NA_REAL),
    _["StarvationRate"] = NumericVector(numCohorts, NA_REAL),
    _["DessicationRate"] = NumericVector(numCohorts, NA_REAL),
    _["MortalityRate"] = NumericVector(numCohorts, NA_REAL)
  );
  
  DataFrame plantBiomassBalance = DataFrame::create(_["InitialDensity"] = NumericVector(numCohorts, 0.0),
                                                    _["InitialSapwoodBiomass"] = NumericVector(numCohorts, 0.0),
                                                    _["InitialStructuralBiomass"] = NumericVector(numCohorts, 0.0),
                                                    _["StructuralBiomassBalance"] = NumericVector(numCohorts, 0.0),
                                                    _["StructuralBiomassChange"] = NumericVector(numCohorts, 0.0),
                                                    _["InitialLabileBiomass"] = NumericVector(numCohorts, 0.0),
                                                    _["LabileBiomassBalance"] = NumericVector(numCohorts, 0.0),
                                                    _["LabileBiomassChange"] = NumericVector(numCohorts, 0.0),
                                                    _["InitialLivingPlantBiomass"] = NumericVector(numCohorts, 0.0),
                                                    _["InitialPlantBiomass"] = NumericVector(numCohorts, 0.0),
                                                    _["PlantBiomassBalance"] = NumericVector(numCohorts, 0.0),
                                                    _["PlantBiomassChange"] = NumericVector(numCohorts, 0.0),
                                                    _["MortalityBiomassLoss"] = NumericVector(numCohorts, 0.0),
                                                    _["InitialCohortBiomass"] = NumericVector(numCohorts, 0.0),
                                                    _["CohortBiomassBalance"] = NumericVector(numCohorts, 0.0),
                                                    _["CohortBiomassChange"] = NumericVector(numCohorts, 0.0));
  
  List l = List::create(_["topography"] = spwbOut["topography"],
                        _["weather"] = spwbOut["weather"],
                        _["WaterBalance"] = spwbOut["WaterBalance"], 
                        _["CarbonBalance"] = standCB,
                        _["Soil"] = spwbOut["Soil"], 
                        _["Stand"] = spwbOut["Stand"], 
                        _["Plants"] = spwbOut["Plants"],
                        _["LabileCarbonBalance"] = labileCarbonBalance,
                        _["PlantBiomassBalance"] = plantBiomassBalance,
                        _["PlantStructure"] = plantStructure,
                        _["GrowthMortality"] = growthMortality);
  l.push_back(spwbOut["FireHazard"], "FireHazard");
  l.attr("class") = CharacterVector::create("growth_day","list");
  return(l);
}

List advancedTranspirationCommunicationOutput(int numCohorts, int nlayers, int ncanlayers, int ntimesteps) {
  NumericMatrix SoilWaterExtract(numCohorts, nlayers);
  std::fill(SoilWaterExtract.begin(), SoilWaterExtract.end(), 0.0);
  
  List ExtractionPools(numCohorts);
  for(int c=0;c<numCohorts;c++) {
    NumericMatrix ExtractionPoolsCoh(numCohorts, nlayers);
    std::fill(ExtractionPoolsCoh.begin(), ExtractionPoolsCoh.end(), 0.0);
    ExtractionPools[c] = ExtractionPoolsCoh;
  }
  
  NumericMatrix soilLayerExtractInst(nlayers, ntimesteps);
  std::fill(soilLayerExtractInst.begin(), soilLayerExtractInst.end(), 0.0);

  NumericMatrix minPsiRhizo(numCohorts, nlayers);
  std::fill(minPsiRhizo.begin(), minPsiRhizo.end(), NA_REAL);

  NumericVector Stand = NumericVector::create(_["LAI"] = NA_REAL,
                                              _["LAIlive"] = NA_REAL, 
                                              _["LAIexpanded"] = NA_REAL, 
                                              _["LAIdead"] = NA_REAL);
  
  // ARRANGE OUTPUT
  NumericVector solarHour(ntimesteps, NA_REAL);
  NumericVector Tatm(ntimesteps, NA_REAL), Tcan(ntimesteps, NA_REAL);
  NumericVector Hcansoil(ntimesteps, NA_REAL);
  NumericVector Ebalsoil(ntimesteps, NA_REAL);
  NumericVector LEVsoil(ntimesteps, NA_REAL);
  NumericVector LEFsnow(ntimesteps, NA_REAL);
  NumericVector abs_SWR_soil(ntimesteps, NA_REAL);
  NumericVector net_LWR_soil(ntimesteps, NA_REAL);
  NumericVector abs_SWR_can(ntimesteps, NA_REAL);
  NumericVector net_LWR_can(ntimesteps, NA_REAL);
  NumericVector LEVcan(ntimesteps, NA_REAL), Hcan_heat(ntimesteps, NA_REAL), Ebal(ntimesteps, NA_REAL);
  
  DataFrame Tinst = DataFrame::create(_["SolarHour"] = solarHour, 
                                      _["Tatm"] = Tatm, _["Tcan"] = Tcan);
  DataFrame CEBinst = DataFrame::create(_["SolarHour"] = solarHour, 
                                        _["SWRcan"] = abs_SWR_can, 
                                        _["LWRcan"] = net_LWR_can,
                                        _["LEVcan"] = LEVcan, 
                                        _["LEFsnow"] = LEFsnow, 
                                        _["Hcan"] = Hcan_heat, 
                                        _["Ebalcan"] = Ebal);
  DataFrame SEBinst = DataFrame::create(_["SolarHour"] = solarHour, 
                                        _["Hcansoil"] = Hcansoil, 
                                        _["LEVsoil"] = LEVsoil, 
                                        _["SWRsoil"] = abs_SWR_soil, 
                                        _["LWRsoil"] = net_LWR_soil,
                                        _["Ebalsoil"] = Ebalsoil);
  NumericMatrix Tcan_mat(ntimesteps, ncanlayers);
  NumericMatrix VPcan_mat(ntimesteps, ncanlayers);
  NumericMatrix Tsoil_mat(ntimesteps, nlayers);
  List EB = List::create(_["Temperature"]=Tinst, 
                         _["SoilTemperature"] = Tsoil_mat,
                         _["CanopyEnergyBalance"] = CEBinst, 
                         _["SoilEnergyBalance"] = SEBinst,
                         _["TemperatureLayers"] = Tcan_mat, 
                         _["VaporPressureLayers"] = VPcan_mat);
  
  NumericVector LAI(numCohorts,0.0);
  NumericVector LAIlive(numCohorts,0.0);
  NumericVector PARcohort(numCohorts,0.0);
  NumericVector SoilExtractCoh(numCohorts,0.0);
  NumericVector DDS(numCohorts, 0.0), LFMC(numCohorts, 0.0);
  NumericVector Eplant(numCohorts, 0.0), Anplant(numCohorts, 0.0), Agplant(numCohorts, 0.0);
  NumericVector minStemPsi(numCohorts, NA_REAL), minRootPsi(numCohorts,NA_REAL); //Minimum potentials experienced
  NumericVector minLeafPsi(numCohorts,NA_REAL), maxLeafPsi(numCohorts,NA_REAL); 
  NumericVector PLClm(numCohorts, NA_REAL), PLCsm(numCohorts, NA_REAL);
  NumericVector dEdPm(numCohorts, NA_REAL), PWB(numCohorts,0.0);
  NumericVector RWCsm(numCohorts, NA_REAL), RWClm(numCohorts, NA_REAL);
  
  DataFrame Plants = DataFrame::create(_["LAI"] = LAI,
                                       _["LAIlive"] = LAIlive,
                                       _["FPAR"] = PARcohort,
                                       _["Extraction"] = SoilExtractCoh,
                                       _["Transpiration"] = Eplant,
                                       _["GrossPhotosynthesis"] = Agplant,
                                       _["NetPhotosynthesis"] = Anplant,
                                       _["RootPsi"] = minRootPsi, 
                                       _["StemPsi"] = minStemPsi, 
                                       _["LeafPLC"] = PLClm, //Average daily leaf PLC
                                       _["StemPLC"] = PLCsm, //Average daily stem PLC
                                       _["LeafPsiMin"] = minLeafPsi, 
                                       _["LeafPsiMax"] = maxLeafPsi, 
                                       _["dEdP"] = dEdPm,//Average daily soilplant conductance
                                       _["DDS"] = DDS, //Daily drought stress is the ratio of average soil plant conductance over its maximum value
                                       _["StemRWC"] = RWCsm,
                                       _["LeafRWC"] = RWClm,
                                       _["LFMC"] = LFMC,
                                       _["WaterBalance"] = PWB);

  NumericVector maxGSW_SL(numCohorts,NA_REAL), maxGSW_SH(numCohorts,NA_REAL); 
  NumericVector minGSW_SL(numCohorts,NA_REAL), minGSW_SH(numCohorts,NA_REAL); 
  NumericVector maxTemp_SL(numCohorts,NA_REAL), maxTemp_SH(numCohorts,NA_REAL); 
  NumericVector minTemp_SL(numCohorts,NA_REAL), minTemp_SH(numCohorts,NA_REAL); 
  NumericVector minLeafPsi_SL(numCohorts,NA_REAL), maxLeafPsi_SL(numCohorts,NA_REAL); 
  NumericVector minLeafPsi_SH(numCohorts,NA_REAL), maxLeafPsi_SH(numCohorts,NA_REAL);
  
  DataFrame Sunlit = DataFrame::create(
    _["LeafPsiMin"] = minLeafPsi_SL, 
    _["LeafPsiMax"] = maxLeafPsi_SL, 
    _["GSWMin"] = minGSW_SL,
    _["GSWMax"] = maxGSW_SL,
    _["TempMin"] = minTemp_SL,
    _["TempMax"] = maxTemp_SL  
  );
  // Rcout<<Sunlit.nrow()<<"\n";
  
  DataFrame Shade = DataFrame::create(
    _["LeafPsiMin"] = minLeafPsi_SH, 
    _["LeafPsiMax"] = maxLeafPsi_SH, 
    _["GSWMin"] = minGSW_SH,
    _["GSWMax"] = maxGSW_SH,
    _["TempMin"] = minTemp_SH,
    _["TempMax"] = maxTemp_SH  
  );

  
  NumericMatrix Einst(numCohorts, ntimesteps);
  NumericMatrix Aninst(numCohorts, ntimesteps), Aginst(numCohorts, ntimesteps);
  NumericMatrix dEdPInst(numCohorts, ntimesteps);
  NumericMatrix LeafPsiInst(numCohorts, ntimesteps), StemPsiInst(numCohorts, ntimesteps);
  NumericMatrix RootPsiInst(numCohorts, ntimesteps);
  NumericMatrix LeafSympPsiInst(numCohorts, ntimesteps), StemSympPsiInst(numCohorts, ntimesteps);
  NumericMatrix StemPLC(numCohorts, ntimesteps), LeafPLC(numCohorts, ntimesteps);
  NumericMatrix LeafRWCInst(numCohorts, ntimesteps), StemRWCInst(numCohorts, ntimesteps);
  NumericMatrix LeafSympRWCInst(numCohorts, ntimesteps), StemSympRWCInst(numCohorts, ntimesteps);
  NumericMatrix PWBinst(numCohorts, ntimesteps);

  List PlantsInst = List::create(
    _["E"]=Einst, _["Ag"]=Aginst, _["An"]=Aninst,
    _["dEdP"] = dEdPInst,
    _["RootPsi"] = RootPsiInst, 
    _["StemPsi"] = StemPsiInst,
    _["LeafPsi"] = LeafPsiInst,
    _["StemSympPsi"] = StemSympPsiInst,
    _["LeafSympPsi"] = LeafSympPsiInst,
    _["StemPLC"] = StemPLC, 
    _["LeafPLC"] = LeafPLC, 
    _["StemRWC"] = StemRWCInst,
    _["LeafRWC"] = LeafRWCInst,
    _["StemSympRWC"] = StemSympRWCInst,
    _["LeafSympRWC"] = LeafSympRWCInst,
    _["PWB"] = PWBinst);
  
  NumericMatrix LAI_SL(numCohorts, ntimesteps);
  NumericMatrix LAI_SH(numCohorts, ntimesteps);
  NumericMatrix Vmax298_SL(numCohorts, ntimesteps);
  NumericMatrix Vmax298_SH(numCohorts, ntimesteps);
  NumericMatrix Jmax298_SL(numCohorts, ntimesteps);
  NumericMatrix Jmax298_SH(numCohorts, ntimesteps);
  NumericMatrix SWR_SL(numCohorts, ntimesteps);
  NumericMatrix SWR_SH(numCohorts, ntimesteps);
  NumericMatrix PAR_SL(numCohorts, ntimesteps);
  NumericMatrix PAR_SH(numCohorts, ntimesteps);
  NumericMatrix LWR_SL(numCohorts, ntimesteps);
  NumericMatrix LWR_SH(numCohorts, ntimesteps);
  NumericMatrix An_SL(numCohorts, ntimesteps), Ag_SL(numCohorts, ntimesteps);
  NumericMatrix An_SH(numCohorts, ntimesteps), Ag_SH(numCohorts, ntimesteps);
  NumericMatrix Ci_SL(numCohorts, ntimesteps);
  NumericMatrix Ci_SH(numCohorts, ntimesteps);
  NumericMatrix E_SL(numCohorts, ntimesteps);
  NumericMatrix E_SH(numCohorts, ntimesteps);
  NumericMatrix GSW_SH(numCohorts, ntimesteps);
  NumericMatrix GSW_SL(numCohorts, ntimesteps);
  NumericMatrix VPD_SH(numCohorts, ntimesteps);
  NumericMatrix VPD_SL(numCohorts, ntimesteps);
  NumericMatrix Temp_SH(numCohorts, ntimesteps);
  NumericMatrix Temp_SL(numCohorts, ntimesteps);
  NumericMatrix Psi_SL(numCohorts, ntimesteps);
  NumericMatrix Psi_SH(numCohorts, ntimesteps);

  List ShadeInst = List::create(
    _["LAI"] = LAI_SH,
    _["Vmax298"] = Vmax298_SH,
    _["Jmax298"] = Jmax298_SH,
    _["Abs_SWR"] = SWR_SH,
    _["Abs_PAR"]=PAR_SH,
    _["Net_LWR"] = LWR_SH,
    _["Ag"] = Ag_SH,
    _["An"] = An_SH,
    _["Ci"] = Ci_SH,
    _["E"] = E_SH,
    _["Gsw"] = GSW_SH,
    _["VPD"] = VPD_SH,
    _["Temp"] = Temp_SH,
    _["Psi"] = Psi_SH);
  List SunlitInst = List::create(
    _["LAI"] = LAI_SL,
    _["Vmax298"] = Vmax298_SL,
    _["Jmax298"] = Jmax298_SL,
    _["Abs_SWR"]=SWR_SL,
    _["Abs_PAR"]=PAR_SL,
    _["Net_LWR"] = LWR_SL,
    _["Ag"] = Ag_SL,
    _["An"] = An_SL,
    _["Ci"] = Ci_SL,
    _["E"] = E_SL,
    _["Gsw"] = GSW_SL,
    _["VPD"] = VPD_SL,
    _["Temp"] = Temp_SL,
    _["Psi"] = Psi_SL);
  
  List lwrExtinctionList(ntimesteps);
  for(int n=0;n<ntimesteps;n++) {
    lwrExtinctionList[n] = communicationLongWaveRadiation(ncanlayers);
  }
  List supply(numCohorts); 
  for(int c=0;c<numCohorts;c++) {
    supply[c] = List::create(); //To be replaced
  }

  List l = List::create(_["EnergyBalance"] = EB,
                        _["Extraction"] = SoilWaterExtract,
                        _["ExtractionPools"] = ExtractionPools,
                        _["RhizoPsi"] = minPsiRhizo,
                        _["Stand"] = Stand,
                        _["Plants"] = Plants,
                        _["SunlitLeaves"] = Sunlit,
                        _["ShadeLeaves"] = Shade,
                        _["ExtractionInst"] = soilLayerExtractInst,
                        _["RadiationInputInst"] = DataFrame::create(), //To be replaced
                        _["PlantsInst"] = PlantsInst,
                        _["SunlitLeavesInst"] = SunlitInst,
                        _["ShadeLeavesInst"] = ShadeInst,
                        _["LightExtinction"] = List::create(), //To be replaced
                        _["LWRExtinction"] = lwrExtinctionList,
                        _["CanopyTurbulence"] = communicationCanopyTurbulence(ncanlayers)); //To be replaced
  
  List outPhotoSunlit(numCohorts);
  List outPhotoShade(numCohorts);
  List outPMSunlit(numCohorts);
  List outPMShade(numCohorts);
  l.push_back(supply, "SupplyFunctions");
  l.push_back(outPhotoSunlit, "PhotoSunlitFunctions");
  l.push_back(outPhotoShade, "PhotoShadeFunctions");
  l.push_back(outPMSunlit, "PMSunlitFunctions");
  l.push_back(outPMShade, "PMShadeFunctions");
  
  l.attr("class") = CharacterVector::create("pwb_day","list");
  return(l);
}
List advancedSPWBCommunicationOutput(List outputTransp, int nlayers) {

  NumericVector topo = NumericVector::create(NA_REAL, NA_REAL, NA_REAL);
  topo.attr("names") = CharacterVector::create("elevation", "slope", "aspect");
  NumericVector meteovec_adv = NumericVector::create(
    Named("tmin") = NA_REAL, 
    Named("tmax") = NA_REAL,
    Named("tminPrev") = NA_REAL, 
    Named("tmaxPrev") = NA_REAL, 
    Named("tminNext") = NA_REAL, 
    Named("prec") = NA_REAL,
    Named("rhmin") = NA_REAL, 
    Named("rhmax") = NA_REAL, 
    Named("rad") = NA_REAL, 
    Named("wind") = NA_REAL, 
    Named("Catm") = NA_REAL,
    Named("Patm") = NA_REAL,
    Named("pet") = NA_REAL,
    Named("rint") = NA_REAL);
  
  NumericVector WaterBalance = NumericVector::create(_["PET"] = NA_REAL, 
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
  
  DataFrame Soil = DataFrame::create(_["Psi"] = NumericVector(nlayers, 0.0),
                                     _["HerbTranspiration"] = NumericVector(nlayers, 0.0),
                                     _["HydraulicInput"] = NumericVector(nlayers, 0.0), //Water that entered into the layer across all time steps
                                     _["HydraulicOutput"] = NumericVector(nlayers, 0.0), //Water that left the layer across all time steps
                                     _["PlantExtraction"] = NumericVector(nlayers, 0.0));
  
  List l = List::create(_["topography"] = topo,
                        _["weather"] = meteovec_adv,
                        _["WaterBalance"] = WaterBalance, 
                        _["EnergyBalance"] = outputTransp["EnergyBalance"],
                        _["Soil"] = Soil, 
                        _["Stand"] = Stand, 
                        _["Plants"] = outputTransp["Plants"],
                        _["RhizoPsi"] = outputTransp["RhizoPsi"],
                        _["SunlitLeaves"] = outputTransp["SunlitLeaves"],
                        _["ShadeLeaves"] = outputTransp["ShadeLeaves"],
                        _["ExtractionInst"] = outputTransp["ExtractionInst"],
                        _["PlantsInst"] = outputTransp["PlantsInst"],
                        _["RadiationInputInst"] = outputTransp["RadiationInputInst"],
                        _["SunlitLeavesInst"] = outputTransp["SunlitLeavesInst"],
                        _["ShadeLeavesInst"] = outputTransp["ShadeLeavesInst"],
                        _["LightExtinction"] = outputTransp["LightExtinction"],
                        _["LWRExtinction"] = outputTransp["LWRExtinction"],
                        _["CanopyTurbulence"] = outputTransp["CanopyTurbulence"]);
  l.push_back(List::create(), "FireHazard");
  l.attr("class") = CharacterVector::create("spwb_day","list");
  return(l);
}
List advancedGROWTHCommunicationOutput(List spwbOut, int numCohorts, int ntimesteps) {

  DataFrame labileCarbonBalance = DataFrame::create(_["GrossPhotosynthesis"] = NumericVector(numCohorts, NA_REAL),
                                                    _["MaintenanceRespiration"] = NumericVector(numCohorts, NA_REAL),
                                                    _["GrowthCosts"] = NumericVector(numCohorts, NA_REAL),
                                                    _["RootExudation"] = NumericVector(numCohorts, NA_REAL),
                                                    _["LabileCarbonBalance"] = NumericVector(numCohorts, NA_REAL),
                                                    _["SugarLeaf"] = NumericVector(numCohorts, NA_REAL),
                                                    _["StarchLeaf"] = NumericVector(numCohorts, NA_REAL),
                                                    _["SugarSapwood"] = NumericVector(numCohorts, NA_REAL),
                                                    _["StarchSapwood"] = NumericVector(numCohorts, NA_REAL),
                                                    _["SugarTransport"] = NumericVector(numCohorts, NA_REAL));

  NumericVector standCB = {NA_REAL, NA_REAL, NA_REAL, NA_REAL};
  standCB.attr("names") = CharacterVector({"GrossPrimaryProduction", "MaintenanceRespiration", "SynthesisRespiration", "NetPrimaryProduction"});
  
  //Final Biomass compartments
  DataFrame plantStructure = DataFrame::create(
    _["LeafBiomass"] = NumericVector(numCohorts, NA_REAL),
    _["SapwoodBiomass"] = NumericVector(numCohorts, NA_REAL),
    _["FineRootBiomass"] = NumericVector(numCohorts, NA_REAL),
    _["LeafArea"] = NumericVector(numCohorts, NA_REAL),
    _["SapwoodArea"] = NumericVector(numCohorts, NA_REAL),
    _["FineRootArea"] = NumericVector(numCohorts, NA_REAL),
    _["HuberValue"] = NumericVector(numCohorts, NA_REAL),
    _["RootAreaLeafArea"] = NumericVector(numCohorts, NA_REAL),
    _["DBH"] = NumericVector(numCohorts, NA_REAL),
    _["Height"] = NumericVector(numCohorts, NA_REAL)
  );
  
  DataFrame growthMortality = DataFrame::create(
    _["SAgrowth"] = NumericVector(numCohorts, NA_REAL),
    _["LAgrowth"] = NumericVector(numCohorts, NA_REAL),
    _["FRAgrowth"] = NumericVector(numCohorts, NA_REAL),
    _["StarvationRate"] = NumericVector(numCohorts, NA_REAL),
    _["DessicationRate"] = NumericVector(numCohorts, NA_REAL),
    _["MortalityRate"] = NumericVector(numCohorts, NA_REAL)
  );

  DataFrame plantBiomassBalance = DataFrame::create(_["InitialDensity"] = NumericVector(numCohorts, 0.0),
                                                    _["InitialSapwoodBiomass"] = NumericVector(numCohorts, 0.0),
                                                    _["InitialStructuralBiomass"] = NumericVector(numCohorts, 0.0),
                                                    _["StructuralBiomassBalance"] = NumericVector(numCohorts, 0.0),
                                                    _["StructuralBiomassChange"] = NumericVector(numCohorts, 0.0),
                                                    _["InitialLabileBiomass"] = NumericVector(numCohorts, 0.0),
                                                    _["LabileBiomassBalance"] = NumericVector(numCohorts, 0.0),
                                                    _["LabileBiomassChange"] = NumericVector(numCohorts, 0.0),
                                                    _["InitialLivingPlantBiomass"] = NumericVector(numCohorts, 0.0),
                                                    _["InitialPlantBiomass"] = NumericVector(numCohorts, 0.0),
                                                    _["PlantBiomassBalance"] = NumericVector(numCohorts, 0.0),
                                                    _["PlantBiomassChange"] = NumericVector(numCohorts, 0.0),
                                                    _["MortalityBiomassLoss"] = NumericVector(numCohorts, 0.0),
                                                    _["InitialCohortBiomass"] = NumericVector(numCohorts, 0.0),
                                                    _["CohortBiomassBalance"] = NumericVector(numCohorts, 0.0),
                                                    _["CohortBiomassChange"] = NumericVector(numCohorts, 0.0));

  List l = List::create(_["topography"] = spwbOut["topography"],
                        _["weather"] = spwbOut["weather"],
                        _["WaterBalance"] = spwbOut["WaterBalance"], 
                        _["EnergyBalance"] = spwbOut["EnergyBalance"],
                        _["CarbonBalance"] = standCB,
                        _["Soil"] = spwbOut["Soil"], 
                        _["Stand"] = spwbOut["Stand"], 
                        _["Plants"] = spwbOut["Plants"],
                        _["LabileCarbonBalance"] = labileCarbonBalance,
                        _["PlantBiomassBalance"] = plantBiomassBalance,
                        _["PlantStructure"] = plantStructure,
                        _["GrowthMortality"] = growthMortality,
                        _["RhizoPsi"] = spwbOut["RhizoPsi"],
                        _["SunlitLeaves"] = spwbOut["SunlitLeaves"],
                        _["ShadeLeaves"] = spwbOut["ShadeLeaves"],
                        _["ExtractionInst"] = spwbOut["ExtractionInst"],
                        _["PlantsInst"] = spwbOut["PlantsInst"],
                        _["SunlitLeavesInst"] = spwbOut["SunlitLeavesInst"],
                        _["ShadeLeavesInst"] = spwbOut["ShadeLeavesInst"]);

    
  //Subdaily output matrices (Sperry/Sureau)
  NumericMatrix LabileCarbonBalanceInst(numCohorts, ntimesteps);  
  NumericMatrix GrossPhotosynthesisInst(numCohorts, ntimesteps);  
  NumericMatrix MaintenanceRespirationInst(numCohorts, ntimesteps);  
  NumericMatrix GrowthCostsInst(numCohorts, ntimesteps);  
  NumericMatrix RootExudationInst(numCohorts, ntimesteps);  
  NumericMatrix PlantSugarTransportInst(numCohorts, ntimesteps);
  NumericMatrix PlantSugarLeafInst(numCohorts, ntimesteps), PlantStarchLeafInst(numCohorts, ntimesteps);
  NumericMatrix PlantSugarSapwoodInst(numCohorts, ntimesteps), PlantStarchSapwoodInst(numCohorts, ntimesteps);
    
  List labileCBInst = List::create(
      _["GrossPhotosynthesis"] = GrossPhotosynthesisInst,
      _["MaintenanceRespiration"] = MaintenanceRespirationInst,
      _["GrowthCosts"] = GrowthCostsInst,
      _["RootExudation"] = RootExudationInst,
      _["LabileCarbonBalance"] = LabileCarbonBalanceInst,
      _["SugarLeaf"] = PlantSugarLeafInst,
      _["StarchLeaf"] = PlantStarchLeafInst,
      _["SugarSapwood"] = PlantSugarSapwoodInst,
      _["StarchSapwood"] = PlantStarchSapwoodInst,
      _["SugarTransport"] = PlantSugarTransportInst
  );
  l.push_back(labileCBInst,"LabileCarbonBalanceInst");
  l.push_back(spwbOut["RadiationInputInst"], "RadiationInputInst");
  l.push_back(spwbOut["LightExtinction"], "LightExtinction");
  l.push_back(spwbOut["LWRExtinction"], "LWRExtinction");
  l.push_back(spwbOut["CanopyTurbulence"], "CanopyTurbulence");
  l.push_back(spwbOut["FireHazard"], "FireHazard");
  l.attr("class") = CharacterVector::create("growth_day","list");
  return(l);
}

DataFrame copyDataFrame(DataFrame comm, int numRows) {
  CharacterVector colnames = comm.attr("names");
  int n = colnames.size();
  List out(n); 
  for(int i = 0;i<n;i++) {
    String nameCol = colnames[i];
    NumericVector vecIn = comm[nameCol];
    NumericVector vecOut(numRows);
    for(int c = 0;c<numRows;c++) {
      vecOut[c] = vecIn[c];
    }
    out[i] = vecOut;
  }
  out.attr("names") = clone(colnames);
  DataFrame dfout(out);
  return(dfout);
}
NumericMatrix copyNumericMatrix(NumericMatrix comm, int rows, int cols) {
  NumericMatrix out(rows, cols);
  for(int r=0;r<rows;r++) {
    for(int c=0;c<cols; c++) {
      out(r, c) = comm(r, c);
    }
  }
  return(out);
}
List copyList(List comm, int n) {
  List out(n);
  for(int i=0;i<n;i++) {
    out[i] = clone(as<List>(comm[i]));
  }
  return(out);
}
List copyCommunicationLongWaveRadiation(List clwr, int ncanlayers) {
  List lwr_struct = List::create(_["LWR_layer"] = copyDataFrame(as<DataFrame>(clwr["LWR_layer"]), ncanlayers),
                                 _["Ldown_ground"] = clwr["Ldown_ground"],
                                 _["Lup_ground"] = clwr["Lup_ground"],
                                 _["Lnet_ground"] = clwr["Lnet_ground"],
                                 _["Ldown_canopy"] = clwr["Ldown_canopy"],
                                 _["Lup_canopy"] = clwr["Lup_canopy"],
                                 _["Lnet_canopy"] = clwr["Lnet_canopy"],
                                 _["Lnet_cohort_layer"] = clwr["Lnet_cohort_layer"]);
  return(lwr_struct);
}
List copyBasicTranspirationOutput(List btc, List x) {
  List control = x["control"];
  String rhizosphereOverlap = control["rhizosphereOverlap"];
  bool plantWaterPools = (rhizosphereOverlap!="total");
  DataFrame cohorts = Rcpp::as<Rcpp::DataFrame>(x["cohorts"]);
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  DataFrame soil = Rcpp::as<Rcpp::DataFrame>(x["soil"]);
  int nlayers = soil.nrow();
  int numCohorts = above.nrow();
  
  NumericMatrix ExtractionComm = btc["Extraction"];
  NumericMatrix Extraction = copyNumericMatrix(ExtractionComm, numCohorts, nlayers); // this is final extraction of each cohort from each layer
  Extraction.attr("dimnames") = List::create(above.attr("row.names"), seq(1,nlayers));
  
  List ExtractionPools(numCohorts);
  List ExtractionPoolsComm = btc["ExtractionPools"];
  if(plantWaterPools) {
    for(int c=0;c<numCohorts;c++) {
      NumericMatrix ExtractionPoolsCohComm = ExtractionPoolsComm[c];
      NumericMatrix ExtractionPoolsCohComm_c = copyNumericMatrix(ExtractionPoolsCohComm, numCohorts, nlayers); // this is final extraction of each cohort from each layer
      ExtractionPoolsCohComm_c.attr("dimnames") = List::create(above.attr("row.names"), seq(1,nlayers));
      ExtractionPools[c] = ExtractionPoolsCohComm_c;
    }
    ExtractionPools.attr("names") = above.attr("row.names");
  }
  
  NumericVector StandComm = btc["Stand"];
  NumericVector Stand = clone(StandComm);
  
  DataFrame PlantsComm = Rcpp::as<Rcpp::DataFrame>(btc["Plants"]);
  DataFrame Plants = copyDataFrame(PlantsComm, numCohorts);
  Plants.attr("row.names") = above.attr("row.names");
  
  List l = List::create(_["cohorts"] = clone(cohorts),
                        _["Stand"] = Stand,
                        _["Plants"] = Plants,
                        _["Extraction"] = Extraction,
                        _["ExtractionPools"] = ExtractionPools);
  return(l);
}
List copyPlantsInstOutput(List PlantsInstComm, List x) {
  List control = x["control"];
  DataFrame cohorts = Rcpp::as<Rcpp::DataFrame>(x["cohorts"]);
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  int numCohorts = cohorts.nrow();
  int ntimesteps = control["ndailysteps"];
  NumericMatrix Einst = copyNumericMatrix(as<NumericMatrix>(PlantsInstComm["E"]), numCohorts, ntimesteps);
  Einst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  NumericMatrix Aginst = copyNumericMatrix(as<NumericMatrix>(PlantsInstComm["Ag"]), numCohorts, ntimesteps);
  Aginst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  NumericMatrix Aninst = copyNumericMatrix(as<NumericMatrix>(PlantsInstComm["An"]), numCohorts, ntimesteps);
  Aninst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  NumericMatrix dEdPInst = copyNumericMatrix(as<NumericMatrix>(PlantsInstComm["dEdP"]), numCohorts, ntimesteps);
  dEdPInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  NumericMatrix LeafPsiInst = copyNumericMatrix(as<NumericMatrix>(PlantsInstComm["LeafPsi"]), numCohorts, ntimesteps);
  LeafPsiInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  NumericMatrix StemPsiInst = copyNumericMatrix(as<NumericMatrix>(PlantsInstComm["StemPsi"]), numCohorts, ntimesteps);
  StemPsiInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  NumericMatrix RootPsiInst = copyNumericMatrix(as<NumericMatrix>(PlantsInstComm["RootPsi"]), numCohorts, ntimesteps);
  RootPsiInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  NumericMatrix LeafSympPsiInst = copyNumericMatrix(as<NumericMatrix>(PlantsInstComm["LeafSympPsi"]), numCohorts, ntimesteps);
  LeafSympPsiInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  NumericMatrix StemSympPsiInst = copyNumericMatrix(as<NumericMatrix>(PlantsInstComm["StemSympPsi"]), numCohorts, ntimesteps);
  StemSympPsiInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  NumericMatrix LeafSympRWCInst = copyNumericMatrix(as<NumericMatrix>(PlantsInstComm["LeafSympRWC"]), numCohorts, ntimesteps);
  LeafSympRWCInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  NumericMatrix StemSympRWCInst = copyNumericMatrix(as<NumericMatrix>(PlantsInstComm["StemSympRWC"]), numCohorts, ntimesteps);
  StemSympRWCInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  NumericMatrix StemPLC = copyNumericMatrix(as<NumericMatrix>(PlantsInstComm["StemPLC"]), numCohorts, ntimesteps);
  StemPLC.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  NumericMatrix LeafPLC = copyNumericMatrix(as<NumericMatrix>(PlantsInstComm["LeafPLC"]), numCohorts, ntimesteps);
  LeafPLC.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  NumericMatrix LeafRWCInst = copyNumericMatrix(as<NumericMatrix>(PlantsInstComm["LeafRWC"]), numCohorts, ntimesteps);
  LeafRWCInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  NumericMatrix StemRWCInst = copyNumericMatrix(as<NumericMatrix>(PlantsInstComm["StemRWC"]), numCohorts, ntimesteps);
  StemRWCInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  NumericMatrix PWBinst = copyNumericMatrix(as<NumericMatrix>(PlantsInstComm["PWB"]), numCohorts, ntimesteps);
  PWBinst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  
  List PlantsInst = List::create(
    _["E"]=Einst, _["Ag"]=Aginst, _["An"]=Aninst,
    _["dEdP"] = dEdPInst,
    _["RootPsi"] = RootPsiInst, 
    _["StemPsi"] = StemPsiInst,
    _["LeafPsi"] = LeafPsiInst,
    _["StemSympPsi"] = StemSympPsiInst,
    _["LeafSympPsi"] = LeafSympPsiInst,
    _["StemPLC"] = StemPLC, 
    _["LeafPLC"] = LeafPLC, 
    _["StemRWC"] = StemRWCInst,
    _["LeafRWC"] = LeafRWCInst,
    _["StemSympRWC"] = StemSympRWCInst,
    _["LeafSympRWC"] = LeafSympRWCInst,
    _["PWB"] = PWBinst);
  return(PlantsInst);
}
List copyLeavesInstOutput(List LeavesInstComm, List x) {
  List control = x["control"];
  DataFrame cohorts = Rcpp::as<Rcpp::DataFrame>(x["cohorts"]);
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  int numCohorts = cohorts.nrow();
  int ntimesteps = control["ndailysteps"];
  
  NumericMatrix LAI = copyNumericMatrix(as<NumericMatrix>(LeavesInstComm["LAI"]), numCohorts, ntimesteps);
  LAI.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  NumericMatrix Vmax298 = copyNumericMatrix(as<NumericMatrix>(LeavesInstComm["Vmax298"]), numCohorts, ntimesteps);
  Vmax298.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  NumericMatrix Jmax298 = copyNumericMatrix(as<NumericMatrix>(LeavesInstComm["Jmax298"]), numCohorts, ntimesteps);
  Jmax298.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  NumericMatrix SWR = copyNumericMatrix(as<NumericMatrix>(LeavesInstComm["Abs_SWR"]), numCohorts, ntimesteps);
  SWR.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  NumericMatrix PAR = copyNumericMatrix(as<NumericMatrix>(LeavesInstComm["Abs_PAR"]), numCohorts, ntimesteps);
  PAR.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  NumericMatrix LWR = copyNumericMatrix(as<NumericMatrix>(LeavesInstComm["Net_LWR"]), numCohorts, ntimesteps);
  LWR.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  NumericMatrix An = copyNumericMatrix(as<NumericMatrix>(LeavesInstComm["An"]), numCohorts, ntimesteps);
  An.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  NumericMatrix Ag = copyNumericMatrix(as<NumericMatrix>(LeavesInstComm["Ag"]), numCohorts, ntimesteps);
  Ag.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  NumericMatrix Ci = copyNumericMatrix(as<NumericMatrix>(LeavesInstComm["Ci"]), numCohorts, ntimesteps);
  Ci.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  NumericMatrix E = copyNumericMatrix(as<NumericMatrix>(LeavesInstComm["E"]), numCohorts, ntimesteps);
  E.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  NumericMatrix GSW = copyNumericMatrix(as<NumericMatrix>(LeavesInstComm["Gsw"]), numCohorts, ntimesteps);
  GSW.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  NumericMatrix VPD = copyNumericMatrix(as<NumericMatrix>(LeavesInstComm["VPD"]), numCohorts, ntimesteps);
  VPD.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  NumericMatrix Temp = copyNumericMatrix(as<NumericMatrix>(LeavesInstComm["Temp"]), numCohorts, ntimesteps);
  Temp.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  NumericMatrix Psi = copyNumericMatrix(as<NumericMatrix>(LeavesInstComm["Psi"]), numCohorts, ntimesteps);
  Psi.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));

  
  List LeavesInst = List::create(
    _["LAI"] = LAI,
    _["Vmax298"] = Vmax298,
    _["Jmax298"] = Jmax298,
    _["Abs_SWR"]=SWR,
    _["Abs_PAR"]=PAR,
    _["Net_LWR"] = LWR,
    _["Ag"] = Ag,
    _["An"] = An,
    _["Ci"] = Ci,
    _["E"] = E,
    _["Gsw"] = GSW,
    _["VPD"] = VPD,
    _["Temp"] = Temp,
    _["Psi"] = Psi);
  
  return(LeavesInst);
}
List copyEnergyBalanceOutput(List EnergyBalanceComm, List x) {
  List control = x["control"];
  int ntimesteps = control["ndailysteps"];
  DataFrame canopyParams = Rcpp::as<Rcpp::DataFrame>(x["canopy"]);
  int ncanlayers = canopyParams.nrow(); //Number of canopy layers
  DataFrame soil = Rcpp::as<Rcpp::DataFrame>(x["soil"]);
  int nlayers = soil.nrow();
  
  DataFrame Tinst = copyDataFrame(as<DataFrame>(EnergyBalanceComm["Temperature"]), ntimesteps);
  DataFrame CEBinst = copyDataFrame(as<DataFrame>(EnergyBalanceComm["CanopyEnergyBalance"]), ntimesteps);
  DataFrame SEBinst = copyDataFrame(as<DataFrame>(EnergyBalanceComm["SoilEnergyBalance"]), ntimesteps);
  NumericMatrix Tcan_mat= copyNumericMatrix(as<NumericMatrix>(EnergyBalanceComm["TemperatureLayers"]), ntimesteps, ncanlayers);
  Tcan_mat.attr("dimnames") = List::create(seq(1,ntimesteps), seq(1,ncanlayers));
  NumericMatrix VPcan_mat= copyNumericMatrix(as<NumericMatrix>(EnergyBalanceComm["VaporPressureLayers"]), ntimesteps, ncanlayers);
  VPcan_mat.attr("dimnames") = List::create(seq(1,ntimesteps), seq(1,ncanlayers));
  NumericMatrix Tsoil_mat= copyNumericMatrix(as<NumericMatrix>(EnergyBalanceComm["SoilTemperature"]), ntimesteps, nlayers);
  Tsoil_mat.attr("dimnames") = List::create(seq(1,ntimesteps), seq(1,nlayers));
  List EnergyBalance = List::create(_["Temperature"]=Tinst, 
                                    _["SoilTemperature"] = Tsoil_mat,
                                    _["CanopyEnergyBalance"] = CEBinst, 
                                    _["SoilEnergyBalance"] = SEBinst,
                                    _["TemperatureLayers"] = Tcan_mat, 
                                    _["VaporPressureLayers"] = VPcan_mat);
  return(EnergyBalance);
}
List copyAdvancedTranspirationOutput(List atc, List x) {
  List control = x["control"];
  int ntimesteps = control["ndailysteps"];
  String transpirationMode = control["transpirationMode"];
  String rhizosphereOverlap = control["rhizosphereOverlap"];
  bool plantWaterPools = (rhizosphereOverlap!="total");
  DataFrame cohorts = Rcpp::as<Rcpp::DataFrame>(x["cohorts"]);
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  DataFrame soil = Rcpp::as<Rcpp::DataFrame>(x["soil"]);
  int nlayers = soil.nrow();
  int numCohorts = above.nrow();
  DataFrame canopyParams = Rcpp::as<Rcpp::DataFrame>(x["canopy"]);
  int ncanlayers = canopyParams.nrow(); //Number of canopy layers
  
  NumericMatrix SoilWaterExtractComm = atc["Extraction"];
  NumericMatrix SoilWaterExtract = copyNumericMatrix(SoilWaterExtractComm, numCohorts, nlayers);
  SoilWaterExtract.attr("dimnames") = List::create(above.attr("row.names"), seq(1,nlayers));
  
  List ExtractionPools(numCohorts);
  List ExtractionPoolsComm = atc["ExtractionPools"];
  if(plantWaterPools) {
    for(int c=0;c<numCohorts;c++) {
      NumericMatrix ExtractionPoolsCohComm = ExtractionPoolsComm[c];
      NumericMatrix ExtractionPoolsCohComm_c = copyNumericMatrix(ExtractionPoolsCohComm, numCohorts, nlayers); // this is final extraction of each cohort from each layer
      ExtractionPoolsCohComm_c.attr("dimnames") = List::create(above.attr("row.names"), seq(1,nlayers));
      ExtractionPools[c] = ExtractionPoolsCohComm_c;
    }
    ExtractionPools.attr("names") = above.attr("row.names");
  }
  
  NumericMatrix soilLayerExtractInstComm = atc["ExtractionInst"];
  NumericMatrix soilLayerExtractInst = copyNumericMatrix(soilLayerExtractInstComm, nlayers, ntimesteps);
  soilLayerExtractInst.attr("dimnames") = List::create(seq(1,nlayers), seq(1,ntimesteps));
  
  NumericMatrix minPsiRhizoComm = atc["RhizoPsi"];
  NumericMatrix minPsiRhizo = copyNumericMatrix(minPsiRhizoComm, numCohorts, nlayers);
  minPsiRhizo.attr("dimnames") = List::create(above.attr("row.names"), seq(1,nlayers));
  
  NumericVector StandComm = atc["Stand"];
  NumericVector Stand = clone(StandComm);
  
  List EnergyBalance = copyEnergyBalanceOutput(as<List>(atc["EnergyBalance"]), x);

  DataFrame Plants = copyDataFrame(as<DataFrame>(atc["Plants"]), numCohorts);
  Plants.attr("row.names") = above.attr("row.names");
  
  DataFrame Sunlit = copyDataFrame(as<DataFrame>(atc["SunlitLeaves"]), numCohorts);
  Sunlit.attr("row.names") = above.attr("row.names");
  DataFrame Shade = copyDataFrame(as<DataFrame>(atc["ShadeLeaves"]), numCohorts);
  Shade.attr("row.names") = above.attr("row.names");
  
  List PlantsInst = copyPlantsInstOutput(as<List>(atc["PlantsInst"]), x);
  
  List SunlitInst = copyLeavesInstOutput(as<List>(atc["SunlitLeavesInst"]), x);
  List ShadeInst = copyLeavesInstOutput(as<List>(atc["ShadeLeavesInst"]), x);
  
  List lwrExtinctionListComm = atc["LWRExtinction"];
  List lwrExtinctionList(ntimesteps);
  for(int n=0;n<ntimesteps;n++) {
    lwrExtinctionList[n] = copyCommunicationLongWaveRadiation(as<List>(lwrExtinctionListComm[n]), ncanlayers);
  }

  List l = List::create(_["cohorts"] = clone(cohorts),
                        _["EnergyBalance"] = EnergyBalance,
                        _["Extraction"] = SoilWaterExtract,
                        _["ExtractionPools"] = ExtractionPools,
                        _["RhizoPsi"] = minPsiRhizo,
                        _["Stand"] = Stand,
                        _["Plants"] = Plants,
                        _["SunlitLeaves"] = Sunlit,
                        _["ShadeLeaves"] = Shade,
                        _["ExtractionInst"] = soilLayerExtractInst,
                        _["RadiationInputInst"] = clone(as<List>(atc["RadiationInputInst"])), 
                        _["PlantsInst"] = PlantsInst,
                        _["SunlitLeavesInst"] = SunlitInst,
                        _["ShadeLeavesInst"] = ShadeInst,
                        _["LightExtinction"] =  clone(as<List>(atc["LightExtinction"])), 
                        _["LWRExtinction"] = lwrExtinctionList,
                        _["CanopyTurbulence"] = copyDataFrame(as<DataFrame>(atc["CanopyTurbulence"]), ncanlayers)); //To be replaced
  
  List supply = copyList(as<List>(atc["SupplyFunctions"]), numCohorts);
  supply.attr("names") = above.attr("row.names");
  List outPhotoSunlit = copyList(as<List>(atc["PhotoSunlitFunctions"]), numCohorts);
  outPhotoSunlit.attr("names") = above.attr("row.names");
  List outPhotoShade = copyList(as<List>(atc["PhotoShadeFunctions"]), numCohorts);
  outPhotoShade.attr("names") = above.attr("row.names");
  List outPMSunlit = copyList(as<List>(atc["PMSunlitFunctions"]), numCohorts);
  outPMSunlit.attr("names") = above.attr("row.names");
  List outPMShade = copyList(as<List>(atc["PMShadeFunctions"]), numCohorts);
  outPMShade.attr("names") = above.attr("row.names");
  l.push_back(supply, "SupplyFunctions");
  l.push_back(outPhotoSunlit, "PhotoSunlitFunctions");
  l.push_back(outPhotoShade, "PhotoShadeFunctions");
  l.push_back(outPMSunlit, "PMSunlitFunctions");
  l.push_back(outPMShade, "PMShadeFunctions");
  
  l.attr("class") = CharacterVector::create("pwb_day","list");
  return(l);
}


List copyBasicSPWBOutput(List boc, List x) {
  List control = x["control"];
  DataFrame cohorts = Rcpp::as<Rcpp::DataFrame>(x["cohorts"]);
  DataFrame soil = Rcpp::as<Rcpp::DataFrame>(x["soil"]);
  int nlayers = soil.nrow();
  int numCohorts = cohorts.nrow();
  
  NumericVector topo = clone(as<NumericVector>(boc["topography"]));
  NumericVector meteovec_bas = clone(as<NumericVector>(boc["weather"]));
  NumericVector WaterBalance = clone(as<NumericVector>(boc["WaterBalance"]));
  
  List l = List::create(_["cohorts"] = clone(cohorts),
                        _["topography"] = topo,
                        _["weather"] = meteovec_bas,
                        _["WaterBalance"] = WaterBalance);

  if(control["soilResults"]) {
    DataFrame Soil = copyDataFrame(as<DataFrame>(boc["Soil"]), nlayers);
    l.push_back(Soil, "Soil");
  }  
  if(control["standResults"]) {
    NumericVector Stand = clone(as<NumericVector>(boc["Stand"]));
    l.push_back(Stand, "Stand");
  }
  if(control["plantResults"]) {
    DataFrame PlantsComm = Rcpp::as<Rcpp::DataFrame>(boc["Plants"]);
    DataFrame Plants = copyDataFrame(PlantsComm, numCohorts);
    Plants.attr("row.names") = cohorts.attr("row.names");
    l.push_back(Plants, "Plants");
  }
  if(control["fireHazardResults"]) {
    l.push_back(clone(as<NumericVector>(boc["FireHazard"])), "FireHazard"); 
  }
  l.attr("class") = CharacterVector::create("spwb_day","list");
  return(l);
}
List copyBasicGROWTHOutput(List boc, List x) {
  List control = x["control"];
  DataFrame cohorts = Rcpp::as<Rcpp::DataFrame>(x["cohorts"]);
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  int numCohorts = cohorts.nrow();
  
  List spwbOut = copyBasicSPWBOutput(boc, x);
  
  
  NumericVector standCB = clone(as<NumericVector>(boc["CarbonBalance"]));
    
  
  List l = List::create(_["cohorts"] = spwbOut["cohorts"],
                        _["topography"] = spwbOut["topography"],
                        _["weather"] = spwbOut["weather"],
                        _["WaterBalance"] = spwbOut["WaterBalance"], 
                        _["CarbonBalance"] = standCB);
  if(control["soilResults"]) {
    l.push_back(spwbOut["Soil"], "Soil");
  }
  if(control["standResults"]) {
    l.push_back(spwbOut["Stand"], "Stand");
  }
  if(control["plantResults"]) {
    l.push_back(spwbOut["Plants"], "Plants");
  }
  if(control["labileCarbonBalanceResults"]) {
    DataFrame labileCarbonBalance = copyDataFrame(as<DataFrame>(boc["LabileCarbonBalance"]), numCohorts);
    labileCarbonBalance.attr("row.names") = above.attr("row.names");
    l.push_back(labileCarbonBalance, "LabileCarbonBalance");
  }
  DataFrame plantBiomassBalance = copyDataFrame(as<DataFrame>(boc["PlantBiomassBalance"]), numCohorts);
  plantBiomassBalance.attr("row.names") = above.attr("row.names");
  l.push_back(plantBiomassBalance, "PlantBiomassBalance");
  
  if(control["plantStructureResults"]) {
    DataFrame plantStructure = copyDataFrame(as<DataFrame>(boc["PlantStructure"]), numCohorts);
    plantStructure.attr("row.names") = above.attr("row.names");
    l.push_back(plantStructure, "PlantStructure");
  }
  
  if(control["growthMortalityResults"]) {
    DataFrame growthMortality  = copyDataFrame(as<DataFrame>(boc["GrowthMortality"]), numCohorts);
    growthMortality.attr("row.names") = above.attr("row.names");
    l.push_back(growthMortality, "GrowthMortality");
  }
  
  if(control["fireHazardResults"]) l.push_back(spwbOut["FireHazard"], "FireHazard");
  l.attr("class") = CharacterVector::create("growth_day","list");
  return(l);
}

List copyAdvancedSPWBOutput(List aoc, List x) {
  List control = x["control"];
  DataFrame cohorts = Rcpp::as<Rcpp::DataFrame>(x["cohorts"]);
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  DataFrame soil = Rcpp::as<Rcpp::DataFrame>(x["soil"]);
  DataFrame canopyParams = Rcpp::as<Rcpp::DataFrame>(x["canopy"]);
  int ncanlayers = canopyParams.nrow(); //Number of canopy layers
  int nlayers = soil.nrow();
  int numCohorts = cohorts.nrow();
  int ntimesteps = control["ndailysteps"];
  
  NumericVector topo = clone(as<NumericVector>(aoc["topography"]));
  NumericVector meteovec_adv = clone(as<NumericVector>(aoc["weather"]));
  
  NumericVector WaterBalance = clone(as<NumericVector>(aoc["WaterBalance"]));
  List EnergyBalance = copyEnergyBalanceOutput(as<List>(aoc["EnergyBalance"]), x);
  List l = List::create(_["cohorts"] = clone(cohorts),
                        _["topography"] = topo,
                        _["weather"] = meteovec_adv,
                        _["WaterBalance"] = WaterBalance,
                        _["EnergyBalance"] = EnergyBalance);

  if(control["soilResults"]){
    DataFrame Soil = copyDataFrame(as<DataFrame>(aoc["Soil"]), nlayers);
    l.push_back(Soil, "Soil");
  }
  if(control["standResults"]){
    NumericVector Stand = clone(as<NumericVector>(aoc["Stand"]));
    l.push_back(Stand, "Stand");
  }
  if(control["plantResults"]){
    DataFrame Plants = copyDataFrame(as<DataFrame>(aoc["Plants"]), numCohorts);
    Plants.attr("row.names") = cohorts.attr("row.names");
    l.push_back(Plants, "Plants"); 
  }
  if(control["leafResults"]){
    DataFrame Sunlit = copyDataFrame(as<DataFrame>(aoc["SunlitLeaves"]), numCohorts);
    Sunlit.attr("row.names") = above.attr("row.names");
    DataFrame Shade = copyDataFrame(as<DataFrame>(aoc["ShadeLeaves"]), numCohorts);
    Shade.attr("row.names") = above.attr("row.names");
    l.push_back(Sunlit, "SunlitLeaves");
    l.push_back(Shade, "ShadeLeaves");
  }

  NumericMatrix minPsiRhizoComm = aoc["RhizoPsi"];
  NumericMatrix minPsiRhizo = copyNumericMatrix(minPsiRhizoComm, numCohorts, nlayers);
  minPsiRhizo.attr("dimnames") = List::create(above.attr("row.names"), seq(1,nlayers));
  l.push_back(minPsiRhizo, "RhizoPsi");

  
  if(control["subdailyResults"]) {
    NumericMatrix soilLayerExtractInstComm = aoc["ExtractionInst"];
    NumericMatrix soilLayerExtractInst = copyNumericMatrix(soilLayerExtractInstComm, nlayers, ntimesteps);
    soilLayerExtractInst.attr("dimnames") = List::create(seq(1,nlayers), seq(1,ntimesteps));
    l.push_back(soilLayerExtractInst, "ExtractionInst");
    List PlantsInst = copyPlantsInstOutput(as<List>(aoc["PlantsInst"]), x);
    l.push_back(PlantsInst, "PlantsInst");
    l.push_back(clone(as<List>(aoc["RadiationInputInst"])), "RadiationInputInst");
    List SunlitInst = copyLeavesInstOutput(as<List>(aoc["SunlitLeavesInst"]), x);
    List ShadeInst = copyLeavesInstOutput(as<List>(aoc["ShadeLeavesInst"]), x);
    l.push_back(SunlitInst, "SunlitLeavesInst");
    l.push_back(ShadeInst, "ShadeLeavesInst");
  }
  List lwrExtinctionListComm = aoc["LWRExtinction"];
  List lwrExtinctionList(ntimesteps);
  for(int n=0;n<ntimesteps;n++) {
    lwrExtinctionList[n] = copyCommunicationLongWaveRadiation(as<List>(lwrExtinctionListComm[n]), ncanlayers);
  }
  l.push_back(clone(as<List>(aoc["LightExtinction"])), "LightExtinction");
  l.push_back(lwrExtinctionList, "LWRExtinction");
  l.push_back(copyDataFrame(as<DataFrame>(aoc["CanopyTurbulence"]), ncanlayers), "CanopyTurbulence");
  l.attr("class") = CharacterVector::create("spwb_day","list");
  return(l);
}

List copyAdvancedGROWTHOutput(List aoc, List x) {
  List control = x["control"];
  DataFrame cohorts = Rcpp::as<Rcpp::DataFrame>(x["cohorts"]);
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  int numCohorts = cohorts.nrow();
  int ntimesteps = control["ndailysteps"];
  bool subdailyCarbonBalance = control["subdailyCarbonBalance"];

  List spwbOut = copyAdvancedSPWBOutput(aoc, x);

  NumericVector standCB = clone(as<NumericVector>(aoc["CarbonBalance"]));
  
  List l = List::create(_["cohorts"] = spwbOut["cohorts"],
                        _["topography"] = spwbOut["topography"],
                        _["weather"] = spwbOut["weather"],
                        _["WaterBalance"] = spwbOut["WaterBalance"], 
                        _["EnergyBalance"] = spwbOut["EnergyBalance"],
                        _["CarbonBalance"] = standCB);
  
  if(control["soilResults"]) {
    l.push_back(spwbOut["Soil"], "Soil");
  }
  if(control["standResults"]){
    l.push_back(spwbOut["Stand"], "Stand");
  }
  if(control["plantResults"]){
    l.push_back(spwbOut["Plants"], "Plants");
  }
  if(control["labileCarbonBalanceResults"]) {
    DataFrame labileCarbonBalance = copyDataFrame(as<DataFrame>(aoc["LabileCarbonBalance"]), numCohorts);
    labileCarbonBalance.attr("row.names") = above.attr("row.names");
    l.push_back(labileCarbonBalance, "LabileCarbonBalance");
  }
  DataFrame plantBiomassBalance = copyDataFrame(as<DataFrame>(aoc["PlantBiomassBalance"]), numCohorts);
  plantBiomassBalance.attr("row.names") = above.attr("row.names");
  l.push_back(plantBiomassBalance, "PlantBiomassBalance");
  
  if(control["plantStructureResults"]) {
    DataFrame plantStructure = copyDataFrame(as<DataFrame>(aoc["PlantStructure"]), numCohorts);
    plantStructure.attr("row.names") = above.attr("row.names");
    l.push_back(plantStructure, "PlantStructure");
  }
  if(control["growthMortalityResults"]) {
    DataFrame growthMortality  = copyDataFrame(as<DataFrame>(aoc["GrowthMortality"]), numCohorts);
    growthMortality.attr("row.names") = above.attr("row.names");
    l.push_back(growthMortality, "GrowthMortality");
  }
  
  l.push_back(spwbOut["RhizoPsi"], "RhizoPsi");
  
  if(control["leafResults"]){
    l.push_back(spwbOut["SunlitLeaves"], "SunlitLeaves");
    l.push_back(spwbOut["ShadeLeaves"], "ShadeLeaves");
  }
  if(control["subdailyResults"]) {
    l.push_back(spwbOut["ExtractionInst"], "ExtractionInst");
    l.push_back(spwbOut["PlantsInst"], "PlantsInst");
    l.push_back(spwbOut["RadiationInputInst"], "RadiationInputInst");
    l.push_back(spwbOut["SunlitLeavesInst"], "SunlitLeavesInst");
    l.push_back(spwbOut["ShadeLeavesInst"], "ShadeLeavesInst");
  }
  if(subdailyCarbonBalance) {
    List labileCBInstComm = aoc["LabileCarbonBalanceInst"];
    
    NumericMatrix LabileCarbonBalanceInst = copyNumericMatrix(as<NumericMatrix>(labileCBInstComm["LabileCarbonBalance"]), numCohorts, ntimesteps);  
    NumericMatrix GrossPhotosynthesisInst = copyNumericMatrix(as<NumericMatrix>(labileCBInstComm["GrossPhotosynthesis"]), numCohorts, ntimesteps);  
    NumericMatrix MaintenanceRespirationInst = copyNumericMatrix(as<NumericMatrix>(labileCBInstComm["MaintenanceRespiration"]), numCohorts, ntimesteps);    
    NumericMatrix GrowthCostsInst = copyNumericMatrix(as<NumericMatrix>(labileCBInstComm["GrowthCosts"]), numCohorts, ntimesteps);    
    NumericMatrix RootExudationInst = copyNumericMatrix(as<NumericMatrix>(labileCBInstComm["RootExudation"]), numCohorts, ntimesteps);    
    NumericMatrix PlantSugarLeafInst = copyNumericMatrix(as<NumericMatrix>(labileCBInstComm["SugarLeaf"]), numCohorts, ntimesteps);  
    NumericMatrix PlantStarchLeafInst = copyNumericMatrix(as<NumericMatrix>(labileCBInstComm["StarchLeaf"]), numCohorts, ntimesteps);  
    NumericMatrix PlantSugarSapwoodInst= copyNumericMatrix(as<NumericMatrix>(labileCBInstComm["SugarSapwood"]), numCohorts, ntimesteps);  
    NumericMatrix PlantStarchSapwoodInst = copyNumericMatrix(as<NumericMatrix>(labileCBInstComm["StarchSapwood"]), numCohorts, ntimesteps);  
    NumericMatrix PlantSugarTransportInst = copyNumericMatrix(as<NumericMatrix>(labileCBInstComm["SugarTransport"]), numCohorts, ntimesteps);  
    
    GrossPhotosynthesisInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
    MaintenanceRespirationInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
    LabileCarbonBalanceInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
    GrowthCostsInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
    PlantSugarLeafInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
    PlantStarchLeafInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
    PlantSugarSapwoodInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
    PlantStarchSapwoodInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
    PlantSugarTransportInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
    RootExudationInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
    
    List labileCBInst = List::create(
      _["GrossPhotosynthesis"] = GrossPhotosynthesisInst,
      _["MaintenanceRespiration"] = MaintenanceRespirationInst,
      _["GrowthCosts"] = GrowthCostsInst,
      _["RootExudation"] = RootExudationInst,
      _["LabileCarbonBalance"] = LabileCarbonBalanceInst,
      _["SugarLeaf"] = PlantSugarLeafInst,
      _["StarchLeaf"] = PlantStarchLeafInst,
      _["SugarSapwood"] = PlantSugarSapwoodInst,
      _["StarchSapwood"] = PlantStarchSapwoodInst,
      _["SugarTransport"] = PlantSugarTransportInst
    );
    l.push_back(labileCBInst,"LabileCarbonBalanceInst");
  }
  l.push_back(spwbOut["LightExtinction"], "LightExtinction");
  l.push_back(spwbOut["CanopyTurbulence"], "CanopyTurbulence");
  if(control["fireHazardResults"]) l.push_back(spwbOut["FireHazard"], "FireHazard");
  l.attr("class") = CharacterVector::create("growth_day","list");
  return(l);
}

List copySPWBOutput(List internalCommunication, List x) {
  List control = x["control"];
  String transpirationMode = control["transpirationMode"];
  List modelOutput;
  if(transpirationMode =="Granier") {
    List modelOutputComm = internalCommunication["basicSPWBOutput"];
    modelOutput = copyBasicSPWBOutput(modelOutputComm, x);
  } else {
    List modelOutputComm = internalCommunication["advancedSPWBOutput"];
    modelOutput = copyAdvancedSPWBOutput(modelOutputComm, x);
  }
  return(modelOutput);
}
List copyGROWTHOutput(List internalCommunication, List x) {
  List control = x["control"];
  String transpirationMode = control["transpirationMode"];
  List modelOutput;
  if(transpirationMode =="Granier") {
    List modelOutputComm = internalCommunication["basicGROWTHOutput"];
    modelOutput = copyBasicGROWTHOutput(modelOutputComm, x);
  } else {
    List modelOutputComm = internalCommunication["advancedGROWTHOutput"];
    modelOutput = copyAdvancedGROWTHOutput(modelOutputComm, x);
  }
  return(modelOutput);
}

DataFrame communicationCarbonCompartments(int numCohorts) {
  DataFrame df = DataFrame::create(
    _["LeafStorageVolume"] = NumericVector(numCohorts, NA_REAL),
    _["SapwoodStorageVolume"] = NumericVector(numCohorts, NA_REAL),
    _["LeafStarchMaximumConcentration"] = NumericVector(numCohorts, NA_REAL),
    _["SapwoodStarchMaximumConcentration"] = NumericVector(numCohorts, NA_REAL),
    _["LeafStarchCapacity"] = NumericVector(numCohorts, NA_REAL),
    _["SapwoodStarchCapacity"] = NumericVector(numCohorts, NA_REAL),
    _["LeafStructuralBiomass"] = NumericVector(numCohorts, NA_REAL),
    _["SapwoodStructuralBiomass"] = NumericVector(numCohorts, NA_REAL),
    _["SapwoodLivingStructuralBiomass"] = NumericVector(numCohorts, NA_REAL),
    _["FineRootBiomass"] = NumericVector(numCohorts, NA_REAL),
    _["StructuralBiomass"] = NumericVector(numCohorts, NA_REAL),
    _["LabileBiomass"] = NumericVector(numCohorts, NA_REAL),
    _["TotalLivingBiomass"] = NumericVector(numCohorts, NA_REAL),
    _["TotalBiomass"] = NumericVector(numCohorts, NA_REAL)
  );
  return(df);
}

List communicationInitialFinalCarbonCompartments(int numCohorts) {
  DataFrame ccFin_g_ind = communicationCarbonCompartments(numCohorts);
  DataFrame ccIni_g_ind = clone(ccFin_g_ind);
  List l = List::create(_["ccIni_g_ind"] = ccIni_g_ind,
                        _["ccFin_g_ind"] = ccFin_g_ind);
  return(l);
}

//' Internal communication
//'
//' Functions for internal communication. Not to be called by users.
//' 
//' @param internalCommunication List for internal communication.
//' @param x An object of class \code{\link{spwbInput}} or \code{\link{growthInput}}.
//' @param date Date as string "yyyy-mm-dd".
//' @param meteovec A named numerical vector with weather data. See variable names in parameter \code{meteo} of \code{\link{spwb}}.
//' @param latitude Latitude (in degrees).
//' @param elevation,slope,aspect Elevation above sea level (in m), slope (in degrees) and aspect (in degrees from North). 
//' @param runon Surface water amount running on the target area from upslope (in mm).
//' @param lateralFlows Lateral source/sink terms for each soil layer (interflow/to from adjacent locations) as mm/day.
//' @param waterTableDepth Water table depth (in mm). When not missing, capillarity rise will be allowed if lower than total soil depth.
//' @param modifyInput Boolean flag to indicate that the input \code{x} object is allowed to be modified during the simulation.
//' @param model String for model, either "spwb" or "growth".
//' 
//' @name communication
//' @keywords internal
// [[Rcpp::export("copy_model_output")]]
List copyModelOutput(List internalCommunication, List x, String model) {
  List out;
  if(model=="spwb") {
    out= copySPWBOutput(internalCommunication, x);
  } else {
    out= copyGROWTHOutput(internalCommunication, x);
  }
  return(out);
}

// General communication structures for many inputs
//' @rdname communication
//' @keywords internal
// [[Rcpp::export("general_communication_structures")]]
List generalCommunicationStructures(int numCohorts, int nlayers, int ncanlayers, int ntimesteps,
                                    String model) {
  List SWBcommunication = communicationSoilWaterBalance(nlayers);
  List SEBcommunication = communicationSoilEnergyBalance(nlayers);
  List basicTranspirationOutput = basicTranspirationCommunicationOutput(numCohorts, nlayers);
  List advancedTranspirationOutput = advancedTranspirationCommunicationOutput(numCohorts, nlayers, ncanlayers, ntimesteps);
  List basicSPWBOutput = basicSPWBCommunicationOutput(basicTranspirationOutput, nlayers);
  List advancedSPWBOutput = advancedSPWBCommunicationOutput(advancedTranspirationOutput, nlayers);
  
  List ic = List::create(_["SWBcommunication"] = SWBcommunication,
                         _["SEBcommunication"] = SEBcommunication,
                         _["basicTranspirationOutput"] = basicTranspirationOutput,
                         _["advancedTranspirationOutput"] = advancedTranspirationOutput,
                         _["basicSPWBOutput"] = basicSPWBOutput,
                         _["advancedSPWBOutput"] = advancedSPWBOutput);
  if(model=="growth") {
    List basicGROWTHOutput = basicGROWTHCommunicationOutput(basicSPWBOutput, numCohorts, nlayers);
    List advancedGROWTHOutput = advancedGROWTHCommunicationOutput(advancedSPWBOutput, numCohorts, ntimesteps);
    List initialFinalCC = communicationInitialFinalCarbonCompartments(numCohorts);
    ic.push_back(basicGROWTHOutput, "basicGROWTHOutput");
    ic.push_back(advancedGROWTHOutput, "advancedGROWTHOutput");
    ic.push_back(initialFinalCC, "initialFinalCC");
  }
  return(ic);
}



// Communication structures fitted to the model input 
//' @rdname communication
//' @keywords internal
// [[Rcpp::export("instance_communication_structures")]]
List instanceCommunicationStructures(List x, String model) {
  List control = x["control"];
  DataFrame cohorts = Rcpp::as<Rcpp::DataFrame>(x["cohorts"]);
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  DataFrame soil = Rcpp::as<Rcpp::DataFrame>(x["soil"]);
  DataFrame canopyParams = Rcpp::as<Rcpp::DataFrame>(x["canopy"]);
  int ncanlayers = canopyParams.nrow(); //Number of canopy layers
  int nlayers = soil.nrow();
  int numCohorts = cohorts.nrow();
  int ntimesteps = control["ndailysteps"];
  return(generalCommunicationStructures(numCohorts, nlayers, ncanlayers, ntimesteps, model));
}

