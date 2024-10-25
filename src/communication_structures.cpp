#include <Rcpp.h>
using namespace Rcpp;

const int MAX_SOIL_LAYERS = 10;
const int MAX_CANOPY_LAYERS = 100;
const int MAX_NUM_COHORTS = 100;
const int MAX_NUM_TIMESTEPS = 100;

List communicationLongWaveRadiation() {
  NumericVector Lup(MAX_CANOPY_LAYERS, NA_REAL), Ldown(MAX_CANOPY_LAYERS, NA_REAL), Lnet(MAX_CANOPY_LAYERS, NA_REAL);
  NumericVector tau(MAX_CANOPY_LAYERS, NA_REAL), sumTauComp(MAX_CANOPY_LAYERS, NA_REAL);
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

DataFrame communicationCanopyTurbulence() {
  DataFrame output = DataFrame::create(Named("zmid") = NumericVector(MAX_CANOPY_LAYERS, NA_REAL),
                                       Named("u") = NumericVector(MAX_CANOPY_LAYERS, NA_REAL),
                                       Named("du") = NumericVector(MAX_CANOPY_LAYERS, NA_REAL),
                                       Named("epsilon") = NumericVector(MAX_CANOPY_LAYERS, NA_REAL),
                                       Named("k") = NumericVector(MAX_CANOPY_LAYERS, NA_REAL),
                                       Named("uw") = NumericVector(MAX_CANOPY_LAYERS, NA_REAL));
  return(output);
}
List basicTranspirationCommunicationOutput() {

  NumericMatrix Extraction(MAX_NUM_COHORTS, MAX_SOIL_LAYERS); // this is final extraction of each cohort from each layer

  List ExtractionPools(MAX_NUM_COHORTS);
  for(int c=0;c<MAX_NUM_COHORTS;c++) {
    NumericMatrix ExtractionPoolsCoh(MAX_NUM_COHORTS, MAX_SOIL_LAYERS);
    std::fill(ExtractionPoolsCoh.begin(), ExtractionPoolsCoh.end(), 0.0);
    ExtractionPools[c] = ExtractionPoolsCoh;
  }
  
  NumericVector Stand = NumericVector::create(_["LAI"] = NA_REAL,
                                              _["LAIlive"] = NA_REAL, 
                                              _["LAIexpanded"] = NA_REAL, 
                                              _["LAIdead"] = NA_REAL);
  
  NumericVector LAI(MAX_NUM_COHORTS, NA_REAL);
  NumericVector LAIlive(MAX_NUM_COHORTS, NA_REAL);
  NumericVector FPAR(MAX_NUM_COHORTS, NA_REAL);
  NumericVector AbsorbedSWRFraction(MAX_NUM_COHORTS, NA_REAL);
  NumericVector ExtractionByPlant(MAX_NUM_COHORTS, NA_REAL);
  NumericVector Transpiration(MAX_NUM_COHORTS, NA_REAL);
  NumericVector GrossPhotosynthesis(MAX_NUM_COHORTS, NA_REAL);
  NumericVector PlantPsi(MAX_NUM_COHORTS, NA_REAL);
  NumericVector DDS(MAX_NUM_COHORTS, NA_REAL);
  NumericVector StemRWC(MAX_NUM_COHORTS, NA_REAL);
  NumericVector LeafRWC(MAX_NUM_COHORTS, NA_REAL);
  NumericVector LFMC(MAX_NUM_COHORTS, NA_REAL);
  NumericVector StemPLC(MAX_NUM_COHORTS, NA_REAL);
  NumericVector LeafPLC(MAX_NUM_COHORTS, NA_REAL);
  NumericVector PWB(MAX_NUM_COHORTS, NA_REAL);
  
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
List basicSPWBCommunicationOutput(List outputTransp) {
  
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
  
  DataFrame Soil = DataFrame::create(_["Psi"] = NumericVector(MAX_SOIL_LAYERS, 0.0),
                                     _["HerbTranspiration"] = NumericVector(MAX_SOIL_LAYERS, 0.0),
                                     _["HydraulicInput"] = NumericVector(MAX_SOIL_LAYERS, 0.0), //Water that entered into the layer across all time steps
                                     _["HydraulicOutput"] = NumericVector(MAX_SOIL_LAYERS, 0.0), //Water that left the layer across all time steps
                                     _["PlantExtraction"] = NumericVector(MAX_SOIL_LAYERS, 0.0));
  
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
List advancedTranspirationCommunicationOutput() {
  NumericMatrix SoilWaterExtract(MAX_NUM_COHORTS, MAX_SOIL_LAYERS);
  std::fill(SoilWaterExtract.begin(), SoilWaterExtract.end(), 0.0);
  
  List ExtractionPools(MAX_NUM_COHORTS);
  for(int c=0;c<MAX_NUM_COHORTS;c++) {
    NumericMatrix ExtractionPoolsCoh(MAX_NUM_COHORTS, MAX_SOIL_LAYERS);
    std::fill(ExtractionPoolsCoh.begin(), ExtractionPoolsCoh.end(), 0.0);
    ExtractionPools[c] = ExtractionPoolsCoh;
  }
  
  NumericMatrix soilLayerExtractInst(MAX_SOIL_LAYERS, MAX_NUM_TIMESTEPS);
  std::fill(soilLayerExtractInst.begin(), soilLayerExtractInst.end(), 0.0);

  NumericMatrix minPsiRhizo(MAX_NUM_COHORTS, MAX_SOIL_LAYERS);
  std::fill(minPsiRhizo.begin(), minPsiRhizo.end(), NA_REAL);

  NumericVector Stand = NumericVector::create(_["LAI"] = NA_REAL,
                                              _["LAIlive"] = NA_REAL, 
                                              _["LAIexpanded"] = NA_REAL, 
                                              _["LAIdead"] = NA_REAL);
  
  // ARRANGE OUTPUT
  NumericVector solarHour(MAX_NUM_TIMESTEPS, NA_REAL);
  NumericVector Tatm(MAX_NUM_TIMESTEPS, NA_REAL), Tcan(MAX_NUM_TIMESTEPS, NA_REAL);
  NumericVector Hcansoil(MAX_NUM_TIMESTEPS, NA_REAL);
  NumericVector Ebalsoil(MAX_NUM_TIMESTEPS, NA_REAL);
  NumericVector LEVsoil(MAX_NUM_TIMESTEPS, NA_REAL);
  NumericVector LEFsnow(MAX_NUM_TIMESTEPS, NA_REAL);
  NumericVector abs_SWR_soil(MAX_NUM_TIMESTEPS, NA_REAL);
  NumericVector net_LWR_soil(MAX_NUM_TIMESTEPS, NA_REAL);
  NumericVector abs_SWR_can(MAX_NUM_TIMESTEPS, NA_REAL);
  NumericVector net_LWR_can(MAX_NUM_TIMESTEPS, NA_REAL);
  NumericVector LEVcan(MAX_NUM_TIMESTEPS, NA_REAL), Hcan_heat(MAX_NUM_TIMESTEPS, NA_REAL), Ebal(MAX_NUM_TIMESTEPS, NA_REAL);
  
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
  NumericMatrix Tcan_mat(MAX_NUM_TIMESTEPS, MAX_CANOPY_LAYERS);
  NumericMatrix VPcan_mat(MAX_NUM_TIMESTEPS, MAX_CANOPY_LAYERS);
  NumericMatrix Tsoil_mat(MAX_NUM_TIMESTEPS, MAX_SOIL_LAYERS);
  List EB = List::create(_["Temperature"]=Tinst, 
                         _["SoilTemperature"] = Tsoil_mat,
                         _["CanopyEnergyBalance"] = CEBinst, 
                         _["SoilEnergyBalance"] = SEBinst,
                         _["TemperatureLayers"] = Tcan_mat, 
                         _["VaporPressureLayers"] = VPcan_mat);
  
  NumericVector LAI(MAX_NUM_COHORTS,0.0);
  NumericVector LAIlive(MAX_NUM_COHORTS,0.0);
  NumericVector PARcohort(MAX_NUM_COHORTS,0.0);
  NumericVector SoilExtractCoh(MAX_NUM_COHORTS,0.0);
  NumericVector DDS(MAX_NUM_COHORTS, 0.0), LFMC(MAX_NUM_COHORTS, 0.0);
  NumericVector Eplant(MAX_NUM_COHORTS, 0.0), Anplant(MAX_NUM_COHORTS, 0.0), Agplant(MAX_NUM_COHORTS, 0.0);
  NumericVector minStemPsi(MAX_NUM_COHORTS, NA_REAL), minRootPsi(MAX_NUM_COHORTS,NA_REAL); //Minimum potentials experienced
  NumericVector minLeafPsi(MAX_NUM_COHORTS,NA_REAL), maxLeafPsi(MAX_NUM_COHORTS,NA_REAL); 
  NumericVector PLClm(MAX_NUM_COHORTS, NA_REAL), PLCsm(MAX_NUM_COHORTS, NA_REAL);
  NumericVector dEdPm(MAX_NUM_COHORTS, NA_REAL), PWB(MAX_NUM_COHORTS,0.0);
  NumericVector RWCsm(MAX_NUM_COHORTS, NA_REAL), RWClm(MAX_NUM_COHORTS, NA_REAL);
  
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

  NumericVector maxGSW_SL(MAX_NUM_COHORTS,NA_REAL), maxGSW_SH(MAX_NUM_COHORTS,NA_REAL); 
  NumericVector minGSW_SL(MAX_NUM_COHORTS,NA_REAL), minGSW_SH(MAX_NUM_COHORTS,NA_REAL); 
  NumericVector maxTemp_SL(MAX_NUM_COHORTS,NA_REAL), maxTemp_SH(MAX_NUM_COHORTS,NA_REAL); 
  NumericVector minTemp_SL(MAX_NUM_COHORTS,NA_REAL), minTemp_SH(MAX_NUM_COHORTS,NA_REAL); 
  NumericVector minLeafPsi_SL(MAX_NUM_COHORTS,NA_REAL), maxLeafPsi_SL(MAX_NUM_COHORTS,NA_REAL); 
  NumericVector minLeafPsi_SH(MAX_NUM_COHORTS,NA_REAL), maxLeafPsi_SH(MAX_NUM_COHORTS,NA_REAL);
  
  DataFrame Sunlit = DataFrame::create(
    _["LeafPsiMin"] = minLeafPsi_SL, 
    _["LeafPsiMax"] = maxLeafPsi_SL, 
    _["GSWMin"] = minGSW_SL,
    _["GSWMax"] = maxGSW_SL,
    _["TempMin"] = minTemp_SL,
    _["TempMax"] = maxTemp_SL  
  );
  DataFrame Shade = DataFrame::create(
    _["LeafPsiMin"] = minLeafPsi_SH, 
    _["LeafPsiMax"] = maxLeafPsi_SH, 
    _["GSWMin"] = minGSW_SH,
    _["GSWMax"] = maxGSW_SH,
    _["TempMin"] = minTemp_SH,
    _["TempMax"] = maxTemp_SH  
  );

  
  NumericMatrix Einst(MAX_NUM_COHORTS, MAX_NUM_TIMESTEPS);
  NumericMatrix Aninst(MAX_NUM_COHORTS, MAX_NUM_TIMESTEPS), Aginst(MAX_NUM_COHORTS, MAX_NUM_TIMESTEPS);
  NumericMatrix dEdPInst(MAX_NUM_COHORTS, MAX_NUM_TIMESTEPS);
  NumericMatrix LeafPsiInst(MAX_NUM_COHORTS, MAX_NUM_TIMESTEPS), StemPsiInst(MAX_NUM_COHORTS, MAX_NUM_TIMESTEPS);
  NumericMatrix RootPsiInst(MAX_NUM_COHORTS, MAX_NUM_TIMESTEPS);
  NumericMatrix LeafSympPsiInst(MAX_NUM_COHORTS, MAX_NUM_TIMESTEPS), StemSympPsiInst(MAX_NUM_COHORTS, MAX_NUM_TIMESTEPS);
  NumericMatrix StemPLC(MAX_NUM_COHORTS, MAX_NUM_TIMESTEPS), LeafPLC(MAX_NUM_COHORTS, MAX_NUM_TIMESTEPS);
  NumericMatrix LeafRWCInst(MAX_NUM_COHORTS, MAX_NUM_TIMESTEPS), StemRWCInst(MAX_NUM_COHORTS, MAX_NUM_TIMESTEPS);
  NumericMatrix LeafSympRWCInst(MAX_NUM_COHORTS, MAX_NUM_TIMESTEPS), StemSympRWCInst(MAX_NUM_COHORTS, MAX_NUM_TIMESTEPS);
  NumericMatrix PWBinst(MAX_NUM_COHORTS, MAX_NUM_TIMESTEPS);

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
  
  NumericMatrix LAI_SL(MAX_NUM_COHORTS, MAX_NUM_TIMESTEPS);
  NumericMatrix LAI_SH(MAX_NUM_COHORTS, MAX_NUM_TIMESTEPS);
  NumericMatrix Vmax298_SL(MAX_NUM_COHORTS, MAX_NUM_TIMESTEPS);
  NumericMatrix Vmax298_SH(MAX_NUM_COHORTS, MAX_NUM_TIMESTEPS);
  NumericMatrix Jmax298_SL(MAX_NUM_COHORTS, MAX_NUM_TIMESTEPS);
  NumericMatrix Jmax298_SH(MAX_NUM_COHORTS, MAX_NUM_TIMESTEPS);
  NumericMatrix SWR_SL(MAX_NUM_COHORTS, MAX_NUM_TIMESTEPS);
  NumericMatrix SWR_SH(MAX_NUM_COHORTS, MAX_NUM_TIMESTEPS);
  NumericMatrix PAR_SL(MAX_NUM_COHORTS, MAX_NUM_TIMESTEPS);
  NumericMatrix PAR_SH(MAX_NUM_COHORTS, MAX_NUM_TIMESTEPS);
  NumericMatrix LWR_SL(MAX_NUM_COHORTS, MAX_NUM_TIMESTEPS);
  NumericMatrix LWR_SH(MAX_NUM_COHORTS, MAX_NUM_TIMESTEPS);
  NumericMatrix An_SL(MAX_NUM_COHORTS, MAX_NUM_TIMESTEPS), Ag_SL(MAX_NUM_COHORTS, MAX_NUM_TIMESTEPS);
  NumericMatrix An_SH(MAX_NUM_COHORTS, MAX_NUM_TIMESTEPS), Ag_SH(MAX_NUM_COHORTS, MAX_NUM_TIMESTEPS);
  NumericMatrix Ci_SL(MAX_NUM_COHORTS, MAX_NUM_TIMESTEPS);
  NumericMatrix Ci_SH(MAX_NUM_COHORTS, MAX_NUM_TIMESTEPS);
  NumericMatrix E_SL(MAX_NUM_COHORTS, MAX_NUM_TIMESTEPS);
  NumericMatrix E_SH(MAX_NUM_COHORTS, MAX_NUM_TIMESTEPS);
  NumericMatrix GSW_SH(MAX_NUM_COHORTS, MAX_NUM_TIMESTEPS);
  NumericMatrix GSW_SL(MAX_NUM_COHORTS, MAX_NUM_TIMESTEPS);
  NumericMatrix VPD_SH(MAX_NUM_COHORTS, MAX_NUM_TIMESTEPS);
  NumericMatrix VPD_SL(MAX_NUM_COHORTS, MAX_NUM_TIMESTEPS);
  NumericMatrix Temp_SH(MAX_NUM_COHORTS, MAX_NUM_TIMESTEPS);
  NumericMatrix Temp_SL(MAX_NUM_COHORTS, MAX_NUM_TIMESTEPS);
  NumericMatrix Psi_SL(MAX_NUM_COHORTS, MAX_NUM_TIMESTEPS);
  NumericMatrix Psi_SH(MAX_NUM_COHORTS, MAX_NUM_TIMESTEPS);

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
  
  List lwrExtinctionList(MAX_NUM_TIMESTEPS);
  for(int n=0;n<MAX_NUM_TIMESTEPS;n++) {
    lwrExtinctionList[n] = communicationLongWaveRadiation();
  }
  List supply(MAX_NUM_COHORTS); 
  for(int c=0;c<MAX_NUM_COHORTS;c++) {
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
                        _["CanopyTurbulence"] = communicationCanopyTurbulence()); //To be replaced
  
  List outPhotoSunlit(MAX_NUM_COHORTS);
  List outPhotoShade(MAX_NUM_COHORTS);
  List outPMSunlit(MAX_NUM_COHORTS);
  List outPMShade(MAX_NUM_COHORTS);
  l.push_back(supply, "SupplyFunctions");
  l.push_back(outPhotoSunlit, "PhotoSunlitFunctions");
  l.push_back(outPhotoShade, "PhotoShadeFunctions");
  l.push_back(outPMSunlit, "PMSunlitFunctions");
  l.push_back(outPMShade, "PMShadeFunctions");
  
  l.attr("class") = CharacterVector::create("pwb_day","list");
  return(l);
}
List advancedSPWBCommunicationOutput(List outputTransp) {

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
  
  DataFrame Soil = DataFrame::create(_["Psi"] = NumericVector(MAX_SOIL_LAYERS, 0.0),
                                     _["HerbTranspiration"] = NumericVector(MAX_SOIL_LAYERS, 0.0),
                                     _["HydraulicInput"] = NumericVector(MAX_SOIL_LAYERS, 0.0), //Water that entered into the layer across all time steps
                                     _["HydraulicOutput"] = NumericVector(MAX_SOIL_LAYERS, 0.0), //Water that left the layer across all time steps
                                     _["PlantExtraction"] = NumericVector(MAX_SOIL_LAYERS, 0.0));
  
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
  l.attr("class") = CharacterVector::create("spwb_day","list");
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
      ExtractionPools[c] = copyNumericMatrix(ExtractionPoolsCohComm, numCohorts, nlayers); // this is final extraction of each cohort from each layer
    }
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
      ExtractionPools[c] = copyNumericMatrix(ExtractionPoolsCohComm, numCohorts, nlayers); // this is final extraction of each cohort from each layer
    }
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
  
  NumericVector Stand = clone(as<NumericVector>(boc["Stand"]));
  
  DataFrame Soil = copyDataFrame(as<DataFrame>(boc["Soil"]), nlayers);
  
  DataFrame PlantsComm = Rcpp::as<Rcpp::DataFrame>(boc["Plants"]);
  DataFrame Plants = copyDataFrame(PlantsComm, numCohorts);
  Plants.attr("row.names") = cohorts.attr("row.names");
  
  List l = List::create(_["cohorts"] = clone(cohorts),
                        _["topography"] = topo,
                        _["weather"] = meteovec_bas,
                        _["WaterBalance"] = WaterBalance, 
                        _["Soil"] = Soil,
                        _["Stand"] = Stand,
                        _["Plants"] = Plants);
  if(control["fireHazardResults"]) {
    l.push_back(clone(as<NumericVector>(boc["FireHazard"])), "FireHazard"); 
  }
  l.attr("class") = CharacterVector::create("spwb_day","list");
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
  
  NumericVector Stand = clone(as<NumericVector>(aoc["Stand"]));
  
  DataFrame Soil = copyDataFrame(as<DataFrame>(aoc["Soil"]), nlayers);
  
  DataFrame Plants = copyDataFrame(as<DataFrame>(aoc["Plants"]), numCohorts);
  Plants.attr("row.names") = cohorts.attr("row.names");

  DataFrame Sunlit = copyDataFrame(as<DataFrame>(aoc["SunlitLeaves"]), numCohorts);
  Sunlit.attr("row.names") = above.attr("row.names");
  DataFrame Shade = copyDataFrame(as<DataFrame>(aoc["ShadeLeaves"]), numCohorts);
  Shade.attr("row.names") = above.attr("row.names");
  
  
  NumericMatrix minPsiRhizoComm = aoc["RhizoPsi"];
  NumericMatrix minPsiRhizo = copyNumericMatrix(minPsiRhizoComm, numCohorts, nlayers);
  minPsiRhizo.attr("dimnames") = List::create(above.attr("row.names"), seq(1,nlayers));
  
  NumericMatrix soilLayerExtractInstComm = aoc["ExtractionInst"];
  NumericMatrix soilLayerExtractInst = copyNumericMatrix(soilLayerExtractInstComm, nlayers, ntimesteps);
  soilLayerExtractInst.attr("dimnames") = List::create(seq(1,nlayers), seq(1,ntimesteps));
  
  List PlantsInst = copyPlantsInstOutput(as<List>(aoc["PlantsInst"]), x);
  List SunlitInst = copyLeavesInstOutput(as<List>(aoc["SunlitLeavesInst"]), x);
  List ShadeInst = copyLeavesInstOutput(as<List>(aoc["ShadeLeavesInst"]), x);
  
  List lwrExtinctionListComm = aoc["LWRExtinction"];
  List lwrExtinctionList(ntimesteps);
  for(int n=0;n<ntimesteps;n++) {
    lwrExtinctionList[n] = copyCommunicationLongWaveRadiation(as<List>(lwrExtinctionListComm[n]), ncanlayers);
  }
  

  List l = List::create(_["cohorts"] = clone(cohorts),
                        _["topography"] = topo,
                        _["weather"] = meteovec_adv,
                        _["WaterBalance"] = WaterBalance, 
                        _["EnergyBalance"] = EnergyBalance,
                        _["Soil"] = Soil, 
                        _["Stand"] = Stand, 
                        _["Plants"] = Plants,
                        _["RhizoPsi"] = minPsiRhizo,
                        _["SunlitLeaves"] = Sunlit,
                        _["ShadeLeaves"] = Shade,
                        _["ExtractionInst"] = soilLayerExtractInst,
                        _["PlantsInst"] = PlantsInst,
                        _["RadiationInputInst"] = clone(as<List>(aoc["RadiationInputInst"])),
                        _["SunlitLeavesInst"] = SunlitInst,
                        _["ShadeLeavesInst"] = ShadeInst,
                        _["LightExtinction"] = clone(as<List>(aoc["LightExtinction"])),
                        _["LWRExtinction"] = lwrExtinctionListComm,
                        _["CanopyTurbulence"] = copyDataFrame(as<DataFrame>(aoc["CanopyTurbulence"]), ncanlayers));
  l.attr("class") = CharacterVector::create("spwb_day","list");
  return(l);
}

List basicGROWTHOutput(List x, List outputTransp) {
  List control = x["control"];
  DataFrame cohorts = Rcpp::as<Rcpp::DataFrame>(x["cohorts"]);
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  int numCohorts = cohorts.nrow();
  List spwbOut = basicSPWBCommunicationOutput(outputTransp);
  
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
  labileCarbonBalance.attr("row.names") = above.attr("row.names");
  
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
  growthMortality.attr("row.names") = above.attr("row.names");
  
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
  plantBiomassBalance.attr("row.names") = above.attr("row.names");
  
  List l = List::create(_["cohorts"] = spwbOut["cohorts"],
                        _["topography"] = spwbOut["topography"],
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
  if(control["fireHazardResults"]) l.push_back(spwbOut["FireHazard"], "FireHazard");
  l.attr("class") = CharacterVector::create("growth_day","list");
  return(l);
}
List advancedGROWTHOutput(List x, List outputTransp) {
  List control = x["control"];
  bool subdailyCarbonBalance = control["subdailyCarbonBalance"];
  DataFrame cohorts = Rcpp::as<Rcpp::DataFrame>(x["cohorts"]);
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  int numCohorts = cohorts.nrow();
  List spwbOut = advancedSPWBCommunicationOutput(outputTransp);
  
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
  labileCarbonBalance.attr("row.names") = above.attr("row.names");
  
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
  growthMortality.attr("row.names") = above.attr("row.names");
  
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
  plantBiomassBalance.attr("row.names") = above.attr("row.names");
  
  
  
  List l = List::create(_["cohorts"] = spwbOut["cohorts"],
                        _["topography"] = spwbOut["topography"],
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
  if(subdailyCarbonBalance) {
    List PlantsInst = spwbOut["PlantsInst"];
    NumericMatrix AgStep  =  Rcpp::as<Rcpp::NumericMatrix>(PlantsInst["Ag"]);
    int numSteps = AgStep.ncol();

    //Subdaily output matrices (Sperry/Sureau)
    NumericMatrix LabileCarbonBalanceInst(numCohorts, numSteps);  
    NumericMatrix GrossPhotosynthesisInst(numCohorts, numSteps);  
    NumericMatrix MaintenanceRespirationInst(numCohorts, numSteps);  
    NumericMatrix GrowthCostsInst(numCohorts, numSteps);  
    NumericMatrix RootExudationInst(numCohorts, numSteps);  
    NumericMatrix PlantSugarTransportInst(numCohorts, numSteps);
    NumericMatrix PlantSugarLeafInst(numCohorts, numSteps), PlantStarchLeafInst(numCohorts, numSteps);
    NumericMatrix PlantSugarSapwoodInst(numCohorts, numSteps), PlantStarchSapwoodInst(numCohorts, numSteps);
    
    GrossPhotosynthesisInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,numSteps));
    MaintenanceRespirationInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,numSteps));
    LabileCarbonBalanceInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,numSteps));
    GrowthCostsInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,numSteps));
    PlantSugarLeafInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,numSteps));
    PlantStarchLeafInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,numSteps));
    PlantSugarSapwoodInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,numSteps));
    PlantStarchSapwoodInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,numSteps));
    PlantSugarTransportInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,numSteps));
    RootExudationInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,numSteps));
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
    if(subdailyCarbonBalance) l.push_back(labileCBInst,"LabileCarbonBalanceInst");
  }
  l.push_back(spwbOut["LightExtinction"], "LightExtinction");
  l.push_back(spwbOut["CanopyTurbulence"], "CanopyTurbulence");
  if(control["fireHazardResults"]) l.push_back(spwbOut["FireHazard"], "FireHazard");
  l.attr("class") = CharacterVector::create("growth_day","list");
  return(l);
}

DataFrame communicationCarbonCompartments() {
  DataFrame df = DataFrame::create(
    _["LeafStorageVolume"] = NumericVector(MAX_NUM_COHORTS, NA_REAL),
    _["SapwoodStorageVolume"] = NumericVector(MAX_NUM_COHORTS, NA_REAL),
    _["LeafStarchMaximumConcentration"] = NumericVector(MAX_NUM_COHORTS, NA_REAL),
    _["SapwoodStarchMaximumConcentration"] = NumericVector(MAX_NUM_COHORTS, NA_REAL),
    _["LeafStarchCapacity"] = NumericVector(MAX_NUM_COHORTS, NA_REAL),
    _["SapwoodStarchCapacity"] = NumericVector(MAX_NUM_COHORTS, NA_REAL),
    _["LeafStructuralBiomass"] = NumericVector(MAX_NUM_COHORTS, NA_REAL),
    _["SapwoodStructuralBiomass"] = NumericVector(MAX_NUM_COHORTS, NA_REAL),
    _["SapwoodLivingStructuralBiomass"] = NumericVector(MAX_NUM_COHORTS, NA_REAL),
    _["FineRootBiomass"] = NumericVector(MAX_NUM_COHORTS, NA_REAL),
    _["StructuralBiomass"] = NumericVector(MAX_NUM_COHORTS, NA_REAL),
    _["LabileBiomass"] = NumericVector(MAX_NUM_COHORTS, NA_REAL),
    _["TotalLivingBiomass"] = NumericVector(MAX_NUM_COHORTS, NA_REAL),
    _["TotalBiomass"] = NumericVector(MAX_NUM_COHORTS, NA_REAL)
  );
  return(df);
}

List communicationInitialFinalCarbonCompartments() {
  DataFrame ccFin_g_ind = communicationCarbonCompartments();
  DataFrame ccIni_g_ind = clone(ccFin_g_ind);
  List l = List::create(_["ccIni_g_ind"] = ccIni_g_ind,
                        _["ccFin_g_ind"] = ccFin_g_ind);
  return(l);
}


// [[Rcpp::export(".instanceCommunicationStructures")]]
List instanceCommunicationStructures() {
  List basicTranspirationOutput = basicTranspirationCommunicationOutput();
  List basicSPWBOutput = basicSPWBCommunicationOutput(basicTranspirationOutput);
  List advancedTranspirationOutput = advancedTranspirationCommunicationOutput();
  List advancedSPWBOutput = advancedSPWBCommunicationOutput(advancedTranspirationOutput);
  List ic = List::create(_["basicTranspirationOutput"] = basicTranspirationOutput,
                         _["advancedTranspirationOutput"] = advancedTranspirationOutput,
                         _["basicSPWBOutput"] = basicSPWBOutput,
                         _["advancedSPWBOutput"] = advancedSPWBOutput);
  // ic.push_back(advancedTranspirationCommunicationOutput(), "advancedTranspirationOutput"); 
  // List communicationOutputTransp = ic["basicTranspirationOutput"];
  // if(transpirationMode=="Granier") {
  //   if(!ic.containsElementNamed("modelOutput")) {
  //     if(model=="spwb") ic.push_back(basicSPWBOutput(x, communicationOutputTransp), "modelOutput"); 
  //     else ic.push_back(basicGROWTHOutput(x, communicationOutputTransp), "modelOutput"); 
  //   } 
  // } else {
  //   DataFrame paramsCanopydf = as<DataFrame>(x["canopy"]);
  //   List outputTransp = ic["transpirationOutput"];
  //   if(!ic.containsElementNamed("modelOutput")) {
  //     if(model=="spwb") ic.push_back(advancedSPWBOutput(x, outputTransp), "modelOutput"); 
  //     else ic.push_back(advancedGROWTHOutput(x, outputTransp), "modelOutput"); 
  //   } 
  // }
  // if(model=="growth") {
  //   if(!ic.containsElementNamed("initialFinalCC")) ic.push_back(communicationInitialFinalCarbonCompartments(), "initialFinalCC");
  // }
  return(ic);
}

// [[Rcpp::export(".clearCommunicationStructures")]]
void clearCommunicationStructures(List x) {
  List control = x["control"];
  x["internalCommunication"] = List::create();
}