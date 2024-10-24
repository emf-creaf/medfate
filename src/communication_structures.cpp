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
      std::fill(ExtractionPoolsCoh.begin(), ExtractionPoolsCoh.end(), 0.0);
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
List advancedTranspirationOutput(List x) {
  List control = x["control"];
  int ntimesteps = control["ndailysteps"];
  String transpirationMode = control["transpirationMode"];
  String rhizosphereOverlap = control["rhizosphereOverlap"];
  bool plantWaterPools = (rhizosphereOverlap!="total");
  DataFrame cohorts = Rcpp::as<Rcpp::DataFrame>(x["cohorts"]);
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  DataFrame soil = Rcpp::as<Rcpp::DataFrame>(x["soil"]);
  int nlayers = Rcpp::as<Rcpp::NumericVector>(soil["widths"]).size();
  int numCohorts = above.nrow();
  DataFrame canopyParams = Rcpp::as<Rcpp::DataFrame>(x["canopy"]);
  NumericVector Tair = canopyParams["Tair"];
  int ncanlayers = Tair.size(); //Number of canopy layers
  
  NumericMatrix SoilWaterExtract(numCohorts, nlayers);
  std::fill(SoilWaterExtract.begin(), SoilWaterExtract.end(), 0.0);
  
  List ExtractionPools(numCohorts);
  if(plantWaterPools) {
    for(int c=0;c<numCohorts;c++) {
      NumericMatrix ExtractionPoolsCoh(numCohorts, nlayers);
      std::fill(ExtractionPoolsCoh.begin(), ExtractionPoolsCoh.end(), 0.0);
      ExtractionPools[c] = ExtractionPoolsCoh;
    }
  }
  
  NumericMatrix soilLayerExtractInst(nlayers, ntimesteps);
  std::fill(soilLayerExtractInst.begin(), soilLayerExtractInst.end(), 0.0);
  soilLayerExtractInst.attr("dimnames") = List::create(seq(1,nlayers), seq(1,ntimesteps));
  
  NumericMatrix minPsiRhizo(numCohorts, nlayers);
  if(numCohorts>0) std::fill(minPsiRhizo.begin(), minPsiRhizo.end(), NA_REAL);
  minPsiRhizo.attr("dimnames") = List::create(above.attr("row.names"), seq(1,nlayers));
  
  
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
  Tcan_mat.attr("dimnames") = List::create(seq(1,ntimesteps), seq(1,ncanlayers));
  VPcan_mat.attr("dimnames") = List::create(seq(1,ntimesteps), seq(1,ncanlayers));
  NumericMatrix Tsoil_mat(ntimesteps, nlayers);
  Tsoil_mat.attr("dimnames") = List::create(seq(1,ntimesteps), seq(1,nlayers));
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
  Plants.attr("row.names") = above.attr("row.names");
  
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
  DataFrame Shade = DataFrame::create(
    _["LeafPsiMin"] = minLeafPsi_SH, 
    _["LeafPsiMax"] = maxLeafPsi_SH, 
    _["GSWMin"] = minGSW_SH,
    _["GSWMax"] = maxGSW_SH,
    _["TempMin"] = minTemp_SH,
    _["TempMax"] = maxTemp_SH  
  );
  Sunlit.attr("row.names") = above.attr("row.names");
  Shade.attr("row.names") = above.attr("row.names");
  
  
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
  Einst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  dEdPInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  LeafPsiInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  StemPsiInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  LeafSympPsiInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  StemSympPsiInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  LeafSympRWCInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  StemSympRWCInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  RootPsiInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  Aginst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  Aninst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  StemPLC.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  LeafPLC.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  LeafRWCInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  StemRWCInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
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
  LAI_SH.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  LAI_SL.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  Vmax298_SH.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  Vmax298_SL.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  Jmax298_SH.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  Jmax298_SL.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  E_SH.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  E_SL.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  Psi_SH.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  Psi_SL.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  An_SH.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  An_SL.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  Ag_SH.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  Ag_SL.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  SWR_SH.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  SWR_SL.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  PAR_SH.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  PAR_SL.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  LWR_SH.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  LWR_SL.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  Ci_SH.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  Ci_SL.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  GSW_SH.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  GSW_SL.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  Temp_SH.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  Temp_SL.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  VPD_SH.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  VPD_SL.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  
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
  supply.attr("names") = above.attr("row.names");
  
  List l = List::create(_["cohorts"] = clone(cohorts),
                        _["EnergyBalance"] = EB,
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
  outPhotoSunlit.attr("names") = above.attr("row.names");
  outPhotoShade.attr("names") = above.attr("row.names");
  outPMSunlit.attr("names") = above.attr("row.names");
  outPMShade.attr("names") = above.attr("row.names");
  l.push_back(supply, "SupplyFunctions");
  l.push_back(outPhotoSunlit, "PhotoSunlitFunctions");
  l.push_back(outPhotoShade, "PhotoShadeFunctions");
  l.push_back(outPMSunlit, "PMSunlitFunctions");
  l.push_back(outPMShade, "PMShadeFunctions");
  
  l.attr("class") = CharacterVector::create("pwb_day","list");
  return(l);
}

List basicSPWBOutput(List x, List outputTransp) {
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
  
  List l = List::create(_["cohorts"] = clone(cohorts),
                        _["topography"] = topo,
                        _["weather"] = meteovec_bas,
                        _["WaterBalance"] = WaterBalance, 
                        _["Soil"] = Soil,
                        _["Stand"] = Stand,
                        _["Plants"] = outputTransp["Plants"]);
  if(control["fireHazardResults"]) l.push_back(List::create(), "FireHazard");
  l.attr("class") = CharacterVector::create("spwb_day","list");
  return(l);
}
List advancedSPWBOutput(List x, List outputTransp) {
  List control = x["control"];
  DataFrame cohorts = Rcpp::as<Rcpp::DataFrame>(x["cohorts"]);
  DataFrame soil = Rcpp::as<Rcpp::DataFrame>(x["soil"]);
  int nlayers = Rcpp::as<Rcpp::NumericVector>(soil["widths"]).size();
  
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
  
  List l = List::create(_["cohorts"] = clone(cohorts),
                        _["topography"] = topo,
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

List basicGROWTHOutput(List x, List outputTransp) {
  List control = x["control"];
  DataFrame cohorts = Rcpp::as<Rcpp::DataFrame>(x["cohorts"]);
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  int numCohorts = cohorts.nrow();
  List spwbOut = basicSPWBOutput(x, outputTransp);
  
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
  List spwbOut = advancedSPWBOutput(x, outputTransp);
  
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

DataFrame communicationCarbonCompartments(DataFrame above) {
  int numCohorts = above.nrow();
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
  df.attr("row.names") = above.attr("row.names");
  return(df);
}

List communicationInitialFinalCarbonCompartments(DataFrame above) {
  DataFrame ccFin_g_ind = communicationCarbonCompartments(above);
  DataFrame ccIni_g_ind = clone(ccFin_g_ind);
  List l = List::create(_["ccIni_g_ind"] = ccIni_g_ind,
                        _["ccFin_g_ind"] = ccFin_g_ind);
  return(l);
}


// [[Rcpp::export(".addCommunicationStructures")]]
void addCommunicationStructures(List x) {
  String model = "spwb";
  if(x.inherits("growthInput")) model = "growth";
  List control = x["control"];
  String transpirationMode = control["transpirationMode"];
  List ic = as<List>(x["internalCommunication"]);
  if(transpirationMode=="Granier") {
    if(!ic.containsElementNamed("transpirationOutput")) ic.push_back(basicTranspirationOutput(x), "transpirationOutput"); 
    List outputTransp = ic["transpirationOutput"];
    if(!ic.containsElementNamed("modelOutput")) {
      if(model=="spwb") ic.push_back(basicSPWBOutput(x, outputTransp), "modelOutput"); 
      else ic.push_back(basicGROWTHOutput(x, outputTransp), "modelOutput"); 
    } 
  } else {
    DataFrame paramsCanopydf = as<DataFrame>(x["canopy"]);
    if(!ic.containsElementNamed("transpirationOutput")) ic.push_back(advancedTranspirationOutput(x), "transpirationOutput"); 
    List outputTransp = ic["transpirationOutput"];
    if(!ic.containsElementNamed("modelOutput")) {
      if(model=="spwb") ic.push_back(advancedSPWBOutput(x, outputTransp), "modelOutput"); 
      else ic.push_back(advancedGROWTHOutput(x, outputTransp), "modelOutput"); 
    } 
  }
  if(model=="growth") {
    DataFrame above = as<DataFrame>(x["above"]);
    if(!ic.containsElementNamed("initialFinalCC")) ic.push_back(communicationInitialFinalCarbonCompartments(above), "initialFinalCC");
  }
  x["internalCommunication"] = ic;
}

// [[Rcpp::export(".clearCommunicationStructures")]]
void clearCommunicationStructures(List x) {
  List control = x["control"];
  x["internalCommunication"] = List::create();
}