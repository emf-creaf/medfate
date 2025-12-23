// [[Rcpp::interfaces(r,cpp)]]
#define STRICT_R_HEADERS
#include <Rcpp.h>
#include <numeric>
#include "biophysicsutils.h"
#include "carbon.h"
#include "communication_structures.h"
#include "forestutils.h"
#include "growth_day.h"
#include "hydrology.h"
#include "modelInput.h"
#include "phenology.h"
#include "spwb.h"
#include "soil.h"
#include "tissuemoisture.h"
#include <meteoland.h>
using namespace Rcpp;


void checkgrowthInput(List x, String transpirationMode, String soilFunctions) {
  
  DataFrame soil = Rcpp::as<Rcpp::DataFrame>(x["soil"]);
  if(!x.containsElementNamed("above")) stop("above missing in growthInput");
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  if(!above.containsElementNamed("LAI_live")) stop("LAI_live missing in growthInput$above");
  if(!above.containsElementNamed("LAI_expanded")) stop("LAI_expanded missing in growthInput$above");
  if(!above.containsElementNamed("LAI_dead")) stop("LAI_dead missing in growthInput$above");
  if(!above.containsElementNamed("SA")) stop("SA missing in growthInput$above");
  if(!above.containsElementNamed("CR")) stop("CR missing in growthInput$above");
  if(!above.containsElementNamed("H")) stop("H missing in growthInput$above");
  if(!above.containsElementNamed("N")) stop("N missing in growthInput$above");
  if(!above.containsElementNamed("DBH")) stop("DBH missing in growthInput$above");
  
  if(!x.containsElementNamed("below")) stop("below missing in growthInput");
  DataFrame below = Rcpp::as<Rcpp::DataFrame>(x["below"]);
  if(!below.containsElementNamed("Z50")) stop("Z50 missing in growthInput$below");
  if(!below.containsElementNamed("Z95")) stop("Z95 missing in growthInput$below");
  if(!x.containsElementNamed("belowLayers")) stop("belowLayers missing in growthInput");
  List belowLayers = Rcpp::as<Rcpp::List>(x["belowLayers"]);
  if(!belowLayers.containsElementNamed("V")) stop("V missing in growthInput$belowLayers");
  if((transpirationMode=="Sperry") || (transpirationMode=="Sureau")) {
    if(!belowLayers.containsElementNamed("VGrhizo_kmax")) stop("VGrhizo_kmax missing in growthInput$below");
    if(!belowLayers.containsElementNamed("VCroot_kmax")) stop("VCroot_kmax missing in growthInput$below");
  }  
  
  if(!x.containsElementNamed("paramsPhenology")) stop("paramsPhenology missing in growthInput");
  DataFrame paramsPhenology = Rcpp::as<Rcpp::DataFrame>(x["paramsPhenology"]);
  if(!paramsPhenology.containsElementNamed("Sgdd")) stop("Sgdd missing in paramsPhenology");
  
  if(!x.containsElementNamed("paramsInterception")) stop("paramsInterception missing in growthInput");
  DataFrame paramsInterception = Rcpp::as<Rcpp::DataFrame>(x["paramsInterception"]);
  if(!paramsInterception.containsElementNamed("kPAR")) stop("kPAR missing in growthInput$paramsInterception");
  if(!paramsInterception.containsElementNamed("g")) stop("g missing in growthInput$paramsInterception");
  
  if(!x.containsElementNamed("paramsGrowth")) stop("paramsGrowth missing in growthInput");
  DataFrame paramsGrowth = Rcpp::as<Rcpp::DataFrame>(x["paramsGrowth"]);
  if(!paramsGrowth.containsElementNamed("WoodC")) stop("WoodC missing in growthInput$paramsGrowth");
  if(!paramsGrowth.containsElementNamed("RGRleafmax")) stop("RGRleafmax missing in growthInput$paramsGrowth");
  if(!paramsGrowth.containsElementNamed("RGRsapwoodmax")) stop("RGRsapwoodmax missing in growthInput$paramsGrowth");
  if(!paramsGrowth.containsElementNamed("RGRfinerootmax")) stop("RGRfinerootmax missing in growthInput$paramsGrowth");
  if(!paramsGrowth.containsElementNamed("RERleaf")) stop("RERleaf missing in growthInput$paramsGrowth");
  if(!paramsGrowth.containsElementNamed("RERsapwood")) stop("RERsapwood missing in growthInput$paramsGrowth");
  if(!paramsGrowth.containsElementNamed("RERfineroot")) stop("RERfineroot missing in growthInput$paramsGrowth");
  
  if(!x.containsElementNamed("paramsAnatomy")) stop("paramsAnatomy missing in growthInput");
  DataFrame paramsAnatomy = Rcpp::as<Rcpp::DataFrame>(x["paramsAnatomy"]);
  if(!paramsAnatomy.containsElementNamed("SLA")) stop("SLA missing in paramsAnatomy$paramsGrowth");
  if(!paramsAnatomy.containsElementNamed("Al2As")) stop("Al2As missing in paramsAnatomy$paramsGrowth");
  if(!paramsAnatomy.containsElementNamed("WoodDensity")) stop("WoodDensity missing in paramsAnatomy$paramsGrowth");
  
  if(!x.containsElementNamed("paramsTranspiration")) stop("paramsTranspiration missing in growthInput");
  DataFrame paramsTranspiration = Rcpp::as<Rcpp::DataFrame>(x["paramsTranspiration"]);
  if(transpirationMode=="Granier") {
    if(!paramsTranspiration.containsElementNamed("Psi_Extract")) stop("Psi_Extract missing in growthInput$paramsTransp");
    if(!paramsTranspiration.containsElementNamed("WUE")) stop("WUE missing in growthInput$paramsTransp");
  } else if((transpirationMode=="Sperry") || (transpirationMode=="Sureau")) {
    if(!soil.containsElementNamed("VG_n")) stop("VG_n missing in soil");
    if(!soil.containsElementNamed("VG_alpha")) stop("VG_alpha missing in soil");
    
    if(!paramsTranspiration.containsElementNamed("VCstem_kmax")) stop("VCstem_kmax missing in growthInput");
    if(!paramsTranspiration.containsElementNamed("VCstem_c")) stop("VCstem_c missing in growthInput$paramsTransp");
    if(!paramsTranspiration.containsElementNamed("VCstem_d")) stop("VCstem_d missing in growthInput$paramsTransp");
    if(!paramsTranspiration.containsElementNamed("VCroot_c")) stop("VCroot_c missing in growthInput$paramsTransp");
    if(!paramsTranspiration.containsElementNamed("VCroot_d")) stop("VCroot_d missing in growthInput$paramsTransp");
    if(!paramsTranspiration.containsElementNamed("Gswmax")) stop("Gswmax missing in growthInput$paramsTransp");
    if(!paramsTranspiration.containsElementNamed("Vmax298")) stop("Vmax298 missing in growthInput$paramsTransp");
    if(!paramsTranspiration.containsElementNamed("Jmax298")) stop("Jmax298 missing in growthInput$paramsTransp");
  }
  if(!soil.containsElementNamed("W")) stop("W missing in soil");
  if(!soil.containsElementNamed("widths")) stop("widths missing in soil");
  if(!soil.containsElementNamed("macro")) stop("macro missing in soil");
  if(soilFunctions=="SX") {
    if(!soil.containsElementNamed("clay")) stop("clay missing in soil");
    if(!soil.containsElementNamed("sand")) stop("sand missing in soil");
  }
  if(soilFunctions=="VG") {
    if(!soil.containsElementNamed("VG_n")) stop("VG_n missing in soil");
    if(!soil.containsElementNamed("VG_alpha")) stop("VG_alpha missing in soil");
    if(!soil.containsElementNamed("VG_theta_res")) stop("VG_theta_res missing in soil");
    if(!soil.containsElementNamed("VG_theta_sat")) stop("VG_theta_sat missing in soil");
  }
}


NumericVector standLevelBiomassBalance(DataFrame biomassBalance) {
  return(NumericVector::create(
      _["StructuralBalance"] = sum(Rcpp::as<Rcpp::NumericVector>(biomassBalance["StructuralBiomassBalance"])),
      _["LabileBalance"] = sum(Rcpp::as<Rcpp::NumericVector>(biomassBalance["LabileBiomassBalance"])),
      _["PlantBalance"] = sum(Rcpp::as<Rcpp::NumericVector>(biomassBalance["PlantBiomassBalance"])),
      _["MortalityLoss"] = sum(Rcpp::as<Rcpp::NumericVector>(biomassBalance["MortalityBiomassLoss"])),
      _["CohortBalance"] = sum(Rcpp::as<Rcpp::NumericVector>(biomassBalance["CohortBiomassBalance"]))
  ));
}

// [[Rcpp::export(".defineGrowthDailyOutput")]]
List defineGrowthDailyOutput(double latitude, double elevation, double slope, double aspect, 
                             CharacterVector dateStrings, List x) {
  
  NumericVector topo = NumericVector::create(elevation, slope, aspect);
  topo.attr("names") = CharacterVector::create("elevation", "slope", "aspect");
  
  List growthInput = clone(x);
  
  //Check if input version is lower than current medfate version. If so, try to complete fields
  if(isLowerVersion(x)) growthInputVersionUpdate(x);
  if(!x.containsElementNamed("version")) {
    x.push_back(medfateVersionString(), "version");
  }
  
  
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  int numCohorts = above.nrow();
  int numDays = dateStrings.size();
  List control = x["control"];
  DataFrame soil = Rcpp::as<Rcpp::DataFrame>(x["soil"]);
  String transpirationMode = control["transpirationMode"];
  
  DataFrame DWB = defineWaterBalanceDailyOutput(dateStrings, transpirationMode);
  List Soil = defineSoilDailyOutput(dateStrings, soil, true);
  DataFrame Snow = defineSnowDailyOutput(dateStrings);
  List sunlitDO = defineSunlitShadeLeavesDailyOutput(dateStrings, above);
  List shadeDO = defineSunlitShadeLeavesDailyOutput(dateStrings, above);
  List plantDWOL = definePlantWaterDailyOutput(dateStrings, above, soil, control);
  DataFrame Stand = defineStandDailyOutput(dateStrings);
  
  //Detailed subday results
  List subdailyRes(numDays);
  subdailyRes.attr("names") = dateStrings ;
  
  //Plant carbon output variables
  NumericMatrix LabileCarbonBalance(numDays, numCohorts);
  NumericMatrix MaintenanceRespiration(numDays, numCohorts);
  NumericMatrix GrowthCosts(numDays, numCohorts);
  NumericMatrix PlantSugarLeaf(numDays, numCohorts);
  NumericMatrix PlantStarchLeaf(numDays, numCohorts);
  NumericMatrix PlantSugarSapwood(numDays, numCohorts);
  NumericMatrix PlantStarchSapwood(numDays, numCohorts);
  NumericMatrix PlantSugarTransport(numDays, numCohorts);
  NumericMatrix SapwoodBiomass(numDays, numCohorts);
  NumericMatrix LeafBiomass(numDays, numCohorts);
  NumericMatrix SapwoodArea(numDays, numCohorts);
  NumericMatrix LeafArea(numDays, numCohorts);
  NumericMatrix FineRootArea(numDays, numCohorts);
  NumericMatrix FineRootBiomass(numDays, numCohorts);
  NumericMatrix HuberValue(numDays, numCohorts);
  NumericMatrix RootAreaLeafArea(numDays, numCohorts);
  NumericMatrix DBH(numDays, numCohorts);
  NumericMatrix Height(numDays, numCohorts);
  NumericMatrix LabileBiomass(numDays, numCohorts);
  NumericMatrix TotalBiomass(numDays, numCohorts);
  NumericMatrix SAgrowth(numDays, numCohorts), LAgrowth(numDays, numCohorts), FRAgrowth(numDays, numCohorts);
  NumericMatrix starvationRate(numDays, numCohorts), dessicationRate(numDays, numCohorts), mortalityRate(numDays, numCohorts);
  NumericMatrix GrossPhotosynthesis(numDays, numCohorts);
  NumericMatrix PlantLAIexpanded(numDays, numCohorts), PlantLAIdead(numDays, numCohorts), PlantLAIlive(numDays, numCohorts);
  NumericMatrix RootExudation(numDays, numCohorts);
  NumericMatrix StructuralBiomassBalance(numDays, numCohorts);
  NumericMatrix LabileBiomassBalance(numDays, numCohorts);
  NumericMatrix PlantBiomassBalance(numDays, numCohorts);
  NumericMatrix MortalityBiomassLoss(numDays, numCohorts);
  NumericMatrix CohortBiomassBalance(numDays, numCohorts);
  NumericMatrix StandBiomassBalance(numDays, 5);
  NumericMatrix StandCarbonBalance(numDays, 6);
  
  
  //Add matrix dimnames
  LabileCarbonBalance.attr("dimnames") = List::create(dateStrings, above.attr("row.names"));
  GrossPhotosynthesis.attr("dimnames") = List::create(dateStrings, above.attr("row.names"));
  MaintenanceRespiration.attr("dimnames") = List::create(dateStrings, above.attr("row.names"));
  GrowthCosts.attr("dimnames") = List::create(dateStrings, above.attr("row.names"));
  PlantSugarLeaf.attr("dimnames") = List::create(dateStrings, above.attr("row.names"));
  PlantStarchLeaf.attr("dimnames") = List::create(dateStrings, above.attr("row.names")) ;
  PlantSugarSapwood.attr("dimnames") = List::create(dateStrings, above.attr("row.names")) ;
  PlantStarchSapwood.attr("dimnames") = List::create(dateStrings, above.attr("row.names")) ;
  PlantSugarTransport.attr("dimnames") = List::create(dateStrings, above.attr("row.names")) ;
  SapwoodBiomass.attr("dimnames") = List::create(dateStrings, above.attr("row.names")) ;
  LeafBiomass.attr("dimnames") = List::create(dateStrings, above.attr("row.names")) ;
  FineRootArea.attr("dimnames") = List::create(dateStrings, above.attr("row.names")) ;
  SapwoodArea.attr("dimnames") = List::create(dateStrings, above.attr("row.names")) ;
  LeafArea.attr("dimnames") = List::create(dateStrings, above.attr("row.names")) ;
  HuberValue.attr("dimnames") = List::create(dateStrings, above.attr("row.names")) ;
  RootAreaLeafArea.attr("dimnames") = List::create(dateStrings, above.attr("row.names")) ;
  FineRootBiomass.attr("dimnames") = List::create(dateStrings, above.attr("row.names"));
  DBH.attr("dimnames") = List::create(dateStrings, above.attr("row.names"));
  Height.attr("dimnames") = List::create(dateStrings, above.attr("row.names"));
  LAgrowth.attr("dimnames") = List::create(dateStrings, above.attr("row.names")) ;
  SAgrowth.attr("dimnames") = List::create(dateStrings, above.attr("row.names")) ;
  FRAgrowth.attr("dimnames") = List::create(dateStrings, above.attr("row.names"));
  dessicationRate.attr("dimnames") = List::create(dateStrings, above.attr("row.names")) ;
  mortalityRate.attr("dimnames") = List::create(dateStrings, above.attr("row.names")) ;
  starvationRate.attr("dimnames") = List::create(dateStrings, above.attr("row.names"));
  RootExudation.attr("dimnames") = List::create(dateStrings, above.attr("row.names"));
  StructuralBiomassBalance.attr("dimnames") = List::create(dateStrings, above.attr("row.names"));
  LabileBiomassBalance.attr("dimnames") = List::create(dateStrings, above.attr("row.names"));
  PlantBiomassBalance.attr("dimnames") = List::create(dateStrings, above.attr("row.names"));
  MortalityBiomassLoss.attr("dimnames") = List::create(dateStrings, above.attr("row.names"));
  CohortBiomassBalance.attr("dimnames") = List::create(dateStrings, above.attr("row.names"));
  
  StandBiomassBalance.attr("dimnames") = List::create(dateStrings, 
                           CharacterVector::create("StructuralBalance", "LabileBalance", "PlantBalance", "MortalityLoss", "CohortBalance"));
  StandCarbonBalance.attr("dimnames") = List::create(dateStrings, 
                          CharacterVector::create("GrossPrimaryProduction", "MaintenanceRespiration", "SynthesisRespiration", "NetPrimaryProduction",
                                                  "HeterotrophicRespiration", "NetEcosystemExchange"));
  // Assemble output
  List labileCarbonBalance = List::create(
    Named("GrossPhotosynthesis") = GrossPhotosynthesis,
    Named("MaintenanceRespiration") = MaintenanceRespiration,
    Named("GrowthCosts") = GrowthCosts,
    Named("RootExudation") = RootExudation,
    Named("LabileCarbonBalance") = LabileCarbonBalance,
    Named("SugarLeaf") = PlantSugarLeaf,
    Named("StarchLeaf") = PlantStarchLeaf,
    Named("SugarSapwood") = PlantSugarSapwood,
    Named("StarchSapwood") = PlantStarchSapwood,
    Named("SugarTransport") = PlantSugarTransport
  );
  List plantBiomassBalance = List::create(_["StructuralBiomassBalance"] = StructuralBiomassBalance,
                                          _["LabileBiomassBalance"] = LabileBiomassBalance,
                                          _["PlantBiomassBalance"] = PlantBiomassBalance,
                                          _["MortalityBiomassLoss"] = MortalityBiomassLoss,
                                          _["CohortBiomassBalance"] = CohortBiomassBalance);

  
  List growthMortality, plantStructure;
  plantStructure = List::create(Named("LeafBiomass")=LeafBiomass,
                                Named("SapwoodBiomass") = SapwoodBiomass,
                                Named("FineRootBiomass") = FineRootBiomass,
                                Named("LeafArea") = LeafArea,
                                Named("SapwoodArea")=SapwoodArea,
                                Named("FineRootArea") = FineRootArea,
                                Named("HuberValue") = HuberValue,
                                Named("RootAreaLeafArea") = RootAreaLeafArea,
                                Named("DBH") = DBH,
                                Named("Height") = Height);
  growthMortality = List::create(Named("LAgrowth") = LAgrowth,
                                 Named("SAgrowth") = SAgrowth,
                                 Named("FRAgrowth") = FRAgrowth,
                                 Named("StarvationRate") = starvationRate,
                                 Named("DessicationRate") = dessicationRate,
                                 Named("MortalityRate") = mortalityRate);

  NumericMatrix DecompositionPools(numDays, 10);
  DecompositionPools.attr("dimnames") = List::create(dateStrings, 
                          CharacterVector::create("SurfaceSnags","SurfaceLitter", "SoilLitter", "SurfaceMetabolic", "SoilMetabolic",
                                                  "SurfaceActive", "SoilActive", "SurfaceSlow", "SoilSlow", "SoilPassive"));
  
  List l;
  if(transpirationMode=="Granier") {
    l = List::create(Named("latitude") = latitude,
                     Named("topography") = topo,
                     Named("weather") = NA_REAL,
                     Named("growthInput") = growthInput,
                     Named("growthOutput") = x,
                     Named("WaterBalance")= DWB, 
                     Named("CarbonBalance")=StandCarbonBalance, 
                     Named("BiomassBalance") = StandBiomassBalance);
    if(control["soilResults"]) l.push_back(Soil, "Soil");
    if(control["snowResults"]) l.push_back(Snow, "Snow");
    if(control["standResults"]) l.push_back(Stand, "Stand");
    if(control["plantResults"]) l.push_back(plantDWOL, "Plants");
    if(control["labileCarbonBalanceResults"]) l.push_back(labileCarbonBalance, "LabileCarbonBalance");
    l.push_back(plantBiomassBalance, "PlantBiomassBalance");
    if(control["plantStructureResults"]) l.push_back(plantStructure, "PlantStructure");
    if(control["growthMortalityResults"]) l.push_back(growthMortality, "GrowthMortality");
    if(control["decompositionPoolResults"]) l.push_back(DecompositionPools, "DecompositionPools"); 
    if(control["fireHazardResults"]) {
      DataFrame fireHazard = defineFireHazardOutput(dateStrings);
      l.push_back(fireHazard, "FireHazard");
    }
  } else {
    
    DataFrame canopy = Rcpp::as<Rcpp::DataFrame>(x["canopy"]);
    DataFrame DEB = defineEnergyBalanceDailyOutput(dateStrings);
    DataFrame DT = defineTemperatureDailyOutput(dateStrings);
    NumericMatrix DLT =  defineTemperatureLayersDailyOutput(dateStrings, canopy);
    l = List::create(Named("latitude") = latitude,
                     Named("topography") = topo,
                     Named("weather") = NA_REAL,
                     Named("growthInput") = growthInput,
                     Named("growthOutput") = x,
                     Named("WaterBalance")= DWB, 
                     Named("CarbonBalance")=StandCarbonBalance, 
                     Named("EnergyBalance") = DEB,
                     Named("BiomassBalance") = StandBiomassBalance);
    if(control["temperatureResults"]) {
      l.push_back(DT, "Temperature");
      if(control["multiLayerBalance"]) l.push_back(DLT,"TemperatureLayers");
    }
    if(control["soilResults"]) l.push_back(Soil, "Soil");
    if(control["snowResults"]) l.push_back(Snow, "Snow");
    if(control["standResults"]) l.push_back(Stand, "Stand");
    if(control["plantResults"]) l.push_back(plantDWOL, "Plants");
    if(control["labileCarbonBalanceResults"]) l.push_back(labileCarbonBalance, "LabileCarbonBalance");
    l.push_back(plantBiomassBalance, "PlantBiomassBalance");
    if(control["plantStructureResults"]) l.push_back(plantStructure, "PlantStructure");
    if(control["growthMortalityResults"]) l.push_back(growthMortality, "GrowthMortality");
    if(control["decompositionPoolResults"]) l.push_back(DecompositionPools, "DecompositionPools"); 
    if(control["leafResults"]) {
      l.push_back(sunlitDO, "SunlitLeaves");
      l.push_back(shadeDO, "ShadeLeaves");
    }
    if(control["fireHazardResults"]) {
      DataFrame fireHazard = defineFireHazardOutput(dateStrings);
      l.push_back(fireHazard, "FireHazard");
    }
  }
  if(control["subdailyResults"]) l.push_back(subdailyRes,"subdaily");
  
  l.attr("class") = CharacterVector::create("growth","list");
  return(l);
}


// [[Rcpp::export(".fillGrowthDailyOutput")]]
void fillGrowthDailyOutput(List l, List x, List sDay, int iday) {
  
  List control = x["control"];
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  int numCohorts = above.nrow();
  int ntimesteps = control["ndailysteps"];
  DataFrame canopyParams = Rcpp::as<Rcpp::DataFrame>(x["canopy"]);
  int ncanlayers = canopyParams.nrow(); //Number of canopy layers
  DataFrame soil = Rcpp::as<Rcpp::DataFrame>(x["soil"]);
  String transpirationMode = control["transpirationMode"];
  
  DataFrame DWB = Rcpp::as<Rcpp::DataFrame>(l["WaterBalance"]);
  int numDays = DWB.nrow();
  fillWaterBalanceDailyOutput(DWB, sDay, iday, transpirationMode);
  
  if(control["soilResults"]) {
    String soilFunctions = control["soilFunctions"];
    List Soil = Rcpp::as<Rcpp::List>(l["Soil"]);
    fillSoilDailyOutput(Soil, soil, sDay, 
                        iday, numDays, soilFunctions,
                        true);
  }
  if(control["snowResults"]) {
    DataFrame Snow = Rcpp::as<Rcpp::DataFrame>(l["Snow"]);
    fillSnowDailyOutput(Snow, x, iday);
  }
  if(control["standResults"]) {
    DataFrame Stand = Rcpp::as<Rcpp::DataFrame>(l["Stand"]);
    fillStandDailyOutput(Stand, sDay,iday); 
  }
  if(control["plantResults"]) {
    List plantDWOL = l["Plants"];
    fillPlantWaterDailyOutput(plantDWOL, sDay, iday, transpirationMode); 
    if(transpirationMode!= "Granier") {
      List sunlitDO = l["SunlitLeaves"];
      List shadeDO = l["ShadeLeaves"];
      fillSunlitShadeLeavesDailyOutput(sunlitDO, shadeDO, sDay, iday, numCohorts); 
    } 
  }
  if(transpirationMode!= "Granier") {
    List DEB = l["EnergyBalance"];
    fillEnergyBalanceDailyOutput(DEB,sDay, iday, ntimesteps);
    if(control["temperatureResults"]) {
      List DT = l["Temperature"];
      fillTemperatureDailyOutput(DT,sDay, iday, ntimesteps);
      if(control["multiLayerBalance"]) {
        NumericMatrix DLT = l["TemperatureLayers"];
        fillTemperatureLayersDailyOutput(DLT,sDay, iday, ncanlayers, ntimesteps);
      }
    }
  } 
  if(control["fireHazardResults"]) {
    DataFrame fireHazard = Rcpp::as<Rcpp::DataFrame>(l["FireHazard"]);
    fillFireHazardOutput(fireHazard, sDay, iday);
  }
  
  if(control["subdailyResults"]) {
    List subdailyRes = Rcpp::as<Rcpp::List>(l["subdaily"]);
    subdailyRes[iday] = copyAdvancedGROWTHOutput(sDay, x); //Clones subdaily results because they are communication structures
  }
  
  List sb = sDay["Soil"];
  List db = sDay["WaterBalance"];
  List Plants = sDay["Plants"];
  DataFrame gm = Rcpp::as<Rcpp::DataFrame>(sDay["GrowthMortality"]);
  DataFrame cb = Rcpp::as<Rcpp::DataFrame>(sDay["LabileCarbonBalance"]);
  DataFrame bb = Rcpp::as<Rcpp::DataFrame>(sDay["PlantBiomassBalance"]);
  DataFrame ps = Rcpp::as<Rcpp::DataFrame>(sDay["PlantStructure"]);
  

  if(control["labileCarbonBalanceResults"]) {
    List labileCarbonBalance = Rcpp::as<Rcpp::List>(l["LabileCarbonBalance"]);
    NumericMatrix LabileCarbonBalance = Rcpp::as<Rcpp::NumericMatrix>(labileCarbonBalance["LabileCarbonBalance"]);
    NumericMatrix GrossPhotosynthesis = Rcpp::as<Rcpp::NumericMatrix>(labileCarbonBalance["GrossPhotosynthesis"]);
    NumericMatrix GrowthCosts = Rcpp::as<Rcpp::NumericMatrix>(labileCarbonBalance["GrowthCosts"]);
    NumericMatrix RootExudation = Rcpp::as<Rcpp::NumericMatrix>(labileCarbonBalance["RootExudation"]);
    NumericMatrix MaintenanceRespiration = Rcpp::as<Rcpp::NumericMatrix>(labileCarbonBalance["MaintenanceRespiration"]);
    NumericMatrix PlantSugarLeaf = Rcpp::as<Rcpp::NumericMatrix>(labileCarbonBalance["SugarLeaf"]);
    NumericMatrix PlantStarchLeaf = Rcpp::as<Rcpp::NumericMatrix>(labileCarbonBalance["StarchLeaf"]);
    NumericMatrix PlantSugarSapwood = Rcpp::as<Rcpp::NumericMatrix>(labileCarbonBalance["SugarSapwood"]);
    NumericMatrix PlantStarchSapwood = Rcpp::as<Rcpp::NumericMatrix>(labileCarbonBalance["StarchSapwood"]);
    NumericMatrix PlantSugarTransport = Rcpp::as<Rcpp::NumericMatrix>(labileCarbonBalance["SugarTransport"]);
    NumericVector LabileCarbonBalanceIN = Rcpp::as<Rcpp::NumericVector>(cb["LabileCarbonBalance"]);
    NumericVector MaintenanceRespirationIN = Rcpp::as<Rcpp::NumericVector>(cb["MaintenanceRespiration"]);
    NumericVector GrowthCostsIN = Rcpp::as<Rcpp::NumericVector>(cb["GrowthCosts"]);
    NumericVector GrossPhotosynthesisIN = Rcpp::as<Rcpp::NumericVector>(cb["GrossPhotosynthesis"]);
    NumericVector PlantSugarLeafIN = Rcpp::as<Rcpp::NumericVector>(cb["SugarLeaf"]);
    NumericVector PlantStarchLeafIN = Rcpp::as<Rcpp::NumericVector>(cb["StarchLeaf"]);
    NumericVector PlantSugarSapwoodIN = Rcpp::as<Rcpp::NumericVector>(cb["SugarSapwood"]);
    NumericVector PlantStarchSapwoodIN = Rcpp::as<Rcpp::NumericVector>(cb["StarchSapwood"]);
    NumericVector PlantSugarTransportIN = Rcpp::as<Rcpp::NumericVector>(cb["SugarTransport"]);
    NumericVector RootExudationIN = Rcpp::as<Rcpp::NumericVector>(cb["RootExudation"]);
    for(int i =0;i<numCohorts;i++) {
      LabileCarbonBalance(iday,i) = LabileCarbonBalanceIN[i];
      MaintenanceRespiration(iday,i) = MaintenanceRespirationIN[i];
      GrowthCosts(iday,i) = GrowthCostsIN[i];
      GrossPhotosynthesis(iday,i) = GrossPhotosynthesisIN[i];
      PlantSugarLeaf(iday,i) = PlantSugarLeafIN[i];
      PlantStarchLeaf(iday,i) = PlantStarchLeafIN[i];
      PlantSugarSapwood(iday,i) = PlantSugarSapwoodIN[i];
      PlantStarchSapwood(iday,i) = PlantStarchSapwoodIN[i];
      PlantSugarTransport(iday,i) = PlantSugarTransportIN[i];
      RootExudation(iday,i) = RootExudationIN[i];
    }
  }
  
  if(control["plantStructureResults"]) {
    List plantStructure = Rcpp::as<Rcpp::List>(l["PlantStructure"]);
    NumericMatrix LeafBiomass = Rcpp::as<Rcpp::NumericMatrix>(plantStructure["LeafBiomass"]);
    NumericMatrix SapwoodBiomass = Rcpp::as<Rcpp::NumericMatrix>(plantStructure["SapwoodBiomass"]);
    NumericMatrix FineRootBiomass = Rcpp::as<Rcpp::NumericMatrix>(plantStructure["FineRootBiomass"]);
    NumericMatrix LeafArea = Rcpp::as<Rcpp::NumericMatrix>(plantStructure["LeafArea"]);
    NumericMatrix SapwoodArea = Rcpp::as<Rcpp::NumericMatrix>(plantStructure["SapwoodArea"]);
    NumericMatrix FineRootArea = Rcpp::as<Rcpp::NumericMatrix>(plantStructure["FineRootArea"]);
    NumericMatrix HuberValue = Rcpp::as<Rcpp::NumericMatrix>(plantStructure["HuberValue"]);
    NumericMatrix RootAreaLeafArea = Rcpp::as<Rcpp::NumericMatrix>(plantStructure["RootAreaLeafArea"]);
    NumericMatrix DBH = Rcpp::as<Rcpp::NumericMatrix>(plantStructure["DBH"]);
    NumericMatrix Height = Rcpp::as<Rcpp::NumericMatrix>(plantStructure["Height"]);
    
    NumericVector SapwoodBiomassIN = Rcpp::as<Rcpp::NumericVector>(ps["SapwoodBiomass"]);
    NumericVector LeafBiomassIN = Rcpp::as<Rcpp::NumericVector>(ps["LeafBiomass"]);
    NumericVector FineRootBiomassIN = Rcpp::as<Rcpp::NumericVector>(ps["FineRootBiomass"]);
    NumericVector SapwoodAreaIN = Rcpp::as<Rcpp::NumericVector>(ps["SapwoodArea"]);
    NumericVector LeafAreaIN = Rcpp::as<Rcpp::NumericVector>(ps["LeafArea"]);
    NumericVector FineRootAreaIN = Rcpp::as<Rcpp::NumericVector>(ps["FineRootArea"]);
    NumericVector HuberValueIN = Rcpp::as<Rcpp::NumericVector>(ps["HuberValue"]);
    NumericVector RootAreaLeafAreaIN = Rcpp::as<Rcpp::NumericVector>(ps["RootAreaLeafArea"]);
    NumericVector DBHIN = Rcpp::as<Rcpp::NumericVector>(ps["DBH"]);
    NumericVector HeightIN = Rcpp::as<Rcpp::NumericVector>(ps["Height"]);
    
    for(int i =0;i<numCohorts;i++) {
      SapwoodBiomass(iday,i) = SapwoodBiomassIN[i];
      LeafBiomass(iday,i) = LeafBiomassIN[i];
      FineRootBiomass(iday,i) = FineRootBiomassIN[i];
      SapwoodArea(iday,i) = SapwoodAreaIN[i];
      LeafArea(iday,i) = LeafAreaIN[i];
      FineRootArea(iday,i) = FineRootAreaIN[i];
      HuberValue(iday,i) = HuberValueIN[i];
      RootAreaLeafArea(iday,i) = RootAreaLeafAreaIN[i];
      DBH(iday,i) = DBHIN[i];
      Height(iday,i) = HeightIN[i];
    }
  }
  
  List plantBiomassBalance = Rcpp::as<Rcpp::List>(l["PlantBiomassBalance"]);
  NumericMatrix StructuralBiomassBalance = Rcpp::as<Rcpp::NumericMatrix>(plantBiomassBalance["StructuralBiomassBalance"]);
  NumericMatrix LabileBiomassBalance = Rcpp::as<Rcpp::NumericMatrix>(plantBiomassBalance["LabileBiomassBalance"]);
  NumericMatrix PlantBiomassBalance = Rcpp::as<Rcpp::NumericMatrix>(plantBiomassBalance["PlantBiomassBalance"]);
  NumericMatrix MortalityBiomassLoss = Rcpp::as<Rcpp::NumericMatrix>(plantBiomassBalance["MortalityBiomassLoss"]);
  NumericMatrix CohortBiomassBalance = Rcpp::as<Rcpp::NumericMatrix>(plantBiomassBalance["CohortBiomassBalance"]);
  
  NumericVector StructuralBiomassBalanceIN = Rcpp::as<Rcpp::NumericVector>(bb["StructuralBiomassBalance"]);
  NumericVector LabileBiomassBalanceIN = Rcpp::as<Rcpp::NumericVector>(bb["LabileBiomassBalance"]);
  NumericVector PlantBiomassBalanceIN = Rcpp::as<Rcpp::NumericVector>(bb["PlantBiomassBalance"]);
  NumericVector MortalityBiomassLossIN = Rcpp::as<Rcpp::NumericVector>(bb["MortalityBiomassLoss"]);
  NumericVector CohortBiomassBalanceIN = Rcpp::as<Rcpp::NumericVector>(bb["CohortBiomassBalance"]);

  for(int i =0;i<numCohorts;i++) {
    StructuralBiomassBalance(iday,i) = StructuralBiomassBalanceIN[i];
    LabileBiomassBalance(iday,i) = LabileBiomassBalanceIN[i];
    PlantBiomassBalance(iday,i) = PlantBiomassBalanceIN[i];
    MortalityBiomassLoss(iday,i) = MortalityBiomassLossIN[i];
    CohortBiomassBalance(iday,i) = CohortBiomassBalanceIN[i];
  }
  
  if(control["growthMortalityResults"]) {
    List growthMortality = Rcpp::as<Rcpp::List>(l["GrowthMortality"]);
    NumericMatrix LAgrowth = Rcpp::as<Rcpp::NumericMatrix>(growthMortality["LAgrowth"]);
    NumericMatrix SAgrowth = Rcpp::as<Rcpp::NumericMatrix>(growthMortality["SAgrowth"]);
    NumericMatrix FRAgrowth = Rcpp::as<Rcpp::NumericMatrix>(growthMortality["FRAgrowth"]);
    NumericMatrix starvationRate = Rcpp::as<Rcpp::NumericMatrix>(growthMortality["StarvationRate"]);
    NumericMatrix dessicationRate = Rcpp::as<Rcpp::NumericMatrix>(growthMortality["DessicationRate"]);
    NumericMatrix mortalityRate = Rcpp::as<Rcpp::NumericMatrix>(growthMortality["MortalityRate"]);
    NumericVector LAgrowthIN = Rcpp::as<Rcpp::NumericVector>(gm["LAgrowth"]);
    NumericVector SAgrowthIN = Rcpp::as<Rcpp::NumericVector>(gm["SAgrowth"]);
    NumericVector FRAgrowthIN = Rcpp::as<Rcpp::NumericVector>(gm["FRAgrowth"]);
    NumericVector starvationRateIN = Rcpp::as<Rcpp::NumericVector>(gm["StarvationRate"]);
    NumericVector dessicationRateIN = Rcpp::as<Rcpp::NumericVector>(gm["DessicationRate"]);
    NumericVector mortalityRateIN = Rcpp::as<Rcpp::NumericVector>(gm["MortalityRate"]);
    
    for(int i =0;i<numCohorts;i++) {
      LAgrowth(iday,i) = LAgrowthIN[i];
      SAgrowth(iday,i) = SAgrowthIN[i];
      FRAgrowth(iday,i) = FRAgrowthIN[i];
      starvationRate(iday,i) = starvationRateIN[i];
      dessicationRate(iday,i) = dessicationRateIN[i];
      mortalityRate(iday,i) = mortalityRateIN[i];
    }
  }
  
  if(control["decompositionPoolResults"]) {
    NumericMatrix decompositionPool = Rcpp::as<Rcpp::NumericMatrix>(l["DecompositionPools"]);
    NumericVector internalSOC = x["internalSOC"];
    DataFrame internalSnags = Rcpp::as<Rcpp::DataFrame>(x["internalSnags"]);
    NumericVector structural_snags_smallbranches = internalSnags["SmallBranches"]; 
    NumericVector structural_snags_largewood = internalSnags["LargeWood"]; 
    DataFrame internalLitter = Rcpp::as<Rcpp::DataFrame>(x["internalLitter"]);
    NumericVector structural_litter_leaves = internalLitter["Leaves"]; 
    NumericVector structural_litter_twigs = internalLitter["Twigs"]; 
    NumericVector structural_litter_smallbranches = internalLitter["SmallBranches"]; 
    NumericVector structural_litter_largewood = internalLitter["LargeWood"]; 
    NumericVector structural_litter_coarseroots = internalLitter["CoarseRoots"]; 
    NumericVector structural_litter_fineroots = internalLitter["FineRoots"]; 
    decompositionPool(iday,0) = sum(structural_snags_smallbranches) + sum(structural_snags_largewood);
    decompositionPool(iday,1) = sum(structural_litter_leaves) + sum(structural_litter_twigs) +sum(structural_litter_smallbranches) + sum(structural_litter_largewood);
    decompositionPool(iday,2) = sum(structural_litter_fineroots)+ sum(structural_litter_coarseroots);
    decompositionPool(iday,3) = internalSOC["SurfaceMetabolic"];
    decompositionPool(iday,4) = internalSOC["SoilMetabolic"];
    decompositionPool(iday,5) = internalSOC["SurfaceActive"];
    decompositionPool(iday,6) = internalSOC["SoilActive"];
    decompositionPool(iday,7) = internalSOC["SurfaceSlow"];
    decompositionPool(iday,8) = internalSOC["SoilSlow"];
    decompositionPool(iday,9) = internalSOC["SoilPassive"];
  }
  
  NumericMatrix StandBiomassBalance = Rcpp::as<Rcpp::NumericMatrix>(l["BiomassBalance"]);
  StandBiomassBalance(iday,_) = standLevelBiomassBalance(bb);
  
  NumericMatrix StandCarbonBalance = Rcpp::as<Rcpp::NumericMatrix>(l["CarbonBalance"]);
  StandCarbonBalance(iday,_) = Rcpp::as<Rcpp::NumericVector>(sDay["CarbonBalance"]);
  
}

//' Forest growth
//' 
//' Function \code{growth} is a process-based model that performs energy, water and carbon balances; 
//' and determines changes in water/carbon pools, functional variables (leaf area, sapwood area, root area) 
//' and structural ones (tree diameter, tree height, shrub cover) for woody plant cohorts in a given forest stand 
//' during a period specified in the input climatic data. 
//' 
//' @param x An object of class \code{\link{growthInput}}.
//' @param meteo A data frame with daily meteorological data series (see \code{\link{spwb}}).
//' @param latitude Latitude (in degrees).
//' @param elevation,slope,aspect Elevation above sea level (in m), slope (in degrees) and aspect (in degrees from North). 
//' @param CO2ByYear A named numeric vector with years as names and atmospheric CO2 concentration (in ppm) as values. Used to specify annual changes in CO2 concentration along the simulation (as an alternative to specifying daily values in \code{meteo}).
//' @param waterTableDepth Water table depth (in mm). When not missing, capillarity rise will be allowed if lower than total soil depth.
//' 
//' @details
//' Detailed model description is available in the medfate book. 
//' 
//' Forest growth simulations allow using different sub-models for bulk soil water flows and different sub-models of transpiration and photosynthesis (see details in \code{\link{spwb}}). 
//' 
//' @return
//' Function \code{growth} returns a list of class 'growth'. Since lists are difficult to handle, we recommend using
//' function \code{\link{extract}} to reshape simulation results (including their units) from those objects. 
//' 
//' List elements are as follows:
//' \itemize{
//'   \item{\code{"latitude"}: Latitude (in degrees) given as input.} 
//'   \item{\code{"topography"}: Vector with elevation, slope and aspect given as input.} 
//'   \item{\code{"weather"}: A copy of the input weather data frame.}
//'   \item{\code{"growthInput"}: A copy of the object \code{x} of class \code{\link{growthInput}} given as input.}
//'   \item{\code{"growthOutput"}: An copy of the final state of the object \code{x} of class \code{\link{growthInput}}.}
//'   \item{\code{"WaterBalance"}: A data frame where different water balance variables (see \code{\link{spwb}}).}
//'   \item{\code{"EnergyBalance"}: A data frame with the daily values of energy balance components for the soil and the canopy (only for \code{transpirationMode = "Sperry"} or \code{transpirationMode = "Sureau"}; see \code{\link{spwb}}).}
//'   \item{\code{"CarbonBalance"}: A data frame where different stand-level carbon balance components (gross primary production, maintenance respiration, synthesis respiration, net primary production, heterotrophic respiration and net ecosystem exchange.), all in g C · m-2.}
//'   \item{\code{"BiomassBalance"}: A data frame with the daily values of stand biomass balance components (in g dry · m-2.}
//'   \item{\code{"Temperature"}: A data frame with the daily values of minimum/mean/maximum temperatures for the atmosphere (input), canopy and soil (only for \code{transpirationMode = "Sperry"} or \code{transpirationMode = "Sureau"}; see \code{\link{spwb}}).}
//'   \item{\code{"Soil"}: A data frame where different soil variables  (see \code{\link{spwb}}).}
//'   \item{\code{"Stand"}: A data frame where different stand-level variables (see \code{\link{spwb}}).}
//'   \item{\code{"Plants"}: A list of daily results for plant cohorts (see \code{\link{spwb}}).}
//'   \item{\code{"SunlitLeaves"} and \code{"ShadeLeaves"}: A list with daily results for sunlit and shade leaves (only for \code{transpirationMode = "Sperry"} or \code{transpirationMode = "Sureau"}; see \code{\link{spwb}}).}
//'   \item{\code{"LabileCarbonBalance"}: A list of daily labile carbon balance results for plant cohorts, with elements:}
//'   \itemize{
//'     \item{\code{"GrossPhotosynthesis"}: Daily gross photosynthesis per dry weight of living biomass (g gluc · g dry-1).}
//'     \item{\code{"MaintentanceRespiration"}: Daily maintenance respiration per dry weight of living biomass (g gluc · g dry-1).}
//'     \item{\code{"GrowthCosts"}: Daily growth costs per dry weight of living biomass (g gluc · g dry-1).}
//'     \item{\code{"RootExudation"}: Root exudation per dry weight of living biomass (g gluc · g dry-1).}    
//'     \item{\code{"LabileCarbonBalance"}: Daily labile carbon balance (photosynthesis - maintenance respiration - growth costs - root exudation) per dry weight of living biomass (g gluc · g dry-1).}
//'     \item{\code{"SugarLeaf"}: Sugar concentration (mol·l-1) in leaves.}
//'     \item{\code{"StarchLeaf"}: Starch concentration (mol·l-1) in leaves.}
//'     \item{\code{"SugarSapwood"}: Sugar concentration (mol·l-1) in sapwood.}
//'     \item{\code{"StarchSapwood"}: Starch concentration (mol·l-1) in sapwood.}
//'     \item{\code{"SugarTransport"}:  Average instantaneous rate of carbon transferred between leaves and stem compartments via floem (mol gluc·s-1).}
//'   }
//'   \item{\code{"PlantBiomassBalance"}: A list of daily plant biomass balance results for plant cohorts, with elements:}
//'   \itemize{
//'     \item{\code{"StructuralBiomassBalance"}: Daily structural biomass balance (g dry · m-2).}
//'     \item{\code{"LabileBiomassBalance"}: Daily labile biomass balance (g dry · m-2).}
//'     \item{\code{"PlantBiomassBalance"}: Daily plant biomass balance, i.e. labile change + structural change (g dry · m-2).}
//'     \item{\code{"MortalityBiomassLoss"}: Biomass loss due to mortality (g dry · m-2).}    
//'     \item{\code{"CohortBiomassBalance"}: Daily cohort biomass balance (including mortality) (g dry · m-2).}
//'   }
//'   \item{\code{"PlantStructure"}: A list of daily area and biomass values for compartments of plant cohorts, with elements:}
//'   \itemize{
//'     \item{\code{"LeafBiomass"}: Daily amount of leaf structural biomass (in g dry) for an average individual of each plant cohort.}
//'     \item{\code{"SapwoodBiomass"}: Daily amount of sapwood structural biomass (in g dry) for an average individual of each plant cohort.}
//'     \item{\code{"FineRootBiomass"}: Daily amount of fine root biomass (in g dry) for an average individual of each plant cohort.}
//'     \item{\code{"LeafArea"}: Daily amount of leaf area (in m2) for an average individual of each plant cohort.}
//'     \item{\code{"SapwoodArea"}: Daily amount of sapwood area (in cm2) for an average individual of each plant cohort.}
//'     \item{\code{"FineRootArea"}: Daily amount of fine root area (in m2) for an average individual of each plant cohort.}
//'     \item{\code{"HuberValue"}: The ratio of sapwood area to (target) leaf area (in cm2/m2).}
//'     \item{\code{"RootAreaLeafArea"}: The ratio of fine root area to (target) leaf area (in m2/m2).}
//'     \item{\code{"DBH"}: Diameter at breast height (in cm) for an average individual of each plant cohort.}
//'     \item{\code{"Height"}: Height (in cm) for an average individual of each plant cohort.}
//'   }
//'   \item{\code{"GrowthMortality"}: A list of daily growth and mortality rates for plant cohorts, with elements:}
//'   \itemize{
//'     \item{\code{"LAgrowth"}: Leaf area growth (in m2·day-1) for an average individual of each plant cohort.}
//'     \item{\code{"SAgrowth"}: Sapwood area growth rate (in cm2·day-1) for an average individual of each plant cohort.}
//'     \item{\code{"FRAgrowth"}: Fine root area growth (in m2·day-1) for an average individual of each plant cohort.}
//'     \item{\code{"StarvationRate"}: Daily mortality rate from starvation (ind/d-1).}
//'     \item{\code{"DessicationRate"}: Daily mortality rate from dessication (ind/d-1).}
//'     \item{\code{"MortalityRate"}: Daily mortality rate (any cause) (ind/d-1).}
//'   }
//'   \item{\code{"DecompositionPools"}: A data frame with the mass of different decomposition carbon pools, all in g C · m-2.}
//'   \item{\code{"subdaily"}: A list of objects of class \code{\link{growth_day}}, one per day simulated (only if required in \code{control} parameters, see \code{\link{defaultControl}}).}
//' }
//' 
//' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
//' 
//' @seealso \code{\link{growthInput}}, \code{\link{growth_day}}, \code{\link{plot.growth}}
//' 
//' @references
//' De Cáceres M, Molowny-Horas R, Cabon A, Martínez-Vilalta J, Mencuccini M, García-Valdés R, Nadal-Sala D, Sabaté S, 
//' Martin-StPaul N, Morin X, D'Adamo F, Batllori E, Améztegui A (2023) MEDFATE 2.9.3: A trait-enabled model to simulate 
//' Mediterranean forest function and dynamics at regional scales. 
//' Geoscientific Model Development 16: 3165-3201 (https://doi.org/10.5194/gmd-16-3165-2023).
//' 
//' @examples
//' \donttest{
//' #Load example daily meteorological data
//' data(examplemeteo)
//' 
//' #Load example plot plant data
//' data(exampleforest)
//' 
//' #Default species parameterization
//' data(SpParamsMED)
//'   
//' #Initialize control parameters
//' control <- defaultControl("Granier")
//'   
//' #Initialize soil with default soil params (4 layers)
//' examplesoil <- defaultSoilParams(4)
//' 
//' #Initialize model input
//' x1 <- growthInput(exampleforest, examplesoil, SpParamsMED, control)
//' 
//' #Call simulation function
//' G1 <- growth(x1, examplemeteo, latitude = 41.82592, elevation = 100)
//'  
//' #Switch to 'Sperry' transpiration mode
//' control <- defaultControl("Sperry")
//' 
//' #Initialize model input
//' x2 <- growthInput(exampleforest,examplesoil, SpParamsMED, control)
//' 
//' #Call simulation function
//' G2 <-growth(x2, examplemeteo, latitude = 41.82592, elevation = 100)
//' 
//' #Switch to 'Sureau' transpiration mode
//' control <- defaultControl("Sureau")
//' 
//' #Initialize model input
//' x3 <- growthInput(exampleforest,examplesoil, SpParamsMED, control)
//' 
//' #Call simulation function
//' G3 <-growth(x3, examplemeteo, latitude = 41.82592, elevation = 100)
//' }
//'       
// [[Rcpp::export("growth")]]
List growth(List x, DataFrame meteo, double latitude, 
            double elevation, double slope = NA_REAL, double aspect = NA_REAL,
            NumericVector CO2ByYear = NumericVector(0), double waterTableDepth = NA_REAL) {

  //Clone input
  x = clone(x);
  
  //Control params 
  List control =x["control"];  
  String transpirationMode = control["transpirationMode"];
  String soilFunctions = control["soilFunctions"];
  bool verbose = control["verbose"];
  bool leafPhenology = control["leafPhenology"];
  bool unlimitedSoilWater = control["unlimitedSoilWater"];
  NumericVector defaultRainfallIntensityPerMonth = control["defaultRainfallIntensityPerMonth"];
  checkgrowthInput(x, transpirationMode, soilFunctions);


  //Cohort info
  DataFrame cohorts = Rcpp::as<Rcpp::DataFrame>(x["cohorts"]);
  IntegerVector SP = Rcpp::as<Rcpp::IntegerVector>(cohorts["SP"]);
  IntegerVector SPunique = uniqueSpp(SP);

  if(NumericVector::is_na(latitude)) stop("Value for 'latitude' should not be missing.");
  double latrad = latitude * (M_PI/180.0);
  if(NumericVector::is_na(aspect)) aspect = 0.0;
  if(NumericVector::is_na(slope)) slope = 0.0;
  double asprad = aspect * (M_PI/180.0);
  double slorad = slope * (M_PI/180.0);
  
  //Meteorological input    
  NumericVector MinTemperature, MaxTemperature;
  NumericVector MinRelativeHumidity, MaxRelativeHumidity;
  NumericVector Radiation;
  if(!meteo.containsElementNamed("Precipitation")) stop("Please include variable 'Precipitation' in weather input.");
  NumericVector Precipitation = meteo["Precipitation"];
  int numDays = Precipitation.size();
  if(NumericVector::is_na(elevation)) stop("Value for 'elevation' should not be missing.");
  if(!meteo.containsElementNamed("MinTemperature")) stop("Please include variable 'MinTemperature' in weather input.");
  MinTemperature = meteo["MinTemperature"];
  if(!meteo.containsElementNamed("MaxTemperature")) stop("Please include variable 'MaxTemperature' in weather input.");
  MaxTemperature = meteo["MaxTemperature"];
  if(!meteo.containsElementNamed("MinRelativeHumidity")) stop("Please include variable 'MinRelativeHumidity' in weather input.");
  MinRelativeHumidity = meteo["MinRelativeHumidity"];
  if(!meteo.containsElementNamed("MaxRelativeHumidity")) stop("Please include variable 'MaxRelativeHumidity' in weather input.");
  MaxRelativeHumidity = meteo["MaxRelativeHumidity"];
  if(!meteo.containsElementNamed("Radiation")) stop("Please include variable 'Radiation' in weather input.");
  Radiation = meteo["Radiation"];
  
  if(any(is_na(Precipitation))) stop("Missing values in 'Precipitation' are not allowed");
  if(any(is_na(MinTemperature))) stop("Missing values in 'MinTemperature' are not allowed");
  if(any(is_na(MaxTemperature))) stop("Missing values in 'MaxTemperature' are not allowed");
  if(any(is_na(MinRelativeHumidity))) warning("Missing values in 'MinRelativeHumidity' were estimated from temperature range");
  if(any(is_na(MaxRelativeHumidity))) warning("Missing values in 'MaxRelativeHumidity' were assumed to be 100");
  if(any(is_na(Radiation))) warning("Missing values in 'Radiation' were estimated");
  
  NumericVector WindSpeed(numDays, NA_REAL);
  if(meteo.containsElementNamed("WindSpeed")) WindSpeed = meteo["WindSpeed"];
  
  NumericVector FireProbability(numDays, 0.0);
  if(meteo.containsElementNamed("FireProbability")) FireProbability = meteo["FireProbability"];
  
  NumericVector PET = NumericVector(numDays, NA_REAL);
  
  NumericVector CO2(Precipitation.length(), NA_REAL);
  if(meteo.containsElementNamed("CO2")) {
    CO2 = meteo["CO2"];
    if(verbose) {
      Rcout<<"CO2 taken from input column 'CO2'\n";
    }
    if(any(is_na(CO2))) stop("Missing values in 'CO2'");
  }
  NumericVector Patm(Precipitation.length(), NA_REAL);
  if(meteo.containsElementNamed("Patm")) {
    Patm = meteo["Patm"];
    if(verbose) {
      Rcout<<"Patm taken from input column 'Patm'\n";
    }
  }
  NumericVector RainfallIntensity(Precipitation.length(), NA_REAL);
  if(meteo.containsElementNamed("RainfallIntensity")) {
    RainfallIntensity = meteo["RainfallIntensity"];
    if(verbose) {
      Rcout<<"Rainfall intensity taken from input column 'RainfallIntensity'\n";
    }
  }
  
  IntegerVector DOY, JulianDay;
  NumericVector Photoperiod;
  bool doy_input = false, photoperiod_input = false, julianday_input = false;
  if(meteo.containsElementNamed("DOY")) {
    DOY = meteo["DOY"];
    doy_input = true;
    if(verbose) {
      Rcout<<"DOY taken from input column 'DOY'\n";
    }
  }
  if(meteo.containsElementNamed("Photoperiod")) {
    Photoperiod = meteo["Photoperiod"];
    photoperiod_input = true;
    if(verbose) {
      Rcout<<"Photoperiod taken from input column 'Photoperiod'\n";
    }
  }
  if(meteo.containsElementNamed("JulianDay")) {
    JulianDay = meteo["JulianDay"];
    julianday_input = true;
    if(verbose) {
      Rcout<<"Julian day taken from input column 'JulianDay'\n";
    }
  }
  
  // Dates
  CharacterVector dateStrings = getWeatherDates(meteo);
  if(!doy_input) DOY = date2doy(dateStrings);
  if(!photoperiod_input) Photoperiod = date2photoperiod(dateStrings, latrad);
  
  //Soil params 
  DataFrame soil = Rcpp::as<Rcpp::DataFrame>(x["soil"]);
  
  //Output list
  List outputList = defineGrowthDailyOutput(latitude, elevation, slope, aspect,
                                            dateStrings, x);
  outputList["weather"] = clone(meteo);

  //Count years (times structural variables will be updated)
  int numYears = 0;
  for(int i=0;i<numDays;i++) {
    if(((DOY[i]==1) && (i>0)) || ((i==(numDays-1)) && (DOY[i]>=365))) numYears = numYears + 1;
  }

  NumericVector initialSoilContent = water(soil, soilFunctions);
  NumericVector initialPlantContent = plantWaterContent(x);
  double initialSnowContent = x["snowpack"];
  
  
  DataFrame ccIni_m2 = carbonCompartments(x, "g_m2");
  double cohortBiomassBalanceSum = 0.0;
  double initialCohortBiomass = sum(Rcpp::as<Rcpp::NumericVector>(ccIni_m2["TotalBiomass"]));
  
  if(verbose) {
    Rcout<<"Initial plant cohort biomass (g/m2): "<<initialCohortBiomass<<"\n";
    Rcout<<"Initial plant water content (mm): "<< sum(initialPlantContent)<<"\n";
    Rcout<<"Initial soil water content (mm): "<< sum(initialSoilContent)<<"\n";
    Rcout<<"Initial snowpack content (mm): "<< initialSnowContent<<"\n";
  }
  
  //Instance communication structures
  List internalCommunication = instanceCommunicationStructures(x, "growth");

  bool error_occurence = false;
  if(verbose) Rcout << "Performing daily simulations\n";
  List s;
  std::string yearString;
  for(int i=0;i<numDays;i++) {
    std::string c = as<std::string>(dateStrings[i]);
    yearString = c.substr(0, 4);
    if(verbose) {
      if(DOY[i]==1 || i==0) {
        Rcout<<"\n Year "<< yearString<< ":";
      } 
      else if(i%30 == 0) Rcout<<".";//<<i;
    } 
    
    double wind = WindSpeed[i];
    if(NumericVector::is_na(wind)) wind = control["defaultWindSpeed"]; //Default 1 m/s -> 10% of fall every day
    if(wind<0.1) wind = 0.1; //Minimum windspeed abovecanopy
    
    double Catm = CO2[i];
    //If missing, use
    if(NumericVector::is_na(Catm)) {
      if(CO2ByYear.attr("names") != R_NilValue) Catm = CO2ByYear[yearString];
    }
    //If still missing, use default control value
    if(NumericVector::is_na(Catm)) {
      Catm = control["defaultCO2"];
    }
    
    double Rint = RainfallIntensity[i];
    if(NumericVector::is_na(Rint)) {
      int month = std::atoi(c.substr(5,2).c_str());
      Rint = rainfallIntensity(month, Precipitation[i], defaultRainfallIntensityPerMonth);
    }
    
    if(unlimitedSoilWater) {
      NumericVector W = soil["W"];
      for(int h=0;h<W.size();h++) W[h] = 1.0;
    }
    

    //Julian day from either input column or date
    int J = NA_INTEGER;
    if(julianday_input) J = JulianDay[i];
    if(IntegerVector::is_na(J)){
      std::string c = as<std::string>(dateStrings[i]);
      J = meteoland::radiation_julianDay(std::atoi(c.substr(0, 4).c_str()),std::atoi(c.substr(5,2).c_str()),std::atoi(c.substr(8,2).c_str())); 
    }
    double delta = meteoland::radiation_solarDeclination(J);
    double solarConstant = meteoland::radiation_solarConstant(J);
    
    double tmin = MinTemperature[i];
    double tmax = MaxTemperature[i];
    if(tmin > tmax) {
      warning("tmin > tmax. Swapping values.");
      double swap = tmin;
      tmin = tmax;
      tmax = swap;
    }
    double prec = Precipitation[i];
    double tday = meteoland::utils_averageDaylightTemperature(tmin, tmax);
    double rhmin = MinRelativeHumidity[i];
    double rhmax = MaxRelativeHumidity[i];
    double rad = Radiation[i];
    if(NumericVector::is_na(rhmax)) {
      rhmax = 100.0;
    }
    if(NumericVector::is_na(rhmin)) {
      double vp_tmin = meteoland::utils_saturationVP(tmin);
      double vp_tmax = meteoland::utils_saturationVP(tmax);
      rhmin = std::min(rhmax, 100.0*(vp_tmin/vp_tmax));
    }
    if(rhmin > rhmax) {
      warning("rhmin > rhmax. Swapping values.");
      double swap = rhmin;
      rhmin = rhmax;
      rhmax = swap;
    }
    if(NumericVector::is_na(rad)) {
      double vpa = meteoland::utils_averageDailyVP(tmin, tmax, rhmin, rhmax);
      rad = meteoland::radiation_solarRadiation(solarConstant, latrad, elevation,
                                                slorad, asprad, delta, tmax -tmin, tmax-tmin,
                                                vpa, prec);
    }
    PET[i] = meteoland::penman(latrad, elevation, slorad, asprad, J, 
                               tmin, tmax, rhmin, rhmax, rad, wind);
    
    //1. Phenology (only leaf fall)
    if(leafPhenology) {
      updatePhenology(x, DOY[i], Photoperiod[i], tday);
      updateLeaves(x, wind, true);
    }
    
    //2. Water balance and photosynthesis
    if(transpirationMode=="Granier") {
      NumericVector meteovec_inner = NumericVector::create(
        Named("tday") = tday, Named("tmax") = tmax, Named("tmin") = tmin,
        Named("prec") = prec, 
        Named("rhmin") = rhmin, Named("rhmax") = rhmax,
        Named("rad") = rad, 
        Named("wind") = wind, 
        Named("pet") = PET[i],
        Named("Catm") = Catm);
      meteovec_inner.push_back(Patm[i], "Patm");
      meteovec_inner.push_back(Rint, "rint");
      meteovec_inner.push_back(FireProbability[i], "pfire"); 
      try{
        growthDay_private(internalCommunication, x, meteovec_inner,  
                           latitude, elevation, slope, aspect,
                           solarConstant, delta, 
                           0.0, R_NilValue, waterTableDepth,
                           false); //No Runon in simulations for a single cell
        fillGrowthDailyOutput(outputList, x, internalCommunication["basicGROWTHOutput"], i);
      } catch(std::exception& ex) {
        Rcerr<< "c++ error: "<< ex.what() <<"\n";
        error_occurence = true;
      }
    } else {
      double tmaxPrev = tmax;
      double tminPrev = tmin;
      double tminNext = tmin;
      if(i>0) {
        tmaxPrev = MaxTemperature[i-1];
        tminPrev = MinTemperature[i-1];
      }
      if(i<(numDays-1)) tminNext = MinTemperature[i+1]; 
      NumericVector meteovec = NumericVector::create(
        Named("tday") = tday,
        Named("tmin") = tmin, 
        Named("tmax") = tmax,
        Named("tminPrev") = tminPrev, 
        Named("tmaxPrev") = tmaxPrev, 
        Named("tminNext") = tminNext, 
        Named("prec") = Precipitation[i],
        Named("rhmin") = rhmin, 
        Named("rhmax") = rhmax, 
        Named("rad") = rad, 
        Named("wind") = wind, 
        Named("Catm") = Catm,
        Named("Patm") = Patm[i],
        Named("pet") = PET[i],
        Named("rint") = Rint);
      meteovec.push_back(FireProbability[i], "pfire"); 
      try{
        growthDay_private(internalCommunication, x, meteovec, 
                           latitude, elevation, slope, aspect,
                           solarConstant, delta, 
                           0.0, R_NilValue, waterTableDepth,
                           verbose);
        fillGrowthDailyOutput(outputList, x, internalCommunication["advancedGROWTHOutput"], i);
      } catch(std::exception& ex) {
        Rcerr<< "c++ error: "<< ex.what() <<"\n";
        error_occurence = true;
      }
    }    
    
    //Add cohort biomass sum
    List plantBiomassBalance = Rcpp::as<Rcpp::List>(outputList["PlantBiomassBalance"]);
    NumericMatrix CohortBiomassBalance = Rcpp::as<Rcpp::NumericMatrix>(plantBiomassBalance["CohortBiomassBalance"]);
    cohortBiomassBalanceSum += sum(CohortBiomassBalance(i,_));
  }
  if(verbose) Rcout << "\n\n";
  
  // Check biomass balance
  DataFrame ccFin_m2 = carbonCompartments(x, "g_m2");
  double finalCohortBiomass = sum(Rcpp::as<Rcpp::NumericVector>(ccFin_m2["TotalBiomass"]));
  if(verbose) {
    NumericMatrix StandBiomassBalance = Rcpp::as<Rcpp::NumericMatrix>(outputList["BiomassBalance"]);
    
    Rcout<<"Final plant cohort biomass (g/m2): "<<finalCohortBiomass<<"\n";
    Rcout<<"Change in plant cohort biomass (g/m2): " << finalCohortBiomass - initialCohortBiomass <<"\n";
    Rcout<<"Plant biomass balance result (g/m2): " <<  cohortBiomassBalanceSum<<"\n";
    Rcout<<"Plant biomass balance components:\n";
    
    Rcout<<"  Structural balance (g/m2) "  <<round(sum(StandBiomassBalance(_,0)))<<" Labile balance (g/m2) "  <<round(sum(StandBiomassBalance(_,1))) <<"\n";
    Rcout<<"  Plant individual balance (g/m2) "  <<round(sum(StandBiomassBalance(_,2)))<<" Mortality loss (g/m2) "  <<round(sum(StandBiomassBalance(_,3))) <<"\n";
    
    printWaterBalanceResult(outputList, x,
                            initialPlantContent, initialSoilContent, initialSnowContent,
                            transpirationMode);
    
    if(error_occurence) {
      Rcout<< " ERROR: Calculations stopped because of numerical error: Revise parameters\n";
    }
  }
  
  return(outputList);
}
