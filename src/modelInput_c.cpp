#include "modelInput_c.h"
#include "control_c.h"
#include <RcppArmadillo.h>

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

  //Canopy parameters
  Rcpp::DataFrame canopyDF = Rcpp::as<Rcpp::DataFrame>(x["canopy"]);
  canopy.zlow = Rcpp::as< std::vector<double> >(canopyDF["zlow"]);
  canopy.zmid = Rcpp::as< std::vector<double> >(canopyDF["zmid"]);
  canopy.zup = Rcpp::as< std::vector<double> >(canopyDF["zup"]);
  canopy.LAIlive = Rcpp::as< std::vector<double> >(canopyDF["LAIlive"]);
  canopy.LAIexpanded = Rcpp::as< std::vector<double> >(canopyDF["LAIexpanded"]);
  canopy.LAIdead = Rcpp::as< std::vector<double> >(canopyDF["LAIdead"]);
  canopy.Tair = Rcpp::as< std::vector<double> >(canopyDF["Tair"]);
  canopy.Cair = Rcpp::as< std::vector<double> >(canopyDF["Cair"]);
  canopy.VPair = Rcpp::as< std::vector<double> >(canopyDF["VPair"]);
  
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

  //Interception parameters
  Rcpp::DataFrame intercDF = Rcpp::as<Rcpp::DataFrame>(x["paramsInterception"]);
  if(intercDF.containsElementNamed("LeafAngle")) paramsInterception.LeafAngle = Rcpp::as< std::vector<double> >(intercDF["LeafAngle"]);
  if(intercDF.containsElementNamed("LeafAngleSD")) paramsInterception.LeafAngleSD = Rcpp::as< std::vector<double> >(intercDF["LeafAngleSD"]);
  if(intercDF.containsElementNamed("Beta_p")) paramsInterception.Beta_p = Rcpp::as< std::vector<double> >(intercDF["Beta_p"]);
  if(intercDF.containsElementNamed("Beta_q")) paramsInterception.Beta_q = Rcpp::as< std::vector<double> >(intercDF["Beta_q"]);
  if(intercDF.containsElementNamed("ClumpingIndex")) paramsInterception.ClumpingIndex = Rcpp::as< std::vector<double> >(intercDF["ClumpingIndex"]);
  if(intercDF.containsElementNamed("kPAR")) paramsInterception.kPAR = Rcpp::as< std::vector<double> >(intercDF["kPAR"]);
  if(intercDF.containsElementNamed("kSWR")) paramsInterception.kSWR = Rcpp::as< std::vector<double> >(intercDF["kSWR"]);
  if(intercDF.containsElementNamed("alphaSWR")) paramsInterception.alphaSWR = Rcpp::as< std::vector<double> >(intercDF["alphaSWR"]);
  if(intercDF.containsElementNamed("alphaLWR")) paramsInterception.alphaLWR = Rcpp::as< std::vector<double> >(intercDF["alphaLWR"]);
  if(intercDF.containsElementNamed("g")) paramsInterception.g = Rcpp::as< std::vector<double> >(intercDF["g"]);
  
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

  //Water storage parameters
  Rcpp::DataFrame waterStorageDF = Rcpp::as<Rcpp::DataFrame>(x["paramsWaterStorage"]);
  paramsWaterStorage.maxFMC = Rcpp::as< std::vector<double> >(waterStorageDF["maxFMC"]);
  if(waterStorageDF.containsElementNamed("maxMCleaf")) paramsWaterStorage.maxMCleaf = Rcpp::as< std::vector<double> >(waterStorageDF["maxMCleaf"]);
  if(waterStorageDF.containsElementNamed("maxMCstem")) paramsWaterStorage.maxMCstem = Rcpp::as< std::vector<double> >(waterStorageDF["maxMCstem"]);
  paramsWaterStorage.LeafPI0 = Rcpp::as< std::vector<double> >(waterStorageDF["LeafPI0"]);
  paramsWaterStorage.LeafEPS = Rcpp::as< std::vector<double> >(waterStorageDF["LeafEPS"]);
  paramsWaterStorage.LeafAF = Rcpp::as< std::vector<double> >(waterStorageDF["LeafAF"]);
  paramsWaterStorage.Vleaf = Rcpp::as< std::vector<double> >(waterStorageDF["Vleaf"]);
  paramsWaterStorage.StemPI0 = Rcpp::as< std::vector<double> >(waterStorageDF["StemPI0"]);
  paramsWaterStorage.StemEPS = Rcpp::as< std::vector<double> >(waterStorageDF["StemEPS"]);
  paramsWaterStorage.StemAF = Rcpp::as< std::vector<double> >(waterStorageDF["StemAF"]);
  paramsWaterStorage.Vsapwood = Rcpp::as< std::vector<double> >(waterStorageDF["Vsapwood"]);
  
  //Transpiration parameters
  Rcpp::DataFrame transpDF = Rcpp::as<Rcpp::DataFrame>(x["paramsTranspiration"]);
  if(waterStorageDF.containsElementNamed("Tmax_LAI")) paramsTranspiration.Tmax_LAI = Rcpp::as< std::vector<double> >(transpDF["Tmax_LAI"]);
  if(waterStorageDF.containsElementNamed("Tmax_LAIsq")) paramsTranspiration.Tmax_LAIsq = Rcpp::as< std::vector<double> >(transpDF["Tmax_LAIsq"]);
  if(waterStorageDF.containsElementNamed("Psi_Extract")) paramsTranspiration.Psi_Extract = Rcpp::as< std::vector<double> >(transpDF["Psi_Extract"]);
  if(waterStorageDF.containsElementNamed("Exp_Extract")) paramsTranspiration.Exp_Extract = Rcpp::as< std::vector<double> >(transpDF["Exp_Extract"]);
  if(waterStorageDF.containsElementNamed("Gswmin")) paramsTranspiration.Gswmin = Rcpp::as< std::vector<double> >(transpDF["Gswmin"]);
  if(waterStorageDF.containsElementNamed("Gswmax")) paramsTranspiration.Gswmax = Rcpp::as< std::vector<double> >(transpDF["Gswmax"]);
  if(waterStorageDF.containsElementNamed("Gs_Toptim")) paramsTranspiration.Gs_Toptim = Rcpp::as< std::vector<double> >(transpDF["Gs_Toptim"]);
  if(waterStorageDF.containsElementNamed("Gs_Tsens")) paramsTranspiration.Gs_Tsens = Rcpp::as< std::vector<double> >(transpDF["Gs_Tsens"]);
  if(waterStorageDF.containsElementNamed("Gsw_AC_slope")) paramsTranspiration.Gsw_AC_slope = Rcpp::as< std::vector<double> >(transpDF["Gsw_AC_slope"]);
  if(waterStorageDF.containsElementNamed("Gs_P50")) paramsTranspiration.Gs_P50 = Rcpp::as< std::vector<double> >(transpDF["Gs_P50"]);
  if(waterStorageDF.containsElementNamed("Gs_slope")) paramsTranspiration.Gs_slope = Rcpp::as< std::vector<double> >(transpDF["Gs_slope"]);
  if(waterStorageDF.containsElementNamed("WUE")) paramsTranspiration.WUE = Rcpp::as< std::vector<double> >(transpDF["WUE"]);
  if(waterStorageDF.containsElementNamed("WUE_par")) paramsTranspiration.WUE_par = Rcpp::as< std::vector<double> >(transpDF["WUE_par"]);
  if(waterStorageDF.containsElementNamed("WUE_co2")) paramsTranspiration.WUE_co2 = Rcpp::as< std::vector<double> >(transpDF["WUE_co2"]);
  if(waterStorageDF.containsElementNamed("WUE_vpd")) paramsTranspiration.WUE_vpd = Rcpp::as< std::vector<double> >(transpDF["WUE_vpd"]);
  if(waterStorageDF.containsElementNamed("Vmax298")) paramsTranspiration.Vmax298 = Rcpp::as< std::vector<double> >(transpDF["Vmax298"]);
  if(waterStorageDF.containsElementNamed("Jmax298")) paramsTranspiration.Jmax298 = Rcpp::as< std::vector<double> >(transpDF["Jmax298"]);
  if(waterStorageDF.containsElementNamed("Kmax_stemxylem")) paramsTranspiration.Kmax_stemxylem = Rcpp::as< std::vector<double> >(transpDF["Kmax_stemxylem"]);
  if(waterStorageDF.containsElementNamed("Kmax_rootxylem")) paramsTranspiration.Kmax_rootxylem = Rcpp::as< std::vector<double> >(transpDF["Kmax_rootxylem"]);
  if(waterStorageDF.containsElementNamed("VCleaf_kmax")) paramsTranspiration.VCleaf_kmax = Rcpp::as< std::vector<double> >(transpDF["VCleaf_kmax"]);
  if(waterStorageDF.containsElementNamed("VCleafapo_kmax")) paramsTranspiration.VCleafapo_kmax = Rcpp::as< std::vector<double> >(transpDF["VCleafapo_kmax"]);
  if(waterStorageDF.containsElementNamed("VCleaf_slope")) paramsTranspiration.VCleaf_slope = Rcpp::as< std::vector<double> >(transpDF["VCleaf_slope"]);
  if(waterStorageDF.containsElementNamed("VCleaf_P50")) paramsTranspiration.VCleaf_P50 = Rcpp::as< std::vector<double> >(transpDF["VCleaf_P50"]);
  if(waterStorageDF.containsElementNamed("VCleaf_c")) paramsTranspiration.VCleaf_c = Rcpp::as< std::vector<double> >(transpDF["VCleaf_c"]);
  if(waterStorageDF.containsElementNamed("VCleaf_d")) paramsTranspiration.VCleaf_d = Rcpp::as< std::vector<double> >(transpDF["VCleaf_d"]);
  if(waterStorageDF.containsElementNamed("kleaf_symp")) paramsTranspiration.kleaf_symp = Rcpp::as< std::vector<double> >(transpDF["kleaf_symp"]);
  if(waterStorageDF.containsElementNamed("VCstem_kmax")) paramsTranspiration.VCstem_kmax = Rcpp::as< std::vector<double> >(transpDF["VCstem_kmax"]);
  if(waterStorageDF.containsElementNamed("VCstem_slope")) paramsTranspiration.VCstem_slope = Rcpp::as< std::vector<double> >(transpDF["VCstem_slope"]);
  if(waterStorageDF.containsElementNamed("VCstem_P50")) paramsTranspiration.VCstem_P50 = Rcpp::as< std::vector<double> >(transpDF["VCstem_P50"]);
  if(waterStorageDF.containsElementNamed("VCstem_c")) paramsTranspiration.VCstem_c = Rcpp::as< std::vector<double> >(transpDF["VCstem_c"]);
  if(waterStorageDF.containsElementNamed("VCstem_d")) paramsTranspiration.VCstem_d = Rcpp::as< std::vector<double> >(transpDF["VCstem_d"]);
  if(waterStorageDF.containsElementNamed("kstem_xylem")) paramsTranspiration.kstem_xylem = Rcpp::as< std::vector<double> >(transpDF["kstem_xylem"]);
  if(waterStorageDF.containsElementNamed("kstem_symp")) paramsTranspiration.kstem_symp = Rcpp::as< std::vector<double> >(transpDF["kstem_symp"]);
  if(waterStorageDF.containsElementNamed("VCroottot_kmax")) paramsTranspiration.VCroottot_kmax = Rcpp::as< std::vector<double> >(transpDF["VCroottot_kmax"]);
  if(waterStorageDF.containsElementNamed("VCroot_slope")) paramsTranspiration.VCroot_slope = Rcpp::as< std::vector<double> >(transpDF["VCroot_slope"]);
  if(waterStorageDF.containsElementNamed("VCroot_P50")) paramsTranspiration.VCroot_P50 = Rcpp::as< std::vector<double> >(transpDF["VCroot_P50"]);
  if(waterStorageDF.containsElementNamed("VCroot_c")) paramsTranspiration.VCroot_c = Rcpp::as< std::vector<double> >(transpDF["VCroot_c"]);
  if(waterStorageDF.containsElementNamed("VCroot_d")) paramsTranspiration.VCroot_d = Rcpp::as< std::vector<double> >(transpDF["VCroot_d"]);
  if(waterStorageDF.containsElementNamed("VGrhizotot_kmax")) paramsTranspiration.VGrhizotot_kmax = Rcpp::as< std::vector<double> >(transpDF["VGrhizotot_kmax"]);
  if(waterStorageDF.containsElementNamed("Plant_kmax")) paramsTranspiration.Plant_kmax = Rcpp::as< std::vector<double> >(transpDF["Plant_kmax"]);
  if(waterStorageDF.containsElementNamed("FR_leaf")) paramsTranspiration.FR_leaf = Rcpp::as< std::vector<double> >(transpDF["FR_leaf"]);
  if(waterStorageDF.containsElementNamed("FR_stem")) paramsTranspiration.FR_stem = Rcpp::as< std::vector<double> >(transpDF["FR_stem"]);
  if(waterStorageDF.containsElementNamed("FR_root")) paramsTranspiration.FR_root = Rcpp::as< std::vector<double> >(transpDF["FR_root"]);
  
  //Internal phenology variables
  Rcpp::DataFrame internalPhenoDF = Rcpp::as<Rcpp::DataFrame>(x["internalPhenology"]);
  internalPhenology.gdd = Rcpp::as< std::vector<double> >(internalPhenoDF["gdd"]);
  internalPhenology.sen = Rcpp::as< std::vector<double> >(internalPhenoDF["sen"]);
  internalPhenology.budFormation = Rcpp::as< std::vector<bool> >(internalPhenoDF["budFormation"]);
  internalPhenology.leafUnfolding = Rcpp::as< std::vector<bool> >(internalPhenoDF["leafUnfolding"]);
  internalPhenology.leafSenescence = Rcpp::as< std::vector<bool> >(internalPhenoDF["leafSenescence"]);
  internalPhenology.leafDormancy = Rcpp::as< std::vector<bool> >(internalPhenoDF["leafDormancy"]);
  internalPhenology.phi = Rcpp::as< std::vector<double> >(internalPhenoDF["phi"]);

  //Internal LAI distribution
  if(x.containsElementNamed("internalLAIDistribution")){
    Rcpp::List intLAIDist = x["internalLAIDistribution"];
    internalLAIDistribution.PrevLAIdead = Rcpp::as< std::vector<double> >(intLAIDist["PrevLAIdead"]);
    internalLAIDistribution.PrevLAIexpanded = Rcpp::as< std::vector<double> >(intLAIDist["PrevLAIexpanded"]);
    internalLAIDistribution.PARcohort = Rcpp::as< std::vector<double> >(intLAIDist["PARcohort"]);
    Rcpp::NumericMatrix liveMat = Rcpp::as<Rcpp::NumericMatrix>(intLAIDist["live"]);
    Rcpp::NumericMatrix expandedMat = Rcpp::as<Rcpp::NumericMatrix>(intLAIDist["expanded"]);
    Rcpp::NumericMatrix deadMat = Rcpp::as<Rcpp::NumericMatrix>(intLAIDist["dead"]);
    internalLAIDistribution.live = Rcpp::as<arma::mat>(liveMat);
    internalLAIDistribution.dead = Rcpp::as<arma::mat>(deadMat);
    internalLAIDistribution.expanded = Rcpp::as<arma::mat>(expandedMat);
  }
  //Internal water variables
  if(x.containsElementNamed("internalWater")) {
    Rcpp::DataFrame internalWaterDF = Rcpp::as<Rcpp::DataFrame>(x["internalWater"]);
    if(internalWaterDF.containsElementNamed("Einst")) internalWater.Einst = Rcpp::as< std::vector<double> >(internalWaterDF["Einst"]);
    if(internalWaterDF.containsElementNamed("Elim")) internalWater.Elim = Rcpp::as< std::vector<double> >(internalWaterDF["Elim"]);
    if(internalWaterDF.containsElementNamed("Emin_L")) internalWater.Emin_L = Rcpp::as< std::vector<double> >(internalWaterDF["Emin_L"]);
    if(internalWaterDF.containsElementNamed("Emin_S")) internalWater.Emin_S = Rcpp::as< std::vector<double> >(internalWaterDF["Emin_S"]);
    if(internalWaterDF.containsElementNamed("PlantPsi")) internalWater.PlantPsi = Rcpp::as< std::vector<double> >(internalWaterDF["PlantPsi"]);
    if(internalWaterDF.containsElementNamed("RootCrownPsi")) internalWater.RootCrownPsi = Rcpp::as< std::vector<double> >(internalWaterDF["RootCrownPsi"]);
    if(internalWaterDF.containsElementNamed("LeafPsi")) internalWater.LeafPsi = Rcpp::as< std::vector<double> >(internalWaterDF["LeafPsi"]);
    if(internalWaterDF.containsElementNamed("StemPsi")) internalWater.StemPsi = Rcpp::as< std::vector<double> >(internalWaterDF["StemPsi"]);
    if(internalWaterDF.containsElementNamed("LeafSympPsi")) internalWater.LeafSympPsi = Rcpp::as< std::vector<double> >(internalWaterDF["LeafSympPsi"]);
    if(internalWaterDF.containsElementNamed("StemSympPsi")) internalWater.StemSympPsi = Rcpp::as< std::vector<double> >(internalWaterDF["StemSympPsi"]);
    if(internalWaterDF.containsElementNamed("LeafPLC")) internalWater.LeafPLC = Rcpp::as< std::vector<double> >(internalWaterDF["LeafPLC"]);
    if(internalWaterDF.containsElementNamed("StemPLC")) internalWater.StemPLC = Rcpp::as< std::vector<double> >(internalWaterDF["StemPLC"]);
  }
  //Internal carbon variables
  if(x.containsElementNamed("internalCarbon")) {
    Rcpp::DataFrame internalCarbonDF = Rcpp::as<Rcpp::DataFrame>(x["internalCarbon"]);
    internalCarbon.sugarLeaf = Rcpp::as< std::vector<double> >(internalCarbonDF["sugarLeaf"]);
    internalCarbon.starchLeaf = Rcpp::as< std::vector<double> >(internalCarbonDF["starchLeaf"]);
    internalCarbon.sugarSapwood = Rcpp::as< std::vector<double> >(internalCarbonDF["sugarSapwood"]);
    internalCarbon.starchSapwood = Rcpp::as< std::vector<double> >(internalCarbonDF["starchSapwood"]);
  }
  
  //Internal mortality variables
  if(x.containsElementNamed("internalMortality")) {
    Rcpp::DataFrame internalMortalityDF = Rcpp::as<Rcpp::DataFrame>(x["internalMortality"]);
    internalMortality.N_dead = Rcpp::as< std::vector<double> >(internalMortalityDF["N_dead"]);
    internalMortality.N_starvation = Rcpp::as< std::vector<double> >(internalMortalityDF["N_starvation"]);
    internalMortality.N_dessication = Rcpp::as< std::vector<double> >(internalMortalityDF["N_dessication"]);
    internalMortality.N_burnt = Rcpp::as< std::vector<double> >(internalMortalityDF["N_burnt"]);
    internalMortality.N_resprouting_stumps = Rcpp::as< std::vector<double> >(internalMortalityDF["N_resprouting_stumps"]);
    internalMortality.Cover_dead = Rcpp::as< std::vector<double> >(internalMortalityDF["Cover_dead"]);
    internalMortality.Cover_starvation = Rcpp::as< std::vector<double> >(internalMortalityDF["Cover_starvation"]);
    internalMortality.Cover_dessication = Rcpp::as< std::vector<double> >(internalMortalityDF["Cover_dessication"]);
    internalMortality.Cover_burnt = Rcpp::as< std::vector<double> >(internalMortalityDF["Cover_burnt"]);
    internalMortality.Cover_resprouting_stumps = Rcpp::as< std::vector<double> >(internalMortalityDF["Cover_resprouting_stumps"]);
    if(internalMortalityDF.containsElementNamed("Snag_smallbranches")) internalMortality.Snag_smallbranches = Rcpp::as< std::vector<double> >(internalMortalityDF["Snag_smallbranches"]);
    if(internalMortalityDF.containsElementNamed("Snag_largewood")) internalMortality.Snag_largewood = Rcpp::as< std::vector<double> >(internalMortalityDF["Snag_largewood"]);
  }
  
  //Internal allocation variables
  if(x.containsElementNamed("internalAllocation")) {
    Rcpp::DataFrame internalAllocationDF = Rcpp::as<Rcpp::DataFrame>(x["internalAllocation"]);
    internalAllocation.allocationTarget = Rcpp::as< std::vector<double> >(internalAllocationDF["allocationTarget"]);
    internalAllocation.leafAreaTarget = Rcpp::as< std::vector<double> >(internalAllocationDF["leafAreaTarget"]);
    internalAllocation.sapwoodAreaTarget = Rcpp::as< std::vector<double> >(internalAllocationDF["sapwoodAreaTarget"]);
    internalAllocation.fineRootBiomassTarget = Rcpp::as< std::vector<double> >(internalAllocationDF["fineRootBiomassTarget"]);
    internalAllocation.crownBudPercent = Rcpp::as< std::vector<double> >(internalAllocationDF["crownBudPercent"]);
  }
  
  //Internal snag variables
  if(x.containsElementNamed("internalSnags")) {
    Rcpp::DataFrame internalSnagDF = Rcpp::as<Rcpp::DataFrame>(x["internalSnags"]);
    internalSnags.Species = Rcpp::as< std::vector<std::string> >(internalSnagDF["Species"]);
    internalSnags.DBH = Rcpp::as< std::vector<double> >(internalSnagDF["DBH"]);
    internalSnags.Height = Rcpp::as< std::vector<double> >(internalSnagDF["Height"]);
    internalSnags.SmallBranches = Rcpp::as< std::vector<double> >(internalSnagDF["SmallBranches"]);
    internalSnags.LargeWood = Rcpp::as< std::vector<double> >(internalSnagDF["LargeWood"]);
  }
  
  //Internal litter variables
  if(x.containsElementNamed("internalLitter")) {
    Rcpp::DataFrame internalLitterDF = Rcpp::as<Rcpp::DataFrame>(x["internalLitter"]);
    internalLitter.Species = Rcpp::as< std::vector<std::string> >(internalLitterDF["Species"]);
    internalLitter.Leaves = Rcpp::as< std::vector<double> >(internalLitterDF["Leaves"]);
    internalLitter.Twigs = Rcpp::as< std::vector<double> >(internalLitterDF["Twigs"]);
    internalLitter.SmallBranches = Rcpp::as< std::vector<double> >(internalLitterDF["SmallBranches"]);
    internalLitter.LargeWood = Rcpp::as< std::vector<double> >(internalLitterDF["LargeWood"]);
    internalLitter.CoarseRoots = Rcpp::as< std::vector<double> >(internalLitterDF["CoarseRoots"]);
    internalLitter.FineRoots = Rcpp::as< std::vector<double> >(internalLitterDF["FineRoots"]);
  }
  //Internal SOC variables
  if(x.containsElementNamed("internalSOC")) {
    Rcpp::NumericVector internalSOCDF = Rcpp::as<Rcpp::NumericVector>(x["internalSOC"]);
    internalSOC.SurfaceMetabolic = Rcpp::as<double>(internalSOCDF["SurfaceMetabolic"]);
    internalSOC.SoilMetabolic = Rcpp::as<double>(internalSOCDF["SoilMetabolic"]);
    internalSOC.SurfaceActive = Rcpp::as<double>(internalSOCDF["SurfaceActive"]);
    internalSOC.SoilActive = Rcpp::as<double>(internalSOCDF["SoilActive"]);
    internalSOC.SurfaceSlow = Rcpp::as<double>(internalSOCDF["SurfaceSlow"]);
    internalSOC.SoilSlow = Rcpp::as<double>(internalSOCDF["SoilSlow"]);
    internalSOC.SoilPassive = Rcpp::as<double>(internalSOCDF["SoilPassive"]);
  }
  
  //Internal FCCS variables
  if(x.containsElementNamed("internalFCCS")) {
    Rcpp::DataFrame fccsDF = Rcpp::as<Rcpp::DataFrame>(x["internalFCCS"]);
    if(fccsDF.nrows()>0) {
      internalFCCS.w = Rcpp::as< std::vector<double> >(fccsDF["w"]);
      internalFCCS.cover = Rcpp::as< std::vector<double> >(fccsDF["cover"]);
      internalFCCS.hbc = Rcpp::as< std::vector<double> >(fccsDF["hbc"]);
      internalFCCS.htc = Rcpp::as< std::vector<double> >(fccsDF["htc"]);
      internalFCCS.habc = Rcpp::as< std::vector<double> >(fccsDF["habc"]);
      internalFCCS.hatc = Rcpp::as< std::vector<double> >(fccsDF["hatc"]);
      internalFCCS.delta = Rcpp::as< std::vector<double> >(fccsDF["delta"]);
      internalFCCS.rhob = Rcpp::as< std::vector<double> >(fccsDF["rhob"]);
      internalFCCS.rhop = Rcpp::as< std::vector<double> >(fccsDF["rhop"]);
      internalFCCS.PV = Rcpp::as< std::vector<double> >(fccsDF["PV"]);
      internalFCCS.beta = Rcpp::as< std::vector<double> >(fccsDF["beta"]);
      internalFCCS.betarel = Rcpp::as< std::vector<double> >(fccsDF["betarel"]);
      internalFCCS.etabetarel = Rcpp::as< std::vector<double> >(fccsDF["etabetarel"]);
      internalFCCS.sigma = Rcpp::as< std::vector<double> >(fccsDF["sigma"]);
      internalFCCS.pDead = Rcpp::as< std::vector<double> >(fccsDF["pDead"]);
      internalFCCS.FAI = Rcpp::as< std::vector<double> >(fccsDF["FAI"]);
      internalFCCS.h = Rcpp::as< std::vector<double> >(fccsDF["h"]);
      internalFCCS.RV = Rcpp::as< std::vector<double> >(fccsDF["RV"]);
      internalFCCS.MinFMC = Rcpp::as< std::vector<double> >(fccsDF["MinFMC"]);
      internalFCCS.MaxFMC = Rcpp::as< std::vector<double> >(fccsDF["MaxFMC"]);
      internalFCCS.ActFMC = Rcpp::as< std::vector<double> >(fccsDF["ActFMC"]);
    }
  }
  
  if(x.containsElementNamed("version")) version = Rcpp::as<std::string>(x["version"]);
}