#include <RcppArmadillo.h>
#include "fuelmoisture_c.h"
#include "communication_structures.h"
#include "decomposition_c.h"
using namespace Rcpp;


//' @param x Soil water pH (0-14)
//' @param pool String indicating the decomposition pool
//' 
//' @rdname decomposition_annualLitterDecompositionRate
// [[Rcpp::export("decomposition_pHEffect")]]
double pHEffect(double x, String pool) {
  std::string pool_c = pool.get_cstring();
  return(pHEffect_c(x, pool_c));
}

// [[Rcpp::export(".decomposition_DAYCENTsnags")]]
NumericVector DAYCENTsnags(DataFrame snags, 
                           NumericVector baseAnnualRates, DataFrame paramsLitterDecomposition,
                           double airTemperature, double airRelativeHumidity, 
                           double tstep = 1.0) {
  
  
  Rcpp::CharacterVector Species_SN = snags["Species"];
  Rcpp::NumericVector DBH_SN = snags["DBH"];
  Rcpp::NumericVector Height_SN = snags["Height"];
  Rcpp::NumericVector DeathAge_SN = snags["DeathAge"];
  Rcpp::NumericVector SmallBranches_SN = snags["SmallBranches"];
  Rcpp::NumericVector LargeWood_SN = snags["LargeWood"];
  
  InternalSnags internalSnags;
  internalSnags.Species = Rcpp::as< std::vector<std::string> >(Species_SN);
  internalSnags.DBH = Rcpp::as< std::vector<double> >(DBH_SN);
  internalSnags.Height = Rcpp::as< std::vector<double> >(Height_SN);
  internalSnags.DeathAge = Rcpp::as< std::vector<double> >(DeathAge_SN);
  internalSnags.SmallBranches = Rcpp::as< std::vector<double> >(SmallBranches_SN);
  internalSnags.LargeWood = Rcpp::as< std::vector<double> >(LargeWood_SN);
  
  LitterDecompositionParams paramsLitterDecomposition_c;
  if(paramsLitterDecomposition.containsElementNamed("Species")) paramsLitterDecomposition_c.Species = Rcpp::as< std::vector<std::string> >(paramsLitterDecomposition["Species"]);
  if(paramsLitterDecomposition.containsElementNamed("LeafLignin")) paramsLitterDecomposition_c.LeafLignin = Rcpp::as< std::vector<double> >(paramsLitterDecomposition["LeafLignin"]);
  if(paramsLitterDecomposition.containsElementNamed("WoodLignin")) paramsLitterDecomposition_c.WoodLignin = Rcpp::as< std::vector<double> >(paramsLitterDecomposition["WoodLignin"]);
  if(paramsLitterDecomposition.containsElementNamed("FineRootLignin")) paramsLitterDecomposition_c.FineRootLignin = Rcpp::as< std::vector<double> >(paramsLitterDecomposition["FineRootLignin"]);
  if(paramsLitterDecomposition.containsElementNamed("Nleaf")) paramsLitterDecomposition_c.Nleaf = Rcpp::as< std::vector<double> >(paramsLitterDecomposition["Nleaf"]);
  if(paramsLitterDecomposition.containsElementNamed("Nsapwood")) paramsLitterDecomposition_c.Nsapwood = Rcpp::as< std::vector<double> >(paramsLitterDecomposition["Nsapwood"]);
  if(paramsLitterDecomposition.containsElementNamed("Nfineroot")) paramsLitterDecomposition_c.Nfineroot = Rcpp::as< std::vector<double> >(paramsLitterDecomposition["Nfineroot"]);
  
  DecompositionAnnualBaseRates baseAnnualRates_c;
  baseAnnualRates_c.SurfaceMetabolic = baseAnnualRates["SurfaceMetabolic"];
  baseAnnualRates_c.SoilMetabolic = baseAnnualRates["SoilMetabolic"];
  baseAnnualRates_c.Leaves = baseAnnualRates["Leaves"];
  baseAnnualRates_c.FineRoots = baseAnnualRates["FineRoots"];
  baseAnnualRates_c.Twigs = baseAnnualRates["Twigs"];
  baseAnnualRates_c.SmallBranches = baseAnnualRates["SmallBranches"];
  baseAnnualRates_c.LargeWood = baseAnnualRates["LargeWood"];
  baseAnnualRates_c.CoarseRoots = baseAnnualRates["CoarseRoots"];
  baseAnnualRates_c.SurfaceActive = baseAnnualRates["SurfaceActive"];
  baseAnnualRates_c.SurfaceSlow = baseAnnualRates["SurfaceSlow"];
  baseAnnualRates_c.SoilActive = baseAnnualRates["SoilActive"];
  baseAnnualRates_c.SoilSlow = baseAnnualRates["SoilSlow"];
  baseAnnualRates_c.SoilPassive = baseAnnualRates["SoilPassive"];  
  SnagDecomposition_COMM sdo;
  DAYCENTsnagsInner_c(sdo,
                      internalSnags, paramsLitterDecomposition_c,
                      baseAnnualRates_c,
                      airTemperature, airRelativeHumidity,
                      tstep);
  
  for(int c = 0;c < Species_SN.size(); c++) {
    Species_SN[c] = internalSnags.Species[c];
    DBH_SN[c] = internalSnags.DBH[c];
    Height_SN[c] = internalSnags.Height[c];
    DeathAge_SN[c] = internalSnags.DeathAge[c];
    SmallBranches_SN[c] = internalSnags.SmallBranches[c];
    LargeWood_SN[c] = internalSnags.LargeWood[c];
  }
  
  NumericVector output = NumericVector::create(_["transfer_surface_active"] = sdo.transfer_surface_active,
                                               _["transfer_surface_slow"] = sdo.transfer_surface_slow,
                                               _["flux_respiration"] = sdo.flux_respiration);
  return(output);
}



// [[Rcpp::export(".decomposition_DAYCENTlitter")]]
NumericVector DAYCENTlitter(DataFrame litter, DataFrame paramsLitterDecomposition,
                            NumericVector baseAnnualRates,
                            double sand, double clay, double soilTemperature, double soilMoisture, double soilPH, 
                            double soilO2 = 1.0, double cultfac = 1.0,
                            double tstep = 1.0) {

  LitterDecompositionParams paramsLitterDecomposition_c;
  if(paramsLitterDecomposition.containsElementNamed("Species")) paramsLitterDecomposition_c.Species = Rcpp::as< std::vector<std::string> >(paramsLitterDecomposition["Species"]);
  if(paramsLitterDecomposition.containsElementNamed("LeafLignin")) paramsLitterDecomposition_c.LeafLignin = Rcpp::as< std::vector<double> >(paramsLitterDecomposition["LeafLignin"]);
  if(paramsLitterDecomposition.containsElementNamed("WoodLignin")) paramsLitterDecomposition_c.WoodLignin = Rcpp::as< std::vector<double> >(paramsLitterDecomposition["WoodLignin"]);
  if(paramsLitterDecomposition.containsElementNamed("FineRootLignin")) paramsLitterDecomposition_c.FineRootLignin = Rcpp::as< std::vector<double> >(paramsLitterDecomposition["FineRootLignin"]);
  if(paramsLitterDecomposition.containsElementNamed("Nleaf")) paramsLitterDecomposition_c.Nleaf = Rcpp::as< std::vector<double> >(paramsLitterDecomposition["Nleaf"]);
  if(paramsLitterDecomposition.containsElementNamed("Nsapwood")) paramsLitterDecomposition_c.Nsapwood = Rcpp::as< std::vector<double> >(paramsLitterDecomposition["Nsapwood"]);
  if(paramsLitterDecomposition.containsElementNamed("Nfineroot")) paramsLitterDecomposition_c.Nfineroot = Rcpp::as< std::vector<double> >(paramsLitterDecomposition["Nfineroot"]);
  
  
  DecompositionAnnualBaseRates baseAnnualRates_c;
  baseAnnualRates_c.SurfaceMetabolic = baseAnnualRates["SurfaceMetabolic"];
  baseAnnualRates_c.SoilMetabolic = baseAnnualRates["SoilMetabolic"];
  baseAnnualRates_c.Leaves = baseAnnualRates["Leaves"];
  baseAnnualRates_c.FineRoots = baseAnnualRates["FineRoots"];
  baseAnnualRates_c.Twigs = baseAnnualRates["Twigs"];
  baseAnnualRates_c.SmallBranches = baseAnnualRates["SmallBranches"];
  baseAnnualRates_c.LargeWood = baseAnnualRates["LargeWood"];
  baseAnnualRates_c.CoarseRoots = baseAnnualRates["CoarseRoots"];
  baseAnnualRates_c.SurfaceActive = baseAnnualRates["SurfaceActive"];
  baseAnnualRates_c.SurfaceSlow = baseAnnualRates["SurfaceSlow"];
  baseAnnualRates_c.SoilActive = baseAnnualRates["SoilActive"];
  baseAnnualRates_c.SoilSlow = baseAnnualRates["SoilSlow"];
  baseAnnualRates_c.SoilPassive = baseAnnualRates["SoilPassive"];  
  
  Rcpp::CharacterVector Species_LI = litter["Species"];
  Rcpp::NumericVector Leaves_LI = litter["Leaves"];
  Rcpp::NumericVector Twigs_LI = litter["Twigs"];
  Rcpp::NumericVector SmallBranches_LI = litter["SmallBranches"];
  Rcpp::NumericVector LargeWood_LI = litter["LargeWood"];
  Rcpp::NumericVector CoarseRoots_LI = litter["CoarseRoots"];
  Rcpp::NumericVector FineRoots_LI = litter["FineRoots"];
  InternalLitter internalLitter;
  internalLitter.Species = Rcpp::as< std::vector<std::string> >(Species_LI);
  internalLitter.Leaves = Rcpp::as< std::vector<double> >(Leaves_LI);
  internalLitter.Twigs = Rcpp::as< std::vector<double> >(Twigs_LI);
  internalLitter.SmallBranches = Rcpp::as< std::vector<double> >(SmallBranches_LI);
  internalLitter.LargeWood = Rcpp::as< std::vector<double> >(LargeWood_LI);
  internalLitter.CoarseRoots = Rcpp::as< std::vector<double> >(CoarseRoots_LI);
  internalLitter.FineRoots = Rcpp::as< std::vector<double> >(FineRoots_LI);
  
  LitterDecomposition_COMM ldo;
  DAYCENTlitterInner_c(ldo, 
                       internalLitter, paramsLitterDecomposition_c, 
                       baseAnnualRates_c,
                       sand, clay, soilTemperature, soilMoisture, soilPH,
                       soilO2, cultfac,
                       tstep);
  NumericVector output = NumericVector::create(_["transfer_surface_active"] = ldo.transfer_surface_active,
                                               _["transfer_surface_slow"] = ldo.transfer_surface_slow,
                                               _["transfer_soil_active"] = ldo.transfer_soil_active,
                                               _["transfer_soil_slow"] = ldo.transfer_soil_slow,
                                               _["flux_respiration"] = ldo.flux_respiration);
  return(output);
}


//' DAYCENT decomposition
//' 
//' Function implementing a modification of the DAYCENT carbon decomposition model (Parton et al. 1988, 1993, 1998), inspired by the description given in Chapter 18 of Bonan (2019).
//' 
//' @param snags A data frame with initial dead standing (snag) cohort information (see \code{\link{growthInput}}).
//' @param litter A data frame with initial aboveground and belowground structural carbon pools corresponding to plant cohorts, in g C/m2  (see \code{\link{growthInput}}).
//' @param SOC A named numeric vector with initial metabolic, active, slow and passive carbon pools for surface and soil, in g C/m2  (see \code{\link{growthInput}}).
//' @param paramsLitterDecomposition A data frame of species-specific litter decomposition parameters (see \code{\link{growthInput}}).
//' @param baseAnnualRates A named vector of annual decomposition rates, in yr-1 (see \code{\link{defaultControl}}).
//' @param annualTurnoverRate Annual turnover rate, in yr-1  (see \code{\link{defaultControl}}).
//' @param sand,clay Soil texture (sand and sand) in percent volume (percent). 
//' @param environmentalConditions A data frame containing environmental conditions for each time step to simulate:
//'   \itemize{
//'     \item{\code{AirTempperature}: Mean air temperature (in Celsius).}
//'     \item{\code{AirRelativeHumidity}: Mean relative humidity (percent).}
//'     \item{\code{SoilTemperature}: Mean soil temperature (in Celsius).}
//'     \item{\code{SoilMoisture}: Mean soil moisture, relative to saturation (0-1).}
//'   }
//' @param litterProduction A data frame containing litter inputs corresponding to time steps:
//'   \itemize{
//'     \item{\code{Step}: Integer indicating the time step where litter is generated, corresponding to a row in \code{environmentalConditions}.}
//'     \item{\code{Species}: Mean relative humidity (percent).}
//'     \item{\code{Leaves}: Leaf litter production (g C·m-2).}
//'     \item{\code{Twigs}: Twig litter production (g C·m-2).}
//'     \item{\code{SmallBranches}: Small branch litter production (g C·m-2).}
//'     \item{\code{LargeWood}: Large wood litter production (g C·m-2).}
//'     \item{\code{CoarseRoots}: Coarse root litter production (g C·m-2).}
//'     \item{\code{FineRoots}: Fine root litter production (g C·m-2).}
//'   }
//' @param soilPH Soil pH (0-14).
//' @param soilO2 Soil oxygen factor (0-1).
//' @param cultfac Cultivation factor (0-1).
//' @param tstep Time step in days. By default, one day. For annual time steps, use \code{tstep = 365.25}.
//' 
//' @author 
//' Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
//' 
//' Roberto Molowny-Horas, CREAF
//' 
//' Inés Delsman Valderrama, CREAF
//' 
//' @details Species names must match between inputs \code{paramsLitterDecomposition}, \code{litter} and \code{litterProduction}.
//' 
//' \emph{IMPORTANT NOTE}: Decomposition functions modify the input data (i.e. \code{snags}, \code{litter} and/or \code{SOC}) according to decomposition rates and carbon transfer rates. When used as part of \code{\link{growth}} simulations,
//' soil physical and chemical parameters correspond to the uppermost soil layer.
//' 
//' @returns A list with two elements:
//'   \itemize{
//'     \item{\code{"HeterotrophicRespiration"}: A numeric vector with the heterotrophic respiration (g C·m-2) corresponding to the decomposition in each time step.}
//'     \item{\code{"DecompositionPools"}: A numeric matrix with the mass of different decomposition carbon pools, all in g C · m-2.}
//'   }
//' 
//' @references
//' Bonan, G. (2019). Climate change and terrestrial ecosystem modeling. Cambridge University Press, Cambridge, UK.
//' 
//' Parton WJ, Steward JWB, Cole CV (1988). Dynamics of C, N, P and S in grassland soils: a model. Biogeochemistry 5: 109-131.
//' 
//' Parton WJ, Scurlock JMO, Ojima DS, Gilmanov TG, Scoles RJ et al. (1993). Observations and modeling of biomass and soil organic matter dynamics for the grassland biome worldwide. Global Biogeochemical Cycles 7: 785-809.
//'  
//' Parton WJ, Hartman M, Ojima DS, Schimel D (1998). DAYCENT and its land surface submodel: Description and testing. Global and Planetary Change, 19: 35-48.
//' 
//' @seealso \code{\link{decomposition_temperatureEffect}}, \code{\link{growthInput}}, \code{\link{growth}}
//' 
//' @keywords internal
// [[Rcpp::export("decomposition_DAYCENT")]]
List DAYCENT(DataFrame snags, DataFrame litter, NumericVector SOC,
             DataFrame paramsLitterDecomposition,
             NumericVector baseAnnualRates, double annualTurnoverRate,
             DataFrame environmentalConditions,
             DataFrame litterProduction, 
             double sand, double clay,  double soilPH, 
             double soilO2 = 1.0, double cultfac = 1.0,
             double tstep = 1.0) {
  
  // Environmental vectors
  NumericVector airTemperature = environmentalConditions["AirTemperature"];
  NumericVector airRelativeHumidity = environmentalConditions["AirRelativeHumidity"];
  NumericVector soilTemperature = environmentalConditions["SoilTemperature"];
  NumericVector soilMoisture = environmentalConditions["SoilMoisture"];

  //Litter production vectors
  IntegerVector stepProduction = litterProduction["Step"];
  CharacterVector speciesProduction = litterProduction["Species"];
  NumericVector leafProduction = litterProduction["Leaves"];
  NumericVector twigProduction = litterProduction["Twigs"];
  NumericVector smallBranchProduction = litterProduction["SmallBranches"];
  NumericVector largeWoodProduction = litterProduction["LargeWood"];
  NumericVector coarseRootProduction = litterProduction["CoarseRoots"];
  NumericVector fineRootProduction = litterProduction["FineRoots"];
  
  // Build structure InternalLitter
  Rcpp::CharacterVector Species_LI = litter["Species"];
  Rcpp::NumericVector Leaves_LI = litter["Leaves"];
  Rcpp::NumericVector Twigs_LI = litter["Twigs"];
  Rcpp::NumericVector SmallBranches_LI = litter["SmallBranches"];
  Rcpp::NumericVector LargeWood_LI = litter["LargeWood"];
  Rcpp::NumericVector CoarseRoots_LI = litter["CoarseRoots"];
  Rcpp::NumericVector FineRoots_LI = litter["FineRoots"];
  InternalLitter internalLitter;
  internalLitter.Species = Rcpp::as< std::vector<std::string> >(Species_LI);
  internalLitter.Leaves = Rcpp::as< std::vector<double> >(Leaves_LI);
  internalLitter.Twigs = Rcpp::as< std::vector<double> >(Twigs_LI);
  internalLitter.SmallBranches = Rcpp::as< std::vector<double> >(SmallBranches_LI);
  internalLitter.LargeWood = Rcpp::as< std::vector<double> >(LargeWood_LI);
  internalLitter.CoarseRoots = Rcpp::as< std::vector<double> >(CoarseRoots_LI);
  internalLitter.FineRoots = Rcpp::as< std::vector<double> >(FineRoots_LI);
  
  // Build structure InternalSnags
  Rcpp::CharacterVector Species_SN = snags["Species"];
  Rcpp::NumericVector DBH_SN = snags["DBH"];
  Rcpp::NumericVector Height_SN = snags["Height"];
  Rcpp::NumericVector DeathAge_SN = snags["DeathAge"];
  Rcpp::NumericVector SmallBranches_SN = snags["SmallBranches"];
  Rcpp::NumericVector LargeWood_SN = snags["LargeWood"];
  
  InternalSnags internalSnags;
  internalSnags.Species = Rcpp::as< std::vector<std::string> >(Species_SN);
  internalSnags.DBH = Rcpp::as< std::vector<double> >(DBH_SN);
  internalSnags.Height = Rcpp::as< std::vector<double> >(Height_SN);
  internalSnags.DeathAge = Rcpp::as< std::vector<double> >(DeathAge_SN);
  internalSnags.SmallBranches = Rcpp::as< std::vector<double> >(SmallBranches_SN);
  internalSnags.LargeWood = Rcpp::as< std::vector<double> >(LargeWood_SN);
  
  // Build structure InternalSOC
  InternalSOC internalSOC;
  internalSOC.SurfaceMetabolic = Rcpp::as<double>(SOC["SurfaceMetabolic"]);
  internalSOC.SoilMetabolic = Rcpp::as<double>(SOC["SoilMetabolic"]);
  internalSOC.SurfaceActive = Rcpp::as<double>(SOC["SurfaceActive"]);
  internalSOC.SoilActive = Rcpp::as<double>(SOC["SoilActive"]);
  internalSOC.SurfaceSlow = Rcpp::as<double>(SOC["SurfaceSlow"]);
  internalSOC.SoilSlow = Rcpp::as<double>(SOC["SoilSlow"]);
  internalSOC.SoilPassive = Rcpp::as<double>(SOC["SoilPassive"]);
  
  
  // Build decomposition annual base rates
  DecompositionAnnualBaseRates baseAnnualRates_c;
  baseAnnualRates_c.SurfaceMetabolic = baseAnnualRates["SurfaceMetabolic"];
  baseAnnualRates_c.SoilMetabolic = baseAnnualRates["SoilMetabolic"];
  baseAnnualRates_c.Leaves = baseAnnualRates["Leaves"];
  baseAnnualRates_c.FineRoots = baseAnnualRates["FineRoots"];
  baseAnnualRates_c.Twigs = baseAnnualRates["Twigs"];
  baseAnnualRates_c.SmallBranches = baseAnnualRates["SmallBranches"];
  baseAnnualRates_c.LargeWood = baseAnnualRates["LargeWood"];
  baseAnnualRates_c.CoarseRoots = baseAnnualRates["CoarseRoots"];
  baseAnnualRates_c.SurfaceActive = baseAnnualRates["SurfaceActive"];
  baseAnnualRates_c.SurfaceSlow = baseAnnualRates["SurfaceSlow"];
  baseAnnualRates_c.SoilActive = baseAnnualRates["SoilActive"];
  baseAnnualRates_c.SoilSlow = baseAnnualRates["SoilSlow"];
  baseAnnualRates_c.SoilPassive = baseAnnualRates["SoilPassive"];  
  
  // Copy species decomposition params
  LitterDecompositionParams paramsLitterDecomposition_c;
  if(paramsLitterDecomposition.containsElementNamed("Species")) paramsLitterDecomposition_c.Species = Rcpp::as< std::vector<std::string> >(paramsLitterDecomposition["Species"]);
  if(paramsLitterDecomposition.containsElementNamed("LeafLignin")) paramsLitterDecomposition_c.LeafLignin = Rcpp::as< std::vector<double> >(paramsLitterDecomposition["LeafLignin"]);
  if(paramsLitterDecomposition.containsElementNamed("WoodLignin")) paramsLitterDecomposition_c.WoodLignin = Rcpp::as< std::vector<double> >(paramsLitterDecomposition["WoodLignin"]);
  if(paramsLitterDecomposition.containsElementNamed("FineRootLignin")) paramsLitterDecomposition_c.FineRootLignin = Rcpp::as< std::vector<double> >(paramsLitterDecomposition["FineRootLignin"]);
  if(paramsLitterDecomposition.containsElementNamed("Nleaf")) paramsLitterDecomposition_c.Nleaf = Rcpp::as< std::vector<double> >(paramsLitterDecomposition["Nleaf"]);
  if(paramsLitterDecomposition.containsElementNamed("Nsapwood")) paramsLitterDecomposition_c.Nsapwood = Rcpp::as< std::vector<double> >(paramsLitterDecomposition["Nsapwood"]);
  if(paramsLitterDecomposition.containsElementNamed("Nfineroot")) paramsLitterDecomposition_c.Nfineroot = Rcpp::as< std::vector<double> >(paramsLitterDecomposition["Nfineroot"]);
  
  Decomposition_COMM dec_com(7);
  int nsteps = environmentalConditions.nrow();
  int nrow_prod = stepProduction.size();
  
  NumericVector heterotrophicRespiration(nsteps, 0.0);
  NumericMatrix decompositionPools(nsteps, 10);
  decompositionPools.attr("dimnames") = List::create(seq(1, nsteps), 
                          CharacterVector::create("SurfaceSnags","SurfaceLitter", "SoilLitter", "SurfaceMetabolic", "SoilMetabolic",
                                                  "SurfaceActive", "SoilActive", "SurfaceSlow", "SoilSlow", "SoilPassive"));
  
  // Rcout<<"\nStarting simulation\n";
  for(int s = 0; s< nsteps; s++ ) {
    //Add Litter production
    for(int i=0;i< nrow_prod;i++) {
      if(stepProduction[i] == (s+1)) {
        std::string species_litter = as<std::string>(speciesProduction[i]);
        addLeafTwigLitter_c(species_litter, 
                            leafProduction[i], twigProduction[i],
                            internalLitter, 
                            paramsLitterDecomposition_c,
                            internalSOC);
        addSmallBranchLitter_c(species_litter, smallBranchProduction[i], 
                               internalLitter);
        addLargeWoodLitter_c(species_litter, largeWoodProduction[i], 
                             internalLitter);
        addCoarseRootLitter_c(species_litter, coarseRootProduction[i], 
                              internalLitter);
        addFineRootLitter_c(species_litter, fineRootProduction[i], 
                            internalLitter, 
                            paramsLitterDecomposition_c,
                            internalSOC);
      }
    }
    // decomposition
    // Rcout<<"Entering decomposition: \n";
    heterotrophicRespiration[s]  = DAYCENTInner_c(dec_com,
                                                  internalSnags, internalLitter, internalSOC,
                                                  paramsLitterDecomposition_c,
                                                  baseAnnualRates_c, annualTurnoverRate,
                                                  airTemperature[s], airRelativeHumidity[s],  
                                                  sand, clay, soilTemperature[s], soilMoisture[s], soilPH, 
                                                  soilO2, cultfac,
                                                  tstep);
    
    // Rcout<<"Storing state \n";
    decompositionPools(s,0) = std::accumulate(internalSnags.SmallBranches.begin(), internalSnags.SmallBranches.end(),0.0) + std::accumulate(internalSnags.LargeWood.begin(), internalSnags.LargeWood.end(),0.0);
    decompositionPools(s,1) = std::accumulate(internalLitter.Leaves.begin(), internalLitter.Leaves.end(),0.0);
    decompositionPools(s,1) += std::accumulate(internalLitter.Twigs.begin(), internalLitter.Twigs.end(),0.0);
    decompositionPools(s,1) += std::accumulate(internalLitter.SmallBranches.begin(), internalLitter.SmallBranches.end(),0.0);
    decompositionPools(s,1) += std::accumulate(internalLitter.LargeWood.begin(), internalLitter.LargeWood.end(),0.0);
    decompositionPools(s,2) = std::accumulate(internalLitter.FineRoots.begin(), internalLitter.FineRoots.end(),0.0) + std::accumulate(internalLitter.CoarseRoots.begin(), internalLitter.CoarseRoots.end(),0.0) ;
    decompositionPools(s,3) = internalSOC.SurfaceMetabolic;
    decompositionPools(s,4) = internalSOC.SoilMetabolic;
    decompositionPools(s,5) = internalSOC.SurfaceActive;
    decompositionPools(s,6) = internalSOC.SoilActive;
    decompositionPools(s,7) = internalSOC.SurfaceSlow;
    decompositionPools(s,8) = internalSOC.SoilSlow;
    decompositionPools(s,9) = internalSOC.SoilPassive;
  }
  
  // Update inputs for snag stock
  for(int c = 0;c < Species_SN.size(); c++) {
    Species_SN[c] = internalSnags.Species[c];
    DBH_SN[c] = internalSnags.DBH[c];
    Height_SN[c] = internalSnags.Height[c];
    DeathAge_SN[c] = internalSnags.DeathAge[c];
    SmallBranches_SN[c] = internalSnags.SmallBranches[c];
    LargeWood_SN[c] = internalSnags.LargeWood[c];
  }
  // Update inputs for litter stock
  for(int c = 0;c < Species_LI.size(); c++) {
    Species_LI[c] = internalLitter.Species[c];
    Leaves_LI[c] = internalLitter.Leaves[c];
    Twigs_LI[c] = internalLitter.Twigs[c];
    SmallBranches_LI[c] = internalLitter.SmallBranches[c];
    LargeWood_LI[c] = internalLitter.LargeWood[c];
    CoarseRoots_LI[c] = internalLitter.CoarseRoots[c];
    FineRoots_LI[c] = internalLitter.FineRoots[c];
  }
  // Update inputs for SOC
  SOC["SurfaceMetabolic"] = internalSOC.SurfaceMetabolic;
  SOC["SoilMetabolic"] = internalSOC.SoilMetabolic;
  SOC["SurfaceActive"] = internalSOC.SurfaceActive;
  SOC["SoilActive"] = internalSOC.SoilActive;
  SOC["SurfaceSlow"] = internalSOC.SurfaceSlow;
  SOC["SoilSlow"] = internalSOC.SoilSlow;
  SOC["SoilPassive"] = internalSOC.SoilPassive;
  
  //Assemble output
  List l = List::create(_["HeterotrophicRespiration"] = heterotrophicRespiration,
                        _["DecompositionPools"] = decompositionPools);
  return(l);
}
